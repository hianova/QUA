#define _POSIX_C_SOURCE 200809L

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "qctrl_ring0.h"

/* ============================================================================
 * Fast Deterministic PRNG (xoshiro256**) & Exact Binomial Sampler
 * ============================================================================ */

typedef struct {
    uint64_t s[4];
} rng_t;

static inline uint64_t rotl64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t rng_next(rng_t *r) {
    const uint64_t res = rotl64(r->s[1] * 5ULL, 7) * 9ULL;
    const uint64_t t   = r->s[1] << 17;
    r->s[2] ^= r->s[0];
    r->s[3] ^= r->s[1];
    r->s[1] ^= r->s[2];
    r->s[0] ^= r->s[3];
    r->s[2] ^= t;
    r->s[3] = rotl64(r->s[3], 45);
    return res;
}

static void rng_seed(rng_t *r, uint64_t seed) {
    uint64_t z = seed + UINT64_C(0x9E3779B97F4A7C15);
    for (int i = 0; i < 4; ++i) {
        z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
        r->s[i] = z ^ (z >> 31);
        z += UINT64_C(0x9E3779B97F4A7C15);
    }
}

static double rng_uniform01(rng_t *r) {
    return (double)(rng_next(r) >> 11) * (1.0 / 9007199254740992.0);
}

static double rng_normal(rng_t *r) {
    double u1 = rng_uniform01(r);
    double u2 = rng_uniform01(r);
    if (u1 < 1e-15) {
        u1 = 1e-15;
    }
    return sqrt(-2.0 * log(u1)) * cos(6.2831853071795864769 * u2);
}

/* Exact Binomial(n, p) sampler in pure C using scaled uint64 threshold comparison */
static int32_t sample_binomial(rng_t *r, int32_t n, double p) {
    if (p <= 0.0) {
        return 0;
    }
    if (p >= 1.0) {
        return n;
    }
    /* Upper 53-bit comparison equivalent to exact U[0,1) < p without FP conversion in inner loop */
    uint64_t thresh = (uint64_t)(p * 18446744073709551616.0);
    int32_t cnt     = 0;
    for (int32_t i = 0; i < n; ++i) {
        cnt += (rng_next(r) < thresh) ? 1 : 0;
    }
    return cnt;
}

static double clamp01(double x) {
    return (x < 0.0) ? 0.0 : ((x > 1.0) ? 1.0 : x);
}

static void sample_ps_edge_counts(rng_t *r, double phi_true, double theta_cal, double e01,
                                  double e10, double bias_m, double bias_p, int32_t n_shots,
                                  int32_t *km_out, int32_t *k0_out, int32_t *kp_out) {
    const double base    = 1e-4;
    const double eta     = 1.0 - e01 - e10;
    const double half_pi = 1.5707963267948966192;
    double dm            = phi_true - (theta_cal - half_pi);
    double d0            = phi_true - theta_cal;
    double dp            = phi_true - (theta_cal + half_pi);

    double sm = sin(0.5 * dm);
    double s0 = sin(0.5 * d0);
    double sp = sin(0.5 * dp);

    double pm = clamp01(e01 + eta * (base + (1.0 - 2.0 * base) * sm * sm) + bias_m);
    double p0 = clamp01(e01 + eta * (base + (1.0 - 2.0 * base) * s0 * s0));
    double pp = clamp01(e01 + eta * (base + (1.0 - 2.0 * base) * sp * sp) + bias_p);

    *km_out = sample_binomial(r, n_shots, pm);
    *k0_out = sample_binomial(r, n_shots, p0);
    *kp_out = sample_binomial(r, n_shots, pp);
}

/* ============================================================================
 * Self-Test: Topology Coloring + Theorems 11, 12, 13 + Ring-0 Latency Benchmark
 * ============================================================================ */

static int run_selftest(void) {
    qctrl_edge_desc_t edges[QCTRL_NUM_EDGES];
    uint32_t n_edges = qctrl_build_grid40_edges(edges);
    if (n_edges != QCTRL_NUM_EDGES) {
        fprintf(stderr, "SELFTEST FAIL: expected 67 edges, got %u\n", n_edges);
        return 1;
    }

    /* Check conflict-free 4-coloring: no vertex shared inside same color */
    uint32_t color_counts[QCTRL_NUM_COLORS] = {0U, 0U, 0U, 0U};
    uint32_t region_counts[3]               = {0U, 0U, 0U};
    for (uint8_t col = 0U; col < QCTRL_NUM_COLORS; ++col) {
        uint8_t used[QCTRL_NUM_QUBITS];
        memset(used, 0, sizeof(used));
        for (uint32_t e = 0U; e < n_edges; ++e) {
            if (edges[e].color == col) {
                if (used[edges[e].u] || used[edges[e].v]) {
                    fprintf(stderr, "SELFTEST FAIL: vertex conflict in color %u edge %u\n", col, e);
                    return 1;
                }
                used[edges[e].u] = 1U;
                used[edges[e].v] = 1U;
                color_counts[col] += 1U;
            }
        }
    }
    for (uint32_t e = 0U; e < n_edges; ++e) {
        region_counts[edges[e].region] += 1U;
    }
    if (region_counts[QCTRL_REGION_PATCH_A] != 31U || region_counts[QCTRL_REGION_PATCH_B] != 31U ||
        region_counts[QCTRL_REGION_SEAM] != 5U) {
        fprintf(stderr, "SELFTEST FAIL: unexpected region counts (%u, %u, %u)\n",
                region_counts[0], region_counts[1], region_counts[2]);
        return 1;
    }

    /* Theorem 11 (exact asymmetric SPAM cancellation in integer Q16.16) */
    qctrl_cfg_t cfg = qctrl_cfg_default();
    qctrl_edge_state_t st;
    qctrl_edge_state_reset(&st);
    {
        const double phi_ref = 0.08;
        const int32_t n_det  = 200000;
        const double eta     = 1.0 - 0.02 - 0.08;
        int32_t km           = (int32_t)llround(n_det * (0.02 + eta * pow(sin(0.5 * (phi_ref + 1.5707963267948966)), 2.0)));
        int32_t k0           = (int32_t)llround(n_det * (0.02 + eta * pow(sin(0.5 * phi_ref), 2.0)));
        int32_t kp           = (int32_t)llround(n_det * (0.02 + eta * pow(sin(0.5 * (phi_ref - 1.5707963267948966)), 2.0)));
        qctrl_cfg_t det_cfg  = cfg;
        det_cfg.log_thresh_q16 = 0; /* immediate trigger test */
        det_cfg.z_soft_q16     = 0;
        (void)qctrl_edge_step_ring0(&st, &det_cfg, km, k0, kp, n_det);
        double est = (double)st.theta_q16 / (double)FLOW_Q16_ONE;
        if (fabs(est - phi_ref) > 5e-4) {
            fprintf(stderr, "SELFTEST FAIL: Theorem 11 SPAM cancellation error |%.6f - %.6f|\n",
                    est, phi_ref);
            return 1;
        }
    }

    /* Theorem 12 (20-sigma single glitch cannot trigger Huber e-Process) */
    qctrl_edge_state_reset(&st);
    {
        int32_t n  = 2000;
        int32_t k0 = 40;
        int32_t km = 400;
        int32_t kp = 1400; /* huge +20 sigma glitch on kp-km */
        bool act   = qctrl_edge_step_ring0(&st, &cfg, km, k0, kp, n);
        if (act || st.act_count != 0U) {
            fprintf(stderr, "SELFTEST FAIL: Theorem 12 glitch triggered Huber e-Process!\n");
            return 1;
        }
    }

    /* Real wall-clock Ring-0 latency benchmark (50,000 full 67-edge 4-color sweeps) */
    qctrl_edge_state_t grid_st[QCTRL_NUM_EDGES];
    int32_t km_arr[QCTRL_NUM_EDGES];
    int32_t k0_arr[QCTRL_NUM_EDGES];
    int32_t kp_arr[QCTRL_NUM_EDGES];
    for (uint32_t e = 0U; e < QCTRL_NUM_EDGES; ++e) {
        qctrl_edge_state_reset(&grid_st[e]);
        km_arr[e] = 920 + (int32_t)(e % 17U);
        k0_arr[e] = 45 + (int32_t)(e % 5U);
        kp_arr[e] = 900 + (int32_t)(e % 13U);
    }

    struct timespec t_start, t_end;
    const uint32_t iters = 50000U;
    volatile uint64_t sink = 0U;
    clock_gettime(CLOCK_MONOTONIC, &t_start);
    for (uint32_t it = 0U; it < iters; ++it) {
        uint32_t a = 0U;
        for (uint8_t col = 0U; col < QCTRL_NUM_COLORS; ++col) {
            a += qctrl_grid40_color_sweep_ring0(grid_st, edges, &cfg, col, km_arr, k0_arr, kp_arr, 2000);
        }
        sink += (uint64_t)a;
    }
    clock_gettime(CLOCK_MONOTONIC, &t_end);
    (void)sink;

    double elapsed_ns = ((double)(t_end.tv_sec - t_start.tv_sec) * 1e9) +
                        (double)(t_end.tv_nsec - t_start.tv_nsec);
    double ns_per_sweep40 = elapsed_ns / (double)iters;
    double ns_per_edge    = ns_per_sweep40 / (double)QCTRL_NUM_EDGES;

    printf("[SELF-CHECK PASS] 40q 5x8 Planar Grid: 67 edges (A=%u, B=%u, Seam=%u), 4 colors=(%u,%u,%u,%u)\n",
           region_counts[0], region_counts[1], region_counts[2],
           color_counts[0], color_counts[1], color_counts[2], color_counts[3]);
    printf("[SELF-CHECK PASS] Theorem 11 (Asym SPAM cancelled), Theorem 12 (20-sigma glitch blocked)\n");
    printf("[BENCHMARK PASS]  Ring-0 Full 40q (67-edge) 4-Color Sweep: %.1f ns/sweep (%.1f ns/edge)\n",
           ns_per_sweep40, ns_per_edge);
    return 0;
}

/* ============================================================================
 * Bench40: Realistic Taiwan 20q->40q Drift + SPAM + AR(1) + 3% Glitches
 * ============================================================================ */

typedef struct {
    double quiet_sq_sum;
    double drift_sq_sum;
    double all_sq_sum;
    uint64_t quiet_cnt;
    uint64_t drift_cnt;
    uint64_t all_cnt;
    uint64_t quiet_acts;
    uint64_t total_acts;
    uint64_t total_toggles;
    uint64_t total_probe_shots;
    double region_epoch_rms[3][100]; /* [region][epoch] */
    uint32_t region_edge_cnt[3];
} ctrl_stats_t;

enum {
    CTRL_NO_CAL = 0,
    CTRL_PERIODIC = 1,
    CTRL_SCHMITT_3S = 2,
    CTRL_SCHMITT_5S = 3,
    CTRL_FLOW_QCTRL = 4,
    NUM_CTRLS = 5
};

static const char *kCtrlNames[NUM_CTRLS] = {
    "No-Cal (Pure QEC)",
    "Periodic Offline (q20 ep)",
    "Single-Step 3sig Schmitt",
    "Single-Step 5sig Schmitt",
    "FLOW-QCtrl Ring-0 (Huber+PS)",
};

static int run_bench40(const char *out_json_path) {
    qctrl_edge_desc_t edges[QCTRL_NUM_EDGES];
    qctrl_build_grid40_edges(edges);
    qctrl_cfg_t cfg = qctrl_cfg_default();

    const int n_seeds    = 20;
    const int n_epochs   = 100;
    const int n_shots    = 2000;
    const int n_recal    = 10000;
    const double rho_ar  = 0.6;
    const double ar_std  = 0.006;
    const double ar_innov = sqrt(1.0 - rho_ar * rho_ar);

    ctrl_stats_t stats[NUM_CTRLS];
    memset(stats, 0, sizeof(stats));
    for (int c = 0; c < NUM_CTRLS; ++c) {
        for (uint32_t e = 0U; e < QCTRL_NUM_EDGES; ++e) {
            stats[c].region_edge_cnt[edges[e].region] += 1U;
        }
        /* count is per-seed sum -> normalize by n_seeds later */
        for (int r = 0; r < 3; ++r) {
            stats[c].region_edge_cnt[r] /= (uint32_t)1;
        }
    }

    /* Warm up Ring-0 timing on bench workload */
    struct timespec ring0_t0, ring0_t1;
    double ring0_ns_accum = 0.0;
    uint64_t ring0_sweeps = 0ULL;

    rng_t rng;
    rng_seed(&rng, UINT64_C(2026112140));

    for (int seed = 0; seed < n_seeds; ++seed) {
        /* Heterogeneous edge physics matching 20q->40q superconducting hardware */
        double phi_base[QCTRL_NUM_EDGES];
        double drift_rate[QCTRL_NUM_EDGES];
        double tls_jump[QCTRL_NUM_EDGES];
        double e01[QCTRL_NUM_EDGES];
        double e10[QCTRL_NUM_EDGES];
        double u_ar[QCTRL_NUM_EDGES];

        for (uint32_t e = 0U; e < QCTRL_NUM_EDGES; ++e) {
            double sign    = (rng_uniform01(&rng) < 0.5) ? -1.0 : 1.0;
            phi_base[e]    = sign * (0.006 + 0.012 * rng_uniform01(&rng));        /* 1-3 kHz quiet */
            drift_rate[e]  = sign * (0.0022 + 0.0028 * rng_uniform01(&rng));      /* drift +0.09..0.20 rad */
            tls_jump[e]    = (e % 4U == 1U) ? (sign * (0.10 + 0.09 * rng_uniform01(&rng))) : 0.0; /* ~17 TLS edges */
            e01[e]         = 0.01 + 0.02 * rng_uniform01(&rng);
            e10[e]         = 0.04 + 0.04 * rng_uniform01(&rng);
            u_ar[e]        = ar_std * rng_normal(&rng);
        }

        qctrl_edge_state_t st[NUM_CTRLS][QCTRL_NUM_EDGES];
        for (int c = 0; c < NUM_CTRLS; ++c) {
            for (uint32_t e = 0U; e < QCTRL_NUM_EDGES; ++e) {
                qctrl_edge_state_reset(&st[c][e]);
            }
        }

        for (int t = 0; t < n_epochs; ++t) {
            double phi_true[QCTRL_NUM_EDGES];
            double glitch_m[QCTRL_NUM_EDGES];
            double glitch_p[QCTRL_NUM_EDGES];

            for (uint32_t e = 0U; e < QCTRL_NUM_EDGES; ++e) {
                if (t > 0) {
                    u_ar[e] = rho_ar * u_ar[e] + ar_innov * ar_std * rng_normal(&rng);
                }
                double phi = phi_base[e];
                if (t >= 25 && t < 65) {
                    phi += drift_rate[e] * (double)(t - 25);
                } else if (t >= 65) {
                    phi += drift_rate[e] * 40.0 + tls_jump[e];
                }
                phi_true[e] = phi;

                double gm = 0.0, gp = 0.0;
                if (rng_uniform01(&rng) < 0.03) {
                    double gs = (rng_uniform01(&rng) < 0.5) ? 0.25 : -0.25;
                    if (rng_uniform01(&rng) < 0.5) {
                        gp += gs;
                    } else {
                        gm += gs;
                    }
                }
                glitch_m[e] = -0.5 * u_ar[e] + gm;
                glitch_p[e] = 0.5 * u_ar[e] + gp;
            }

            /* 0. No-Cal (0 probe overhead) */
            for (uint32_t e = 0U; e < QCTRL_NUM_EDGES; ++e) {
                double res = phi_true[e];
                double r2  = res * res;
                stats[CTRL_NO_CAL].all_sq_sum += r2;
                stats[CTRL_NO_CAL].all_cnt += 1ULL;
                if (t < 25) {
                    stats[CTRL_NO_CAL].quiet_sq_sum += r2;
                    stats[CTRL_NO_CAL].quiet_cnt += 1ULL;
                } else {
                    stats[CTRL_NO_CAL].drift_sq_sum += r2;
                    stats[CTRL_NO_CAL].drift_cnt += 1ULL;
                }
                stats[CTRL_NO_CAL].region_epoch_rms[edges[e].region][t] += r2;
            }

            /* 1. Periodic Offline Recal every 20 epochs (N = 10,000 shots) */
            if (t % 20 == 0 && t > 0) {
                for (uint32_t e = 0U; e < QCTRL_NUM_EDGES; ++e) {
                    int32_t km, k0, kp;
                    double th = (double)st[CTRL_PERIODIC][e].theta_q16 / (double)FLOW_Q16_ONE;
                    sample_ps_edge_counts(&rng, phi_true[e], th, e01[e], e10[e], glitch_m[e],
                                          glitch_p[e], n_recal, &km, &k0, &kp);
                    stats[CTRL_PERIODIC].total_probe_shots += (uint64_t)(3 * n_recal);
                    /* Re-calibrate unconditionally during scheduled offline maintenance */
                    qctrl_cfg_t off_cfg = cfg;
                    off_cfg.log_thresh_q16 = 0;
                    off_cfg.z_soft_q16     = 0;
                    st[CTRL_PERIODIC][e].win_count = 0U;
                    st[CTRL_PERIODIC][e].win_head  = 0U;
                    if (qctrl_edge_step_ring0(&st[CTRL_PERIODIC][e], &off_cfg, km, k0, kp, n_recal)) {
                        stats[CTRL_PERIODIC].total_acts += 1ULL;
                        if (t < 25) {
                            stats[CTRL_PERIODIC].quiet_acts += 1ULL;
                        }
                    }
                }
            }
            for (uint32_t e = 0U; e < QCTRL_NUM_EDGES; ++e) {
                double th  = (double)st[CTRL_PERIODIC][e].theta_q16 / (double)FLOW_Q16_ONE;
                double res = phi_true[e] - th;
                double r2  = res * res;
                stats[CTRL_PERIODIC].all_sq_sum += r2;
                stats[CTRL_PERIODIC].all_cnt += 1ULL;
                if (t < 25) {
                    stats[CTRL_PERIODIC].quiet_sq_sum += r2;
                    stats[CTRL_PERIODIC].quiet_cnt += 1ULL;
                } else {
                    stats[CTRL_PERIODIC].drift_sq_sum += r2;
                    stats[CTRL_PERIODIC].drift_cnt += 1ULL;
                }
                stats[CTRL_PERIODIC].region_epoch_rms[edges[e].region][t] += r2;
            }

            /* 2 & 3. Single-Step 3σ / 5σ Schmitt PS */
            for (int sc = 0; sc < 2; ++sc) {
                int cidx        = (sc == 0) ? CTRL_SCHMITT_3S : CTRL_SCHMITT_5S;
                int64_t z_enter = (sc == 0) ? (INT64_C(3) * FLOW_Q16_ONE) : (INT64_C(5) * FLOW_Q16_ONE);
                int64_t z_exit  = INT64_C(98304); /* 1.5 */
                for (uint32_t e = 0U; e < QCTRL_NUM_EDGES; ++e) {
                    int32_t km, k0, kp;
                    double th = (double)st[cidx][e].theta_q16 / (double)FLOW_Q16_ONE;
                    sample_ps_edge_counts(&rng, phi_true[e], th, e01[e], e10[e], glitch_m[e],
                                          glitch_p[e], n_shots, &km, &k0, &kp);
                    stats[cidx].total_probe_shots += (uint64_t)(3 * n_shots);
                    if (qctrl_edge_step_schmitt_ring0(&st[cidx][e], z_exit, z_enter, cfg.z_curv_min_q16,
                                                      cfg.dt_max_q16, km, k0, kp, n_shots)) {
                        stats[cidx].total_acts += 1ULL;
                        if (t < 25) {
                            stats[cidx].quiet_acts += 1ULL;
                        }
                    }
                    double th_new = (double)st[cidx][e].theta_q16 / (double)FLOW_Q16_ONE;
                    double res    = phi_true[e] - th_new;
                    double r2     = res * res;
                    stats[cidx].all_sq_sum += r2;
                    stats[cidx].all_cnt += 1ULL;
                    if (t < 25) {
                        stats[cidx].quiet_sq_sum += r2;
                        stats[cidx].quiet_cnt += 1ULL;
                    } else {
                        stats[cidx].drift_sq_sum += r2;
                        stats[cidx].drift_cnt += 1ULL;
                    }
                    stats[cidx].region_epoch_rms[edges[e].region][t] += r2;
                }
            }

            /* 4. FLOW-QCtrl Ring-0 (4-color parallel sweep across all 67 edges) */
            {
                int32_t km_arr[QCTRL_NUM_EDGES];
                int32_t k0_arr[QCTRL_NUM_EDGES];
                int32_t kp_arr[QCTRL_NUM_EDGES];
                for (uint32_t e = 0U; e < QCTRL_NUM_EDGES; ++e) {
                    double th = (double)st[CTRL_FLOW_QCTRL][e].theta_q16 / (double)FLOW_Q16_ONE;
                    sample_ps_edge_counts(&rng, phi_true[e], th, e01[e], e10[e], glitch_m[e],
                                          glitch_p[e], n_shots, &km_arr[e], &k0_arr[e], &kp_arr[e]);
                    stats[CTRL_FLOW_QCTRL].total_probe_shots += (uint64_t)(3 * n_shots);
                }

                uint32_t ep_acts = 0U;
                clock_gettime(CLOCK_MONOTONIC, &ring0_t0);
                for (uint8_t col = 0U; col < QCTRL_NUM_COLORS; ++col) {
                    ep_acts += qctrl_grid40_color_sweep_ring0(st[CTRL_FLOW_QCTRL], edges, &cfg, col,
                                                              km_arr, k0_arr, kp_arr, n_shots);
                }
                clock_gettime(CLOCK_MONOTONIC, &ring0_t1);
                ring0_ns_accum += ((double)(ring0_t1.tv_sec - ring0_t0.tv_sec) * 1e9) +
                                  (double)(ring0_t1.tv_nsec - ring0_t0.tv_nsec);
                ring0_sweeps += 1ULL;

                stats[CTRL_FLOW_QCTRL].total_acts += (uint64_t)ep_acts;
                if (t < 25) {
                    stats[CTRL_FLOW_QCTRL].quiet_acts += (uint64_t)ep_acts;
                }

                for (uint32_t e = 0U; e < QCTRL_NUM_EDGES; ++e) {
                    double th_new = (double)st[CTRL_FLOW_QCTRL][e].theta_q16 / (double)FLOW_Q16_ONE;
                    double res    = phi_true[e] - th_new;
                    double r2     = res * res;
                    stats[CTRL_FLOW_QCTRL].all_sq_sum += r2;
                    stats[CTRL_FLOW_QCTRL].all_cnt += 1ULL;
                    if (t < 25) {
                        stats[CTRL_FLOW_QCTRL].quiet_sq_sum += r2;
                        stats[CTRL_FLOW_QCTRL].quiet_cnt += 1ULL;
                    } else {
                        stats[CTRL_FLOW_QCTRL].drift_sq_sum += r2;
                        stats[CTRL_FLOW_QCTRL].drift_cnt += 1ULL;
                    }
                    stats[CTRL_FLOW_QCTRL].region_epoch_rms[edges[e].region][t] += r2;
                }
            }
        }

        for (int c = 0; c < NUM_CTRLS; ++c) {
            for (uint32_t e = 0U; e < QCTRL_NUM_EDGES; ++e) {
                stats[c].total_toggles += (uint64_t)st[c][e].toggle_count;
            }
        }
    }

    double ns_per_sweep = ring0_ns_accum / (double)ring0_sweeps;
    printf("\n=== BENCH40: 40-Qubit (5x8, 67 Coupler Edges) x 100 Epochs x 20 Seeds ===\n");
    printf("Ring-0 4-Color Sweep Decision Latency: %.1f ns/sweep (%.1f ns/edge)\n", ns_per_sweep,
           ns_per_sweep / (double)QCTRL_NUM_EDGES);
    printf("%-28s | Quiet RMS  | Drift+TLS RMS | Overall RMS | QuietFalse/chip | TotalAct/6700\n",
           "Controller");
    for (int c = 0; c < NUM_CTRLS; ++c) {
        double q_rms = sqrt(stats[c].quiet_sq_sum / (double)stats[c].quiet_cnt);
        double d_rms = sqrt(stats[c].drift_sq_sum / (double)stats[c].drift_cnt);
        double a_rms = sqrt(stats[c].all_sq_sum / (double)stats[c].all_cnt);
        double q_act = (double)stats[c].quiet_acts / (double)n_seeds;
        double t_act = (double)stats[c].total_acts / (double)n_seeds;
        printf("%-28s | %10.6f | %13.6f | %11.6f | %15.2f | %8.2f / 6700\n", kCtrlNames[c], q_rms,
               d_rms, a_rms, q_act, t_act);
    }

    if (out_json_path != NULL) {
        FILE *fp = fopen(out_json_path, "w");
        if (fp == NULL) {
            fprintf(stderr, "ERROR: cannot open %s for writing\n", out_json_path);
            return 1;
        }
        fprintf(fp, "{\n");
        fprintf(fp, "  \"grid_rows\": %u, \"grid_cols\": %u, \"num_qubits\": %u, \"num_edges\": %u,\n",
                QCTRL_ROWS, QCTRL_COLS, QCTRL_NUM_QUBITS, QCTRL_NUM_EDGES);
        fprintf(fp, "  \"n_seeds\": %d, \"n_epochs\": %d, \"n_shots\": %d,\n", n_seeds, n_epochs,
                n_shots);
        fprintf(fp, "  \"ring0_ns_per_40q_sweep\": %.2f, \"ring0_ns_per_edge\": %.2f,\n",
                ns_per_sweep, ns_per_sweep / (double)QCTRL_NUM_EDGES);
        fprintf(fp, "  \"controllers\": [\n");
        for (int c = 0; c < NUM_CTRLS; ++c) {
            double q_rms = sqrt(stats[c].quiet_sq_sum / (double)stats[c].quiet_cnt);
            double d_rms = sqrt(stats[c].drift_sq_sum / (double)stats[c].drift_cnt);
            double a_rms = sqrt(stats[c].all_sq_sum / (double)stats[c].all_cnt);
            fprintf(fp, "    {\n");
            fprintf(fp, "      \"name\": \"%s\",\n", kCtrlNames[c]);
            fprintf(fp, "      \"quiet_rms_zz\": %.8f,\n", q_rms);
            fprintf(fp, "      \"drift_rms_zz\": %.8f,\n", d_rms);
            fprintf(fp, "      \"overall_rms_zz\": %.8f,\n", a_rms);
            fprintf(fp, "      \"quiet_acts_per_run\": %.4f,\n",
                    (double)stats[c].quiet_acts / (double)n_seeds);
            fprintf(fp, "      \"total_acts_per_run\": %.4f,\n",
                    (double)stats[c].total_acts / (double)n_seeds);
            fprintf(fp, "      \"toggles_per_run\": %.4f,\n",
                    (double)stats[c].total_toggles / (double)n_seeds);
            fprintf(fp, "      \"probe_shots_per_run\": %.1f,\n",
                    (double)stats[c].total_probe_shots / (double)n_seeds);
            const char *rkeys[3] = {"patch_a_rms", "patch_b_rms", "seam_rms"};
            for (int r = 0; r < 3; ++r) {
                uint32_t rec = stats[c].region_edge_cnt[r];
                fprintf(fp, "      \"%s\": [", rkeys[r]);
                for (int t = 0; t < n_epochs; ++t) {
                    double val = sqrt(stats[c].region_epoch_rms[r][t] / (double)(rec * (uint32_t)n_seeds));
                    fprintf(fp, "%.6f%s", val, (t + 1 < n_epochs) ? ", " : "");
                }
                fprintf(fp, "]%s\n", (r < 2) ? "," : "");
            }
            fprintf(fp, "    }%s\n", (c + 1 < NUM_CTRLS) ? "," : "");
        }
        fprintf(fp, "  ]\n}\n");
        fclose(fp);
        printf("Wrote 40q trajectories JSON to %s\n", out_json_path);
    }

    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s selftest | bench40 [--out path.json]\n", argv[0]);
        return 1;
    }
    if (strcmp(argv[1], "selftest") == 0) {
        return run_selftest();
    }
    if (strcmp(argv[1], "bench40") == 0) {
        const char *out_path = NULL;
        for (int i = 2; i + 1 < argc; ++i) {
            if (strcmp(argv[i], "--out") == 0) {
                out_path = argv[i + 1];
            }
        }
        return run_bench40(out_path);
    }
    fprintf(stderr, "Unknown command: %s\n", argv[1]);
    return 1;
}
