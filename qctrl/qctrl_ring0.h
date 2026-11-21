#ifndef QCTRL_RING0_H
#define QCTRL_RING0_H

/* ============================================================================
 * FLOW-QCtrl Ring-0 Autonomous 40-Qubit Calibration Engine (C17 / C++20)
 * ============================================================================
 * Pure integer Q16.16 / 128-bit Guard-Bit kernel built directly on FLOW v8.0.0
 * (`../FLOW/core/flow.h`):
 *   - Zero FPU / SIMD instructions (-mgeneral-regs-only clean)
 *   - Zero udiv / sdiv / __divti3 (uses flow_sym_div128_clamped + flow_isqrt)
 *   - Zero heap allocation (static POD topology & state arrays)
 *   - 4-color edge parallel schedule across 5x8 = 40 qubits (67 coupler edges)
 *   - Enforces FLOW Bundle Theorems 9, 11, 12, 13 (Schmitt hysteresis, affine
 *     SPAM cancellation, Huber glitch immunity, Soft-Threshold de-biasing)
 * ============================================================================ */

#include "flow.h"

#ifdef __cplusplus
extern "C" {
#endif

#define QCTRL_ROWS        5U
#define QCTRL_COLS        8U
#define QCTRL_NUM_QUBITS  (QCTRL_ROWS * QCTRL_COLS) /* 40 qubits */
#define QCTRL_NUM_EDGES   67U                       /* 5*7 horiz + 4*8 vert */
#define QCTRL_NUM_COLORS  4U                        /* 2 horiz + 2 vert colors */
#define QCTRL_NUM_LAMBDAS 4U
#define QCTRL_WIN_CAP     8U

enum qctrl_region_e {
    QCTRL_REGION_PATCH_A = 0, /* cols 0..3 -> 17q Patch A rotated surface code */
    QCTRL_REGION_PATCH_B = 1, /* cols 4..7 -> 17q Patch B rotated surface code */
    QCTRL_REGION_SEAM    = 2  /* horizontal couplers (r,3)-(r,4) for lattice surgery */
};

typedef struct {
    uint8_t u;
    uint8_t v;
    uint8_t color;  /* 0..3 conflict-free parallel coloring */
    uint8_t region; /* enum qctrl_region_e */
} qctrl_edge_desc_t;

typedef struct {
    int64_t theta_q16;                      /* calibrated ZZ compensation angle */
    int64_t log_e_pos_q16[QCTRL_NUM_LAMBDAS]; /* signed log-martingale (+) */
    int64_t log_e_neg_q16[QCTRL_NUM_LAMBDAS]; /* signed log-martingale (-) */
    int64_t win_slope_cnt[QCTRL_WIN_CAP];   /* recent kp - km counts */
    int64_t win_curv_cnt[QCTRL_WIN_CAP];    /* recent kp + km - 2*k0 counts */
    int64_t win_se1_q16[QCTRL_WIN_CAP];     /* Q16.16 SE of kp - km */
    uint32_t win_count;
    uint32_t win_head;
    uint32_t act_count;
    uint32_t toggle_count;
    uint8_t schmitt_active;
    uint8_t reserved[7];
} qctrl_edge_state_t;

typedef struct {
    int64_t huber_bound_q16; /* Q16.16 clip bound B (default 3.0 * 65536 = 196608) */
    int64_t inv_lrv_q16;     /* Q16.16 1 / LRV_eff scale (default 0.74 * 65536 = 48497) */
    int64_t lrv_q16;         /* Q16.16 LRV_eff scale (default 1.35 * 65536 = 88474) */
    int64_t log_thresh_q16;  /* Q16.16 ln(800) = 438083 => E_t >= 100 strict alpha<=0.01 */
    int64_t z_soft_q16;      /* Q16.16 soft-threshold multiplier (default 1.0 * 65536) */
    int64_t z_curv_min_q16;  /* Q16.16 minimum curvature SNR (default 2.0 * 65536) */
    int64_t dt_max_q16;      /* Q16.16 trust-region max step (default 0.8 * 65536 = 52429) */
} qctrl_cfg_t;

FLOW_INL qctrl_cfg_t qctrl_cfg_default(void) {
    qctrl_cfg_t c;
    c.huber_bound_q16 = INT64_C(3) * FLOW_Q16_ONE;
    c.inv_lrv_q16     = INT64_C(48497);  /* ~1 / 1.35 in Q16.16 */
    c.lrv_q16         = INT64_C(88474);  /* ~1.35 in Q16.16 */
    c.log_thresh_q16  = INT64_C(438083); /* ln(800) in Q16.16 */
    c.z_soft_q16      = FLOW_Q16_ONE;    /* 1.0 sigma soft de-biasing */
    c.z_curv_min_q16  = INT64_C(2) * FLOW_Q16_ONE;
    c.dt_max_q16      = INT64_C(52429);  /* 0.8 rad trust-region */
    return c;
}

FLOW_INL void qctrl_edge_state_reset(qctrl_edge_state_t *st) {
    uint32_t i;
    st->theta_q16 = 0;
    for (i = 0U; i < QCTRL_NUM_LAMBDAS; ++i) {
        st->log_e_pos_q16[i] = 0;
        st->log_e_neg_q16[i] = 0;
    }
    for (i = 0U; i < QCTRL_WIN_CAP; ++i) {
        st->win_slope_cnt[i] = 0;
        st->win_curv_cnt[i]  = 0;
        st->win_se1_q16[i]   = 0;
    }
    st->win_count      = 0U;
    st->win_head       = 0U;
    st->act_count      = 0U;
    st->toggle_count   = 0U;
    st->schmitt_active = 0U;
    for (i = 0U; i < 7U; ++i) {
        st->reserved[i] = 0U;
    }
}

/* Build 5x8 = 40-qubit planar grid with 67 coupler edges and conflict-free 4-coloring */
FLOW_INL uint32_t qctrl_build_grid40_edges(qctrl_edge_desc_t edges[QCTRL_NUM_EDGES]) {
    uint32_t idx = 0U;
    uint32_t r;
    uint32_t c;
    for (r = 0U; r < QCTRL_ROWS; ++r) {
        for (c = 0U; c + 1U < QCTRL_COLS; ++c) {
            edges[idx].u     = (uint8_t)(r * QCTRL_COLS + c);
            edges[idx].v     = (uint8_t)(r * QCTRL_COLS + c + 1U);
            edges[idx].color = (uint8_t)(c & 1U); /* colors 0, 1: parallel horiz */
            edges[idx].region =
                (uint8_t)((c == 3U) ? QCTRL_REGION_SEAM
                                    : ((c < 3U) ? QCTRL_REGION_PATCH_A : QCTRL_REGION_PATCH_B));
            ++idx;
        }
    }
    for (r = 0U; r + 1U < QCTRL_ROWS; ++r) {
        for (c = 0U; c < QCTRL_COLS; ++c) {
            edges[idx].u      = (uint8_t)(r * QCTRL_COLS + c);
            edges[idx].v      = (uint8_t)((r + 1U) * QCTRL_COLS + c);
            edges[idx].color  = (uint8_t)(2U + (r & 1U)); /* colors 2, 3: parallel vert */
            edges[idx].region = (uint8_t)((c < 4U) ? QCTRL_REGION_PATCH_A : QCTRL_REGION_PATCH_B);
            ++idx;
        }
    }
    return idx;
}

/* Branchless median for <= 8 int64_t values (zero heap, zero qsort, constant-time pass) */
FLOW_INL int64_t qctrl__median_i64(const int64_t *vals, uint32_t n) {
    int64_t a[QCTRL_WIN_CAP];
    uint32_t i;
    uint32_t j;
    if (n == 0U) {
        return 0;
    }
    if (n > QCTRL_WIN_CAP) {
        n = QCTRL_WIN_CAP;
    }
    for (i = 0U; i < n; ++i) {
        a[i] = vals[i];
    }
    for (i = 0U; i < n; ++i) {
        for (j = i + 1U; j < n; ++j) {
            bool sw   = (a[j] < a[i]);
            int64_t lo = flow_sel_i64(sw, a[j], a[i]);
            int64_t hi = flow_sel_i64(sw, a[i], a[j]);
            a[i]       = lo;
            a[j]       = hi;
        }
    }
    if ((n & 1U) != 0U) {
        return a[n >> 1U];
    }
    return flow_ashr_i64(a[(n >> 1U) - 1U] + a[n >> 1U], 1U);
}

/* Shared integer probe jet calculator from raw readout shot counts (0 FPU, 0 udiv) */
FLOW_INL void qctrl__probe_counts_ring0(int32_t km, int32_t k0, int32_t kp, int32_t n_shots,
                                        int64_t *out_slope_cnt, int64_t *out_curv_cnt,
                                        int64_t *out_se1_q16, int64_t *out_z1_q16,
                                        int64_t *out_snr2_q16) {
    int64_t im    = (int64_t)km;
    int64_t i0    = (int64_t)k0;
    int64_t ip    = (int64_t)kp;
    int64_t n     = (int64_t)n_shots;
    int64_t d1    = ip - im;
    int64_t d2    = ip + im - (2 * i0);
    flow_i128 u_p = (flow_i128)(2 * ip + 1) * (flow_i128)(2 * (n - ip) + 1);
    flow_i128 u_0 = (flow_i128)(2 * i0 + 1) * (flow_i128)(2 * (n - i0) + 1);
    flow_i128 u_m = (flow_i128)(2 * im + 1) * (flow_i128)(2 * (n - im) + 1);
    flow_i128 var1_q32 =
        (flow_i128)flow_sym_div128_clamped((u_p + u_m) << 30U, (flow_i128)(n + 1), INT64_C(1) << 60);
    flow_i128 var2_q32 = (flow_i128)flow_sym_div128_clamped(
        (u_p + (4 * u_0) + u_m) << 30U, (flow_i128)(n + 1), INT64_C(1) << 60);
    int64_t se1_q16 = (int64_t)flow_isqrt_round_u128((flow_u128)flow_sel_i128(var1_q32 > 1, var1_q32, 1));
    int64_t se2_q16 = (int64_t)flow_isqrt_round_u128((flow_u128)flow_sel_i128(var2_q32 > 1, var2_q32, 1));
    se1_q16         = flow_sel_i64(se1_q16 > 0, se1_q16, 1);
    se2_q16         = flow_sel_i64(se2_q16 > 0, se2_q16, 1);

    *out_slope_cnt = d1;
    *out_curv_cnt  = d2;
    *out_se1_q16   = se1_q16;
    *out_z1_q16    = flow_sym_div128_clamped(((flow_i128)(-d1)) << 32U, (flow_i128)se1_q16,
                                             INT64_C(32) * FLOW_Q16_ONE);
    *out_snr2_q16  = flow_sym_div128_clamped(((flow_i128)d2) << 32U, (flow_i128)se2_q16,
                                             INT64_C(64) * FLOW_Q16_ONE);
}

/* Single-edge Ring-0 update: Robust Huber-Clipped LRV e-Process + Median/Soft PS 2-Jet */
FLOW_INL bool qctrl_edge_step_ring0(qctrl_edge_state_t *st, const qctrl_cfg_t *cfg, int32_t km,
                                    int32_t k0, int32_t kp, int32_t n_shots) {
    /* λ ∈ {0.25, 0.50, 1.00, 2.00} in Q16.16 and 0.5*λ^2 in Q16.16 */
    static const int64_t kLamQ16[QCTRL_NUM_LAMBDAS]      = {16384, 32768, 65536, 131072};
    static const int64_t kHalfLam2Q16[QCTRL_NUM_LAMBDAS] = {2048, 8192, 32768, 131072};
    /* 1 / sqrt(m) in Q16.16 for m = 0..8 */
    static const int64_t kInvSqrtMQ16[QCTRL_WIN_CAP + 1U] = {65536, 65536, 46341, 37837, 32768,
                                                             29309, 26755, 24770, 23170};
    int64_t d1 = 0;
    int64_t d2 = 0;
    int64_t se1_q16 = 0;
    int64_t z_raw_q16 = 0;
    int64_t snr2_q16 = 0;
    int64_t z_clip_q16;
    int64_t z_eff_q16;
    int64_t max_log_q16;
    uint32_t slot;
    uint32_t i;

    qctrl__probe_counts_ring0(km, k0, kp, n_shots, &d1, &d2, &se1_q16, &z_raw_q16, &snr2_q16);
    (void)snr2_q16;

    slot                     = st->win_head;
    st->win_slope_cnt[slot]  = d1;
    st->win_curv_cnt[slot]   = d2;
    st->win_se1_q16[slot]    = se1_q16;
    st->win_head             = (slot + 1U >= QCTRL_WIN_CAP) ? 0U : (slot + 1U);
    if (st->win_count < QCTRL_WIN_CAP) {
        st->win_count += 1U;
    }

    z_clip_q16 = flow_clamp_sym(z_raw_q16, cfg->huber_bound_q16);
    z_eff_q16  = (int64_t)(((flow_i128)z_clip_q16 * (flow_i128)cfg->inv_lrv_q16) >> FLOW_Q16_SHIFT);

    max_log_q16 = INT64_MIN;
    for (i = 0U; i < QCTRL_NUM_LAMBDAS; ++i) {
        int64_t lz_q16 = (int64_t)(((flow_i128)kLamQ16[i] * (flow_i128)z_eff_q16) >> FLOW_Q16_SHIFT);
        int64_t lp     = st->log_e_pos_q16[i] + lz_q16 - kHalfLam2Q16[i];
        int64_t ln     = st->log_e_neg_q16[i] - lz_q16 - kHalfLam2Q16[i];
        lp             = flow_sel_i64(lp < -INT64_C(200) * FLOW_Q16_ONE, -INT64_C(200) * FLOW_Q16_ONE, lp);
        ln             = flow_sel_i64(ln < -INT64_C(200) * FLOW_Q16_ONE, -INT64_C(200) * FLOW_Q16_ONE, ln);
        st->log_e_pos_q16[i] = lp;
        st->log_e_neg_q16[i] = ln;
        max_log_q16          = flow_sel_i64(lp > max_log_q16, lp, max_log_q16);
        max_log_q16          = flow_sel_i64(ln > max_log_q16, ln, max_log_q16);
    }

    if (max_log_q16 < cfg->log_thresh_q16) {
        return false;
    }

    {
        uint32_t m          = st->win_count;
        int64_t med_d1      = qctrl__median_i64(st->win_slope_cnt, m);
        int64_t med_d2      = qctrl__median_i64(st->win_curv_cnt, m);
        int64_t med_se_q16  = qctrl__median_i64(st->win_se1_q16, m);
        int64_t inv_sqrt_m  = kInvSqrtMQ16[m];
        flow_i128 se_eff    = (((flow_i128)med_se_q16 * (flow_i128)inv_sqrt_m) >> FLOW_Q16_SHIFT);
        flow_i128 se_med    = (se_eff * (flow_i128)cfg->lrv_q16) >> FLOW_Q16_SHIFT;
        int64_t bias_q16    = (int64_t)((se_med * (flow_i128)cfg->z_soft_q16) >> FLOW_Q16_SHIFT);
        int64_t abs_d1_q16  = (int64_t)(flow_abs_i64(med_d1) << FLOW_Q16_SHIFT);
        int64_t deb_mag_q16 = abs_d1_q16 - bias_q16;
        bool actuated       = false;

        if (deb_mag_q16 > 0 && med_d2 > 0) {
            flow_i128 num_q16 = (med_d1 < 0) ? (flow_i128)deb_mag_q16 : -(flow_i128)deb_mag_q16;
            int64_t step_q16  = flow_sym_div128_clamped(num_q16, (flow_i128)med_d2, cfg->dt_max_q16);
            if (step_q16 != 0) {
                st->theta_q16 += step_q16;
                st->act_count += 1U;
                actuated = true;
            }
        }
        for (i = 0U; i < QCTRL_NUM_LAMBDAS; ++i) {
            st->log_e_pos_q16[i] = 0;
            st->log_e_neg_q16[i] = 0;
        }
        st->win_count = 0U;
        st->win_head  = 0U;
        return actuated;
    }
}

/* Baseline single-step SNR-Schmitt PS controller in identical pure Ring-0 integer math */
FLOW_INL bool qctrl_edge_step_schmitt_ring0(qctrl_edge_state_t *st, int64_t z_exit_q16,
                                            int64_t z_enter_q16, int64_t z_curv_min_q16,
                                            int64_t dt_max_q16, int32_t km, int32_t k0, int32_t kp,
                                            int32_t n_shots) {
    int64_t d1 = 0;
    int64_t d2 = 0;
    int64_t se1_q16 = 0;
    int64_t z_raw_q16 = 0;
    int64_t snr2_q16 = 0;
    int64_t abs_z_q16;
    uint8_t next_active;

    qctrl__probe_counts_ring0(km, k0, kp, n_shots, &d1, &d2, &se1_q16, &z_raw_q16, &snr2_q16);
    (void)se1_q16;

    abs_z_q16   = (int64_t)flow_abs_i64(z_raw_q16);
    next_active = st->schmitt_active ? (uint8_t)(abs_z_q16 >= z_exit_q16)
                                     : (uint8_t)(abs_z_q16 >= z_enter_q16);
    if (next_active != st->schmitt_active) {
        st->toggle_count += 1U;
    }
    st->schmitt_active = next_active;
    if (next_active == 0U || d2 <= 0 || snr2_q16 < z_curv_min_q16) {
        return false;
    }

    {
        /* James-Stein shrinkage: 1 - 1 / snr^2 */
        flow_i128 snr2_sq = ((flow_i128)abs_z_q16 * (flow_i128)abs_z_q16) >> FLOW_Q16_SHIFT; /* Q16 */
        int64_t inv_sq_q16 =
            flow_sym_div128_clamped((flow_i128)FLOW_Q16_ONE << FLOW_Q16_SHIFT, snr2_sq, FLOW_Q16_ONE);
        int64_t shrink_q16 = FLOW_Q16_ONE - inv_sq_q16;
        flow_i128 num_q16  = ((flow_i128)(-d1) * (flow_i128)shrink_q16);
        int64_t step_q16   = flow_sym_div128_clamped(num_q16, (flow_i128)d2, dt_max_q16);
        if (step_q16 != 0) {
            st->theta_q16 += step_q16;
            st->act_count += 1U;
            return true;
        }
    }
    return false;
}

#ifdef QCTRL_EXPORT_RING0
uint32_t qctrl_grid40_color_sweep_ring0(qctrl_edge_state_t states[QCTRL_NUM_EDGES],
                                        const qctrl_edge_desc_t edges[QCTRL_NUM_EDGES],
                                        const qctrl_cfg_t *cfg, uint8_t target_color,
                                        const int32_t km[QCTRL_NUM_EDGES],
                                        const int32_t k0[QCTRL_NUM_EDGES],
                                        const int32_t kp[QCTRL_NUM_EDGES], int32_t n_shots);
bool qctrl_edge_step_ring0_export(qctrl_edge_state_t *st, const qctrl_cfg_t *cfg, int32_t km,
                                  int32_t k0, int32_t kp, int32_t n_shots);
#else
FLOW_INL uint32_t qctrl_grid40_color_sweep_ring0(qctrl_edge_state_t states[QCTRL_NUM_EDGES],
                                                 const qctrl_edge_desc_t edges[QCTRL_NUM_EDGES],
                                                 const qctrl_cfg_t *cfg, uint8_t target_color,
                                                 const int32_t km[QCTRL_NUM_EDGES],
                                                 const int32_t k0[QCTRL_NUM_EDGES],
                                                 const int32_t kp[QCTRL_NUM_EDGES],
                                                 int32_t n_shots) {
    uint32_t acts = 0U;
    uint32_t e;
    for (e = 0U; e < QCTRL_NUM_EDGES; ++e) {
        if (edges[e].color == target_color) {
            if (qctrl_edge_step_ring0(&states[e], cfg, km[e], k0[e], kp[e], n_shots)) {
                acts += 1U;
            }
        }
    }
    return acts;
}
#endif

#ifdef __cplusplus
}
#endif

#endif /* QCTRL_RING0_H */
