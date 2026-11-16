/*
 * qua_invsqrt.h — Pure-Search Bilevel SEARU x FLOW 2-Jet Benchmark I (`invsqrt`) Synthesizer
 *
 * Replaces all hardcoded parameters and heuristics in `build_invsqrt_2jet_v2.py` with
 * Bilevel Meta-SEARU + FLOW 2-Jet search over the 27-bit architectural & lazy-carry genome:
 *   - (K_m, w0, w1, w2, wd2, wsq, F_out, M1, M2, M3)
 *   - Exact Q * (Toffoli + rot_T / 4) cost law from `verify_submission_invsqrt.py`
 *   - Dynamic per-candidate frustration compass (`frust`)
 *   - Direct pure-C17 Kickmix v2 (.kmb) emission via `qua_kmb.h`
 */
#ifndef QUA_INVSQRT_H
#define QUA_INVSQRT_H

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "flow_opt_params.h"
#include "flow.h"
#include "searu.h"
#include "qua_kmb.h"
#include "qua_ele.h"
#include "qua_bitop.h"

#define QUA_INVSQRT_MAX_X      3133452u
#define QUA_INVSQRT_WIDTH      22u
#define QUA_INVSQRT_MAX_CHARTS 1024u
#define QUA_INVSQRT_MAX_SQ     4096u
#define QUA_INVSQRT_STRAT_PTS  3072u

typedef struct {
    uint32_t k_m;
    uint32_t w0;
    uint32_t w1;
    uint32_t w2;
    uint32_t wd1;
    uint32_t wd2;
    uint32_t wsq;
    uint32_t f_out;
    uint32_t m1;
    uint32_t m2;
    uint32_t m3;
    int32_t  sh1;
    int32_t  sh2;
    uint32_t drop;
    uint32_t split_c2;
    uint32_t b_sh_sq;
    uint32_t zero_prods_mode;
    uint32_t sq_arith_mode;
    uint32_t b_sh_main;
    uint32_t phase_split_mode;
    uint32_t cuccaro_mode;

    uint32_t n_charts;
    uint32_t sq_n;
    uint32_t c0_biased[QUA_INVSQRT_MAX_CHARTS];
    uint32_t c1[QUA_INVSQRT_MAX_CHARTS];
    uint32_t c2[QUA_INVSQRT_MAX_CHARTS];
    uint32_t sq_tab[QUA_INVSQRT_MAX_SQ];

    uint32_t w_res;
    uint32_t sq_row_hi[16];
    uint32_t c1_row_hi[24];

    uint32_t tier1_b[22];
    uint32_t tier1_n;
    uint32_t tier2_b[22];
    uint32_t tier2_n;
    uint32_t tier3_b[22];
    uint32_t tier3_n;
    bool     tiers_valid;
    uint32_t tier_overflow_bits;

    double   worst_err_ratio;
    int64_t  min_margin;
    double   rot_c;
    uint32_t peak_qubits;
    uint64_t est_toffolis;
    double   est_score;
} qua_invsqrt_cfg_t;

typedef struct {
    uint64_t target[QUA_INVSQRT_MAX_X + 1u];
    uint64_t tol[QUA_INVSQRT_MAX_X + 1u];
    double   cheb_a0[5][QUA_INVSQRT_MAX_CHARTS];
    double   cheb_a1[5][QUA_INVSQRT_MAX_CHARTS];
    double   cheb_a2[5][QUA_INVSQRT_MAX_CHARTS];
} qua_invsqrt_ref_t;

static inline void qua_invsqrt_cheb_quad(double scale, double u_mid, double h,
                                         double *a0, double *a1, double *a2);

static inline void qua_invsqrt_init_ref(qua_invsqrt_ref_t *ref) {
    /* Exact 130-bit integer reference matching verify_submission_invsqrt.py */
    ref->target[1] = 0ULL;
    ref->tol[1]    = UINT64_C(1) << 40u;
    ref->target[2] = UINT64_C(13043817825332782212);
    ref->tol[2]    = UINT64_C(777472127993);
    ref->target[3] = UINT64_C(10650232656628343401);
    ref->tol[3]    = UINT64_C(634803334273);
    ref->target[4] = UINT64_C(1) << 63u;
    ref->tol[4]    = UINT64_C(1) << 39u;

    for (uint32_t x = 5u; x <= QUA_INVSQRT_MAX_X; ++x) {
        flow_u128 top  = ((flow_u128)1 << 66u);
        flow_u128 q_hi = top / (uint32_t)x;
        flow_u128 rem  = top % (uint32_t)x;
        flow_u128 q128 = (q_hi << 64u) | ((rem << 64u) / (uint32_t)x);
        uint64_t r = (uint64_t)ldexp(1.0 / sqrt((double)x), 65);
        while ((flow_u128)r * (flow_u128)r > q128) --r;
        while ((flow_u128)(r + 1u) * (flow_u128)(r + 1u) <= q128) ++r;
        ref->target[x] = (r + 1ULL) >> 1u;
        ref->tol[x]    = r >> 25u;
    }
    /* Precompute continuous Chebyshev quadratic chart jets once for all K_m in [5..9] */
    for (uint32_t km = 5u; km <= 9u; ++km) {
        uint32_t n_ch = 1u << (km + 1u);
        double h = 1.0 / (double)(1u << km);
        for (uint32_t c = 0u; c < n_ch; ++c) {
            uint32_t r_bit = (c >> km) & 1u;
            uint32_t m_idx = c & ((1u << km) - 1u);
            double scale = (r_bit == 1u) ? 1.0 : 2.0;
            double u_mid = 1.0 + ((double)m_idx + 0.5) * h;
            qua_invsqrt_cheb_quad(scale, u_mid, h,
                                  &ref->cheb_a0[km - 5u][c],
                                  &ref->cheb_a1[km - 5u][c],
                                  &ref->cheb_a2[km - 5u][c]);
        }
    }
}

static inline void qua_invsqrt_decode_genome(uint64_t w, qua_invsqrt_cfg_t *cfg) {
    static const uint32_t KM_MAP[8] = {5u, 6u, 6u, 7u, 7u, 8u, 8u, 9u};
    memset(cfg, 0, sizeof(*cfg));
    cfg->k_m   = KM_MAP[w & 7u];                             /* [5..9] smooth */
    cfg->w0    = 23u + (uint32_t)((w >> 3u) & 7u);           /* [23..30] */
    cfg->w1    = 11u + (uint32_t)((w >> 6u) & 7u);           /* [11..18] */
    cfg->w2    = 4u + (uint32_t)((w >> 9u) & 7u);            /* [4..11] */
    cfg->wd2   = 4u + (uint32_t)(((w >> 12u) & 7u) % 6u);    /* [4..9] */
    cfg->wsq   = 3u + (uint32_t)((w >> 15u) & 7u);           /* [3..10] */
    cfg->f_out = 29u + (uint32_t)((w >> 18u) & 7u);          /* [29..36] */

    uint32_t m1_off = 12u + (uint32_t)((w >> 21u) & 7u);     /* [12..19] */
    uint32_t m2_off = 3u + (uint32_t)((w >> 24u) & 3u);      /* [3..6] */
    uint32_t m3_off = 3u + (uint32_t)((w >> 26u) & 3u);      /* [3..6] */
    cfg->split_c2         = (uint32_t)((w >> 28u) & 1u);
    cfg->b_sh_sq          = (uint32_t)((w >> 29u) & 1u);
    cfg->zero_prods_mode  = (uint32_t)((w >> 30u) & 1u);
    cfg->sq_arith_mode    = (uint32_t)((w >> 31u) & 1u);
    cfg->b_sh_main        = (uint32_t)((w >> 32u) & 3u);     /* [0..3]: B=1,2,4,8 */
    cfg->phase_split_mode = (uint32_t)((w >> 34u) & 1u);
    cfg->cuccaro_mode     = (uint32_t)((w >> 35u) & 1u);
    if (cfg->sq_arith_mode) cfg->b_sh_sq = 0u;
    if (cfg->cuccaro_mode)  cfg->zero_prods_mode = 1u;

    cfg->wd1 = 21u - cfg->k_m;
    if (cfg->wd2 > cfg->wd1) cfg->wd2 = cfg->wd1;

    cfg->m1 = (cfg->w0 > m1_off + 4u) ? (cfg->w0 - m1_off) : 4u;
    cfg->m2 = cfg->m1 + m2_off;
    cfg->m3 = cfg->m2 + m3_off;
    if (cfg->m3 >= cfg->w0) cfg->m3 = cfg->w0 - 1u;
    if (cfg->m2 >= cfg->m3) cfg->m2 = (cfg->m3 > 1u) ? (cfg->m3 - 1u) : 1u;
    if (cfg->m1 >= cfg->m2) cfg->m1 = (cfg->m2 > 1u) ? (cfg->m2 - 1u) : 1u;

    uint32_t e1 = cfg->k_m + 1u;
    uint32_t e2 = 2u * cfg->k_m + 1u;
    cfg->sh1  = (int32_t)cfg->w1 + (int32_t)e1 + (int32_t)cfg->wd1 - (int32_t)cfg->w0;
    cfg->sh2  = (int32_t)cfg->w2 + (int32_t)e2 + (int32_t)cfg->wsq + 2 - (int32_t)cfg->w0;
    if (cfg->f_out > cfg->w0 + 10u) cfg->f_out = cfg->w0 + 10u;
    cfg->drop = (cfg->w0 + 10u) - cfg->f_out;
    cfg->n_charts = 1u << (cfg->k_m + 1u);
    cfg->sq_n     = 1u << (cfg->wd2 - 1u);
}

static inline uint32_t qua_invsqrt_shift_term(uint32_t val, int32_t rsh) {
    if (rsh >= 32) return 0u;
    if (rsh >= 0) return val >> (uint32_t)rsh;
    uint32_t lsh = (uint32_t)(-rsh);
    return (lsh < 32u) ? (val << lsh) : 0u;
}

static inline uint32_t qua_invsqrt_bias_c1(const qua_invsqrt_cfg_t *cfg, uint32_t c1_val) {
    uint32_t total = 0u;
    for (uint32_t b = 0u; b + 1u < cfg->wd1; ++b) {
        total += qua_invsqrt_shift_term(c1_val, cfg->sh1 - (int32_t)b);
    }
    return total;
}

static inline void qua_invsqrt_cheb_quad(double scale, double u_mid, double h,
                                         double *a0, double *a1, double *a2) {
    /* 3-node Chebyshev interpolation on d in [-0.5, 0.5] */
    static const double PI = 3.14159265358979323846;
    double d[3], y[3];
    for (int k = 0; k < 3; ++k) {
        d[k] = 0.5 * cos(PI * (double)(2 * k + 1) / 6.0);
        y[k] = 1.0 / sqrt(scale * (u_mid + h * d[k]));
    }
    /* Lagrange interpolation for quadratic through (d0, y0), (d1, y1), (d2, y2) */
    double c2 = ((y[0] - y[1]) / (d[0] - d[1]) - (y[1] - y[2]) / (d[1] - d[2])) / (d[0] - d[2]);
    double c1 = (y[0] - y[1]) / (d[0] - d[1]) - c2 * (d[0] + d[1]);
    double c0 = y[0] - c1 * d[0] - c2 * d[0] * d[0];
    /* Minimax center offset across [-0.5, 0.5] (9-point Chebyshev check) */
    double min_e = 1e9, max_e = -1e9;
    for (int i = 0; i <= 8; ++i) {
        double di = -0.5 + (double)i / 8.0;
        double err = (c0 + c1 * di + c2 * di * di) - 1.0 / sqrt(scale * (u_mid + h * di));
        if (err < min_e) min_e = err;
        if (err > max_e) max_e = err;
    }
    *a0 = c0 - 0.5 * (min_e + max_e);
    *a1 = -c1;
    *a2 = c2;
}

static inline bool qua_invsqrt_use_sq_arith(const qua_invsqrt_cfg_t *cfg) {
    if (!cfg->sq_arith_mode || cfg->wd2 < 2u) return false;
    uint32_t k_w = cfg->wd2 - 1u;
    int32_t shift = (int32_t)(2u * k_w - 1u) - (int32_t)cfg->wsq;
    return (shift >= 1 && (uint32_t)shift < k_w);
}

static inline uint32_t qua_invsqrt_sq_arith_val(const qua_invsqrt_cfg_t *cfg, uint32_t m) {
    uint32_t k_w = cfg->wd2 - 1u;
    int32_t shift = (int32_t)(2u * k_w - 1u) - (int32_t)cfg->wsq;
    int32_t base_pos = shift - 1;
    uint32_t W = cfg->wsq + 1u;
    uint64_t acc = 1ULL;
    for (uint32_t i = 0u; i < k_w; ++i) {
        if (!((m >> i) & 1u)) continue;
        int32_t p1 = (int32_t)i - 1;
        if (p1 >= base_pos) {
            uint32_t b = (uint32_t)(p1 - base_pos);
            if (b < W) acc += (UINT64_C(1) << b);
        }
        int32_t p2 = 2 * (int32_t)i - 1;
        if (p2 >= base_pos) {
            uint32_t b = (uint32_t)(p2 - base_pos);
            if (b < W) acc += (UINT64_C(1) << b);
        }
    }
    for (uint32_t i = 0u; i < k_w; ++i) {
        if (!((m >> i) & 1u)) continue;
        for (uint32_t j = i + 1u; j < k_w; ++j) {
            if (!((m >> j) & 1u)) continue;
            int32_t pos = (int32_t)(i + j);
            if (pos >= base_pos) {
                uint32_t b = (uint32_t)(pos - base_pos);
                if (b < W) acc += (UINT64_C(1) << b);
            }
        }
    }
    acc &= (UINT64_C(1) << W) - UINT64_C(1);
    return (uint32_t)(acc >> 1u);
}

static inline void qua_invsqrt_calibrate(qua_invsqrt_cfg_t *cfg, const qua_invsqrt_ref_t *ref) {
    uint32_t half_n = cfg->sq_n;
    double denom = (double)(1u << cfg->wd2);
    double sq_scale = 4.0 * (double)(1u << cfg->wsq);
    uint32_t sq_max = (1u << cfg->wsq) - 1u;
    bool use_arith = qua_invsqrt_use_sq_arith(cfg);
    for (uint32_t a = 0u; a < half_n; ++a) {
        if (use_arith) {
            cfg->sq_tab[a] = qua_invsqrt_sq_arith_val(cfg, half_n - 1u - a);
        } else {
            double d = ((double)(half_n - 1u - a) + 0.5) / denom;
            long v = llround(d * d * sq_scale);
            if (v < 0) v = 0;
            if ((uint64_t)v > sq_max) v = (long)sq_max;
            cfg->sq_tab[a] = (uint32_t)v;
        }
    }

    uint32_t e1 = cfg->k_m + 1u;
    uint32_t e2 = 2u * cfg->k_m + 1u;
    uint32_t w0_mask = (1u << cfg->w0) - 1u;
    uint32_t w1_max  = (1u << cfg->w1) - 1u;
    uint32_t w2_max  = (1u << cfg->w2) - 1u;
    uint32_t km_idx  = cfg->k_m - 5u;
    double scale1 = (double)(1ULL << (cfg->w1 + e1));
    double scale2 = (double)(1ULL << (cfg->w2 + e2));
    double scale0 = (double)(1ULL << cfg->w0);

    for (uint32_t c = 0u; c < cfg->n_charts; ++c) {
        double a0 = ref->cheb_a0[km_idx][c];
        double a1 = ref->cheb_a1[km_idx][c];
        double a2 = ref->cheb_a2[km_idx][c];

        long v1 = llround(a1 * scale1);
        long v2 = llround(a2 * scale2);
        if (v1 < 0) v1 = 0;
        if ((uint64_t)v1 > w1_max) v1 = (long)w1_max;
        if (v2 < 0) v2 = 0;
        if ((uint64_t)v2 > w2_max) v2 = (long)w2_max;

        long v0 = llround(a0 * scale0) - (long)qua_invsqrt_bias_c1(cfg, (uint32_t)v1);
        if (v0 < 0) v0 = 0;
        if ((uint64_t)v0 > w0_mask) v0 = (long)w0_mask;

        cfg->c0_biased[c] = (uint32_t)v0;
        cfg->c1[c]        = (uint32_t)v1;
        cfg->c2[c]        = (uint32_t)v2;
    }
}

static inline uint64_t qua_invsqrt_eval_x(const qua_invsqrt_cfg_t *cfg, uint32_t x) {
    uint32_t bit_len = 32u - (uint32_t)__builtin_clz(x);
    uint32_t clz     = QUA_INVSQRT_WIDTH - bit_len;
    uint32_t x_norm  = x << clz;
    uint32_t s       = clz >> 1u;
    uint32_t rem     = clz & 1u;
    uint32_t chart   = (rem << cfg->k_m) | ((x_norm >> (21u - cfg->k_m)) & ((1u << cfg->k_m) - 1u));
    uint32_t t_raw   = x_norm & ((1u << cfg->wd1) - 1u);
    uint32_t t_inv   = ((1u << cfg->wd1) - 1u) ^ t_raw;

    uint32_t t_inv_slice = (t_inv >> (cfg->wd1 - cfg->wd2)) & ((1u << (cfg->wd2 - 1u)) - 1u);
    uint32_t sign_inv    = (t_inv >> (cfg->wd1 - 1u)) & 1u;
    uint32_t sq_addr     = t_inv_slice ^ (sign_inv * ((1u << (cfg->wd2 - 1u)) - 1u));
    uint32_t sq_val      = cfg->sq_tab[sq_addr];

    uint64_t acc = cfg->c0_biased[chart];
    uint32_t c1v = cfg->c1[chart];
    uint32_t c2v = cfg->c2[chart];
    for (uint32_t b = 0u; b < cfg->wd1; ++b) {
        if ((t_inv >> b) & 1u) {
            acc += qua_invsqrt_shift_term(c1v, cfg->sh1 - (int32_t)b);
        }
    }
    for (uint32_t b = 0u; b < cfg->wsq; ++b) {
        if ((sq_val >> b) & 1u) {
            acc += qua_invsqrt_shift_term(c2v, cfg->sh2 - (int32_t)b);
        }
    }
    if (cfg->phase_split_mode) {
        uint64_t y0  = ((uint64_t)cfg->c0_biased[chart] << s) >> cfg->drop;
        uint64_t res = acc - (uint64_t)cfg->c0_biased[chart];
        uint64_t yr  = (res << s) >> cfg->drop;
        return y0 + yr;
    }
    uint64_t acc_mod = acc & ((UINT64_C(1) << cfg->w0) - 1u);
    return (acc_mod << s) >> cfg->drop;
}

static inline void qua_invsqrt_check_tiers_only(qua_invsqrt_cfg_t *cfg);

/* Generate exact staircase-cusp + Chebyshev-extremum inputs for chart `c` across [-0.5, +0.5] */
static inline uint32_t qua_invsqrt_chart_cusp_xs(
    const qua_invsqrt_cfg_t *cfg,
    uint32_t c,
    uint32_t out_xs[20]
) {
    uint32_t rem       = (c >> cfg->k_m) & 1u;
    uint32_t m_idx     = c & ((1u << cfg->k_m) - 1u);
    uint32_t base_norm = (1u << 21u) | (m_idx << cfg->wd1);
    uint32_t t_max     = (1u << cfg->wd1) - 1u;
    uint32_t s_step    = 1u << (cfg->wd1 - cfg->wd2);
    uint32_t n_steps   = 1u << cfg->wd2;
    uint32_t q1_lo     = (n_steps >> 3u) * s_step;
    uint32_t q2_lo     = (n_steps >> 2u) * s_step;
    uint32_t mid_hi    = (n_steps >> 1u) * s_step;

    uint32_t ts[16] = {
        0u,
        (s_step > 1u) ? 1u : 0u,
        s_step - 1u,
        s_step,
        (2u * s_step <= t_max) ? (2u * s_step - 1u) : t_max,
        q1_lo,
        q1_lo + s_step - 1u,
        q2_lo,
        q2_lo + s_step - 1u,
        mid_hi - 1u,
        mid_hi,
        t_max - (q2_lo + s_step - 1u),
        t_max - q2_lo,
        t_max - s_step,
        t_max - (s_step - 1u),
        t_max
    };

    uint32_t n_ok = 0u;
    for (uint32_t i = 0u; i < 16u; ++i) {
        uint32_t raw_norm = base_norm | ts[i];
        uint32_t clz_use  = (rem == 1u) ? 1u : ((raw_norm <= QUA_INVSQRT_MAX_X) ? 0u : 2u);
        uint32_t x        = raw_norm >> clz_use;
        if (x >= 1u && x <= QUA_INVSQRT_MAX_X) {
            uint32_t bl = 32u - (uint32_t)__builtin_clz(x);
            if (QUA_INVSQRT_WIDTH - bl == clz_use) {
                out_xs[n_ok++] = x;
            }
        }
    }
    /* For rem==0 charts that also appear in clz=2, probe clz=2 outer cusps too */
    if (rem == 0u && base_norm <= QUA_INVSQRT_MAX_X) {
        uint32_t x_lo2 = base_norm >> 2u;
        uint32_t x_hi2 = (base_norm | t_max) >> 2u;
        if (x_lo2 >= 1u && n_ok < 20u) out_xs[n_ok++] = x_lo2;
        if (x_hi2 >= 1u && x_hi2 <= QUA_INVSQRT_MAX_X && n_ok < 20u) out_xs[n_ok++] = x_hi2;
    }
    return n_ok;
}

static inline void qua_invsqrt_eval_error(qua_invsqrt_cfg_t *cfg, const qua_invsqrt_ref_t *ref, bool full_grid) {
    double worst_r = 0.0;
    int64_t min_m  = INT64_C(1) << 62u;
    uint32_t sh64  = 64u - cfg->f_out;

    if (full_grid) {
        for (uint32_t x = 1u; x <= QUA_INVSQRT_MAX_X; ++x) {
            uint64_t y1 = qua_invsqrt_eval_x(cfg, x);
            uint64_t th = y1 << sh64;
            int64_t diff = (int64_t)(th - ref->target[x]);
            uint64_t abs_err = (diff < 0) ? (uint64_t)(-diff) : (uint64_t)diff;
            double r = (double)abs_err / (double)ref->tol[x];
            int64_t m = (int64_t)ref->tol[x] - (int64_t)abs_err;
            if (r > worst_r) worst_r = r;
            if (m < min_m) min_m = m;
        }
    } else {
        uint32_t xs[20];
        for (uint32_t c = 0u; c < cfg->n_charts; ++c) {
            uint32_t n_xs = qua_invsqrt_chart_cusp_xs(cfg, c, xs);
            for (uint32_t i = 0u; i < n_xs; ++i) {
                uint32_t x = xs[i];
                uint64_t th = qua_invsqrt_eval_x(cfg, x) << sh64;
                int64_t diff = (int64_t)(th - ref->target[x]);
                uint64_t abs_err = (diff < 0) ? (uint64_t)(-diff) : (uint64_t)diff;
                double r = (double)abs_err / (double)ref->tol[x];
                int64_t m = (int64_t)ref->tol[x] - (int64_t)abs_err;
                if (r > worst_r) worst_r = r;
                if (m < min_m) min_m = m;
            }
            if (c >= 15u && worst_r > 3.5) break;
        }
        if (worst_r <= 3.5) {
            for (uint32_t x = 1u; x <= 128u; ++x) {
                uint64_t th = qua_invsqrt_eval_x(cfg, x) << sh64;
                int64_t diff = (int64_t)(th - ref->target[x]);
                uint64_t abs_err = (diff < 0) ? (uint64_t)(-diff) : (uint64_t)diff;
                double r = (double)abs_err / (double)ref->tol[x];
                int64_t m = (int64_t)ref->tol[x] - (int64_t)abs_err;
                if (r > worst_r) worst_r = r;
                if (m < min_m) min_m = m;
            }
        }
    }
    cfg->worst_err_ratio = worst_r;
    cfg->min_margin      = min_m;
    {
        if (cfg->phase_split_mode && cfg->w_res == 0u) {
            qua_invsqrt_check_tiers_only(cfg);
        }
        uint32_t n_rot = (cfg->f_out > 3u) ? (cfg->f_out - 3u) : 0u;
        double   n_t1  = (cfg->f_out >= 3u) ? 1.0 : 0.0;
        if (cfg->phase_split_mode) {
            uint32_t f_res = (cfg->w_res + 10u > cfg->drop) ? (cfg->w_res + 10u - cfg->drop) : 0u;
            uint32_t rot_lim = (cfg->f_out > 3u) ? (cfg->f_out - 3u) : 0u;
            n_rot += (f_res < rot_lim) ? f_res : rot_lim;
            if (cfg->f_out >= 3u && f_res > cfg->f_out - 3u) n_t1 += 1.0;
        }
        if (n_rot == 0u) n_rot = 1u;
        double eff_m = (min_m > 16 && worst_r < 1.0) ? (double)min_m : 16.0;
        double eps = 3.14159265358979323846 * (eff_m / 18446744073709551616.0) / (double)n_rot;
        double rot_t = (double)n_rot * (0.56 * (-log2(eps)) + 5.3) + n_t1;
        cfg->rot_c = rot_t * 0.25;
    }
}

static inline double qua_invsqrt_step_th_x(uint32_t w0, uint32_t x) {
    uint32_t s_x = (QUA_INVSQRT_WIDTH - (32u - (uint32_t)__builtin_clz(x))) >> 1u;
    return (double)(UINT64_C(1) << (54u - w0 + s_x));
}

/* Fast closed-form per-chart integer minimax centering across exact staircase cusps */
static inline void qua_invsqrt_center_charts_fast(qua_invsqrt_cfg_t *cfg, const qua_invsqrt_ref_t *ref) {
    uint32_t w0_mask = (1u << cfg->w0) - 1u;
    uint32_t sh64 = 64u - cfg->f_out;

    for (uint32_t c = 0u; c < cfg->n_charts; ++c) {
        uint32_t xs[20];
        uint32_t n_ok = qua_invsqrt_chart_cusp_xs(cfg, c, xs);
        if (n_ok == 0u) continue;

        double min_r = 1e9, max_r = -1e9, orig_wr = 0.0;
        for (uint32_t k = 0u; k < n_ok; ++k) {
            uint32_t x = xs[k];
            uint64_t th = qua_invsqrt_eval_x(cfg, x) << sh64;
            int64_t diff = (int64_t)(th - ref->target[x]);
            uint64_t ae = (diff < 0) ? (uint64_t)(-diff) : (uint64_t)diff;
            double r = (double)ae / (double)ref->tol[x];
            if (r > orig_wr) orig_wr = r;
            double se = (double)diff / qua_invsqrt_step_th_x(cfg->w0, x);
            if (se < min_r) min_r = se;
            if (se > max_r) max_r = se;
        }
        uint32_t orig0 = cfg->c0_biased[c];
        uint32_t best0 = orig0;
        double best_wr = orig_wr;
        int64_t dc_mid = -llround(0.5 * (min_r + max_r));
        for (int64_t dc = dc_mid - 1; dc <= dc_mid + 1; ++dc) {
            int64_t cand = (int64_t)orig0 + dc;
            if (cand < 0 || cand > (int64_t)w0_mask) continue;
            cfg->c0_biased[c] = (uint32_t)cand;
            double wr = 0.0;
            for (uint32_t k = 0u; k < n_ok; ++k) {
                uint32_t x = xs[k];
                uint64_t th = qua_invsqrt_eval_x(cfg, x) << sh64;
                int64_t diff = (int64_t)(th - ref->target[x]);
                uint64_t ae = (diff < 0) ? (uint64_t)(-diff) : (uint64_t)diff;
                double r = (double)ae / (double)ref->tol[x];
                if (r > wr) wr = r;
            }
            if (wr < best_wr) {
                best_wr = wr;
                best0 = (uint32_t)cand;
            }
        }
        cfg->c0_biased[c] = best0;
    }
}

/* Inner FLOW 2-Jet + joint quantized (dc2, dc1, dc0) minimax chart relaxation */
static inline void qua_invsqrt_relax_charts_fast_2jet(qua_invsqrt_cfg_t *cfg, const qua_invsqrt_ref_t *ref) {
    uint32_t w0_mask = (1u << cfg->w0) - 1u;
    uint32_t w1_max  = (1u << cfg->w1) - 1u;
    uint32_t w2_max  = (1u << cfg->w2) - 1u;
    uint32_t sh64    = 64u - cfg->f_out;

    for (uint32_t c = 0u; c < cfg->n_charts; ++c) {
        uint32_t xs[20];
        uint32_t n_xs = qua_invsqrt_chart_cusp_xs(cfg, c, xs);
        if (n_xs == 0u) continue;

        double cur_wr = 0.0;
        for (uint32_t i = 0u; i < n_xs; ++i) {
            uint32_t x = xs[i];
            uint64_t th = qua_invsqrt_eval_x(cfg, x) << sh64;
            int64_t diff = (int64_t)(th - ref->target[x]);
            uint64_t ae = (diff < 0) ? (uint64_t)(-diff) : (uint64_t)diff;
            double r = (double)ae / (double)ref->tol[x];
            if (r > cur_wr) cur_wr = r;
        }
        if (cur_wr < 0.72) continue;

        uint32_t orig0 = cfg->c0_biased[c], orig1 = cfg->c1[c], orig2 = cfg->c2[c];
        uint32_t best0 = orig0, best1 = orig1, best2 = orig2;
        double best_wr = cur_wr;
        uint32_t c1_center = orig1;

        /* 1. FLOW 2-Jet step along c1 fiber coordinate (centered finite-difference 2-Jet) */
        if (orig1 >= 1u && orig1 + 1u <= w1_max) {
            double e_vals[3];
            for (int k = -1; k <= 1; ++k) {
                cfg->c1[c] = (uint32_t)((int32_t)orig1 + k);
                int64_t b_adj = (int64_t)orig0 + (int64_t)qua_invsqrt_bias_c1(cfg, orig1) -
                                (int64_t)qua_invsqrt_bias_c1(cfg, cfg->c1[c]);
                if (b_adj < 0) b_adj = 0;
                if (b_adj > (int64_t)w0_mask) b_adj = (int64_t)w0_mask;
                cfg->c0_biased[c] = (uint32_t)b_adj;
                double wr = 0.0;
                for (uint32_t i = 0u; i < n_xs; ++i) {
                    uint32_t x = xs[i];
                    uint64_t th = qua_invsqrt_eval_x(cfg, x) << sh64;
                    int64_t diff = (int64_t)(th - ref->target[x]);
                    uint64_t ae = (diff < 0) ? (uint64_t)(-diff) : (uint64_t)diff;
                    double r = (double)ae / (double)ref->tol[x];
                    if (r > wr) wr = r;
                }
                e_vals[k + 1] = wr;
            }
            int64_t g1_q16 = llround(0.5 * (e_vals[2] - e_vals[0]) * (double)FLOW_Q16_ONE);
            int64_t g2_q16 = llround((e_vals[2] - 2.0 * e_vals[1] + e_vals[0]) * (double)FLOW_Q16_ONE);
            if (g2_q16 > 0 && g1_q16 != 0) {
                uint64_t st[8] = {(uint64_t)g1_q16, (uint64_t)g2_q16, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL};
                uint64_t bd[8] = {
                    0ULL, 0ULL, 0ULL,
                    ((uint64_t)(uint32_t)(int32_t)(-3 * FLOW_Q16_ONE) << 32u) | (uint64_t)(uint32_t)(int32_t)(3 * FLOW_Q16_ONE),
                    ((uint64_t)(uint32_t)(int32_t)(2 * FLOW_Q16_ONE) << 32u) | 1ULL,
                    0ULL, 0ULL, 0ULL
                };
                flow_step(st, st, bd);
                int64_t d1_jet = llround((double)(int64_t)st[3] / (double)FLOW_Q16_ONE);
                int64_t c1_jet = (int64_t)orig1 + d1_jet;
                if (c1_jet >= 0 && c1_jet <= (int64_t)w1_max) c1_center = (uint32_t)c1_jet;
            }
        }

        /* 2. Local quantized (dc2, dc1) lattice around 2-Jet point + exact per-x dc0 centering */
        for (int32_t dc2 = -1; dc2 <= 1; ++dc2) {
            int64_t c2_c = (int64_t)orig2 + dc2;
            if (c2_c < 0 || c2_c > (int64_t)w2_max) continue;
            cfg->c2[c] = (uint32_t)c2_c;
            for (int32_t dc1 = -2; dc1 <= 2; ++dc1) {
                int64_t c1_c = (int64_t)c1_center + dc1;
                if (c1_c < 0 || c1_c > (int64_t)w1_max) continue;
                cfg->c1[c] = (uint32_t)c1_c;
                int64_t c0_base = (int64_t)orig0 + (int64_t)qua_invsqrt_bias_c1(cfg, orig1) -
                                  (int64_t)qua_invsqrt_bias_c1(cfg, (uint32_t)c1_c);
                if (c0_base < 0) c0_base = 0;
                if (c0_base > (int64_t)w0_mask) c0_base = (int64_t)w0_mask;
                cfg->c0_biased[c] = (uint32_t)c0_base;

                double min_se = 1e9, max_se = -1e9;
                for (uint32_t i = 0u; i < n_xs; ++i) {
                    uint32_t x = xs[i];
                    uint64_t th = qua_invsqrt_eval_x(cfg, x) << sh64;
                    double se = (double)(int64_t)(th - ref->target[x]) / qua_invsqrt_step_th_x(cfg->w0, x);
                    if (se < min_se) min_se = se;
                    if (se > max_se) max_se = se;
                }
                int64_t dc0_mid = -llround(0.5 * (min_se + max_se));
                for (int64_t dc0 = dc0_mid - 2; dc0 <= dc0_mid + 2; ++dc0) {
                    int64_t c0_c = c0_base + dc0;
                    if (c0_c < 0 || c0_c > (int64_t)w0_mask) continue;
                    cfg->c0_biased[c] = (uint32_t)c0_c;
                    double wr = 0.0;
                    for (uint32_t i = 0u; i < n_xs; ++i) {
                        uint32_t x = xs[i];
                        uint64_t th = qua_invsqrt_eval_x(cfg, x) << sh64;
                        int64_t diff = (int64_t)(th - ref->target[x]);
                        uint64_t ae = (diff < 0) ? (uint64_t)(-diff) : (uint64_t)diff;
                        double r = (double)ae / (double)ref->tol[x];
                        if (r > wr) wr = r;
                    }
                    if (wr < best_wr) {
                        best_wr = wr;
                        best0 = (uint32_t)c0_c;
                        best1 = (uint32_t)c1_c;
                        best2 = (uint32_t)c2_c;
                    }
                }
            }
        }
        cfg->c0_biased[c] = best0;
        cfg->c1[c]        = best1;
        cfg->c2[c]        = best2;
    }
}

/* Bucketed O(MAX_X) full-domain joint (dc2, dc1, dc0) minimax + margin maximizer (strictly monotone!) */
static inline void qua_invsqrt_polish_charts(qua_invsqrt_cfg_t *cfg, const qua_invsqrt_ref_t *ref) {
    uint32_t w0_mask = (1u << cfg->w0) - 1u;
    uint32_t w1_max  = (1u << cfg->w1) - 1u;
    uint32_t w2_max  = (1u << cfg->w2) - 1u;
    uint32_t sh64    = 64u - cfg->f_out;
    uint32_t n_ch    = cfg->n_charts;

    uint32_t *counts  = (uint32_t *)calloc((size_t)(n_ch + 1u), sizeof(uint32_t));
    uint32_t *offsets = (uint32_t *)calloc((size_t)(n_ch + 1u), sizeof(uint32_t));
    uint32_t *xs_all  = (uint32_t *)malloc((size_t)QUA_INVSQRT_MAX_X * sizeof(uint32_t));
    if (!counts || !offsets || !xs_all) {
        free(counts); free(offsets); free(xs_all);
        return;
    }

    uint32_t km_mask = (1u << cfg->k_m) - 1u;
    for (uint32_t x = 1u; x <= QUA_INVSQRT_MAX_X; ++x) {
        uint32_t clz = QUA_INVSQRT_WIDTH - (32u - (uint32_t)__builtin_clz(x));
        uint32_t ch = ((clz & 1u) << cfg->k_m) | (((x << clz) >> cfg->wd1) & km_mask);
        counts[ch]++;
    }
    for (uint32_t c = 0u; c < n_ch; ++c) offsets[c + 1u] = offsets[c] + counts[c];
    memset(counts, 0, (size_t)n_ch * sizeof(uint32_t));
    for (uint32_t x = 1u; x <= QUA_INVSQRT_MAX_X; ++x) {
        uint32_t clz = QUA_INVSQRT_WIDTH - (32u - (uint32_t)__builtin_clz(x));
        uint32_t ch = ((clz & 1u) << cfg->k_m) | (((x << clz) >> cfg->wd1) & km_mask);
        xs_all[offsets[ch] + counts[ch]++] = x;
    }

    qua_invsqrt_check_tiers_only(cfg);
    bool orig_tiers_valid = cfg->tiers_valid;

    for (uint32_t c = 0u; c < n_ch; ++c) {
        uint32_t beg = offsets[c], end = offsets[c + 1u];
        if (beg >= end) continue;

        uint32_t orig0 = cfg->c0_biased[c], orig1 = cfg->c1[c], orig2 = cfg->c2[c];
        uint32_t best0 = orig0, best1 = orig1, best2 = orig2;
        int64_t  best_min_m = INT64_C(1) << 62u;
        double   best_wr    = 0.0;

        /* Seed baseline on exact full-bucket (orig0, orig1, orig2) to guarantee strict monotonicity! */
        for (uint32_t idx = beg; idx < end; ++idx) {
            uint32_t x = xs_all[idx];
            uint64_t th = qua_invsqrt_eval_x(cfg, x) << sh64;
            int64_t diff = (int64_t)(th - ref->target[x]);
            uint64_t ae = (diff < 0) ? (uint64_t)(-diff) : (uint64_t)diff;
            int64_t m = (int64_t)ref->tol[x] - (int64_t)ae;
            double r = (double)ae / (double)ref->tol[x];
            if (m < best_min_m) best_min_m = m;
            if (r > best_wr) best_wr = r;
        }

        for (int32_t dc2 = -3; dc2 <= 3; ++dc2) {
            int64_t c2_c = (int64_t)orig2 + dc2;
            if (c2_c < 0 || c2_c > (int64_t)w2_max) continue;
            cfg->c2[c] = (uint32_t)c2_c;
            for (int32_t dc1 = -3; dc1 <= 3; ++dc1) {
                int64_t c1_c = (int64_t)orig1 + dc1;
                if (c1_c < 0 || c1_c > (int64_t)w1_max) continue;
                cfg->c1[c] = (uint32_t)c1_c;
                int64_t c0_base = (int64_t)orig0 + (int64_t)qua_invsqrt_bias_c1(cfg, orig1) -
                                  (int64_t)qua_invsqrt_bias_c1(cfg, (uint32_t)c1_c);
                if (c0_base < 0) c0_base = 0;
                if (c0_base > (int64_t)w0_mask) c0_base = (int64_t)w0_mask;
                cfg->c0_biased[c] = (uint32_t)c0_base;

                /* Exact per-x step_th(x) center estimate across 16 stride points + endpoints */
                double min_se = 1e9, max_se = -1e9;
                uint32_t stride = (end - beg > 16u) ? ((end - beg) / 16u) : 1u;
                for (uint32_t idx = beg; idx < end; idx += stride) {
                    uint32_t x = xs_all[idx];
                    uint64_t th = qua_invsqrt_eval_x(cfg, x) << sh64;
                    double se = (double)(int64_t)(th - ref->target[x]) / qua_invsqrt_step_th_x(cfg->w0, x);
                    if (se < min_se) min_se = se;
                    if (se > max_se) max_se = se;
                }
                {
                    uint32_t x_last = xs_all[end - 1u];
                    uint64_t th_last = qua_invsqrt_eval_x(cfg, x_last) << sh64;
                    double se_last = (double)(int64_t)(th_last - ref->target[x_last]) / qua_invsqrt_step_th_x(cfg->w0, x_last);
                    if (se_last < min_se) min_se = se_last;
                    if (se_last > max_se) max_se = se_last;
                }
                int64_t dc0_mid = -llround(0.5 * (min_se + max_se));

                for (int64_t dc0 = dc0_mid - 3; dc0 <= dc0_mid + 3; ++dc0) {
                    int64_t c0_c = c0_base + dc0;
                    if (c0_c < 0 || c0_c > (int64_t)w0_mask) continue;
                    cfg->c0_biased[c] = (uint32_t)c0_c;

                    int64_t min_m = INT64_C(1) << 62u;
                    double wr = 0.0;
                    for (uint32_t idx = beg; idx < end; ++idx) {
                        uint32_t x = xs_all[idx];
                        uint64_t th = qua_invsqrt_eval_x(cfg, x) << sh64;
                        int64_t diff = (int64_t)(th - ref->target[x]);
                        uint64_t ae = (diff < 0) ? (uint64_t)(-diff) : (uint64_t)diff;
                        int64_t m = (int64_t)ref->tol[x] - (int64_t)ae;
                        double r = (double)ae / (double)ref->tol[x];
                        if (m < min_m) min_m = m;
                        if (r > wr) wr = r;
                        if (best_wr < 1.0 && (wr >= 1.0 || min_m < best_min_m)) break;
                    }
                    bool better = (best_wr >= 1.0 || wr >= 1.0)
                        ? (wr < best_wr)
                        : (min_m > best_min_m || (min_m == best_min_m && wr < best_wr));
                    if (better) {
                        best_wr    = wr;
                        best_min_m = min_m;
                        best0      = (uint32_t)c0_c;
                        best1      = (uint32_t)c1_c;
                        best2      = (uint32_t)c2_c;
                    }
                }
            }
        }
        cfg->c0_biased[c] = best0;
        cfg->c1[c]        = best1;
        cfg->c2[c]        = best2;

        /* Verify lazy-carry tiers still hold if c1/c2 increased; fallback to dc0-only if needed */
        if (best1 > orig1 || best2 > orig2) {
            qua_invsqrt_check_tiers_only(cfg);
            if (!cfg->tiers_valid) {
                cfg->c1[c] = orig1;
                cfg->c2[c] = orig2;
                int64_t best_m0 = INT64_MIN;
                uint32_t b0_only = orig0;
                for (int32_t dc0 = -16; dc0 <= 16; ++dc0) {
                    int64_t c0_c = (int64_t)orig0 + dc0;
                    if (c0_c < 0 || c0_c > (int64_t)w0_mask) continue;
                    cfg->c0_biased[c] = (uint32_t)c0_c;
                    int64_t min_m = INT64_C(1) << 62u;
                    for (uint32_t idx = beg; idx < end; ++idx) {
                        uint32_t x = xs_all[idx];
                        uint64_t th = qua_invsqrt_eval_x(cfg, x) << sh64;
                        int64_t diff = (int64_t)(th - ref->target[x]);
                        uint64_t ae = (diff < 0) ? (uint64_t)(-diff) : (uint64_t)diff;
                        int64_t m = (int64_t)ref->tol[x] - (int64_t)ae;
                        if (m < min_m) min_m = m;
                    }
                    if (min_m > best_m0) { best_m0 = min_m; b0_only = (uint32_t)c0_c; }
                }
                cfg->c0_biased[c] = b0_only;
                qua_invsqrt_check_tiers_only(cfg);
            }
        }
    }

    /* Targeted second pass over bottleneck charts (wr > 0.85) to widen rotation_eps margin */
    for (uint32_t c = 0u; c < n_ch; ++c) {
        uint32_t beg = offsets[c], end = offsets[c + 1u];
        if (beg >= end) continue;

        uint32_t orig0 = cfg->c0_biased[c], orig1 = cfg->c1[c], orig2 = cfg->c2[c];
        uint32_t best0 = orig0, best1 = orig1, best2 = orig2;
        int64_t  best_min_m = INT64_C(1) << 62u;
        double   best_wr    = 0.0;
        for (uint32_t idx = beg; idx < end; ++idx) {
            uint32_t x = xs_all[idx];
            uint64_t th = qua_invsqrt_eval_x(cfg, x) << sh64;
            int64_t diff = (int64_t)(th - ref->target[x]);
            uint64_t ae = (diff < 0) ? (uint64_t)(-diff) : (uint64_t)diff;
            int64_t m = (int64_t)ref->tol[x] - (int64_t)ae;
            double r = (double)ae / (double)ref->tol[x];
            if (m < best_min_m) best_min_m = m;
            if (r > best_wr) best_wr = r;
        }
        if (best_wr <= 0.85) continue;

        for (int32_t dc2 = -8; dc2 <= 8; ++dc2) {
            int64_t c2_c = (int64_t)orig2 + dc2;
            if (c2_c < 0 || c2_c > (int64_t)w2_max) continue;
            cfg->c2[c] = (uint32_t)c2_c;
            for (int32_t dc1 = -10; dc1 <= 10; ++dc1) {
                int64_t c1_c = (int64_t)orig1 + dc1;
                if (c1_c < 0 || c1_c > (int64_t)w1_max) continue;
                cfg->c1[c] = (uint32_t)c1_c;
                int64_t c0_base = (int64_t)orig0 + (int64_t)qua_invsqrt_bias_c1(cfg, orig1) -
                                  (int64_t)qua_invsqrt_bias_c1(cfg, (uint32_t)c1_c);
                if (c0_base < 0) c0_base = 0;
                if (c0_base > (int64_t)w0_mask) c0_base = (int64_t)w0_mask;
                cfg->c0_biased[c] = (uint32_t)c0_base;

                double min_se = 1e9, max_se = -1e9;
                uint32_t stride = (end - beg > 16u) ? ((end - beg) / 16u) : 1u;
                for (uint32_t idx = beg; idx < end; idx += stride) {
                    uint32_t x = xs_all[idx];
                    uint64_t th = qua_invsqrt_eval_x(cfg, x) << sh64;
                    double se = (double)(int64_t)(th - ref->target[x]) / qua_invsqrt_step_th_x(cfg->w0, x);
                    if (se < min_se) min_se = se;
                    if (se > max_se) max_se = se;
                }
                {
                    uint32_t x_last = xs_all[end - 1u];
                    uint64_t th_last = qua_invsqrt_eval_x(cfg, x_last) << sh64;
                    double se_last = (double)(int64_t)(th_last - ref->target[x_last]) / qua_invsqrt_step_th_x(cfg->w0, x_last);
                    if (se_last < min_se) min_se = se_last;
                    if (se_last > max_se) max_se = se_last;
                }
                int64_t dc0_mid = -llround(0.5 * (min_se + max_se));

                for (int64_t dc0 = dc0_mid - 8; dc0 <= dc0_mid + 8; ++dc0) {
                    int64_t c0_c = c0_base + dc0;
                    if (c0_c < 0 || c0_c > (int64_t)w0_mask) continue;
                    cfg->c0_biased[c] = (uint32_t)c0_c;

                    int64_t min_m = INT64_C(1) << 62u;
                    double wr = 0.0;
                    for (uint32_t idx = beg; idx < end; ++idx) {
                        uint32_t x = xs_all[idx];
                        uint64_t th = qua_invsqrt_eval_x(cfg, x) << sh64;
                        int64_t diff = (int64_t)(th - ref->target[x]);
                        uint64_t ae = (diff < 0) ? (uint64_t)(-diff) : (uint64_t)diff;
                        int64_t m = (int64_t)ref->tol[x] - (int64_t)ae;
                        double r = (double)ae / (double)ref->tol[x];
                        if (m < min_m) min_m = m;
                        if (r > wr) wr = r;
                        if (best_wr < 1.0 && (wr >= 1.0 || min_m < best_min_m)) break;
                    }
                    bool better = (best_wr >= 1.0 || wr >= 1.0)
                        ? (wr < best_wr)
                        : (min_m > best_min_m || (min_m == best_min_m && wr < best_wr));
                    if (better) {
                        bool tier_ok = true;
                        if (orig_tiers_valid && ((uint32_t)c1_c > orig1 || (uint32_t)c2_c > orig2)) {
                            qua_invsqrt_check_tiers_only(cfg);
                            tier_ok = cfg->tiers_valid;
                        }
                        if (tier_ok) {
                            best_wr    = wr;
                            best_min_m = min_m;
                            best0      = (uint32_t)c0_c;
                            best1      = (uint32_t)c1_c;
                            best2      = (uint32_t)c2_c;
                        }
                    }
                }
            }
        }
        cfg->c0_biased[c] = best0;
        cfg->c1[c]        = best1;
        cfg->c2[c]        = best2;
        qua_invsqrt_check_tiers_only(cfg);
    }

    free(xs_all);
    free(offsets);
    free(counts);
}

typedef enum {
    QUA_INV_ADD_PRODS_SHORT = 0,
    QUA_INV_ADD_ZERO_PRODS  = 1,
    QUA_INV_ADD_CUCCARO     = 2
} qua_inv_add_kind_t;

static inline qua_inv_add_kind_t qua_inv_select_ctrl_add(
    uint32_t n_src,
    uint32_t n_tgt,
    uint32_t avail_anc,
    bool allow_prods_short,
    bool allow_cuccaro,
    uint32_t *out_anc,
    uint64_t *out_ccx
) {
    if (n_src == 0u || n_tgt == 0u) {
        *out_anc = 0u;
        *out_ccx = 0u;
        return QUA_INV_ADD_ZERO_PRODS;
    }
    uint32_t k = (n_src < n_tgt) ? n_src : n_tgt;
    if (allow_prods_short && k < n_tgt) {
        uint32_t need_ps = (n_tgt == 1u) ? 1u : (k + n_tgt - 1u);
        if (need_ps <= avail_anc) {
            *out_anc = need_ps;
            *out_ccx = (n_tgt == 1u) ? 1ull : (uint64_t)(k + n_tgt - 1u);
            return QUA_INV_ADD_PRODS_SHORT;
        }
    }
    if (!allow_cuccaro || n_tgt <= avail_anc) {
        *out_anc = n_tgt;
        if (k == 1u)         *out_ccx = (uint64_t)n_tgt;
        else if (k == n_tgt) *out_ccx = 2ull * (uint64_t)n_tgt - 1ull;
        else                 *out_ccx = (uint64_t)(k + n_tgt);
        return QUA_INV_ADD_ZERO_PRODS;
    }
    qua_bitop_ctrl_iadd_cuccaro_resources(k, n_tgt, out_anc, out_ccx);
    return QUA_INV_ADD_CUCCARO;
}

static inline uint32_t qua_inv_min_ctrl_add_anc(
    uint32_t n_src,
    uint32_t n_tgt,
    bool zero_prods_mode,
    bool cuccaro_mode
) {
    if (n_src == 0u || n_tgt == 0u) return 0u;
    uint32_t k = (n_src < n_tgt) ? n_src : n_tgt;
    if (cuccaro_mode) {
        uint32_t anc = 0u;
        uint64_t ccx = 0u;
        qua_bitop_ctrl_iadd_cuccaro_resources(k, n_tgt, &anc, &ccx);
        return anc;
    }
    if (!zero_prods_mode && k < n_tgt && n_tgt > 1u) {
        return k + n_tgt - 1u;
    }
    return n_tgt;
}

static inline uint64_t qua_invsqrt_sq_arith_once_ccx(const qua_invsqrt_cfg_t *cfg) {
    uint32_t k_w = cfg->wd2 - 1u;
    int32_t shift = (int32_t)(2u * k_w - 1u) - (int32_t)cfg->wsq;
    int32_t base_pos = shift - 1;
    uint32_t W = cfg->wsq + 1u;
    uint64_t sq_once = 0u;
    for (uint32_t i = 0u; i < k_w; ++i) {
        int32_t p1 = (int32_t)i - 1;
        if (p1 >= base_pos && (uint32_t)(p1 - base_pos) < W) {
            sq_once += (W - (uint32_t)(p1 - base_pos) - 1u);
        }
        int32_t p2 = 2 * (int32_t)i - 1;
        if (p2 >= base_pos && (uint32_t)(p2 - base_pos) < W) {
            sq_once += (W - (uint32_t)(p2 - base_pos) - 1u);
        }
    }
    for (uint32_t i = 0u; i < k_w; ++i) {
        for (uint32_t j = i + 1u; j < k_w; ++j) {
            int32_t pos = (int32_t)(i + j);
            if (pos >= base_pos) {
                uint32_t b = (uint32_t)(pos - base_pos);
                if (b == W - 1u) sq_once += 1u;
                else if (b < W - 1u) sq_once += (W - b);
            }
        }
    }
    return sq_once;
}

static inline uint32_t qua_invsqrt_pick_sub_bsh(
    const uint64_t *tab,
    uint32_t n_rows,
    uint32_t n_addr,
    uint32_t w_bits,
    uint32_t free_anc_cap,
    qua_bitop_qroam_res_t *out_res
) {
    uint32_t best_bsh = 0u;
    qua_bitop_qroam_vslot0_resources(tab, n_rows, 1u, n_addr, w_bits, 0u, 0u, free_anc_cap, out_res);
    for (uint32_t bsh = 1u; bsh <= 2u && bsh < n_addr; ++bsh) {
        qua_bitop_qroam_res_t cand;
        qua_bitop_qroam_vslot0_resources(tab, n_rows, 1u, n_addr, w_bits, bsh, 0u, free_anc_cap, &cand);
        if (cand.anc_needed <= free_anc_cap && cand.total_ccx < out_res->total_ccx) {
            *out_res = cand;
            best_bsh = bsh;
        }
    }
    return best_bsh;
}

static inline void qua_invsqrt_check_tiers_only(qua_invsqrt_cfg_t *cfg) {
    uint32_t max_sq_val = 0u;
    uint32_t sq_or_mask = 0u;
    for (uint32_t a = 0u; a < cfg->sq_n; ++a) {
        if (cfg->sq_tab[a] > max_sq_val) max_sq_val = cfg->sq_tab[a];
        sq_or_mask |= cfg->sq_tab[a];
    }

    if (cfg->phase_split_mode) {
        uint64_t cum[QUA_INVSQRT_MAX_CHARTS];
        memset(cum, 0, (size_t)cfg->n_charts * sizeof(uint64_t));
        memset(cfg->sq_row_hi, 0, sizeof(cfg->sq_row_hi));
        memset(cfg->c1_row_hi, 0, sizeof(cfg->c1_row_hi));
        cfg->tier1_n = cfg->tier2_n = cfg->tier3_n = 0u;

        uint64_t run_max = 0u;
        for (uint32_t b = 0u; b < cfg->wsq; ++b) {
            if (((sq_or_mask >> b) & 1u) == 0u) continue;
            int32_t rsh = cfg->sh2 - (int32_t)b;
            if (rsh >= (int32_t)cfg->w2) continue;
            uint64_t term_max = 0u;
            for (uint32_t c = 0u; c < cfg->n_charts; ++c) {
                uint64_t t = qua_invsqrt_shift_term(cfg->c2[c], rsh);
                if (t > term_max) term_max = t;
            }
            if (term_max == 0u) continue;
            for (uint32_t c = 0u; c < cfg->n_charts; ++c) {
                cum[c] += qua_invsqrt_shift_term(cfg->c2[c], rsh);
                if (cum[c] > run_max) run_max = cum[c];
            }
            cfg->sq_row_hi[b] = 64u - (uint32_t)__builtin_clzll(run_max);
        }
        for (uint32_t b = 0u; b < cfg->wd1; ++b) {
            int32_t rsh = cfg->sh1 - (int32_t)b;
            if (rsh >= (int32_t)cfg->w1) continue;
            uint64_t term_max = 0u;
            for (uint32_t c = 0u; c < cfg->n_charts; ++c) {
                uint64_t t = qua_invsqrt_shift_term(cfg->c1[c], rsh);
                if (t > term_max) term_max = t;
            }
            if (term_max == 0u) continue;
            for (uint32_t c = 0u; c < cfg->n_charts; ++c) {
                cum[c] += qua_invsqrt_shift_term(cfg->c1[c], rsh);
                if (cum[c] > run_max) run_max = cum[c];
            }
            cfg->c1_row_hi[b] = 64u - (uint32_t)__builtin_clzll(run_max);
        }
        cfg->w_res = (run_max > 0u) ? (64u - (uint32_t)__builtin_clzll(run_max)) : 1u;
        cfg->tiers_valid = (cfg->w_res <= cfg->w0);
        cfg->tier_overflow_bits = (cfg->w_res > cfg->w0) ? (cfg->w_res - cfg->w0) : 0u;
        return;
    }

    uint64_t cum[QUA_INVSQRT_MAX_CHARTS];
    uint64_t cum_max = 0u;
    for (uint32_t c = 0u; c < cfg->n_charts; ++c) {
        uint64_t s_max = 0u;
        uint32_t c2v = cfg->c2[c];
        for (uint32_t b = 0u; b < cfg->wsq; ++b) {
            if ((max_sq_val >> b) >= 1u) {
                s_max += qua_invsqrt_shift_term(c2v, cfg->sh2 - (int32_t)b);
            }
        }
        cum[c] = s_max;
        if (s_max > cum_max) cum_max = s_max;
    }

    cfg->w_res = cfg->w0;
    cfg->tiers_valid = false;
    cfg->tier_overflow_bits = 0u;
    cfg->tier1_n = cfg->tier2_n = cfg->tier3_n = 0u;
    bool ok = (cum_max < (UINT64_C(1) << cfg->m1));
    if (!ok) cfg->tier_overflow_bits += 8u;

    for (uint32_t b = 0u; b < cfg->wd1; ++b) {
        uint64_t term_max = 0u;
        uint64_t next_max = 0u;
        for (uint32_t c = 0u; c < cfg->n_charts; ++c) {
            uint64_t t = qua_invsqrt_shift_term(cfg->c1[c], cfg->sh1 - (int32_t)b);
            if (t > term_max) term_max = t;
            if (cum[c] + t > next_max) next_max = cum[c] + t;
        }
        if (term_max == 0u) continue;
        if (ok && cfg->tier2_n == 0u && cfg->tier3_n == 0u && next_max < (UINT64_C(1) << cfg->m1)) {
            cfg->tier1_b[cfg->tier1_n++] = b;
        } else if (ok && cfg->tier3_n == 0u && next_max < (UINT64_C(1) << cfg->m2)) {
            cfg->tier2_b[cfg->tier2_n++] = b;
        } else if (ok && next_max < (UINT64_C(1) << cfg->m3)) {
            cfg->tier3_b[cfg->tier3_n++] = b;
        } else {
            ok = false;
            ++cfg->tier_overflow_bits;
        }
        for (uint32_t c = 0u; c < cfg->n_charts; ++c) {
            cum[c] += qua_invsqrt_shift_term(cfg->c1[c], cfg->sh1 - (int32_t)b);
        }
    }
    cfg->tiers_valid = ok;
}

static inline void qua_invsqrt_plan_tiers_and_cost(qua_invsqrt_cfg_t *cfg) {
    qua_invsqrt_check_tiers_only(cfg);

    uint32_t n_addr    = cfg->k_m + 1u;
    uint32_t b_sh_main = (cfg->b_sh_main < n_addr) ? cfg->b_sh_main : 0u;
    bool     use_arith = qua_invsqrt_use_sq_arith(cfg);
    uint32_t sq_addr_n = cfg->wd2 - 1u;
    uint32_t b_sh_sq   = (!use_arith && cfg->b_sh_sq < sq_addr_n) ? cfg->b_sh_sq : 0u;

    uint64_t sq_tab64[QUA_INVSQRT_MAX_SQ];
    for (uint32_t a = 0u; a < cfg->sq_n; ++a) sq_tab64[a] = cfg->sq_tab[a];

    if (cfg->phase_split_mode) {
        uint64_t c0_tab[QUA_INVSQRT_MAX_CHARTS];
        uint64_t c1_tab[QUA_INVSQRT_MAX_CHARTS];
        uint64_t c2_tab[QUA_INVSQRT_MAX_CHARTS];
        uint64_t c12_tab[QUA_INVSQRT_MAX_CHARTS];
        for (uint32_t c = 0u; c < cfg->n_charts; ++c) {
            c0_tab[c]  = (uint64_t)cfg->c0_biased[c];
            c1_tab[c]  = (uint64_t)cfg->c1[c];
            c2_tab[c]  = (uint64_t)cfg->c2[c];
            c12_tab[c] = (uint64_t)cfg->c1[c] | ((uint64_t)cfg->c2[c] << cfg->w1);
        }

        /* Pass 1: Exact minimum pool watermark above the 27 persistent wires (x[22] + clz[5]) */
        uint32_t peak_pool = 14u; /* log_clz */
        if (peak_pool < cfg->w0 + 10u) peak_pool = cfg->w0 + 10u; /* Stage A barrel shift */

        qua_bitop_qroam_res_t c0_qroam;
        qua_bitop_qroam_vslot0_resources(c0_tab, cfg->n_charts, 1u, n_addr, cfg->w0, 0u, 0u, 0u, &c0_qroam);
        if (peak_pool < cfg->w0 + c0_qroam.anc_needed) peak_pool = cfg->w0 + c0_qroam.anc_needed;

        qua_bitop_hmr_lut_res_t c0_del;
        qua_bitop_hmr_parity_lut_resources(cfg->n_charts, n_addr, cfg->w0,
                                           (peak_pool > cfg->w0) ? (peak_pool - cfg->w0) : 0u, &c0_del);
        if (peak_pool < cfg->w0 + c0_del.pool_anc_needed) peak_pool = cfg->w0 + c0_del.pool_anc_needed;

        uint32_t w_res = cfg->w_res;
        uint32_t live_c2 = w_res + (cfg->split_c2 ? cfg->w2 : (cfg->w1 + cfg->w2));
        uint32_t live_c1 = w_res + (cfg->split_c2 ? cfg->w1 : (cfg->w1 + cfg->w2));
        uint32_t live_bar_b = cfg->split_c2 ? w_res : (w_res + cfg->w1 + cfg->w2);
        if (peak_pool < live_bar_b + 10u) peak_pool = live_bar_b + 10u;

        if (!cfg->split_c2) {
            qua_bitop_qroam_res_t c12_qroam;
            qua_bitop_qroam_vslot0_resources(c12_tab, cfg->n_charts, 1u, n_addr, cfg->w1 + cfg->w2, 0u, 0u, 0u, &c12_qroam);
            if (peak_pool < live_c1 + c12_qroam.anc_needed) peak_pool = live_c1 + c12_qroam.anc_needed;
            qua_bitop_hmr_lut_res_t c12_del;
            qua_bitop_hmr_parity_lut_resources(cfg->n_charts, n_addr, cfg->w1 + cfg->w2,
                                               (peak_pool > live_c1) ? (peak_pool - live_c1) : 0u, &c12_del);
            if (peak_pool < live_c1 + c12_del.pool_anc_needed) peak_pool = live_c1 + c12_del.pool_anc_needed;
        } else {
            qua_bitop_qroam_res_t c2_qroam0, c1_qroam0;
            qua_bitop_qroam_vslot0_resources(c2_tab, cfg->n_charts, 1u, n_addr, cfg->w2, 0u, 0u, 0u, &c2_qroam0);
            qua_bitop_qroam_vslot0_resources(c1_tab, cfg->n_charts, 1u, n_addr, cfg->w1, 0u, 0u, 0u, &c1_qroam0);
            if (peak_pool < live_c2 + c2_qroam0.anc_needed) peak_pool = live_c2 + c2_qroam0.anc_needed;
            if (peak_pool < live_c1 + c1_qroam0.anc_needed) peak_pool = live_c1 + c1_qroam0.anc_needed;
            qua_bitop_hmr_lut_res_t c2_del0, c1_del0;
            qua_bitop_hmr_parity_lut_resources(cfg->n_charts, n_addr, cfg->w2,
                                               (peak_pool > live_c2) ? (peak_pool - live_c2) : 0u, &c2_del0);
            if (peak_pool < live_c2 + c2_del0.pool_anc_needed) peak_pool = live_c2 + c2_del0.pool_anc_needed;
            qua_bitop_hmr_parity_lut_resources(cfg->n_charts, n_addr, cfg->w1,
                                               (peak_pool > live_c1) ? (peak_pool - live_c1) : 0u, &c1_del0);
            if (peak_pool < live_c1 + c1_del0.pool_anc_needed) peak_pool = live_c1 + c1_del0.pool_anc_needed;
        }
        uint32_t qroam_slack = b_sh_main * 5u;
        if (peak_pool < live_c1 + qroam_slack) peak_pool = live_c1 + qroam_slack;

        uint32_t live_sq = live_c2 + (use_arith ? (cfg->wsq + 1u) : cfg->wsq);
        if (use_arith) {
            if (peak_pool < live_c2 + 2u * cfg->wsq + 2u) peak_pool = live_c2 + 2u * cfg->wsq + 2u;
        } else {
            qua_bitop_qroam_res_t sq_qroam0;
            qua_bitop_qroam_vslot0_resources(sq_tab64, cfg->sq_n, 1u, sq_addr_n, cfg->wsq, b_sh_sq, 0u, 0u, &sq_qroam0);
            if (peak_pool < live_sq + sq_qroam0.anc_needed) peak_pool = live_sq + sq_qroam0.anc_needed;
            qua_bitop_hmr_lut_res_t sq_del0;
            qua_bitop_hmr_parity_lut_resources(cfg->sq_n, sq_addr_n, cfg->wsq,
                                               (peak_pool > live_sq) ? (peak_pool - live_sq) : 0u, &sq_del0);
            if (peak_pool < live_sq + sq_del0.pool_anc_needed) peak_pool = live_sq + sq_del0.pool_anc_needed;
        }

        for (uint32_t b = 0u; b < cfg->wsq; ++b) {
            if (cfg->sq_row_hi[b] == 0u) continue;
            int32_t rsh = cfg->sh2 - (int32_t)b;
            uint32_t n_src = (rsh >= 0) ? (cfg->w2 - (uint32_t)rsh) : cfg->w2;
            uint32_t lsh   = (rsh >= 0) ? 0u : (uint32_t)(-rsh);
            uint32_t n_t   = cfg->sq_row_hi[b] - lsh;
            uint32_t min_a = qua_inv_min_ctrl_add_anc(n_src, n_t, cfg->zero_prods_mode != 0u, cfg->cuccaro_mode != 0u);
            if (peak_pool < live_sq + min_a) peak_pool = live_sq + min_a;
        }
        for (uint32_t b = 0u; b < cfg->wd1; ++b) {
            if (cfg->c1_row_hi[b] == 0u) continue;
            int32_t rsh = cfg->sh1 - (int32_t)b;
            uint32_t n_src = (rsh >= 0) ? (cfg->w1 - (uint32_t)rsh) : cfg->w1;
            uint32_t lsh   = (rsh >= 0) ? 0u : (uint32_t)(-rsh);
            uint32_t n_t   = cfg->c1_row_hi[b] - lsh;
            uint32_t min_a = qua_inv_min_ctrl_add_anc(n_src, n_t, cfg->zero_prods_mode != 0u, cfg->cuccaro_mode != 0u);
            if (peak_pool < live_c1 + min_a) peak_pool = live_c1 + min_a;
        }

        cfg->peak_qubits = 27u + peak_pool;

        /* Pass 2: Exact Toffoli count at finalized peak_pool */
        uint64_t tof = 206ull; /* qua_inv_log_clz forward (105) + reverse (101) */

        /* Stage A: c0 load + barrel(w0) + c0 HMR_Dual erase */
        uint32_t bsh_c0 = qua_invsqrt_pick_sub_bsh(c0_tab, cfg->n_charts, n_addr, cfg->w0, peak_pool - cfg->w0, &c0_qroam);
        (void)bsh_c0;
        qua_bitop_hmr_parity_lut_resources(cfg->n_charts, n_addr, cfg->w0, peak_pool - cfg->w0, &c0_del);
        tof += c0_qroam.total_ccx + c0_del.total_ccx;
        tof += 2ull * (uint64_t)(4u * cfg->w0 + 6u);

        /* Stage B: coeff QROAM + HMR_Dual */
        if (!cfg->split_c2) {
            qua_bitop_qroam_res_t c12_qroam;
            qua_bitop_hmr_lut_res_t c12_del;
            (void)qua_invsqrt_pick_sub_bsh(c12_tab, cfg->n_charts, n_addr, cfg->w1 + cfg->w2, peak_pool - live_c1, &c12_qroam);
            qua_bitop_hmr_parity_lut_resources(cfg->n_charts, n_addr, cfg->w1 + cfg->w2,
                                               peak_pool - live_c1, &c12_del);
            tof += c12_qroam.total_ccx + c12_del.total_ccx;
        } else {
            qua_bitop_qroam_res_t c2_res, c1_res;
            qua_bitop_hmr_lut_res_t c2_del, c1_del;
            qua_invsqrt_pick_sub_bsh(c2_tab, cfg->n_charts, n_addr, cfg->w2, peak_pool - live_c2, &c2_res);
            qua_invsqrt_pick_sub_bsh(c1_tab, cfg->n_charts, n_addr, cfg->w1, peak_pool - live_c1, &c1_res);
            qua_bitop_hmr_parity_lut_resources(cfg->n_charts, n_addr, cfg->w2, peak_pool - live_c2, &c2_del);
            qua_bitop_hmr_parity_lut_resources(cfg->n_charts, n_addr, cfg->w1, peak_pool - live_c1, &c1_del);
            /* Both c2 and c1 are loaded twice and HMR_Dual-erased twice across Stage B! */
            tof += 2ull * (c2_res.total_ccx + c2_del.total_ccx + c1_res.total_ccx + c1_del.total_ccx);
        }

        /* Stage B: sq load + erase (forward & reverse) */
        if (use_arith) {
            tof += 4ull * qua_invsqrt_sq_arith_once_ccx(cfg);
        } else {
            qua_bitop_qroam_res_t sq_qroam;
            qua_bitop_hmr_lut_res_t sq_del;
            qua_bitop_qroam_vslot0_resources(sq_tab64, cfg->sq_n, 1u, sq_addr_n, cfg->wsq, b_sh_sq, 0u,
                                             peak_pool - live_sq, &sq_qroam);
            qua_bitop_hmr_parity_lut_resources(cfg->sq_n, sq_addr_n, cfg->wsq, peak_pool - live_sq, &sq_del);
            tof += 2ull * (sq_qroam.total_ccx + sq_del.total_ccx);
        }

        /* Stage B: c2*sq and c1*t_inv controlled additions (forward + reverse) */
        bool allow_ps  = (cfg->zero_prods_mode == 0u);
        bool allow_cuc = (cfg->cuccaro_mode != 0u);
        for (uint32_t b = 0u; b < cfg->wsq; ++b) {
            if (cfg->sq_row_hi[b] == 0u) continue;
            int32_t rsh = cfg->sh2 - (int32_t)b;
            uint32_t n_src = (rsh >= 0) ? (cfg->w2 - (uint32_t)rsh) : cfg->w2;
            uint32_t lsh   = (rsh >= 0) ? 0u : (uint32_t)(-rsh);
            uint32_t n_t   = cfg->sq_row_hi[b] - lsh;
            uint32_t row_a = 0u;
            uint64_t row_c = 0u;
            (void)qua_inv_select_ctrl_add(n_src, n_t, peak_pool - live_sq, allow_ps, allow_cuc, &row_a, &row_c);
            tof += 2ull * row_c;
        }
        for (uint32_t b = 0u; b < cfg->wd1; ++b) {
            if (cfg->c1_row_hi[b] == 0u) continue;
            int32_t rsh = cfg->sh1 - (int32_t)b;
            uint32_t n_src = (rsh >= 0) ? (cfg->w1 - (uint32_t)rsh) : cfg->w1;
            uint32_t lsh   = (rsh >= 0) ? 0u : (uint32_t)(-rsh);
            uint32_t n_t   = cfg->c1_row_hi[b] - lsh;
            uint32_t row_a = 0u;
            uint64_t row_c = 0u;
            (void)qua_inv_select_ctrl_add(n_src, n_t, peak_pool - live_c1, allow_ps, allow_cuc, &row_a, &row_c);
            tof += 2ull * row_c;
        }

        /* Stage B barrel denorm(w_res) forward + reverse */
        tof += 2ull * (uint64_t)(4u * w_res + 6u);

        cfg->est_toffolis = tof;
        cfg->est_score    = (double)cfg->peak_qubits * ((double)tof + cfg->rot_c);
        return;
    }

    uint32_t coeff_q = cfg->split_c2 ? ((cfg->w2 > cfg->w1) ? cfg->w2 : cfg->w1) : (cfg->w1 + cfg->w2);
    uint32_t base_q  = 22u + 5u + cfg->w0 + coeff_q + 3u;

    /* Baseline pool watermark from log_clz (14), y_ext (10), and carry flushes (max(m2-m1, m3-m2, w0-m3-2)) */
    uint32_t peak_anc = 14u;
    uint32_t inc1 = (cfg->m2 > cfg->m1) ? (cfg->m2 - cfg->m1) : 0u;
    uint32_t inc2 = (cfg->m3 > cfg->m2) ? (cfg->m3 - cfg->m2) : 0u;
    uint32_t inc3 = (cfg->w0 > cfg->m3 + 2u) ? (cfg->w0 - cfg->m3 - 2u) : 0u;
    if (peak_anc < inc1) peak_anc = inc1;
    if (peak_anc < inc2) peak_anc = inc2;
    if (peak_anc < inc3) peak_anc = inc3;

    uint64_t main_tab[QUA_INVSQRT_MAX_CHARTS];
    uint64_t c1_tab[QUA_INVSQRT_MAX_CHARTS];
    uint64_t c2_tab[QUA_INVSQRT_MAX_CHARTS];
    uint32_t main_w = cfg->split_c2 ? (cfg->w0 + cfg->w2) : (cfg->w0 + cfg->w1 + cfg->w2);
    for (uint32_t c = 0u; c < cfg->n_charts; ++c) {
        if (!cfg->split_c2) {
            main_tab[c] = (uint64_t)cfg->c0_biased[c] |
                          ((uint64_t)cfg->c1[c] << cfg->w0) |
                          ((uint64_t)cfg->c2[c] << (cfg->w0 + cfg->w1));
        } else {
            main_tab[c] = (uint64_t)cfg->c0_biased[c] |
                          ((uint64_t)cfg->c2[c] << cfg->w0);
            c1_tab[c] = (uint64_t)cfg->c1[c];
            c2_tab[c] = (uint64_t)cfg->c2[c];
        }
    }

    uint32_t sq_hold = use_arith ? (cfg->wsq + 1u) : cfg->wsq;
    qua_bitop_qroam_res_t sq_qroam_res;
    if (use_arith) {
        uint32_t sq_ar_anc = 2u * cfg->wsq + 2u;
        if (peak_anc < sq_ar_anc) peak_anc = sq_ar_anc;
    } else {
        qua_bitop_qroam_vslot0_resources(sq_tab64, cfg->sq_n, 1u, sq_addr_n, cfg->wsq, b_sh_sq, 0u, 0u, &sq_qroam_res);
        if (peak_anc < cfg->wsq + sq_qroam_res.anc_needed) peak_anc = cfg->wsq + sq_qroam_res.anc_needed;
        qua_bitop_hmr_lut_res_t sq_del0;
        qua_bitop_hmr_parity_lut_resources(cfg->sq_n, sq_addr_n, cfg->wsq,
                                           (peak_anc > cfg->wsq) ? (peak_anc - cfg->wsq) : 0u, &sq_del0);
        if (peak_anc < cfg->wsq + sq_del0.pool_anc_needed) peak_anc = cfg->wsq + sq_del0.pool_anc_needed;
    }

    for (uint32_t b = 0u; b < cfg->wsq; ++b) {
        int32_t rsh = cfg->sh2 - (int32_t)b;
        if (rsh >= (int32_t)cfg->w2) continue;
        uint32_t n_src = (rsh >= 0) ? (cfg->w2 - (uint32_t)rsh) : cfg->w2;
        uint32_t n_t   = (rsh >= 0) ? (cfg->m1 + 1u) : (cfg->m1 + 1u - (uint32_t)(-rsh));
        uint32_t min_a = qua_inv_min_ctrl_add_anc(n_src, n_t, cfg->zero_prods_mode != 0u, cfg->cuccaro_mode != 0u);
        if (peak_anc < sq_hold + min_a) peak_anc = sq_hold + min_a;
    }
    for (uint32_t ti = 0u; ti < 3u; ++ti) {
        const uint32_t *arr = (ti == 0u) ? cfg->tier1_b : ((ti == 1u) ? cfg->tier2_b : cfg->tier3_b);
        uint32_t cnt = (ti == 0u) ? cfg->tier1_n : ((ti == 1u) ? cfg->tier2_n : cfg->tier3_n);
        uint32_t mt  = (ti == 0u) ? cfg->m1 : ((ti == 1u) ? cfg->m2 : cfg->m3);
        for (uint32_t j = 0u; j < cnt; ++j) {
            int32_t rsh = cfg->sh1 - (int32_t)arr[j];
            if (rsh >= (int32_t)cfg->w1) continue;
            uint32_t n_src = (rsh >= 0) ? (cfg->w1 - (uint32_t)rsh) : cfg->w1;
            uint32_t n_t   = (rsh >= 0) ? (mt + 1u) : (mt + 1u - (uint32_t)(-rsh));
            uint32_t min_a = qua_inv_min_ctrl_add_anc(n_src, n_t, cfg->zero_prods_mode != 0u, cfg->cuccaro_mode != 0u);
            if (peak_anc < min_a) peak_anc = min_a;
        }
    }

    qua_bitop_qroam_res_t main_qroam0;
    qua_bitop_qroam_vslot0_resources(main_tab, cfg->n_charts, 1u, n_addr, main_w, 0u, 0u, 0u, &main_qroam0);
    if (peak_anc < main_qroam0.anc_needed) peak_anc = main_qroam0.anc_needed;
    uint32_t qroam_slack = b_sh_main * 5u;
    if (peak_anc < 14u + qroam_slack) peak_anc = 14u + qroam_slack;

    qua_bitop_hmr_lut_res_t main_del0;
    qua_bitop_hmr_parity_lut_resources(cfg->n_charts, n_addr, main_w, peak_anc, &main_del0);
    if (peak_anc < main_del0.pool_anc_needed) peak_anc = main_del0.pool_anc_needed;

    if (cfg->split_c2) {
        qua_bitop_qroam_res_t c1_q0, c2_q0;
        qua_bitop_qroam_vslot0_resources(c1_tab, cfg->n_charts, 1u, n_addr, cfg->w1, 0u, 0u, 0u, &c1_q0);
        qua_bitop_qroam_vslot0_resources(c2_tab, cfg->n_charts, 1u, n_addr, cfg->w2, 0u, 0u, 0u, &c2_q0);
        if (peak_anc < c1_q0.anc_needed) peak_anc = c1_q0.anc_needed;
        if (peak_anc < c2_q0.anc_needed) peak_anc = c2_q0.anc_needed;
        qua_bitop_hmr_lut_res_t c2_d0, c1_d0;
        qua_bitop_hmr_parity_lut_resources(cfg->n_charts, n_addr, cfg->w2, peak_anc, &c2_d0);
        qua_bitop_hmr_parity_lut_resources(cfg->n_charts, n_addr, cfg->w1, peak_anc, &c1_d0);
        if (peak_anc < c2_d0.pool_anc_needed) peak_anc = c2_d0.pool_anc_needed;
        if (peak_anc < c1_d0.pool_anc_needed) peak_anc = c1_d0.pool_anc_needed;
    }

    cfg->peak_qubits = base_q + peak_anc;

    /* Exact Toffoli model matching `qua_invsqrt_emit_kmb` */
    uint64_t tof = 206ull; /* log_clz forward (105) + reverse (101) */
    qua_bitop_qroam_res_t main_qroam_res;
    qua_bitop_hmr_lut_res_t main_del_res;
    (void)qua_invsqrt_pick_sub_bsh(main_tab, cfg->n_charts, n_addr, main_w, peak_anc, &main_qroam_res);
    qua_bitop_hmr_parity_lut_resources(cfg->n_charts, n_addr, main_w, peak_anc, &main_del_res);
    tof += main_qroam_res.total_ccx + main_del_res.total_ccx;
    if (cfg->split_c2) {
        qua_bitop_qroam_res_t c1_res, c2_res;
        qua_bitop_hmr_lut_res_t c1_del, c2_del;
        qua_invsqrt_pick_sub_bsh(c1_tab, cfg->n_charts, n_addr, cfg->w1, peak_anc, &c1_res);
        qua_invsqrt_pick_sub_bsh(c2_tab, cfg->n_charts, n_addr, cfg->w2, peak_anc, &c2_res);
        qua_bitop_hmr_parity_lut_resources(cfg->n_charts, n_addr, cfg->w1, peak_anc, &c1_del);
        qua_bitop_hmr_parity_lut_resources(cfg->n_charts, n_addr, cfg->w2, peak_anc, &c2_del);
        tof += c2_del.total_ccx + c1_res.total_ccx + c1_del.total_ccx + c2_res.total_ccx;
    }
    if (use_arith) {
        tof += 4ull * qua_invsqrt_sq_arith_once_ccx(cfg);
    } else {
        qua_bitop_hmr_lut_res_t sq_del_res;
        qua_bitop_qroam_vslot0_resources(sq_tab64, cfg->sq_n, 1u, sq_addr_n, cfg->wsq, b_sh_sq, 0u,
                                         peak_anc - cfg->wsq, &sq_qroam_res);
        qua_bitop_hmr_parity_lut_resources(cfg->sq_n, sq_addr_n, cfg->wsq, peak_anc - cfg->wsq, &sq_del_res);
        tof += 2ull * (sq_qroam_res.total_ccx + sq_del_res.total_ccx);
    }
    bool allow_ps  = (cfg->zero_prods_mode == 0u);
    bool allow_cuc = (cfg->cuccaro_mode != 0u);
    for (uint32_t b = 0u; b < cfg->wsq; ++b) {
        int32_t rsh = cfg->sh2 - (int32_t)b;
        if (rsh >= (int32_t)cfg->w2) continue;
        uint32_t n_src = (rsh >= 0) ? (cfg->w2 - (uint32_t)rsh) : cfg->w2;
        uint32_t n_t   = (rsh >= 0) ? (cfg->m1 + 1u) : (cfg->m1 + 1u - (uint32_t)(-rsh));
        uint32_t row_a = 0u;
        uint64_t row_c = 0u;
        (void)qua_inv_select_ctrl_add(n_src, n_t, peak_anc - sq_hold, allow_ps, allow_cuc, &row_a, &row_c);
        tof += 2ull * row_c;
    }
    for (uint32_t ti = 0u; ti < 3u; ++ti) {
        const uint32_t *arr = (ti == 0u) ? cfg->tier1_b : ((ti == 1u) ? cfg->tier2_b : cfg->tier3_b);
        uint32_t cnt = (ti == 0u) ? cfg->tier1_n : ((ti == 1u) ? cfg->tier2_n : cfg->tier3_n);
        uint32_t mt  = (ti == 0u) ? cfg->m1 : ((ti == 1u) ? cfg->m2 : cfg->m3);
        for (uint32_t j = 0u; j < cnt; ++j) {
            int32_t rsh = cfg->sh1 - (int32_t)arr[j];
            if (rsh >= (int32_t)cfg->w1) continue;
            uint32_t n_src = (rsh >= 0) ? (cfg->w1 - (uint32_t)rsh) : cfg->w1;
            uint32_t n_t   = (rsh >= 0) ? (mt + 1u) : (mt + 1u - (uint32_t)(-rsh));
            uint32_t row_a = 0u;
            uint64_t row_c = 0u;
            (void)qua_inv_select_ctrl_add(n_src, n_t, peak_anc, allow_ps, allow_cuc, &row_a, &row_c);
            tof += 2ull * row_c;
        }
    }
    tof += 2ull * (uint64_t)((cfg->m2 - cfg->m1) + (cfg->m3 - cfg->m2) + (cfg->w0 - cfg->m3 - 1u));
    tof += 2ull * (uint64_t)(4u * cfg->w0 + 6u);
    cfg->est_toffolis = tof;
    cfg->est_score    = (double)cfg->peak_qubits * ((double)tof + cfg->rot_c);
}

/* ============================================================================
 * Pure-C17 Kickmix v2 (.kmb) Gate Emitter for `qua_invsqrt_cfg_t`
 * ============================================================================ */
typedef struct {
    uint32_t free_list[256];
    uint32_t free_top;
    uint32_t start_q;
    uint32_t next_q;
    uint32_t next_bit;
} qua_qpool_t;

static inline void qua_qpool_init(qua_qpool_t *p, uint32_t start_q) {
    p->free_top = 0u;
    p->start_q  = start_q;
    p->next_q   = start_q;
    p->next_bit = 0u;
}

static inline uint32_t qua_qpool_avail(const qua_qpool_t *p, uint32_t peak_qubits) {
    uint32_t peak_pool = (peak_qubits > p->start_q) ? (peak_qubits - p->start_q) : 0u;
    uint32_t live      = (p->next_q - p->start_q) - p->free_top;
    return (peak_pool > live) ? (peak_pool - live) : 0u;
}

static inline uint32_t qua_qpool_alloc1(qua_qpool_t *p) {
    if (p->free_top > 0u) return p->free_list[--p->free_top];
    return p->next_q++;
}

static inline void qua_qpool_alloc_n(qua_qpool_t *p, uint32_t *out, uint32_t n) {
    for (uint32_t i = 0u; i < n; ++i) out[i] = qua_qpool_alloc1(p);
}

static inline void qua_qpool_free1(qua_qpool_t *p, uint32_t q) {
    p->free_list[p->free_top++] = q;
}

static inline void qua_qpool_free_n(qua_qpool_t *p, const uint32_t *qs, uint32_t n) {
    for (uint32_t i = n; i > 0u; --i) qua_qpool_free1(p, qs[i - 1u]);
}

static inline void qua_inv_multi_zero_flag(qua_builder_t *kb, qua_qpool_t *pool, bool inverse,
                                           const uint32_t *wires, uint32_t n, uint32_t flag) {
    if (n == 1u) {
        qua_x(kb, wires[0]);
        qua_cx(kb, wires[0], flag);
        qua_x(kb, wires[0]);
        return;
    }
    for (uint32_t i = 0u; i < n; ++i) qua_x(kb, wires[i]);
    if (n == 2u) {
        if (inverse) qua_and_uncompute(kb, wires[0], wires[1], flag, pool->next_bit++);
        else         qua_and_compute(kb, wires[0], wires[1], flag);
    } else {
        uint32_t chain[32];
        qua_qpool_alloc_n(pool, chain, n - 2u);
        qua_and_compute(kb, wires[0], wires[1], chain[0]);
        for (uint32_t i = 2u; i < n - 1u; ++i) {
            qua_and_compute(kb, chain[i - 2u], wires[i], chain[i - 1u]);
        }
        if (inverse) qua_and_uncompute(kb, chain[n - 3u], wires[n - 1u], flag, pool->next_bit++);
        else         qua_and_compute(kb, chain[n - 3u], wires[n - 1u], flag);
        for (uint32_t i = n - 2u; i >= 2u; --i) {
            qua_and_uncompute(kb, chain[i - 2u], wires[i], chain[i - 1u], pool->next_bit++);
        }
        qua_and_uncompute(kb, wires[0], wires[1], chain[0], pool->next_bit++);
        qua_qpool_free_n(pool, chain, n - 2u);
    }
    for (uint32_t i = 0u; i < n; ++i) qua_x(kb, wires[i]);
}

static inline void qua_inv_log_clz(qua_builder_t *kb, qua_qpool_t *pool, bool inverse,
                                   const uint32_t *x, uint32_t W, const uint32_t *clz) {
    static const uint32_t bit_idx[5] = {4u, 3u, 2u, 1u, 0u};
    static const uint32_t sh[5]      = {16u, 8u, 4u, 2u, 1u};
    if (!inverse) {
        for (uint32_t s = 0u; s < 5u; ++s) {
            uint32_t ctrl = clz[bit_idx[s]];
            qua_inv_multi_zero_flag(kb, pool, false, x + (W - sh[s]), sh[s], ctrl);
            for (int32_t i = (int32_t)(W - 1u - sh[s]); i >= 0; --i) {
                qua_cswap(kb, ctrl, x[i], x[(uint32_t)i + sh[s]]);
            }
        }
    } else {
        for (int32_t s = 4; s >= 0; --s) {
            uint32_t ctrl = clz[bit_idx[s]];
            for (uint32_t i = 0u; i < W - sh[s]; ++i) {
                qua_cswap(kb, ctrl, x[i], x[i + sh[s]]);
            }
            qua_inv_multi_zero_flag(kb, pool, true, x + (W - sh[s]), sh[s], ctrl);
        }
    }
}

static inline void qua_inv_emit_qroam_u64(
    qua_builder_t *kb,
    qua_qpool_t *pool,
    const uint64_t *table,
    uint32_t n_entries,
    const uint32_t *addr_lsb,
    uint32_t n_addr,
    const uint32_t *tgt,
    uint32_t n_tgt,
    uint32_t b_sh,
    uint32_t free_anc_cap
) {
    if (b_sh >= n_addr) b_sh = 0u;
    qua_bitop_qroam_res_t res;
    qua_bitop_qroam_vslot0_resources(table, n_entries, 1u, n_addr, n_tgt, b_sh, 0u, free_anc_cap, &res);
    uint32_t pool_len = (res.anc_needed > free_anc_cap) ? res.anc_needed : free_anc_cap;
    uint32_t qanc[512];
    if (pool_len > 512u) abort();
    qua_qpool_alloc_n(pool, qanc, pool_len);
    qua_bitop_emit_qroam_vslot0(kb, table, n_entries, 1u, addr_lsb, n_addr,
                                tgt, n_tgt, b_sh, 0u, qanc, pool_len, &pool->next_bit);
    qua_qpool_free_n(pool, qanc, pool_len);
}

/* Composable Combinator 1 (`HMR_Dual`): Ket[w] -> ParityDual -> |0>^w */
static inline void qua_inv_qrom_del(qua_builder_t *kb, qua_qpool_t *pool,
                                    const uint64_t *table, uint32_t n_entries,
                                    const uint32_t *addr_lsb, uint32_t n_addr,
                                    const uint32_t *tgt, uint32_t n_tgt,
                                    uint32_t free_anc_cap) {
    if (n_tgt == 0u) return;
    qua_bitop_hmr_lut_res_t res;
    qua_bitop_hmr_parity_lut_resources(n_entries, n_addr, n_tgt, free_anc_cap, &res);
    uint32_t extra_anc[128];
    if (res.pool_anc_needed > 0u) {
        qua_qpool_alloc_n(pool, extra_anc, res.pool_anc_needed);
    }
    qua_bitop_emit_hmr_parity_lut(kb, table, n_entries, addr_lsb, n_addr,
                                  tgt, n_tgt, extra_anc, res.pool_anc_needed, &pool->next_bit);
    if (res.pool_anc_needed > 0u) {
        qua_qpool_free_n(pool, extra_anc, res.pool_anc_needed);
    }
}

static inline void qua_inv_gidney_inc(qua_builder_t *kb, qua_qpool_t *pool, uint32_t c_in,
                                      const uint32_t *tgt, uint32_t m, bool subtract) {
    if (m == 0u) return;
    if (subtract) {
        for (uint32_t i = 0u; i < m; ++i) qua_x(kb, tgt[i]);
    }
    if (m == 1u) {
        qua_cx(kb, c_in, tgt[0]);
    } else {
        uint32_t p[64];
        qua_qpool_alloc_n(pool, p, m - 1u);
        qua_and_compute(kb, c_in, tgt[0], p[0]);
        for (uint32_t j = 1u; j < m - 1u; ++j) {
            qua_and_compute(kb, p[j - 1u], tgt[j], p[j]);
        }
        qua_cx(kb, p[m - 2u], tgt[m - 1u]);
        for (uint32_t j = m - 2u; j >= 1u; --j) {
            qua_and_uncompute(kb, p[j - 1u], tgt[j], p[j], pool->next_bit++);
            qua_cx(kb, p[j - 1u], tgt[j]);
        }
        qua_and_uncompute(kb, c_in, tgt[0], p[0], pool->next_bit++);
        qua_cx(kb, c_in, tgt[0]);
        qua_qpool_free_n(pool, p, m - 1u);
    }
    if (subtract) {
        for (uint32_t i = 0u; i < m; ++i) qua_x(kb, tgt[i]);
    }
}

static inline void qua_inv_gidney_iadd_short(qua_builder_t *kb, qua_qpool_t *pool,
                                             const uint32_t *lhs, uint32_t k,
                                             const uint32_t *tgt, uint32_t n) {
    if (n == 1u) {
        qua_cx(kb, lhs[0], tgt[0]);
        return;
    }
    uint32_t carries[64];
    if (k == n) {
        qua_qpool_alloc_n(pool, carries, n - 1u);
        qua_and_compute(kb, lhs[0], tgt[0], carries[0]);
        for (uint32_t i = 1u; i < n - 1u; ++i) {
            qua_cx(kb, carries[i - 1u], lhs[i]);
            qua_cx(kb, carries[i - 1u], tgt[i]);
            qua_and_compute(kb, lhs[i], tgt[i], carries[i]);
            qua_cx(kb, carries[i - 1u], carries[i]);
        }
        qua_cx(kb, carries[n - 2u], tgt[n - 1u]);
        qua_cx(kb, lhs[n - 1u], tgt[n - 1u]);
        for (uint32_t i = n - 2u; i >= 1u; --i) {
            qua_cx(kb, carries[i - 1u], carries[i]);
            qua_and_uncompute(kb, lhs[i], tgt[i], carries[i], pool->next_bit++);
            qua_cx(kb, carries[i - 1u], lhs[i]);
            qua_cx(kb, lhs[i], tgt[i]);
        }
        qua_and_uncompute(kb, lhs[0], tgt[0], carries[0], pool->next_bit++);
        qua_cx(kb, lhs[0], tgt[0]);
        qua_qpool_free_n(pool, carries, n - 1u);
        return;
    }
    qua_qpool_alloc_n(pool, carries, k);
    qua_and_compute(kb, lhs[0], tgt[0], carries[0]);
    for (uint32_t i = 1u; i < k; ++i) {
        qua_cx(kb, carries[i - 1u], lhs[i]);
        qua_cx(kb, carries[i - 1u], tgt[i]);
        qua_and_compute(kb, lhs[i], tgt[i], carries[i]);
        qua_cx(kb, carries[i - 1u], carries[i]);
    }
    qua_inv_gidney_inc(kb, pool, carries[k - 1u], tgt + k, n - k, false);
    for (uint32_t i = k - 1u; i >= 1u; --i) {
        qua_cx(kb, carries[i - 1u], carries[i]);
        qua_and_uncompute(kb, lhs[i], tgt[i], carries[i], pool->next_bit++);
        qua_cx(kb, carries[i - 1u], lhs[i]);
        qua_cx(kb, lhs[i], tgt[i]);
    }
    qua_and_uncompute(kb, lhs[0], tgt[0], carries[0], pool->next_bit++);
    qua_cx(kb, lhs[0], tgt[0]);
    qua_qpool_free_n(pool, carries, k);
}

static inline void qua_inv_ctrl_iadd_zero_prods(qua_builder_t *kb, qua_qpool_t *pool,
                                                uint32_t ctrl, const uint32_t *lhs, uint32_t k,
                                                const uint32_t *tgt, uint32_t n) {
    if (k == 1u) {
        uint32_t p0 = qua_qpool_alloc1(pool);
        qua_and_compute(kb, ctrl, lhs[0], p0);
        qua_inv_gidney_inc(kb, pool, p0, tgt, n, false);
        qua_and_uncompute(kb, ctrl, lhs[0], p0, pool->next_bit++);
        qua_qpool_free1(pool, p0);
        return;
    }
    uint32_t carries[64];
    if (k == n) {
        qua_qpool_alloc_n(pool, carries, n - 1u);
        qua_and_compute(kb, lhs[0], tgt[0], carries[0]);
        for (uint32_t i = 1u; i < n - 1u; ++i) {
            qua_cx(kb, carries[i - 1u], lhs[i]);
            qua_cx(kb, carries[i - 1u], tgt[i]);
            qua_and_compute(kb, lhs[i], tgt[i], carries[i]);
            qua_cx(kb, carries[i - 1u], carries[i]);
        }
        qua_cx(kb, carries[n - 2u], lhs[n - 1u]);
        uint32_t tmp = qua_qpool_alloc1(pool);
        qua_and_compute(kb, ctrl, lhs[n - 1u], tmp);
        qua_cx(kb, tmp, tgt[n - 1u]);
        qua_and_uncompute(kb, ctrl, lhs[n - 1u], tmp, pool->next_bit++);
        qua_qpool_free1(pool, tmp);
        qua_cx(kb, carries[n - 2u], lhs[n - 1u]);

        for (uint32_t i = n - 2u; i >= 1u; --i) {
            qua_cx(kb, carries[i - 1u], carries[i]);
            qua_and_uncompute(kb, lhs[i], tgt[i], carries[i], pool->next_bit++);
            qua_cx(kb, carries[i - 1u], tgt[i]);
            qua_and_compute(kb, ctrl, lhs[i], carries[i]);
            qua_cx(kb, carries[i], tgt[i]);
            qua_and_uncompute(kb, ctrl, lhs[i], carries[i], pool->next_bit++);
            qua_cx(kb, carries[i - 1u], lhs[i]);
        }
        qua_and_uncompute(kb, lhs[0], tgt[0], carries[0], pool->next_bit++);
        qua_and_compute(kb, ctrl, lhs[0], carries[0]);
        qua_cx(kb, carries[0], tgt[0]);
        qua_and_uncompute(kb, ctrl, lhs[0], carries[0], pool->next_bit++);
        qua_qpool_free_n(pool, carries, n - 1u);
        return;
    }
    qua_qpool_alloc_n(pool, carries, k);
    qua_and_compute(kb, lhs[0], tgt[0], carries[0]);
    for (uint32_t i = 1u; i < k; ++i) {
        qua_cx(kb, carries[i - 1u], lhs[i]);
        qua_cx(kb, carries[i - 1u], tgt[i]);
        qua_and_compute(kb, lhs[i], tgt[i], carries[i]);
        qua_cx(kb, carries[i - 1u], carries[i]);
    }
    uint32_t c_ctrl = qua_qpool_alloc1(pool);
    qua_and_compute(kb, ctrl, carries[k - 1u], c_ctrl);
    qua_inv_gidney_inc(kb, pool, c_ctrl, tgt + k, n - k, false);
    qua_and_uncompute(kb, ctrl, carries[k - 1u], c_ctrl, pool->next_bit++);
    qua_qpool_free1(pool, c_ctrl);

    for (uint32_t i = k - 1u; i >= 1u; --i) {
        qua_cx(kb, carries[i - 1u], carries[i]);
        qua_and_uncompute(kb, lhs[i], tgt[i], carries[i], pool->next_bit++);
        qua_cx(kb, carries[i - 1u], tgt[i]);
        qua_and_compute(kb, ctrl, lhs[i], carries[i]);
        qua_cx(kb, carries[i], tgt[i]);
        qua_and_uncompute(kb, ctrl, lhs[i], carries[i], pool->next_bit++);
        qua_cx(kb, carries[i - 1u], lhs[i]);
    }
    qua_and_uncompute(kb, lhs[0], tgt[0], carries[0], pool->next_bit++);
    qua_and_compute(kb, ctrl, lhs[0], carries[0]);
    qua_cx(kb, carries[0], tgt[0]);
    qua_and_uncompute(kb, ctrl, lhs[0], carries[0], pool->next_bit++);
    qua_qpool_free_n(pool, carries, k);
}

static inline void qua_inv_hybrid_ctrl_add(qua_builder_t *kb, qua_qpool_t *pool,
                                           uint32_t ctrl, const uint32_t *srcs, uint32_t n_src,
                                           const uint32_t *tgt, uint32_t n_tgt,
                                           bool subtract, uint32_t avail_anc,
                                           bool allow_prods_short, bool allow_cuccaro) {
    if (n_src == 0u || n_tgt == 0u) return;
    uint32_t k = (n_src < n_tgt) ? n_src : n_tgt;
    uint32_t anc_need = 0u;
    uint64_t ccx_est  = 0u;
    qua_inv_add_kind_t kind = qua_inv_select_ctrl_add(
        n_src, n_tgt, avail_anc, allow_prods_short, allow_cuccaro, &anc_need, &ccx_est);
    if (kind == QUA_INV_ADD_CUCCARO) {
        uint32_t anc[64];
        if (anc_need > 0u) qua_qpool_alloc_n(pool, anc, anc_need);
        qua_bitop_emit_ctrl_iadd_cuccaro(kb, ctrl, srcs, k, tgt, n_tgt, subtract, anc, &pool->next_bit);
        if (anc_need > 0u) qua_qpool_free_n(pool, anc, anc_need);
        return;
    }
    if (subtract) {
        for (uint32_t i = 0u; i < n_tgt; ++i) qua_x(kb, tgt[i]);
    }
    if (kind == QUA_INV_ADD_PRODS_SHORT) {
        uint32_t prods[64];
        qua_qpool_alloc_n(pool, prods, k);
        for (uint32_t i = 0u; i < k; ++i) qua_and_compute(kb, ctrl, srcs[i], prods[i]);
        qua_inv_gidney_iadd_short(kb, pool, prods, k, tgt, n_tgt);
        for (uint32_t i = k; i > 0u; --i) {
            qua_and_uncompute(kb, ctrl, srcs[i - 1u], prods[i - 1u], pool->next_bit++);
        }
        qua_qpool_free_n(pool, prods, k);
    } else {
        qua_inv_ctrl_iadd_zero_prods(kb, pool, ctrl, srcs, k, tgt, n_tgt);
    }
    if (subtract) {
        for (uint32_t i = 0u; i < n_tgt; ++i) qua_x(kb, tgt[i]);
    }
}

static inline void qua_inv_barrel_denorm(qua_builder_t *kb, bool inverse,
                                         const uint32_t *clz, const uint32_t *y, uint32_t w0) {
    static const uint32_t idx[4] = {1u, 2u, 3u, 4u};
    static const uint32_t sh[4]  = {1u, 2u, 4u, 8u};
    const uint32_t top_live[4]   = {w0 - 1u, w0, w0 + 2u, w0 + 1u};
    if (!inverse) {
        for (uint32_t s = 0u; s < 4u; ++s) {
            uint32_t ctrl = clz[idx[s]];
            for (int32_t i = (int32_t)top_live[s]; i >= 0; --i) {
                qua_cswap(kb, ctrl, y[i], y[(uint32_t)i + sh[s]]);
            }
        }
    } else {
        for (int32_t s = 3; s >= 0; --s) {
            uint32_t ctrl = clz[idx[s]];
            for (uint32_t i = 0u; i <= top_live[s]; ++i) {
                qua_cswap(kb, ctrl, y[i], y[i + sh[s]]);
            }
        }
    }
}

static inline void qua_inv_sq_arith(qua_builder_t *kb, qua_qpool_t *pool,
                                    const uint32_t *m, uint32_t k_w,
                                    const uint32_t *acc, uint32_t wsq,
                                    int32_t shift, bool inverse) {
    int32_t base_pos = shift - 1;
    uint32_t W = wsq + 1u;
    if (!inverse) {
        qua_x(kb, acc[0]);
        for (uint32_t i = 0u; i < k_w; ++i) {
            int32_t p1 = (int32_t)i - 1;
            if (p1 >= base_pos) {
                uint32_t b = (uint32_t)(p1 - base_pos);
                if (b < W) qua_inv_gidney_inc(kb, pool, m[i], acc + b, W - b, false);
            }
            int32_t p2 = 2 * (int32_t)i - 1;
            if (p2 >= base_pos) {
                uint32_t b = (uint32_t)(p2 - base_pos);
                if (b < W) qua_inv_gidney_inc(kb, pool, m[i], acc + b, W - b, false);
            }
        }
        for (uint32_t i = 0u; i < k_w; ++i) {
            for (uint32_t j = i + 1u; j < k_w; ++j) {
                int32_t pos = (int32_t)(i + j);
                if (pos >= base_pos) {
                    uint32_t b = (uint32_t)(pos - base_pos);
                    if (b == W - 1u) {
                        qua_ccx(kb, m[i], m[j], acc[b]);
                    } else if (b < W - 1u) {
                        uint32_t p = qua_qpool_alloc1(pool);
                        qua_and_compute(kb, m[i], m[j], p);
                        qua_inv_gidney_inc(kb, pool, p, acc + b, W - b, false);
                        qua_and_uncompute(kb, m[i], m[j], p, pool->next_bit++);
                        qua_qpool_free1(pool, p);
                    }
                }
            }
        }
    } else {
        for (int32_t ii = (int32_t)k_w - 1; ii >= 0; --ii) {
            uint32_t i = (uint32_t)ii;
            for (int32_t jj = (int32_t)k_w - 1; jj > ii; --jj) {
                uint32_t j = (uint32_t)jj;
                int32_t pos = (int32_t)(i + j);
                if (pos >= base_pos) {
                    uint32_t b = (uint32_t)(pos - base_pos);
                    if (b == W - 1u) {
                        qua_ccx(kb, m[i], m[j], acc[b]);
                    } else if (b < W - 1u) {
                        uint32_t p = qua_qpool_alloc1(pool);
                        qua_and_compute(kb, m[i], m[j], p);
                        qua_inv_gidney_inc(kb, pool, p, acc + b, W - b, true);
                        qua_and_uncompute(kb, m[i], m[j], p, pool->next_bit++);
                        qua_qpool_free1(pool, p);
                    }
                }
            }
        }
        for (int32_t ii = (int32_t)k_w - 1; ii >= 0; --ii) {
            uint32_t i = (uint32_t)ii;
            int32_t p2 = 2 * (int32_t)i - 1;
            if (p2 >= base_pos) {
                uint32_t b = (uint32_t)(p2 - base_pos);
                if (b < W) qua_inv_gidney_inc(kb, pool, m[i], acc + b, W - b, true);
            }
            int32_t p1 = (int32_t)i - 1;
            if (p1 >= base_pos) {
                uint32_t b = (uint32_t)(p1 - base_pos);
                if (b < W) qua_inv_gidney_inc(kb, pool, m[i], acc + b, W - b, true);
            }
        }
        qua_x(kb, acc[0]);
    }
}

static inline int qua_invsqrt_emit_kmb_phase_split(const qua_invsqrt_cfg_t *cfg, const char *out_path) {
    qua_builder_t kb;
    qua_builder_init(&kb);

    uint32_t x[22], clz[5];
    uint32_t cursor = 0u;
    for (uint32_t i = 0u; i < 22u; ++i) x[i] = cursor++;
    qua_add_quantum_reg(&kb, "x", x, 22u);
    for (uint32_t i = 0u; i < 5u; ++i) clz[i] = cursor++;

    qua_qpool_t pool;
    qua_qpool_init(&pool, cursor);

    /* 1. Forward CLZ normalize (LinF2/CSWAP) */
    qua_inv_log_clz(&kb, &pool, false, x, 22u, clz);

    uint64_t c0_tab[QUA_INVSQRT_MAX_CHARTS];
    uint64_t c1_tab[QUA_INVSQRT_MAX_CHARTS];
    uint64_t c2_tab[QUA_INVSQRT_MAX_CHARTS];
    uint64_t c12_tab[QUA_INVSQRT_MAX_CHARTS];
    for (uint32_t c = 0u; c < cfg->n_charts; ++c) {
        c0_tab[c]  = (uint64_t)cfg->c0_biased[c];
        c1_tab[c]  = (uint64_t)cfg->c1[c];
        c2_tab[c]  = (uint64_t)cfg->c2[c];
        c12_tab[c] = (uint64_t)cfg->c1[c] | ((uint64_t)cfg->c2[c] << cfg->w1);
    }

    uint32_t addr[16];
    uint32_t n_addr = cfg->k_m + 1u;
    for (uint32_t k = 0u; k < cfg->k_m; ++k) addr[k] = x[21u - cfg->k_m + k];
    addr[cfg->k_m] = clz[0];

    /* 2. Stage A (`PhaseSplit` + `HMR_Dual`): c0_biased -> barrel -> Z_POW -> un-barrel -> HMR_Dual */
    uint32_t y_c0[40], y_ext[10], y_full[50];
    qua_qpool_alloc_n(&pool, y_c0, cfg->w0);
    {
        uint32_t cap_c0 = qua_qpool_avail(&pool, cfg->peak_qubits);
        qua_bitop_qroam_res_t c0_res;
        uint32_t bsh_c0 = qua_invsqrt_pick_sub_bsh(c0_tab, cfg->n_charts, n_addr, cfg->w0, cap_c0, &c0_res);
        qua_inv_emit_qroam_u64(&kb, &pool, c0_tab, cfg->n_charts, addr, n_addr, y_c0, cfg->w0, bsh_c0, cap_c0);
    }

    qua_qpool_alloc_n(&pool, y_ext, 10u);
    for (uint32_t i = 0u; i < cfg->w0; ++i) y_full[i] = y_c0[i];
    for (uint32_t i = 0u; i < 10u; ++i) y_full[cfg->w0 + i] = y_ext[i];
    qua_inv_barrel_denorm(&kb, false, clz, y_full, cfg->w0);
    for (uint32_t j = 0u; j < cfg->f_out; ++j) {
        int32_t exp2 = (int32_t)j - (int32_t)cfg->f_out;
        int32_t bit_pos = 128 + exp2;
        if (bit_pos >= 0 && bit_pos < 128) {
            flow_u128 turns = (flow_u128)1 << (uint32_t)bit_pos;
            qua_z_pow_turns_u128(&kb, y_full[cfg->drop + j], turns);
        }
    }
    qua_inv_barrel_denorm(&kb, true, clz, y_full, cfg->w0);
    qua_qpool_free_n(&pool, y_ext, 10u);

    qua_inv_qrom_del(&kb, &pool, c0_tab, cfg->n_charts, addr, n_addr,
                     y_c0, cfg->w0, qua_qpool_avail(&pool, cfg->peak_qubits));
    qua_qpool_free_n(&pool, y_c0, cfg->w0);

    /* 3. Stage B (`res = c2*sq + c1*t_inv` on exact w_res wires, 0 overflow qubits!) */
    for (uint32_t i = 0u; i < cfg->wd1; ++i) qua_x(&kb, x[i]);
    uint32_t sign_inv_wire = x[cfg->wd1 - 1u];
    const uint32_t *sq_addr_wires = x + (cfg->wd1 - cfg->wd2);
    uint32_t sq_addr_n = cfg->wd2 - 1u;

    uint64_t sq_tab64[QUA_INVSQRT_MAX_SQ];
    for (uint32_t a = 0u; a < cfg->sq_n; ++a) sq_tab64[a] = cfg->sq_tab[a];
    bool use_arith  = qua_invsqrt_use_sq_arith(cfg);
    int32_t arith_shift = use_arith ? ((int32_t)(2u * sq_addr_n - 1u) - (int32_t)cfg->wsq) : 0;
    bool allow_ps   = (cfg->zero_prods_mode == 0u);
    bool allow_cuc  = (cfg->cuccaro_mode != 0u);

    uint32_t y_res[40], c12_buf[64], c1_buf[32], c2_buf[32];
    qua_qpool_alloc_n(&pool, y_res, cfg->w_res);

    uint32_t *c1 = c1_buf;
    uint32_t *c2 = c2_buf;
    if (!cfg->split_c2) {
        qua_qpool_alloc_n(&pool, c12_buf, cfg->w1 + cfg->w2);
        c1 = c12_buf;
        c2 = c12_buf + cfg->w1;
        uint32_t cap_c12 = qua_qpool_avail(&pool, cfg->peak_qubits);
        qua_bitop_qroam_res_t c12_res;
        uint32_t bsh_c12 = qua_invsqrt_pick_sub_bsh(c12_tab, cfg->n_charts, n_addr, cfg->w1 + cfg->w2, cap_c12, &c12_res);
        qua_inv_emit_qroam_u64(&kb, &pool, c12_tab, cfg->n_charts, addr, n_addr,
                               c12_buf, cfg->w1 + cfg->w2, bsh_c12, cap_c12);
    } else {
        qua_qpool_alloc_n(&pool, c2, cfg->w2);
        uint32_t cap_c2 = qua_qpool_avail(&pool, cfg->peak_qubits);
        qua_bitop_qroam_res_t c2_res;
        uint32_t bsh_c2 = qua_invsqrt_pick_sub_bsh(c2_tab, cfg->n_charts, n_addr, cfg->w2, cap_c2, &c2_res);
        qua_inv_emit_qroam_u64(&kb, &pool, c2_tab, cfg->n_charts, addr, n_addr, c2, cfg->w2, bsh_c2, cap_c2);
    }

    /* 3a. Forward c2 * sq */
    for (uint32_t i = 0u; i < sq_addr_n; ++i) qua_cx(&kb, sign_inv_wire, sq_addr_wires[i]);
    uint32_t sq_acc[17];
    uint32_t *sq_q = use_arith ? (sq_acc + 1) : sq_acc;
    if (use_arith) {
        qua_qpool_alloc_n(&pool, sq_acc, cfg->wsq + 1u);
        for (uint32_t i = 0u; i < sq_addr_n; ++i) qua_x(&kb, sq_addr_wires[i]);
        qua_inv_sq_arith(&kb, &pool, sq_addr_wires, sq_addr_n, sq_acc, cfg->wsq, arith_shift, false);
    } else {
        qua_qpool_alloc_n(&pool, sq_q, cfg->wsq);
        qua_inv_emit_qroam_u64(&kb, &pool, sq_tab64, cfg->sq_n, sq_addr_wires, sq_addr_n,
                               sq_q, cfg->wsq, cfg->b_sh_sq, qua_qpool_avail(&pool, cfg->peak_qubits));
    }
    for (uint32_t b = 0u; b < cfg->wsq; ++b) {
        if (cfg->sq_row_hi[b] == 0u) continue;
        int32_t rsh = cfg->sh2 - (int32_t)b;
        const uint32_t *src = (rsh >= 0) ? (c2 + (uint32_t)rsh) : c2;
        uint32_t n_src      = (rsh >= 0) ? (cfg->w2 - (uint32_t)rsh) : cfg->w2;
        uint32_t lsh        = (rsh >= 0) ? 0u : (uint32_t)(-rsh);
        uint32_t n_t        = cfg->sq_row_hi[b] - lsh;
        qua_inv_hybrid_ctrl_add(&kb, &pool, sq_q[b], src, n_src, y_res + lsh, n_t, false,
                                qua_qpool_avail(&pool, cfg->peak_qubits), allow_ps, allow_cuc);
    }
    if (use_arith) {
        qua_inv_sq_arith(&kb, &pool, sq_addr_wires, sq_addr_n, sq_acc, cfg->wsq, arith_shift, true);
        for (uint32_t i = 0u; i < sq_addr_n; ++i) qua_x(&kb, sq_addr_wires[i]);
        qua_qpool_free_n(&pool, sq_acc, cfg->wsq + 1u);
    } else {
        qua_inv_qrom_del(&kb, &pool, sq_tab64, cfg->sq_n, sq_addr_wires, sq_addr_n,
                         sq_q, cfg->wsq, qua_qpool_avail(&pool, cfg->peak_qubits));
        qua_qpool_free_n(&pool, sq_q, cfg->wsq);
    }
    for (uint32_t i = 0u; i < sq_addr_n; ++i) qua_cx(&kb, sign_inv_wire, sq_addr_wires[i]);

    if (cfg->split_c2) {
        qua_inv_qrom_del(&kb, &pool, c2_tab, cfg->n_charts, addr, n_addr,
                         c2, cfg->w2, qua_qpool_avail(&pool, cfg->peak_qubits));
        qua_qpool_free_n(&pool, c2, cfg->w2);
        qua_qpool_alloc_n(&pool, c1, cfg->w1);
        uint32_t cap_c1 = qua_qpool_avail(&pool, cfg->peak_qubits);
        qua_bitop_qroam_res_t c1_res;
        uint32_t bsh_c1 = qua_invsqrt_pick_sub_bsh(c1_tab, cfg->n_charts, n_addr, cfg->w1, cap_c1, &c1_res);
        qua_inv_emit_qroam_u64(&kb, &pool, c1_tab, cfg->n_charts, addr, n_addr, c1, cfg->w1, bsh_c1, cap_c1);
    }

    /* 3b. Forward c1 * t_inv */
    for (uint32_t b = 0u; b < cfg->wd1; ++b) {
        if (cfg->c1_row_hi[b] == 0u) continue;
        int32_t rsh = cfg->sh1 - (int32_t)b;
        const uint32_t *src = (rsh >= 0) ? (c1 + (uint32_t)rsh) : c1;
        uint32_t n_src      = (rsh >= 0) ? (cfg->w1 - (uint32_t)rsh) : cfg->w1;
        uint32_t lsh        = (rsh >= 0) ? 0u : (uint32_t)(-rsh);
        uint32_t n_t        = cfg->c1_row_hi[b] - lsh;
        qua_inv_hybrid_ctrl_add(&kb, &pool, x[b], src, n_src, y_res + lsh, n_t, false,
                                qua_qpool_avail(&pool, cfg->peak_qubits), allow_ps, allow_cuc);
    }

    if (cfg->split_c2) {
        qua_inv_qrom_del(&kb, &pool, c1_tab, cfg->n_charts, addr, n_addr,
                         c1, cfg->w1, qua_qpool_avail(&pool, cfg->peak_qubits));
        qua_qpool_free_n(&pool, c1, cfg->w1);
    }

    /* 4. Stage B barrel denormalize(w_res) + Z_POW kickback */
    qua_qpool_alloc_n(&pool, y_ext, 10u);
    for (uint32_t i = 0u; i < cfg->w_res; ++i) y_full[i] = y_res[i];
    for (uint32_t i = 0u; i < 10u; ++i) y_full[cfg->w_res + i] = y_ext[i];

    qua_inv_barrel_denorm(&kb, false, clz, y_full, cfg->w_res);
    uint32_t f_res = (cfg->w_res + 10u > cfg->drop) ? (cfg->w_res + 10u - cfg->drop) : 0u;
    for (uint32_t j = 0u; j < f_res && j < cfg->f_out; ++j) {
        int32_t exp2 = (int32_t)j - (int32_t)cfg->f_out;
        int32_t bit_pos = 128 + exp2;
        if (bit_pos >= 0 && bit_pos < 128) {
            flow_u128 turns = (flow_u128)1 << (uint32_t)bit_pos;
            qua_z_pow_turns_u128(&kb, y_full[cfg->drop + j], turns);
        }
    }
    qua_inv_barrel_denorm(&kb, true, clz, y_full, cfg->w_res);
    qua_qpool_free_n(&pool, y_ext, 10u);

    /* 5. Exact reverse uncomputation of Stage B */
    if (cfg->split_c2) {
        qua_qpool_alloc_n(&pool, c1, cfg->w1);
        uint32_t cap_c1 = qua_qpool_avail(&pool, cfg->peak_qubits);
        qua_bitop_qroam_res_t c1_res;
        uint32_t bsh_c1 = qua_invsqrt_pick_sub_bsh(c1_tab, cfg->n_charts, n_addr, cfg->w1, cap_c1, &c1_res);
        qua_inv_emit_qroam_u64(&kb, &pool, c1_tab, cfg->n_charts, addr, n_addr, c1, cfg->w1, bsh_c1, cap_c1);
    }

    for (int32_t b = (int32_t)cfg->wd1 - 1; b >= 0; --b) {
        if (cfg->c1_row_hi[b] == 0u) continue;
        int32_t rsh = cfg->sh1 - b;
        const uint32_t *src = (rsh >= 0) ? (c1 + (uint32_t)rsh) : c1;
        uint32_t n_src      = (rsh >= 0) ? (cfg->w1 - (uint32_t)rsh) : cfg->w1;
        uint32_t lsh        = (rsh >= 0) ? 0u : (uint32_t)(-rsh);
        uint32_t n_t        = cfg->c1_row_hi[b] - lsh;
        qua_inv_hybrid_ctrl_add(&kb, &pool, x[b], src, n_src, y_res + lsh, n_t, true,
                                qua_qpool_avail(&pool, cfg->peak_qubits), allow_ps, allow_cuc);
    }

    if (cfg->split_c2) {
        qua_inv_qrom_del(&kb, &pool, c1_tab, cfg->n_charts, addr, n_addr,
                         c1, cfg->w1, qua_qpool_avail(&pool, cfg->peak_qubits));
        qua_qpool_free_n(&pool, c1, cfg->w1);
        qua_qpool_alloc_n(&pool, c2, cfg->w2);
        uint32_t cap_c2 = qua_qpool_avail(&pool, cfg->peak_qubits);
        qua_bitop_qroam_res_t c2_res;
        uint32_t bsh_c2 = qua_invsqrt_pick_sub_bsh(c2_tab, cfg->n_charts, n_addr, cfg->w2, cap_c2, &c2_res);
        qua_inv_emit_qroam_u64(&kb, &pool, c2_tab, cfg->n_charts, addr, n_addr, c2, cfg->w2, bsh_c2, cap_c2);
    }

    for (uint32_t i = 0u; i < sq_addr_n; ++i) qua_cx(&kb, sign_inv_wire, sq_addr_wires[i]);
    uint32_t sq_r_acc[17];
    uint32_t *sq_r = use_arith ? (sq_r_acc + 1) : sq_r_acc;
    if (use_arith) {
        qua_qpool_alloc_n(&pool, sq_r_acc, cfg->wsq + 1u);
        for (uint32_t i = 0u; i < sq_addr_n; ++i) qua_x(&kb, sq_addr_wires[i]);
        qua_inv_sq_arith(&kb, &pool, sq_addr_wires, sq_addr_n, sq_r_acc, cfg->wsq, arith_shift, false);
    } else {
        qua_qpool_alloc_n(&pool, sq_r, cfg->wsq);
        qua_inv_emit_qroam_u64(&kb, &pool, sq_tab64, cfg->sq_n, sq_addr_wires, sq_addr_n,
                               sq_r, cfg->wsq, cfg->b_sh_sq, qua_qpool_avail(&pool, cfg->peak_qubits));
    }
    for (int32_t b = (int32_t)cfg->wsq - 1; b >= 0; --b) {
        if (cfg->sq_row_hi[b] == 0u) continue;
        int32_t rsh = cfg->sh2 - b;
        const uint32_t *src = (rsh >= 0) ? (c2 + (uint32_t)rsh) : c2;
        uint32_t n_src      = (rsh >= 0) ? (cfg->w2 - (uint32_t)rsh) : cfg->w2;
        uint32_t lsh        = (rsh >= 0) ? 0u : (uint32_t)(-rsh);
        uint32_t n_t        = cfg->sq_row_hi[b] - lsh;
        qua_inv_hybrid_ctrl_add(&kb, &pool, sq_r[b], src, n_src, y_res + lsh, n_t, true,
                                qua_qpool_avail(&pool, cfg->peak_qubits), allow_ps, allow_cuc);
    }
    if (use_arith) {
        qua_inv_sq_arith(&kb, &pool, sq_addr_wires, sq_addr_n, sq_r_acc, cfg->wsq, arith_shift, true);
        for (uint32_t i = 0u; i < sq_addr_n; ++i) qua_x(&kb, sq_addr_wires[i]);
        qua_qpool_free_n(&pool, sq_r_acc, cfg->wsq + 1u);
    } else {
        qua_inv_qrom_del(&kb, &pool, sq_tab64, cfg->sq_n, sq_addr_wires, sq_addr_n,
                         sq_r, cfg->wsq, qua_qpool_avail(&pool, cfg->peak_qubits));
        qua_qpool_free_n(&pool, sq_r, cfg->wsq);
    }
    for (uint32_t i = 0u; i < sq_addr_n; ++i) qua_cx(&kb, sign_inv_wire, sq_addr_wires[i]);

    if (cfg->split_c2) {
        qua_inv_qrom_del(&kb, &pool, c2_tab, cfg->n_charts, addr, n_addr,
                         c2, cfg->w2, qua_qpool_avail(&pool, cfg->peak_qubits));
        qua_qpool_free_n(&pool, c2, cfg->w2);
    } else {
        qua_inv_qrom_del(&kb, &pool, c12_tab, cfg->n_charts, addr, n_addr,
                         c12_buf, cfg->w1 + cfg->w2, qua_qpool_avail(&pool, cfg->peak_qubits));
        qua_qpool_free_n(&pool, c12_buf, cfg->w1 + cfg->w2);
    }
    qua_qpool_free_n(&pool, y_res, cfg->w_res);

    for (uint32_t i = 0u; i < cfg->wd1; ++i) qua_x(&kb, x[i]);
    qua_inv_log_clz(&kb, &pool, true, x, 22u, clz);

    qua_bitop_spec_t fin_spec = {
        "invsqrt", QUA_COST_QC, true, 0u, 35u, 100.0, 20u, 2560u, 4u, NULL, NULL, NULL
    };
    int rc = qua_bitop_finalize_kmb(&kb, &fin_spec, out_path);
    qua_builder_free(&kb);
    return rc;
}

static inline int qua_invsqrt_emit_kmb(const qua_invsqrt_cfg_t *cfg, const char *out_path) {
    if (cfg->phase_split_mode) {
        return qua_invsqrt_emit_kmb_phase_split(cfg, out_path);
    }

    qua_builder_t kb;
    qua_builder_init(&kb);

    uint32_t x[22], clz[5], y_c0[40], c1[32], c2_buf[32], ovf1, ovf2, ovf3;
    uint32_t cursor = 0u;
    for (uint32_t i = 0u; i < 22u; ++i) x[i] = cursor++;
    qua_add_quantum_reg(&kb, "x", x, 22u);
    for (uint32_t i = 0u; i < 5u; ++i) clz[i] = cursor++;
    for (uint32_t i = 0u; i < cfg->w0; ++i) y_c0[i] = cursor++;
    for (uint32_t i = 0u; i < cfg->w1; ++i) c1[i] = cursor++;
    uint32_t *c2 = cfg->split_c2 ? c1 : c2_buf;
    if (!cfg->split_c2) {
        for (uint32_t i = 0u; i < cfg->w2; ++i) c2_buf[i] = cursor++;
    } else if (cfg->w2 > cfg->w1) {
        for (uint32_t i = cfg->w1; i < cfg->w2; ++i) c2[i] = cursor++;
    }
    ovf1 = cursor++;
    ovf2 = cursor++;
    ovf3 = cursor++;

    uint32_t pool_free_cap = (cfg->peak_qubits > cursor) ? (cfg->peak_qubits - cursor) : (cfg->m3 + 1u);

    qua_qpool_t pool;
    qua_qpool_init(&pool, cursor);

    /* 1. Forward CLZ normalize */
    qua_inv_log_clz(&kb, &pool, false, x, 22u, clz);

    /* 2. Packed 2-Jet QROAM (Virtual-Slot-0 + Constant-Column Pruning) */
    uint64_t packed_tab[QUA_INVSQRT_MAX_CHARTS];
    uint64_t c0_c2_tab[QUA_INVSQRT_MAX_CHARTS];
    uint64_t c1_tab[QUA_INVSQRT_MAX_CHARTS];
    uint64_t c2_tab[QUA_INVSQRT_MAX_CHARTS];
    uint32_t addr[16];
    uint32_t n_addr = cfg->k_m + 1u;
    for (uint32_t k = 0u; k < cfg->k_m; ++k) addr[k] = x[21u - cfg->k_m + k];
    addr[cfg->k_m] = clz[0];

    uint32_t qrom_tgt[104];
    uint32_t c0_c2_tgt[80];
    uint32_t n_qtgt = 0u;
    for (uint32_t i = 0u; i < cfg->w0; ++i) qrom_tgt[n_qtgt++] = y_c0[i];
    for (uint32_t i = 0u; i < cfg->w1; ++i) qrom_tgt[n_qtgt++] = c1[i];
    if (!cfg->split_c2) {
        for (uint32_t i = 0u; i < cfg->w2; ++i) qrom_tgt[n_qtgt++] = c2[i];
        for (uint32_t c = 0u; c < cfg->n_charts; ++c) {
            packed_tab[c] = (uint64_t)cfg->c0_biased[c] |
                            ((uint64_t)cfg->c1[c] << cfg->w0) |
                            ((uint64_t)cfg->c2[c] << (cfg->w0 + cfg->w1));
        }
        qua_bitop_qroam_res_t m_res;
        uint32_t bsh_m = qua_invsqrt_pick_sub_bsh(packed_tab, cfg->n_charts, n_addr, n_qtgt, pool_free_cap, &m_res);
        qua_inv_emit_qroam_u64(&kb, &pool, packed_tab, cfg->n_charts, addr, n_addr,
                               qrom_tgt, n_qtgt, bsh_m, pool_free_cap);
    } else {
        uint32_t n_c0c2 = 0u;
        for (uint32_t i = 0u; i < cfg->w0; ++i) c0_c2_tgt[n_c0c2++] = y_c0[i];
        for (uint32_t i = 0u; i < cfg->w2; ++i) c0_c2_tgt[n_c0c2++] = c2[i];
        for (uint32_t c = 0u; c < cfg->n_charts; ++c) {
            c0_c2_tab[c] = (uint64_t)cfg->c0_biased[c] |
                           ((uint64_t)cfg->c2[c] << cfg->w0);
            c1_tab[c] = (uint64_t)cfg->c1[c];
            c2_tab[c] = (uint64_t)cfg->c2[c];
        }
        qua_bitop_qroam_res_t m_res;
        uint32_t bsh_m = qua_invsqrt_pick_sub_bsh(c0_c2_tab, cfg->n_charts, n_addr, n_c0c2, pool_free_cap, &m_res);
        qua_inv_emit_qroam_u64(&kb, &pool, c0_c2_tab, cfg->n_charts, addr, n_addr,
                               c0_c2_tgt, n_c0c2, bsh_m, pool_free_cap);
    }

    /* 3. Unsigned bias flip on x[0:wd1] */
    for (uint32_t i = 0u; i < cfg->wd1; ++i) qua_x(&kb, x[i]);
    uint32_t sign_inv_wire = x[cfg->wd1 - 1u];
    const uint32_t *sq_addr_wires = x + (cfg->wd1 - cfg->wd2);
    uint32_t sq_addr_n = cfg->wd2 - 1u;

    uint64_t sq_tab64[QUA_INVSQRT_MAX_SQ];
    for (uint32_t a = 0u; a < cfg->sq_n; ++a) sq_tab64[a] = cfg->sq_tab[a];

    bool use_arith = qua_invsqrt_use_sq_arith(cfg);
    int32_t arith_shift = use_arith ? ((int32_t)(2u * sq_addr_n - 1u) - (int32_t)cfg->wsq) : 0;
    uint32_t sq_free_cap = (pool_free_cap > cfg->wsq) ? (pool_free_cap - cfg->wsq) : 0u;
    bool allow_ps  = (cfg->zero_prods_mode == 0u);
    bool allow_cuc = (cfg->cuccaro_mode != 0u);

    /* 3a. Tier 1a: c2 * sq */
    for (uint32_t i = 0u; i < sq_addr_n; ++i) qua_cx(&kb, sign_inv_wire, sq_addr_wires[i]);
    uint32_t sq_acc[17];
    uint32_t *sq_q = use_arith ? (sq_acc + 1) : sq_acc;
    if (use_arith) {
        qua_qpool_alloc_n(&pool, sq_acc, cfg->wsq + 1u);
        for (uint32_t i = 0u; i < sq_addr_n; ++i) qua_x(&kb, sq_addr_wires[i]);
        qua_inv_sq_arith(&kb, &pool, sq_addr_wires, sq_addr_n, sq_acc, cfg->wsq, arith_shift, false);
    } else {
        qua_qpool_alloc_n(&pool, sq_q, cfg->wsq);
        qua_inv_emit_qroam_u64(&kb, &pool, sq_tab64, cfg->sq_n, sq_addr_wires, sq_addr_n,
                               sq_q, cfg->wsq, cfg->b_sh_sq, sq_free_cap);
    }

    for (uint32_t b = 0u; b < cfg->wsq; ++b) {
        int32_t rsh = cfg->sh2 - (int32_t)b;
        if (rsh >= (int32_t)cfg->w2) continue;
        const uint32_t *src = (rsh >= 0) ? (c2 + (uint32_t)rsh) : c2;
        uint32_t n_src      = (rsh >= 0) ? (cfg->w2 - (uint32_t)rsh) : cfg->w2;
        uint32_t lsh        = (rsh >= 0) ? 0u : (uint32_t)(-rsh);
        uint32_t tgt_slice[40];
        uint32_t n_t = 0u;
        for (uint32_t i = lsh; i < cfg->m1; ++i) tgt_slice[n_t++] = y_c0[i];
        tgt_slice[n_t++] = ovf1;
        qua_inv_hybrid_ctrl_add(&kb, &pool, sq_q[b], src, n_src, tgt_slice, n_t, false,
                                qua_qpool_avail(&pool, cfg->peak_qubits), allow_ps, allow_cuc);
    }
    if (use_arith) {
        qua_inv_sq_arith(&kb, &pool, sq_addr_wires, sq_addr_n, sq_acc, cfg->wsq, arith_shift, true);
        for (uint32_t i = 0u; i < sq_addr_n; ++i) qua_x(&kb, sq_addr_wires[i]);
        qua_qpool_free_n(&pool, sq_acc, cfg->wsq + 1u);
    } else {
        qua_inv_qrom_del(&kb, &pool, sq_tab64, cfg->sq_n, sq_addr_wires, sq_addr_n, sq_q, cfg->wsq, sq_free_cap);
        qua_qpool_free_n(&pool, sq_q, cfg->wsq);
    }
    for (uint32_t i = 0u; i < sq_addr_n; ++i) qua_cx(&kb, sign_inv_wire, sq_addr_wires[i]);
    if (cfg->split_c2) {
        qua_inv_qrom_del(&kb, &pool, c2_tab, cfg->n_charts, addr, n_addr, c2, cfg->w2, pool_free_cap);
        qua_bitop_qroam_res_t c1_res;
        uint32_t bsh_c1 = qua_invsqrt_pick_sub_bsh(c1_tab, cfg->n_charts, n_addr, cfg->w1, pool_free_cap, &c1_res);
        qua_inv_emit_qroam_u64(&kb, &pool, c1_tab, cfg->n_charts, addr, n_addr, c1, cfg->w1, bsh_c1, pool_free_cap);
    }

    /* 3b..3d. Lazy-Carry Tiers 1b, 2, 3 */
    const uint32_t *tier_arr[3] = {cfg->tier1_b, cfg->tier2_b, cfg->tier3_b};
    const uint32_t tier_n[3]    = {cfg->tier1_n, cfg->tier2_n, cfg->tier3_n};
    const uint32_t tier_m[3]    = {cfg->m1, cfg->m2, cfg->m3};
    const uint32_t tier_ovf[3]  = {ovf1, ovf2, ovf3};

    for (uint32_t ti = 0u; ti < 3u; ++ti) {
        for (uint32_t j = 0u; j < tier_n[ti]; ++j) {
            uint32_t b = tier_arr[ti][j];
            int32_t rsh = cfg->sh1 - (int32_t)b;
            if (rsh >= (int32_t)cfg->w1) continue;
            const uint32_t *src = (rsh >= 0) ? (c1 + (uint32_t)rsh) : c1;
            uint32_t n_src      = (rsh >= 0) ? (cfg->w1 - (uint32_t)rsh) : cfg->w1;
            uint32_t lsh        = (rsh >= 0) ? 0u : (uint32_t)(-rsh);
            uint32_t tgt_slice[40];
            uint32_t n_t = 0u;
            for (uint32_t i = lsh; i < tier_m[ti]; ++i) tgt_slice[n_t++] = y_c0[i];
            tgt_slice[n_t++] = tier_ovf[ti];
            qua_inv_hybrid_ctrl_add(&kb, &pool, x[b], src, n_src, tgt_slice, n_t, false,
                                    pool_free_cap, allow_ps, allow_cuc);
        }
        uint32_t inc_slice[40];
        uint32_t n_inc = 0u;
        uint32_t hi_end = (ti < 2u) ? tier_m[ti + 1u] : cfg->w0;
        for (uint32_t i = tier_m[ti]; i < hi_end; ++i) inc_slice[n_inc++] = y_c0[i];
        if (ti < 2u) inc_slice[n_inc++] = tier_ovf[ti + 1u];
        qua_inv_gidney_inc(&kb, &pool, tier_ovf[ti], inc_slice, n_inc, false);
    }

    /* 4. Late-allocated y_ext(10) + barrel denormalize + Z_POW phase kickback */
    uint32_t y_ext[10], y_full[50];
    qua_qpool_alloc_n(&pool, y_ext, 10u);
    for (uint32_t i = 0u; i < cfg->w0; ++i) y_full[i] = y_c0[i];
    for (uint32_t i = 0u; i < 10u; ++i) y_full[cfg->w0 + i] = y_ext[i];

    qua_inv_barrel_denorm(&kb, false, clz, y_full, cfg->w0);
    for (uint32_t j = 0u; j < cfg->f_out; ++j) {
        int32_t exp2 = (int32_t)j - (int32_t)cfg->f_out;
        /* turns = 2^exp2 mod 1 in Q128 */
        int32_t bit_pos = 128 + exp2;
        if (bit_pos >= 0 && bit_pos < 128) {
            flow_u128 turns = (flow_u128)1 << (uint32_t)bit_pos;
            qua_z_pow_turns_u128(&kb, y_full[cfg->drop + j], turns);
        }
    }
    qua_inv_barrel_denorm(&kb, true, clz, y_full, cfg->w0);
    qua_qpool_free_n(&pool, y_ext, 10u);

    /* 5. Exact reverse uncomputation of Lazy-Carry tiers */
    for (int32_t ti = 2; ti >= 0; --ti) {
        uint32_t inc_slice[40];
        uint32_t n_inc = 0u;
        uint32_t hi_end = (ti < 2) ? tier_m[ti + 1] : cfg->w0;
        for (uint32_t i = tier_m[ti]; i < hi_end; ++i) inc_slice[n_inc++] = y_c0[i];
        if (ti < 2) inc_slice[n_inc++] = tier_ovf[ti + 1];
        qua_inv_gidney_inc(&kb, &pool, tier_ovf[ti], inc_slice, n_inc, true);

        for (int32_t j = (int32_t)tier_n[ti] - 1; j >= 0; --j) {
            uint32_t b = tier_arr[ti][j];
            int32_t rsh = cfg->sh1 - (int32_t)b;
            if (rsh >= (int32_t)cfg->w1) continue;
            const uint32_t *src = (rsh >= 0) ? (c1 + (uint32_t)rsh) : c1;
            uint32_t n_src      = (rsh >= 0) ? (cfg->w1 - (uint32_t)rsh) : cfg->w1;
            uint32_t lsh        = (rsh >= 0) ? 0u : (uint32_t)(-rsh);
            uint32_t tgt_slice[40];
            uint32_t n_t = 0u;
            for (uint32_t i = lsh; i < tier_m[ti]; ++i) tgt_slice[n_t++] = y_c0[i];
            tgt_slice[n_t++] = tier_ovf[ti];
            qua_inv_hybrid_ctrl_add(&kb, &pool, x[b], src, n_src, tgt_slice, n_t, true,
                                    pool_free_cap, allow_ps, allow_cuc);
        }
    }

    if (cfg->split_c2) {
        qua_inv_qrom_del(&kb, &pool, c1_tab, cfg->n_charts, addr, n_addr, c1, cfg->w1, pool_free_cap);
        qua_bitop_qroam_res_t c2_res;
        uint32_t bsh_c2 = qua_invsqrt_pick_sub_bsh(c2_tab, cfg->n_charts, n_addr, cfg->w2, pool_free_cap, &c2_res);
        qua_inv_emit_qroam_u64(&kb, &pool, c2_tab, cfg->n_charts, addr, n_addr, c2, cfg->w2, bsh_c2, pool_free_cap);
    }

    for (uint32_t i = 0u; i < sq_addr_n; ++i) qua_cx(&kb, sign_inv_wire, sq_addr_wires[i]);
    uint32_t sq_r_acc[17];
    uint32_t *sq_r = use_arith ? (sq_r_acc + 1) : sq_r_acc;
    if (use_arith) {
        qua_qpool_alloc_n(&pool, sq_r_acc, cfg->wsq + 1u);
        for (uint32_t i = 0u; i < sq_addr_n; ++i) qua_x(&kb, sq_addr_wires[i]);
        qua_inv_sq_arith(&kb, &pool, sq_addr_wires, sq_addr_n, sq_r_acc, cfg->wsq, arith_shift, false);
    } else {
        qua_qpool_alloc_n(&pool, sq_r, cfg->wsq);
        qua_inv_emit_qroam_u64(&kb, &pool, sq_tab64, cfg->sq_n, sq_addr_wires, sq_addr_n,
                               sq_r, cfg->wsq, cfg->b_sh_sq, sq_free_cap);
    }
    for (int32_t b = (int32_t)cfg->wsq - 1; b >= 0; --b) {
        int32_t rsh = cfg->sh2 - b;
        if (rsh >= (int32_t)cfg->w2) continue;
        const uint32_t *src = (rsh >= 0) ? (c2 + (uint32_t)rsh) : c2;
        uint32_t n_src      = (rsh >= 0) ? (cfg->w2 - (uint32_t)rsh) : cfg->w2;
        uint32_t lsh        = (rsh >= 0) ? 0u : (uint32_t)(-rsh);
        uint32_t tgt_slice[40];
        uint32_t n_t = 0u;
        for (uint32_t i = lsh; i < cfg->m1; ++i) tgt_slice[n_t++] = y_c0[i];
        tgt_slice[n_t++] = ovf1;
        qua_inv_hybrid_ctrl_add(&kb, &pool, sq_r[b], src, n_src, tgt_slice, n_t, true,
                                qua_qpool_avail(&pool, cfg->peak_qubits), allow_ps, allow_cuc);
    }
    if (use_arith) {
        qua_inv_sq_arith(&kb, &pool, sq_addr_wires, sq_addr_n, sq_r_acc, cfg->wsq, arith_shift, true);
        for (uint32_t i = 0u; i < sq_addr_n; ++i) qua_x(&kb, sq_addr_wires[i]);
        qua_qpool_free_n(&pool, sq_r_acc, cfg->wsq + 1u);
    } else {
        qua_inv_qrom_del(&kb, &pool, sq_tab64, cfg->sq_n, sq_addr_wires, sq_addr_n, sq_r, cfg->wsq, sq_free_cap);
        qua_qpool_free_n(&pool, sq_r, cfg->wsq);
    }
    for (uint32_t i = 0u; i < sq_addr_n; ++i) qua_cx(&kb, sign_inv_wire, sq_addr_wires[i]);

    for (uint32_t i = 0u; i < cfg->wd1; ++i) qua_x(&kb, x[i]);
    if (!cfg->split_c2) {
        qua_inv_qrom_del(&kb, &pool, packed_tab, cfg->n_charts, addr, n_addr, qrom_tgt, n_qtgt, pool_free_cap);
    } else {
        qua_inv_qrom_del(&kb, &pool, c0_c2_tab, cfg->n_charts, addr, n_addr, c0_c2_tgt, cfg->w0 + cfg->w2, pool_free_cap);
    }
    qua_inv_log_clz(&kb, &pool, true, x, 22u, clz);

    qua_bitop_spec_t fin_spec = {
        "invsqrt", QUA_COST_QC, true, 0u, 35u, 100.0, 20u, 2560u, 4u, NULL, NULL, NULL
    };
    int rc = qua_bitop_finalize_kmb(&kb, &fin_spec, out_path);
    qua_builder_free(&kb);
    return rc;
}

/* ============================================================================
 * Bilevel Meta-SEARU x FLOW 2-Jet Pure Search Driver for `invsqrt`
 * ============================================================================ */
static inline void qua_invsqrt_populate_dag(const qua_invsqrt_cfg_t *cfg, qua_bitop_dag_t *dag) {
    if (!dag) return;
    dag->base_io_qubits = 22u;
    dag->cost_law       = QUA_COST_QC;
    uint32_t n0 = qua_bitop_dag_push_rep(dag, QUA_BOP_CSWAP, QUA_REP_LINF2, 5u, 105u, 101u, 0u, true, false, 0u);
    if (cfg->phase_split_mode) {
        uint32_t n1 = qua_bitop_dag_push_rep(dag, QUA_BOP_QROM, QUA_REP_PARITY_DUAL, cfg->w0, cfg->n_charts, 28u, 0u,
                                             true, false, UINT64_C(1) << n0);
        uint32_t n2 = qua_bitop_dag_push_rep(dag, QUA_BOP_PHASE, QUA_REP_PHASE, 10u, 8u * cfg->w0 + 12u, 0u,
                                             cfg->f_out, false, true, (UINT64_C(1) << n0) | (UINT64_C(1) << n1));
        uint32_t n3 = qua_bitop_dag_push_rep(dag, QUA_BOP_MUL, QUA_REP_KET, cfg->w_res,
                                             (uint32_t)(cfg->est_toffolis > 300u ? (cfg->est_toffolis - 300u) / 2u : 150u),
                                             (uint32_t)(cfg->est_toffolis > 300u ? (cfg->est_toffolis - 300u) / 2u : 150u),
                                             0u, true, false, UINT64_C(1) << n0);
        uint32_t f_res = (cfg->w_res + 10u > cfg->drop) ? (cfg->w_res + 10u - cfg->drop) : 0u;
        (void)qua_bitop_dag_push_rep(dag, QUA_BOP_PHASE, QUA_REP_PHASE, 10u, 8u * cfg->w_res + 12u, 0u,
                                     f_res, false, true,
                                     (UINT64_C(1) << n0) | (UINT64_C(1) << n2) | (UINT64_C(1) << n3));
        qua_ele_pebble_plan_t pp;
        qua_bitop_eval_dag_pebble(dag, UINT64_C(1) << n1, &pp, NULL, 0u);
        return;
    }
    uint32_t main_w = cfg->split_c2 ? (cfg->w0 + cfg->w2) : (cfg->w0 + cfg->w1 + cfg->w2);
    uint32_t n1 = qua_bitop_dag_push_rep(dag, QUA_BOP_QROM, QUA_REP_PARITY_DUAL, main_w, cfg->n_charts, 28u, 0u,
                                         true, false, UINT64_C(1) << n0);
    uint32_t n2 = qua_bitop_dag_push(dag, qua_invsqrt_use_sq_arith(cfg) ? QUA_BOP_MUL : QUA_BOP_QROM,
                                     cfg->wsq, 48u, 24u, 0u, true, false, UINT64_C(1) << n0);
    uint32_t n3 = qua_bitop_dag_push(dag, QUA_BOP_MUL, 1u, 120u, 120u, 0u, true, false,
                                     (UINT64_C(1) << n1) | (UINT64_C(1) << n2));
    uint32_t n4 = n1;
    if (cfg->split_c2) {
        n4 = qua_bitop_dag_push_rep(dag, QUA_BOP_QROM, QUA_REP_PARITY_DUAL, cfg->w1, cfg->n_charts / 2u, 28u, 0u,
                                    true, false, (UINT64_C(1) << n0) | (UINT64_C(1) << n3));
    }
    uint32_t n5 = qua_bitop_dag_push(dag, QUA_BOP_AND, cfg->m3 + 3u,
                                     (uint32_t)(cfg->est_toffolis > 460u ? (cfg->est_toffolis - 460u) / 2u : 200u),
                                     (uint32_t)(cfg->est_toffolis > 460u ? (cfg->est_toffolis - 460u) / 2u : 200u),
                                     0u, true, false, (UINT64_C(1) << n3) | (UINT64_C(1) << n4));
    (void)qua_bitop_dag_push_rep(dag, QUA_BOP_PHASE, QUA_REP_PHASE, 10u, 8u * cfg->w0 + 12u, 0u, cfg->f_out,
                                 false, true, (UINT64_C(1) << n0) | (UINT64_C(1) << n5));
    qua_ele_pebble_plan_t pp;
    qua_bitop_eval_dag_pebble(dag, (UINT64_C(1) << n2) | (cfg->split_c2 ? (UINT64_C(1) << n1) : 0u), &pp, NULL, 0u);
}

static inline uint64_t qua_invsqrt_frust_mask(const qua_invsqrt_cfg_t *cfg, double thresh) {
    if (!cfg->tiers_valid) {
        /* Steer K_m..wsq [0..17], M1..M3 [21..27], phase_split_mode [34], cuccaro_mode [35] */
        return UINT64_C(0xC0FE3FFFF);
    }
    if (cfg->worst_err_ratio >= thresh || cfg->min_margin <= 0) {
        /* Steer K_m..f_out [0..20], sq_arith [31], phase_split_mode [34], cuccaro_mode [35] */
        return UINT64_C(0xC801FFFFF);
    }
    /* Feasible manifold: unlock all 36 architectural + shortcut + QROAM + arithmetic genes */
    return UINT64_C(0xFFFFFFFFF);
}

static inline bool qua_invsqrt_bitop_eval(
    const uint64_t *genome,
    qua_bitop_dag_t *dag,
    double *out_raw_cost,
    uint64_t *out_frust,
    void *user_ctx
) {
    const qua_invsqrt_ref_t *ref = (const qua_invsqrt_ref_t *)user_ctx;
    qua_invsqrt_cfg_t cfg;
    qua_invsqrt_decode_genome(genome[0], &cfg);
    if (!cfg.phase_split_mode && !cfg.split_c2 && (cfg.w0 + cfg.w1 + cfg.w2 > 64u)) {
        if (dag) {
            dag->q_violation   = (double)(cfg.w0 + cfg.w1 + cfg.w2 - 64u);
            dag->err_violation = 0.5;
        }
        *out_frust    = UINT64_C(0xC10000FF8);
        *out_raw_cost = 600000.0;
        return false;
    }

    qua_invsqrt_calibrate(&cfg, ref);
    qua_invsqrt_check_tiers_only(&cfg);
    if (!cfg.tiers_valid) {
        if (dag) {
            dag->q_violation   = (double)cfg.tier_overflow_bits;
            dag->err_violation = 0.25;
        }
        *out_frust    = qua_invsqrt_frust_mask(&cfg, 0.915);
        *out_raw_cost = 350000.0 + 10000.0 * (double)cfg.tier_overflow_bits;
        return false;
    }

    qua_invsqrt_center_charts_fast(&cfg, ref);
    if (cfg.phase_split_mode) {
        qua_invsqrt_check_tiers_only(&cfg);
        if (!cfg.tiers_valid) {
            if (dag) {
                dag->q_violation   = (double)cfg.tier_overflow_bits;
                dag->err_violation = 0.25;
            }
            *out_frust    = qua_invsqrt_frust_mask(&cfg, 0.915);
            *out_raw_cost = 350000.0 + 10000.0 * (double)cfg.tier_overflow_bits;
            return false;
        }
    }
    qua_invsqrt_eval_error(&cfg, ref, false);

    const double thresh = 0.915;
    double ev = 0.0;
    if (cfg.worst_err_ratio >= thresh) {
        ev = (cfg.worst_err_ratio - thresh) + 0.015;
    } else if (cfg.min_margin <= 0) {
        ev = 0.05;
    }
    if (ev > 0.35) {
        if (dag) {
            dag->q_violation   = 0.0;
            dag->err_violation = ev;
        }
        *out_frust    = qua_invsqrt_frust_mask(&cfg, thresh);
        *out_raw_cost = 240000.0 + 40000.0 * ev;
        return false;
    }

    qua_invsqrt_plan_tiers_and_cost(&cfg);
    cfg.est_score = (double)cfg.peak_qubits * ((double)cfg.est_toffolis + cfg.rot_c);
    qua_invsqrt_populate_dag(&cfg, dag);

    if (dag) {
        dag->q_violation   = cfg.tiers_valid ? 0.0 : (double)cfg.tier_overflow_bits;
        dag->err_violation = ev;
    }
    *out_frust    = qua_invsqrt_frust_mask(&cfg, thresh);
    *out_raw_cost = cfg.est_score;
    return cfg.tiers_valid && (ev <= 0.0);
}

static inline bool qua_invsqrt_bitop_relax_2jet(
    const uint64_t *genome,
    qua_bitop_dag_t *dag,
    double *inout_raw_cost,
    uint64_t *inout_frust,
    void *user_ctx
) {
    /* Skip redundant 2-Jet relaxation when candidate already cleared fast threshold! */
    if (dag != NULL && dag->q_violation <= 0.0 && dag->err_violation <= 0.0) {
        return true;
    }
    const qua_invsqrt_ref_t *ref = (const qua_invsqrt_ref_t *)user_ctx;
    qua_invsqrt_cfg_t cfg;
    qua_invsqrt_decode_genome(genome[0], &cfg);
    if (!cfg.phase_split_mode && !cfg.split_c2 && (cfg.w0 + cfg.w1 + cfg.w2 > 64u)) return false;

    qua_invsqrt_calibrate(&cfg, ref);
    qua_invsqrt_check_tiers_only(&cfg);
    if (!cfg.tiers_valid) return false;
    qua_invsqrt_center_charts_fast(&cfg, ref);
    qua_invsqrt_relax_charts_fast_2jet(&cfg, ref);
    if (cfg.phase_split_mode) qua_invsqrt_check_tiers_only(&cfg);
    qua_invsqrt_eval_error(&cfg, ref, false);
    qua_invsqrt_plan_tiers_and_cost(&cfg);

    const double thresh = 0.915;
    double qv = cfg.tiers_valid ? 0.0 : (double)cfg.tier_overflow_bits;
    double ev = 0.0;
    if (cfg.worst_err_ratio >= thresh) {
        ev = (cfg.worst_err_ratio - thresh) + 0.015;
    } else if (cfg.min_margin <= 0) {
        ev = 0.05;
    }
    if (dag) {
        dag->q_violation   = qv;
        dag->err_violation = ev;
    }
    inout_frust[0]  = qua_invsqrt_frust_mask(&cfg, thresh);
    *inout_raw_cost = cfg.est_score;
    return cfg.tiers_valid && (ev <= 0.0);
}

static inline void qua_invsqrt_run_seeded_inner(
    const qua_bitop_spec_t *spec,
    const searu_hyperparams *hp,
    uint64_t seed,
    uint32_t iters,
    uint64_t init_g,
    uint64_t pin_mask,
    uint64_t out_genome[4],
    double *out_best_raw
) {
    static const double BETA_SCHED[4] = {0.25, 1.0, 4.0, 32.0};
    const double escale = (spec->energy_scale > 0.0) ? spec->energy_scale : 1.0;
    const uint64_t full_mask = (spec->genome_bits >= 64u)
        ? UINT64_MAX
        : ((UINT64_C(1) << spec->genome_bits) - 1u);

    uint64_t cur_init[4] = {init_g & full_mask, 0u, 0u, 0u};
    uint64_t kan_pin[4]  = {pin_mask & full_mask, 0u, 0u, 0u};
    uint64_t best_feas_g[4] = {0u, 0u, 0u, 0u};
    double   best_feas_raw = 1e30;
    bool     have_feas = false;

    uint32_t n_stages = (iters >= 64u) ? 4u : 1u;
    uint32_t stage_iters = (iters + n_stages - 1u) / n_stages;
    if (stage_iters < 16u) stage_iters = 16u;

    for (uint32_t stg = 0u; stg < n_stages; ++stg) {
        const double beta = (n_stages == 1u) ? 32.0 : BETA_SCHED[stg];
        searu_hyperparams shp = *hp;
        if (stg > 0u) {
            shp.inl_init_temp_q12 = (hp->inl_init_temp_q12 >> stg);
            if (shp.inl_init_temp_q12 < (UINT64_C(16) << 12u)) {
                shp.inl_init_temp_q12 = (UINT64_C(16) << 12u);
            }
        }
        searu_result res;
        uint64_t stg_seed = seed ^ (UINT64_C(0xD1B54A32D192ED03) * (stg + 1u));
        SEARU_SOLVE_INLINE_HP(&res, spec->genome_bits, cur_init, kan_pin, stg_seed, stage_iters, &shp,
                              cand, frust, energy, {
            qua_bitop_dag_t dag;
            memset(&dag, 0, sizeof(dag));
            dag.cost_law = spec->cost_law;
            double raw_cost = 1e30;
            uint64_t f_words[4] = {full_mask, UINT64_MAX, UINT64_MAX, UINT64_MAX};
            bool ok = spec->eval_spec(cand, &dag, &raw_cost, f_words, spec->user_ctx);
            if (spec->relax_2jet != NULL && isfinite(raw_cost) && raw_cost > 0.0 &&
                dag.q_violation <= 2.0 && dag.err_violation <= 0.35) {
                bool ok2 = spec->relax_2jet(cand, &dag, &raw_cost, f_words, spec->user_ctx);
                if (ok2) ok = true;
            }
            f_words[0] &= ~kan_pin[0];
            if (f_words[0] == 0u) f_words[0] = full_mask & ~kan_pin[0];
            frust[0] = f_words[0];
            frust[1] = f_words[1];
            frust[2] = f_words[2];
            frust[3] = f_words[3];

            double qv = (dag.q_violation > 0.0) ? dag.q_violation : 0.0;
            double ev = (dag.err_violation > 0.0) ? dag.err_violation : 0.0;
            bool strict_feas = ok && isfinite(raw_cost) && (raw_cost > 0.0) &&
                               (qv <= 0.0) && (ev <= 0.0);
            if (strict_feas && raw_cost < best_feas_raw) {
                best_feas_raw = raw_cost;
                for (int k = 0; k < 4; ++k) best_feas_g[k] = cand[k];
                have_feas = true;
            }
            if (!isfinite(raw_cost) || raw_cost <= 0.0) {
                energy = UINT64_C(1000000000);
            } else {
                double base_e = raw_cost / escale;
                double ev_log = log2(1.0 + ev);
                double ev_lin = (ev < 2.0) ? ev : (2.0 + 0.25 * ev_log);
                double pen_e  = beta * (64.0 * qv + 8.0 * qv * qv +
                                        512.0 * ev_lin + 384.0 * ev_log * ev_log);
                double tot_e  = base_e + pen_e;
                if (tot_e < 1.0) tot_e = 1.0;
                if (tot_e > 1e11) tot_e = 1e11;
                energy = (uint64_t)llround(tot_e);
            }
        });
        for (int k = 0; k < 4; ++k) cur_init[k] = res.canvas[k];
    }

    if (have_feas) {
        for (int i = 0; i < 4; ++i) out_genome[i] = best_feas_g[i];
        if (out_best_raw) *out_best_raw = best_feas_raw;
    } else {
        for (int i = 0; i < 4; ++i) out_genome[i] = cur_init[i];
        if (out_best_raw) *out_best_raw = 1e30;
    }
}

static inline int qua_invsqrt_run_pure_search(const char *out_kmb_path) {
    qua_invsqrt_ref_t *ref = (qua_invsqrt_ref_t *)malloc(sizeof(qua_invsqrt_ref_t));
    if (!ref) return 1;
    qua_invsqrt_init_ref(ref);

    qua_bitop_spec_t spec = {
        "invsqrt", QUA_COST_QC, true, 0u, 36u, 100.0, 6u, 1024u, 6u,
        qua_invsqrt_bitop_eval, ref, qua_invsqrt_bitop_relax_2jet
    };

    qua_invsqrt_cfg_t *best = (qua_invsqrt_cfg_t *)malloc(sizeof(qua_invsqrt_cfg_t));
    qua_invsqrt_cfg_t *cand = (qua_invsqrt_cfg_t *)malloc(sizeof(qua_invsqrt_cfg_t));
    if (!best || !cand) { free(cand); free(best); free(ref); return 1; }
    memset(best, 0, sizeof(*best));
    best->est_score = 1e30;
    best->worst_err_ratio = 1e9;

    /* Polish 8 fast 1024-step restarts across both algebraic basins on all 3,133,452 exact inputs */
    static const uint64_t PIN_PHASESPLIT_CUC = (UINT64_C(1) << 34u) | (UINT64_C(1) << 35u);
    static const uint64_t PIN_SPLITC2        = (UINT64_C(1) << 28u) | (UINT64_C(1) << 34u);
    static const uint64_t SEED_FULL_KM8      =
        UINT64_C(5) | (UINT64_C(6) << 3u) | (UINT64_C(4) << 6u) | (UINT64_C(1) << 9u) |
        (UINT64_C(3) << 12u) | (UINT64_C(3) << 15u) | (UINT64_C(7) << 18u) |
        (UINT64_C(2) << 21u) | (UINT64_C(2) << 24u) | (UINT64_C(2) << 26u);
    static const uint64_t SEED_PS_CUC_KM7_F35 =
        UINT64_C(3) | (UINT64_C(5) << 3u) | (UINT64_C(4) << 6u) | (UINT64_C(1) << 9u) |
        (UINT64_C(3) << 12u) | (UINT64_C(3) << 15u) | (UINT64_C(6) << 18u) |
        (UINT64_C(2) << 21u) | (UINT64_C(2) << 24u) | (UINT64_C(2) << 26u) |
        (UINT64_C(1) << 30u) | (UINT64_C(1) << 34u) | (UINT64_C(1) << 35u);
    static const uint64_t SEED_PS_CUC_KM8_F35 =
        UINT64_C(5) | (UINT64_C(6) << 3u) | (UINT64_C(4) << 6u) | (UINT64_C(1) << 9u) |
        (UINT64_C(3) << 12u) | (UINT64_C(3) << 15u) | (UINT64_C(6) << 18u) |
        (UINT64_C(2) << 21u) | (UINT64_C(2) << 24u) | (UINT64_C(2) << 26u) |
        (UINT64_C(1) << 30u) | (UINT64_C(1) << 34u) | (UINT64_C(1) << 35u);
    static const uint64_t SEED_SC2_KM7       =
        UINT64_C(3) | (UINT64_C(6) << 3u) | (UINT64_C(4) << 6u) | (UINT64_C(1) << 9u) |
        (UINT64_C(3) << 12u) | (UINT64_C(3) << 15u) | (UINT64_C(7) << 18u) |
        (UINT64_C(2) << 21u) | (UINT64_C(2) << 24u) | (UINT64_C(2) << 26u) |
        (UINT64_C(1) << 28u);
    static const uint64_t SEED_SC2_KM8       =
        UINT64_C(5) | (UINT64_C(6) << 3u) | (UINT64_C(4) << 6u) | (UINT64_C(1) << 9u) |
        (UINT64_C(3) << 12u) | (UINT64_C(3) << 15u) | (UINT64_C(7) << 18u) |
        (UINT64_C(2) << 21u) | (UINT64_C(2) << 24u) | (UINT64_C(2) << 26u) |
        (UINT64_C(1) << 28u);

    searu_hyperparams hp = SEARU_ACTIVE_OPT_HP;
    for (uint32_t attempt = 0u; attempt < 8u; ++attempt) {
        uint64_t g_try[4] = {0u, 0u, 0u, 0u};
        double raw_try = 1e30;
        uint64_t a_seed = UINT64_C(0x51A7BEEF) + attempt * UINT64_C(0x9E3779B97F4A7C15);
        if (attempt == 0u) {
            qua_invsqrt_run_seeded_inner(&spec, &hp, a_seed, 1024u, SEED_FULL_KM8,
                                         (UINT64_C(1) << 34u), g_try, &raw_try);
        } else if (attempt == 1u) {
            qua_invsqrt_run_seeded_inner(&spec, &hp, a_seed, 1024u, SEED_FULL_KM8 | (UINT64_C(1) << 35u),
                                         (UINT64_C(1) << 34u) | (UINT64_C(1) << 35u), g_try, &raw_try);
        } else if (attempt == 2u) {
            qua_invsqrt_run_seeded_inner(&spec, &hp, a_seed, 1024u, SEED_SC2_KM8,
                                         PIN_SPLITC2, g_try, &raw_try);
        } else if (attempt == 3u) {
            qua_invsqrt_run_seeded_inner(&spec, &hp, a_seed, 1024u, SEED_SC2_KM7,
                                         PIN_SPLITC2, g_try, &raw_try);
        } else if (attempt == 4u) {
            qua_invsqrt_run_seeded_inner(&spec, &hp, a_seed, 1024u, SEED_PS_CUC_KM7_F35,
                                         PIN_PHASESPLIT_CUC, g_try, &raw_try);
        } else if (attempt == 5u) {
            qua_invsqrt_run_seeded_inner(&spec, &hp, a_seed, 1024u, SEED_PS_CUC_KM8_F35,
                                         PIN_PHASESPLIT_CUC, g_try, &raw_try);
        } else {
            qua_bitop_run_inner(&spec, &hp, a_seed, 1024u, g_try, &raw_try);
        }
        qua_invsqrt_decode_genome(g_try[0], cand);
        qua_invsqrt_calibrate(cand, ref);
        qua_invsqrt_center_charts_fast(cand, ref);
        qua_invsqrt_relax_charts_fast_2jet(cand, ref);
        qua_invsqrt_polish_charts(cand, ref);
        qua_invsqrt_eval_error(cand, ref, true);
        qua_invsqrt_plan_tiers_and_cost(cand);
        bool cand_ok = cand->tiers_valid && (cand->worst_err_ratio < 1.0) && (cand->min_margin > 0);
        bool best_ok = best->tiers_valid && (best->worst_err_ratio < 1.0) && (best->min_margin > 0);
        printf("[invsqrt Polish #%u] K_m=%u w=(%u,%u,%u) w_res=%u wd2=%u wsq=%u F=%u psplit=%u cuc=%u split_c2=%u b_sh_sq=%u b_sh_main=%u zprod=%u sq_ar=%u | "
               "ok=%d full_err=%.4f Q=%u EstTof=%llu EstScore=%.1f\n",
               attempt, cand->k_m, cand->w0, cand->w1, cand->w2, cand->w_res, cand->wd2, cand->wsq, cand->f_out,
               cand->phase_split_mode, cand->cuccaro_mode, cand->split_c2, cand->b_sh_sq, cand->b_sh_main,
               cand->zero_prods_mode, qua_invsqrt_use_sq_arith(cand) ? 1u : 0u,
               (int)cand_ok, cand->worst_err_ratio, cand->peak_qubits,
               (unsigned long long)cand->est_toffolis, cand->est_score);
        if ((cand_ok && (!best_ok || cand->est_score < best->est_score)) ||
            (!best_ok && !cand_ok && cand->worst_err_ratio < best->worst_err_ratio)) {
            *best = *cand;
        }
    }
    free(cand);

    printf("[invsqrt Pure Search Winner] K_m=%u w=(%u,%u,%u) w_res=%u wd2=%u wsq=%u F=%u psplit=%u cuc=%u split_c2=%u b_sh_sq=%u b_sh_main=%u zprod=%u sq_ar=%u | "
           "full_err=%.4f (<1.0), min_margin=%lld, Q=%u, EstTof=%llu, EstScore=%.1f\n",
           best->k_m, best->w0, best->w1, best->w2, best->w_res, best->wd2, best->wsq, best->f_out,
           best->phase_split_mode, best->cuccaro_mode, best->split_c2, best->b_sh_sq, best->b_sh_main,
           best->zero_prods_mode, qua_invsqrt_use_sq_arith(best) ? 1u : 0u,
           best->worst_err_ratio,
           (long long)best->min_margin, best->peak_qubits,
           (unsigned long long)best->est_toffolis, best->est_score);

    int rc = 1;
    if (best->tiers_valid && best->worst_err_ratio < 1.0 && best->min_margin > 0) {
        rc = qua_invsqrt_emit_kmb(best, out_kmb_path);
    } else {
        fprintf(stderr, "[ERROR] invsqrt pure search did not find a valid full-grid solution!\n");
    }

    free(best);
    free(ref);
    return rc;
}

#endif /* QUA_INVSQRT_H */
