/*
 * qua_dqi.h — Pure-Search Bilevel SEARU x FLOW 2-Jet Benchmark II (`dqi`) Synthesizer
 *
 * Replaces all hardcoded parameters and heuristics in `build_dqi_v4.py` with
 * Bilevel Meta-SEARU + FLOW 2-Jet search over the Half-EEA + Toeplitz + Vandermonde genome:
 *   - (R1, S1, K_MIN, halt_start, toeplitz_base, vand_mode)
 *   - Exact R(eps_plus) * Q^2 * T cost law from `verify_submission_dqi.py`
 *   - Dynamic per-candidate frustration compass (`frust`)
 *   - SEARU Vandermonde evaluation-point ordering
 *   - Direct pure-C17 Kickmix v2 (.kmb) emission via `qua_kmb.h` & `qua_ele.h`
 */
#ifndef QUA_DQI_H
#define QUA_DQI_H

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

#define QUA_DQI_M            4095u
#define QUA_DQI_N            70u
#define QUA_DQI_DEG          12u
#define QUA_DQI_R            2016u
#define QUA_DQI_L            35u
#define QUA_DQI_MOD_NATIVE   0x1053u
#define QUA_DQI_MOD_GF8      0xBu
#define QUA_DQI_MAX_PTS      80u
#define QUA_DQI_POOL_CAP     4096u

/* ============================================================================
 * 1. Qubit & Classical Bit Pool + GF(2^12) / Tower Algebra Tables
 * ============================================================================ */

typedef struct {
    qua_builder_t *kb;
    uint32_t next_q;
    uint32_t peak_q;
    uint32_t free_q[QUA_DQI_POOL_CAP];
    uint32_t num_free_q;
    uint32_t next_b;
    uint32_t free_b[QUA_DQI_POOL_CAP];
    uint32_t num_free_b;
} qua_dqi_ctx_t;

static inline void qua_dqi_ctx_init(qua_dqi_ctx_t *ctx, qua_builder_t *kb, uint32_t base_q, uint32_t base_b) {
    memset(ctx, 0, sizeof(*ctx));
    ctx->kb     = kb;
    ctx->next_q = base_q;
    ctx->peak_q = base_q;
    ctx->next_b = base_b;
}

static inline uint32_t qua_dqi_alloc_q1(qua_dqi_ctx_t *ctx) {
    if (ctx->num_free_q > 0u) {
        return ctx->free_q[--ctx->num_free_q];
    }
    uint32_t q = ctx->next_q++;
    if (ctx->next_q > ctx->peak_q) ctx->peak_q = ctx->next_q;
    return q;
}

static inline void qua_dqi_alloc_q(qua_dqi_ctx_t *ctx, uint32_t *out, uint32_t n) {
    for (uint32_t i = 0u; i < n; ++i) out[i] = qua_dqi_alloc_q1(ctx);
}

static inline void qua_dqi_free_q1(qua_dqi_ctx_t *ctx, uint32_t q) {
    if (ctx->num_free_q < QUA_DQI_POOL_CAP) {
        ctx->free_q[ctx->num_free_q++] = q;
    }
}

static inline void qua_dqi_free_q(qua_dqi_ctx_t *ctx, const uint32_t *qs, uint32_t n) {
    for (uint32_t i = n; i > 0u; --i) {
        qua_dqi_free_q1(ctx, qs[i - 1u]);
    }
}

static inline uint32_t qua_dqi_alloc_b1(qua_dqi_ctx_t *ctx) {
    if (ctx->num_free_b > 0u) {
        return ctx->free_b[--ctx->num_free_b];
    }
    return ctx->next_b++;
}

static inline void qua_dqi_alloc_b(qua_dqi_ctx_t *ctx, uint32_t *out, uint32_t n) {
    for (uint32_t i = 0u; i < n; ++i) out[i] = qua_dqi_alloc_b1(ctx);
}

static inline void qua_dqi_free_b(qua_dqi_ctx_t *ctx, const uint32_t *bs, uint32_t n) {
    for (uint32_t i = n; i > 0u; --i) {
        if (ctx->num_free_b < QUA_DQI_POOL_CAP) {
            ctx->free_b[ctx->num_free_b++] = bs[i - 1u];
        }
    }
}

typedef struct {
    double   weight_probs[QUA_DQI_L + 1u];
    uint16_t gf12_exp[4096];
    uint16_t gf12_log[4096];
    qua_ele_gf12_rank42_t mul42;
    qua_ele_gf64_rank15_t mul64;
    double   cp_eps[32];

    /* Tower GF(((2^3)^2)^2) isomorphism & inversion tables */
    uint8_t  beta_gf64;
    uint16_t t_nat2tow[12];  /* 12 rows, each 12-bit mask */
    uint16_t t_tow2nat[12];
    uint16_t m_delta_lin[6]; /* 6 rows, each 12-bit mask */
    uint8_t  m_d_lin[3];     /* 3 rows, each 6-bit mask */
    uint8_t  l_gf8_inv[3];   /* 3 rows, each 3-bit mask */
    uint8_t  q_gf8_inv[3];   /* 3 rows, each 3-bit mask */

    /* GF(8) rank-6 bilinear plan (6 terms) */
    uint8_t  gf8_u[6];
    uint8_t  gf8_v[6];
    uint8_t  gf8_w[6];

    /* Precomputed GF(64) phase table: gf64_phase[k][i] is 6-bit mask on j */
    uint8_t  gf64_phase[6][6];
    /* Precomputed GF(8) phase table: gf8_phase[k][i] is 3-bit mask on j */
    uint8_t  gf8_phase[3][3];
    /* Precomputed GF(2^12) phase table: gf12_phase[k][i] is 12-bit mask on j */
    uint16_t gf12_phase[12][12];

    /* Good evaluation points in GF(2^12) for Vandermonde */
    uint16_t eval_pts[QUA_DQI_MAX_PTS];
} qua_dqi_tables_t;

static inline uint8_t qua_dqi_gf64_mul(uint8_t a, uint8_t b) {
    uint8_t a0 = a & 7u, a1 = (a >> 3u) & 7u;
    uint8_t b0 = b & 7u, b1 = (b >> 3u) & 7u;
    uint8_t p0 = qua_gf8_mul(a0, b0);
    uint8_t p1 = qua_gf8_mul(a1, b1);
    uint8_t pm = qua_gf8_mul((uint8_t)(a0 ^ a1), (uint8_t)(b0 ^ b1));
    return (uint8_t)((p0 ^ p1) | ((pm ^ p0) << 3u));
}

static inline uint16_t qua_dqi_tower_mul(uint16_t a, uint16_t b, uint8_t beta) {
    uint8_t a0 = (uint8_t)(a & 63u), a1 = (uint8_t)((a >> 6u) & 63u);
    uint8_t b0 = (uint8_t)(b & 63u), b1 = (uint8_t)((b >> 6u) & 63u);
    uint8_t p0 = qua_dqi_gf64_mul(a0, b0);
    uint8_t p1 = qua_dqi_gf64_mul(a1, b1);
    uint8_t pm = qua_dqi_gf64_mul((uint8_t)(a0 ^ a1), (uint8_t)(b0 ^ b1));
    uint8_t c0 = p0 ^ qua_dqi_gf64_mul(beta, p1);
    uint8_t c1 = pm ^ p0;
    return (uint16_t)(c0 | ((uint16_t)c1 << 6u));
}

static inline uint16_t qua_dqi_tower_pow(uint16_t a, uint32_t e, uint8_t beta) {
    uint16_t res = 1u;
    uint16_t base = a;
    while (e > 0u) {
        if (e & 1u) res = qua_dqi_tower_mul(res, base, beta);
        base = qua_dqi_tower_mul(base, base, beta);
        e >>= 1u;
    }
    return res;
}

static inline uint16_t qua_dqi_gf12_pow(const qua_dqi_tables_t *tab, uint16_t a, uint32_t e) {
    if (a == 0u) return (e == 0u) ? 1u : 0u;
    uint32_t lg = tab->gf12_log[a];
    uint32_t idx = (uint32_t)(((uint64_t)lg * (uint64_t)e) % 4095ull);
    return tab->gf12_exp[idx];
}

static inline bool qua_dqi_invert_f2_12(const uint16_t in_rows[12], uint16_t out_rows[12]) {
    uint16_t mat[12], inv[12];
    for (uint32_t r = 0u; r < 12u; ++r) {
        mat[r] = in_rows[r] & 0xFFFu;
        inv[r] = (uint16_t)(1u << r);
    }
    for (uint32_t col = 0u; col < 12u; ++col) {
        uint32_t piv = col;
        while (piv < 12u && (((mat[piv] >> col) & 1u) == 0u)) ++piv;
        if (piv == 12u) return false;
        if (piv != col) {
            uint16_t tm = mat[col]; mat[col] = mat[piv]; mat[piv] = tm;
            uint16_t ti = inv[col]; inv[col] = inv[piv]; inv[piv] = ti;
        }
        for (uint32_t r = 0u; r < 12u; ++r) {
            if (r != col && ((mat[r] >> col) & 1u)) {
                mat[r] ^= mat[col];
                inv[r] ^= inv[col];
            }
        }
    }
    for (uint32_t r = 0u; r < 12u; ++r) out_rows[r] = inv[r];
    return true;
}

/* Compute principal eigenvector of (L+1)x(L+1) symmetric tridiagonal DQI matrix via Jacobi */
static inline void qua_dqi_compute_weights(double probs[QUA_DQI_L + 1u]) {
    const uint32_t N = QUA_DQI_L + 1u; /* 36 */
    double A[36][36];
    double V[36][36];
    memset(A, 0, sizeof(A));
    memset(V, 0, sizeof(V));

    double q = 4096.0;
    double r = (double)QUA_DQI_R;
    double d = (q - 2.0 * r) / sqrt(r * (q - r));
    for (uint32_t i = 0u; i < N; ++i) {
        V[i][i] = 1.0;
        A[i][i] = (double)i * d;
        if (i + 1u < N) {
            double k = (double)(i + 1u);
            double off = sqrt(k * ((double)QUA_DQI_M - k + 1.0));
            A[i][i + 1u] = off;
            A[i + 1u][i] = off;
        }
    }

    for (uint32_t sweep = 0u; sweep < 40u; ++sweep) {
        double max_off = 0.0;
        for (uint32_t p = 0u; p < N; ++p) {
            for (uint32_t q_i = p + 1u; q_i < N; ++q_i) {
                double apq = A[p][q_i];
                if (fabs(apq) < 1e-13) continue;
                if (fabs(apq) > max_off) max_off = fabs(apq);
                double app = A[p][p];
                double aqq = A[q_i][q_i];
                double tau = (aqq - app) / (2.0 * apq);
                double t = (tau >= 0.0)
                    ? (1.0 / (tau + sqrt(1.0 + tau * tau)))
                    : (-1.0 / (-tau + sqrt(1.0 + tau * tau)));
                double c = 1.0 / sqrt(1.0 + t * t);
                double s = t * c;

                A[p][p] = app - t * apq;
                A[q_i][q_i] = aqq + t * apq;
                A[p][q_i] = 0.0;
                A[q_i][p] = 0.0;
                for (uint32_t k = 0u; k < N; ++k) {
                    if (k != p && k != q_i) {
                        double akp = A[k][p];
                        double akq = A[k][q_i];
                        A[k][p] = c * akp - s * akq;
                        A[p][k] = A[k][p];
                        A[k][q_i] = s * akp + c * akq;
                        A[q_i][k] = A[k][q_i];
                    }
                    double vkp = V[k][p];
                    double vkq = V[k][q_i];
                    V[k][p] = c * vkp - s * vkq;
                    V[k][q_i] = s * vkp + c * vkq;
                }
            }
        }
        if (max_off < 1e-13) break;
    }

    uint32_t best_col = 0u;
    for (uint32_t i = 1u; i < N; ++i) {
        if (A[i][i] > A[best_col][best_col]) best_col = i;
    }
    double norm2 = 0.0;
    for (uint32_t i = 0u; i < N; ++i) {
        probs[i] = V[i][best_col] * V[i][best_col];
        norm2 += probs[i];
    }
    for (uint32_t i = 0u; i < N; ++i) probs[i] /= norm2;
}

static inline void qua_dqi_init_tables(qua_dqi_tables_t *tab) {
    memset(tab, 0, sizeof(*tab));
    qua_dqi_compute_weights(tab->weight_probs);

    /* Log / Antilog for GF(2^12) with prim = 2 */
    uint16_t cur = 1u;
    for (uint32_t i = 0u; i < 4095u; ++i) {
        tab->gf12_exp[i] = cur;
        tab->gf12_log[cur] = (uint16_t)i;
        cur = qua_gf12_mul_ref(cur, 2u, QUA_DQI_MOD_NATIVE);
    }
    qua_ele_build_gf12_rank42(QUA_DQI_MOD_NATIVE, &tab->mul42);
    qua_ele_build_gf64_rank15(qua_dqi_gf64_mul, &tab->mul64);

    /* Exact Clopper-Pearson upper bounds eps^+(k; N=102400, delta=2^-128) for k = 0..31 */
    {
        const double n_sh = 102400.0;
        const double log_delta = -128.0 * log(2.0);
        const double lg_n1 = lgamma(n_sh + 1.0);
        tab->cp_eps[0] = -expm1(log_delta / n_sh);
        for (uint32_t k = 1u; k < 32u; ++k) {
            double lo = (double)k / n_sh;
            double hi = 0.05;
            for (uint32_t it = 0u; it < 52u; ++it) {
                double mid = 0.5 * (lo + hi);
                double lp = log(mid), l1p = log1p(-mid);
                double max_lt = -1e300;
                double lterms[32];
                for (uint32_t j = 0u; j <= k; ++j) {
                    double lt = lg_n1 - lgamma((double)j + 1.0) - lgamma(n_sh - (double)j + 1.0)
                                + (double)j * lp + (n_sh - (double)j) * l1p;
                    lterms[j] = lt;
                    if (lt > max_lt) max_lt = lt;
                }
                double sum_e = 0.0;
                for (uint32_t j = 0u; j <= k; ++j) sum_e += exp(lterms[j] - max_lt);
                double log_cdf = max_lt + log(sum_e);
                if (log_cdf < log_delta) hi = mid;
                else lo = mid;
            }
            tab->cp_eps[k] = 0.5 * (lo + hi);
        }
    }

    /* Irreducible beta in GF(64) for w^2 + w + beta */
    for (uint8_t b = 1u; b < 64u; ++b) {
        bool has_root = false;
        for (uint8_t x = 0u; x < 64u; ++x) {
            if ((qua_dqi_gf64_mul(x, x) ^ x ^ b) == 0u) {
                has_root = true;
                break;
            }
        }
        if (!has_root) { tab->beta_gf64 = b; break; }
    }

    /* Find root r in Tower of x^12 + x^6 + x^4 + x + 1 = 0 -> T_TOW2NAT / T_NAT2TOW */
    for (uint16_t r = 2u; r < 4096u; ++r) {
        uint16_t r2 = qua_dqi_tower_mul(r, r, tab->beta_gf64);
        uint16_t r4 = qua_dqi_tower_mul(r2, r2, tab->beta_gf64);
        uint16_t r6 = qua_dqi_tower_mul(r4, r2, tab->beta_gf64);
        uint16_t r12 = qua_dqi_tower_mul(r6, r6, tab->beta_gf64);
        if ((r12 ^ r6 ^ r4 ^ r ^ 1u) == 0u) {
            uint16_t cols[12];
            uint16_t pw = 1u;
            for (uint32_t k = 0u; k < 12u; ++k) {
                cols[k] = pw;
                pw = qua_dqi_tower_mul(pw, r, tab->beta_gf64);
            }
            uint16_t rows[12];
            for (uint32_t row = 0u; row < 12u; ++row) {
                uint16_t m = 0u;
                for (uint32_t col = 0u; col < 12u; ++col) {
                    if ((cols[col] >> row) & 1u) m |= (uint16_t)(1u << col);
                }
                rows[row] = m;
            }
            if (qua_dqi_invert_f2_12(rows, tab->t_tow2nat)) {
                for (uint32_t row = 0u; row < 12u; ++row) tab->t_nat2tow[row] = rows[row];
                break;
            }
        }
    }

    /* M_DELTA_LIN (6x12): a^2 + beta*b^2 */
    for (uint32_t bit = 0u; bit < 6u; ++bit) {
        uint8_t a_sq = qua_dqi_gf64_mul((uint8_t)(1u << bit), (uint8_t)(1u << bit));
        uint8_t b_sq = qua_dqi_gf64_mul(tab->beta_gf64, a_sq);
        for (uint32_t r = 0u; r < 6u; ++r) {
            if ((a_sq >> r) & 1u) tab->m_delta_lin[r] |= (uint16_t)(1u << bit);
            if ((b_sq >> r) & 1u) tab->m_delta_lin[r] |= (uint16_t)(1u << (6u + bit));
        }
    }
    /* M_D_LIN (3x6): u0^2 + u1^2 */
    for (uint32_t bit = 0u; bit < 3u; ++bit) {
        uint8_t sq = qua_gf8_mul((uint8_t)(1u << bit), (uint8_t)(1u << bit));
        for (uint32_t r = 0u; r < 3u; ++r) {
            if ((sq >> r) & 1u) {
                tab->m_d_lin[r] |= (uint8_t)(1u << bit);
                tab->m_d_lin[r] |= (uint8_t)(1u << (3u + bit));
            }
        }
    }
    /* GF(8) inversion: d^{-1} = d^6 = L*d + Q*[d0d1, d0d2, d1d2] */
    uint8_t lin_vals[3];
    for (uint32_t b = 0u; b < 3u; ++b) {
        lin_vals[b] = qua_gf8_inv((uint8_t)(1u << b));
        for (uint32_t r = 0u; r < 3u; ++r) {
            if ((lin_vals[b] >> r) & 1u) tab->l_gf8_inv[r] |= (uint8_t)(1u << b);
        }
    }
    uint8_t q_pairs[3] = {3u, 5u, 6u};
    uint8_t q_b0[3] = {1u, 1u, 2u};
    uint8_t q_b1[3] = {2u, 4u, 4u};
    for (uint32_t p = 0u; p < 3u; ++p) {
        uint8_t inv_p = qua_gf8_inv(q_pairs[p]);
        uint8_t lin_p = qua_gf8_inv(q_b0[p]) ^ qua_gf8_inv(q_b1[p]);
        uint8_t quad  = inv_p ^ lin_p;
        for (uint32_t r = 0u; r < 3u; ++r) {
            if ((quad >> r) & 1u) tab->q_gf8_inv[r] |= (uint8_t)(1u << p);
        }
    }

    /* GF(8) 6-term Karatsuba plan */
    static const uint8_t GF8_U[6] = {1u, 2u, 4u, 3u, 5u, 6u};
    static const uint8_t GF8_W[3] = {
        /* z0 = m0+m1+m2+m12 */ 0x27u,
        /* z1 = m0+m01+m12   */ 0x29u,
        /* z2 = m0+m1+m02    */ 0x13u
    };
    for (uint32_t t = 0u; t < 6u; ++t) {
        tab->gf8_u[t] = GF8_U[t];
        tab->gf8_v[t] = GF8_U[t];
        uint8_t w = 0u;
        for (uint32_t r = 0u; r < 3u; ++r) {
            if ((GF8_W[r] >> t) & 1u) w |= (uint8_t)(1u << r);
        }
        tab->gf8_w[t] = w;
    }

    /* Precompute bitwise phase-kickback tables for 0-CCX HMR uncomputation over
     * GF(64), GF(8), and GF(2^12) (paired with Rank-15 mul64 and Rank-42 mul42). */
    for (uint32_t i = 0u; i < 6u; ++i) {
        for (uint32_t j = 0u; j < 6u; ++j) {
            uint8_t p = qua_dqi_gf64_mul((uint8_t)(1u << i), (uint8_t)(1u << j));
            for (uint32_t k = 0u; k < 6u; ++k) {
                if ((p >> k) & 1u) tab->gf64_phase[k][i] |= (uint8_t)(1u << j);
            }
        }
    }
    for (uint32_t i = 0u; i < 3u; ++i) {
        for (uint32_t j = 0u; j < 3u; ++j) {
            uint8_t p = qua_gf8_mul((uint8_t)(1u << i), (uint8_t)(1u << j));
            for (uint32_t k = 0u; k < 3u; ++k) {
                if ((p >> k) & 1u) tab->gf8_phase[k][i] |= (uint8_t)(1u << j);
            }
        }
    }
    for (uint32_t i = 0u; i < 12u; ++i) {
        for (uint32_t j = 0u; j < 12u; ++j) {
            uint16_t p = qua_gf12_mul_ref((uint16_t)(1u << i), (uint16_t)(1u << j), QUA_DQI_MOD_NATIVE);
            for (uint32_t k = 0u; k < 12u; ++k) {
                if ((p >> k) & 1u) tab->gf12_phase[k][i] |= (uint16_t)(1u << j);
            }
        }
    }

    /* Select 80 low-CNOT evaluation points in GF(2^12): {0, 1, 2, 4, ...} */
    tab->eval_pts[0] = 0u;
    uint32_t cnt = 1u;
    for (uint32_t w = 1u; w <= 12u && cnt < QUA_DQI_MAX_PTS; ++w) {
        for (uint16_t x = 1u; x < 4096u && cnt < QUA_DQI_MAX_PTS; ++x) {
            if ((uint32_t)__builtin_popcount((unsigned)x) == w) {
                tab->eval_pts[cnt++] = x;
            }
        }
    }
}

/* ============================================================================
 * 2. Pure-Search Architectural Evaluation & Bilevel SEARU x 2-Jet Optimizer
 * ============================================================================ */

typedef struct {
    uint32_t r1;
    uint32_t s1;
    uint32_t k_min;
    uint32_t halt_start;
    uint32_t toeplitz_base;
    bool     pebble_scaff;
    bool     mbuc_q0;
    bool     three_stage;
    uint32_t r0;
    uint32_t s0;
    uint32_t peak_qubits;
    uint64_t est_toffolis;
    double   p_trunc;
    double   eps_plus;
    double   rep_factor;
    double   est_score;
    bool     valid;
} qua_dqi_cfg_t;

static inline uint64_t qua_dqi_est_toeplitz_ccx(uint32_t m, uint32_t t_base) {
    if (m <= t_base) {
        /* Per step k=0..m-1: 1 mul42 + 12 cswap + k mul42 (0-CCX HMR_Dual erasure) */
        return (uint64_t)m * (42ull + 12ull) + ((uint64_t)m * (uint64_t)(m - 1u) / 2ull) * 42ull;
    }
    uint32_t h = m >> 1u;
    uint32_t slice_m = (m - 2u) + h;
    return qua_dqi_est_toeplitz_ccx(h, t_base)
         + (uint64_t)slice_m * 42ull
         + qua_dqi_est_toeplitz_ccx(m - h, t_base);
}

static inline uint64_t qua_dqi_est_round_pair_ccx(
    uint32_t r,
    uint32_t lo,
    uint32_t h_start,
    uint64_t inv_pair_ccx,
    bool mbuc_q0
) {
    uint64_t flags_ccx = (r == 0u)
        ? 23ull
        : ((r <= h_start) ? 24ull : ((r < 34u) ? 26ull : 12ull));
    if (r == 0u) {
        uint64_t mac_quo_r0 = (uint64_t)(76u - lo) * 42ull; /* (71 - lo) poly_mac + 5 quo/tail */
        uint64_t csw_r0     = 12ull * (uint64_t)(71u - lo);
        return 2ull * (flags_ccx + 2ull * inv_pair_ccx + mac_quo_r0 + csw_r0);
    }
    uint64_t mac_r  = (uint64_t)(74u - lo) * 42ull;
    uint64_t q2_q1  = 192ull; /* 54 (q2 fwd) + 84 (q1 fwd) + 42 (q1 unquo1) + 12 (q2 hmr) */
    uint64_t q0_fwd = ((mbuc_q0 ? 3ull : 5ull) + ((r == 1u) ? 1ull : 2ull)) * 42ull;
    uint64_t q0_rev = (((r == 1u) ? (mbuc_q0 ? 2ull : 3ull) : (mbuc_q0 ? 3ull : 5ull)) + 2ull) * 42ull;
    uint64_t sh_csw = 1ull + 24ull * (2ull + ((r == 1u) ? 1ull : ((r == 2u) ? 2ull : 3ull)))
                    + 12ull * (uint64_t)(71u - lo);
    return 2ull * (flags_ccx + 2ull * inv_pair_ccx + mac_r + q2_q1 + sh_csw) + q0_fwd + q0_rev;
}

static inline void qua_dqi_eval_cfg(
    const qua_dqi_tables_t *tab,
    uint64_t word0,
    qua_dqi_cfg_t *cfg
) {
    memset(cfg, 0, sizeof(*cfg));
    uint32_t r1_raw          = (uint32_t)(word0 & 15u);          /* [0..3]:   16..30 */
    uint32_t hs_raw          = (uint32_t)((word0 >> 4u) & 31u);  /* [4..8]:   6..28  */
    uint32_t tb_raw          = (uint32_t)((word0 >> 9u) & 3u);   /* [9..10]:  2..5   */
    /* [11..14] pinned off: enforce k_off=0, d_s1=0 so Phase-1 HGCD never drops low syndromes */
    uint32_t k_off           = 0u;
    uint32_t d_s1            = 0u;
    bool     peb_sc          = ((word0 >> 15u) & 1u) != 0u;      /* [15]:     0..1   */
    bool     mbuc_q0         = ((word0 >> 16u) & 1u) != 0u;      /* [16]:     0..1   */
    bool     three_stage_req = ((word0 >> 17u) & 1u) != 0u;      /* [17]:     0..1   */
    uint32_t r0_raw          = (uint32_t)((word0 >> 18u) & 15u); /* [18..21]: 0..15  */

    uint32_t r1 = 16u + (r1_raw % 15u); /* 16..30 */
    int32_t  s1_i = (int32_t)(70u - 2u * r1) - (int32_t)d_s1;
    uint32_t s1 = (s1_i >= 8) ? (uint32_t)s1_i : 0u;
    uint32_t k_min = (s1 < r1) ? s1 : r1;
    if (k_min > k_off + 8u) k_min -= k_off;
    if (k_min > s1) k_min = s1;
    if (k_min > r1) k_min = r1;

    uint32_t h_start = 6u + (hs_raw % 23u); /* 6..28 */
    uint32_t t_base  = 2u + tb_raw;         /* 2..5  */

    cfg->r1            = r1;
    cfg->s1            = s1;
    cfg->k_min         = k_min;
    cfg->halt_start    = h_start;
    cfg->toeplitz_base = t_base;
    cfg->pebble_scaff  = peb_sc;
    cfg->mbuc_q0       = mbuc_q0;
    cfg->three_stage   = false;
    cfg->r0            = 0u;
    cfg->s0            = s1;

    if (s1 < 8u || r1 < 16u || r1 + s1 > 71u) {
        cfg->valid = false;
        cfg->est_score = 1e30;
        return;
    }

    if (three_stage_req && r1 >= 12u && s1 < 35u) {
        uint32_t r0 = 8u + (r0_raw % (r1 - 8u));
        uint32_t s0 = (2u * r0 < 70u) ? (70u - 2u * r0) : 0u;
        if (s0 > 35u) s0 = 35u;
        if (s0 > s1) {
            uint32_t r_bounds[3] = {r0, r1, 35u};
            uint32_t s_floors[3] = {s0, s1, k_min};
            if (qua_bitop_verify_hgcd_bands(70u, 35u, 3u, r_bounds, s_floors)) {
                cfg->three_stage = true;
                cfg->r0          = r0;
                cfg->s0          = s0;
            }
        }
    }

    /* Exact Peak Qubits matching `qua_dqi_emit_kmb`:
     * Persistent: 840 (syn) + 12*(35 - k_min) [B_lo] + 12 [b70] + 12*(71 - k_min) [A_pool]
     *           + (35 - h_start) [do_div_hi] + 34 [jump_at]
     * Transient EEA scratch peak:
     *   pebble_scaff=1: 53 qubits (mbuc_q0=0) or 65 qubits (mbuc_q0=1)
     *   pebble_scaff=0: 71 qubits (mbuc_q0=0) or 83 qubits (mbuc_q0=1) */
    uint32_t b_lo_regs   = 35u - k_min;
    uint32_t a_work_regs = 71u - k_min;
    uint32_t scratch_q   = peb_sc ? (mbuc_q0 ? 65u : 53u) : (mbuc_q0 ? 83u : 71u);
    uint32_t eea_q = 840u + 12u * (b_lo_regs + 1u + a_work_regs) + (35u - h_start) + 34u + scratch_q;
    cfg->peak_qubits = eea_q;

    /* Exact Toffoli count matching `qua_dqi_emit_kmb` */
    uint64_t inv_pair_ccx = peb_sc ? (69ull + 36ull) : 69ull;
    uint64_t tof_p01 = 0ull;
    for (uint32_t r = 0u; r < r1; ++r) {
        uint32_t lo = (cfg->three_stage && r < cfg->r0) ? cfg->s0 : s1;
        tof_p01 += qua_dqi_est_round_pair_ccx(r, lo, h_start, inv_pair_ccx, mbuc_q0);
    }

    uint64_t tof_cp0 = 0ull;
    if (cfg->three_stage) {
        uint32_t len_c0 = cfg->s0 - s1;
        uint64_t tof_cp0_once = (uint64_t)(cfg->r0 + len_c0) * 42ull + 12ull
                              + (uint64_t)(cfg->r0 + len_c0 - 1u) * 42ull;
        tof_cp0 = 2ull * tof_cp0_once;
    }

    /* Midpoint Vandermonde (applied in both forward and reverse!) */
    uint64_t tof_mid_once = (uint64_t)(r1 + s1) * 42ull + 12ull + (uint64_t)(r1 + s1 - 1u) * 42ull;
    uint64_t tof_mid = 2ull * tof_mid_once;

    uint64_t tof_p2 = 0ull;
    for (uint32_t r = r1; r < 35u; ++r) {
        uint32_t lo_r = (k_min > (r - 1u)) ? k_min : (r - 1u);
        tof_p2 += qua_dqi_est_round_pair_ccx(r, lo_r, h_start, inv_pair_ccx, mbuc_q0);
    }

    /* Toeplitz (fwd + rev, each wraps 11 is_zero + tower_inv [69 or 105 if peb_sc]) */
    uint64_t tof_toep = 2ull * (80ull + (peb_sc ? 36ull : 0ull) + qua_dqi_est_toeplitz_ccx(35u, t_base));
    /* Chien-Forney: 4095 positions * (11 is_zero + 69 Rank-15 GF64 tower_inv = 80 CCX) */
    uint64_t tof_cf   = 4095ull * 80ull;

    uint64_t total_tof = tof_p01 + tof_cp0 + tof_mid + tof_p2 + tof_toep + tof_cf;
    cfg->est_toffolis = total_tof;

    /* Continuous Poisson-expected Clopper-Pearson repetition factor (zero llround cliff!) */
    double p_trunc = 0.0;
    for (uint32_t k = 0u; k < h_start && k <= QUA_DQI_L; ++k) {
        p_trunc += tab->weight_probs[k];
    }
    cfg->p_trunc = p_trunc;
    const double n_sh = 102400.0;
    double lam = p_trunc * n_sh;
    if (lam > 18.0) {
        cfg->eps_plus   = (p_trunc < 0.9) ? (p_trunc + 4.5 * sqrt(p_trunc / n_sh)) : 1.0;
        cfg->valid      = false;
        cfg->rep_factor = 1e6;
        cfg->est_score  = 1e30;
        return;
    }

    double prob_sum = 0.0;
    double mean_r   = 0.0;
    double mean_r2  = 0.0;
    double mean_eps = 0.0;
    double p_k = exp(-lam);
    for (uint32_t k = 0u; k < 32u; ++k) {
        double ep = tab->cp_eps[k];
        double rk = log(2.0) / -log(0.5 + 2.0 * sqrt(ep));
        prob_sum += p_k;
        mean_eps += p_k * ep;
        mean_r   += p_k * rk;
        mean_r2  += p_k * rk * rk;
        p_k *= lam / (double)(k + 1u);
    }
    if (prob_sum > 0.0) {
        mean_eps /= prob_sum;
        mean_r   /= prob_sum;
        mean_r2  /= prob_sum;
    }
    double var_r = mean_r2 - mean_r * mean_r;
    double std_r = (var_r > 0.0) ? sqrt(var_r) : 0.0;
    double eff_rep = mean_r + 2.5 * std_r; /* conservative 2.5-sigma upper-confidence R */

    cfg->eps_plus   = mean_eps;
    cfg->rep_factor = eff_rep;
    cfg->valid      = (mean_eps < 1.0 / 16.0) && (lam <= 0.25);
    cfg->est_score  = eff_rep * (double)eea_q * (double)eea_q * (double)total_tof;
}

/* ============================================================================
 * 3. Linear & Bilinear F2 Quantum Circuit Primitives
 * ============================================================================ */

static inline void qua_dqi_emit_inplace_f2(
    qua_builder_t *kb,
    const uint32_t *q,
    uint32_t n,
    const uint16_t *mat_rows /* n rows, each n bits */
) {
    if (n > 12u) abort();
    uint16_t A[12];
    for (uint32_t r = 0u; r < n; ++r) A[r] = mat_rows[r];

    typedef struct { uint8_t kind, a, b; } op_t;
    op_t ops[256];
    uint32_t nops = 0u;

    for (uint32_t col = 0u; col < n; ++col) {
        /* Minimum-weight pivot row selection to minimize CNOT fanout */
        uint32_t best_piv = col;
        uint32_t best_cost = UINT32_MAX;
        for (uint32_t r = col; r < n; ++r) {
            if (((A[r] >> col) & 1u) == 0u) continue;
            uint32_t cost = (uint32_t)__builtin_popcount((unsigned)A[r]) + ((r != col) ? 2u : 0u);
            if (cost < best_cost) {
                best_cost = cost;
                best_piv  = r;
            }
        }
        if (best_piv != col) {
            if (nops >= 256u) abort();
            uint16_t tmp = A[col]; A[col] = A[best_piv]; A[best_piv] = tmp;
            ops[nops++] = (op_t){1u, (uint8_t)col, (uint8_t)best_piv};
        }
        for (uint32_t r = 0u; r < n; ++r) {
            if (r != col && ((A[r] >> col) & 1u)) {
                if (nops >= 256u) abort();
                A[r] ^= A[col];
                ops[nops++] = (op_t){0u, (uint8_t)r, (uint8_t)col};
            }
        }
    }
    for (uint32_t k = nops; k > 0u; --k) {
        op_t op = ops[k - 1u];
        if (op.kind == 1u) {
            qua_cx(kb, q[op.a], q[op.b]);
            qua_cx(kb, q[op.b], q[op.a]);
            qua_cx(kb, q[op.a], q[op.b]);
        } else {
            qua_cx(kb, q[op.b], q[op.a]);
        }
    }
}

static inline void qua_dqi_gf12_imul_const(
    qua_builder_t *kb,
    const uint32_t q[12],
    uint16_t c_val
) {
    if (c_val == 1u) return;
    uint16_t rows[12];
    memset(rows, 0, sizeof(rows));
    for (uint32_t b = 0u; b < 12u; ++b) {
        uint16_t p = qua_gf12_mul_ref(c_val, (uint16_t)(1u << b), QUA_DQI_MOD_NATIVE);
        for (uint32_t r = 0u; r < 12u; ++r) {
            if ((p >> r) & 1u) rows[r] |= (uint16_t)(1u << b);
        }
    }
    qua_dqi_emit_inplace_f2(kb, q, 12u, rows);
}

static inline void qua_dqi_mul42(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t qa[12],
    const uint32_t qb[12],
    const uint32_t qout[12]
) {
    uint32_t anc = qua_dqi_alloc_q1(ctx);
    uint32_t hb  = qua_dqi_alloc_b1(ctx);
    qua_ele_emit_gf2_12_rank42_mul(ctx->kb, qa, qb, qout, anc, hb, &tab->mul42);
    qua_dqi_free_b(ctx, &hb, 1u);
    qua_dqi_free_q1(ctx, anc);
}

static inline void qua_dqi_phase_by_prod12(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t qa[12],
    const uint32_t qb[12],
    const uint32_t mask_bits[12]
) {
    for (uint32_t k = 0u; k < 12u; ++k) {
        for (uint32_t i = 0u; i < 12u; ++i) {
            uint16_t m = tab->gf12_phase[k][i];
            while (m) {
                uint32_t j = (uint32_t)__builtin_ctz((unsigned)m);
                m &= (uint16_t)(m - 1u);
                qua_cz_if(ctx->kb, qa[i], qb[j], mask_bits[k]);
            }
        }
    }
}

static inline void qua_dqi_unmul_gf12_hmr(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t tgt[12],
    const uint32_t lhs[12],
    const uint32_t rhs[12]
) {
    uint32_t mb[12];
    qua_dqi_alloc_b(ctx, mb, 12u);
    for (uint32_t k = 0u; k < 12u; ++k) qua_hmr(ctx->kb, tgt[k], mb[k]);
    qua_dqi_phase_by_prod12(ctx, tab, lhs, rhs, mb);
    qua_dqi_free_b(ctx, mb, 12u);
}

/* ============================================================================
 * 4. 12-bit Zero Test + Composite Tower GF(2^12) Inverter (HMR Uncompute)
 * ============================================================================ */

static inline void qua_dqi_init_is_zero_12(
    qua_dqi_ctx_t *ctx,
    const uint32_t reg[12],
    uint32_t ladder[11]
) {
    qua_dqi_alloc_q(ctx, ladder, 11u);
    for (uint32_t i = 0u; i < 12u; ++i) qua_x(ctx->kb, reg[i]);
    for (uint32_t i = 0u; i < 6u; ++i)  qua_and_compute(ctx->kb, reg[2u * i], reg[2u * i + 1u], ladder[i]);
    for (uint32_t i = 0u; i < 3u; ++i)  qua_and_compute(ctx->kb, ladder[2u * i], ladder[2u * i + 1u], ladder[6u + i]);
    qua_and_compute(ctx->kb, ladder[6], ladder[7], ladder[9]);
    qua_and_compute(ctx->kb, ladder[8], ladder[9], ladder[10]);
    for (uint32_t i = 0u; i < 12u; ++i) qua_x(ctx->kb, reg[i]);
}

static inline void qua_dqi_del_is_zero_12(
    qua_dqi_ctx_t *ctx,
    const uint32_t reg[12],
    const uint32_t ladder[11]
) {
    uint32_t sb = qua_dqi_alloc_b1(ctx);
    for (uint32_t i = 0u; i < 12u; ++i) qua_x(ctx->kb, reg[i]);
    qua_and_uncompute(ctx->kb, ladder[8], ladder[9], ladder[10], sb);
    qua_and_uncompute(ctx->kb, ladder[6], ladder[7], ladder[9], sb);
    for (uint32_t i = 3u; i > 0u; --i) {
        uint32_t k = i - 1u;
        qua_and_uncompute(ctx->kb, ladder[2u * k], ladder[2u * k + 1u], ladder[6u + k], sb);
    }
    for (uint32_t i = 6u; i > 0u; --i) {
        uint32_t k = i - 1u;
        qua_and_uncompute(ctx->kb, reg[2u * k], reg[2u * k + 1u], ladder[k], sb);
    }
    for (uint32_t i = 0u; i < 12u; ++i) qua_x(ctx->kb, reg[i]);
    qua_dqi_free_b(ctx, &sb, 1u);
    qua_dqi_free_q(ctx, ladder, 11u);
}

static inline void qua_dqi_gf8_mul6_emit(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t qa[3],
    const uint32_t qb[3],
    const uint32_t qout[3]
) {
    uint32_t anc = qua_dqi_alloc_q1(ctx);
    uint32_t hb  = qua_dqi_alloc_b1(ctx);
    for (uint32_t t = 0u; t < 6u; ++t) {
        uint8_t lm = tab->gf8_u[t], rm = tab->gf8_v[t], wm = tab->gf8_w[t];
        uint32_t pa = (uint32_t)__builtin_ctz((unsigned)lm);
        uint32_t pb = (uint32_t)__builtin_ctz((unsigned)rm);
        for (uint32_t i = pa + 1u; i < 3u; ++i) if ((lm >> i) & 1u) qua_cx(ctx->kb, qa[i], qa[pa]);
        for (uint32_t j = pb + 1u; j < 3u; ++j) if ((rm >> j) & 1u) qua_cx(ctx->kb, qb[j], qb[pb]);
        qua_and_compute(ctx->kb, qa[pa], qb[pb], anc);
        for (uint32_t k = 0u; k < 3u; ++k) if ((wm >> k) & 1u) qua_cx(ctx->kb, anc, qout[k]);
        qua_and_uncompute(ctx->kb, qa[pa], qb[pb], anc, hb);
        for (uint32_t j = 2u; j > pb; --j) if ((rm >> j) & 1u) qua_cx(ctx->kb, qb[j], qb[pb]);
        for (uint32_t i = 2u; i > pa; --i) if ((lm >> i) & 1u) qua_cx(ctx->kb, qa[i], qa[pa]);
    }
    qua_dqi_free_b(ctx, &hb, 1u);
    qua_dqi_free_q1(ctx, anc);
}

static inline void qua_dqi_gf64_mul_emit(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t qa[6],
    const uint32_t qb[6],
    const uint32_t qout[6]
) {
    /* Exact Rank-15 Chudnovsky P^1(GF(4)) GF(64) bilinear multiplier (15 CCX vs 18 Karatsuba!) */
    uint32_t anc = qua_dqi_alloc_q1(ctx);
    uint32_t hb  = qua_dqi_alloc_b1(ctx);
    qua_ele_emit_gf64_rank15_mul(ctx->kb, qa, qb, qout, anc, hb, &tab->mul64);
    qua_dqi_free_b(ctx, &hb, 1u);
    qua_dqi_free_q1(ctx, anc);
}

static inline void qua_dqi_phase_by_prod64(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t lhs[6],
    const uint32_t rhs[6],
    const uint32_t mask_bits[6]
) {
    for (uint32_t k = 0u; k < 6u; ++k) {
        for (uint32_t i = 0u; i < 6u; ++i) {
            uint8_t m = tab->gf64_phase[k][i];
            while (m) {
                uint32_t j = (uint32_t)__builtin_ctz((unsigned)m);
                m &= (uint8_t)(m - 1u);
                qua_cz_if(ctx->kb, lhs[i], rhs[j], mask_bits[k]);
            }
        }
    }
}

static inline void qua_dqi_unmul_gf64_hmr(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t tgt[6],
    const uint32_t lhs[6],
    const uint32_t rhs[6]
) {
    uint32_t mb[6];
    qua_dqi_alloc_b(ctx, mb, 6u);
    for (uint32_t k = 0u; k < 6u; ++k) qua_hmr(ctx->kb, tgt[k], mb[k]);
    qua_dqi_phase_by_prod64(ctx, tab, lhs, rhs, mb);
    qua_dqi_free_b(ctx, mb, 6u);
}

static inline void qua_dqi_unmul_gf8_hmr(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t tgt[3],
    const uint32_t lhs[3],
    const uint32_t rhs[3]
) {
    uint32_t mb[3];
    qua_dqi_alloc_b(ctx, mb, 3u);
    for (uint32_t k = 0u; k < 3u; ++k) qua_hmr(ctx->kb, tgt[k], mb[k]);
    for (uint32_t k = 0u; k < 3u; ++k) {
        for (uint32_t i = 0u; i < 3u; ++i) {
            uint8_t m = tab->gf8_phase[k][i];
            while (m) {
                uint32_t j = (uint32_t)__builtin_ctz((unsigned)m);
                m &= (uint8_t)(m - 1u);
                qua_cz_if(ctx->kb, lhs[i], rhs[j], mb[k]);
            }
        }
    }
    qua_dqi_free_b(ctx, mb, 3u);
}

static inline void qua_dqi_erase_d_inv_and_d(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t Delta[6],
    const uint32_t d[3],
    const uint32_t d_inv[3]
) {
    for (uint32_t r = 0u; r < 3u; ++r) {
        for (uint32_t c = 0u; c < 3u; ++c) {
            if ((tab->l_gf8_inv[r] >> c) & 1u) qua_cx(ctx->kb, d[c], d_inv[r]);
        }
    }
    uint16_t q8_inv[12], q8_rows[12];
    memset(q8_rows, 0, sizeof(q8_rows));
    for (uint32_t r = 0u; r < 3u; ++r) q8_rows[r] = tab->q_gf8_inv[r];
    for (uint32_t r = 3u; r < 12u; ++r) q8_rows[r] = (uint16_t)(1u << r);
    qua_dqi_invert_f2_12(q8_rows, q8_inv);
    qua_dqi_emit_inplace_f2(ctx->kb, d_inv, 3u, q8_inv);

    uint32_t mb3[3];
    qua_dqi_alloc_b(ctx, mb3, 3u);
    for (uint32_t i = 0u; i < 3u; ++i) qua_hmr(ctx->kb, d_inv[i], mb3[i]);
    qua_cz_if(ctx->kb, d[0], d[1], mb3[0]);
    qua_cz_if(ctx->kb, d[0], d[2], mb3[1]);
    qua_cz_if(ctx->kb, d[1], d[2], mb3[2]);
    qua_dqi_free_b(ctx, mb3, 3u);

    for (uint32_t r = 0u; r < 3u; ++r) {
        for (uint32_t c = 0u; c < 6u; ++c) {
            if ((tab->m_d_lin[r] >> c) & 1u) qua_cx(ctx->kb, Delta[c], d[r]);
        }
    }
    qua_dqi_unmul_gf8_hmr(ctx, tab, d, Delta, Delta + 3u);
}

static inline void qua_dqi_init_tower_inv_ex(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t x_nat[12],
    const uint32_t target_nat[12],
    uint32_t gate_ctrl, /* UINT32_MAX if ungated */
    uint32_t scaff[21],
    bool pebble_scaff
) {
    if (!pebble_scaff) {
        qua_dqi_emit_inplace_f2(ctx->kb, x_nat, 12u, tab->t_nat2tow);
        const uint32_t *a = x_nat, *b = x_nat + 6u;

        uint32_t *Delta = scaff;
        uint32_t *d = scaff + 6u;
        uint32_t *d_inv = scaff + 9u;
        uint32_t *d_inv_eff = (gate_ctrl != UINT32_MAX) ? (scaff + 12u) : d_inv;
        uint32_t *Delta_inv = (gate_ctrl != UINT32_MAX) ? (scaff + 15u) : (scaff + 12u);
        qua_dqi_alloc_q(ctx, scaff, (gate_ctrl != UINT32_MAX) ? 21u : 18u);

        qua_dqi_gf64_mul_emit(ctx, tab, a, b, Delta);
        for (uint32_t r = 0u; r < 6u; ++r) {
            for (uint32_t c = 0u; c < 12u; ++c) {
                if ((tab->m_delta_lin[r] >> c) & 1u) qua_cx(ctx->kb, x_nat[c], Delta[r]);
            }
        }

        const uint32_t *u0 = Delta, *u1 = Delta + 3u;
        for (uint32_t r = 0u; r < 3u; ++r) {
            for (uint32_t c = 0u; c < 6u; ++c) {
                if ((tab->m_d_lin[r] >> c) & 1u) qua_cx(ctx->kb, Delta[c], d[r]);
            }
        }
        qua_dqi_gf8_mul6_emit(ctx, tab, u0, u1, d);

        qua_ccx(ctx->kb, d[0], d[1], d_inv[0]);
        qua_ccx(ctx->kb, d[0], d[2], d_inv[1]);
        qua_ccx(ctx->kb, d[1], d[2], d_inv[2]);
        uint16_t q8[12];
        memset(q8, 0, sizeof(q8));
        for (uint32_t r = 0u; r < 3u; ++r) q8[r] = tab->q_gf8_inv[r];
        qua_dqi_emit_inplace_f2(ctx->kb, d_inv, 3u, q8);
        for (uint32_t r = 0u; r < 3u; ++r) {
            for (uint32_t c = 0u; c < 3u; ++c) {
                if ((tab->l_gf8_inv[r] >> c) & 1u) qua_cx(ctx->kb, d[c], d_inv[r]);
            }
        }

        if (gate_ctrl != UINT32_MAX) {
            for (uint32_t i = 0u; i < 3u; ++i) qua_and_compute(ctx->kb, gate_ctrl, d_inv[i], d_inv_eff[i]);
        }

        for (uint32_t i = 0u; i < 3u; ++i) qua_cx(ctx->kb, u1[i], u0[i]);
        qua_dqi_gf8_mul6_emit(ctx, tab, u0, d_inv_eff, Delta_inv);
        for (uint32_t i = 0u; i < 3u; ++i) qua_cx(ctx->kb, u1[i], u0[i]);
        qua_dqi_gf8_mul6_emit(ctx, tab, u1, d_inv_eff, Delta_inv + 3u);

        for (uint32_t i = 0u; i < 6u; ++i) qua_cx(ctx->kb, b[i], a[i]);
        qua_dqi_gf64_mul_emit(ctx, tab, a, Delta_inv, target_nat);
        for (uint32_t i = 0u; i < 6u; ++i) qua_cx(ctx->kb, b[i], a[i]);
        qua_dqi_gf64_mul_emit(ctx, tab, b, Delta_inv, target_nat + 6u);

        qua_dqi_emit_inplace_f2(ctx->kb, target_nat, 12u, tab->t_tow2nat);
        qua_dqi_emit_inplace_f2(ctx->kb, x_nat, 12u, tab->t_tow2nat);
        return;
    }

    /* Pebbled Scaffold Mode: holds ONLY d_inv_eff[3] (3 qubits vs 21!) live across EEA round! */
    qua_dqi_emit_inplace_f2(ctx->kb, x_nat, 12u, tab->t_nat2tow);
    const uint32_t *a = x_nat, *b = x_nat + 6u;

    uint32_t *d_inv_eff = scaff;
    qua_dqi_alloc_q(ctx, d_inv_eff, 3u);

    uint32_t tmp12[12];
    qua_dqi_alloc_q(ctx, tmp12, 12u);
    uint32_t *Delta = tmp12;
    uint32_t *d     = tmp12 + 6u;
    uint32_t *d_inv = tmp12 + 9u;
    const uint32_t *u0 = Delta, *u1 = Delta + 3u;

    qua_dqi_gf64_mul_emit(ctx, tab, a, b, Delta);
    for (uint32_t r = 0u; r < 6u; ++r) {
        for (uint32_t c = 0u; c < 12u; ++c) {
            if ((tab->m_delta_lin[r] >> c) & 1u) qua_cx(ctx->kb, x_nat[c], Delta[r]);
        }
    }

    for (uint32_t r = 0u; r < 3u; ++r) {
        for (uint32_t c = 0u; c < 6u; ++c) {
            if ((tab->m_d_lin[r] >> c) & 1u) qua_cx(ctx->kb, Delta[c], d[r]);
        }
    }
    qua_dqi_gf8_mul6_emit(ctx, tab, u0, u1, d);

    qua_ccx(ctx->kb, d[0], d[1], d_inv[0]);
    qua_ccx(ctx->kb, d[0], d[2], d_inv[1]);
    qua_ccx(ctx->kb, d[1], d[2], d_inv[2]);
    uint16_t q8[12];
    memset(q8, 0, sizeof(q8));
    for (uint32_t r = 0u; r < 3u; ++r) q8[r] = tab->q_gf8_inv[r];
    qua_dqi_emit_inplace_f2(ctx->kb, d_inv, 3u, q8);
    for (uint32_t r = 0u; r < 3u; ++r) {
        for (uint32_t c = 0u; c < 3u; ++c) {
            if ((tab->l_gf8_inv[r] >> c) & 1u) qua_cx(ctx->kb, d[c], d_inv[r]);
        }
    }

    if (gate_ctrl != UINT32_MAX) {
        for (uint32_t i = 0u; i < 3u; ++i) qua_and_compute(ctx->kb, gate_ctrl, d_inv[i], d_inv_eff[i]);
    } else {
        for (uint32_t i = 0u; i < 3u; ++i) qua_cx(ctx->kb, d_inv[i], d_inv_eff[i]);
    }

    /* Immediately MBUC-erase d_inv[3] and d[3] (0 CCX!) and recycle those 6 wires for Delta_inv[6]! */
    qua_dqi_erase_d_inv_and_d(ctx, tab, Delta, d, d_inv);
    qua_dqi_free_q(ctx, tmp12 + 6u, 6u);

    uint32_t Delta_inv[6];
    qua_dqi_alloc_q(ctx, Delta_inv, 6u);
    for (uint32_t i = 0u; i < 3u; ++i) qua_cx(ctx->kb, u1[i], u0[i]);
    qua_dqi_gf8_mul6_emit(ctx, tab, u0, d_inv_eff, Delta_inv);
    for (uint32_t i = 0u; i < 3u; ++i) qua_cx(ctx->kb, u1[i], u0[i]);
    qua_dqi_gf8_mul6_emit(ctx, tab, u1, d_inv_eff, Delta_inv + 3u);

    for (uint32_t i = 0u; i < 6u; ++i) qua_cx(ctx->kb, b[i], a[i]);
    qua_dqi_gf64_mul_emit(ctx, tab, a, Delta_inv, target_nat);
    for (uint32_t i = 0u; i < 6u; ++i) qua_cx(ctx->kb, b[i], a[i]);
    qua_dqi_gf64_mul_emit(ctx, tab, b, Delta_inv, target_nat + 6u);

    /* Immediately MBUC-erase Delta_inv[6] and Delta[6] (0 CCX!), leaving only d_inv_eff[3] live! */
    qua_dqi_unmul_gf8_hmr(ctx, tab, Delta_inv + 3u, u1, d_inv_eff);
    for (uint32_t i = 0u; i < 3u; ++i) qua_cx(ctx->kb, u1[i], u0[i]);
    qua_dqi_unmul_gf8_hmr(ctx, tab, Delta_inv, u0, d_inv_eff);
    for (uint32_t i = 0u; i < 3u; ++i) qua_cx(ctx->kb, u1[i], u0[i]);
    qua_dqi_free_q(ctx, Delta_inv, 6u);

    for (uint32_t r = 0u; r < 6u; ++r) {
        for (uint32_t c = 0u; c < 12u; ++c) {
            if ((tab->m_delta_lin[r] >> c) & 1u) qua_cx(ctx->kb, x_nat[c], Delta[r]);
        }
    }
    qua_dqi_unmul_gf64_hmr(ctx, tab, Delta, a, b);
    qua_dqi_free_q(ctx, tmp12, 6u);

    qua_dqi_emit_inplace_f2(ctx->kb, target_nat, 12u, tab->t_tow2nat);
    qua_dqi_emit_inplace_f2(ctx->kb, x_nat, 12u, tab->t_tow2nat);
}

static inline void qua_dqi_init_tower_inv(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t x_nat[12],
    const uint32_t target_nat[12],
    uint32_t gate_ctrl,
    uint32_t scaff[21]
) {
    qua_dqi_init_tower_inv_ex(ctx, tab, x_nat, target_nat, gate_ctrl, scaff, false);
}

static inline void qua_dqi_del_tower_inv_ex(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t x_nat[12],
    const uint32_t target_nat[12],
    uint32_t gate_ctrl,
    const uint32_t scaff[21],
    bool pebble_scaff
) {
    if (!pebble_scaff) {
        qua_dqi_emit_inplace_f2(ctx->kb, x_nat, 12u, tab->t_nat2tow);
        qua_dqi_emit_inplace_f2(ctx->kb, target_nat, 12u, tab->t_nat2tow);
        const uint32_t *a = x_nat, *b = x_nat + 6u;
        const uint32_t *Delta = scaff;
        const uint32_t *d = scaff + 6u;
        const uint32_t *d_inv = scaff + 9u;
        const uint32_t *d_inv_eff = (gate_ctrl != UINT32_MAX) ? (scaff + 12u) : d_inv;
        const uint32_t *Delta_inv = (gate_ctrl != UINT32_MAX) ? (scaff + 15u) : (scaff + 12u);
        const uint32_t *u0 = Delta, *u1 = Delta + 3u;

        qua_dqi_unmul_gf64_hmr(ctx, tab, target_nat + 6u, b, Delta_inv);
        for (uint32_t i = 0u; i < 6u; ++i) qua_cx(ctx->kb, b[i], a[i]);
        qua_dqi_unmul_gf64_hmr(ctx, tab, target_nat, a, Delta_inv);
        for (uint32_t i = 0u; i < 6u; ++i) qua_cx(ctx->kb, b[i], a[i]);

        qua_dqi_unmul_gf8_hmr(ctx, tab, Delta_inv + 3u, u1, d_inv_eff);
        for (uint32_t i = 0u; i < 3u; ++i) qua_cx(ctx->kb, u1[i], u0[i]);
        qua_dqi_unmul_gf8_hmr(ctx, tab, Delta_inv, u0, d_inv_eff);
        for (uint32_t i = 0u; i < 3u; ++i) qua_cx(ctx->kb, u1[i], u0[i]);

        if (gate_ctrl != UINT32_MAX) {
            uint32_t sb = qua_dqi_alloc_b1(ctx);
            for (uint32_t i = 3u; i > 0u; --i) {
                qua_and_uncompute(ctx->kb, gate_ctrl, d_inv[i - 1u], d_inv_eff[i - 1u], sb);
            }
            qua_dqi_free_b(ctx, &sb, 1u);
        }

        qua_dqi_erase_d_inv_and_d(ctx, tab, Delta, d, d_inv);

        for (uint32_t r = 0u; r < 6u; ++r) {
            for (uint32_t c = 0u; c < 12u; ++c) {
                if ((tab->m_delta_lin[r] >> c) & 1u) qua_cx(ctx->kb, x_nat[c], Delta[r]);
            }
        }
        qua_dqi_unmul_gf64_hmr(ctx, tab, Delta, a, b);

        qua_dqi_emit_inplace_f2(ctx->kb, x_nat, 12u, tab->t_tow2nat);
        qua_dqi_free_q(ctx, scaff, (gate_ctrl != UINT32_MAX) ? 21u : 18u);
        return;
    }

    /* Pebbled Scaffold Uncomputation with HMR_Dual target_nat reuse (peak 4 qubits incl. scaff[3]+sb, 36 CCX) */
    qua_dqi_emit_inplace_f2(ctx->kb, x_nat, 12u, tab->t_nat2tow);
    qua_dqi_emit_inplace_f2(ctx->kb, target_nat, 12u, tab->t_nat2tow);

    uint32_t mb_tgt[12];
    qua_dqi_alloc_b(ctx, mb_tgt, 12u);
    for (uint32_t k = 0u; k < 12u; ++k) qua_hmr(ctx->kb, target_nat[k], mb_tgt[k]);

    const uint32_t *a = x_nat, *b = x_nat + 6u;
    const uint32_t *d_inv_eff = scaff;
    const uint32_t *Delta     = target_nat;
    const uint32_t *Delta_inv = target_nat + 6u;
    const uint32_t *u0 = Delta, *u1 = Delta + 3u;

    qua_dqi_gf64_mul_emit(ctx, tab, a, b, Delta);
    for (uint32_t r = 0u; r < 6u; ++r) {
        for (uint32_t c = 0u; c < 12u; ++c) {
            if ((tab->m_delta_lin[r] >> c) & 1u) qua_cx(ctx->kb, x_nat[c], Delta[r]);
        }
    }
    for (uint32_t i = 0u; i < 3u; ++i) qua_cx(ctx->kb, u1[i], u0[i]);
    qua_dqi_gf8_mul6_emit(ctx, tab, u0, d_inv_eff, Delta_inv);
    for (uint32_t i = 0u; i < 3u; ++i) qua_cx(ctx->kb, u1[i], u0[i]);
    qua_dqi_gf8_mul6_emit(ctx, tab, u1, d_inv_eff, Delta_inv + 3u);

    qua_dqi_phase_by_prod64(ctx, tab, b, Delta_inv, mb_tgt + 6u);
    for (uint32_t i = 0u; i < 6u; ++i) qua_cx(ctx->kb, b[i], a[i]);
    qua_dqi_phase_by_prod64(ctx, tab, a, Delta_inv, mb_tgt);
    for (uint32_t i = 0u; i < 6u; ++i) qua_cx(ctx->kb, b[i], a[i]);
    qua_dqi_free_b(ctx, mb_tgt, 12u);

    qua_dqi_unmul_gf8_hmr(ctx, tab, Delta_inv + 3u, u1, d_inv_eff);
    for (uint32_t i = 0u; i < 3u; ++i) qua_cx(ctx->kb, u1[i], u0[i]);
    qua_dqi_unmul_gf8_hmr(ctx, tab, Delta_inv, u0, d_inv_eff);
    for (uint32_t i = 0u; i < 3u; ++i) qua_cx(ctx->kb, u1[i], u0[i]);

    const uint32_t *d     = target_nat + 6u;
    const uint32_t *d_inv = target_nat + 9u;
    for (uint32_t r = 0u; r < 3u; ++r) {
        for (uint32_t c = 0u; c < 6u; ++c) {
            if ((tab->m_d_lin[r] >> c) & 1u) qua_cx(ctx->kb, Delta[c], d[r]);
        }
    }
    qua_dqi_gf8_mul6_emit(ctx, tab, u0, u1, d);
    qua_ccx(ctx->kb, d[0], d[1], d_inv[0]);
    qua_ccx(ctx->kb, d[0], d[2], d_inv[1]);
    qua_ccx(ctx->kb, d[1], d[2], d_inv[2]);
    uint16_t q8[12];
    memset(q8, 0, sizeof(q8));
    for (uint32_t r = 0u; r < 3u; ++r) q8[r] = tab->q_gf8_inv[r];
    qua_dqi_emit_inplace_f2(ctx->kb, d_inv, 3u, q8);
    for (uint32_t r = 0u; r < 3u; ++r) {
        for (uint32_t c = 0u; c < 3u; ++c) {
            if ((tab->l_gf8_inv[r] >> c) & 1u) qua_cx(ctx->kb, d[c], d_inv[r]);
        }
    }

    if (gate_ctrl != UINT32_MAX) {
        uint32_t sb = qua_dqi_alloc_b1(ctx);
        for (uint32_t i = 3u; i > 0u; --i) {
            qua_and_uncompute(ctx->kb, gate_ctrl, d_inv[i - 1u], d_inv_eff[i - 1u], sb);
        }
        qua_dqi_free_b(ctx, &sb, 1u);
    } else {
        for (uint32_t i = 0u; i < 3u; ++i) qua_cx(ctx->kb, d_inv[i], d_inv_eff[i]);
    }

    qua_dqi_erase_d_inv_and_d(ctx, tab, Delta, d, d_inv);

    for (uint32_t r = 0u; r < 6u; ++r) {
        for (uint32_t c = 0u; c < 12u; ++c) {
            if ((tab->m_delta_lin[r] >> c) & 1u) qua_cx(ctx->kb, x_nat[c], Delta[r]);
        }
    }
    qua_dqi_unmul_gf64_hmr(ctx, tab, Delta, a, b);
    qua_dqi_free_q(ctx, scaff, 3u);

    qua_dqi_emit_inplace_f2(ctx->kb, x_nat, 12u, tab->t_tow2nat);
}

static inline void qua_dqi_del_tower_inv(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t x_nat[12],
    const uint32_t target_nat[12],
    uint32_t gate_ctrl,
    const uint32_t scaff[21]
) {
    qua_dqi_del_tower_inv_ex(ctx, tab, x_nat, target_nat, gate_ctrl, scaff, false);
}

/* ============================================================================
 * 5. In-Place Vandermonde Polynomial Slice & MAC over GF(2^12)
 * ============================================================================ */

static inline bool qua_dqi_invert_f2_gen(uint8_t *M, uint8_t *Inv, uint32_t N) {
    /* M and Inv are N x N bit-matrices stored as uint8_t (0 or 1), row-major */
    for (uint32_t r = 0u; r < N; ++r) {
        for (uint32_t c = 0u; c < N; ++c) Inv[r * N + c] = (r == c) ? 1u : 0u;
    }
    for (uint32_t col = 0u; col < N; ++col) {
        uint32_t piv = col;
        while (piv < N && M[piv * N + col] == 0u) ++piv;
        if (piv == N) return false;
        if (piv != col) {
            for (uint32_t c = 0u; c < N; ++c) {
                uint8_t tm = M[col * N + c]; M[col * N + c] = M[piv * N + c]; M[piv * N + c] = tm;
                uint8_t ti = Inv[col * N + c]; Inv[col * N + c] = Inv[piv * N + c]; Inv[piv * N + c] = ti;
            }
        }
        for (uint32_t r = 0u; r < N; ++r) {
            if (r != col && M[r * N + col]) {
                for (uint32_t c = 0u; c < N; ++c) {
                    M[r * N + c] ^= M[col * N + c];
                    Inv[r * N + c] ^= Inv[col * N + c];
                }
            }
        }
    }
    return true;
}

static inline void qua_dqi_emit_poly_slice(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t *const *q_regs, uint32_t len_q,
    const uint32_t *const *c_regs, uint32_t len_c,
    const uint32_t *const *d_regs,
    uint32_t M,
    uint32_t out_lo,
    uint32_t out_hi,
    const uint32_t prod12[12]
) {
    if (len_q * len_c <= M + 1u) {
        for (uint32_t qi = 0u; qi < len_q; ++qi) {
            for (uint32_t ci = 0u; ci < len_c; ++ci) {
                uint32_t deg = qi + ci;
                if (deg >= out_lo && deg < out_hi) {
                    qua_dqi_mul42(ctx, tab, q_regs[qi], c_regs[ci], d_regs[deg - out_lo]);
                }
            }
        }
        return;
    }

    uint32_t dim = 12u * M;
    uint8_t *Vand = (uint8_t *)calloc((size_t)dim * (size_t)dim, 1u);
    uint8_t *Inv  = (uint8_t *)calloc((size_t)dim * (size_t)dim, 1u);
    if (!Vand || !Inv) abort();
    uint32_t deg_q = len_q - 1u;
    uint32_t deg_c = len_c - 1u;
    uint32_t qt_cols = 12u * deg_q;
    uint32_t ct_cols = 12u * deg_c;
    uint8_t *q_tail = (deg_q > 0u) ? (uint8_t *)calloc((size_t)M * 12u * (size_t)qt_cols, 1u) : NULL;
    uint8_t *c_tail = (deg_c > 0u) ? (uint8_t *)calloc((size_t)M * 12u * (size_t)ct_cols, 1u) : NULL;
    if ((deg_q > 0u && !q_tail) || (deg_c > 0u && !c_tail)) abort();

    for (uint32_t p = 0u; p < M; ++p) {
        uint16_t gamma = tab->eval_pts[p];
        uint16_t g_pow = 1u;
        for (uint32_t k = 0u; k < M; ++k) {
            for (uint32_t b = 0u; b < 12u; ++b) {
                uint16_t col_val = qua_gf12_mul_ref(g_pow, (uint16_t)(1u << b), QUA_DQI_MOD_NATIVE);
                for (uint32_t r = 0u; r < 12u; ++r) {
                    uint8_t bit = (uint8_t)((col_val >> r) & 1u);
                    if (bit) {
                        Vand[(12u * p + r) * dim + (12u * k + b)] = 1u;
                        if (k >= 1u && k <= deg_q) {
                            q_tail[(p * 12u + r) * qt_cols + (12u * (k - 1u) + b)] = 1u;
                        }
                        if (k >= 1u && k <= deg_c) {
                            c_tail[(p * 12u + r) * ct_cols + (12u * (k - 1u) + b)] = 1u;
                        }
                    }
                }
            }
            g_pow = qua_gf12_mul_ref(g_pow, gamma, QUA_DQI_MOD_NATIVE);
        }
    }
    qua_dqi_invert_f2_gen(Vand, Inv, dim);

    /* Nearest-neighbor + 2-opt TSP evaluation point order starting & ending at 0 */
    uint32_t dist[QUA_DQI_MAX_PTS][QUA_DQI_MAX_PTS];
    memset(dist, 0, sizeof(dist));
    for (uint32_t a_i = 0u; a_i < M; ++a_i) {
        for (uint32_t b_i = a_i + 1u; b_i < M; ++b_i) {
            uint32_t w = 0u;
            if (q_tail) {
                const uint8_t *a = q_tail + a_i * 12u * qt_cols;
                const uint8_t *b = q_tail + b_i * 12u * qt_cols;
                for (uint32_t i = 0u; i < 12u * qt_cols; ++i) w += (a[i] ^ b[i]);
            }
            if (c_tail) {
                const uint8_t *a = c_tail + a_i * 12u * ct_cols;
                const uint8_t *b = c_tail + b_i * 12u * ct_cols;
                for (uint32_t i = 0u; i < 12u * ct_cols; ++i) w += (a[i] ^ b[i]);
            }
            dist[a_i][b_i] = w;
            dist[b_i][a_i] = w;
        }
    }

    uint32_t order[QUA_DQI_MAX_PTS];
    bool used[QUA_DQI_MAX_PTS];
    memset(used, 0, sizeof(used));
    order[0] = 0u;
    used[0] = true;
    for (uint32_t step = 1u; step < M; ++step) {
        uint32_t prev = order[step - 1u];
        uint32_t best_p = 0u;
        uint32_t best_w = UINT32_MAX;
        for (uint32_t p = 1u; p < M; ++p) {
            if (used[p]) continue;
            if (dist[prev][p] < best_w) { best_w = dist[prev][p]; best_p = p; }
        }
        used[best_p] = true;
        order[step] = best_p;
    }

    /* Closed-tour 2-opt local search (note: p=0 has all-zero tail, so dist[last][0] = return cost!) */
    if (M >= 4u) {
        for (uint32_t pass = 0u; pass < 8u; ++pass) {
            bool improved = false;
            for (uint32_t i = 1u; i + 1u < M; ++i) {
                for (uint32_t k = i + 1u; k < M; ++k) {
                    uint32_t u0 = order[i - 1u], u1 = order[i];
                    uint32_t v0 = order[k], v1 = (k + 1u < M) ? order[k + 1u] : 0u;
                    uint32_t old_d = dist[u0][u1] + dist[v0][v1];
                    uint32_t new_d = dist[u0][v0] + dist[u1][v1];
                    if (new_d < old_d) {
                        uint32_t lo = i, hi = k;
                        while (lo < hi) {
                            uint32_t t = order[lo]; order[lo] = order[hi]; order[hi] = t;
                            ++lo; --hi;
                        }
                        improved = true;
                    }
                }
            }
            if (!improved) break;
        }
    }

    const uint32_t *q0 = q_regs[0], *c0 = c_regs[0];
    uint32_t mbits[12];
    qua_dqi_alloc_b(ctx, mbits, 12u);

    uint8_t *prev_qt = q_tail ? (uint8_t *)calloc(12u * qt_cols, 1u) : NULL;
    uint8_t *prev_ct = c_tail ? (uint8_t *)calloc(12u * ct_cols, 1u) : NULL;
    if ((q_tail && !prev_qt) || (c_tail && !prev_ct)) abort();

    for (uint32_t step = 0u; step < M; ++step) {
        uint32_t p = order[step];
        if (q_tail) {
            const uint8_t *cur_qt = q_tail + p * 12u * qt_cols;
            for (uint32_t r = 0u; r < 12u; ++r) {
                for (uint32_t idx = 0u; idx < qt_cols; ++idx) {
                    if (prev_qt[r * qt_cols + idx] ^ cur_qt[r * qt_cols + idx]) {
                        qua_cx(ctx->kb, q_regs[1u + idx / 12u][idx % 12u], q0[r]);
                    }
                }
            }
            memcpy(prev_qt, cur_qt, 12u * qt_cols);
        }
        if (c_tail) {
            const uint8_t *cur_ct = c_tail + p * 12u * ct_cols;
            for (uint32_t r = 0u; r < 12u; ++r) {
                for (uint32_t idx = 0u; idx < ct_cols; ++idx) {
                    if (prev_ct[r * ct_cols + idx] ^ cur_ct[r * ct_cols + idx]) {
                        qua_cx(ctx->kb, c_regs[1u + idx / 12u][idx % 12u], c0[r]);
                    }
                }
            }
            memcpy(prev_ct, cur_ct, 12u * ct_cols);
        }
        qua_dqi_mul42(ctx, tab, q0, c0, prod12);
        for (uint32_t d_i = out_lo; d_i < out_hi; ++d_i) {
            for (uint32_t r = 0u; r < 12u; ++r) {
                const uint8_t *irow = Inv + (12u * d_i + r) * dim + 12u * p;
                for (uint32_t b = 0u; b < 12u; ++b) {
                    if (irow[b]) qua_cx(ctx->kb, prod12[b], d_regs[d_i - out_lo][r]);
                }
            }
        }
        for (uint32_t k = 0u; k < 12u; ++k) qua_hmr(ctx->kb, prod12[k], mbits[k]);
        qua_dqi_phase_by_prod12(ctx, tab, q0, c0, mbits);
    }

    if (c_tail) {
        for (uint32_t r = 0u; r < 12u; ++r) {
            for (uint32_t idx = 0u; idx < ct_cols; ++idx) {
                if (prev_ct[r * ct_cols + idx]) qua_cx(ctx->kb, c_regs[1u + idx / 12u][idx % 12u], c0[r]);
            }
        }
    }
    if (q_tail) {
        for (uint32_t r = 0u; r < 12u; ++r) {
            for (uint32_t idx = 0u; idx < qt_cols; ++idx) {
                if (prev_qt[r * qt_cols + idx]) qua_cx(ctx->kb, q_regs[1u + idx / 12u][idx % 12u], q0[r]);
            }
        }
    }

    qua_dqi_free_b(ctx, mbits, 12u);
    free(prev_qt);
    free(prev_ct);
    free(q_tail);
    free(c_tail);
    free(Vand);
    free(Inv);
}

static inline void qua_dqi_emit_poly_mac(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t *const *q_regs, uint32_t len_q,
    const uint32_t *const *c_regs, uint32_t len_c,
    const uint32_t *const *d_regs,
    bool top_vanishes,
    const uint32_t prod12[12]
) {
    uint32_t out_len = len_c + len_q - 1u - (top_vanishes ? 1u : 0u);
    qua_dqi_emit_poly_slice(ctx, tab, q_regs, len_q, c_regs, len_c, d_regs, out_len, 0u, out_len, prod12);
}

static inline void qua_dqi_mul_xor_temp_hmr(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t lhs[12],
    const uint32_t rhs[12],
    const uint32_t acc[12],
    const uint32_t inv_reg[12],
    const uint32_t target_q[12]
) {
    uint32_t t_prod[12], mbits[12];
    qua_dqi_alloc_q(ctx, t_prod, 12u);
    qua_dqi_alloc_b(ctx, mbits, 12u);
    qua_dqi_mul42(ctx, tab, lhs, rhs, t_prod);
    for (uint32_t i = 0u; i < 12u; ++i) qua_cx(ctx->kb, t_prod[i], acc[i]);
    qua_dqi_mul42(ctx, tab, acc, inv_reg, target_q);
    for (uint32_t i = 0u; i < 12u; ++i) qua_cx(ctx->kb, t_prod[i], acc[i]);
    for (uint32_t k = 0u; k < 12u; ++k) qua_hmr(ctx->kb, t_prod[k], mbits[k]);
    qua_dqi_phase_by_prod12(ctx, tab, lhs, rhs, mbits);
    qua_dqi_free_b(ctx, mbits, 12u);
    qua_dqi_free_q(ctx, t_prod, 12u);
}

/* Ket-exact 3-mul42 (or 2-mul42 if rhs2==NULL) MBUC double-product quotient helper:
 * target_q ^= (acc ^ lhs1*rhs1 ^ lhs2*rhs2) * inv_reg, saving 84 CCX per call! */
static inline void qua_dqi_mul2_xor_temp_hmr(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t lhs1[12],
    const uint32_t rhs1[12],
    const uint32_t lhs2[12],
    const uint32_t rhs2[12], /* NULL if absent (e.g. r == 1) */
    const uint32_t acc[12],
    const uint32_t inv_reg[12],
    const uint32_t target_q[12]
) {
    uint32_t t_sum[12], mbits[12];
    qua_dqi_alloc_q(ctx, t_sum, 12u);
    qua_dqi_alloc_b(ctx, mbits, 12u);
    qua_dqi_mul42(ctx, tab, lhs1, rhs1, t_sum);
    if (rhs2 != NULL) {
        qua_dqi_mul42(ctx, tab, lhs2, rhs2, t_sum);
    }
    for (uint32_t i = 0u; i < 12u; ++i) qua_cx(ctx->kb, t_sum[i], acc[i]);
    qua_dqi_mul42(ctx, tab, acc, inv_reg, target_q);
    for (uint32_t i = 0u; i < 12u; ++i) qua_cx(ctx->kb, t_sum[i], acc[i]);
    for (uint32_t k = 0u; k < 12u; ++k) qua_hmr(ctx->kb, t_sum[k], mbits[k]);
    qua_dqi_phase_by_prod12(ctx, tab, lhs1, rhs1, mbits);
    if (rhs2 != NULL) {
        qua_dqi_phase_by_prod12(ctx, tab, lhs2, rhs2, mbits);
    }
    qua_dqi_free_b(ctx, mbits, 12u);
    qua_dqi_free_q(ctx, t_sum, 12u);
}

/* HMR_Dual quotient uncomputation: mb_q already measured & freed from target_q.
 * Uncomputes target_q == (acc ^ lhs1*rhs1 ^ (rhs2 ? lhs2*rhs2 : 0)) * inv_reg
 * using 2*42 CCX (or 1*42 if rhs2==NULL) instead of 3*42..5*42. */
static inline void qua_dqi_unquo2_hmr(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t mb_q[12],
    const uint32_t lhs1[12],
    const uint32_t rhs1[12],
    const uint32_t lhs2[12],
    const uint32_t rhs2[12],
    const uint32_t acc[12],
    const uint32_t inv_reg[12]
) {
    uint32_t t_sum[12], mb_t[12];
    qua_dqi_alloc_q(ctx, t_sum, 12u);
    qua_dqi_alloc_b(ctx, mb_t, 12u);
    qua_dqi_mul42(ctx, tab, lhs1, rhs1, t_sum);
    if (rhs2 != NULL) {
        qua_dqi_mul42(ctx, tab, lhs2, rhs2, t_sum);
    }
    for (uint32_t i = 0u; i < 12u; ++i) qua_cx(ctx->kb, t_sum[i], acc[i]);
    qua_dqi_phase_by_prod12(ctx, tab, acc, inv_reg, mb_q);
    for (uint32_t i = 0u; i < 12u; ++i) qua_cx(ctx->kb, t_sum[i], acc[i]);
    for (uint32_t k = 0u; k < 12u; ++k) qua_hmr(ctx->kb, t_sum[k], mb_t[k]);
    qua_dqi_phase_by_prod12(ctx, tab, lhs1, rhs1, mb_t);
    if (rhs2 != NULL) {
        qua_dqi_phase_by_prod12(ctx, tab, lhs2, rhs2, mb_t);
    }
    qua_dqi_free_b(ctx, mb_t, 12u);
    qua_dqi_free_q(ctx, t_sum, 12u);
}

static inline void qua_dqi_unquo1_hmr(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    const uint32_t mb_q[12],
    const uint32_t lhs[12],
    const uint32_t rhs[12],
    const uint32_t acc[12],
    const uint32_t inv_reg[12]
) {
    qua_dqi_unquo2_hmr(ctx, tab, mb_q, lhs, rhs, NULL, NULL, acc, inv_reg);
}

/* ============================================================================
 * 6. Half-EEA + Toeplitz + Chien-Forney Decoder Emission
 * ============================================================================ */

static inline void qua_dqi_toggle_round_flags(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    uint32_t *B_regs[71],
    uint32_t r,
    uint32_t hs,
    const uint32_t *do_div_hi,
    const uint32_t *jump_at
) {
    (void)tab;
    uint32_t lad_lead[11];
    qua_dqi_init_is_zero_12(ctx, B_regs[69u - r], lad_lead);
    uint32_t lead_z = lad_lead[10];

    if (r == 0u) {
        uint32_t lad_next[11];
        qua_dqi_init_is_zero_12(ctx, B_regs[68], lad_next);
        uint32_t next_z = lad_next[10];
        qua_x(ctx->kb, next_z);
        qua_ccx(ctx->kb, lead_z, next_z, jump_at[0]);
        qua_x(ctx->kb, next_z);
        qua_dqi_del_is_zero_12(ctx, B_regs[68], lad_next);
    } else if (r <= hs) {
        if (r == hs) {
            qua_x(ctx->kb, lead_z);
            qua_cx(ctx->kb, lead_z, do_div_hi[0]);
            qua_x(ctx->kb, lead_z);
        }
        uint32_t lad_next[11];
        qua_dqi_init_is_zero_12(ctx, B_regs[68u - r], lad_next);
        uint32_t next_z = lad_next[10];
        uint32_t t1 = qua_dqi_alloc_q1(ctx);
        uint32_t sb = qua_dqi_alloc_b1(ctx);
        qua_x(ctx->kb, next_z);
        qua_x(ctx->kb, jump_at[r - 1u]);
        qua_and_compute(ctx->kb, lead_z, jump_at[r - 1u], t1);
        qua_x(ctx->kb, jump_at[r - 1u]);
        qua_ccx(ctx->kb, t1, next_z, jump_at[r]);
        qua_x(ctx->kb, jump_at[r - 1u]);
        qua_and_uncompute(ctx->kb, lead_z, jump_at[r - 1u], t1, sb);
        qua_x(ctx->kb, jump_at[r - 1u]);
        qua_x(ctx->kb, next_z);
        qua_dqi_free_b(ctx, &sb, 1u);
        qua_dqi_free_q1(ctx, t1);
        qua_dqi_del_is_zero_12(ctx, B_regs[68u - r], lad_next);
    } else {
        uint32_t alive_q = qua_dqi_alloc_q1(ctx);
        qua_cx(ctx->kb, do_div_hi[r - 1u - hs], alive_q);
        qua_cx(ctx->kb, jump_at[r - 1u], alive_q);

        qua_x(ctx->kb, lead_z);
        qua_ccx(ctx->kb, alive_q, lead_z, do_div_hi[r - hs]);
        qua_x(ctx->kb, lead_z);

        if (r < 34u) {
            uint32_t lad_next[11];
            qua_dqi_init_is_zero_12(ctx, B_regs[68u - r], lad_next);
            uint32_t next_z = lad_next[10];
            uint32_t t2[2], sb = qua_dqi_alloc_b1(ctx);
            qua_dqi_alloc_q(ctx, t2, 2u);
            qua_x(ctx->kb, next_z);
            qua_and_compute(ctx->kb, alive_q, lead_z, t2[0]);
            qua_x(ctx->kb, jump_at[r - 1u]);
            qua_and_compute(ctx->kb, jump_at[r - 1u], next_z, t2[1]);
            qua_x(ctx->kb, jump_at[r - 1u]);
            qua_ccx(ctx->kb, t2[0], t2[1], jump_at[r]);
            qua_x(ctx->kb, jump_at[r - 1u]);
            qua_and_uncompute(ctx->kb, jump_at[r - 1u], next_z, t2[1], sb);
            qua_x(ctx->kb, jump_at[r - 1u]);
            qua_and_uncompute(ctx->kb, alive_q, lead_z, t2[0], sb);
            qua_x(ctx->kb, next_z);
            qua_dqi_free_b(ctx, &sb, 1u);
            qua_dqi_free_q(ctx, t2, 2u);
            qua_dqi_del_is_zero_12(ctx, B_regs[68u - r], lad_next);
        }

        qua_cx(ctx->kb, jump_at[r - 1u], alive_q);
        qua_cx(ctx->kb, do_div_hi[r - 1u - hs], alive_q);
        qua_dqi_free_q1(ctx, alive_q);
    }

    qua_dqi_del_is_zero_12(ctx, B_regs[69u - r], lad_lead);
}

static inline void qua_dqi_cswap_slice(
    qua_builder_t *kb,
    uint32_t ctrl,
    uint32_t *const *A,
    uint32_t *const *B,
    uint32_t lo,
    uint32_t hi
) {
    for (uint32_t idx = lo; idx < hi; ++idx) {
        for (uint32_t b = 0u; b < 12u; ++b) {
            qua_cswap(kb, ctrl, A[idx][b], B[idx][b]);
        }
    }
}

static inline void qua_dqi_eea_round_body(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    uint32_t *A_regs[71],
    uint32_t *B_regs[71],
    uint32_t r,
    uint32_t lo_idx,
    uint32_t act_q,
    const uint32_t *jump_at,
    bool is_adjoint,
    bool pebble_scaff,
    bool mbuc_q0
) {
    if (r == 0u) {
        if (!is_adjoint) {
            uint32_t q0[12], q1[12], inv_b[12], sc_b[21];
            qua_dqi_alloc_q(ctx, q0, 12u);
            qua_dqi_alloc_q(ctx, q1, 12u);

            qua_dqi_alloc_q(ctx, inv_b, 12u);
            qua_dqi_init_tower_inv_ex(ctx, tab, B_regs[69], inv_b, act_q, sc_b, pebble_scaff);
            qua_dqi_mul42(ctx, tab, A_regs[70], inv_b, q1);
            qua_dqi_mul_xor_temp_hmr(ctx, tab, q1, B_regs[68], A_regs[69], inv_b, q0);
            qua_dqi_del_tower_inv_ex(ctx, tab, B_regs[69], inv_b, act_q, sc_b, pebble_scaff);
            qua_dqi_free_q(ctx, inv_b, 12u);

            uint32_t prod12[12];
            qua_dqi_alloc_q(ctx, prod12, 12u);
            const uint32_t *qs[2] = {q0, q1};
            qua_dqi_emit_poly_mac(ctx, tab, qs, 2u, (const uint32_t *const *)&B_regs[lo_idx], 70u - lo_idx,
                                  (const uint32_t *const *)&A_regs[lo_idx], false, prod12);
            qua_dqi_free_q(ctx, prod12, 12u);
            qua_dqi_mul42(ctx, tab, q0, B_regs[70], A_regs[70]);
            qua_dqi_mul42(ctx, tab, q1, B_regs[70], A_regs[69]);

            uint32_t mb_q0[12], mb_q1[12], inv_u[12], sc_u[21];
            qua_dqi_alloc_b(ctx, mb_q0, 12u);
            for (uint32_t k = 0u; k < 12u; ++k) qua_hmr(ctx->kb, q0[k], mb_q0[k]);
            qua_dqi_free_q(ctx, q0, 12u);

            qua_dqi_alloc_b(ctx, mb_q1, 12u);
            for (uint32_t k = 0u; k < 12u; ++k) qua_hmr(ctx->kb, q1[k], mb_q1[k]);
            qua_dqi_free_q(ctx, q1, 12u);

            qua_dqi_alloc_q(ctx, inv_u, 12u);
            qua_dqi_init_tower_inv_ex(ctx, tab, B_regs[70], inv_u, act_q, sc_u, pebble_scaff);
            qua_dqi_phase_by_prod12(ctx, tab, A_regs[70], inv_u, mb_q0);
            qua_dqi_phase_by_prod12(ctx, tab, A_regs[69], inv_u, mb_q1);
            qua_dqi_del_tower_inv_ex(ctx, tab, B_regs[70], inv_u, act_q, sc_u, pebble_scaff);
            qua_dqi_free_q(ctx, inv_u, 12u);
            qua_dqi_free_b(ctx, mb_q1, 12u);
            qua_dqi_free_b(ctx, mb_q0, 12u);

            qua_dqi_cswap_slice(ctx->kb, act_q, A_regs, B_regs, lo_idx, 71u);
        } else {
            qua_dqi_cswap_slice(ctx->kb, act_q, A_regs, B_regs, lo_idx, 71u);

            uint32_t q0[12], q1[12], inv_u[12], sc_u[21];
            qua_dqi_alloc_q(ctx, q0, 12u);
            qua_dqi_alloc_q(ctx, q1, 12u);

            qua_dqi_alloc_q(ctx, inv_u, 12u);
            qua_dqi_init_tower_inv_ex(ctx, tab, B_regs[70], inv_u, act_q, sc_u, pebble_scaff);
            qua_dqi_mul42(ctx, tab, A_regs[70], inv_u, q0);
            qua_dqi_mul42(ctx, tab, A_regs[69], inv_u, q1);
            qua_dqi_del_tower_inv_ex(ctx, tab, B_regs[70], inv_u, act_q, sc_u, pebble_scaff);
            qua_dqi_free_q(ctx, inv_u, 12u);

            uint32_t prod12[12];
            qua_dqi_alloc_q(ctx, prod12, 12u);
            const uint32_t *qs[2] = {q0, q1};
            qua_dqi_emit_poly_mac(ctx, tab, qs, 2u, (const uint32_t *const *)&B_regs[lo_idx], 70u - lo_idx,
                                  (const uint32_t *const *)&A_regs[lo_idx], false, prod12);
            qua_dqi_free_q(ctx, prod12, 12u);
            qua_dqi_mul42(ctx, tab, q0, B_regs[70], A_regs[70]);
            qua_dqi_mul42(ctx, tab, q1, B_regs[70], A_regs[69]);

            uint32_t mb_q0[12], inv_b[12], sc_b[21];
            qua_dqi_alloc_b(ctx, mb_q0, 12u);
            for (uint32_t k = 0u; k < 12u; ++k) qua_hmr(ctx->kb, q0[k], mb_q0[k]);
            qua_dqi_free_q(ctx, q0, 12u);

            qua_dqi_alloc_q(ctx, inv_b, 12u);
            qua_dqi_init_tower_inv_ex(ctx, tab, B_regs[69], inv_b, act_q, sc_b, pebble_scaff);
            qua_dqi_unquo2_hmr(ctx, tab, mb_q0, q1, B_regs[68], NULL, NULL, A_regs[69], inv_b);
            qua_dqi_free_b(ctx, mb_q0, 12u);

            qua_dqi_unmul_gf12_hmr(ctx, tab, q1, A_regs[70], inv_b);
            qua_dqi_free_q(ctx, q1, 12u);
            qua_dqi_del_tower_inv_ex(ctx, tab, B_regs[69], inv_b, act_q, sc_b, pebble_scaff);
            qua_dqi_free_q(ctx, inv_b, 12u);
        }
        return;
    }

    uint32_t prev_jump_q = jump_at[r - 1u];
    if (is_adjoint) qua_dqi_cswap_slice(ctx->kb, act_q, A_regs, B_regs, lo_idx, 71u);

    if (!is_adjoint) {
        uint32_t q2[12], inv_b[12], sc_b[21];
        qua_dqi_alloc_q(ctx, inv_b, 12u);
        qua_dqi_init_tower_inv_ex(ctx, tab, B_regs[69u - r], inv_b, act_q, sc_b, pebble_scaff);
        qua_dqi_alloc_q(ctx, q2, 12u);

        uint32_t a_jmp[12], sb = qua_dqi_alloc_b1(ctx);
        qua_dqi_alloc_q(ctx, a_jmp, 12u);
        for (uint32_t i = 0u; i < 12u; ++i) qua_and_compute(ctx->kb, prev_jump_q, A_regs[71u - r][i], a_jmp[i]);
        qua_dqi_mul42(ctx, tab, a_jmp, inv_b, q2);
        for (uint32_t i = 12u; i > 0u; --i) qua_and_uncompute(ctx->kb, prev_jump_q, A_regs[71u - r][i - 1u], a_jmp[i - 1u], sb);
        qua_dqi_free_b(ctx, &sb, 1u);
        qua_dqi_free_q(ctx, a_jmp, 12u);

        uint32_t q1[12], q0[12];
        qua_dqi_alloc_q(ctx, q1, 12u);
        qua_dqi_mul_xor_temp_hmr(ctx, tab, q2, B_regs[68u - r], A_regs[70u - r], inv_b, q1);

        qua_dqi_alloc_q(ctx, q0, 12u);
        if (mbuc_q0) {
            qua_dqi_mul2_xor_temp_hmr(ctx, tab, q1, B_regs[68u - r], q2, B_regs[67u - r], A_regs[69u - r], inv_b, q0);
        } else {
            qua_dqi_mul42(ctx, tab, q2, B_regs[67u - r], A_regs[69u - r]);
            qua_dqi_mul42(ctx, tab, q1, B_regs[68u - r], A_regs[69u - r]);
            qua_dqi_mul42(ctx, tab, A_regs[69u - r], inv_b, q0);
            qua_dqi_mul42(ctx, tab, q1, B_regs[68u - r], A_regs[69u - r]);
            qua_dqi_mul42(ctx, tab, q2, B_regs[67u - r], A_regs[69u - r]);
        }

        qua_dqi_del_tower_inv_ex(ctx, tab, B_regs[69u - r], inv_b, act_q, sc_b, pebble_scaff);
        qua_dqi_free_q(ctx, inv_b, 12u);

        uint32_t prod12[12];
        qua_dqi_alloc_q(ctx, prod12, 12u);
        const uint32_t *qs[3] = {q0, q1, q2};
        qua_dqi_emit_poly_mac(ctx, tab, qs, 3u, (const uint32_t *const *)&B_regs[lo_idx], (70u - r) - lo_idx,
                              (const uint32_t *const *)&A_regs[lo_idx], false, prod12);
        const uint32_t *u_b[40], *u_a[40];
        for (uint32_t j = 0u; j <= r; ++j) u_b[j] = B_regs[70u - j];
        for (uint32_t j = 0u; j <= r + 1u; ++j) u_a[j] = A_regs[70u - j];
        qua_dqi_emit_poly_mac(ctx, tab, qs, 3u, u_b, r + 1u, u_a, true, prod12);
        qua_dqi_free_q(ctx, prod12, 12u);

        uint32_t jmp = qua_dqi_alloc_q1(ctx);
        qua_and_compute(ctx->kb, act_q, prev_jump_q, jmp);

        for (uint32_t b = 0u; b < 12u; ++b) {
            qua_cswap(ctx->kb, jmp, A_regs[69u - r][b], A_regs[70u - r][b]);
            qua_cswap(ctx->kb, jmp, A_regs[70u - r][b], A_regs[71u - r][b]);
        }
        uint32_t *X0 = A_regs[70u - r], *X1 = A_regs[69u - r], *X2_raw = A_regs[71u - r];
        uint32_t *B0, *B1, *B2 = NULL;
        if (r == 1u) {
            for (uint32_t b = 0u; b < 12u; ++b) qua_cswap(ctx->kb, jmp, B_regs[69][b], B_regs[70][b]);
            B0 = B_regs[69]; B1 = B_regs[70];
        } else if (r == 2u) {
            for (uint32_t b = 0u; b < 12u; ++b) {
                qua_cswap(ctx->kb, jmp, B_regs[68][b], B_regs[69][b]);
                qua_cswap(ctx->kb, jmp, B_regs[69][b], B_regs[70][b]);
            }
            B0 = B_regs[68]; B1 = B_regs[69]; B2 = B_regs[70];
        } else {
            for (uint32_t b = 0u; b < 12u; ++b) {
                qua_cswap(ctx->kb, jmp, B_regs[70u - r][b], B_regs[71u - r][b]);
                qua_cswap(ctx->kb, jmp, B_regs[71u - r][b], B_regs[72u - r][b]);
                qua_cswap(ctx->kb, jmp, B_regs[72u - r][b], B_regs[73u - r][b]);
            }
            B0 = B_regs[70u - r]; B1 = B_regs[71u - r]; B2 = B_regs[72u - r];
        }

        uint32_t mb_q0[12];
        qua_dqi_alloc_b(ctx, mb_q0, 12u);
        for (uint32_t k = 0u; k < 12u; ++k) qua_hmr(ctx->kb, q0[k], mb_q0[k]);
        qua_dqi_free_q(ctx, q0, 12u);

        uint32_t inv_u[12], sc_u[21];
        qua_dqi_alloc_q(ctx, inv_u, 12u);
        qua_dqi_init_tower_inv_ex(ctx, tab, B0, inv_u, act_q, sc_u, pebble_scaff);

        qua_dqi_unquo2_hmr(ctx, tab, mb_q0, q1, B1, q2, B2, X0, inv_u);
        qua_dqi_free_b(ctx, mb_q0, 12u);

        uint32_t mb_q1[12];
        qua_dqi_alloc_b(ctx, mb_q1, 12u);
        for (uint32_t k = 0u; k < 12u; ++k) qua_hmr(ctx->kb, q1[k], mb_q1[k]);
        qua_dqi_free_q(ctx, q1, 12u);
        qua_dqi_unquo1_hmr(ctx, tab, mb_q1, q2, B1, X1, inv_u);
        qua_dqi_free_b(ctx, mb_q1, 12u);

        uint32_t mb_q2[12], x2_eff[12];
        qua_dqi_alloc_b(ctx, mb_q2, 12u);
        for (uint32_t k = 0u; k < 12u; ++k) qua_hmr(ctx->kb, q2[k], mb_q2[k]);
        qua_dqi_free_q(ctx, q2, 12u);

        sb = qua_dqi_alloc_b1(ctx);
        qua_dqi_alloc_q(ctx, x2_eff, 12u);
        for (uint32_t i = 0u; i < 12u; ++i) qua_and_compute(ctx->kb, jmp, X2_raw[i], x2_eff[i]);
        qua_dqi_phase_by_prod12(ctx, tab, x2_eff, inv_u, mb_q2);
        for (uint32_t i = 12u; i > 0u; --i) qua_and_uncompute(ctx->kb, jmp, X2_raw[i - 1u], x2_eff[i - 1u], sb);
        qua_dqi_free_b(ctx, &sb, 1u);
        qua_dqi_free_q(ctx, x2_eff, 12u);
        qua_dqi_free_b(ctx, mb_q2, 12u);

        qua_dqi_del_tower_inv_ex(ctx, tab, B0, inv_u, act_q, sc_u, pebble_scaff);
        qua_dqi_free_q(ctx, inv_u, 12u);

        if (r == 1u) {
            for (uint32_t b = 0u; b < 12u; ++b) qua_cswap(ctx->kb, jmp, B_regs[69][b], B_regs[70][b]);
        } else if (r == 2u) {
            for (uint32_t b = 0u; b < 12u; ++b) {
                qua_cswap(ctx->kb, jmp, B_regs[69][b], B_regs[70][b]);
                qua_cswap(ctx->kb, jmp, B_regs[68][b], B_regs[69][b]);
            }
        } else {
            for (uint32_t b = 0u; b < 12u; ++b) {
                qua_cswap(ctx->kb, jmp, B_regs[72u - r][b], B_regs[73u - r][b]);
                qua_cswap(ctx->kb, jmp, B_regs[71u - r][b], B_regs[72u - r][b]);
                qua_cswap(ctx->kb, jmp, B_regs[70u - r][b], B_regs[71u - r][b]);
            }
        }
        for (uint32_t b = 0u; b < 12u; ++b) {
            qua_cswap(ctx->kb, jmp, A_regs[70u - r][b], A_regs[71u - r][b]);
            qua_cswap(ctx->kb, jmp, A_regs[69u - r][b], A_regs[70u - r][b]);
        }
        sb = qua_dqi_alloc_b1(ctx);
        qua_and_uncompute(ctx->kb, act_q, prev_jump_q, jmp, sb);
        qua_dqi_free_b(ctx, &sb, 1u);
        qua_dqi_free_q1(ctx, jmp);

        qua_dqi_cswap_slice(ctx->kb, act_q, A_regs, B_regs, lo_idx, 71u);
    } else {
        uint32_t jmp = qua_dqi_alloc_q1(ctx);
        qua_and_compute(ctx->kb, act_q, prev_jump_q, jmp);

        for (uint32_t b = 0u; b < 12u; ++b) {
            qua_cswap(ctx->kb, jmp, A_regs[69u - r][b], A_regs[70u - r][b]);
            qua_cswap(ctx->kb, jmp, A_regs[70u - r][b], A_regs[71u - r][b]);
        }
        uint32_t *X0 = A_regs[70u - r], *X1 = A_regs[69u - r], *X2_raw = A_regs[71u - r];
        uint32_t *B0, *B1, *B2 = NULL;
        if (r == 1u) {
            for (uint32_t b = 0u; b < 12u; ++b) qua_cswap(ctx->kb, jmp, B_regs[69][b], B_regs[70][b]);
            B0 = B_regs[69]; B1 = B_regs[70];
        } else if (r == 2u) {
            for (uint32_t b = 0u; b < 12u; ++b) {
                qua_cswap(ctx->kb, jmp, B_regs[68][b], B_regs[69][b]);
                qua_cswap(ctx->kb, jmp, B_regs[69][b], B_regs[70][b]);
            }
            B0 = B_regs[68]; B1 = B_regs[69]; B2 = B_regs[70];
        } else {
            for (uint32_t b = 0u; b < 12u; ++b) {
                qua_cswap(ctx->kb, jmp, B_regs[70u - r][b], B_regs[71u - r][b]);
                qua_cswap(ctx->kb, jmp, B_regs[71u - r][b], B_regs[72u - r][b]);
                qua_cswap(ctx->kb, jmp, B_regs[72u - r][b], B_regs[73u - r][b]);
            }
            B0 = B_regs[70u - r]; B1 = B_regs[71u - r]; B2 = B_regs[72u - r];
        }

        uint32_t inv_u[12], sc_u[21];
        qua_dqi_alloc_q(ctx, inv_u, 12u);
        qua_dqi_init_tower_inv_ex(ctx, tab, B0, inv_u, act_q, sc_u, pebble_scaff);

        uint32_t q2[12], x2_eff[12], sb = qua_dqi_alloc_b1(ctx);
        qua_dqi_alloc_q(ctx, q2, 12u);
        qua_dqi_alloc_q(ctx, x2_eff, 12u);
        for (uint32_t i = 0u; i < 12u; ++i) qua_and_compute(ctx->kb, jmp, X2_raw[i], x2_eff[i]);
        qua_dqi_mul42(ctx, tab, x2_eff, inv_u, q2);
        for (uint32_t i = 12u; i > 0u; --i) qua_and_uncompute(ctx->kb, jmp, X2_raw[i - 1u], x2_eff[i - 1u], sb);
        qua_dqi_free_b(ctx, &sb, 1u);
        qua_dqi_free_q(ctx, x2_eff, 12u);

        uint32_t q1[12], q0[12];
        qua_dqi_alloc_q(ctx, q1, 12u);
        qua_dqi_mul_xor_temp_hmr(ctx, tab, q2, B1, X1, inv_u, q1);

        qua_dqi_alloc_q(ctx, q0, 12u);
        if (mbuc_q0) {
            qua_dqi_mul2_xor_temp_hmr(ctx, tab, q1, B1, q2, B2, X0, inv_u, q0);
        } else {
            qua_dqi_mul42(ctx, tab, q1, B1, X0);
            if (B2) qua_dqi_mul42(ctx, tab, q2, B2, X0);
            qua_dqi_mul42(ctx, tab, X0, inv_u, q0);
            if (B2) qua_dqi_mul42(ctx, tab, q2, B2, X0);
            qua_dqi_mul42(ctx, tab, q1, B1, X0);
        }

        qua_dqi_del_tower_inv_ex(ctx, tab, B0, inv_u, act_q, sc_u, pebble_scaff);
        qua_dqi_free_q(ctx, inv_u, 12u);

        if (r == 1u) {
            for (uint32_t b = 0u; b < 12u; ++b) qua_cswap(ctx->kb, jmp, B_regs[69][b], B_regs[70][b]);
        } else if (r == 2u) {
            for (uint32_t b = 0u; b < 12u; ++b) {
                qua_cswap(ctx->kb, jmp, B_regs[69][b], B_regs[70][b]);
                qua_cswap(ctx->kb, jmp, B_regs[68][b], B_regs[69][b]);
            }
        } else {
            for (uint32_t b = 0u; b < 12u; ++b) {
                qua_cswap(ctx->kb, jmp, B_regs[72u - r][b], B_regs[73u - r][b]);
                qua_cswap(ctx->kb, jmp, B_regs[71u - r][b], B_regs[72u - r][b]);
                qua_cswap(ctx->kb, jmp, B_regs[70u - r][b], B_regs[71u - r][b]);
            }
        }
        for (uint32_t b = 0u; b < 12u; ++b) {
            qua_cswap(ctx->kb, jmp, A_regs[70u - r][b], A_regs[71u - r][b]);
            qua_cswap(ctx->kb, jmp, A_regs[69u - r][b], A_regs[70u - r][b]);
        }
        sb = qua_dqi_alloc_b1(ctx);
        qua_and_uncompute(ctx->kb, act_q, prev_jump_q, jmp, sb);
        qua_dqi_free_b(ctx, &sb, 1u);
        qua_dqi_free_q1(ctx, jmp);

        uint32_t prod12[12];
        qua_dqi_alloc_q(ctx, prod12, 12u);
        const uint32_t *qs[3] = {q0, q1, q2};
        qua_dqi_emit_poly_mac(ctx, tab, qs, 3u, (const uint32_t *const *)&B_regs[lo_idx], (70u - r) - lo_idx,
                              (const uint32_t *const *)&A_regs[lo_idx], false, prod12);
        const uint32_t *u_b[40], *u_a[40];
        for (uint32_t j = 0u; j <= r; ++j) u_b[j] = B_regs[70u - j];
        for (uint32_t j = 0u; j <= r + 1u; ++j) u_a[j] = A_regs[70u - j];
        qua_dqi_emit_poly_mac(ctx, tab, qs, 3u, u_b, r + 1u, u_a, true, prod12);
        qua_dqi_free_q(ctx, prod12, 12u);

        uint32_t mb_q0[12];
        qua_dqi_alloc_b(ctx, mb_q0, 12u);
        for (uint32_t k = 0u; k < 12u; ++k) qua_hmr(ctx->kb, q0[k], mb_q0[k]);
        qua_dqi_free_q(ctx, q0, 12u);

        uint32_t inv_b[12], sc_b[21];
        qua_dqi_alloc_q(ctx, inv_b, 12u);
        qua_dqi_init_tower_inv_ex(ctx, tab, B_regs[69u - r], inv_b, act_q, sc_b, pebble_scaff);

        qua_dqi_unquo2_hmr(ctx, tab, mb_q0, q1, B_regs[68u - r], q2, B_regs[67u - r], A_regs[69u - r], inv_b);
        qua_dqi_free_b(ctx, mb_q0, 12u);

        uint32_t mb_q1[12];
        qua_dqi_alloc_b(ctx, mb_q1, 12u);
        for (uint32_t k = 0u; k < 12u; ++k) qua_hmr(ctx->kb, q1[k], mb_q1[k]);
        qua_dqi_free_q(ctx, q1, 12u);
        qua_dqi_unquo1_hmr(ctx, tab, mb_q1, q2, B_regs[68u - r], A_regs[70u - r], inv_b);
        qua_dqi_free_b(ctx, mb_q1, 12u);

        uint32_t mb_q2[12], a_jmp[12];
        qua_dqi_alloc_b(ctx, mb_q2, 12u);
        for (uint32_t k = 0u; k < 12u; ++k) qua_hmr(ctx->kb, q2[k], mb_q2[k]);
        qua_dqi_free_q(ctx, q2, 12u);

        sb = qua_dqi_alloc_b1(ctx);
        qua_dqi_alloc_q(ctx, a_jmp, 12u);
        for (uint32_t i = 0u; i < 12u; ++i) qua_and_compute(ctx->kb, prev_jump_q, A_regs[71u - r][i], a_jmp[i]);
        qua_dqi_phase_by_prod12(ctx, tab, a_jmp, inv_b, mb_q2);
        for (uint32_t i = 12u; i > 0u; --i) qua_and_uncompute(ctx->kb, prev_jump_q, A_regs[71u - r][i - 1u], a_jmp[i - 1u], sb);
        qua_dqi_free_b(ctx, &sb, 1u);
        qua_dqi_free_q(ctx, a_jmp, 12u);
        qua_dqi_free_b(ctx, mb_q2, 12u);

        qua_dqi_del_tower_inv_ex(ctx, tab, B_regs[69u - r], inv_b, act_q, sc_b, pebble_scaff);
        qua_dqi_free_q(ctx, inv_b, 12u);
    }
}

static inline void qua_dqi_toeplitz_rec(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    uint32_t *const *s_slice,
    uint32_t m,
    uint32_t t_base,
    const uint32_t sigma0[12],
    const uint32_t inv_s0[12],
    uint32_t s0_nz,
    uint32_t *const sigma_off[34],
    const uint32_t t_work[12],
    bool rev
) {
    if (m <= t_base) {
        for (uint32_t ii = 0u; ii < m; ++ii) {
            uint32_t k = rev ? ii : (m - 1u - ii);
            if (!rev) {
                qua_dqi_mul42(ctx, tab, sigma0, s_slice[k], t_work);
                for (uint32_t b = 0u; b < 12u; ++b) qua_cswap(ctx->kb, s0_nz, s_slice[k][b], t_work[b]);
                qua_dqi_unmul_gf12_hmr(ctx, tab, t_work, inv_s0, s_slice[k]);
                for (uint32_t j = 1u; j <= k; ++j) {
                    qua_dqi_mul42(ctx, tab, sigma_off[j - 1u], s_slice[k - j], s_slice[k]);
                }
            } else {
                for (uint32_t j = k; j >= 1u; --j) {
                    qua_dqi_mul42(ctx, tab, sigma_off[j - 1u], s_slice[k - j], s_slice[k]);
                }
                qua_dqi_mul42(ctx, tab, inv_s0, s_slice[k], t_work);
                for (uint32_t b = 0u; b < 12u; ++b) qua_cswap(ctx->kb, s0_nz, s_slice[k][b], t_work[b]);
                qua_dqi_unmul_gf12_hmr(ctx, tab, t_work, sigma0, s_slice[k]);
            }
        }
        return;
    }
    uint32_t h = m / 2u;
    if (!rev) {
        qua_dqi_toeplitz_rec(ctx, tab, s_slice + h, m - h, t_base, sigma0, inv_s0, s0_nz, sigma_off, t_work, false);
        qua_dqi_emit_poly_slice(ctx, tab, (const uint32_t *const *)sigma_off, m - 1u,
                                (const uint32_t *const *)s_slice, h, (const uint32_t *const *)(s_slice + h),
                                (m - 2u) + h, h - 1u, m - 1u, t_work);
        qua_dqi_toeplitz_rec(ctx, tab, s_slice, h, t_base, sigma0, inv_s0, s0_nz, sigma_off, t_work, false);
    } else {
        qua_dqi_toeplitz_rec(ctx, tab, s_slice, h, t_base, sigma0, inv_s0, s0_nz, sigma_off, t_work, true);
        qua_dqi_emit_poly_slice(ctx, tab, (const uint32_t *const *)sigma_off, m - 1u,
                                (const uint32_t *const *)s_slice, h, (const uint32_t *const *)(s_slice + h),
                                (m - 2u) + h, h - 1u, m - 1u, t_work);
        qua_dqi_toeplitz_rec(ctx, tab, s_slice + h, m - h, t_base, sigma0, inv_s0, s0_nz, sigma_off, t_work, true);
    }
}

static inline void qua_dqi_chien_forney(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    uint32_t *B[71],
    uint32_t mbuc_bits[QUA_DQI_M][12]
) {
    uint32_t *sigma_eval = B[70];
    uint32_t *sigma_prime_x = B[69];
    uint32_t *omega_x = B[0];

    uint32_t sigma_prime_inv[12];
    qua_dqi_alloc_q(ctx, sigma_prime_inv, 12u);
    uint32_t is_error = qua_dqi_alloc_b1(ctx);

    uint16_t mult[70];
    uint32_t mult_idx[70];
    uint32_t n_mult = 0u;
    for (uint32_t k = 1u; k <= 35u; ++k) {
        mult_idx[n_mult] = 70u - k;
        mult[n_mult++] = qua_dqi_gf12_pow(tab, 2u, (4095u - k) % 4095u);
    }
    for (uint32_t k = 0u; k < 35u; ++k) {
        mult_idx[n_mult] = k;
        mult[n_mult++] = qua_dqi_gf12_pow(tab, 2u, (4095u - (k + 1u)) % 4095u);
    }

    uint32_t *sum_src[70], *sum_dst[70];
    uint32_t n_sum = 0u;
    for (uint32_t k = 1u; k < 35u; ++k) { sum_src[n_sum] = B[k]; sum_dst[n_sum++] = omega_x; }
    for (uint32_t k = 3u; k < 36u; k += 2u) { sum_src[n_sum] = B[70u - k]; sum_dst[n_sum++] = sigma_prime_x; }
    for (uint32_t k = 2u; k < 36u; k += 2u) { sum_src[n_sum] = B[70u - k]; sum_dst[n_sum++] = sigma_eval; }
    sum_src[n_sum] = sigma_prime_x; sum_dst[n_sum++] = sigma_eval;

    /* Precompute 70 constant-multiplication F2 gate sequences once (avoids 286,510 Gaussian eliminations!) */
    typedef struct { uint8_t kind, a, b; } cf_op_t;
    typedef struct { cf_op_t ops[144]; uint8_t nops; } cf_seq_t;
    cf_seq_t *mult_seqs = (cf_seq_t *)calloc(n_mult, sizeof(cf_seq_t));
    if (!mult_seqs) abort();
    for (uint32_t m_i = 0u; m_i < n_mult; ++m_i) {
        uint16_t c_val = mult[m_i];
        if (c_val == 1u) continue;
        uint16_t A[12];
        memset(A, 0, sizeof(A));
        for (uint32_t b = 0u; b < 12u; ++b) {
            uint16_t p = qua_gf12_mul_ref(c_val, (uint16_t)(1u << b), QUA_DQI_MOD_NATIVE);
            for (uint32_t r = 0u; r < 12u; ++r) {
                if ((p >> r) & 1u) A[r] |= (uint16_t)(1u << b);
            }
        }
        uint8_t nops = 0u;
        for (uint32_t col = 0u; col < 12u; ++col) {
            uint32_t best_piv = col, best_cost = UINT32_MAX;
            for (uint32_t r = col; r < 12u; ++r) {
                if (((A[r] >> col) & 1u) == 0u) continue;
                uint32_t cost = (uint32_t)__builtin_popcount((unsigned)A[r]) + ((r != col) ? 2u : 0u);
                if (cost < best_cost) { best_cost = cost; best_piv = r; }
            }
            if (best_piv != col) {
                if (nops >= 144u) abort();
                uint16_t tmp = A[col]; A[col] = A[best_piv]; A[best_piv] = tmp;
                mult_seqs[m_i].ops[nops++] = (cf_op_t){1u, (uint8_t)col, (uint8_t)best_piv};
            }
            for (uint32_t r = 0u; r < 12u; ++r) {
                if (r != col && ((A[r] >> col) & 1u)) {
                    if (nops >= 144u) abort();
                    A[r] ^= A[col];
                    mult_seqs[m_i].ops[nops++] = (cf_op_t){0u, (uint8_t)r, (uint8_t)col};
                }
            }
        }
        mult_seqs[m_i].nops = nops;
    }

    for (uint32_t j = 0u; j < QUA_DQI_M; ++j) {
        for (uint32_t b = 0u; b < 12u; ++b) qua_bit_store1_if(ctx->kb, is_error, mbuc_bits[j][b]);
        qua_push_condition(ctx->kb, is_error);

        for (uint32_t p = 0u; p < n_sum; ++p) {
            for (uint32_t b = 0u; b < 12u; ++b) qua_cx(ctx->kb, sum_src[p][b], sum_dst[p][b]);
        }
        uint32_t zero_lad[11];
        qua_dqi_init_is_zero_12(ctx, sigma_eval, zero_lad);
        uint32_t is_sigma_zero = zero_lad[10];

        uint32_t inv_sc[21];
        qua_dqi_init_tower_inv(ctx, tab, sigma_prime_x, sigma_prime_inv, is_sigma_zero, inv_sc);
        qua_dqi_phase_by_prod12(ctx, tab, omega_x, sigma_prime_inv, mbuc_bits[j]);
        qua_dqi_del_tower_inv(ctx, tab, sigma_prime_x, sigma_prime_inv, is_sigma_zero, inv_sc);
        qua_dqi_del_is_zero_12(ctx, sigma_eval, zero_lad);

        for (uint32_t p = n_sum; p > 0u; --p) {
            for (uint32_t b = 0u; b < 12u; ++b) qua_cx(ctx->kb, sum_src[p - 1u][b], sum_dst[p - 1u][b]);
        }

        qua_pop_condition(ctx->kb);
        qua_bit_store0(ctx->kb, is_error);

        if (j + 1u < QUA_DQI_M) {
            for (uint32_t m_i = 0u; m_i < n_mult; ++m_i) {
                const uint32_t *qreg = B[mult_idx[m_i]];
                const cf_seq_t *sq = &mult_seqs[m_i];
                for (uint32_t k = sq->nops; k > 0u; --k) {
                    cf_op_t op = sq->ops[k - 1u];
                    if (op.kind == 1u) {
                        qua_cx(ctx->kb, qreg[op.a], qreg[op.b]);
                        qua_cx(ctx->kb, qreg[op.b], qreg[op.a]);
                        qua_cx(ctx->kb, qreg[op.a], qreg[op.b]);
                    } else {
                        qua_cx(ctx->kb, qreg[op.b], qreg[op.a]);
                    }
                }
            }
        }
    }
    free(mult_seqs);

    for (uint32_t m_i = 0u; m_i < n_mult; ++m_i) {
        uint16_t inv_m = qua_dqi_gf12_pow(tab, qua_dqi_gf12_pow(tab, mult[m_i], QUA_DQI_M - 1u), 4094u);
        qua_dqi_gf12_imul_const(ctx->kb, B[mult_idx[m_i]], inv_m);
    }

    qua_dqi_free_b(ctx, &is_error, 1u);
    qua_dqi_free_q(ctx, sigma_prime_inv, 12u);
}

static inline void qua_dqi_emit_hgcd_checkpoint(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    uint32_t *const A_regs[71],
    uint32_t *const B_regs[71],
    uint32_t syn_q[][12],
    const uint32_t jump_at[34],
    uint32_t r_split,
    uint32_t syn_off,
    uint32_t syn_len,
    uint32_t dst_lo
) {
    uint32_t prod12[12];
    qua_dqi_alloc_q(ctx, prod12, 12u);
    uint32_t slice_lo = dst_lo - syn_off;
    uint32_t out_hi = r_split + syn_len;
    const uint32_t *ub[40], *b_dst[60], *s_seg[40];
    for (uint32_t j = 0u; j <= r_split; ++j) ub[j] = B_regs[70u - j];
    for (uint32_t i = 0u; i < syn_len; ++i) s_seg[i] = syn_q[syn_off + i];
    for (uint32_t d = slice_lo; d < out_hi; ++d) b_dst[d - slice_lo] = B_regs[syn_off + d];
    qua_dqi_emit_poly_slice(ctx, tab, ub, r_split + 1u, s_seg, syn_len, b_dst, out_hi, slice_lo, out_hi, prod12);

    uint32_t ua_last[12], sb = qua_dqi_alloc_b1(ctx);
    qua_dqi_alloc_q(ctx, ua_last, 12u);
    qua_x(ctx->kb, jump_at[r_split - 1u]);
    for (uint32_t i = 0u; i < 12u; ++i) qua_and_compute(ctx->kb, jump_at[r_split - 1u], A_regs[71u - r_split][i], ua_last[i]);
    qua_x(ctx->kb, jump_at[r_split - 1u]);

    const uint32_t *ua[40], *a_dst[60];
    uint32_t a_out_hi = r_split + syn_len - 1u;
    for (uint32_t j = 0u; j < r_split - 1u; ++j) ua[j] = A_regs[70u - j];
    ua[r_split - 1u] = ua_last;
    for (uint32_t d = slice_lo; d < a_out_hi; ++d) a_dst[d - slice_lo] = A_regs[syn_off + d];
    qua_dqi_emit_poly_slice(ctx, tab, ua, r_split, s_seg, syn_len, a_dst, a_out_hi, slice_lo, a_out_hi, prod12);

    qua_x(ctx->kb, jump_at[r_split - 1u]);
    for (uint32_t i = 12u; i > 0u; --i) qua_and_uncompute(ctx->kb, jump_at[r_split - 1u], A_regs[71u - r_split][i - 1u], ua_last[i - 1u], sb);
    qua_x(ctx->kb, jump_at[r_split - 1u]);
    qua_dqi_free_b(ctx, &sb, 1u);
    qua_dqi_free_q(ctx, ua_last, 12u);
    qua_dqi_free_q(ctx, prod12, 12u);
}

static inline void qua_dqi_emit_hgcd_phase(
    qua_dqi_ctx_t *ctx,
    const qua_dqi_tables_t *tab,
    uint32_t *A_regs[71],
    uint32_t *B_regs[71],
    uint32_t r_begin,
    uint32_t r_end,
    uint32_t lo_idx,
    bool use_lo_r,
    uint32_t k_min,
    uint32_t hs,
    const uint32_t *do_div_hi,
    const uint32_t *jump_at,
    bool is_adjoint,
    bool pebble_scaff,
    bool mbuc_q0
) {
    if (!is_adjoint) {
        for (uint32_t r = r_begin; r < r_end; ++r) {
            uint32_t lo_r = use_lo_r ? ((k_min > (r - 1u)) ? k_min : (r - 1u)) : lo_idx;
            qua_dqi_toggle_round_flags(ctx, tab, B_regs, r, hs, do_div_hi, jump_at);
            if (r < hs) {
                qua_x(ctx->kb, jump_at[r]);
                qua_dqi_eea_round_body(ctx, tab, A_regs, B_regs, r, lo_r, jump_at[r], jump_at, false, pebble_scaff, mbuc_q0);
                qua_x(ctx->kb, jump_at[r]);
            } else {
                qua_dqi_eea_round_body(ctx, tab, A_regs, B_regs, r, lo_r, do_div_hi[r - hs], jump_at, false, pebble_scaff, mbuc_q0);
            }
        }
    } else {
        for (uint32_t rr = r_end; rr > r_begin; --rr) {
            uint32_t r = rr - 1u;
            uint32_t lo_r = use_lo_r ? ((k_min > (r - 1u)) ? k_min : (r - 1u)) : lo_idx;
            if (r < hs) {
                qua_x(ctx->kb, jump_at[r]);
                qua_dqi_eea_round_body(ctx, tab, A_regs, B_regs, r, lo_r, jump_at[r], jump_at, true, pebble_scaff, mbuc_q0);
                qua_x(ctx->kb, jump_at[r]);
            } else {
                qua_dqi_eea_round_body(ctx, tab, A_regs, B_regs, r, lo_r, do_div_hi[r - hs], jump_at, true, pebble_scaff, mbuc_q0);
            }
            qua_dqi_toggle_round_flags(ctx, tab, B_regs, r, hs, do_div_hi, jump_at);
        }
    }
}

static inline int qua_dqi_emit_kmb(
    const qua_dqi_tables_t *tab,
    const qua_dqi_cfg_t *cfg,
    const char *out_path
) {
    qua_builder_t kb;
    qua_builder_init(&kb);

    static uint32_t syn_q[QUA_DQI_N][12];
    static uint32_t mbuc_b[QUA_DQI_M][12];
    for (uint32_t i = 0u; i < QUA_DQI_N; ++i) {
        for (uint32_t b = 0u; b < 12u; ++b) syn_q[i][b] = i * 12u + b;
        char rname[16];
        snprintf(rname, sizeof(rname), "S%u", i);
        qua_add_quantum_reg(&kb, rname, syn_q[i], 12u);
    }
    for (uint32_t j = 0u; j < QUA_DQI_M; ++j) {
        for (uint32_t b = 0u; b < 12u; ++b) mbuc_b[j][b] = j * 12u + b;
        char rname[16];
        snprintf(rname, sizeof(rname), "E%u", j);
        qua_add_classical_reg(&kb, rname, mbuc_b[j], 12u);
    }

    qua_dqi_ctx_t ctx;
    qua_dqi_ctx_init(&ctx, &kb, QUA_DQI_N * 12u, QUA_DQI_M * 12u);

    uint32_t k_min = cfg->k_min, s1 = cfg->s1, r1 = cfg->r1, hs = cfg->halt_start;
    static uint32_t B_pool[40][12], A_pool[71][12];
    uint32_t b70[12];
    uint32_t *B_regs[71], *A_regs[71];
    memset(B_regs, 0, sizeof(B_regs));
    memset(A_regs, 0, sizeof(A_regs));

    for (uint32_t idx = k_min; idx < 35u; ++idx) {
        qua_dqi_alloc_q(&ctx, B_pool[idx - k_min], 12u);
        B_regs[idx] = B_pool[idx - k_min];
    }
    for (uint32_t idx = 35u; idx < 70u; ++idx) B_regs[idx] = syn_q[idx];
    qua_dqi_alloc_q(&ctx, b70, 12u);
    B_regs[70] = b70;

    for (uint32_t idx = k_min; idx <= 70u; ++idx) {
        qua_dqi_alloc_q(&ctx, A_pool[idx - k_min], 12u);
        A_regs[idx] = A_pool[idx - k_min];
    }

    uint32_t do_div_hi[35], jump_at[34];
    qua_dqi_alloc_q(&ctx, do_div_hi, 35u - hs);
    qua_dqi_alloc_q(&ctx, jump_at, 34u);

    uint32_t init_s = cfg->three_stage ? cfg->s0 : s1;
    for (uint32_t i = init_s; i < 35u; ++i) {
        for (uint32_t b = 0u; b < 12u; ++b) qua_cx(&kb, syn_q[i][b], B_regs[i][b]);
    }
    qua_x(&kb, A_regs[70][0]);
    qua_x(&kb, B_regs[70][0]);

    if (cfg->three_stage) {
        /* Forward Phase 0: 0 .. r0 - 1 (band lo = s0) */
        qua_dqi_emit_hgcd_phase(&ctx, tab, A_regs, B_regs, 0u, cfg->r0, cfg->s0, false, k_min, hs, do_div_hi, jump_at, false, cfg->pebble_scaff, cfg->mbuc_q0);
        /* Forward Checkpoint 0: lift [s1 .. s0) */
        qua_dqi_emit_hgcd_checkpoint(&ctx, tab, A_regs, B_regs, syn_q, jump_at, cfg->r0, s1, cfg->s0 - s1, s1);
        /* Forward Phase 1: r0 .. r1 - 1 (band lo = s1) */
        qua_dqi_emit_hgcd_phase(&ctx, tab, A_regs, B_regs, cfg->r0, r1, s1, false, k_min, hs, do_div_hi, jump_at, false, cfg->pebble_scaff, cfg->mbuc_q0);
    } else {
        /* Forward Phase 1: 0 .. r1 - 1 */
        qua_dqi_emit_hgcd_phase(&ctx, tab, A_regs, B_regs, 0u, r1, s1, false, k_min, hs, do_div_hi, jump_at, false, cfg->pebble_scaff, cfg->mbuc_q0);
    }

    /* Self-Inverse Midpoint Vandermonde Convolution (Checkpoint 1) */
    qua_dqi_emit_hgcd_checkpoint(&ctx, tab, A_regs, B_regs, syn_q, jump_at, r1, 0u, s1, k_min);

    /* Forward Phase 2: r1 .. 34 */
    qua_dqi_emit_hgcd_phase(&ctx, tab, A_regs, B_regs, r1, 35u, 0u, true, k_min, hs, do_div_hi, jump_at, false, cfg->pebble_scaff, cfg->mbuc_q0);

    /* Toeplitz + Chien-Forney + Reverse Toeplitz */
    {
        uint32_t *sigma0 = B_regs[70];
        uint32_t *sigma_off[34];
        for (uint32_t j = 0u; j < 34u; ++j) sigma_off[j] = B_regs[69u - j];
        uint32_t *syn_ptrs[35];
        for (uint32_t i = 0u; i < 35u; ++i) syn_ptrs[i] = syn_q[i];

        uint32_t lad_s0[11], inv_s0[12], sc_s0[21], t_work[12];
        qua_dqi_init_is_zero_12(&ctx, sigma0, lad_s0);
        uint32_t s0_nz = lad_s0[10];
        qua_x(&kb, s0_nz);
        qua_dqi_alloc_q(&ctx, inv_s0, 12u);
        qua_dqi_init_tower_inv_ex(&ctx, tab, sigma0, inv_s0, s0_nz, sc_s0, cfg->pebble_scaff);
        qua_dqi_alloc_q(&ctx, t_work, 12u);

        qua_dqi_toeplitz_rec(&ctx, tab, syn_ptrs, 35u, cfg->toeplitz_base, sigma0, inv_s0, s0_nz, sigma_off, t_work, false);

        qua_dqi_free_q(&ctx, t_work, 12u);
        qua_dqi_del_tower_inv_ex(&ctx, tab, sigma0, inv_s0, s0_nz, sc_s0, cfg->pebble_scaff);
        qua_dqi_free_q(&ctx, inv_s0, 12u);
        qua_x(&kb, s0_nz);
        qua_dqi_del_is_zero_12(&ctx, sigma0, lad_s0);

        uint32_t *B_chien[71];
        for (uint32_t i = 0u; i < 35u; ++i) B_chien[i] = syn_q[i];
        for (uint32_t i = 35u; i < 71u; ++i) B_chien[i] = B_regs[i];
        qua_dqi_chien_forney(&ctx, tab, B_chien, mbuc_b);

        qua_dqi_init_is_zero_12(&ctx, sigma0, lad_s0);
        s0_nz = lad_s0[10];
        qua_x(&kb, s0_nz);
        qua_dqi_alloc_q(&ctx, inv_s0, 12u);
        qua_dqi_init_tower_inv_ex(&ctx, tab, sigma0, inv_s0, s0_nz, sc_s0, cfg->pebble_scaff);
        qua_dqi_alloc_q(&ctx, t_work, 12u);

        qua_dqi_toeplitz_rec(&ctx, tab, syn_ptrs, 35u, cfg->toeplitz_base, sigma0, inv_s0, s0_nz, sigma_off, t_work, true);

        qua_dqi_free_q(&ctx, t_work, 12u);
        qua_dqi_del_tower_inv_ex(&ctx, tab, sigma0, inv_s0, s0_nz, sc_s0, cfg->pebble_scaff);
        qua_dqi_free_q(&ctx, inv_s0, 12u);
        qua_x(&kb, s0_nz);
        qua_dqi_del_is_zero_12(&ctx, sigma0, lad_s0);
    }

    /* Reverse Phase 2: 34 .. r1 */
    qua_dqi_emit_hgcd_phase(&ctx, tab, A_regs, B_regs, r1, 35u, 0u, true, k_min, hs, do_div_hi, jump_at, true, cfg->pebble_scaff, cfg->mbuc_q0);

    /* Reverse Midpoint Vandermonde Convolution (Checkpoint 1) */
    qua_dqi_emit_hgcd_checkpoint(&ctx, tab, A_regs, B_regs, syn_q, jump_at, r1, 0u, s1, k_min);

    if (cfg->three_stage) {
        /* Reverse Phase 1: r1 - 1 .. r0 (band lo = s1) */
        qua_dqi_emit_hgcd_phase(&ctx, tab, A_regs, B_regs, cfg->r0, r1, s1, false, k_min, hs, do_div_hi, jump_at, true, cfg->pebble_scaff, cfg->mbuc_q0);
        /* Reverse Checkpoint 0 */
        qua_dqi_emit_hgcd_checkpoint(&ctx, tab, A_regs, B_regs, syn_q, jump_at, cfg->r0, s1, cfg->s0 - s1, s1);
        /* Reverse Phase 0: r0 - 1 .. 0 (band lo = s0) */
        qua_dqi_emit_hgcd_phase(&ctx, tab, A_regs, B_regs, 0u, cfg->r0, cfg->s0, false, k_min, hs, do_div_hi, jump_at, true, cfg->pebble_scaff, cfg->mbuc_q0);
    } else {
        /* Reverse Phase 1: r1 - 1 .. 0 */
        qua_dqi_emit_hgcd_phase(&ctx, tab, A_regs, B_regs, 0u, r1, s1, false, k_min, hs, do_div_hi, jump_at, true, cfg->pebble_scaff, cfg->mbuc_q0);
    }

    qua_x(&kb, A_regs[70][0]);
    qua_x(&kb, B_regs[70][0]);
    for (uint32_t i = init_s; i < 35u; ++i) {
        for (uint32_t b = 0u; b < 12u; ++b) qua_cx(&kb, syn_q[i][b], B_regs[i][b]);
    }

    qua_bitop_spec_t fin_spec = {
        "dqi", QUA_COST_Q2T, false, 0u, 22u, 1e8, 64u, 2048u, 8u, NULL, NULL, NULL
    };
    int rc = qua_bitop_finalize_kmb(&kb, &fin_spec, out_path);
    qua_builder_free(&kb);
    return rc;
}

static inline void qua_dqi_populate_dag(const qua_dqi_cfg_t *c, qua_bitop_dag_t *dag) {
    if (!dag) return;
    dag->num_nodes      = 0u;
    dag->base_io_qubits = QUA_DQI_N * 12u;
    dag->cost_law       = QUA_COST_Q2T;
    uint32_t b_lo_q   = 12u * (35u - c->k_min) + 12u;
    uint32_t a_pool_q = 12u * (71u - c->k_min);
    uint32_t flags_q  = (35u - c->halt_start) + 34u;
    uint32_t sc_q     = c->pebble_scaff ? (c->mbuc_q0 ? 65u : 53u) : (c->mbuc_q0 ? 83u : 71u);
    uint32_t init_s   = c->three_stage ? c->s0 : c->s1;
    uint32_t n_b = qua_bitop_dag_push_rep(dag, QUA_BOP_XOR,   QUA_REP_LINF2,       b_lo_q,   12u * (35u - init_s), 0u, 0u, false, false, 0u);
    uint32_t n_a = qua_bitop_dag_push_rep(dag, QUA_BOP_IO,    QUA_REP_KET,         a_pool_q, 1u,                   0u, 0u, false, false, 0u);
    uint32_t n_f = qua_bitop_dag_push_rep(dag, QUA_BOP_AND,   QUA_REP_KET,         flags_q,  35u * 20u,            0u, 0u, false, false, (UINT64_C(1) << n_b));
    uint64_t pers_mask = (UINT64_C(1) << n_b) | (UINT64_C(1) << n_a) | (UINT64_C(1) << n_f);
    uint32_t n_p1 = qua_bitop_dag_push_rep(dag, QUA_BOP_MUL,  QUA_REP_PARITY_DUAL, sc_q,     (uint32_t)(c->est_toffolis / 5ull), 36u, 0u, true, false, pers_mask);
    uint32_t mv_ccx = ((c->r1 + c->s1) + (c->three_stage ? (c->r0 + (c->s0 - c->s1)) : 0u)) * 84u;
    uint32_t n_mv = qua_bitop_dag_push_rep(dag, QUA_BOP_MUL,  QUA_REP_PARITY_DUAL, 24u,      mv_ccx,               0u, 0u, true,  false, pers_mask | (UINT64_C(1) << n_p1));
    uint32_t n_p2 = qua_bitop_dag_push_rep(dag, QUA_BOP_MUL,  QUA_REP_PARITY_DUAL, sc_q,     (uint32_t)(c->est_toffolis / 5ull), 36u, 0u, true, false, pers_mask | (UINT64_C(1) << n_mv));
    uint32_t tp_q = c->pebble_scaff ? 39u : 57u;
    uint32_t n_tp = qua_bitop_dag_push_rep(dag, QUA_BOP_MUL,  QUA_REP_PARITY_DUAL, tp_q,     (uint32_t)qua_dqi_est_toeplitz_ccx(35u, c->toeplitz_base), c->pebble_scaff ? 36u : 0u, 0u, true, false, pers_mask | (UINT64_C(1) << n_p2));
    (void)qua_bitop_dag_push_rep(dag,           QUA_BOP_PHASE, QUA_REP_PHASE,       45u,     4095u * 80u,          0u, 0u, false, true,  pers_mask | (UINT64_C(1) << n_tp));
    qua_ele_pebble_plan_t pp;
    uint64_t eager = c->pebble_scaff ? ((UINT64_C(1) << n_p1) | (UINT64_C(1) << n_mv) | (UINT64_C(1) << n_p2) | (UINT64_C(1) << n_tp)) : 0u;
    qua_bitop_eval_dag_pebble(dag, eager, &pp, NULL, 0u);
}

static inline bool qua_dqi_eval_spec(
    const uint64_t *genome,
    qua_bitop_dag_t *dag,
    double *out_raw_cost,
    uint64_t *out_frust,
    void *user_ctx
) {
    const qua_dqi_tables_t *tab = (const qua_dqi_tables_t *)user_ctx;
    qua_dqi_cfg_t c;
    qua_dqi_eval_cfg(tab, genome[0], &c);
    if (!c.valid) {
        if (dag) {
            double lam = c.p_trunc * 102400.0;
            if (c.s1 < 8u || c.r1 + c.s1 > 71u) {
                dag->err_violation = 10.0 + (double)((c.s1 < 8u) ? (8u - c.s1) : (c.r1 + c.s1 - 71u));
                /* Steer r1 [0..3] */
                *out_frust = UINT64_C(0x000F);
            } else {
                dag->err_violation = (lam > 0.25) ? (lam - 0.25) : 1.0;
                /* Truncation rate too high: steer halt_start [4..8] */
                *out_frust = UINT64_C(0x01F0);
            }
        } else {
            *out_frust = UINT64_C(0x3F87FF);
        }
        return false;
    }
    qua_dqi_populate_dag(&c, dag);
    *out_raw_cost = c.est_score;
    /* Co-search active 18 genes (excluding pinned bits [11..14]) */
    *out_frust = UINT64_C(0x3F87FF);
    return true;
}

static inline int qua_dqi_run_pure_search(const char *out_kmb_path) {
    qua_dqi_tables_t *tab = (qua_dqi_tables_t *)malloc(sizeof(qua_dqi_tables_t));
    if (!tab) return 1;
    qua_dqi_init_tables(tab);

    qua_bitop_spec_t spec = {
        "dqi", QUA_COST_Q2T, false, 0u, 22u, 1e8, 64u, 2048u, 8u, qua_dqi_eval_spec, tab, NULL
    };
    uint64_t best_g[4] = {0u, 0u, 0u, 0u};
    double best_raw = 1e30;
    qua_bitop_solve(&spec, UINT64_C(0xD012026), best_g, &best_raw);

    qua_dqi_cfg_t best;
    qua_dqi_eval_cfg(tab, best_g[0], &best);

    printf("[DQI Pure-Search Winner] R1=%u, S1=%u, K_MIN=%u, HALT_START=%u, toep_base=%u, peb_sc=%d, mbuc_q0=%d, 3stage=%d, r0=%u, s0=%u | est_Q=%u, est_T=%llu, p_trunc=%.3e, E[R]=%.4f -> est_score=%.6e\n",
           best.r1, best.s1, best.k_min, best.halt_start, best.toeplitz_base,
           (int)best.pebble_scaff, (int)best.mbuc_q0,
           (int)best.three_stage, best.r0, best.s0,
           best.peak_qubits, (unsigned long long)best.est_toffolis,
           best.p_trunc, best.rep_factor, best.est_score);

    int rc = qua_dqi_emit_kmb(tab, &best, out_kmb_path);
    free(tab);
    return rc;
}

#endif /* QUA_DQI_H */
