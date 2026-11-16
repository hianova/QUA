/* SPDX-License-Identifier: MIT
 * dft_bootstrap.c — Self-Hosting SEARU x FLOW v7.0.1 Symbiotic DF-THC Compiler (v2)
 *
 * Strictly scoped inside qua/kaggle/dfthcsos_opt/.
 * Includes frozen FLOW headers (flow.h, searu.h) read-only (0 bytes modified).
 *
 * 4-Pass Self-Hosting Architecture:
 *   Pass 1   : Progressive Homotopy 4-Tensor Spectral Agglomeration + Apple-AMX
 *              Multi-Threaded Low-Rank `Y_target` Two-Loop L-BFGS coupled with
 *              `flow_step_atlas` 2-Jet Citardauq Strong-Wolfe Line Search +
 *              Exact `T_packed` Block-QR-Eigen Variable Projection (`w_r`).
 *   Pass 2   : First-Principles Orbital Marginal Variance + `SEARU_SOLVE_INLINE`
 *              Multi-Chain MSB-0 Givens Spanning-Tree Synthesis (`tree_edges[75][2]`),
 *              forcing 67..70 of 75 edges to `MSB == 0` across ALL QROM rows.
 *   Pass 2.5 : Multi-Threaded Reverse-Mode Givens-Tree Angle-Space QAT Homotopy
 *              Two-Loop L-BFGS + `FLOW` 2-Jet Line Search (`(alpha, w)` joint).
 *   Pass 3   : Exact Live-Residual 4-Tensor Prefix-Subtree 1-Edge/2-Edge +
 *              `SEARU` CLZ Zipf Sparse-k Coupled Lattice Crucible interleaved
 *              with exact `T_packed` Block-QR-Eigen `w` projection.
 */

#ifndef ACCELERATE_NEW_LAPACK
#define ACCELERATE_NEW_LAPACK
#endif
#include <Accelerate/Accelerate.h>
#include <dispatch/dispatch.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef __APPLE__
typedef __LAPACK_int dft_lapack_int;
#else
typedef int dft_lapack_int;
#endif

#include "flow.h"
#include "searu.h"

#define DFT_N  76
#define DFT_E  75
#define DFT_K  2926 /* 76 * 77 / 2 */

static const double kPi    = 3.14159265358979323846264338327950288;
static const double kSqrt2 = 1.41421356237309504880168872420969808;

static double wall_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

typedef struct {
    int      N;
    int      K;
    int      electrons;
    int      R0, B0, C0;
    int      K_pos;
    double   beta_1;
    double   ecore;
    double   energy_upper;
    double   Y_norm2;
    double   T_norm2;
    double  *T_packed;        /* K x K */
    double  *Y_target;        /* K x K_pos row-major */
    double  *H1_base;         /* N x N */
    double  *H1_base_packed;  /* K */
    double  *shift_1b_packed; /* K */
    double  *U0;              /* R0 x B0 x N */
    double  *w0;              /* R0 x B0 x C0 */
    double  *W0;              /* R0 x C0 */
} dft_ref_t;

static void load_dft_ref(const char *path, dft_ref_t *ref) {
    FILE *fp = fopen(path, "rb");
    if (!fp) {
        fprintf(stderr, "FATAL: cannot open %s\n", path);
        exit(1);
    }
    uint32_t hdr[8];
    if (fread(hdr, sizeof(uint32_t), 8, fp) != 8 || hdr[0] != 0x44465432U) {
        fprintf(stderr, "FATAL: bad magic in %s (expected v2 0x44465432)\n", path);
        exit(1);
    }
    ref->N         = (int)hdr[1];
    ref->K         = (int)hdr[2];
    ref->electrons = (int)hdr[3];
    ref->R0        = (int)hdr[4];
    ref->B0        = (int)hdr[5];
    ref->C0        = (int)hdr[6];
    ref->K_pos     = (int)hdr[7];

    double d4[4];
    if (fread(d4, sizeof(double), 4, fp) != 4) {
        exit(1);
    }
    ref->beta_1       = d4[0];
    ref->ecore        = d4[1];
    ref->energy_upper = d4[2];
    ref->Y_norm2      = d4[3];

    size_t kk = (size_t)DFT_K * (size_t)DFT_K;
    size_t ky = (size_t)DFT_K * (size_t)ref->K_pos;
    ref->T_packed        = (double *)malloc(kk * sizeof(double));
    ref->Y_target        = (double *)malloc(ky * sizeof(double));
    ref->H1_base         = (double *)malloc((size_t)(DFT_N * DFT_N) * sizeof(double));
    ref->H1_base_packed  = (double *)malloc((size_t)DFT_K * sizeof(double));
    ref->shift_1b_packed = (double *)malloc((size_t)DFT_K * sizeof(double));
    ref->U0 = (double *)malloc((size_t)(ref->R0 * ref->B0 * DFT_N) * sizeof(double));
    ref->w0 = (double *)malloc((size_t)(ref->R0 * ref->B0 * ref->C0) * sizeof(double));
    ref->W0 = (double *)malloc((size_t)(ref->R0 * ref->C0) * sizeof(double));
    if (!ref->T_packed || !ref->Y_target || !ref->H1_base || !ref->H1_base_packed ||
        !ref->shift_1b_packed || !ref->U0 || !ref->w0 || !ref->W0) {
        fprintf(stderr, "FATAL: OOM in load_dft_ref\n");
        abort();
    }

    if (fread(ref->T_packed, sizeof(double), kk, fp) != kk ||
        fread(ref->Y_target, sizeof(double), ky, fp) != ky ||
        fread(ref->H1_base, sizeof(double), (size_t)(DFT_N * DFT_N), fp) != (size_t)(DFT_N * DFT_N) ||
        fread(ref->H1_base_packed, sizeof(double), DFT_K, fp) != DFT_K ||
        fread(ref->shift_1b_packed, sizeof(double), DFT_K, fp) != DFT_K ||
        fread(ref->U0, sizeof(double), (size_t)(ref->R0 * ref->B0 * DFT_N), fp) !=
            (size_t)(ref->R0 * ref->B0 * DFT_N) ||
        fread(ref->w0, sizeof(double), (size_t)(ref->R0 * ref->B0 * ref->C0), fp) !=
            (size_t)(ref->R0 * ref->B0 * ref->C0) ||
        fread(ref->W0, sizeof(double), (size_t)(ref->R0 * ref->C0), fp) !=
            (size_t)(ref->R0 * ref->C0)) {
        fprintf(stderr, "FATAL: truncated %s\n", path);
        exit(1);
    }
    fclose(fp);

    double tn2 = 0.0;
    for (size_t i = 0; i < kk; ++i) {
        tn2 += ref->T_packed[i] * ref->T_packed[i];
    }
    ref->T_norm2 = tn2;
}

static inline void pack_u_vec(const double *__restrict__ u, double *__restrict__ p) {
    int k = 0;
    for (int r = 0; r < DFT_N; ++r) {
        const double ur    = u[r];
        const double ur_s2 = ur * kSqrt2;
        for (int c = 0; c < r; ++c) {
            p[k++] = ur_s2 * u[c];
        }
        p[k++] = ur * ur;
    }
}

static inline void normalize_u_vec(double *u) {
    double s2 = 0.0;
    for (int i = 0; i < DFT_N; ++i) s2 += u[i] * u[i];
    double inv = (s2 > 1e-30) ? (1.0 / sqrt(s2)) : 1.0;
    for (int i = 0; i < DFT_N; ++i) u[i] *= inv;
}

static inline void grad_pack_u_vec(const double *__restrict__ gp, const double *__restrict__ u,
                                   double *__restrict__ gu) {
    for (int i = 0; i < DFT_N; ++i) gu[i] = 0.0;
    int k = 0;
    for (int r = 0; r < DFT_N; ++r) {
        const double ur    = u[r];
        const double ur_s2 = ur * kSqrt2;
        double       acc   = 0.0;
        for (int c = 0; c < r; ++c) {
            double g = gp[k++];
            acc += g * u[c];
            gu[c] += g * ur_s2;
        }
        double gd = gp[k++];
        gu[r] += acc * kSqrt2 + 2.0 * gd * ur;
    }
}

static void compute_block_X(const double *U_r, const double *w_r, int B, int C, double *P_r,
                            double *X_r) {
    for (int b = 0; b < B; ++b) {
        pack_u_vec(U_r + (size_t)b * DFT_N, P_r + (size_t)b * DFT_K);
    }
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, C, DFT_K, B, 1.0, w_r, C, P_r, DFT_K, 0.0,
                X_r, DFT_K);
}

static double frob_norm_sym_upper(const double *M, int n) {
    double diag_s = 0.0;
    double off_s  = 0.0;
    for (int j = 0; j < n; ++j) {
        diag_s += M[(size_t)j * n + j] * M[(size_t)j * n + j];
        for (int i = j + 1; i < n; ++i) {
            double v = M[(size_t)j * n + i];
            off_s += v * v;
        }
    }
    return sqrt(diag_s + 2.0 * off_s);
}

static double exact_err2_from_uw(const dft_ref_t *ref, const double *U, const double *w, int R,
                                 int B, int C) {
    int     M     = R * B;
    int     RC    = R * C;
    double *P_all = (double *)malloc((size_t)M * DFT_K * sizeof(double));
    double *X_all = (double *)malloc((size_t)RC * DFT_K * sizeof(double));
    for (int m = 0; m < M; ++m) {
        double u_tmp[DFT_N];
        memcpy(u_tmp, U + (size_t)m * DFT_N, DFT_N * sizeof(double));
        normalize_u_vec(u_tmp);
        pack_u_vec(u_tmp, P_all + (size_t)m * DFT_K);
    }
    for (int r = 0; r < R; ++r) {
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, C, DFT_K, B, 1.0,
                    w + (size_t)r * B * C, C, P_all + (size_t)r * B * DFT_K, DFT_K, 0.0,
                    X_all + (size_t)r * C * DFT_K, DFT_K);
    }
    /* Exact ||X^T X - T_packed||_F^2 = ||G_XX||_F^2 - 2 Tr(X T_packed X^T) + ||T_packed||_F^2 */
    double *G_XX = (double *)malloc((size_t)RC * RC * sizeof(double));
    double *TX   = (double *)malloc((size_t)RC * DFT_K * sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, RC, RC, DFT_K, 1.0, X_all, DFT_K, X_all,
                DFT_K, 0.0, G_XX, RC);
    cblas_dsymm(CblasRowMajor, CblasRight, CblasUpper, RC, DFT_K, 1.0, ref->T_packed, DFT_K, X_all,
                DFT_K, 0.0, TX, DFT_K);
    double g2 = 0.0, xtx = 0.0;
    for (size_t i = 0; i < (size_t)RC * RC; ++i) g2 += G_XX[i] * G_XX[i];
    for (size_t i = 0; i < (size_t)RC * DFT_K; ++i) xtx += X_all[i] * TX[i];
    double err2 = sqrt(fmax(0.0, g2 - 2.0 * xtx + ref->T_norm2));
    free(TX);
    free(G_XX);
    free(X_all);
    free(P_all);
    return err2;
}

/* ============================================================================
 * Generic Self-Hosting Two-Loop L-BFGS + FLOW 2-Jet Citardauq Line Search
 * ============================================================================ */

typedef double (*lbfgs_fg_fn)(const double *x, double *g, void *ctx);

/* Computes FLOW 2-Jet Citardauq stationary contraction u* in [0.10, 0.65] for 1D ray phi(u * alpha) */
static double flow_2jet_line_contract(double f0, double dphi0_scaled, double f1) {
    double curv = 2.0 * (f1 - f0 - dphi0_scaled);
    double u_st = -dphi0_scaled / fmax(1e-24, curv);
    if (curv > 1e-14 && dphi0_scaled < 0.0) {
        double scale   = 4096.0 / fmax(1e-12, fabs(dphi0_scaled) + fabs(curv));
        int64_t b0_q16 = (int64_t)llround(dphi0_scaled * scale * (double)FLOW_Q16_ONE);
        int64_t b1_q16 = (int64_t)llround(curv * scale * (double)FLOW_Q16_ONE);
        if (b1_q16 > 0) {
            FLOW_ALIGN64 flow_cell arena[1];
            const uint64_t chart[8] = {
                (uint64_t)flow_clamp_sym(b0_q16, FLOW_MAX_Q16),
                (uint64_t)flow_clamp_sym(b1_q16, FLOW_MAX_Q16),
                UINT64_C(0),
                ((uint64_t)(uint32_t)(2 * FLOW_Q16_ONE) << 32U) | (uint64_t)(uint32_t)0,
                ((uint64_t)(uint32_t)(2 * FLOW_Q16_ONE) << 32U) | (uint64_t)(uint32_t)(FLOW_Q16_ONE / 2),
                0x0001ULL,
                UINT64_C(1) << 32U,
                0ULL,
            };
            flow_arena_publish_bound(arena, 0U, 0U, 1U, chart);
            flow_cell st = {0ULL, 0ULL, 0ULL, (uint64_t)(FLOW_Q16_ONE / 2), 0ULL, 0x0001ULL, 0ULL, 0ULL};
            for (int s = 0; s < 6; ++s) {
                flow_step_atlas(st, st, arena, 0ULL, 0LL, 0x0001U);
                if (st[6] & FLOW_STATUS_CONVERGED) break;
            }
            double uf = (double)(int64_t)st[3] / (double)FLOW_Q16_ONE;
            if (uf >= 0.08 && uf <= 0.75) {
                u_st = uf;
            }
        }
    }
    return fmin(0.65, fmax(0.12, u_st));
}

static double run_lbfgs_flow2jet(int n_dim, double *x, int max_iter, int m_hist, lbfgs_fg_fn fg,
                                 void *ctx, const char *tag, int log_every) {
    double *g      = (double *)malloc((size_t)n_dim * sizeof(double));
    double *g_prev = (double *)malloc((size_t)n_dim * sizeof(double));
    double *x_prev = (double *)malloc((size_t)n_dim * sizeof(double));
    double *x_try  = (double *)malloc((size_t)n_dim * sizeof(double));
    double *g_try  = (double *)malloc((size_t)n_dim * sizeof(double));
    double *d      = (double *)malloc((size_t)n_dim * sizeof(double));
    double *q      = (double *)malloc((size_t)n_dim * sizeof(double));
    double *S_hist = (double *)malloc((size_t)m_hist * (size_t)n_dim * sizeof(double));
    double *Y_hist = (double *)malloc((size_t)m_hist * (size_t)n_dim * sizeof(double));
    double  rho[64], alpha_h[64];

    double t0 = wall_sec();
    double f  = fg(x, g, ctx);

    int h_cnt  = 0;
    int h_head = 0; /* oldest index when full, next write slot is (h_head + h_cnt) % m_hist */

    for (int it = 0; it < max_iter; ++it) {
        double gnorm = cblas_dnrm2(n_dim, g, 1);
        if (gnorm < 1e-10) break;

        /* Two-loop L-BFGS recursion: d = -H_k * g */
        memcpy(q, g, (size_t)n_dim * sizeof(double));
        for (int c = h_cnt - 1; c >= 0; --c) {
            int           idx = (h_head + c) % m_hist;
            const double *s_i = S_hist + (size_t)idx * n_dim;
            const double *y_i = Y_hist + (size_t)idx * n_dim;
            double        a_i = rho[idx] * cblas_ddot(n_dim, s_i, 1, q, 1);
            alpha_h[c]        = a_i;
            cblas_daxpy(n_dim, -a_i, y_i, 1, q, 1);
        }

        double gamma_0 = 1.0;
        if (h_cnt > 0) {
            int           last = (h_head + h_cnt - 1) % m_hist;
            const double *s_l  = S_hist + (size_t)last * n_dim;
            const double *y_l  = Y_hist + (size_t)last * n_dim;
            double        sy   = cblas_ddot(n_dim, s_l, 1, y_l, 1);
            double        yy   = cblas_ddot(n_dim, y_l, 1, y_l, 1);
            if (yy > 1e-30 && sy > 0.0) gamma_0 = sy / yy;
        }
        cblas_dscal(n_dim, gamma_0, q, 1);

        for (int c = 0; c < h_cnt; ++c) {
            int           idx = (h_head + c) % m_hist;
            const double *s_i = S_hist + (size_t)idx * n_dim;
            const double *y_i = Y_hist + (size_t)idx * n_dim;
            double        b_i = rho[idx] * cblas_ddot(n_dim, y_i, 1, q, 1);
            cblas_daxpy(n_dim, alpha_h[c] - b_i, s_i, 1, q, 1);
        }

        for (int i = 0; i < n_dim; ++i) d[i] = -q[i];

        double dphi0 = cblas_ddot(n_dim, g, 1, d, 1);
        if (dphi0 >= -1e-20) {
            /* Reset L-BFGS history if not a descent direction */
            h_cnt  = 0;
            h_head = 0;
            for (int i = 0; i < n_dim; ++i) d[i] = -g[i];
            dphi0 = -gnorm * gnorm;
        }

        /* FLOW 2-Jet Strong/Weak-Wolfe Line Search */
        double step     = (h_cnt == 0) ? (1.0 / fmax(1.0, gnorm)) : 1.0;
        bool   accepted = false;
        double f_best   = f;

        for (int ls = 0; ls < 10; ++ls) {
            memcpy(x_try, x, (size_t)n_dim * sizeof(double));
            cblas_daxpy(n_dim, step, d, 1, x_try, 1);
            double f_try    = fg(x_try, g_try, ctx);
            double dphi_try = cblas_ddot(n_dim, g_try, 1, d, 1);

            if (f_try <= f + 1e-4 * step * dphi0) {
                if (dphi_try >= 0.9 * dphi0 || ls >= 3) {
                    f_best   = f_try;
                    accepted = true;
                    break;
                }
                /* Slope still very negative -> expand step once or twice! */
                double step_ext = 2.5 * step;
                double *x_ext   = q; /* reuse scratch buffer q */
                memcpy(x_ext, x, (size_t)n_dim * sizeof(double));
                cblas_daxpy(n_dim, step_ext, d, 1, x_ext, 1);
                double *g_ext = (double *)malloc((size_t)n_dim * sizeof(double));
                double  f_ext = fg(x_ext, g_ext, ctx);
                if (f_ext < f_try) {
                    memcpy(x_try, x_ext, (size_t)n_dim * sizeof(double));
                    memcpy(g_try, g_ext, (size_t)n_dim * sizeof(double));
                    f_best = f_ext;
                    step   = step_ext;
                } else {
                    f_best = f_try;
                }
                free(g_ext);
                accepted = true;
                break;
            }
            /* Overshot: contract via FLOW 2-Jet Citardauq stationary step! */
            double u_fac = flow_2jet_line_contract(f, step * dphi0, f_try);
            step *= u_fac;
        }

        if (!accepted) {
            if (h_cnt > 0) {
                /* Reset stale history and retry steepest descent on next iter */
                h_cnt  = 0;
                h_head = 0;
                continue;
            }
            break;
        }

        memcpy(x_prev, x, (size_t)n_dim * sizeof(double));
        memcpy(g_prev, g, (size_t)n_dim * sizeof(double));
        memcpy(x, x_try, (size_t)n_dim * sizeof(double));
        memcpy(g, g_try, (size_t)n_dim * sizeof(double));
        double df_rel = (f - f_best) / fmax(1e-12, fabs(f));
        f             = f_best;

        /* Update L-BFGS history pair (s = x - x_prev, y = g - g_prev) */
        int     slot = (h_cnt < m_hist) ? ((h_head + h_cnt) % m_hist) : h_head;
        double *s_sl = S_hist + (size_t)slot * n_dim;
        double *y_sl = Y_hist + (size_t)slot * n_dim;
        for (int i = 0; i < n_dim; ++i) {
            s_sl[i] = x[i] - x_prev[i];
            y_sl[i] = g[i] - g_prev[i];
        }
        double sy = cblas_ddot(n_dim, s_sl, 1, y_sl, 1);
        double ss = cblas_ddot(n_dim, s_sl, 1, s_sl, 1);
        double yy = cblas_ddot(n_dim, y_sl, 1, y_sl, 1);
        if (sy > 1e-13 * sqrt(ss * yy + 1e-60)) {
            rho[slot] = 1.0 / sy;
            if (h_cnt < m_hist) {
                ++h_cnt;
            } else {
                h_head = (h_head + 1) % m_hist;
            }
        }

        if (log_every > 0 && ((it + 1) % log_every == 0 || it == 0)) {
            printf("    [%s it=%4d/%d] f=%.7e (approx_err2=%.3f mHa) |g|=%.2e [%.2fs]\n", tag,
                   it + 1, max_iter, f, 1e3 * sqrt(fmax(0.0, 4.0 * f)), gnorm, wall_sec() - t0);
            fflush(stdout);
        }
        if (it > 25 && df_rel >= 0.0 && df_rel < 1e-12) break;
    }

    free(Y_hist);
    free(S_hist);
    free(q);
    free(d);
    free(g_try);
    free(x_try);
    free(x_prev);
    free(g_prev);
    free(g);
    return f;
}

/* ============================================================================
 * Pass 1: Continuous Low-Rank Y_target L-BFGS + Progressive Agglomeration
 * ============================================================================ */

typedef struct {
    const dft_ref_t *ref;
    int              R, B, C;
    double           w_scale;
    double           l1_reg;
    double          *Un_buf;     /* M x N */
    double          *norms_buf;  /* M */
    double          *P_buf;      /* M x K */
    double          *X_buf;      /* RC x K */
    double          *G_XX_buf;   /* RC x RC */
    double          *XY_buf;     /* RC x K_pos */
    double          *grad_X_buf; /* RC x K */
    double          *grad_P_buf; /* M x K */
} cont_uw_ctx_t;

static double fg_continuous_uw(const double *x, double *g, void *vctx) {
    cont_uw_ctx_t   *c    = (cont_uw_ctx_t *)vctx;
    const dft_ref_t *ref  = c->ref;
    const int        R    = c->R;
    const int        B    = c->B;
    const int        C    = c->C;
    const int        M    = R * B;
    const int        RC   = R * C;
    const int        u_sz = M * DFT_N;
    const double     ws   = c->w_scale;

    const double *U_in    = x;
    const double *wt_in   = x + u_sz;
    double       *Un      = c->Un_buf;
    double       *norms   = c->norms_buf;
    double       *P       = c->P_buf;
    double       *X_flat  = c->X_buf;
    double       *G_XX    = c->G_XX_buf;
    double       *XY      = c->XY_buf;
    double       *grad_X  = c->grad_X_buf;
    double       *grad_P  = c->grad_P_buf;
    double       *grad_U  = g;
    double       *grad_wt = g + u_sz;

    /* 1. Multi-threaded row normalization + outer-product packing */
    const int chunk    = 64;
    const int n_chunks = (M + chunk - 1) / chunk;
    dispatch_apply((size_t)n_chunks, dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0),
                   ^(size_t cid) {
                       int lo = (int)cid * chunk;
                       int hi = (lo + chunk < M) ? (lo + chunk) : M;
                       for (int m = lo; m < hi; ++m) {
                           const double *ui = U_in + (size_t)m * DFT_N;
                           double       *un = Un + (size_t)m * DFT_N;
                           double        s2 = 0.0;
                           for (int j = 0; j < DFT_N; ++j) s2 += ui[j] * ui[j];
                           double nm  = sqrt(fmax(1e-30, s2));
                           double inv = 1.0 / nm;
                           norms[m]   = nm;
                           for (int j = 0; j < DFT_N; ++j) un[j] = ui[j] * inv;
                           pack_u_vec(un, P + (size_t)m * DFT_K);
                       }
                   });

    /* 2. Block contraction: X_flat[r] = (ws * wt_in[r])^T * P[r] */
    for (int r = 0; r < R; ++r) {
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, C, DFT_K, B, ws,
                    wt_in + (size_t)r * B * C, C, P + (size_t)r * B * DFT_K, DFT_K, 0.0,
                    X_flat + (size_t)r * C * DFT_K, DFT_K);
    }

    /* 3. Fast Global Gram G_XX (RC x RC) & Low-Rank Target XY (RC x K_pos) */
    cblas_dsyrk(CblasRowMajor, CblasUpper, CblasNoTrans, RC, DFT_K, 1.0, X_flat, DFT_K, 0.0, G_XX,
                RC);
    for (int i = 0; i < RC; ++i) {
        for (int j = 0; j < i; ++j) {
            G_XX[(size_t)i * RC + j] = G_XX[(size_t)j * RC + i];
        }
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, RC, ref->K_pos, DFT_K, 1.0, X_flat,
                DFT_K, ref->Y_target, ref->K_pos, 0.0, XY, ref->K_pos);

    double g2 = 0.0, xy2 = 0.0;
    for (size_t i = 0; i < (size_t)RC * RC; ++i) g2 += G_XX[i] * G_XX[i];
    for (size_t i = 0; i < (size_t)RC * ref->K_pos; ++i) xy2 += XY[i] * XY[i];

    double err2_sq = fmax(0.0, g2 - 2.0 * xy2 + ref->Y_norm2);
    double loss    = 0.25 * err2_sq;

    /* 4. Exact gradient grad_X = G_XX * X_flat - XY * Y_target^T */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, RC, DFT_K, RC, 1.0, G_XX, RC, X_flat,
                DFT_K, 0.0, grad_X, DFT_K);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, RC, DFT_K, ref->K_pos, -1.0, XY,
                ref->K_pos, ref->Y_target, ref->K_pos, 1.0, grad_X, DFT_K);

    /* 5. Adjoint wrt wt and P */
    for (int r = 0; r < R; ++r) {
        /* grad_wt[r] = ws * P[r] * grad_X[r]^T */
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, B, C, DFT_K, ws,
                    P + (size_t)r * B * DFT_K, DFT_K, grad_X + (size_t)r * C * DFT_K, DFT_K, 0.0,
                    grad_wt + (size_t)r * B * C, C);
        /* grad_P[r] = (ws * wt_in[r]) * grad_X[r] */
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, B, DFT_K, C, ws,
                    wt_in + (size_t)r * B * C, C, grad_X + (size_t)r * C * DFT_K, DFT_K, 0.0,
                    grad_P + (size_t)r * B * DFT_K, DFT_K);
    }

    if (c->l1_reg > 0.0) {
        const double eps2 = 1e-12;
        for (int r = 0; r < R; ++r) {
            for (int col = 0; col < C; ++col) {
                double L_rc = 0.0;
                for (int b = 0; b < B; ++b) {
                    double wv = ws * wt_in[(r * B + b) * C + col];
                    L_rc += sqrt(wv * wv + eps2);
                }
                loss += c->l1_reg * L_rc * L_rc;
                double fac = 2.0 * c->l1_reg * L_rc * ws;
                for (int b = 0; b < B; ++b) {
                    int    idx = (r * B + b) * C + col;
                    double wv  = ws * wt_in[idx];
                    double aw  = sqrt(wv * wv + eps2);
                    grad_wt[idx] += fac * (wv / aw);
                }
            }
        }
    }

    /* 6. Multi-threaded adjoint grad_P -> grad_U (projected onto sphere tangent / norms) */
    dispatch_apply((size_t)n_chunks, dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0),
                   ^(size_t cid) {
                       int    lo = (int)cid * chunk;
                       int    hi = (lo + chunk < M) ? (lo + chunk) : M;
                       double gu[DFT_N];
                       for (int m = lo; m < hi; ++m) {
                           const double *un = Un + (size_t)m * DFT_N;
                           double       *go = grad_U + (size_t)m * DFT_N;
                           grad_pack_u_vec(grad_P + (size_t)m * DFT_K, un, gu);
                           double dot = 0.0;
                           for (int j = 0; j < DFT_N; ++j) dot += gu[j] * un[j];
                           double inv_nm = 1.0 / norms[m];
                           for (int j = 0; j < DFT_N; ++j) go[j] = (gu[j] - dot * un[j]) * inv_nm;
                       }
                   });

    return loss;
}

static void run_continuous_lbfgs_stage(const dft_ref_t *ref, int R, int B, int C, int max_iter,
                                       double l1_reg, double *U, double *w) {
    int           M    = R * B;
    int           RC   = R * C;
    int           u_sz = M * DFT_N;
    int           w_sz = M * C;
    int           ndim = u_sz + w_sz;
    const double  ws   = 0.05;

    cont_uw_ctx_t ctx = {
        .ref        = ref,
        .R          = R,
        .B          = B,
        .C          = C,
        .w_scale    = ws,
        .l1_reg     = l1_reg,
        .Un_buf     = (double *)malloc((size_t)M * DFT_N * sizeof(double)),
        .norms_buf  = (double *)malloc((size_t)M * sizeof(double)),
        .P_buf      = (double *)malloc((size_t)M * DFT_K * sizeof(double)),
        .X_buf      = (double *)malloc((size_t)RC * DFT_K * sizeof(double)),
        .G_XX_buf   = (double *)malloc((size_t)RC * RC * sizeof(double)),
        .XY_buf     = (double *)malloc((size_t)RC * ref->K_pos * sizeof(double)),
        .grad_X_buf = (double *)malloc((size_t)RC * DFT_K * sizeof(double)),
        .grad_P_buf = (double *)malloc((size_t)M * DFT_K * sizeof(double)),
    };

    double *x = (double *)malloc((size_t)ndim * sizeof(double));
    for (int m = 0; m < M; ++m) {
        memcpy(x + (size_t)m * DFT_N, U + (size_t)m * DFT_N, DFT_N * sizeof(double));
        normalize_u_vec(x + (size_t)m * DFT_N);
    }
    for (int i = 0; i < w_sz; ++i) x[u_sz + i] = w[i] / ws;

    char tag[64];
    snprintf(tag, sizeof(tag), "Cont-LBFGS (%d,%d,%d)", R, B, C);
    run_lbfgs_flow2jet(ndim, x, max_iter, 20, fg_continuous_uw, &ctx, tag, 300);

    for (int m = 0; m < M; ++m) {
        memcpy(U + (size_t)m * DFT_N, x + (size_t)m * DFT_N, DFT_N * sizeof(double));
        normalize_u_vec(U + (size_t)m * DFT_N);
    }
    for (int i = 0; i < w_sz; ++i) w[i] = x[u_sz + i] * ws;

    free(x);
    free(ctx.grad_P_buf);
    free(ctx.grad_X_buf);
    free(ctx.XY_buf);
    free(ctx.G_XX_buf);
    free(ctx.X_buf);
    free(ctx.P_buf);
    free(ctx.norms_buf);
    free(ctx.Un_buf);
}

/* Exact closed-form Block-QR-Eigen `w_r` against full `T_packed` residual `Res` */
static void solve_block_w_qr_eigen(const double *P_r, const double *Res, int B, int C,
                                   double sv_floor, double *w_r_out, double *HQ_buf,
                                   double *Q_col_buf) {
    memcpy(Q_col_buf, P_r, (size_t)B * DFT_K * sizeof(double));
    dft_lapack_int m = DFT_K, n = B, lda = DFT_K, info = 0;
    double         tau[128], qwork = 0.0;
    dft_lapack_int lwork = -1;
    dgeqrf_(&m, &n, Q_col_buf, &lda, tau, &qwork, &lwork, &info);
    lwork        = (dft_lapack_int)qwork + 64;
    double *work = (double *)malloc((size_t)lwork * sizeof(double));
    dgeqrf_(&m, &n, Q_col_buf, &lda, tau, work, &lwork, &info);

    double R_col[128 * 128];
    memset(R_col, 0, sizeof(R_col));
    for (int j = 0; j < B; ++j) {
        for (int i = 0; i <= j; ++i) {
            R_col[j * B + i] = Q_col_buf[(size_t)j * DFT_K + i];
        }
    }
    dorgqr_(&m, &n, &n, Q_col_buf, &lda, tau, work, &lwork, &info);

    cblas_dsymm(CblasColMajor, CblasLeft, CblasLower, DFT_K, B, 1.0, Res, DFT_K, Q_col_buf, DFT_K,
                0.0, HQ_buf, DFT_K);

    double SB[128 * 128];
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, B, B, DFT_K, 1.0, Q_col_buf, DFT_K, HQ_buf,
                DFT_K, 0.0, SB, B);

    double         ev[128];
    dft_lapack_int nb   = B;
    dft_lapack_int lw_e = 256 * B;
    double        *wk_e = (double *)malloc((size_t)lw_e * sizeof(double));
    char           jobz = 'V', uplo = 'U';
    dsyev_(&jobz, &uplo, &nb, SB, &nb, ev, wk_e, &lw_e, &info);

    double Y_col[128 * 128];
    memset(Y_col, 0, sizeof(Y_col));
    int keep = (C < B) ? C : B;
    for (int c = 0; c < keep; ++c) {
        int    idx = B - 1 - c;
        double lam = ev[idx];
        if (lam > 1e-13) {
            double rt = sqrt(lam);
            for (int i = 0; i < B; ++i) Y_col[c * B + i] = SB[idx * B + i] * rt;
        }
    }

    double         U_svd[128 * 128], Vt_svd[128 * 128], s_svd[128];
    char           jobu = 'A', jobvt = 'A';
    dft_lapack_int lw_s = 256 * B;
    double        *wk_s = (double *)malloc((size_t)lw_s * sizeof(double));
    dgesvd_(&jobu, &jobvt, &nb, &nb, R_col, &nb, s_svd, U_svd, &nb, Vt_svd, &nb, wk_s, &lw_s,
            &info);

    double tmp1[128 * 128];
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, B, C, B, 1.0, U_svd, B, Y_col, B, 0.0,
                tmp1, B);
    double reg2 = sv_floor * sv_floor;
    for (int c = 0; c < C; ++c) {
        for (int i = 0; i < B; ++i) {
            double si       = s_svd[i];
            double inv_damp = (si > 1e-12) ? (si / (si * si + reg2)) : 0.0;
            tmp1[c * B + i] *= inv_damp;
        }
    }
    double W_col[128 * 128];
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, B, C, B, 1.0, Vt_svd, B, tmp1, B, 0.0,
                W_col, B);

    for (int b = 0; b < B; ++b) {
        for (int c = 0; c < C; ++c) w_r_out[b * C + c] = W_col[c * B + b];
    }

    free(wk_s);
    free(wk_e);
    free(work);
}

static void polish_w_exact_varpro(const dft_ref_t *ref, int R, int B, int C, int n_sweeps,
                                  double sv_floor, const double *U, double *w) {
    size_t  kk  = (size_t)DFT_K * (size_t)DFT_K;
    double *Res = (double *)malloc(kk * sizeof(double));
    double *P_a = (double *)malloc((size_t)R * B * DFT_K * sizeof(double));
    double *X_a = (double *)malloc((size_t)R * C * DFT_K * sizeof(double));
    double *HQ  = (double *)malloc((size_t)B * DFT_K * sizeof(double));
    double *Qb  = (double *)malloc((size_t)B * DFT_K * sizeof(double));

    memcpy(Res, ref->T_packed, kk * sizeof(double));
    for (int r = 0; r < R; ++r) {
        compute_block_X(U + (size_t)r * B * DFT_N, w + (size_t)r * B * C, B, C,
                        P_a + (size_t)r * B * DFT_K, X_a + (size_t)r * C * DFT_K);
    }
    cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, DFT_K, R * C, -1.0, X_a, DFT_K, 1.0, Res,
                DFT_K);

    for (int sw = 0; sw < n_sweeps; ++sw) {
        for (int r = 0; r < R; ++r) {
            double *Xr = X_a + (size_t)r * C * DFT_K;
            double *wr = w + (size_t)r * B * C;
            double *Pr = P_a + (size_t)r * B * DFT_K;
            cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, DFT_K, C, 1.0, Xr, DFT_K, 1.0, Res,
                        DFT_K);
            solve_block_w_qr_eigen(Pr, Res, B, C, sv_floor, wr, HQ, Qb);
            cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, C, DFT_K, B, 1.0, wr, C, Pr, DFT_K,
                        0.0, Xr, DFT_K);
            cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, DFT_K, C, -1.0, Xr, DFT_K, 1.0,
                        Res, DFT_K);
        }
    }
    free(Qb);
    free(HQ);
    free(X_a);
    free(P_a);
    free(Res);
}

static void prune_and_fit_block_to_X(const double *U_in, int B_in, int B_out, int C,
                                     const double *X_target, double *U_out, double *w_out) {
    double *P_in = (double *)malloc((size_t)B_in * DFT_K * sizeof(double));
    for (int b = 0; b < B_in; ++b) {
        double u_tmp[DFT_N];
        memcpy(u_tmp, U_in + (size_t)b * DFT_N, DFT_N * sizeof(double));
        normalize_u_vec(u_tmp);
        pack_u_vec(u_tmp, P_in + (size_t)b * DFT_K);
    }

    int active[128];
    int m = B_in;
    for (int i = 0; i < m; ++i) active[i] = i;

    double G_full[160 * 160];
    double PX_full[160 * 128];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, B_in, B_in, DFT_K, 1.0, P_in, DFT_K, P_in,
                DFT_K, 0.0, G_full, B_in);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, B_in, C, DFT_K, 1.0, P_in, DFT_K, X_target,
                DFT_K, 0.0, PX_full, C);

    while (m > B_out) {
        double G_inv[160 * 160];
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < m; ++i) {
                G_inv[j * m + i] = G_full[active[i] * B_in + active[j]] + ((i == j) ? 1e-7 : 0.0);
            }
        }
        char           uplo = 'U';
        dft_lapack_int nm = m, info = 0;
        dpotrf_(&uplo, &nm, G_inv, &nm, &info);
        dpotri_(&uplo, &nm, G_inv, &nm, &info);
        for (int j = 0; j < m; ++j) {
            for (int i = j + 1; i < m; ++i) G_inv[j * m + i] = G_inv[i * m + j];
        }

        int    best_drop  = 0;
        double best_score = 1e300;
        for (int i = 0; i < m; ++i) {
            double wn2 = 0.0;
            for (int c = 0; c < C; ++c) {
                double wc = 0.0;
                for (int k = 0; k < m; ++k) wc += G_inv[k * m + i] * PX_full[active[k] * C + c];
                wn2 += wc * wc;
            }
            double diag = fmax(1e-12, G_inv[i * m + i]);
            double sc   = wn2 / diag;
            if (sc < best_score) {
                best_score = sc;
                best_drop  = i;
            }
        }
        for (int i = best_drop; i + 1 < m; ++i) active[i] = active[i + 1];
        --m;
    }

    for (int i = 0; i < B_out; ++i) {
        memcpy(U_out + (size_t)i * DFT_N, U_in + (size_t)active[i] * DFT_N, DFT_N * sizeof(double));
        normalize_u_vec(U_out + (size_t)i * DFT_N);
    }

    double G_b[128 * 128], RHS_col[128 * 128];
    for (int j = 0; j < B_out; ++j) {
        for (int i = 0; i < B_out; ++i) {
            G_b[j * B_out + i] = G_full[active[i] * B_in + active[j]] + ((i == j) ? 1e-6 : 0.0);
        }
    }
    for (int c = 0; c < C; ++c) {
        for (int i = 0; i < B_out; ++i) RHS_col[c * B_out + i] = PX_full[active[i] * C + c];
    }
    char           uplo = 'U';
    dft_lapack_int nb = B_out, nc = C, info = 0;
    dpotrf_(&uplo, &nb, G_b, &nb, &info);
    dpotrs_(&uplo, &nb, &nc, G_b, &nb, RHS_col, &nb, &info);
    for (int b = 0; b < B_out; ++b) {
        for (int c = 0; c < C; ++c) w_out[b * C + c] = RHS_col[c * B_out + b];
    }
    free(P_in);
}

/* Arbitrary-stage 4-Tensor Spectral Agglomeration: (R_src, B_src, C_src) -> (R_dst, B_dst, C_dst) */
static void compress_blocks_4tensor(const double *U_src, const double *w_src, int R_src, int B_src,
                                    int C_src, int R_dst, int B_dst, int C_dst, double *U_dst,
                                    double *w_dst) {
    int cur_R = R_src;
    int B_max = (B_src > B_dst) ? B_src : B_dst;
    if (B_max > 76) B_max = 76;
    int C_cur = (C_dst < C_src) ? C_dst : C_src;

    double *U_list = (double *)malloc((size_t)cur_R * B_max * DFT_N * sizeof(double));
    double *X_list = (double *)malloc((size_t)cur_R * C_cur * DFT_K * sizeof(double));
    int    *B_cnt  = (int *)malloc((size_t)cur_R * sizeof(int));
    double *P_tmp  = (double *)malloc((size_t)(2 * B_max) * DFT_K * sizeof(double));

    for (int r = 0; r < cur_R; ++r) {
        const double *Ur   = U_src + (size_t)r * B_src * DFT_N;
        const double *wr   = w_src + (size_t)r * B_src * C_src;
        double       *Xr_s = (double *)malloc((size_t)C_src * DFT_K * sizeof(double));
        compute_block_X(Ur, wr, B_src, C_src, P_tmp, Xr_s);

        if (C_cur < C_src) {
            double G_s[64 * 64];
            cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, C_src, C_src, DFT_K, 1.0, Xr_s,
                        DFT_K, Xr_s, DFT_K, 0.0, G_s, C_src);
            double         ev[64], wk[2048];
            dft_lapack_int ns = C_src, lw = 2048, info = 0;
            char           jz = 'V', up = 'U';
            dsyev_(&jz, &up, &ns, G_s, &ns, ev, wk, &lw, &info);
            double V_top[64 * 64];
            for (int c = 0; c < C_cur; ++c) {
                int idx = C_src - 1 - c;
                for (int i = 0; i < C_src; ++i) V_top[c * C_src + i] = G_s[idx * C_src + i];
            }
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, DFT_K, C_cur, C_src, 1.0, Xr_s,
                        DFT_K, V_top, C_src, 0.0, X_list + (size_t)r * C_cur * DFT_K, DFT_K);
        } else {
            memcpy(X_list + (size_t)r * C_cur * DFT_K, Xr_s,
                   (size_t)C_cur * DFT_K * sizeof(double));
        }
        free(Xr_s);

        int b_keep = B_src;
        if (b_keep > B_max) b_keep = B_max;
        memcpy(U_list + (size_t)r * B_max * DFT_N, Ur, (size_t)b_keep * DFT_N * sizeof(double));
        B_cnt[r] = b_keep;
    }

    while (cur_R > R_dst) {
        double norms4[65];
        for (int r = 0; r < cur_R; ++r) {
            const double *Xr = X_list + (size_t)r * C_cur * DFT_K;
            double        G[64 * 64];
            cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, C_cur, C_cur, DFT_K, 1.0, Xr,
                        DFT_K, Xr, DFT_K, 0.0, G, C_cur);
            double s2 = 0.0;
            for (int i = 0; i < C_cur * C_cur; ++i) s2 += G[i] * G[i];
            norms4[r] = s2;
        }

        int ord[65];
        for (int r = 0; r < cur_R; ++r) ord[r] = r;
        for (int i = 0; i < cur_R; ++i) {
            for (int j = i + 1; j < cur_R; ++j) {
                if (norms4[ord[j]] < norms4[ord[i]]) {
                    int t  = ord[i];
                    ord[i] = ord[j];
                    ord[j] = t;
                }
            }
        }
        int pool = (cur_R < 28) ? cur_R : 28;

        int     bi = ord[0], bj = ord[1];
        double  best_loss4 = 1e300;
        double *X_cat      = (double *)malloc(128 * DFT_K * sizeof(double));
        double  G2[128 * 128], ev2[128], wk2[4096];

        for (int a = 0; a < pool; ++a) {
            int i = ord[a];
            for (int b = a + 1; b < pool; ++b) {
                int j = ord[b];
                memcpy(X_cat, X_list + (size_t)i * C_cur * DFT_K,
                       (size_t)C_cur * DFT_K * sizeof(double));
                memcpy(X_cat + (size_t)C_cur * DFT_K, X_list + (size_t)j * C_cur * DFT_K,
                       (size_t)C_cur * DFT_K * sizeof(double));
                int c2 = 2 * C_cur;
                cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, c2, c2, DFT_K, 1.0, X_cat,
                            DFT_K, X_cat, DFT_K, 0.0, G2, c2);
                dft_lapack_int nc2 = c2, lw2 = 4096, info = 0;
                char           jz = 'N', up = 'U';
                dsyev_(&jz, &up, &nc2, G2, &nc2, ev2, wk2, &lw2, &info);
                double tail_sq = 0.0;
                for (int k = 0; k < C_cur; ++k) {
                    if (ev2[k] > 0.0) tail_sq += ev2[k] * ev2[k];
                }
                if (tail_sq < best_loss4) {
                    best_loss4 = tail_sq;
                    bi         = (i < j) ? i : j;
                    bj         = (i < j) ? j : i;
                }
            }
        }

        int c2 = 2 * C_cur;
        memcpy(X_cat, X_list + (size_t)bi * C_cur * DFT_K, (size_t)C_cur * DFT_K * sizeof(double));
        memcpy(X_cat + (size_t)C_cur * DFT_K, X_list + (size_t)bj * C_cur * DFT_K,
               (size_t)C_cur * DFT_K * sizeof(double));
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, c2, c2, DFT_K, 1.0, X_cat, DFT_K,
                    X_cat, DFT_K, 0.0, G2, c2);
        dft_lapack_int nc2 = c2, lw2 = 4096, info = 0;
        char           jz = 'V', up = 'U';
        dsyev_(&jz, &up, &nc2, G2, &nc2, ev2, wk2, &lw2, &info);

        double V_top[128 * 64];
        for (int c = 0; c < C_cur; ++c) {
            int idx = c2 - 1 - c;
            for (int l = 0; l < c2; ++l) V_top[c * c2 + l] = G2[idx * c2 + l];
        }
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, DFT_K, C_cur, c2, 1.0, X_cat, DFT_K,
                    V_top, c2, 0.0, X_list + (size_t)bi * C_cur * DFT_K, DFT_K);

        int    btot = B_cnt[bi] + B_cnt[bj];
        double U_cat[128 * DFT_N];
        memcpy(U_cat, U_list + (size_t)bi * B_max * DFT_N,
               (size_t)B_cnt[bi] * DFT_N * sizeof(double));
        memcpy(U_cat + (size_t)B_cnt[bi] * DFT_N, U_list + (size_t)bj * B_max * DFT_N,
               (size_t)B_cnt[bj] * DFT_N * sizeof(double));
        int b_keep = (B_dst + 12 < B_max) ? (B_dst + 12) : B_max;
        if (b_keep > btot) b_keep = btot;
        double w_dummy[128 * 128];
        prune_and_fit_block_to_X(U_cat, btot, b_keep, C_cur, X_list + (size_t)bi * C_cur * DFT_K,
                                 U_list + (size_t)bi * B_max * DFT_N, w_dummy);
        B_cnt[bi] = b_keep;

        for (int r = bj; r + 1 < cur_R; ++r) {
            memcpy(U_list + (size_t)r * B_max * DFT_N, U_list + (size_t)(r + 1) * B_max * DFT_N,
                   (size_t)B_max * DFT_N * sizeof(double));
            memcpy(X_list + (size_t)r * C_cur * DFT_K, X_list + (size_t)(r + 1) * C_cur * DFT_K,
                   (size_t)C_cur * DFT_K * sizeof(double));
            B_cnt[r] = B_cnt[r + 1];
        }
        free(X_cat);
        --cur_R;
    }

    for (int r = 0; r < R_dst; ++r) {
        memset(w_dst + (size_t)r * B_dst * C_dst, 0, (size_t)B_dst * C_dst * sizeof(double));
        double w_fit[128 * 128];
        prune_and_fit_block_to_X(U_list + (size_t)r * B_max * DFT_N, B_cnt[r], B_dst, C_cur,
                                 X_list + (size_t)r * C_cur * DFT_K,
                                 U_dst + (size_t)r * B_dst * DFT_N, w_fit);
        for (int b = 0; b < B_dst; ++b) {
            for (int c = 0; c < C_cur; ++c) {
                w_dst[(r * B_dst + b) * C_dst + c] = w_fit[b * C_cur + c];
            }
        }
    }

    free(P_tmp);
    free(B_cnt);
    free(X_list);
    free(U_list);
}

/* ============================================================================
 * Pass 2: SEARU Multi-Chain MSB-0 Givens Tree Synthesis
 * ============================================================================ */

static void compute_msb0_mask(const int64_t edges[DFT_E][2], int msb0_mask[DFT_E]) {
    for (int t = 0; t < DFT_E; ++t) {
        int ch  = (int)edges[t][1];
        int has = 0;
        for (int s = t + 1; s < DFT_E; ++s) {
            if (edges[s][0] == ch) {
                has = 1;
                break;
            }
        }
        msb0_mask[t] = has;
    }
    /* Edge DFT_E - 1 is ALWAYS a terminal leaf edge: global sign gauge u -> -u forces its angle into [0, pi)! */
    msb0_mask[DFT_E - 1] = 1;
}

static void encode_u_on_tree_msb0(const double *u_in, const int64_t edges[DFT_E][2],
                                  double alpha_out[DFT_E]) {
    double v[DFT_N];
    memcpy(v, u_in, DFT_N * sizeof(double));
    normalize_u_vec(v);

    /* Choose global sign u -> +-u so final edge (p_last, ch_last) has atan2(v[ch_last], v[p_last]) in [0, pi)! */
    int p_last  = (int)edges[DFT_E - 1][0];
    int ch_last = (int)edges[DFT_E - 1][1];
    if (v[ch_last] < 0.0 || (v[ch_last] == 0.0 && v[p_last] < 0.0)) {
        for (int i = 0; i < DFT_N; ++i) v[i] = -v[i];
    }

    for (int t = DFT_E - 1; t >= 0; --t) {
        int    p  = (int)edges[t][0];
        int    ch = (int)edges[t][1];
        double vp = v[p];
        double vc = v[ch];
        double r  = hypot(vp, vc);
        if (r < 1e-15) {
            alpha_out[t] = 0.0;
            v[p]         = 0.0;
            v[ch]        = 0.0;
        } else {
            double th    = atan2(vc, vp);
            alpha_out[t] = th / kPi; /* in [0, 1] whenever msb0_mask[t] == 1! */
            v[p]         = r;
            v[ch]        = 0.0;
        }
    }
}

static void decode_words_on_tree(const int64_t words[DFT_E], const int64_t edges[DFT_E][2],
                                 const double *cos_lut, const double *sin_lut, double u_out[DFT_N]) {
    u_out[0] = 1.0;
    for (int i = 1; i < DFT_N; ++i) u_out[i] = 0.0;
    for (int t = 0; t < DFT_E; ++t) {
        int     p  = (int)edges[t][0];
        int     ch = (int)edges[t][1];
        int64_t w  = words[t];
        double  vp = u_out[p];
        u_out[p]   = cos_lut[w] * vp;
        u_out[ch]  = sin_lut[w] * vp;
    }
}

/* Build K_ch-chain DFS tree (`K_ch in {6,7,8,9}`) with round-robin depth placement from SEARU-permuted orbitals!
 * High-variance orbitals sit near root `0` (depths 1..2 across all chains), AND `75 - K_ch + 1` edges
 * (`67..70` out of `75`!) have `MSB == 0` across ALL rows (`CleanPool` free |0> qubits in `synth_ele_128`)! */
static void build_multichain_tree_from_spec(const int orb_order[DFT_E], const uint64_t canvas[4],
                                            int64_t edges_out[DFT_E][2]) {
    int perm[DFT_E];
    memcpy(perm, orb_order, DFT_E * sizeof(int));
    for (int b = 0; b < 60; ++b) {
        if ((canvas[0] >> b) & 1ULL) {
            int i = b;
            if (i + 1 < DFT_E) {
                int tmp     = perm[i];
                perm[i]     = perm[i + 1];
                perm[i + 1] = tmp;
            }
        }
    }
    const int k_choices[4] = {6, 7, 8, 9};
    int       n_chains     = k_choices[canvas[1] & 3ULL];

    /* Round-robin assign `perm[0..74]` across `n_chains` so heaviest orbitals are at depth 1 of each chain! */
    int chains[12][16];
    int ch_len[12];
    memset(ch_len, 0, sizeof(ch_len));
    for (int idx = 0; idx < DFT_E; ++idx) {
        int c                  = idx % n_chains;
        chains[c][ch_len[c]++] = perm[idx];
    }

    int pos = 0;
    for (int c = 0; c < n_chains; ++c) {
        int prev = 0;
        for (int d = 0; d < ch_len[c]; ++d) {
            int cur           = chains[c][d];
            edges_out[pos][0] = prev;
            edges_out[pos][1] = cur;
            prev              = cur;
            ++pos;
        }
    }
}

static void synthesize_givens_tree_searu(const double *U, const double *w, int R, int B, int C,
                                         int brot, int64_t best_edges_out[DFT_E][2]) {
    double m_orb[DFT_N];
    memset(m_orb, 0, sizeof(m_orb));
    int    M      = R * B;
    double wn2[4096];
    double wn_tot = 0.0;
    for (int rb = 0; rb < M; ++rb) {
        const double *wb = w + (size_t)rb * C;
        double        s2 = 0.0;
        for (int c = 0; c < C; ++c) s2 += wb[c] * wb[c];
        wn2[rb] = s2;
        wn_tot += s2;
        const double *ub = U + (size_t)rb * DFT_N;
        for (int i = 0; i < DFT_N; ++i) m_orb[i] += s2 * ub[i] * ub[i];
    }

    /* Note: orbital 0 is root. Sort 1..75 descending by marginal energy m_orb[i]! */
    int orb_order[DFT_E];
    for (int i = 0; i < DFT_E; ++i) orb_order[i] = i + 1;
    for (int i = 0; i < DFT_E; ++i) {
        for (int j = i + 1; j < DFT_E; ++j) {
            if (m_orb[orb_order[j]] > m_orb[orb_order[i]]) {
                int t        = orb_order[i];
                orb_order[i] = orb_order[j];
                orb_order[j] = t;
            }
        }
    }

    int     mod     = 1 << brot;
    int     half    = 1 << (brot - 1);
    double  scale   = (double)half;
    double *cos_lut = (double *)malloc((size_t)mod * sizeof(double));
    double *sin_lut = (double *)malloc((size_t)mod * sizeof(double));
    for (int k = 0; k < mod; ++k) {
        double th  = kPi * (double)k / scale;
        cos_lut[k] = cos(th);
        sin_lut[k] = sin(th);
    }

    int top_rb[4096], sample_idx[128];
    for (int i = 0; i < M; ++i) top_rb[i] = i;
    int n_samp = (M < 128) ? M : 128;
    for (int i = 0; i < n_samp; ++i) {
        for (int j = i + 1; j < M; ++j) {
            if (wn2[top_rb[j]] > wn2[top_rb[i]]) {
                int t     = top_rb[i];
                top_rb[i] = top_rb[j];
                top_rb[j] = t;
            }
        }
        sample_idx[i] = top_rb[i];
    }

    searu_result sr;
    uint64_t     init_c[4] = {0ULL, 2ULL, 0ULL, 0ULL}; /* start at 8-chain */
    uint64_t     pin_c[4]  = {0ULL, ~UINT64_C(3), ~UINT64_C(0), ~UINT64_C(0)};

    SEARU_SOLVE_INLINE(&sr, 66U, init_c, pin_c, UINT64_C(0xDF7C505), 300U, cand_c, frust_c, e_out, {
        (void)frust_c;
        int64_t cand_edges[DFT_E][2];
        int     msb0[DFT_E];
        build_multichain_tree_from_spec(orb_order, cand_c, cand_edges);
        compute_msb0_mask(cand_edges, msb0);
        double err_acc = 0.0;
        for (int s = 0; s < n_samp; ++s) {
            int           rb = sample_idx[s];
            const double *u  = U + (size_t)rb * DFT_N;
            double        ang[DFT_E];
            int64_t       wd[DFT_E];
            double        uq[DFT_N];
            encode_u_on_tree_msb0(u, cand_edges, ang);
            for (int t = 0; t < DFT_E; ++t) {
                if (msb0[t]) {
                    int64_t k = (int64_t)llround(ang[t] * scale);
                    if (k < 0) k = 0;
                    if (k >= half) k = half - 1;
                    wd[t] = k;
                } else {
                    int64_t k = (int64_t)llround(ang[t] * scale) % mod;
                    if (k < 0) k += mod;
                    wd[t] = k;
                }
            }
            decode_words_on_tree(wd, cand_edges, cos_lut, sin_lut, uq);
            double dot = 0.0;
            for (int i = 0; i < DFT_N; ++i) dot += u[i] * uq[i];
            err_acc += (wn2[rb] / wn_tot) * (1.0 - dot * dot);
        }
        e_out = (uint64_t)llround(err_acc * 1e11) + UINT64_C(1);
    });

    build_multichain_tree_from_spec(orb_order, sr.canvas, best_edges_out);
    int msb0[DFT_E], n_msb0 = 0;
    compute_msb0_mask(best_edges_out, msb0);
    for (int t = 0; t < DFT_E; ++t) n_msb0 += msb0[t];
    printf("  [Pass2 SEARU Tree brot=%d] init_e=%llu -> best_e=%llu | msb0_edges=%d/75\n", brot,
           (unsigned long long)sr.init_energy, (unsigned long long)sr.best_energy, n_msb0);
    fflush(stdout);

    free(sin_lut);
    free(cos_lut);
}

/* ============================================================================
 * Pass 2.5: Multi-Threaded Angle-Space QAT Homotopy L-BFGS + FLOW 2-Jet
 * ============================================================================ */

typedef struct {
    const dft_ref_t *ref;
    const int64_t (*edges)[2];
    const int *msb0_mask;
    int        R, B, C, brot;
    double     w_scale;
    double     gamma_grid;
    double     msb_bar_weight;
    double    *U_buf;     /* M x N */
    double    *P_buf;     /* M x K */
    double    *vp_buf;    /* M x E */
    double    *gg_buf;    /* M x E */
    double    *X_buf;     /* RC x K */
    double    *G_XX_buf;  /* RC x RC */
    double    *XY_buf;    /* RC x K_pos */
    double    *grad_X;    /* RC x K */
    double    *grad_P;    /* M x K */
    double    *wn2_fixed; /* M */
} qat_ctx_t;

static double fg_qat_alpha_w(const double *x, double *g, void *vctx) {
    qat_ctx_t       *c        = (qat_ctx_t *)vctx;
    const dft_ref_t *ref      = c->ref;
    const int        R        = c->R;
    const int        B        = c->B;
    const int        C        = c->C;
    const int        M        = R * B;
    const int        RC       = R * C;
    const int        a_sz     = M * DFT_E;
    const double     ws       = c->w_scale;
    const double     scale    = (double)(1 << (c->brot - 1));
    const double     pi_scale = kPi * scale;
    const double     hi_msb0  = (scale - 1.0) / scale;

    const double *alpha   = x;
    const double *wt_in   = x + a_sz;
    double       *grad_a  = g;
    double       *grad_wt = g + a_sz;

    const int chunk      = 64;
    const int n_chunks   = (M + chunk - 1) / chunk;
    double   *chunk_pen  = (double *)calloc((size_t)n_chunks, sizeof(double));

    /* 1. Multi-threaded forward Givens tree + grid-well penalty + MSB0 box barrier */
    dispatch_apply((size_t)n_chunks, dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0),
                   ^(size_t cid) {
                       int    lo  = (int)cid * chunk;
                       int    hi  = (lo + chunk < M) ? (lo + chunk) : M;
                       double loc = 0.0;
                       for (int m = lo; m < hi; ++m) {
                           const double *a  = alpha + (size_t)m * DFT_E;
                           double       *u  = c->U_buf + (size_t)m * DFT_N;
                           double       *vp = c->vp_buf + (size_t)m * DFT_E;
                           double       *gg = c->gg_buf + (size_t)m * DFT_E;

                           u[0] = 1.0;
                           for (int j = 1; j < DFT_N; ++j) u[j] = 0.0;

                           double wg = c->gamma_grid * c->wn2_fixed[m];
                           for (int t = 0; t < DFT_E; ++t) {
                               int    p      = (int)c->edges[t][0];
                               int    ch     = (int)c->edges[t][1];
                               double th     = kPi * a[t];
                               double cs     = cos(th);
                               double sn     = sin(th);
                               double v_prev = u[p];
                               vp[t]         = v_prev;
                               u[p]          = cs * v_prev;
                               u[ch]         = sn * v_prev;

                               double g_pen = 0.0;
                               if (wg > 0.0) {
                                   double z  = pi_scale * a[t];
                                   double sz = sin(z);
                                   loc += wg * sz * sz;
                                   g_pen += wg * pi_scale * sin(2.0 * z);
                               }
                               if (c->msb0_mask[t]) {
                                   if (a[t] < 0.0) {
                                       double d0 = a[t];
                                       loc += c->msb_bar_weight * d0 * d0;
                                       g_pen += 2.0 * c->msb_bar_weight * d0;
                                   } else if (a[t] > hi_msb0) {
                                       double d1 = a[t] - hi_msb0;
                                       loc += c->msb_bar_weight * d1 * d1;
                                       g_pen += 2.0 * c->msb_bar_weight * d1;
                                   }
                               }
                               gg[t] = g_pen;
                           }
                           pack_u_vec(u, c->P_buf + (size_t)m * DFT_K);
                       }
                       chunk_pen[cid] = loc;
                   });

    double pen_tot = 0.0;
    for (int i = 0; i < n_chunks; ++i) pen_tot += chunk_pen[i];
    free(chunk_pen);

    /* 2. Low-Rank 4-Tensor Contraction */
    for (int r = 0; r < R; ++r) {
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, C, DFT_K, B, ws,
                    wt_in + (size_t)r * B * C, C, c->P_buf + (size_t)r * B * DFT_K, DFT_K, 0.0,
                    c->X_buf + (size_t)r * C * DFT_K, DFT_K);
    }
    cblas_dsyrk(CblasRowMajor, CblasUpper, CblasNoTrans, RC, DFT_K, 1.0, c->X_buf, DFT_K, 0.0,
                c->G_XX_buf, RC);
    for (int i = 0; i < RC; ++i) {
        for (int j = 0; j < i; ++j) c->G_XX_buf[(size_t)i * RC + j] = c->G_XX_buf[(size_t)j * RC + i];
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, RC, ref->K_pos, DFT_K, 1.0, c->X_buf,
                DFT_K, ref->Y_target, ref->K_pos, 0.0, c->XY_buf, ref->K_pos);

    double g2 = 0.0, xy2 = 0.0;
    for (size_t i = 0; i < (size_t)RC * RC; ++i) g2 += c->G_XX_buf[i] * c->G_XX_buf[i];
    for (size_t i = 0; i < (size_t)RC * ref->K_pos; ++i) xy2 += c->XY_buf[i] * c->XY_buf[i];
    double loss = 0.25 * fmax(0.0, g2 - 2.0 * xy2 + ref->Y_norm2) + pen_tot;

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, RC, DFT_K, RC, 1.0, c->G_XX_buf, RC,
                c->X_buf, DFT_K, 0.0, c->grad_X, DFT_K);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, RC, DFT_K, ref->K_pos, -1.0, c->XY_buf,
                ref->K_pos, ref->Y_target, ref->K_pos, 1.0, c->grad_X, DFT_K);

    for (int r = 0; r < R; ++r) {
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, B, C, DFT_K, ws,
                    c->P_buf + (size_t)r * B * DFT_K, DFT_K, c->grad_X + (size_t)r * C * DFT_K,
                    DFT_K, 0.0, grad_wt + (size_t)r * B * C, C);
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, B, DFT_K, C, ws,
                    wt_in + (size_t)r * B * C, C, c->grad_X + (size_t)r * C * DFT_K, DFT_K, 0.0,
                    c->grad_P + (size_t)r * B * DFT_K, DFT_K);
    }

    /* 3. Multi-threaded exact reverse-mode Givens tree adjoint */
    dispatch_apply((size_t)n_chunks, dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0),
                   ^(size_t cid) {
                       int    lo = (int)cid * chunk;
                       int    hi = (lo + chunk < M) ? (lo + chunk) : M;
                       double gu[DFT_N];
                       for (int m = lo; m < hi; ++m) {
                           grad_pack_u_vec(c->grad_P + (size_t)m * DFT_K,
                                           c->U_buf + (size_t)m * DFT_N, gu);
                           const double *a  = alpha + (size_t)m * DFT_E;
                           const double *vp = c->vp_buf + (size_t)m * DFT_E;
                           const double *gg = c->gg_buf + (size_t)m * DFT_E;
                           double       *ga = grad_a + (size_t)m * DFT_E;

                           for (int t = DFT_E - 1; t >= 0; --t) {
                               int    p  = (int)c->edges[t][0];
                               int    ch = (int)c->edges[t][1];
                               double th = kPi * a[t];
                               double cs = cos(th);
                               double sn = sin(th);
                               double gp = gu[p];
                               double gc = gu[ch];
                               ga[t]     = kPi * vp[t] * (-sn * gp + cs * gc) + gg[t];
                               gu[p]     = cs * gp + sn * gc;
                               gu[ch]    = 0.0;
                           }
                       }
                   });

    return loss;
}

static void run_qat_angle_lbfgs_c17(const dft_ref_t *ref, int R, int B, int C, int brot,
                                    const int64_t edges[DFT_E][2], const int msb0_mask[DFT_E],
                                    double *U, double *w, double *alpha_out) {
    const int    M     = R * B;
    const int    RC    = R * C;
    const int    a_sz  = M * DFT_E;
    const int    w_sz  = M * C;
    const int    ndim  = a_sz + w_sz;
    const double ws    = 0.05;
    const double scale = (double)(1 << (brot - 1));
    const double hi_m0 = (scale - 1.0) / scale;

    double *x = (double *)malloc((size_t)ndim * sizeof(double));
    for (int m = 0; m < M; ++m) {
        double *am = x + (size_t)m * DFT_E;
        encode_u_on_tree_msb0(U + (size_t)m * DFT_N, edges, am);
        for (int t = 0; t < DFT_E; ++t) {
            if (msb0_mask[t]) {
                if (am[t] < 0.0) am[t] = 0.0;
                if (am[t] > hi_m0) am[t] = hi_m0;
            }
        }
    }
    for (int i = 0; i < w_sz; ++i) x[a_sz + i] = w[i] / ws;

    qat_ctx_t ctx = {
        .ref            = ref,
        .edges          = edges,
        .msb0_mask      = msb0_mask,
        .R              = R,
        .B              = B,
        .C              = C,
        .brot           = brot,
        .w_scale        = ws,
        .gamma_grid     = 0.0,
        .msb_bar_weight = 2.0,
        .U_buf          = (double *)malloc((size_t)M * DFT_N * sizeof(double)),
        .P_buf          = (double *)malloc((size_t)M * DFT_K * sizeof(double)),
        .vp_buf         = (double *)malloc((size_t)M * DFT_E * sizeof(double)),
        .gg_buf         = (double *)malloc((size_t)M * DFT_E * sizeof(double)),
        .X_buf          = (double *)malloc((size_t)RC * DFT_K * sizeof(double)),
        .G_XX_buf       = (double *)malloc((size_t)RC * RC * sizeof(double)),
        .XY_buf         = (double *)malloc((size_t)RC * ref->K_pos * sizeof(double)),
        .grad_X         = (double *)malloc((size_t)RC * DFT_K * sizeof(double)),
        .grad_P         = (double *)malloc((size_t)M * DFT_K * sizeof(double)),
        .wn2_fixed      = (double *)malloc((size_t)M * sizeof(double)),
    };

    const double gammas[6] = {0.0, 2e-5, 1e-4, 5e-4, 2.5e-3, 1.0e-2};
    const int    iters[6]  = {180, 180, 200, 220, 220, 220};

    for (int st = 0; st < 6; ++st) {
        /* Refresh per-vector weight norm ||w_m||^2 */
        for (int m = 0; m < M; ++m) {
            double s2 = 0.0;
            for (int c = 0; c < C; ++c) {
                double wv = x[a_sz + m * C + c] * ws;
                s2 += wv * wv;
            }
            ctx.wn2_fixed[m] = s2;
        }
        ctx.gamma_grid = gammas[st];
        char tag[64];
        snprintf(tag, sizeof(tag), "QAT-b%d g=%.0e", brot, gammas[st]);
        run_lbfgs_flow2jet(ndim, x, iters[st], 16, fg_qat_alpha_w, &ctx, tag, 0);
    }

    memcpy(alpha_out, x, (size_t)a_sz * sizeof(double));
    for (int i = 0; i < w_sz; ++i) w[i] = x[a_sz + i] * ws;

    free(ctx.wn2_fixed);
    free(ctx.grad_P);
    free(ctx.grad_X);
    free(ctx.XY_buf);
    free(ctx.G_XX_buf);
    free(ctx.X_buf);
    free(ctx.gg_buf);
    free(ctx.vp_buf);
    free(ctx.P_buf);
    free(ctx.U_buf);
    free(x);
}

/* ============================================================================
 * Pass 3: Live-Residual 4-Tensor Prefix-CD + SEARU Sparse-k Coupled Crucible
 * ============================================================================ */

static inline void eval_two_packed(const double *__restrict__ h_pack,
                                   const double *__restrict__ v_pack, const double *__restrict__ u,
                                   double *out_h, double *out_v) {
    double tot_h = 0.0, tot_v = 0.0;
    int    k = 0;
    for (int r = 0; r < DFT_N; ++r) {
        double acc_h = 0.0, acc_v = 0.0;
        for (int c = 0; c < r; ++c) {
            double uc = u[c];
            acc_h += h_pack[k] * uc;
            acc_v += v_pack[k] * uc;
            ++k;
        }
        double ur    = u[r];
        double ur_s2 = ur * kSqrt2;
        double ur2   = ur * ur;
        tot_h += acc_h * ur_s2 + h_pack[k] * ur2;
        tot_v += acc_v * ur_s2 + v_pack[k] * ur2;
        ++k;
    }
    *out_h = tot_h;
    *out_v = tot_v;
}

static inline double eval_4t_surrogate(const double *__restrict__ u_cand,
                                       const double *__restrict__ u0,
                                       const double *__restrict__ h_pack,
                                       const double *__restrict__ v_pack, double h0, double beta0,
                                       double v_norm2, double gamma) {
    double dot0 = 0.0;
    for (int j = 0; j < DFT_N; ++j) dot0 += u0[j] * u_cand[j];
    double alpha   = dot0 * dot0;
    double one_m_a = 1.0 - alpha;
    double hp, beta;
    eval_two_packed(h_pack, v_pack, u_cand, &hp, &beta);
    double db = beta - beta0;
    return (hp - h0) + 4.0 * v_norm2 * one_m_a + 2.0 * db * db +
           2.0 * gamma * gamma * (1.0 - alpha * alpha) + 4.0 * gamma * (beta + beta0) * one_m_a;
}

static inline int64_t step_angle_word(int64_t w, int delta, int is_msb0, int half, int mod) {
    if (is_msb0) {
        int64_t nw = w + delta;
        if (nw < 0) nw = 0;
        if (nw >= half) nw = half - 1;
        return nw;
    }
    int64_t nw = (w + delta) % mod;
    if (nw < 0) nw += mod;
    return nw;
}

static bool optimize_vec_searu_4tensor(int64_t wb_words[DFT_E], const double u0[DFT_N],
                                       const double h_pack[DFT_K], const double v_pack[DFT_K],
                                       double v_norm2, double gamma,
                                       const int64_t edges[DFT_E][2], const int msb0_mask[DFT_E],
                                       const double *cos_lut, const double *sin_lut, int brot,
                                       int rad, int searu_steps, uint64_t *rng_io,
                                       double u_out[DFT_N], double p_out[DFT_K]) {
    const int mod  = 1 << brot;
    const int half = 1 << (brot - 1);
    double    h0, beta0;
    eval_two_packed(h_pack, v_pack, u0, &h0, &beta0);

    double best_df = 0.0;
    bool   changed = false;
    double u_pref[DFT_N], u_buf[DFT_N];
    double vp_sens[DFT_E];

    for (int pass = 0; pass < 2; ++pass) {
        bool imp  = false;
        u_pref[0] = 1.0;
        for (int j = 1; j < DFT_N; ++j) u_pref[j] = 0.0;

        for (int t = 0; t < DFT_E; ++t) {
            int     m0   = msb0_mask[t];
            int64_t orig = wb_words[t];
            int64_t bw   = orig;
            vp_sens[t]   = fabs(u_pref[edges[t][0]]);

            if (vp_sens[t] > 1e-7) {
                int r_eff = (rad > half) ? half : rad;
                for (int d = -r_eff; d <= r_eff; ++d) {
                    if (d == 0) continue;
                    int64_t cand = step_angle_word(orig, d, m0, half, mod);
                    if (cand == orig) continue;
                    memcpy(u_buf, u_pref, DFT_N * sizeof(double));
                    int    p0 = (int)edges[t][0], c0 = (int)edges[t][1];
                    double vp = u_buf[p0];
                    u_buf[p0] = cos_lut[cand] * vp;
                    u_buf[c0] = sin_lut[cand] * vp;
                    for (int s = t + 1; s < DFT_E; ++s) {
                        int     ps = (int)edges[s][0], cs = (int)edges[s][1];
                        int64_t ws = wb_words[s];
                        double  vs = u_buf[ps];
                        u_buf[ps]  = cos_lut[ws] * vs;
                        u_buf[cs]  = sin_lut[ws] * vs;
                    }
                    double df =
                        eval_4t_surrogate(u_buf, u0, h_pack, v_pack, h0, beta0, v_norm2, gamma);
                    if (df < best_df - 1e-18) {
                        best_df = df;
                        bw      = cand;
                        imp     = true;
                        changed = true;
                    }
                }
                wb_words[t] = bw;

                int p_t = (int)edges[t][0], ch_t = (int)edges[t][1];
                int part[2] = {-1, -1};
                for (int s = t + 1; s < DFT_E; ++s) {
                    if (part[0] < 0 && edges[s][0] == ch_t) part[0] = s;
                    if (part[1] < 0 && edges[s][0] == p_t) part[1] = s;
                    if (part[0] >= 0 && part[1] >= 0) break;
                }
                for (int pi = 0; pi < 2; ++pi) {
                    int t2 = part[pi];
                    if (t2 < 0) continue;
                    int     m2  = msb0_mask[t2];
                    int64_t w1c = wb_words[t], w2c = wb_words[t2];
                    int64_t bw1 = w1c, bw2 = w2c;
                    for (int d1 = -2; d1 <= 2; ++d1) {
                        if (d1 == 0) continue;
                        int64_t c1 = step_angle_word(w1c, d1, m0, half, mod);
                        if (c1 == w1c) continue;
                        for (int d2 = -2; d2 <= 2; ++d2) {
                            if (d2 == 0) continue;
                            int64_t c2 = step_angle_word(w2c, d2, m2, half, mod);
                            if (c2 == w2c) continue;
                            wb_words[t]  = c1;
                            wb_words[t2] = c2;
                            memcpy(u_buf, u_pref, DFT_N * sizeof(double));
                            for (int s = t; s < DFT_E; ++s) {
                                int     ps = (int)edges[s][0], cs = (int)edges[s][1];
                                int64_t ws = wb_words[s];
                                double  vs = u_buf[ps];
                                u_buf[ps]  = cos_lut[ws] * vs;
                                u_buf[cs]  = sin_lut[ws] * vs;
                            }
                            double df = eval_4t_surrogate(u_buf, u0, h_pack, v_pack, h0, beta0,
                                                          v_norm2, gamma);
                            if (df < best_df - 1e-18) {
                                best_df = df;
                                bw1     = c1;
                                bw2     = c2;
                                imp     = true;
                                changed = true;
                            }
                        }
                    }
                    wb_words[t]  = bw1;
                    wb_words[t2] = bw2;
                }
            }

            int     pt = (int)edges[t][0], ct = (int)edges[t][1];
            int64_t wt = wb_words[t];
            double  vt = u_pref[pt];
            u_pref[pt] = cos_lut[wt] * vt;
            u_pref[ct] = sin_lut[wt] * vt;
        }
        if (!imp) break;
    }

    if (searu_steps > 0) {
        int active_e[DFT_E];
        int n_act = 0;
        for (int t = 0; t < DFT_E; ++t) {
            if (vp_sens[t] > 0.03) active_e[n_act++] = t;
        }
        if (n_act >= 2) {
            int64_t cur_w[DFT_E], best_w[DFT_E], cand_w[DFT_E];
            memcpy(cur_w, wb_words, sizeof(cur_w));
            memcpy(best_w, wb_words, sizeof(best_w));
            double   cur_df   = best_df;
            uint32_t stag     = 0U;
            uint64_t temp_q12 = UINT64_C(2048);

            for (int it = 0; it < searu_steps; ++it) {
                uint32_t chaos = searu_adapt_chaos_q8(stag, UINT64_C(0));
                uint64_t swan  = 0ULL;
                uint64_t jump  = searu_zipf_jump(rng_io, chaos, &swan);
                uint32_t k_mut = (swan != 0ULL) ? (3U + (searu_clz64(jump) & 1U)) : 2U;

                memcpy(cand_w, cur_w, sizeof(cand_w));
                for (uint32_t f = 0; f < k_mut; ++f) {
                    uint64_t r64  = searu_rng_next(rng_io);
                    int      t    = active_e[searu_Bounded_u32(r64, (uint32_t)n_act)];
                    int      step = 1 + (int)((r64 >> 32U) & 1ULL);
                    if ((r64 >> 33U) & 1ULL) step = -step;
                    cand_w[t] = step_angle_word(cand_w[t], step, msb0_mask[t], half, mod);
                }
                decode_words_on_tree(cand_w, edges, cos_lut, sin_lut, u_buf);
                double cand_df =
                    eval_4t_surrogate(u_buf, u0, h_pack, v_pack, h0, beta0, v_norm2, gamma);

                if (cand_df < best_df - 1e-18) {
                    best_df = cand_df;
                    cur_df  = cand_df;
                    memcpy(best_w, cand_w, sizeof(best_w));
                    memcpy(cur_w, cand_w, sizeof(cur_w));
                    changed = true;
                    stag    = 0U;
                } else {
                    uint64_t c_e = (uint64_t)llround(fmax(0.0, (cand_df - best_df) * 1e9));
                    uint64_t u_e = (uint64_t)llround(fmax(0.0, (cur_df - best_df) * 1e9));
                    bool     epi = false;
                    if (searu_crucible_accept(c_e, u_e, 0ULL, &temp_q12, UINT64_C(2048), swan,
                                              rng_io, &epi)) {
                        cur_df = cand_df;
                        memcpy(cur_w, cand_w, sizeof(cur_w));
                    }
                    ++stag;
                    temp_q12 = (temp_q12 * SEARU_OPT_INL_COOL_Q16) >> 16U;
                }
            }
            memcpy(wb_words, best_w, sizeof(best_w));
        }
    }

    decode_words_on_tree(wb_words, edges, cos_lut, sin_lut, u_out);
    pack_u_vec(u_out, p_out);
    return changed;
}

static void run_searu_4tensor_crucible(const dft_ref_t *ref, int R, int B, int C, int brot,
                                       const int64_t edges[DFT_E][2], const int msb0_mask[DFT_E],
                                       const double *alpha_qat, double *U, double *w,
                                       int64_t *words_out) {
    const int    M     = R * B;
    const int    mod   = 1 << brot;
    const int    half  = 1 << (brot - 1);
    const double scale = (double)half;

    size_t  kk  = (size_t)DFT_K * (size_t)DFT_K;
    double *Res = (double *)malloc(kk * sizeof(double));
    double *X_a = (double *)malloc((size_t)R * C * DFT_K * sizeof(double));
    double *P_a = (double *)malloc((size_t)M * DFT_K * sizeof(double));
    double *Q_b = (double *)malloc((size_t)B * DFT_K * sizeof(double));
    double *HQ  = (double *)malloc((size_t)B * DFT_K * sizeof(double));

    double *cos_lut = (double *)malloc((size_t)mod * sizeof(double));
    double *sin_lut = (double *)malloc((size_t)mod * sizeof(double));
    for (int k = 0; k < mod; ++k) {
        double th  = kPi * (double)k / scale;
        cos_lut[k] = cos(th);
        sin_lut[k] = sin(th);
    }

    for (int rb = 0; rb < M; ++rb) {
        for (int t = 0; t < DFT_E; ++t) {
            double a = alpha_qat[(size_t)rb * DFT_E + t];
            if (msb0_mask[t]) {
                int64_t k = (int64_t)llround(a * scale);
                if (k < 0) k = 0;
                if (k >= half) k = half - 1;
                words_out[(size_t)rb * DFT_E + t] = k;
            } else {
                int64_t k = (int64_t)llround(a * scale) % mod;
                if (k < 0) k += mod;
                words_out[(size_t)rb * DFT_E + t] = k;
            }
        }
        decode_words_on_tree(words_out + (size_t)rb * DFT_E, edges, cos_lut, sin_lut,
                             U + (size_t)rb * DFT_N);
        pack_u_vec(U + (size_t)rb * DFT_N, P_a + (size_t)rb * DFT_K);
    }

    memcpy(Res, ref->T_packed, kk * sizeof(double));
    for (int r = 0; r < R; ++r) {
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, C, DFT_K, B, 1.0,
                    w + (size_t)r * B * C, C, P_a + (size_t)r * B * DFT_K, DFT_K, 0.0,
                    X_a + (size_t)r * C * DFT_K, DFT_K);
    }
    cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, DFT_K, R * C, -1.0, X_a, DFT_K, 1.0, Res,
                DFT_K);

    for (int sw = 0; sw < 3; ++sw) {
        for (int r = 0; r < R; ++r) {
            double *Xr = X_a + (size_t)r * C * DFT_K;
            double *wr = w + (size_t)r * B * C;
            double *Pr = P_a + (size_t)r * B * DFT_K;
            cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, DFT_K, C, 1.0, Xr, DFT_K, 1.0, Res,
                        DFT_K);
            solve_block_w_qr_eigen(Pr, Res, B, C, 2.0e-3, wr, HQ, Q_b);
            cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, C, DFT_K, B, 1.0, wr, C, Pr, DFT_K,
                        0.0, Xr, DFT_K);
            cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, DFT_K, C, -1.0, Xr, DFT_K, 1.0,
                        Res, DFT_K);
        }
    }
    printf("  [Pass3 QAT-Snap+VarPro brot=%d] pre-SEARU err2 = %.4f mHa\n", brot,
           1e3 * frob_norm_sym_upper(Res, DFT_K));
    fflush(stdout);

    const int radii[6]       = {10, 6, 4, 3, 3, 3};
    const int searu_iters[6] = {32, 48, 64, 64, 80, 80};
    uint64_t  rng            = UINT64_C(0x5EA20DF7C5050002);
    int      *rb_ord         = (int *)malloc((size_t)M * sizeof(int));
    double   *wnorm          = (double *)malloc((size_t)M * sizeof(double));
    double    v_pack[DFT_K], q_vec[DFT_K], R4_q[DFT_K], h_pack[DFT_K];
    double    u_new[DFT_N], p_new[DFT_K], dp[DFT_K], R4_dp[DFT_K];

    for (int sw = 0; sw < 6; ++sw) {
        double t_sw = wall_sec();
        for (int rb = 0; rb < M; ++rb) {
            rb_ord[rb]       = rb;
            const double *wb = w + (size_t)rb * C;
            double        s2 = 0.0;
            for (int c = 0; c < C; ++c) s2 += wb[c] * wb[c];
            wnorm[rb] = s2;
        }
        if ((sw & 1) == 0) {
            for (int i = 0; i < M; ++i) {
                for (int j = i + 1; j < M; ++j) {
                    if (wnorm[rb_ord[j]] > wnorm[rb_ord[i]]) {
                        int t     = rb_ord[i];
                        rb_ord[i] = rb_ord[j];
                        rb_ord[j] = t;
                    }
                }
            }
        } else {
            for (int i = M - 1; i > 0; --i) {
                uint32_t j  = searu_Bounded_u32(searu_rng_next(&rng), (uint32_t)(i + 1));
                int      t  = rb_ord[i];
                rb_ord[i]   = rb_ord[j];
                rb_ord[j]   = t;
            }
        }

        int n_acc = 0;
        for (int idx = 0; idx < M; ++idx) {
            int           rb    = rb_ord[idx];
            int           r     = rb / B;
            const double *wb    = w + (size_t)rb * C;
            double        gamma = wnorm[rb];
            if (gamma < 1e-18) continue;

            double *u0 = U + (size_t)rb * DFT_N;
            double *p0 = P_a + (size_t)rb * DFT_K;
            double *Xr = X_a + (size_t)r * C * DFT_K;

            cblas_dgemv(CblasRowMajor, CblasTrans, C, DFT_K, 1.0, Xr, DFT_K, wb, 1, 0.0, v_pack, 1);
            double v_norm2 = 0.0;
            for (int k = 0; k < DFT_K; ++k) {
                v_pack[k] -= gamma * p0[k];
                v_norm2 += v_pack[k] * v_pack[k];
                q_vec[k] = v_pack[k] + gamma * p0[k];
            }

            cblas_dsymv(CblasColMajor, CblasLower, DFT_K, -1.0, Res, DFT_K, q_vec, 1, 0.0, R4_q, 1);
            for (int k = 0; k < DFT_K; ++k) h_pack[k] = 4.0 * R4_q[k];

            int64_t saved_words[DFT_E];
            memcpy(saved_words, words_out + (size_t)rb * DFT_E, sizeof(saved_words));

            bool chg = optimize_vec_searu_4tensor(
                words_out + (size_t)rb * DFT_E, u0, h_pack, v_pack, v_norm2, gamma, edges,
                msb0_mask, cos_lut, sin_lut, brot, radii[sw], searu_iters[sw], &rng, u_new, p_new);

            if (chg) {
                for (int k = 0; k < DFT_K; ++k) dp[k] = p_new[k] - p0[k];
                cblas_dsymv(CblasColMajor, CblasLower, DFT_K, -1.0, Res, DFT_K, dp, 1, 0.0, R4_dp,
                            1);
                double h0_v, b0_v;
                eval_two_packed(h_pack, v_pack, u0, &h0_v, &b0_v);
                double sur_df =
                    eval_4t_surrogate(u_new, u0, h_pack, v_pack, h0_v, b0_v, v_norm2, gamma);
                double dp_R4_dp = 0.0;
                for (int k = 0; k < DFT_K; ++k) dp_R4_dp += dp[k] * R4_dp[k];
                double exact_df = sur_df + 2.0 * gamma * dp_R4_dp;

                if (exact_df < -1e-16) {
                    cblas_dsyr2(CblasColMajor, CblasLower, DFT_K, -1.0, v_pack, 1, dp, 1, Res,
                                DFT_K);
                    cblas_dsyr(CblasColMajor, CblasLower, DFT_K, -gamma, p_new, 1, Res, DFT_K);
                    cblas_dsyr(CblasColMajor, CblasLower, DFT_K, +gamma, p0, 1, Res, DFT_K);
                    cblas_dger(CblasRowMajor, C, DFT_K, 1.0, wb, 1, dp, 1, Xr, DFT_K);
                    memcpy(u0, u_new, DFT_N * sizeof(double));
                    memcpy(p0, p_new, DFT_K * sizeof(double));
                    ++n_acc;
                } else {
                    memcpy(words_out + (size_t)rb * DFT_E, saved_words, sizeof(saved_words));
                }
            }
        }

        for (int gs = 0; gs < 2; ++gs) {
            for (int r = 0; r < R; ++r) {
                double *Xr = X_a + (size_t)r * C * DFT_K;
                double *wr = w + (size_t)r * B * C;
                double *Pr = P_a + (size_t)r * B * DFT_K;
                cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, DFT_K, C, 1.0, Xr, DFT_K, 1.0,
                            Res, DFT_K);
                solve_block_w_qr_eigen(Pr, Res, B, C, 1.8e-3, wr, HQ, Q_b);
                cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, C, DFT_K, B, 1.0, wr, C, Pr,
                            DFT_K, 0.0, Xr, DFT_K);
                cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, DFT_K, C, -1.0, Xr, DFT_K, 1.0,
                            Res, DFT_K);
            }
        }

        double e2 = frob_norm_sym_upper(Res, DFT_K);
        printf("  [Pass3 SEARU-4T brot=%d sw=%d/6] err2 = %.4f mHa (accepted=%d/%d) [%.2fs]\n",
               brot, sw + 1, 1e3 * e2, n_acc, M, wall_sec() - t_sw);
        fflush(stdout);
        if (e2 <= 0.0193 && sw >= 1) break;
    }

    free(wnorm);
    free(rb_ord);
    free(sin_lut);
    free(cos_lut);
    free(HQ);
    free(Q_b);
    free(P_a);
    free(X_a);
    free(Res);
}

/* CLI: ./dft_bootstrap <dft_ref.bin> <out_sol.bin> <R> <B> <C> <brot> [cont_in.bin] [cont_out.bin] */
int main(int argc, char **argv) {
    if (argc < 7) {
        fprintf(stderr,
                "Usage: %s <dft_ref.bin> <out_sol.bin> <R> <B> <C> <brot> [cont_in.bin] "
                "[cont_out.bin]\n",
                argv[0]);
        return 1;
    }
    const char *ref_path = argv[1];
    const char *out_path = argv[2];
    int         R        = atoi(argv[3]);
    int         B        = atoi(argv[4]);
    int         C        = atoi(argv[5]);
    int         brot     = atoi(argv[6]);
    const char *cont_in  =
        (argc >= 8 && strlen(argv[7]) > 0 && strcmp(argv[7], "-") != 0) ? argv[7] : NULL;
    const char *cont_out =
        (argc >= 9 && strlen(argv[8]) > 0 && strcmp(argv[8], "-") != 0) ? argv[8] : NULL;

    double t_start = wall_sec();
    printf("=== SEARU x FLOW v%s DF-THC Bootstrap v2: target (%d, %d, %d) @ brot=%d ===\n",
           FLOW_VERSION_STRING, R, B, C, brot);
    fflush(stdout);

    dft_ref_t ref;
    load_dft_ref(ref_path, &ref);

    double *U = (double *)malloc((size_t)R * B * DFT_N * sizeof(double));
    double *w = (double *)malloc((size_t)R * B * C * sizeof(double));
    if (!U || !w) abort();

    int     cur_R = ref.R0, cur_B = ref.B0, cur_C = ref.C0;
    double *cur_U = (double *)malloc((size_t)cur_R * cur_B * DFT_N * sizeof(double));
    double *cur_w = (double *)malloc((size_t)cur_R * cur_B * cur_C * sizeof(double));
    if (!cur_U || !cur_w) abort();
    memcpy(cur_U, ref.U0, (size_t)cur_R * cur_B * DFT_N * sizeof(double));
    memcpy(cur_w, ref.w0, (size_t)cur_R * cur_B * cur_C * sizeof(double));

    if (cont_in != NULL) {
        FILE *fc = fopen(cont_in, "rb");
        if (!fc) {
            fprintf(stderr, "FATAL: cannot open cont_in %s\n", cont_in);
            return 1;
        }
        int hdr_c[3];
        if (fread(hdr_c, sizeof(int), 3, fc) != 3) return 1;
        free(cur_U);
        free(cur_w);
        cur_R = hdr_c[0];
        cur_B = hdr_c[1];
        cur_C = hdr_c[2];
        cur_U = (double *)malloc((size_t)cur_R * cur_B * DFT_N * sizeof(double));
        cur_w = (double *)malloc((size_t)cur_R * cur_B * cur_C * sizeof(double));
        if (!cur_U || !cur_w) abort();
        if (fread(cur_U, sizeof(double), (size_t)cur_R * cur_B * DFT_N, fc) !=
                (size_t)cur_R * cur_B * DFT_N ||
            fread(cur_w, sizeof(double), (size_t)cur_R * cur_B * cur_C, fc) !=
                (size_t)cur_R * cur_B * cur_C) {
            return 1;
        }
        fclose(fc);
        printf("  Loaded continuous checkpoint %s (%d,%d,%d) err2=%.4f mHa\n", cont_in, cur_R,
               cur_B, cur_C, 1e3 * exact_err2_from_uw(&ref, cur_U, cur_w, cur_R, cur_B, cur_C));
    }

    /* Progressive Homotopy Compression Ladder to (R, B, C) */
    const int ladder[5][3] = {
        {50, 48, 25},
        {38, 38, 22},
        {30, 34, 20},
        {25, 26, 20},
        {R, B, C},
    };
    for (int st = 0; st < 5; ++st) {
        int nxt_R = ladder[st][0];
        int nxt_B = ladder[st][1];
        int nxt_C = ladder[st][2];
        if (st < 4) {
            if (nxt_R <= R || nxt_B < B || nxt_C < C) continue;
            if (nxt_R >= cur_R && nxt_B >= cur_B && nxt_C >= cur_C) continue;
            if (nxt_B < B) nxt_B = B;
            if (nxt_C < C) nxt_C = C;
        } else {
            nxt_R = R;
            nxt_B = B;
            nxt_C = C;
        }

        if (nxt_R != cur_R || nxt_B != cur_B || nxt_C != cur_C) {
            double *nxt_U = (double *)malloc((size_t)nxt_R * nxt_B * DFT_N * sizeof(double));
            double *nxt_w = (double *)malloc((size_t)nxt_R * nxt_B * nxt_C * sizeof(double));
            compress_blocks_4tensor(cur_U, cur_w, cur_R, cur_B, cur_C, nxt_R, nxt_B, nxt_C, nxt_U,
                                    nxt_w);
            free(cur_U);
            free(cur_w);
            cur_R = nxt_R;
            cur_B = nxt_B;
            cur_C = nxt_C;
            cur_U = nxt_U;
            cur_w = nxt_w;
        }

        int n_it = (st == 4) ? 3600 : 1200;
        run_continuous_lbfgs_stage(&ref, cur_R, cur_B, cur_C, n_it, 2e-9, cur_U, cur_w);
        polish_w_exact_varpro(&ref, cur_R, cur_B, cur_C, 3, 1.2e-3, cur_U, cur_w);
        if (st == 4) {
            run_continuous_lbfgs_stage(&ref, cur_R, cur_B, cur_C, 2400, 1e-9, cur_U, cur_w);
            polish_w_exact_varpro(&ref, cur_R, cur_B, cur_C, 3, 1.2e-3, cur_U, cur_w);
        }
        double e2 = exact_err2_from_uw(&ref, cur_U, cur_w, cur_R, cur_B, cur_C);
        printf("  [Pass1 Stage (%d,%d,%d)] EXACT cont_err2 = %.4f mHa [%.2fs]\n", cur_R, cur_B,
               cur_C, 1e3 * e2, wall_sec() - t_start);
        fflush(stdout);
    }

    memcpy(U, cur_U, (size_t)R * B * DFT_N * sizeof(double));
    memcpy(w, cur_w, (size_t)R * B * C * sizeof(double));
    free(cur_U);
    free(cur_w);

    if (cont_out != NULL) {
        FILE *fco = fopen(cont_out, "wb");
        if (fco) {
            int hdr_c[3] = {R, B, C};
            fwrite(hdr_c, sizeof(int), 3, fco);
            fwrite(U, sizeof(double), (size_t)R * B * DFT_N, fco);
            fwrite(w, sizeof(double), (size_t)R * B * C, fco);
            fclose(fco);
            printf("  Saved continuous checkpoint %s\n", cont_out);
        }
    }

    /* If brot <= 0, user requested continuous compression checkpoint only! */
    if (brot <= 0) {
        printf("=== Continuous checkpoint complete in %.2fs ===\n", wall_sec() - t_start);
        return 0;
    }

    /* Pass 2: SEARU Multi-Chain MSB-0 Givens Tree */
    int64_t edges[DFT_E][2];
    int     msb0_mask[DFT_E];
    synthesize_givens_tree_searu(U, w, R, B, C, brot, edges);
    compute_msb0_mask(edges, msb0_mask);

    /* Pass 2.5: Multi-Threaded Angle-Space QAT Homotopy L-BFGS */
    double *alpha_qat = (double *)malloc((size_t)R * B * DFT_E * sizeof(double));
    run_qat_angle_lbfgs_c17(&ref, R, B, C, brot, edges, msb0_mask, U, w, alpha_qat);

    /* Pass 3: Live-Residual 4-Tensor SEARU Sparse-k Lattice Crucible */
    int64_t *words_2b = (int64_t *)malloc((size_t)R * B * DFT_E * sizeof(int64_t));
    run_searu_4tensor_crucible(&ref, R, B, C, brot, edges, msb0_mask, alpha_qat, U, w, words_2b);

    FILE *fo = fopen(out_path, "wb");
    if (!fo) {
        fprintf(stderr, "FATAL: cannot write %s\n", out_path);
        return 1;
    }
    int32_t meta[4] = {R, B, C, brot};
    fwrite(meta, sizeof(int32_t), 4, fo);
    fwrite(edges, sizeof(int64_t), DFT_E * 2, fo);
    fwrite(words_2b, sizeof(int64_t), (size_t)R * B * DFT_E, fo);
    fwrite(U, sizeof(double), (size_t)R * B * DFT_N, fo);
    fwrite(w, sizeof(double), (size_t)R * B * C, fo);
    fclose(fo);

    free(words_2b);
    free(alpha_qat);
    free(w);
    free(U);
    printf("=== Completed %s in %.2fs ===\n", out_path, wall_sec() - t_start);
    return 0;
}
