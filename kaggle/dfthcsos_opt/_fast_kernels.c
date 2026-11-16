#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

void pack_outer_f64(const double* __restrict__ Un, double* __restrict__ P, int M, int N) {
    const double sqrt2 = 1.41421356237309504880;
    const int K = (N * (N + 1)) / 2;
    for (int i = 0; i < M; ++i) {
        const double* u = Un + (size_t)i * N;
        double* p = P + (size_t)i * K;
        int k = 0;
        for (int r = 0; r < N; ++r) {
            const double ur_s2 = u[r] * sqrt2;
            for (int c = 0; c < r; ++c) {
                p[k++] = ur_s2 * u[c];
            }
            p[k++] = u[r] * u[r];
        }
    }
}

void grad_outer_f64(
    const double* __restrict__ grad_P,
    const double* __restrict__ Un,
    double* __restrict__ grad_Un,
    int M,
    int N
) {
    const double sqrt2 = 1.41421356237309504880;
    const int K = (N * (N + 1)) / 2;
    for (int i = 0; i < M; ++i) {
        const double* gp = grad_P + (size_t)i * K;
        const double* u = Un + (size_t)i * N;
        double* gu = grad_Un + (size_t)i * N;
        for (int r = 0; r < N; ++r) {
            gu[r] = 0.0;
        }
        int k = 0;
        for (int r = 0; r < N; ++r) {
            const double ur_s2 = u[r] * sqrt2;
            double acc_r = 0.0;
            for (int c = 0; c < r; ++c) {
                const double g = gp[k++];
                acc_r += g * u[c];
                gu[c] += g * ur_s2;
            }
            const double g_diag = gp[k++];
            gu[r] += acc_r * sqrt2 + 2.0 * g_diag * u[r];
        }
    }
}

static inline void decode_one_lut(
    const int64_t* words,
    const int64_t* edges,
    int N,
    const double* cos_lut,
    const double* sin_lut,
    double* u
) {
    u[0] = 1.0;
    for (int j = 1; j < N; ++j) u[j] = 0.0;
    for (int t = 0; t < N - 1; ++t) {
        int64_t p = edges[2 * t];
        int64_t ch = edges[2 * t + 1];
        int64_t wt = words[t];
        double cs = cos_lut[wt];
        double sn = sin_lut[wt];
        double vp = u[p];
        // Note: child ch has not been visited yet on a tree! Wait: check parent-before-child tree invariant:
        // In reconstruct_u_vectors, child is unvisited so u[ch] == 0 before edge t!
        u[p] = cs * vp;
        u[ch] = sn * vp;
    }
}

static void build_trig_lut(int brot, double* cos_lut, double* sin_lut) {
    const int mod = 1 << brot;
    const double pi = 3.14159265358979323846;
    const double inv_scale_pi = pi / (double)(1 << (brot - 1));
    for (int w = 0; w < mod; ++w) {
        double th = (double)w * inv_scale_pi;
        cos_lut[w] = cos(th);
        sin_lut[w] = sin(th);
    }
}

void cd_quantize_overlap(
    int64_t* __restrict__ words,
    const double* __restrict__ u_target,
    const int64_t* __restrict__ edges,
    const int64_t* __restrict__ mod_per_edge,
    int M,
    int N,
    int brot,
    int n_sweeps
) {
    const int E = N - 1;
    const int mod = 1 << brot;
    double* cos_lut = (double*)malloc((size_t)mod * sizeof(double));
    double* sin_lut = (double*)malloc((size_t)mod * sizeof(double));
    build_trig_lut(brot, cos_lut, sin_lut);
    double u_buf[128];

    for (int i = 0; i < M; ++i) {
        int64_t* w = words + (size_t)i * E;
        const double* ut = u_target + (size_t)i * N;
        decode_one_lut(w, edges, N, cos_lut, sin_lut, u_buf);
        double best_ov = 0.0;
        for (int j = 0; j < N; ++j) best_ov += u_buf[j] * ut[j];
        best_ov = fabs(best_ov);

        for (int sw = 0; sw < n_sweeps; ++sw) {
            int improved = 0;
            for (int t = 0; t < E; ++t) {
                int64_t orig = w[t];
                int64_t mod_t = mod_per_edge[t];
                int64_t best_w = orig;
                for (int delta = -3; delta <= 3; ++delta) {
                    if (delta == 0) continue;
                    int64_t cand = (orig + delta) % mod_t;
                    if (cand < 0) cand += mod_t;
                    w[t] = cand;
                    decode_one_lut(w, edges, N, cos_lut, sin_lut, u_buf);
                    double ov = 0.0;
                    for (int j = 0; j < N; ++j) ov += u_buf[j] * ut[j];
                    ov = fabs(ov);
                    if (ov > best_ov + 1e-15) {
                        best_ov = ov;
                        best_w = cand;
                        improved = 1;
                    }
                }
                w[t] = best_w;
            }
            if (!improved) break;
        }
    }
    free(cos_lut);
    free(sin_lut);
}

// Evaluate u^T M u where M is stored in packed lower-triangular form `m_pack` (len K = N*(N+1)/2)
// with sqrt(2) off-diagonal weights: <m_pack, pack_outer(u)> == sum_r m_rr u_r^2 + sum_{r>c} m_rc * sqrt(2) * u_r * u_c!
static inline double eval_quad_packed(const double* __restrict__ m_pack, const double* __restrict__ u, int N) {
    const double sqrt2 = 1.41421356237309504880;
    double total = 0.0;
    int k = 0;
    for (int r = 0; r < N; ++r) {
        double acc = 0.0;
        for (int c = 0; c < r; ++c) {
            acc += m_pack[k++] * u[c];
        }
        total += acc * (u[r] * sqrt2) + m_pack[k++] * (u[r] * u[r]);
    }
    return total;
}

static inline void pack_one(const double* __restrict__ u, double* __restrict__ p, int N) {
    const double sqrt2 = 1.41421356237309504880;
    int k = 0;
    for (int r = 0; r < N; ++r) {
        double ur_s2 = u[r] * sqrt2;
        for (int c = 0; c < r; ++c) {
            p[k++] = ur_s2 * u[c];
        }
        p[k++] = u[r] * u[r];
    }
}

// Residual-cancelling Block Coordinate Descent for a tensor block X_target (C, K):
// Approximates X_target[c, k] ~= sum_{b=0}^{B-1} w[b, c] * P[b, k]
// If `update_w == 1`, also updates w[b, :] = E_{-b} @ P[b] (since ||P[b]||_2 == 1!) at each step!
// Note: for 1-body h1_prime, C = 1, w[b, 0] = lambda_b, and `update_w == 0` (since eigenvalues lambda_b are fixed)!
void bcd_quantize_block_packed(
    int64_t* __restrict__ words,          // (B, E)
    double* __restrict__ w_mat,           // (B, C)
    const double* __restrict__ X_target,  // (C, K)
    const int64_t* __restrict__ edges,    // (E, 2)
    const int64_t* __restrict__ mod_per_edge, // (E,)
    const int64_t* __restrict__ b_order,  // (B,) order in which to visit basis vectors
    int B,
    int C,
    int N,
    int brot,
    int n_sweeps,
    int search_radius,
    int update_w
) {
    const int E = N - 1;
    const int K = (N * (N + 1)) / 2;
    const int mod = 1 << brot;
    double* cos_lut = (double*)malloc((size_t)mod * sizeof(double));
    double* sin_lut = (double*)malloc((size_t)mod * sizeof(double));
    build_trig_lut(brot, cos_lut, sin_lut);

    // Precompute current P (B, K) and residual R = X_target - w^T P (C, K)
    double* P = (double*)malloc((size_t)B * K * sizeof(double));
    double* R = (double*)malloc((size_t)C * K * sizeof(double));
    double* m_pack = (double*)malloc((size_t)K * sizeof(double));
    double u_buf[128];

    memcpy(R, X_target, (size_t)C * K * sizeof(double));
    for (int b = 0; b < B; ++b) {
        decode_one_lut(words + (size_t)b * E, edges, N, cos_lut, sin_lut, u_buf);
        double* pb = P + (size_t)b * K;
        pack_one(u_buf, pb, N);
        const double* wb = w_mat + (size_t)b * C;
        for (int c = 0; c < C; ++c) {
            double wbc = wb[c];
            double* rc = R + (size_t)c * K;
            for (int k = 0; k < K; ++k) {
                rc[k] -= wbc * pb[k];
            }
        }
    }

    for (int sw = 0; sw < n_sweeps; ++sw) {
        int any_change = 0;
        for (int idx = 0; idx < B; ++idx) {
            int b = (int)b_order[idx];
            int64_t* wb_words = words + (size_t)b * E;
            double* wb = w_mat + (size_t)b * C;
            double* pb = P + (size_t)b * K;

            // Add back b's contribution so R becomes E_{-b} = X_target - sum_{b' != b} w_{b'} P_{b'}
            double w_norm2 = 0.0;
            for (int c = 0; c < C; ++c) {
                double wbc = wb[c];
                w_norm2 += wbc * wbc;
                double* rc = R + (size_t)c * K;
                for (int k = 0; k < K; ++k) {
                    rc[k] += wbc * pb[k];
                }
            }

            if (w_norm2 > 1e-24) {
                // Form m_pack[k] = sum_c wb[c] * E_{-b}[c, k]
                for (int k = 0; k < K; ++k) m_pack[k] = 0.0;
                for (int c = 0; c < C; ++c) {
                    double wbc = wb[c];
                    const double* rc = R + (size_t)c * K;
                    for (int k = 0; k < K; ++k) {
                        m_pack[k] += wbc * rc[k];
                    }
                }

                // Maximize <m_pack, pack(u)> over wb_words!
                double best_score = 0.0;
                for (int k = 0; k < K; ++k) best_score += m_pack[k] * pb[k];

                for (int inner = 0; inner < 2; ++inner) {
                    int inner_imp = 0;
                    for (int t = 0; t < E; ++t) {
                        int64_t orig = wb_words[t];
                        int64_t mod_t = mod_per_edge[t];
                        int64_t best_w = orig;
                        for (int delta = -search_radius; delta <= search_radius; ++delta) {
                            if (delta == 0) continue;
                            int64_t cand = (orig + delta) % mod_t;
                            if (cand < 0) cand += mod_t;
                            wb_words[t] = cand;
                            decode_one_lut(wb_words, edges, N, cos_lut, sin_lut, u_buf);
                            double sc = eval_quad_packed(m_pack, u_buf, N);
                            if (sc > best_score + 1e-15) {
                                best_score = sc;
                                best_w = cand;
                                inner_imp = 1;
                                any_change = 1;
                            }
                        }
                        wb_words[t] = best_w;
                    }
                    if (!inner_imp) break;
                }

                decode_one_lut(wb_words, edges, N, cos_lut, sin_lut, u_buf);
                pack_one(u_buf, pb, N);

                if (update_w) {
                    // Exact optimal wb[c] = <E_{-b}[c, :], pb> since ||pb||_2 == 1!
                    for (int c = 0; c < C; ++c) {
                        const double* rc = R + (size_t)c * K;
                        double dot = 0.0;
                        for (int k = 0; k < K; ++k) dot += rc[k] * pb[k];
                        wb[c] = dot;
                    }
                }
            }

            // Subtract b's updated contribution from R
            for (int c = 0; c < C; ++c) {
                double wbc = wb[c];
                double* rc = R + (size_t)c * K;
                for (int k = 0; k < K; ++k) {
                    rc[k] -= wbc * pb[k];
                }
            }
        }
        if (!any_change) break;
    }

    free(cos_lut);
    free(sin_lut);
    free(P);
    free(R);
    free(m_pack);
}
