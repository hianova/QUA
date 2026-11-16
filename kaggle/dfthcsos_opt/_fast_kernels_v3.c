#include <dispatch/dispatch.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static const double kPi = 3.14159265358979323846264338327950288;
static const double kSqrt2 = 1.41421356237309504880168872420969808;

static inline void pack_one_u(const double* __restrict__ u, double* __restrict__ p, int N) {
    int k = 0;
    for (int r = 0; r < N; ++r) {
        const double ur_s2 = u[r] * kSqrt2;
        for (int c = 0; c < r; ++c) {
            p[k++] = ur_s2 * u[c];
        }
        p[k++] = u[r] * u[r];
    }
}

static inline void grad_one_u(
    const double* __restrict__ gp,
    const double* __restrict__ u,
    double* __restrict__ gu,
    int N
) {
    for (int r = 0; r < N; ++r) gu[r] = 0.0;
    int k = 0;
    for (int r = 0; r < N; ++r) {
        const double ur_s2 = u[r] * kSqrt2;
        double acc_r = 0.0;
        for (int c = 0; c < r; ++c) {
            const double g = gp[k++];
            acc_r += g * u[c];
            gu[c] += g * ur_s2;
        }
        const double g_diag = gp[k++];
        gu[r] += acc_r * kSqrt2 + 2.0 * g_diag * u[r];
    }
}

// Multi-threaded forward Givens angle decoding + outer-product packing + periodic quantization-well loss!
// alpha: (M, E) in units of theta / pi
// edges: (E, 2) int64
// w_norm2: (M,) float64 weight per vector for grid penalty
// U_out: (M, N), P_out: (M, K), vp_out: (M, E), grad_grid_out: (M, E)
// Returns total grid penalty sum_m w_norm2[m] * gamma_grid * sum_t sin^2(pi * scale * alpha[m, t]).
double forward_angles_pack_mt(
    const double* __restrict__ alpha,
    const int64_t* __restrict__ edges,
    const double* __restrict__ w_norm2,
    int M,
    int N,
    int brot,
    double gamma_grid,
    double* __restrict__ U_out,
    double* __restrict__ P_out,
    double* __restrict__ vp_out,
    double* __restrict__ grad_grid_out
) {
    const int E = N - 1;
    const int K = (N * (N + 1)) / 2;
    const double scale = (double)(1 << (brot - 1));
    const double pi_scale = kPi * scale;

    const int chunk = 64;
    const int n_chunks = (M + chunk - 1) / chunk;
    double* chunk_loss = (double*)calloc((size_t)n_chunks, sizeof(double));

    dispatch_apply((size_t)n_chunks, dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0), ^(size_t cid) {
        int i_lo = (int)cid * chunk;
        int i_hi = i_lo + chunk;
        if (i_hi > M) i_hi = M;
        double loc_pen = 0.0;

        for (int i = i_lo; i < i_hi; ++i) {
            const double* a = alpha + (size_t)i * E;
            double* u = U_out + (size_t)i * N;
            double* vp = vp_out + (size_t)i * E;
            double* gg = grad_grid_out + (size_t)i * E;

            u[0] = 1.0;
            for (int j = 1; j < N; ++j) u[j] = 0.0;

            const double wg = gamma_grid * w_norm2[i];
            for (int t = 0; t < E; ++t) {
                int64_t p = edges[2 * t];
                int64_t ch = edges[2 * t + 1];
                double th = kPi * a[t];
                double cs = cos(th);
                double sn = sin(th);
                double v_prev = u[p];
                vp[t] = v_prev;
                u[p] = cs * v_prev;
                u[ch] = sn * v_prev;

                if (wg > 0.0) {
                    double z = pi_scale * a[t];
                    double sz = sin(z);
                    loc_pen += wg * (sz * sz);
                    gg[t] = wg * pi_scale * sin(2.0 * z);
                } else {
                    gg[t] = 0.0;
                }
            }
            pack_one_u(u, P_out + (size_t)i * K, N);
        }
        chunk_loss[cid] = loc_pen;
    });

    double total_pen = 0.0;
    for (int c = 0; c < n_chunks; ++c) total_pen += chunk_loss[c];
    free(chunk_loss);
    return total_pen;
}

// Multi-threaded exact reverse-mode tree adjoint: grad_P -> grad_U -> grad_alpha + grad_grid!
void backward_angles_mt(
    const double* __restrict__ grad_P,
    const double* __restrict__ U_in,
    const double* __restrict__ vp_in,
    const double* __restrict__ alpha,
    const int64_t* __restrict__ edges,
    const double* __restrict__ grad_grid,
    int M,
    int N,
    double* __restrict__ grad_alpha_out
) {
    const int E = N - 1;
    const int K = (N * (N + 1)) / 2;
    const int chunk = 64;
    const int n_chunks = (M + chunk - 1) / chunk;

    dispatch_apply((size_t)n_chunks, dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0), ^(size_t cid) {
        int i_lo = (int)cid * chunk;
        int i_hi = i_lo + chunk;
        if (i_hi > M) i_hi = M;
        double gu[128];

        for (int i = i_lo; i < i_hi; ++i) {
            grad_one_u(grad_P + (size_t)i * K, U_in + (size_t)i * N, gu, N);
            const double* a = alpha + (size_t)i * E;
            const double* vp = vp_in + (size_t)i * E;
            const double* gg = grad_grid + (size_t)i * E;
            double* ga = grad_alpha_out + (size_t)i * E;

            for (int t = E - 1; t >= 0; --t) {
                int64_t p = edges[2 * t];
                int64_t ch = edges[2 * t + 1];
                double th = kPi * a[t];
                double cs = cos(th);
                double sn = sin(th);
                double gp = gu[p];
                double gc = gu[ch];
                ga[t] = kPi * vp[t] * (-sn * gp + cs * gc) + gg[t];
                gu[p] = cs * gp + sn * gc;
                gu[ch] = 0.0;
            }
        }
    });
}

void build_trig_lut_v3(int brot, double* cos_lut, double* sin_lut) {
    const int mod = 1 << brot;
    const double inv_scale_pi = kPi / (double)(1 << (brot - 1));
    for (int w = 0; w < mod; ++w) {
        double th = (double)w * inv_scale_pi;
        cos_lut[w] = cos(th);
        sin_lut[w] = sin(th);
    }
}

static inline void eval_two_quad_packed_v3(
    const double* __restrict__ h_pack,
    const double* __restrict__ v_pack,
    const double* __restrict__ u,
    int N,
    double* out_h,
    double* out_v
) {
    double tot_h = 0.0;
    double tot_v = 0.0;
    int k = 0;
    for (int r = 0; r < N; ++r) {
        double acc_h = 0.0;
        double acc_v = 0.0;
        for (int c = 0; c < r; ++c) {
            double uc = u[c];
            acc_h += h_pack[k] * uc;
            acc_v += v_pack[k] * uc;
            ++k;
        }
        double ur = u[r];
        double ur_s2 = ur * kSqrt2;
        double ur2 = ur * ur;
        tot_h += acc_h * ur_s2 + h_pack[k] * ur2;
        tot_v += acc_v * ur_s2 + v_pack[k] * ur2;
        ++k;
    }
    *out_h = tot_h;
    *out_v = tot_v;
}

static inline double eval_4tensor_df(
    const double* __restrict__ u_cand,
    const double* __restrict__ u0,
    const double* __restrict__ h_pack,
    const double* __restrict__ v_pack,
    double h0,
    double beta0,
    double v_norm2,
    double gamma,
    int N
) {
    double dot0 = 0.0;
    for (int j = 0; j < N; ++j) dot0 += u0[j] * u_cand[j];
    double alpha = dot0 * dot0;
    double one_m_a = 1.0 - alpha;

    double hp, beta;
    eval_two_quad_packed_v3(h_pack, v_pack, u_cand, N, &hp, &beta);
    double db = beta - beta0;

    return (hp - h0)
        + 4.0 * v_norm2 * one_m_a
        + 2.0 * db * db
        + 2.0 * gamma * gamma * (1.0 - alpha * alpha)
        + 4.0 * gamma * (beta + beta0) * one_m_a;
}

// Prefix-cached wide-radius + coupled parent-child 2-edge 4-tensor coordinate descent on ONE vector!
int cd_4tensor_vec_v3(
    int64_t* __restrict__ wb_words,
    const double* __restrict__ u0,
    const double* __restrict__ h_pack,
    const double* __restrict__ v_pack,
    double v_norm2,
    double gamma,
    const int64_t* __restrict__ edges,
    const int64_t* __restrict__ mod_per_edge,
    const double* __restrict__ cos_lut,
    const double* __restrict__ sin_lut,
    int N,
    int search_radius,
    int n_inner,
    int do_coupled_pairs,
    double* __restrict__ u_out,
    double* __restrict__ p_out
) {
    const int E = N - 1;
    double h0, beta0;
    eval_two_quad_packed_v3(h_pack, v_pack, u0, N, &h0, &beta0);

    double best_df = 0.0;
    int changed = 0;
    double u_pref[128];
    double u_buf[128];

    for (int inner = 0; inner < n_inner; ++inner) {
        int imp = 0;
        u_pref[0] = 1.0;
        for (int j = 1; j < N; ++j) u_pref[j] = 0.0;

        for (int t = 0; t < E; ++t) {
            int64_t orig = wb_words[t];
            int64_t mod_t = mod_per_edge[t];
            int64_t best_w = orig;

            int rad = search_radius;
            if (2 * rad + 1 >= (int)mod_t) {
                // Full discrete grid search on edge t!
                for (int64_t cand = 0; cand < mod_t; ++cand) {
                    if (cand == orig) continue;
                    memcpy(u_buf, u_pref, (size_t)N * sizeof(double));
                    int64_t p0 = edges[2 * t], c0 = edges[2 * t + 1];
                    double vp0 = u_buf[p0];
                    u_buf[p0] = cos_lut[cand] * vp0;
                    u_buf[c0] = sin_lut[cand] * vp0;
                    for (int s = t + 1; s < E; ++s) {
                        int64_t ps = edges[2 * s], cs = edges[2 * s + 1];
                        int64_t ws = wb_words[s];
                        double vps = u_buf[ps];
                        u_buf[ps] = cos_lut[ws] * vps;
                        u_buf[cs] = sin_lut[ws] * vps;
                    }
                    double df = eval_4tensor_df(u_buf, u0, h_pack, v_pack, h0, beta0, v_norm2, gamma, N);
                    if (df < best_df - 1e-18) {
                        best_df = df;
                        best_w = cand;
                        imp = 1;
                        changed = 1;
                    }
                }
            } else {
                for (int delta = -rad; delta <= rad; ++delta) {
                    if (delta == 0) continue;
                    int64_t cand = (orig + delta) % mod_t;
                    if (cand < 0) cand += mod_t;
                    memcpy(u_buf, u_pref, (size_t)N * sizeof(double));
                    int64_t p0 = edges[2 * t], c0 = edges[2 * t + 1];
                    double vp0 = u_buf[p0];
                    u_buf[p0] = cos_lut[cand] * vp0;
                    u_buf[c0] = sin_lut[cand] * vp0;
                    for (int s = t + 1; s < E; ++s) {
                        int64_t ps = edges[2 * s], cs = edges[2 * s + 1];
                        int64_t ws = wb_words[s];
                        double vps = u_buf[ps];
                        u_buf[ps] = cos_lut[ws] * vps;
                        u_buf[cs] = sin_lut[ws] * vps;
                    }
                    double df = eval_4tensor_df(u_buf, u0, h_pack, v_pack, h0, beta0, v_norm2, gamma, N);
                    if (df < best_df - 1e-18) {
                        best_df = df;
                        best_w = cand;
                        imp = 1;
                        changed = 1;
                    }
                }
            }
            wb_words[t] = best_w;

            // Coupled 2-edge search with BOTH first child (`edges[2*s] == ch_t`) and next sibling (`edges[2*s] == p_t`)!
            if (do_coupled_pairs) {
                int64_t p_t = edges[2 * t];
                int64_t ch_t = edges[2 * t + 1];
                int partners[2] = {-1, -1};
                for (int s = t + 1; s < E; ++s) {
                    if (partners[0] < 0 && edges[2 * s] == ch_t) partners[0] = s;
                    if (partners[1] < 0 && edges[2 * s] == p_t) partners[1] = s;
                    if (partners[0] >= 0 && partners[1] >= 0) break;
                }
                for (int pidx = 0; pidx < 2; ++pidx) {
                    int t2 = partners[pidx];
                    if (t2 < 0) continue;
                    int64_t w1_cur = wb_words[t];
                    int64_t w2_cur = wb_words[t2];
                    int64_t m2 = mod_per_edge[t2];
                    int64_t bw1 = w1_cur, bw2 = w2_cur;
                    for (int d1 = -2; d1 <= 2; ++d1) {
                        if (d1 == 0) continue;
                        int64_t c1 = (w1_cur + d1) % mod_t;
                        if (c1 < 0) c1 += mod_t;
                        for (int d2 = -2; d2 <= 2; ++d2) {
                            if (d2 == 0) continue;
                            int64_t c2 = (w2_cur + d2) % m2;
                            if (c2 < 0) c2 += m2;
                            wb_words[t] = c1;
                            wb_words[t2] = c2;
                            memcpy(u_buf, u_pref, (size_t)N * sizeof(double));
                            for (int s = t; s < E; ++s) {
                                int64_t ps = edges[2 * s], cs = edges[2 * s + 1];
                                int64_t ws = wb_words[s];
                                double vps = u_buf[ps];
                                u_buf[ps] = cos_lut[ws] * vps;
                                u_buf[cs] = sin_lut[ws] * vps;
                            }
                            double df = eval_4tensor_df(u_buf, u0, h_pack, v_pack, h0, beta0, v_norm2, gamma, N);
                            if (df < best_df - 1e-18) {
                                best_df = df;
                                bw1 = c1;
                                bw2 = c2;
                                imp = 1;
                                changed = 1;
                            }
                        }
                    }
                    wb_words[t] = bw1;
                    wb_words[t2] = bw2;
                }
            }

            // Advance prefix `u_pref` across edge t using updated `wb_words[t]`!
            int64_t pt = edges[2 * t], ct = edges[2 * t + 1];
            int64_t wt = wb_words[t];
            double vpt = u_pref[pt];
            u_pref[pt] = cos_lut[wt] * vpt;
            u_pref[ct] = sin_lut[wt] * vpt;
        }
        if (!imp) break;
    }

    u_out[0] = 1.0;
    for (int j = 1; j < N; ++j) u_out[j] = 0.0;
    for (int t = 0; t < E; ++t) {
        int64_t pt = edges[2 * t], ct = edges[2 * t + 1];
        int64_t wt = wb_words[t];
        double vpt = u_out[pt];
        u_out[pt] = cos_lut[wt] * vpt;
        u_out[ct] = sin_lut[wt] * vpt;
    }
    pack_one_u(u_out, p_out, N);
    return changed;
}
