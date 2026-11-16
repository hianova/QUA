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
        for (int r = 0; r < N; ++r) gu[r] = 0.0;
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
        double vp = u[p];
        u[p] = cos_lut[wt] * vp;
        u[ch] = sin_lut[wt] * vp;
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

// Fused evaluation of BOTH <h_pack, pack(u)> AND <v_pack, pack(u)> in a single L1-cache pass!
static inline void eval_two_quad_packed(
    const double* __restrict__ h_pack,
    const double* __restrict__ v_pack,
    const double* __restrict__ u,
    int N,
    double* out_h,
    double* out_v
) {
    const double sqrt2 = 1.41421356237309504880;
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
        double ur_s2 = ur * sqrt2;
        double ur2 = ur * ur;
        tot_h += acc_h * ur_s2 + h_pack[k] * ur2;
        tot_v += acc_v * ur_s2 + v_pack[k] * ur2;
        ++k;
    }
    *out_h = tot_h;
    *out_v = tot_v;
}

// Coordinate descent on a SINGLE vector `wb_words` (len N-1) minimizing the EXACT 4-tensor delta ||M(p) - T||_F^2!
// Returns 1 if `wb_words` changed, 0 otherwise, and writes new decoded `u_out` (len N) and `p_out` (len K).
int cd_4tensor_single_vec(
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
    double* __restrict__ u_out,
    double* __restrict__ p_out
) {
    const int E = N - 1;
    double h0, beta0;
    eval_two_quad_packed(h_pack, v_pack, u0, N, &h0, &beta0);

    double best_df = 0.0;
    int changed = 0;
    double u_buf[128];

    for (int inner = 0; inner < n_inner; ++inner) {
        int imp = 0;
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

                double dot0 = 0.0;
                for (int j = 0; j < N; ++j) dot0 += u0[j] * u_buf[j];
                double alpha = dot0 * dot0;
                double one_m_a = 1.0 - alpha;

                double hp, beta;
                eval_two_quad_packed(h_pack, v_pack, u_buf, N, &hp, &beta);
                double db = beta - beta0;

                double df = (hp - h0)
                    + 4.0 * v_norm2 * one_m_a
                    + 2.0 * db * db
                    + 2.0 * gamma * gamma * (1.0 - alpha * alpha)
                    + 4.0 * gamma * (beta + beta0) * one_m_a;

                if (df < best_df - 1e-18) {
                    best_df = df;
                    best_w = cand;
                    imp = 1;
                    changed = 1;
                }
            }
            wb_words[t] = best_w;
        }
        if (!imp) break;
    }

    decode_one_lut(wb_words, edges, N, cos_lut, sin_lut, u_out);
    pack_outer_f64(u_out, p_out, 1, N);
    return changed;
}

void build_trig_lut_ext(int brot, double* cos_lut, double* sin_lut) {
    build_trig_lut(brot, cos_lut, sin_lut);
}
