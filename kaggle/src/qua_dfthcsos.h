/*
 * qua_dfthcsos.h — Pure-Search Bilevel SEARU x FLOW 2-Jet Benchmark III (`dfthcsos`) Synthesizer
 *
 * Universal FTQC synthesis for FeMoCo DF-THC-SOS ANGLE_k XOR-oracle circuits:
 *   - Exact `SpinFreeData_V1` + `resource_costs` cost law matching `dfthcsos/verify.py`
 *     (including persistent `J * brot` angle registers, shared temporary workspace, and
 *      forward + reverse `U_T` + `U_T^dagger` Toffoli accounting)
 *   - Virtual-Slot-0 Select-Swap QROAM (`qua_bitop_emit_qroam_vslot0`) with constant-0 column
 *     pruning, collapsed unary spine, 0-CCX MBUC phase fixup, and provably clean `|0>` register
 *     column borrowing (`(w_prev_or | w_cur_or) == 0`)
 *   - Per-batch adaptive `(b_sh_k, collapse_k)` synthesis within the hard `Q <= 128` cap
 *   - Direct pure-C17 Kickmix v2 (.kmb) emission (`angle_k_0.kmb .. angle_k_{K-1}.kmb`)
 */
#ifndef QUA_DFTHCSOS_H
#define QUA_DFTHCSOS_H

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

#define QUA_DFTHC_NUM_ANGLES 75u
#define QUA_DFTHC_QUBIT_CAP  128u
#define QUA_DFTHC_MAX_J      16u
#define QUA_DFTHC_MAX_BSH    5u
#define QUA_DFTHC_NUM_COL    4u

static const uint32_t QUA_DFTHC_COLLAPSE_TAB[QUA_DFTHC_NUM_COL] = {0u, 3u, 4u, 5u};

typedef struct {
    uint32_t R;
    uint32_t B;
    uint32_t C;
    uint32_t brot;
    uint32_t num_rows;
    uint32_t addr_bits;
    uint16_t *words; /* num_rows * 75 */
    uint16_t col_or_mask[QUA_DFTHC_NUM_ANGLES + 1u];
    uint16_t pre_clean[QUA_DFTHC_MAX_J][QUA_DFTHC_NUM_ANGLES];
    uint16_t pre_anc[QUA_DFTHC_MAX_J][QUA_DFTHC_NUM_ANGLES][QUA_DFTHC_MAX_BSH][QUA_DFTHC_NUM_COL];
    uint32_t pre_tof[QUA_DFTHC_MAX_J][QUA_DFTHC_NUM_ANGLES][QUA_DFTHC_MAX_BSH][QUA_DFTHC_NUM_COL];
} qua_dfthc_data_t;

static inline uint16_t qua_dfthc_fk_word(
    const qua_dfthc_data_t *d,
    uint32_t row,
    uint32_t k,
    uint32_t J,
    uint32_t j_off
) {
    if (row >= d->num_rows) return 0u;
    uint32_t idx_cur = k * J + j_off;
    uint16_t w_cur = (idx_cur < QUA_DFTHC_NUM_ANGLES) ? d->words[row * QUA_DFTHC_NUM_ANGLES + idx_cur] : 0u;
    if (k == 0u) return w_cur;
    uint32_t idx_prev = (k - 1u) * J + j_off;
    uint16_t w_prev = (idx_prev < QUA_DFTHC_NUM_ANGLES) ? d->words[row * QUA_DFTHC_NUM_ANGLES + idx_prev] : 0u;
    return (uint16_t)(w_cur ^ w_prev);
}

static inline void qua_dfthc_precompute_tables(qua_dfthc_data_t *d) {
    for (uint32_t a = 0u; a < QUA_DFTHC_NUM_ANGLES; ++a) {
        uint16_t om = 0u;
        for (uint32_t r = 0u; r < d->num_rows; ++r) {
            om |= d->words[r * QUA_DFTHC_NUM_ANGLES + a];
        }
        d->col_or_mask[a] = om;
    }
    d->col_or_mask[QUA_DFTHC_NUM_ANGLES] = 0u;

    memset(d->pre_anc, 0xFF, sizeof(d->pre_anc));
    memset(d->pre_tof, 0xFF, sizeof(d->pre_tof));
    memset(d->pre_clean, 0, sizeof(d->pre_clean));

    for (uint32_t J = 1u; J <= QUA_DFTHC_MAX_J; ++J) {
        uint32_t j_idx = J - 1u;
        uint32_t w_bits = J * d->brot;
        if (d->addr_bits + w_bits > QUA_DFTHC_QUBIT_CAP) continue;

        uint32_t K = (QUA_DFTHC_NUM_ANGLES + J - 1u) / J;
        uint32_t stride_u64 = (w_bits + 63u) / 64u;
        if (stride_u64 < 1u) stride_u64 = 1u;
        uint64_t *tab = (uint64_t *)malloc((size_t)d->num_rows * stride_u64 * sizeof(uint64_t));
        if (!tab) continue;

        for (uint32_t k = 0u; k < K; ++k) {
            memset(tab, 0, (size_t)d->num_rows * stride_u64 * sizeof(uint64_t));
            uint16_t clean_cnt = 0u;
            uint32_t w_act = 0u;
            for (uint32_t j = 0u; j < J; ++j) {
                uint32_t idx_cur = k * J + j;
                uint32_t idx_prev = (k > 0u) ? ((k - 1u) * J + j) : QUA_DFTHC_NUM_ANGLES;
                uint16_t or_cur  = (idx_cur < QUA_DFTHC_NUM_ANGLES) ? d->col_or_mask[idx_cur] : 0u;
                uint16_t or_prev = (idx_prev < QUA_DFTHC_NUM_ANGLES) ? d->col_or_mask[idx_prev] : 0u;
                uint16_t both_or = (uint16_t)(or_cur | or_prev);
                for (uint32_t b = 0u; b < d->brot; ++b) {
                    if (((both_or >> b) & 1u) == 0u) ++clean_cnt;
                }
            }
            d->pre_clean[j_idx][k] = clean_cnt;

            for (uint32_t r = 0u; r < d->num_rows; ++r) {
                for (uint32_t j = 0u; j < J; ++j) {
                    uint16_t val = qua_dfthc_fk_word(d, r, k, J, j);
                    for (uint32_t b = 0u; b < d->brot; ++b) {
                        if ((val >> b) & 1u) {
                            uint32_t bit_idx = j * d->brot + b;
                            tab[(size_t)r * stride_u64 + (bit_idx >> 6u)] |= (UINT64_C(1) << (bit_idx & 63u));
                        }
                    }
                }
            }
            for (uint32_t bit_idx = 0u; bit_idx < w_bits; ++bit_idx) {
                uint32_t wi = bit_idx >> 6u, bi = bit_idx & 63u;
                for (uint32_t r = 0u; r < d->num_rows; ++r) {
                    if ((tab[(size_t)r * stride_u64 + wi] >> bi) & 1ull) {
                        ++w_act;
                        break;
                    }
                }
            }

            for (uint32_t b_sh = 0u; b_sh < QUA_DFTHC_MAX_BSH; ++b_sh) {
                if (b_sh > d->addr_bits) continue;
                if (w_act == 0u && b_sh > 0u) continue;
                uint32_t min_gather_q = (b_sh >= 2u) ? ((1u << b_sh) - 1u) : 0u;
                if (d->addr_bits + w_bits + min_gather_q > QUA_DFTHC_QUBIT_CAP + (uint32_t)clean_cnt) {
                    continue;
                }
                for (uint32_t ci = 0u; ci < QUA_DFTHC_NUM_COL; ++ci) {
                    qua_bitop_qroam_res_t qres;
                    qua_bitop_qroam_vslot0_resources(
                        tab, d->num_rows, stride_u64, d->addr_bits, w_bits,
                        b_sh, QUA_DFTHC_COLLAPSE_TAB[ci], 0u, &qres
                    );
                    uint32_t anc_q = qres.anc_needed;
                    uint64_t est_tof = qres.total_ccx;
                    d->pre_anc[j_idx][k][b_sh][ci] = (anc_q > 65534u) ? 65534u : (uint16_t)anc_q;
                    d->pre_tof[j_idx][k][b_sh][ci] = (est_tof > 0xFFFFFFFEull) ? 0xFFFFFFFEu : (uint32_t)est_tof;
                }
            }
        }
        free(tab);
    }
}

static inline bool qua_dfthc_load_words(const char *path, qua_dfthc_data_t *out) {
    memset(out, 0, sizeof(*out));
    FILE *f = fopen(path, "rb");
    if (!f) return false;
    int32_t hdr[4];
    if (fread(hdr, sizeof(int32_t), 4, f) != 4) { fclose(f); return false; }
    out->R = (uint32_t)hdr[0];
    out->B = (uint32_t)hdr[1];
    out->C = (uint32_t)hdr[2];
    out->brot = (uint32_t)hdr[3];
    out->num_rows = 76u + out->R * out->B;
    uint32_t ab = 1u;
    while ((1u << ab) < out->num_rows) ++ab;
    out->addr_bits = ab;
    size_t cnt = (size_t)out->num_rows * QUA_DFTHC_NUM_ANGLES;
    out->words = (uint16_t *)malloc(cnt * sizeof(uint16_t));
    if (!out->words || fread(out->words, sizeof(uint16_t), cnt, f) != cnt) {
        free(out->words);
        out->words = NULL;
        fclose(f);
        return false;
    }
    fclose(f);
    qua_dfthc_precompute_tables(out);
    return true;
}

static inline void qua_dfthc_free_data(qua_dfthc_data_t *d) {
    free(d->words);
    memset(d, 0, sizeof(*d));
}

typedef struct {
    uint32_t J;
    uint32_t K;
    uint32_t b_cap;
    uint32_t col_cap_idx;
    bool     borrow_clean;
    bool     adaptive_k;
    uint8_t  k_bsh[QUA_DFTHC_NUM_ANGLES];
    uint8_t  k_cidx[QUA_DFTHC_NUM_ANGLES];
    uint32_t q_per_circ;
    uint64_t est_qrom_tof_2pass;
    uint64_t est_heur_tof;
    uint64_t est_total_tof;
    uint64_t est_total_q;
    double   est_score;
    bool     valid;
} qua_dfthc_cfg_t;

static inline uint32_t qua_dfthc_ceil_log2(uint64_t x) {
    uint32_t r = 0u;
    uint64_t p = 1u;
    while (p < x) { p <<= 1u; ++r; }
    return r;
}

static inline void qua_dfthc_qr(uint64_t L, uint64_t M1, uint32_t *out_k, uint64_t *out_cost) {
    double k_f = 0.5 * log2((double)L / (double)M1);
    if (k_f < 0.0) k_f = 0.0;
    uint32_t k0 = (uint32_t)floor(k_f);
    uint32_t k1 = (uint32_t)ceil(k_f);
    double v0 = (double)L / (double)(1ULL << k0) + (double)M1 * (double)((1ULL << k0) - 1ULL);
    double v1 = (double)L / (double)(1ULL << k1) + (double)M1 * (double)((1ULL << k1) - 1ULL);
    uint32_t k_opt = (v1 < v0) ? k1 : k0;
    double v_opt = (v1 < v0) ? v1 : v0;
    *out_k = k_opt;
    *out_cost = (uint64_t)ceil(v_opt);
}

static inline uint64_t qua_dfthc_qi_cost(uint64_t L) {
    double k_f = 0.5 * log2((double)L);
    if (k_f < 0.0) k_f = 0.0;
    uint32_t k0 = (uint32_t)floor(k_f);
    uint32_t k1 = (uint32_t)ceil(k_f);
    double v0 = (double)L / (double)(1ULL << k0) + (double)((1ULL << k0) - 1ULL);
    double v1 = (double)L / (double)(1ULL << k1) + (double)((1ULL << k1) - 1ULL);
    return (uint64_t)ceil((v1 < v0) ? v1 : v0);
}

/* Exact SpinFreeData_V1 + resource_costs matching `dfthcsos/verify.py` */
static inline void qua_dfthc_spinfree_costs(
    uint32_t R, uint32_t B, uint32_t C, uint32_t N, uint32_t chi, uint32_t beta,
    uint32_t J, uint32_t addr_bits, uint32_t max_q_angles,
    uint64_t *out_tof, uint64_t *out_q
) {
    uint64_t X_o = (uint64_t)N + (uint64_t)R * (uint64_t)C;
    uint32_t lo = qua_dfthc_ceil_log2(X_o);
    uint32_t b1 = lo + chi;
    uint32_t k_o = 0u; uint64_t qr_o = 0u;
    qua_dfthc_qr(X_o, b1, &k_o, &qr_o);
    uint64_t outer_t = 8ull * lo + (qr_o + b1) + (qua_dfthc_qi_cost(X_o) + b1 + 1ull);

    uint32_t b2 = qua_dfthc_ceil_log2(B) + chi + 1u;
    uint32_t k2_q = 0u; uint64_t qr_i = 0u;
    qua_dfthc_qr((uint64_t)R * (uint64_t)C * (uint64_t)B, b2, &k2_q, &qr_i);
    for (uint32_t j = 0u; j <= k2_q; ++j) {
        uint32_t kn = k2_q - j;
        uint32_t b_blk = (B + (1u << kn) - 1u) >> kn;
        uint32_t qtmp = qua_dfthc_ceil_log2((uint64_t)R * (uint64_t)C * b_blk) + ((1u << kn) - 1u) * b2;
        if (qtmp < beta * (N - 1u)) { k2_q = kn; break; }
    }
    uint32_t k2_opt = k2_q;
    while (B < (1u << k2_opt) && k2_opt > 0u) --k2_opt;

    uint64_t RC = (uint64_t)R * (uint64_t)C;
    uint64_t for_qroam = RC * ((B + 1u + (1u << k2_opt) - 1u) >> k2_opt) + (uint64_t)(1u << k2_opt) * b2;
    uint32_t k4 = 0u; uint64_t qi_i = 0u;
    qua_dfthc_qr(RC, B + 1u, &k4, &qi_i);
    uint64_t for_qroam_dag = qi_i + (B + 1u) + N + b2;
    uint64_t unif_b = 4ull * qua_dfthc_ceil_log2(B + 1u);
    uint64_t inner_t = for_qroam + for_qroam_dag + 2ull * unif_b;
    uint64_t givens_t = 2ull * (N - 1u) * (2ull * beta - 1ull);
    *out_tof = outer_t + inner_t + givens_t;

    uint64_t outer_pers = (uint64_t)(b1 + chi) + (uint64_t)(lo + 2u);
    uint64_t outer_temp = (uint64_t)qua_dfthc_ceil_log2((X_o + (1ULL << k_o) - 1ULL) >> k_o) +
                          (1ULL << k_o) * (uint64_t)b1 + (uint64_t)(lo > 2u ? lo - 2u : 0u);
    uint64_t inner_pers = (uint64_t)(b2 + chi + qua_dfthc_ceil_log2(R) + 2u) +
                          (uint64_t)(qua_dfthc_ceil_log2(B + 1u) + 2u);
    uint64_t b_blk_q = (B + 1u + (1u << k2_q) - 1u) >> k2_q;
    uint64_t inner_temp = (uint64_t)(qua_dfthc_ceil_log2(B + 1u) > 2u ? qua_dfthc_ceil_log2(B + 1u) - 2u : 0u) +
                          (uint64_t)qua_dfthc_ceil_log2(RC * b_blk_q) + ((1ULL << k2_q) - 1ULL) * (uint64_t)b2;

    uint64_t angles_pers = (uint64_t)J * (uint64_t)beta;
    uint64_t pub_q = (uint64_t)addr_bits + angles_pers;
    uint64_t angles_temp = ((uint64_t)max_q_angles > pub_q) ? ((uint64_t)max_q_angles - pub_q) : 0ull;

    uint64_t max_temp = outer_temp;
    if (inner_temp > max_temp) max_temp = inner_temp;
    if (angles_temp > max_temp) max_temp = angles_temp;
    *out_q = outer_pers + inner_pers + angles_pers + max_temp;
}

static inline void qua_dfthc_eval_cfg(
    const qua_dfthc_data_t *d,
    uint64_t word0,
    qua_dfthc_cfg_t *cfg
) {
    memset(cfg, 0, sizeof(*cfg));
    uint32_t j_raw       = (uint32_t)(word0 & 15u);          /* [0..3]: 1..16 */
    uint32_t b_raw       = (uint32_t)((word0 >> 4u) & 7u);   /* [4..6]: 0..4  */
    uint32_t c_raw       = (uint32_t)((word0 >> 7u) & 3u);   /* [7..8]: 0..3  */
    bool     borrow      = ((word0 >> 9u) & 1u) != 0u;       /* [9]           */
    bool     adaptive    = ((word0 >> 10u) & 1u) != 0u;      /* [10]          */

    uint32_t J = 1u + j_raw;
    uint32_t b_cap = b_raw % QUA_DFTHC_MAX_BSH;
    uint32_t col_cap_idx = c_raw;
    uint32_t K = (QUA_DFTHC_NUM_ANGLES + J - 1u) / J;
    uint32_t w_bits = J * d->brot;

    cfg->J = J;
    cfg->K = K;
    cfg->b_cap = b_cap;
    cfg->col_cap_idx = col_cap_idx;
    cfg->borrow_clean = borrow;
    cfg->adaptive_k = adaptive;

    if (d->addr_bits + w_bits > QUA_DFTHC_QUBIT_CAP) {
        cfg->q_per_circ = d->addr_bits + w_bits;
        cfg->valid = false;
        cfg->est_score = 1e30;
        return;
    }

    uint32_t j_idx = J - 1u;
    uint32_t max_q_circ = 0u;
    uint64_t sum_1pass_tof = 0u;
    bool all_ok = true;

    for (uint32_t k = 0u; k < K; ++k) {
        uint32_t clean_free = borrow ? (uint32_t)d->pre_clean[j_idx][k] : 0u;
        uint32_t best_tof = 0xFFFFFFFFu;
        uint32_t best_qc  = 0xFFFFFFFFu;
        uint32_t best_b   = 0u;
        uint32_t best_ci  = 0u;
        uint32_t min_qc_any = 0xFFFFFFFFu;

        uint32_t b_lo  = adaptive ? 0u : b_cap;
        uint32_t ci_lo = adaptive ? 0u : col_cap_idx;
        for (uint32_t b_sh = b_lo; b_sh <= b_cap; ++b_sh) {
            for (uint32_t ci = ci_lo; ci <= col_cap_idx; ++ci) {
                uint16_t anc_q = d->pre_anc[j_idx][k][b_sh][ci];
                uint32_t tof_k = d->pre_tof[j_idx][k][b_sh][ci];
                if (anc_q == 0xFFFFu) continue;
                uint32_t extra_anc = ((uint32_t)anc_q > clean_free) ? ((uint32_t)anc_q - clean_free) : 0u;
                uint32_t qc = d->addr_bits + w_bits + extra_anc;
                if (qc < min_qc_any) min_qc_any = qc;
                if (qc <= QUA_DFTHC_QUBIT_CAP) {
                    if (tof_k < best_tof || (tof_k == best_tof && qc < best_qc)) {
                        best_tof = tof_k;
                        best_qc  = qc;
                        best_b   = b_sh;
                        best_ci  = ci;
                    }
                }
            }
        }
        if (best_tof == 0xFFFFFFFFu) {
            all_ok = false;
            uint32_t fallback_qc = (min_qc_any < 0xFFFFFFFFu) ? min_qc_any : (QUA_DFTHC_QUBIT_CAP + 16u);
            if (fallback_qc > max_q_circ) max_q_circ = fallback_qc;
        } else {
            cfg->k_bsh[k]  = (uint8_t)best_b;
            cfg->k_cidx[k] = (uint8_t)best_ci;
            if (best_qc > max_q_circ) max_q_circ = best_qc;
            sum_1pass_tof += (uint64_t)best_tof;
        }
    }

    cfg->q_per_circ = max_q_circ;
    if (!all_ok || max_q_circ > QUA_DFTHC_QUBIT_CAP) {
        cfg->valid = false;
        cfg->est_score = 1e30;
        return;
    }

    uint64_t heur_tof = 0u, total_q = 0u;
    qua_dfthc_spinfree_costs(d->R, d->B, d->C, 76u, 10u, d->brot, J, d->addr_bits, max_q_circ, &heur_tof, &total_q);

    cfg->est_qrom_tof_2pass = 2ull * sum_1pass_tof;
    cfg->est_heur_tof       = heur_tof;
    cfg->est_total_tof      = cfg->est_qrom_tof_2pass + heur_tof;
    cfg->est_total_q        = total_q;
    cfg->est_score          = (double)cfg->est_total_tof * (double)cfg->est_total_q;
    cfg->valid              = true;
}

static inline int qua_dfthc_emit_family(
    const qua_dfthc_data_t *d,
    const qua_dfthc_cfg_t *cfg,
    const char *out_dir
) {
    uint64_t total_ccx_1pass = 0u;
    uint32_t max_q = 0u;
    for (uint32_t old_k = 0u; old_k < QUA_DFTHC_NUM_ANGLES; ++old_k) {
        char old_path[512];
        snprintf(old_path, sizeof(old_path), "%s/angle_k_%u.kmb", out_dir, old_k);
        remove(old_path);
    }
    uint32_t j_idx = cfg->J - 1u;
    for (uint32_t k = 0u; k < cfg->K; ++k) {
        qua_builder_t kb;
        qua_builder_init(&kb);

        uint32_t sel_q[16];
        for (uint32_t i = 0u; i < d->addr_bits; ++i) sel_q[i] = i;
        qua_add_quantum_reg(&kb, "selection", sel_q, d->addr_bits);

        uint32_t word_qs[QUA_DFTHC_NUM_ANGLES * 16u];
        uint32_t q_cursor = d->addr_bits;
        for (uint32_t j = 0u; j < cfg->J; ++j) {
            uint32_t wq[16];
            for (uint32_t b = 0u; b < d->brot; ++b) {
                wq[b] = q_cursor;
                word_qs[j * d->brot + b] = q_cursor++;
            }
            char wname[32];
            snprintf(wname, sizeof(wname), "word_%u", j);
            qua_add_quantum_reg(&kb, wname, wq, d->brot);
        }

        uint32_t w_bits = cfg->J * d->brot;
        uint32_t stride_u64 = (w_bits + 63u) / 64u;
        if (stride_u64 < 1u) stride_u64 = 1u;
        uint64_t *tab = (uint64_t *)calloc((size_t)d->num_rows * stride_u64, sizeof(uint64_t));
        if (!tab) abort();
        for (uint32_t r = 0u; r < d->num_rows; ++r) {
            for (uint32_t j = 0u; j < cfg->J; ++j) {
                uint16_t val = qua_dfthc_fk_word(d, r, k, cfg->J, j);
                for (uint32_t b = 0u; b < d->brot; ++b) {
                    if ((val >> b) & 1u) {
                        uint32_t bit_idx = j * d->brot + b;
                        tab[(size_t)r * stride_u64 + (bit_idx >> 6u)] |= (UINT64_C(1) << (bit_idx & 63u));
                    }
                }
            }
        }

        uint32_t b_sh_k = cfg->k_bsh[k];
        uint32_t c_idx_k = cfg->k_cidx[k];
        uint32_t collapse_k = QUA_DFTHC_COLLAPSE_TAB[c_idx_k];
        uint32_t anc_need = d->pre_anc[j_idx][k][b_sh_k][c_idx_k];

        uint32_t anc_pool[512];
        uint32_t anc_have = 0u;
        if (cfg->borrow_clean && anc_need > 0u) {
            for (uint32_t j = 0u; j < cfg->J && anc_have < anc_need; ++j) {
                uint32_t idx_cur = k * cfg->J + j;
                uint32_t idx_prev = (k > 0u) ? ((k - 1u) * cfg->J + j) : QUA_DFTHC_NUM_ANGLES;
                uint16_t or_cur  = (idx_cur < QUA_DFTHC_NUM_ANGLES) ? d->col_or_mask[idx_cur] : 0u;
                uint16_t or_prev = (idx_prev < QUA_DFTHC_NUM_ANGLES) ? d->col_or_mask[idx_prev] : 0u;
                uint16_t both_or = (uint16_t)(or_cur | or_prev);
                for (uint32_t b = 0u; b < d->brot && anc_have < anc_need; ++b) {
                    if (((both_or >> b) & 1u) == 0u) {
                        anc_pool[anc_have++] = word_qs[j * d->brot + b];
                    }
                }
            }
        }
        while (anc_have < anc_need) anc_pool[anc_have++] = q_cursor++;

        uint32_t next_cbit = 0u;
        qua_bitop_emit_qroam_vslot0(
            &kb, tab, d->num_rows, stride_u64, sel_q, d->addr_bits,
            word_qs, w_bits, b_sh_k, collapse_k, anc_pool, anc_need, &next_cbit
        );
        free(tab);

        char path[512];
        snprintf(path, sizeof(path), "%s/angle_k_%u.kmb", out_dir, k);
        qua_bitop_spec_t fin_spec = {
            "dfthcsos", QUA_COST_QT, false, QUA_DFTHC_QUBIT_CAP, 11u, 1e5, 64u, 1024u, 4u, NULL, NULL, NULL
        };
        if (qua_bitop_finalize_kmb(&kb, &fin_spec, path) != 0) {
            qua_builder_free(&kb);
            return 1;
        }
        if (kb.num_qubits > max_q) max_q = kb.num_qubits;
        total_ccx_1pass += kb.count_toffoli;
        qua_builder_free(&kb);
    }
    printf("[DF-THC Unified C17] Emitted %u circuits into %s | max_Q=%u<=128, 1pass_CCX=%llu, 2pass_CCX=%llu\n",
           cfg->K, out_dir, max_q,
           (unsigned long long)total_ccx_1pass,
           (unsigned long long)(2ull * total_ccx_1pass));
    return 0;
}

static inline void qua_dfthc_populate_dag(
    const qua_dfthc_data_t *d,
    const qua_dfthc_cfg_t *c,
    qua_bitop_dag_t *dag
) {
    if (!dag) return;
    uint32_t pub_q = d->addr_bits + c->J * d->brot;
    dag->num_nodes      = 0u;
    dag->base_io_qubits = pub_q;
    dag->cost_law       = QUA_COST_QT;
    uint32_t ang_tmp = (c->q_per_circ > pub_q) ? (c->q_per_circ - pub_q) : 0u;
    uint32_t spine_q = (ang_tmp > 8u) ? 8u : ang_tmp;
    uint32_t bus_q   = (ang_tmp > spine_q) ? (ang_tmp - spine_q) : 0u;
    uint32_t n_io = qua_bitop_dag_push_rep(dag, QUA_BOP_IO,   QUA_REP_KET,         0u,      0u, 0u, 0u, false, true,  0u);
    uint32_t n_sp = qua_bitop_dag_push_rep(dag, QUA_BOP_AND,  QUA_REP_KET,         spine_q, (uint32_t)(c->est_qrom_tof_2pass / 4ull), 0u, 0u, true, false, (UINT64_C(1) << n_io));
    uint32_t n_qr = qua_bitop_dag_push_rep(dag, QUA_BOP_QROM, QUA_REP_PARITY_DUAL, bus_q,   (uint32_t)(c->est_qrom_tof_2pass * 3ull / 4ull), 16u, 0u, true, false, (UINT64_C(1) << n_sp));
    (void)qua_bitop_dag_push_rep(dag, QUA_BOP_XOR, QUA_REP_KET, 0u, 1u, 0u, 0u, false, true, (UINT64_C(1) << n_io) | (UINT64_C(1) << n_qr));
    qua_ele_pebble_plan_t pp;
    qua_bitop_eval_dag_pebble(dag, (UINT64_C(1) << n_sp) | (UINT64_C(1) << n_qr), &pp, NULL, 0u);
}

static inline bool qua_dfthc_eval_spec(
    const uint64_t *genome,
    qua_bitop_dag_t *dag,
    double *out_raw_cost,
    uint64_t *out_frust,
    void *user_ctx
) {
    const qua_dfthc_data_t *d = (const qua_dfthc_data_t *)user_ctx;
    qua_dfthc_cfg_t c;
    qua_dfthc_eval_cfg(d, genome[0], &c);
    if (!c.valid) {
        if (dag) {
            uint32_t over = (c.q_per_circ > QUA_DFTHC_QUBIT_CAP) ? (c.q_per_circ - QUA_DFTHC_QUBIT_CAP) : 1u;
            dag->q_violation = (double)over;
        }
        /* Exceeded 128-qubit cap: steer J [0..3], b_cap [4..6], borrow [9], adaptive [10] */
        *out_frust = UINT64_C(0x67F);
        return false;
    }
    qua_dfthc_populate_dag(d, &c, dag);
    *out_raw_cost = c.est_score;
    *out_frust = UINT64_C(0x7FF);
    return true;
}

static inline int qua_dfthc_run_pure_search(const char *words_bin_path, const char *out_dir) {
    qua_dfthc_data_t d;
    if (!qua_dfthc_load_words(words_bin_path, &d) &&
        !qua_dfthc_load_words("dfthcsos_opt/words.bin", &d)) {
        fprintf(stderr, "[ERROR] DF-THC words file not found at %s (or dfthcsos_opt/words.bin)\n",
                words_bin_path);
        return 1;
    }

    qua_bitop_spec_t spec = {
        "dfthcsos", QUA_COST_QT, false, QUA_DFTHC_QUBIT_CAP, 11u, 1e5, 64u, 1024u, 4u, qua_dfthc_eval_spec, &d, NULL
    };
    uint64_t best_g[4] = {0u, 0u, 0u, 0u};
    double best_raw = 1e30;
    qua_bitop_solve(&spec, UINT64_C(0xDF72026), best_g, &best_raw);

    qua_dfthc_cfg_t best;
    qua_dfthc_eval_cfg(&d, best_g[0], &best);

    printf("[DF-THC Pure-Search Winner] J=%u, K=%u, B_cap=%u, col=%u, borrow=%d, adapt=%d | Q_circ=%u<=128, Q_tot=%llu, 2pass_QROM_T=%llu, T_tot=%llu -> est_Q*T=%.6e\n",
           best.J, best.K, 1u << best.b_cap, QUA_DFTHC_COLLAPSE_TAB[best.col_cap_idx],
           (int)best.borrow_clean, (int)best.adaptive_k,
           best.q_per_circ, (unsigned long long)best.est_total_q,
           (unsigned long long)best.est_qrom_tof_2pass,
           (unsigned long long)best.est_total_tof, best.est_score);

    int rc = qua_dfthc_emit_family(&d, &best, out_dir);
    qua_dfthc_free_data(&d);
    return rc;
}

#endif /* QUA_DFTHCSOS_H */
