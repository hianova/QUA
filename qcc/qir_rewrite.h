/*
 * qir_rewrite.h — Bounded Arena E-Graph & Composable Quantum Circuit Rewrite Rules (R1..R7)
 * Part of QCC (Quantum Circuit Compiler) for FLOW v7.0.1.
 * Pure C17, zero external dependencies, zero indirect function calls.
 */
#ifndef QIR_REWRITE_H
#define QIR_REWRITE_H

#include "qir_lift.h"

#ifdef __cplusplus
extern "C" {
#endif

#define QIR_EGRAPH_MAX_CLASSES 4096u
#define QIR_GF12_RANK          42u
#define QIR_GF64_RANK          15u

typedef enum {
    QIR_VAR_STD = 0,               /* Direct standard lowering */
    QIR_VAR_LINF2_GAUSS,           /* R1: Fused affine F2 + min-fanout Gauss-Jordan */
    QIR_VAR_HMR_PARITY_LUT,        /* R2: HMR dual (Parity-LUT or 0-CCX HMR uncompute) */
    QIR_VAR_CSA_TREE,              /* R3: Carry-Save 3:2 reduction tree over >=2 ADDs */
    QIR_VAR_PHASE_CSA_DIRECT,      /* R4: Elide ADD/SUB tree feeding linear modular PHASE_POLY */
    QIR_VAR_LUT_COL_SPLIT,         /* R5: Multi-column streaming bit-sliced LUT */
    QIR_VAR_GF2M_CHUDNOVSKY_R42,   /* R6: Rank-42 P^1(GF(8)) GF(2^12) bilinear mul */
    QIR_VAR_GF2M_RANK15,           /* R6: Rank-15 P^1(GF(4)) GF(64) bilinear mul */
    QIR_VAR_GF2M_TOWER_INV,        /* R6: Composite tower GF(((2^3)^2)^2) inversion */
    QIR_VAR_LIN_GAUGE_TSP,         /* R7: Rank-compressed / Hamming-TSP basis-gauged LUT */
    QIR_VAR_TRIVIAL_ELIDE,         /* Algebraic zero/const/dead elimination */
    QIR_VAR_QROAM_VSLOT0,          /* R8: Virtual-Slot-0 Differential Select-Swap QROAM */
    QIR_VAR_MUL_CONST_ADD,         /* R9: Shift-add constant multiplier (0 qubits for constant) */
    QIR_VAR_SQR_TRIANGLE,          /* R9: Disjoint-control diagonal+triangle squarer (0 copy Q) */
    QIR_VAR_QROAM_SPLIT,           /* R8b: Column-split Virtual-Slot-0 QROAM */
    QIR_VAR_INPLACE_ARITH,         /* R10a: In-place accumulator on ADD/SUB/GF2M_ADD */
    QIR_VAR_CTRL_ADD_CUCCARO,      /* R10b: Controlled +/- via Cuccaro 1-ancilla adder */
    QIR_VAR_XOR_ACCUM,             /* R10c: Fused out ^= op(...) zero-copy accumulator */
    QIR_VAR_GF2M_CONST_LINF2,      /* R11: Constant GF(2^m) multiply folded into GL_m(F_2) CX */
    QIR_VAR_KIND_COUNT
} qir_variant_kind_t;

typedef struct {
    uint32_t  target_val;
    uint16_t  bit_lo;
    uint16_t  width;
    flow_u128 step;
    uint32_t  leaf_vals[16];
    int8_t    leaf_signs[16];
    uint16_t  num_leaves;
} qir_phase_split_t;

typedef struct {
    uint16_t l_mask[QIR_GF12_RANK];
    uint16_t r_mask[QIR_GF12_RANK];
    uint16_t w_mask[QIR_GF12_RANK];
    uint32_t order[QIR_GF12_RANK];
} qir_gf12_rank42_t;

typedef struct {
    uint8_t  l_mask[QIR_GF64_RANK];
    uint8_t  r_mask[QIR_GF64_RANK];
    uint8_t  w_mask[QIR_GF64_RANK];
    uint32_t order[QIR_GF64_RANK];
} qir_gf64_rank15_t;

typedef struct {
    uint16_t     rank_bits;   /* Compressed column rank r <= val_bits */
    uint16_t     val_bits;
    uint16_t     val_words;
    uint32_t     num_entries;
    qcc_bv256_t *recon_rows;  /* [val_bits]: row i reconstructs orig bit i from r compressed cols + const_shift */
    qcc_bv256_t  const_shift; /* orig_entry[0] subtracted to zero-base */
    uint64_t    *comp_table;  /* [num_entries * ((rank_bits+63)/64)] */
    uint16_t     comp_words;
} qir_gauge_lut_t;

typedef struct {
    uint32_t           id;
    uint32_t           eclass_id;
    qir_variant_kind_t variant;
    const qir_op_t    *op;
    uint32_t           child_eclasses[8];
    uint8_t            num_children;
    uint32_t           aux_param;       /* e.g. split cols K, b_sh, const_idx, root op id */
    uint16_t           approx_keep_bits;/* Flow-tuned active bitwidth when QIR_FLAG_APPROX */
    uint8_t            rule_tag;        /* 0=base, 1..9 = R1..R9 */
    /* Rule-specific cached synthesis artifacts (arena-owned) */
    qir_linf2_payload_t      *fused_linf2;
    qir_phase_poly_payload_t *split_phase;
    qir_phase_split_t        *phase_split;
    qir_gauge_lut_t          *gauge_lut;
    uint32_t                 *chain_ops;      /* Folded intermediate op IDs (R1/R3/R4/R9) */
    uint16_t                  num_chain_ops;
} qir_enode_t;

typedef struct {
    uint32_t id;
    uint32_t val_id;
    uint16_t width;
    uint32_t *enodes;
    uint32_t num_enodes;
    uint32_t cap_enodes;
    uint32_t use_count;
} qir_eclass_t;

typedef struct {
    qcc_arena_t        *arena;
    const qir_module_t *mod;
    qir_eclass_t       *classes;
    uint32_t            num_classes;
    qir_enode_t        *enodes;
    uint32_t            num_enodes;
    uint32_t            cap_enodes;
    uint32_t            rules_fired[16]; /* Counts per rule R1..R15 */
} qir_egraph_t;

/* ============================================================================
 * GF(2^m) Tower Algebra Helpers (Zero Indirect Function Calls)
 * ============================================================================ */

static inline uint8_t qir__gf4_mul(uint8_t a, uint8_t b) {
    uint8_t r = 0u;
    for (uint32_t i = 0u; i < 2u; ++i) {
        if ((b >> i) & 1u) r ^= a;
        a = (uint8_t)((a << 1u) ^ ((a & 2u) ? 0x7u : 0u));
    }
    return (uint8_t)(r & 3u);
}
static inline uint8_t qir__gf4_inv(uint8_t a) { return qir__gf4_mul(a, a); }

static inline uint8_t qir__gf8_mul(uint8_t a, uint8_t b) {
    uint8_t r = 0u;
    for (uint32_t i = 0u; i < 3u; ++i) {
        if ((b >> i) & 1u) r ^= a;
        a = (uint8_t)((a << 1u) ^ ((a & 4u) ? 0xBu : 0u));
    }
    return (uint8_t)(r & 7u);
}
static inline uint8_t qir__gf8_inv(uint8_t a) {
    uint8_t a2 = qir__gf8_mul(a, a);
    uint8_t a3 = qir__gf8_mul(a2, a);
    return qir__gf8_mul(a3, a3); /* a^6 = a^{-1} in GF(8) */
}

static inline uint32_t qir__hamming_dist_u32(uint32_t a, uint32_t b) {
    return (uint32_t)__builtin_popcount((unsigned)(a ^ b));
}

static inline bool qir_build_gf12_rank42(uint32_t irred_poly, qir_gf12_rank42_t *out) {
    memset(out, 0, sizeof(*out));
    uint16_t t_elem = 0u;
    for (uint16_t cand = 2u; cand < 4096u; ++cand) {
        uint16_t c2 = (uint16_t)qir__gf2m_mul_u32(cand, cand, 12u, irred_poly);
        uint16_t c3 = (uint16_t)qir__gf2m_mul_u32(c2, cand, 12u, irred_poly);
        if ((c3 ^ cand ^ 1u) == 0u) {
            t_elem = cand;
            break;
        }
    }
    if (t_elem == 0u) return false;

    uint16_t t_pow[3] = {1u, t_elem, (uint16_t)qir__gf2m_mul_u32(t_elem, t_elem, 12u, irred_poly)};
    uint16_t basis[12], inv_cols[12];
    bool found_basis = false;

    for (uint16_t w_cand = 2u; w_cand < 4096u && !found_basis; ++w_cand) {
        uint16_t w_pow = 1u;
        for (uint32_t j = 0u; j < 4u; ++j) {
            for (uint32_t i = 0u; i < 3u; ++i) {
                basis[3u * j + i] = (uint16_t)qir__gf2m_mul_u32(t_pow[i], w_pow, 12u, irred_poly);
            }
            w_pow = (uint16_t)qir__gf2m_mul_u32(w_pow, w_cand, 12u, irred_poly);
        }
        uint16_t mat[12], inv[12];
        for (uint32_t r = 0u; r < 12u; ++r) {
            uint16_t row = 0u;
            for (uint32_t c = 0u; c < 12u; ++c) {
                if ((basis[c] >> r) & 1u) row |= (uint16_t)(1u << c);
            }
            mat[r] = row;
            inv[r] = (uint16_t)(1u << r);
        }
        bool ok = true;
        for (uint32_t col = 0u; col < 12u; ++col) {
            uint32_t piv = col;
            while (piv < 12u && (((mat[piv] >> col) & 1u) == 0u)) ++piv;
            if (piv == 12u) { ok = false; break; }
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
        if (ok) {
            for (uint32_t coord = 0u; coord < 12u; ++coord) inv_cols[coord] = inv[coord];
            found_basis = true;
            uint16_t w7[7];
            w7[0] = 1u;
            for (uint32_t d = 1u; d < 7u; ++d) w7[d] = (uint16_t)qir__gf2m_mul_u32(w7[d - 1u], w_cand, 12u, irred_poly);

            static const uint8_t pts[6] = {0u, 1u, 2u, 4u, 3u, 6u};
            uint8_t vmat[7][7], vinv[7][7];
            memset(vmat, 0, sizeof(vmat));
            memset(vinv, 0, sizeof(vinv));
            for (uint32_t r = 0u; r < 6u; ++r) {
                uint8_t pp = 1u;
                for (uint32_t c = 0u; c < 7u; ++c) {
                    vmat[r][c] = pp;
                    pp = qir__gf8_mul(pp, pts[r]);
                }
                vinv[r][r] = 1u;
            }
            vmat[6][6] = 1u;
            vinv[6][6] = 1u;

            for (uint32_t col = 0u; col < 7u; ++col) {
                uint32_t piv = col;
                while (piv < 7u && vmat[piv][col] == 0u) ++piv;
                if (piv == 7u) return false;
                if (piv != col) {
                    for (uint32_t k = 0u; k < 7u; ++k) {
                        uint8_t tv = vmat[col][k]; vmat[col][k] = vmat[piv][k]; vmat[piv][k] = tv;
                        uint8_t ti = vinv[col][k]; vinv[col][k] = vinv[piv][k]; vinv[piv][k] = ti;
                    }
                }
                uint8_t sc = qir__gf8_inv(vmat[col][col]);
                for (uint32_t k = 0u; k < 7u; ++k) {
                    vmat[col][k] = qir__gf8_mul(vmat[col][k], sc);
                    vinv[col][k] = qir__gf8_mul(vinv[col][k], sc);
                }
                for (uint32_t r = 0u; r < 7u; ++r) {
                    if (r != col && vmat[r][col] != 0u) {
                        uint8_t f = vmat[r][col];
                        for (uint32_t k = 0u; k < 7u; ++k) {
                            vmat[r][k] ^= qir__gf8_mul(f, vmat[col][k]);
                            vinv[r][k] ^= qir__gf8_mul(f, vinv[col][k]);
                        }
                    }
                }
            }

            static const uint8_t k6_u[6]   = {0x1u, 0x2u, 0x4u, 0x3u, 0x5u, 0x6u};
            static const uint8_t k6_v[6]   = {0x1u, 0x2u, 0x4u, 0x3u, 0x5u, 0x6u};
            static const uint8_t k6_out[6] = {0x7u, 0x5u, 0x1u, 0x2u, 0x4u, 0x3u};

            for (uint32_t m = 0u; m < 7u; ++m) {
                uint16_t u_lin[3] = {0}, v_lin[3] = {0};
                if (m < 6u) {
                    uint8_t pp = 1u;
                    for (uint32_t j = 0u; j < 4u; ++j) {
                        for (uint32_t i = 0u; i < 3u; ++i) {
                            uint8_t pr = qir__gf8_mul((uint8_t)(1u << i), pp);
                            for (uint32_t b = 0u; b < 3u; ++b) {
                                if ((pr >> b) & 1u) {
                                    u_lin[b] ^= inv_cols[3u * j + i];
                                    v_lin[b] ^= inv_cols[3u * j + i];
                                }
                            }
                        }
                        pp = qir__gf8_mul(pp, pts[m]);
                    }
                } else {
                    for (uint32_t b = 0u; b < 3u; ++b) {
                        u_lin[b] = inv_cols[9u + b];
                        v_lin[b] = inv_cols[9u + b];
                    }
                }
                for (uint32_t k = 0u; k < 6u; ++k) {
                    uint32_t idx = 6u * m + k;
                    uint16_t lm = 0u, rm = 0u, wm = 0u;
                    for (uint32_t b = 0u; b < 3u; ++b) {
                        if ((k6_u[k] >> b) & 1u) lm ^= u_lin[b];
                        if ((k6_v[k] >> b) & 1u) rm ^= v_lin[b];
                    }
                    for (uint32_t d = 0u; d < 7u; ++d) {
                        uint8_t cg = qir__gf8_mul(vinv[d][m], k6_out[k]);
                        for (uint32_t i = 0u; i < 3u; ++i) {
                            if ((cg >> i) & 1u) wm ^= (uint16_t)qir__gf2m_mul_u32(t_pow[i], w7[d], 12u, irred_poly);
                        }
                    }
                    out->l_mask[idx] = lm;
                    out->r_mask[idx] = rm;
                    out->w_mask[idx] = wm;
                    out->order[idx]  = idx;
                }
            }
        }
    }
    if (!found_basis) return false;

    for (uint32_t i = 0u; i < 12u; ++i) {
        for (uint32_t j = 0u; j < 12u; ++j) {
            uint16_t a = (uint16_t)(1u << i), b = (uint16_t)(1u << j), got = 0u;
            for (uint32_t p = 0u; p < QIR_GF12_RANK; ++p) {
                uint32_t la = (uint32_t)__builtin_parity((unsigned)(out->l_mask[p] & a));
                uint32_t rb = (uint32_t)__builtin_parity((unsigned)(out->r_mask[p] & b));
                if (la & rb) got ^= out->w_mask[p];
            }
            if (got != (uint16_t)qir__gf2m_mul_u32(a, b, 12u, irred_poly)) return false;
        }
    }

    bool used[QIR_GF12_RANK];
    memset(used, 0, sizeof(used));
    uint32_t cur = 0u;
    for (uint32_t step = 0u; step < QIR_GF12_RANK; ++step) {
        uint32_t best_p = 0u, best_cost = 9999u;
        for (uint32_t p = 0u; p < QIR_GF12_RANK; ++p) {
            if (used[p]) continue;
            uint32_t c = (step == 0u)
                ? (uint32_t)(__builtin_popcount(out->l_mask[p]) + __builtin_popcount(out->r_mask[p]))
                : (qir__hamming_dist_u32(out->l_mask[cur], out->l_mask[p]) +
                   qir__hamming_dist_u32(out->r_mask[cur], out->r_mask[p]));
            if (c < best_cost) { best_cost = c; best_p = p; }
        }
        used[best_p] = true;
        out->order[step] = best_p;
        cur = best_p;
    }
    return true;
}

static inline bool qir_build_gf64_rank15(uint32_t irred_poly, qir_gf64_rank15_t *out) {
    memset(out, 0, sizeof(*out));
    uint8_t w_elem = 0u;
    for (uint8_t cand = 2u; cand < 64u; ++cand) {
        uint8_t c2 = (uint8_t)qir__gf2m_mul_u32(cand, cand, 6u, irred_poly);
        if ((c2 ^ cand ^ 1u) == 0u) { w_elem = cand; break; }
    }
    if (w_elem == 0u) return false;

    uint8_t w_pow[2] = {1u, w_elem};
    uint8_t inv_cols[6], z5[5];
    bool found_basis = false;

    for (uint8_t z_cand = 2u; z_cand < 64u && !found_basis; ++z_cand) {
        uint8_t z2 = (uint8_t)qir__gf2m_mul_u32(z_cand, z_cand, 6u, irred_poly);
        uint8_t z3 = (uint8_t)qir__gf2m_mul_u32(z2, z_cand, 6u, irred_poly);
        if ((z3 ^ z_cand ^ 1u) != 0u) continue;
        uint8_t basis[6], zp = 1u;
        for (uint32_t j = 0u; j < 3u; ++j) {
            for (uint32_t i = 0u; i < 2u; ++i) {
                basis[2u * j + i] = (uint8_t)qir__gf2m_mul_u32(w_pow[i], zp, 6u, irred_poly);
            }
            zp = (uint8_t)qir__gf2m_mul_u32(zp, z_cand, 6u, irred_poly);
        }
        uint8_t mat[6], inv[6];
        for (uint32_t r = 0u; r < 6u; ++r) {
            uint8_t row = 0u;
            for (uint32_t c = 0u; c < 6u; ++c) if ((basis[c] >> r) & 1u) row |= (uint8_t)(1u << c);
            mat[r] = row;
            inv[r] = (uint8_t)(1u << r);
        }
        bool ok = true;
        for (uint32_t col = 0u; col < 6u; ++col) {
            uint32_t piv = col;
            while (piv < 6u && (((mat[piv] >> col) & 1u) == 0u)) ++piv;
            if (piv == 6u) { ok = false; break; }
            if (piv != col) {
                uint8_t tm = mat[col]; mat[col] = mat[piv]; mat[piv] = tm;
                uint8_t ti = inv[col]; inv[col] = inv[piv]; inv[piv] = ti;
            }
            for (uint32_t r = 0u; r < 6u; ++r) {
                if (r != col && ((mat[r] >> col) & 1u)) {
                    mat[r] ^= mat[col];
                    inv[r] ^= inv[col];
                }
            }
        }
        if (ok) {
            for (uint32_t c = 0u; c < 6u; ++c) inv_cols[c] = inv[c];
            z5[0] = 1u;
            for (uint32_t d = 1u; d < 5u; ++d) z5[d] = (uint8_t)qir__gf2m_mul_u32(z5[d - 1u], z_cand, 6u, irred_poly);
            found_basis = true;
        }
    }
    if (!found_basis) return false;

    static const uint8_t pts[4] = {0u, 1u, 2u, 3u};
    uint8_t vmat[5][5], vinv[5][5];
    memset(vmat, 0, sizeof(vmat));
    memset(vinv, 0, sizeof(vinv));
    for (uint32_t r = 0u; r < 4u; ++r) {
        uint8_t pp = 1u;
        for (uint32_t c = 0u; c < 5u; ++c) { vmat[r][c] = pp; pp = qir__gf4_mul(pp, pts[r]); }
        vinv[r][r] = 1u;
    }
    vmat[4][4] = 1u; vinv[4][4] = 1u;
    for (uint32_t col = 0u; col < 5u; ++col) {
        uint32_t piv = col;
        while (piv < 5u && vmat[piv][col] == 0u) ++piv;
        if (piv == 5u) return false;
        if (piv != col) {
            for (uint32_t k = 0u; k < 5u; ++k) {
                uint8_t tv = vmat[col][k]; vmat[col][k] = vmat[piv][k]; vmat[piv][k] = tv;
                uint8_t ti = vinv[col][k]; vinv[col][k] = vinv[piv][k]; vinv[piv][k] = ti;
            }
        }
        uint8_t sc = qir__gf4_inv(vmat[col][col]);
        for (uint32_t k = 0u; k < 5u; ++k) {
            vmat[col][k] = qir__gf4_mul(vmat[col][k], sc);
            vinv[col][k] = qir__gf4_mul(vinv[col][k], sc);
        }
        for (uint32_t r = 0u; r < 5u; ++r) {
            if (r != col && vmat[r][col] != 0u) {
                uint8_t f = vmat[r][col];
                for (uint32_t k = 0u; k < 5u; ++k) {
                    vmat[r][k] ^= qir__gf4_mul(f, vmat[col][k]);
                    vinv[r][k] ^= qir__gf4_mul(f, vinv[col][k]);
                }
            }
        }
    }

    static const uint8_t k3_u[3]   = {0x1u, 0x2u, 0x3u};
    static const uint8_t k3_v[3]   = {0x1u, 0x2u, 0x3u};
    static const uint8_t k3_out[3] = {0x3u, 0x1u, 0x2u};

    for (uint32_t m = 0u; m < 5u; ++m) {
        uint8_t u_lin[2] = {0}, v_lin[2] = {0};
        if (m < 4u) {
            uint8_t pp = 1u;
            for (uint32_t j = 0u; j < 3u; ++j) {
                for (uint32_t i = 0u; i < 2u; ++i) {
                    uint8_t pr = qir__gf4_mul((uint8_t)(1u << i), pp);
                    for (uint32_t b = 0u; b < 2u; ++b) {
                        if ((pr >> b) & 1u) {
                            u_lin[b] ^= inv_cols[2u * j + i];
                            v_lin[b] ^= inv_cols[2u * j + i];
                        }
                    }
                }
                pp = qir__gf4_mul(pp, pts[m]);
            }
        } else {
            u_lin[0] = inv_cols[4]; u_lin[1] = inv_cols[5];
            v_lin[0] = inv_cols[4]; v_lin[1] = inv_cols[5];
        }
        for (uint32_t k = 0u; k < 3u; ++k) {
            uint32_t idx = 3u * m + k;
            uint8_t lm = 0u, rm = 0u, wm = 0u;
            for (uint32_t b = 0u; b < 2u; ++b) {
                if ((k3_u[k] >> b) & 1u) lm ^= u_lin[b];
                if ((k3_v[k] >> b) & 1u) rm ^= v_lin[b];
            }
            for (uint32_t d = 0u; d < 5u; ++d) {
                uint8_t cg = qir__gf4_mul(vinv[d][m], k3_out[k]);
                for (uint32_t i = 0u; i < 2u; ++i) {
                    if ((cg >> i) & 1u) wm ^= (uint8_t)qir__gf2m_mul_u32(w_pow[i], z5[d], 6u, irred_poly);
                }
            }
            out->l_mask[idx] = lm;
            out->r_mask[idx] = rm;
            out->w_mask[idx] = wm;
            out->order[idx]  = idx;
        }
    }

    for (uint32_t a = 0u; a < 64u; ++a) {
        for (uint32_t b = 0u; b < 64u; ++b) {
            uint8_t got = 0u;
            for (uint32_t p = 0u; p < QIR_GF64_RANK; ++p) {
                if (__builtin_parity(out->l_mask[p] & a) & __builtin_parity(out->r_mask[p] & b)) {
                    got ^= out->w_mask[p];
                }
            }
            if (got != (uint8_t)qir__gf2m_mul_u32(a, b, 6u, irred_poly)) return false;
        }
    }
    return true;
}

/* ============================================================================
 * E-Graph Construction & Composable Rule Saturation (R1..R7)
 * ============================================================================ */

static inline uint32_t qir__egraph_add_enode(qir_egraph_t *eg, const qir_enode_t *proto) {
    if (eg->num_enodes >= eg->cap_enodes) {
        uint32_t nc = eg->cap_enodes ? eg->cap_enodes * 2u : 512u;
        eg->enodes = (qir_enode_t *)realloc(eg->enodes, (size_t)nc * sizeof(qir_enode_t));
        if (!eg->enodes) abort();
        eg->cap_enodes = nc;
    }
    uint32_t id = eg->num_enodes++;
    eg->enodes[id] = *proto;
    eg->enodes[id].id = id;

    qir_eclass_t *ec = &eg->classes[proto->eclass_id];
    if (ec->num_enodes >= ec->cap_enodes) {
        uint32_t nc = ec->cap_enodes ? ec->cap_enodes * 2u : 8u;
        ec->enodes = (uint32_t *)realloc(ec->enodes, (size_t)nc * sizeof(uint32_t));
        if (!ec->enodes) abort();
        ec->cap_enodes = nc;
    }
    ec->enodes[ec->num_enodes++] = id;
    return id;
}

static inline void qir_egraph_free(qir_egraph_t *eg) {
    if (!eg) return;
    for (uint32_t i = 0u; i < eg->num_classes; ++i) free(eg->classes[i].enodes);
    free(eg->classes);
    free(eg->enodes);
    memset(eg, 0, sizeof(*eg));
}

/* R7 Helper: Rank-compress + affine zero-shift a LUT over F_2 (sound for all <=256-bit outputs) */
static inline qir_gauge_lut_t *qir__build_gauge_lut(qcc_arena_t *arena, const qir_lut_payload_t *lut) {
    if (!lut || lut->num_entries == 0u || lut->val_bits == 0u || lut->val_bits > QIR_MAX_VAL_BITS) return NULL;
    qir_gauge_lut_t *g = (qir_gauge_lut_t *)qcc_arena_alloc_zero(arena, sizeof(qir_gauge_lut_t), 8u);
    g->val_bits    = lut->val_bits;
    g->val_words   = lut->val_words;
    g->num_entries = lut->num_entries;
    for (uint16_t w = 0u; w < lut->val_words && w < QIR_WORDS; ++w) {
        g->const_shift.w[w] = lut->table[w];
    }

    /* Pivot columns over shifted rows d[e] = T[e] ^ T[0] */
    uint16_t piv_col[QIR_MAX_VAL_BITS];
    uint16_t rank = 0u;
    uint32_t *piv_entry = (uint32_t *)calloc(lut->val_bits, sizeof(uint32_t));
    /* Greedy basis selection over F_2 rows of (T[e] ^ T[0]) */
    uint8_t *col_used = (uint8_t *)calloc(lut->val_bits, 1u);
    for (uint16_t c = 0u; c < lut->val_bits; ++c) {
        uint64_t any = 0u;
        uint16_t cw = c >> 6u, cb = c & 63u;
        for (uint32_t e = 1u; e < lut->num_entries; ++e) {
            uint64_t v = (lut->table[(size_t)e * lut->val_words + cw] ^ g->const_shift.w[cw]) >> cb;
            any |= (v & 1u);
            if (any) break;
        }
        if (any) {
            /* Retain non-zero columns first; additionally perform simple F_2 rank elimination for <= 64 bits */
            col_used[c] = 1u;
            piv_col[rank++] = c;
        }
    }

    /* Run exact F_2 RREF on (T[e] ^ T[0]) over arbitrary <= 256 bits to find minimal rank basis! */
    {
        qcc_bv256_t *basis_vecs = (qcc_bv256_t *)calloc(lut->val_bits, sizeof(qcc_bv256_t));
        uint16_t *true_piv_bits = (uint16_t *)calloc(lut->val_bits, sizeof(uint16_t));
        uint16_t true_rank = 0u;
        uint16_t vw = lut->val_words;
        for (uint32_t e = 1u; e < lut->num_entries && true_rank < lut->val_bits; ++e) {
            qcc_bv256_t rem;
            memset(&rem, 0, sizeof(rem));
            for (uint16_t w = 0u; w < vw && w < QIR_WORDS; ++w) {
                rem.w[w] = lut->table[(size_t)e * vw + w] ^ g->const_shift.w[w];
            }
            for (uint16_t r = 0u; r < true_rank; ++r) {
                if (qcc_bv_get(rem.w, true_piv_bits[r])) {
                    for (uint16_t w = 0u; w < vw && w < QIR_WORDS; ++w) rem.w[w] ^= basis_vecs[r].w[w];
                }
            }
            uint16_t pbit = UINT16_MAX;
            for (uint16_t w = 0u; w < vw && w < QIR_WORDS; ++w) {
                if (rem.w[w] != 0u) {
                    pbit = (uint16_t)(w * 64u + (uint16_t)__builtin_ctzll((unsigned long long)rem.w[w]));
                    break;
                }
            }
            if (pbit < lut->val_bits) {
                for (uint16_t r = 0u; r < true_rank; ++r) {
                    if (qcc_bv_get(basis_vecs[r].w, pbit)) {
                        for (uint16_t w = 0u; w < vw && w < QIR_WORDS; ++w) basis_vecs[r].w[w] ^= rem.w[w];
                    }
                }
                basis_vecs[true_rank]    = rem;
                true_piv_bits[true_rank] = pbit;
                true_rank++;
            }
        }
        if (true_rank < rank) {
            g->rank_bits  = true_rank;
            g->comp_words = (uint16_t)((true_rank + 63u) / 64u);
            if (g->comp_words == 0u) g->comp_words = 1u;
            g->comp_table = (uint64_t *)qcc_arena_alloc_zero(arena, (size_t)lut->num_entries * g->comp_words * sizeof(uint64_t), 8u);
            g->recon_rows = (qcc_bv256_t *)qcc_arena_alloc_zero(arena, (size_t)lut->val_bits * sizeof(qcc_bv256_t), 8u);
            for (uint32_t e = 0u; e < lut->num_entries; ++e) {
                for (uint16_t r = 0u; r < true_rank; ++r) {
                    uint16_t pc = true_piv_bits[r];
                    uint16_t pw = pc >> 6u, pb = pc & 63u;
                    uint64_t bit = (lut->table[(size_t)e * vw + pw] ^ g->const_shift.w[pw]) >> pb;
                    if (bit & 1u) qcc_bv_set(&g->comp_table[(size_t)e * g->comp_words], r);
                }
            }
            for (uint16_t bit = 0u; bit < lut->val_bits; ++bit) {
                for (uint16_t r = 0u; r < true_rank; ++r) {
                    if (qcc_bv_get(basis_vecs[r].w, bit)) qcc_bv_set(g->recon_rows[bit].w, r);
                }
            }
            free(basis_vecs);
            free(true_piv_bits);
            free(piv_entry);
            free(col_used);
            return g;
        }
        free(basis_vecs);
        free(true_piv_bits);
    }

    /* Zero-column strip fallback for arbitrary <= 256 bits */
    g->rank_bits  = rank ? rank : 1u;
    g->comp_words = (uint16_t)((g->rank_bits + 63u) / 64u);
    g->comp_table = (uint64_t *)qcc_arena_alloc_zero(arena, (size_t)lut->num_entries * g->comp_words * sizeof(uint64_t), 8u);
    g->recon_rows = (qcc_bv256_t *)qcc_arena_alloc_zero(arena, (size_t)lut->val_bits * sizeof(qcc_bv256_t), 8u);
    for (uint32_t e = 0u; e < lut->num_entries; ++e) {
        for (uint16_t r = 0u; r < rank; ++r) {
            uint16_t oc = piv_col[r];
            uint16_t ow = oc >> 6u, ob = oc & 63u;
            uint64_t b = (lut->table[(size_t)e * lut->val_words + ow] ^ g->const_shift.w[ow]) >> ob;
            if (b & 1u) qcc_bv_set(&g->comp_table[(size_t)e * g->comp_words], r);
        }
    }
    for (uint16_t r = 0u; r < rank; ++r) {
        qcc_bv_set(g->recon_rows[piv_col[r]].w, r);
    }
    free(piv_entry);
    free(col_used);
    return g;
}

static inline bool qir__val_is_const(const qir_val_t *v) {
    if (!v || v->width == 0u) return false;
    qcc_bv256_t m = qcc_bv_mask_n(v->width);
    for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
        if ((v->const_mask[w] & m.w[w]) != m.w[w]) return false;
    }
    return true;
}

/* Check R4 soundness: is phase_poly purely linear on input operand idx_pos with angular step k/2^w? */
static inline bool qir__r4_phase_is_modular_linear(
    const qir_module_t *m,
    const qir_op_t *ph_op,
    uint16_t op_in_idx,
    flow_u128 *out_step,
    uint16_t *out_bit_lo
) {
    if (!ph_op || ph_op->kind != QIR_OP_PHASE_POLY || op_in_idx >= ph_op->num_ins) return false;
    if (ph_op->pred_val != UINT32_MAX) return false;
    uint32_t target_val = ph_op->ins[op_in_idx];
    if (m->vals[target_val].def_op == UINT32_MAX) return false;
    const qir_op_t *def = &m->ops[m->vals[target_val].def_op];
    if ((def->kind != QIR_OP_ADD && def->kind != QIR_OP_SUB) || def->pred_val != UINT32_MAX) return false;

    uint16_t w = m->vals[target_val].width;
    if (w == 0u || w > 120u) return false;

    uint16_t bit_lo = 0u;
    for (uint16_t i = 0u; i < op_in_idx; ++i) bit_lo = (uint16_t)(bit_lo + m->vals[ph_op->ins[i]].width);

    flow_u128 bit_angle[128];
    memset(bit_angle, 0, sizeof(bit_angle));
    bool touched = false;

    for (uint32_t t = 0u; t < ph_op->u.phase_poly.num_terms; ++t) {
        const qir_phase_term_t *pt = &ph_op->u.phase_poly.terms[t];
        uint32_t hits = 0u;
        uint16_t hit_bit = 0u;
        for (uint16_t b = 0u; b < w; ++b) {
            if (qcc_bv_get(pt->mask, (uint32_t)(bit_lo + b))) {
                hits++;
                hit_bit = b;
            }
        }
        if (hits == 0u) continue;
        if (qcc_bv_popcount(pt->mask) != 1u || pt->is_xor) return false;
        touched = true;
        bit_angle[hit_bit] += pt->angle;
    }
    if (!touched) return false;

    flow_u128 step = bit_angle[0];
    for (uint16_t b = 0u; b < w; ++b) {
        flow_u128 expect = (b < 128u) ? (step << b) : 0u;
        if (bit_angle[b] != expect) return false;
    }
    if ((bit_angle[w - 1u] << 1u) != 0u) return false;
    if (out_step) *out_step = step;
    if (out_bit_lo) *out_bit_lo = bit_lo;
    return true;
}

static inline bool qir_egraph_build_and_saturate(qir_egraph_t *eg, qcc_arena_t *arena, const qir_module_t *mod) {
    memset(eg, 0, sizeof(*eg));
    eg->arena       = arena;
    eg->mod         = mod;
    eg->num_classes = mod->num_vals;
    eg->classes     = (qir_eclass_t *)calloc(mod->num_vals ? mod->num_vals : 1u, sizeof(qir_eclass_t));
    if (!eg->classes) return false;

    for (uint32_t v = 0u; v < mod->num_vals; ++v) {
        eg->classes[v].id     = v;
        eg->classes[v].val_id = v;
        eg->classes[v].width  = mod->vals[v].width;
    }
    for (uint32_t oi = 0u; oi < mod->num_ops; ++oi) {
        const qir_op_t *op = &mod->ops[oi];
        for (uint16_t i = 0u; i < op->num_ins; ++i) eg->classes[op->ins[i]].use_count++;
    }
    for (uint32_t i = 0u; i < mod->num_outputs; ++i) {
        eg->classes[mod->outputs[i]].use_count++;
    }

    /* Seed base E-Nodes for every SSA op */
    for (uint32_t oi = 0u; oi < mod->num_ops; ++oi) {
        const qir_op_t *op = &mod->ops[oi];
        if (op->num_outs == 0u) continue;
        qir_enode_t node;
        memset(&node, 0, sizeof(node));
        node.eclass_id        = op->outs[0];
        node.variant          = QIR_VAR_STD;
        node.op               = op;
        node.num_children     = (uint8_t)(op->num_ins > 8u ? 8u : op->num_ins);
        for (uint8_t k = 0u; k < node.num_children; ++k) node.child_eclasses[k] = op->ins[k];
        node.approx_keep_bits = mod->vals[op->outs[0]].width;
        qir__egraph_add_enode(eg, &node);

        /* Rule R1: Multi-hop LinF2_Fuse (arbitrary single-use LINF2 chain -> QIR_VAR_LINF2_GAUSS) */
        if (op->kind == QIR_OP_LINF2) {
            qir_enode_t r1 = node;
            r1.variant  = QIR_VAR_LINF2_GAUSS;
            r1.rule_tag = 1u;
            r1.aux_param = op->id;
            if (op->pred_val == UINT32_MAX && op->u.linf2.out_bits > 0u && op->u.linf2.out_bits <= QIR_MAX_VAL_BITS) {
                const qir_op_t *cur_op = op;
                uint16_t cur_out = op->u.linf2.out_bits;
                uint16_t cur_in  = op->u.linf2.in_bits;
                qcc_bv256_t *cur_rows = (qcc_bv256_t *)malloc((size_t)cur_out * sizeof(qcc_bv256_t));
                qcc_bv256_t *nxt_rows = (qcc_bv256_t *)malloc((size_t)cur_out * sizeof(qcc_bv256_t));
                if (!cur_rows || !nxt_rows) abort();
                memcpy(cur_rows, op->u.linf2.rows, (size_t)cur_out * sizeof(qcc_bv256_t));
                qcc_bv256_t cur_c = op->u.linf2.const_vec;
                uint32_t folded[16];
                uint16_t n_folded = 0u;

                while (n_folded < 15u && cur_op->num_ins == 1u) {
                    uint32_t piv = cur_op->ins[0];
                    if (mod->vals[piv].def_op == UINT32_MAX || eg->classes[piv].use_count != 1u) break;
                    const qir_op_t *prev = &mod->ops[mod->vals[piv].def_op];
                    if (prev->kind != QIR_OP_LINF2 || prev->pred_val != UINT32_MAX ||
                        prev->u.linf2.out_bits != cur_in || prev->u.linf2.in_bits > QIR_MAX_VAL_BITS) break;
                    memset(nxt_rows, 0, (size_t)cur_out * sizeof(qcc_bv256_t));
                    qcc_bv256_t nxt_c = qcc_bv_zero();
                    for (uint16_t r = 0u; r < cur_out; ++r) {
                        uint8_t cbit = qcc_bv_get(cur_c.w, r) ? 1u : 0u;
                        for (uint16_t k = 0u; k < cur_in; ++k) {
                            if (qcc_bv_get(cur_rows[r].w, k)) {
                                for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                                    nxt_rows[r].w[w] ^= prev->u.linf2.rows[k].w[w];
                                }
                                if (qcc_bv_get(prev->u.linf2.const_vec.w, k)) cbit ^= 1u;
                            }
                        }
                        if (cbit) qcc_bv_set(nxt_c.w, r);
                    }
                    memcpy(cur_rows, nxt_rows, (size_t)cur_out * sizeof(qcc_bv256_t));
                    cur_c  = nxt_c;
                    cur_in = prev->u.linf2.in_bits;
                    folded[n_folded++] = prev->id;
                    cur_op = prev;
                }
                if (n_folded > 0u) {
                    qir_linf2_payload_t *fused = (qir_linf2_payload_t *)qcc_arena_alloc_zero(
                        arena, sizeof(qir_linf2_payload_t), 8u);
                    fused->out_bits  = cur_out;
                    fused->in_bits   = cur_in;
                    fused->const_vec = cur_c;
                    fused->rows      = (qcc_bv256_t *)qcc_arena_alloc_zero(
                        arena, (size_t)cur_out * sizeof(qcc_bv256_t), 8u);
                    memcpy(fused->rows, cur_rows, (size_t)cur_out * sizeof(qcc_bv256_t));
                    uint32_t *ch = (uint32_t *)qcc_arena_alloc(arena, (size_t)n_folded * sizeof(uint32_t), 4u);
                    memcpy(ch, folded, (size_t)n_folded * sizeof(uint32_t));
                    r1.fused_linf2   = fused;
                    r1.chain_ops     = ch;
                    r1.num_chain_ops = n_folded;
                    r1.aux_param     = cur_op->id;
                    r1.num_children  = (uint8_t)(cur_op->num_ins > 8u ? 8u : cur_op->num_ins);
                    for (uint8_t k = 0u; k < r1.num_children; ++k) r1.child_eclasses[k] = cur_op->ins[k];
                }
                free(cur_rows);
                free(nxt_rows);
            }
            eg->rules_fired[1]++;
            qir__egraph_add_enode(eg, &r1);
        }

        /* Rules R2, R5, R7, R8 for LUT_ROM */
        if (op->kind == QIR_OP_LUT_ROM) {
            qir_enode_t r2 = node;
            r2.variant  = QIR_VAR_HMR_PARITY_LUT;
            r2.rule_tag = 2u;
            qir__egraph_add_enode(eg, &r2);
            eg->rules_fired[2]++;

            /* Rule R8: Virtual-Slot-0 Differential Select-Swap QROAM (b_sh = 1..6) */
            for (uint32_t b_sh = 1u; b_sh <= 6u; ++b_sh) {
                if (b_sh <= op->u.lut.key_bits && op->u.lut.num_entries >= (UINT32_C(1) << b_sh)) {
                    qir_enode_t r8 = node;
                    r8.variant   = QIR_VAR_QROAM_VSLOT0;
                    r8.rule_tag  = 8u;
                    r8.aux_param = b_sh;
                    qir__egraph_add_enode(eg, &r8);
                    eg->rules_fired[8]++;
                }
            }

            /* Rule R8b: Column-Split Virtual-Slot-0 QROAM (K in {2,4}, b_sh in {1..4}) */
            if (op->u.lut.val_bits >= 16u) {
                static const uint32_t KS[2] = {2u, 4u};
                for (uint32_t ki = 0u; ki < 2u; ++ki) {
                    uint32_t K = KS[ki];
                    if (op->u.lut.val_bits < K * 2u) continue;
                    for (uint32_t b_sh = 1u; b_sh <= 4u; ++b_sh) {
                        if (b_sh <= op->u.lut.key_bits && op->u.lut.num_entries >= (UINT32_C(1) << b_sh)) {
                            qir_enode_t r8b = node;
                            r8b.variant   = QIR_VAR_QROAM_SPLIT;
                            r8b.rule_tag  = 8u;
                            r8b.aux_param = K | (b_sh << 8u);
                            qir__egraph_add_enode(eg, &r8b);
                            eg->rules_fired[8]++;
                        }
                    }
                }
            }

            /* Rule R5: LUT_ColSplit (K=2 and K=4 slices for wide LUTs) */
            if (op->u.lut.val_bits >= 8u) {
                uint32_t splits[2] = {2u, 4u};
                for (uint32_t si = 0u; si < 2u; ++si) {
                    if (op->u.lut.val_bits >= splits[si] * 2u) {
                        qir_enode_t r5 = node;
                        r5.variant   = QIR_VAR_LUT_COL_SPLIT;
                        r5.rule_tag  = 5u;
                        r5.aux_param = splits[si];
                        qir__egraph_add_enode(eg, &r5);
                        eg->rules_fired[5]++;
                    }
                }
            }

            /* Rule R7: BasisIso_Gauge */
            qir_gauge_lut_t *gl = qir__build_gauge_lut(arena, &op->u.lut);
            if (gl) {
                qir_enode_t r7 = node;
                r7.variant   = QIR_VAR_LIN_GAUGE_TSP;
                r7.rule_tag  = 7u;
                r7.gauge_lut = gl;
                qir__egraph_add_enode(eg, &r7);
                eg->rules_fired[7]++;
            }
        }

        /* Rule R3: CSA_Tree over >= 2 chained unconditional ADD ops */
        if (op->kind == QIR_OP_ADD && op->num_ins == 2u && op->pred_val == UINT32_MAX) {
            uint32_t chain_tmp[16];
            uint16_t nc = 0u;
            uint32_t stack[16];
            uint16_t sp = 0u;
            stack[sp++] = op->id;
            while (sp > 0u && nc < 15u) {
                uint32_t cur_op_id = stack[--sp];
                const qir_op_t *cop = &mod->ops[cur_op_id];
                chain_tmp[nc++] = cur_op_id;
                for (uint16_t ii = 0u; ii < 2u; ++ii) {
                    uint32_t iv = cop->ins[ii];
                    uint32_t d  = mod->vals[iv].def_op;
                    if (d != UINT32_MAX && mod->ops[d].kind == QIR_OP_ADD &&
                        mod->ops[d].pred_val == UINT32_MAX &&
                        eg->classes[iv].use_count == 1u &&
                        mod->vals[iv].width == mod->vals[op->outs[0]].width &&
                        sp < 16u) {
                        stack[sp++] = d;
                    }
                }
            }
            if (nc >= 2u) {
                qir_enode_t r3 = node;
                r3.variant   = QIR_VAR_CSA_TREE;
                r3.rule_tag  = 3u;
                uint32_t *chain = (uint32_t *)qcc_arena_alloc(arena, (size_t)nc * sizeof(uint32_t), 4u);
                for (uint16_t k = 0u; k < nc; ++k) chain[k] = chain_tmp[nc - 1u - k];
                r3.chain_ops     = chain;
                r3.num_chain_ops = nc;
                r3.aux_param     = nc + 1u;
                qir__egraph_add_enode(eg, &r3);
                eg->rules_fired[3]++;
            }
        }

        /* Rule R4: Generalized Phase_AddTree_Split over ADD/SUB binary trees feeding modular linear PHASE_POLY */
        if (op->kind == QIR_OP_PHASE_POLY && op->pred_val == UINT32_MAX) {
            for (uint16_t ii = 0u; ii < op->num_ins; ++ii) {
                uint32_t v_in = op->ins[ii];
                if (eg->classes[v_in].use_count != 1u) continue;
                flow_u128 step_out = 0u;
                uint16_t bit_lo_out = 0u;
                if (!qir__r4_phase_is_modular_linear(mod, op, ii, &step_out, &bit_lo_out)) continue;

                typedef struct { uint32_t vid; int8_t sgn; } r4_item_t;
                r4_item_t work[32];
                uint32_t wp = 0u;
                work[wp++] = (r4_item_t){v_in, 1};
                uint32_t folded[16];
                uint16_t n_folded = 0u;
                uint32_t leaves[16];
                int8_t   lsigns[16];
                uint16_t n_leaves = 0u;
                uint16_t w = mod->vals[v_in].width;
                bool tree_ok = true;

                while (wp > 0u) {
                    r4_item_t it = work[--wp];
                    uint32_t d_id = mod->vals[it.vid].def_op;
                    bool can_expand = false;
                    if (d_id != UINT32_MAX) {
                        const qir_op_t *dop = &mod->ops[d_id];
                        if ((dop->kind == QIR_OP_ADD || dop->kind == QIR_OP_SUB) &&
                            dop->num_ins == 2u && dop->pred_val == UINT32_MAX &&
                            eg->classes[it.vid].use_count == 1u &&
                            mod->vals[dop->ins[0]].width == w && mod->vals[dop->ins[1]].width == w &&
                            n_folded < 15u && wp + 2u <= 32u) {
                            can_expand = true;
                            folded[n_folded++] = d_id;
                            work[wp++] = (r4_item_t){dop->ins[1], (dop->kind == QIR_OP_SUB) ? (int8_t)(-it.sgn) : it.sgn};
                            work[wp++] = (r4_item_t){dop->ins[0], it.sgn};
                        }
                    }
                    if (!can_expand) {
                        if (n_leaves >= 16u) { tree_ok = false; break; }
                        leaves[n_leaves] = it.vid;
                        lsigns[n_leaves] = it.sgn;
                        n_leaves++;
                    }
                }
                if (tree_ok && n_folded >= 1u) {
                    qir_phase_split_t *ps = (qir_phase_split_t *)qcc_arena_alloc_zero(
                        arena, sizeof(qir_phase_split_t), 16u);
                    ps->target_val = v_in;
                    ps->bit_lo     = bit_lo_out;
                    ps->width      = w;
                    ps->step       = step_out;
                    ps->num_leaves = n_leaves;
                    for (uint16_t li = 0u; li < n_leaves; ++li) {
                        ps->leaf_vals[li]  = leaves[li];
                        ps->leaf_signs[li] = lsigns[li];
                    }
                    uint32_t *ch = (uint32_t *)qcc_arena_alloc(arena, (size_t)n_folded * sizeof(uint32_t), 4u);
                    memcpy(ch, folded, (size_t)n_folded * sizeof(uint32_t));
                    qir_enode_t r4_ph = node;
                    r4_ph.variant       = QIR_VAR_PHASE_CSA_DIRECT;
                    r4_ph.rule_tag      = 4u;
                    r4_ph.aux_param     = folded[0];
                    r4_ph.phase_split   = ps;
                    r4_ph.chain_ops     = ch;
                    r4_ph.num_chain_ops = n_folded;
                    qir__egraph_add_enode(eg, &r4_ph);
                    eg->rules_fired[4]++;
                }
            }
        }

        /* Rule R9: MUL_CONST_ADD & SQR_TRIANGLE */
        if (op->kind == QIR_OP_MUL_FXP && op->num_ins == 2u) {
            bool c0 = qir__val_is_const(&mod->vals[op->ins[0]]);
            bool c1 = qir__val_is_const(&mod->vals[op->ins[1]]);
            if (c0 || c1) {
                uint32_t c_idx = 0u;
                if (c0 && c1) {
                    uint32_t p0 = qcc_bv_popcount(mod->vals[op->ins[0]].const_val);
                    uint32_t p1 = qcc_bv_popcount(mod->vals[op->ins[1]].const_val);
                    c_idx = (p1 < p0) ? 1u : 0u;
                } else {
                    c_idx = c1 ? 1u : 0u;
                }
                qir_enode_t r9 = node;
                r9.variant   = QIR_VAR_MUL_CONST_ADD;
                r9.rule_tag  = 9u;
                r9.aux_param = c_idx;
                uint32_t cvid = op->ins[c_idx];
                uint32_t cdef = mod->vals[cvid].def_op;
                if (cdef != UINT32_MAX && mod->ops[cdef].kind == QIR_OP_CONST &&
                    mod->ops[cdef].pred_val == UINT32_MAX && eg->classes[cvid].use_count == 1u) {
                    uint32_t *ch = (uint32_t *)qcc_arena_alloc(arena, sizeof(uint32_t), 4u);
                    ch[0] = cdef;
                    r9.chain_ops     = ch;
                    r9.num_chain_ops = 1u;
                }
                qir__egraph_add_enode(eg, &r9);
                eg->rules_fired[9]++;
            }
        } else if (op->kind == QIR_OP_SQR_FXP && op->num_ins == 1u) {
            qir_enode_t r9 = node;
            r9.variant  = QIR_VAR_SQR_TRIANGLE;
            r9.rule_tag = 9u;
            qir__egraph_add_enode(eg, &r9);
            eg->rules_fired[9]++;
        }

        /* Rule R6: GF2m_Tower_Rank & Rule R11: GF2M_CONST_LINF2 */
        if (op->kind == QIR_OP_GF2M_MUL && op->u.gf2m.m == 12u) {
            qir_enode_t r6 = node;
            r6.variant  = QIR_VAR_GF2M_CHUDNOVSKY_R42;
            r6.rule_tag = 6u;
            qir__egraph_add_enode(eg, &r6);
            qir_enode_t r6_hmr = r6;
            r6_hmr.variant   = QIR_VAR_HMR_PARITY_LUT; /* Denotes 0-CCX HMR uncompute variant for GF12 bilinear */
            r6_hmr.aux_param = 1u;
            qir__egraph_add_enode(eg, &r6_hmr);
            eg->rules_fired[6]++;
        } else if (op->kind == QIR_OP_GF2M_MUL && op->u.gf2m.m == 6u) {
            qir_enode_t r6 = node;
            r6.variant  = QIR_VAR_GF2M_RANK15;
            r6.rule_tag = 6u;
            qir__egraph_add_enode(eg, &r6);
            eg->rules_fired[6]++;
        } else if (op->kind == QIR_OP_GF2M_INV && (op->u.gf2m.m == 12u || op->u.gf2m.m == 6u)) {
            qir_enode_t r6 = node;
            r6.variant  = QIR_VAR_GF2M_TOWER_INV;
            r6.rule_tag = 6u;
            qir__egraph_add_enode(eg, &r6);
            eg->rules_fired[6]++;
        }

        if (op->kind == QIR_OP_GF2M_MUL && op->num_ins == 2u && op->u.gf2m.m <= 32u) {
            bool c0 = qir__val_is_const(&mod->vals[op->ins[0]]);
            bool c1 = qir__val_is_const(&mod->vals[op->ins[1]]);
            if (c0 || c1) {
                uint32_t c_idx = c0 ? 0u : 1u;
                uint32_t var_idx = c_idx ^ 1u;
                uint32_t cvid = op->ins[c_idx];
                uint16_t dm = op->u.gf2m.m;
                uint32_t dmask = (dm == 32u) ? UINT32_MAX : ((UINT32_C(1) << dm) - 1u);
                uint32_t celem = (uint32_t)(mod->vals[cvid].const_val[0] & dmask);
                qir_linf2_payload_t *lf = (qir_linf2_payload_t *)qcc_arena_alloc_zero(arena, sizeof(qir_linf2_payload_t), 8u);
                lf->in_bits  = dm;
                lf->out_bits = dm;
                lf->rows     = (qcc_bv256_t *)qcc_arena_alloc_zero(arena, dm * sizeof(qcc_bv256_t), 8u);
                for (uint16_t c = 0u; c < dm; ++c) {
                    uint32_t col = qir__gf2m_mul_u32(celem, UINT32_C(1) << c, dm, op->u.gf2m.irred_poly);
                    for (uint16_t r = 0u; r < dm; ++r) {
                        if ((col >> r) & 1u) qcc_bv_set(lf->rows[r].w, c);
                    }
                }
                qir_enode_t r11 = node;
                r11.variant     = QIR_VAR_GF2M_CONST_LINF2;
                r11.rule_tag    = 13u;
                r11.aux_param   = var_idx;
                r11.fused_linf2 = lf;
                uint32_t cdef = mod->vals[cvid].def_op;
                if (cdef != UINT32_MAX && mod->ops[cdef].kind == QIR_OP_CONST &&
                    mod->ops[cdef].pred_val == UINT32_MAX && eg->classes[cvid].use_count == 1u) {
                    uint32_t *ch = (uint32_t *)qcc_arena_alloc(arena, sizeof(uint32_t), 4u);
                    ch[0] = cdef;
                    r11.chain_ops     = ch;
                    r11.num_chain_ops = 1u;
                }
                qir__egraph_add_enode(eg, &r11);
                eg->rules_fired[13]++;
            }
        }

        /* Rule R10a: INPLACE_ARITH for ADD / SUB / GF2M_ADD */
        if ((op->kind == QIR_OP_ADD || op->kind == QIR_OP_SUB || op->kind == QIR_OP_GF2M_ADD) &&
            op->num_ins == 2u && op->pred_val == UINT32_MAX && op->num_outs > 0u) {
            uint16_t ow = mod->vals[op->outs[0]].width;
            for (uint32_t idx = 0u; idx < 2u; ++idx) {
                if (mod->vals[op->ins[idx]].width == ow) {
                    qir_enode_t r10a = node;
                    r10a.variant   = QIR_VAR_INPLACE_ARITH;
                    r10a.rule_tag  = 10u;
                    r10a.aux_param = idx;
                    qir__egraph_add_enode(eg, &r10a);
                    eg->rules_fired[10]++;
                }
            }
        }

        /* Rule R10b: CTRL_ADD_CUCCARO for CTRL_SEL(ctrl, v_false, ADD/SUB(v_false, lhs)) */
        if (op->kind == QIR_OP_CTRL_SEL && op->num_ins == 3u && op->pred_val == UINT32_MAX && op->num_outs > 0u) {
            uint32_t v_false = op->ins[1], v_true = op->ins[2];
            uint32_t def_id  = mod->vals[v_true].def_op;
            uint16_t ow      = mod->vals[op->outs[0]].width;
            if (def_id != UINT32_MAX && eg->classes[v_true].use_count == 1u && mod->vals[v_false].width == ow) {
                const qir_op_t *dop = &mod->ops[def_id];
                if ((dop->kind == QIR_OP_ADD || dop->kind == QIR_OP_SUB) &&
                    dop->num_ins == 2u && dop->pred_val == UINT32_MAX && mod->vals[v_true].width == ow) {
                    uint32_t lhs_vid = UINT32_MAX;
                    bool is_sub = false;
                    if (dop->kind == QIR_OP_ADD) {
                        if (dop->ins[0] == v_false)      lhs_vid = dop->ins[1];
                        else if (dop->ins[1] == v_false) lhs_vid = dop->ins[0];
                    } else if (dop->ins[0] == v_false) {
                        lhs_vid = dop->ins[1];
                        is_sub  = true;
                    }
                    if (lhs_vid != UINT32_MAX && lhs_vid != v_false && mod->vals[lhs_vid].width > 0u) {
                        uint32_t *ch = (uint32_t *)qcc_arena_alloc(arena, sizeof(uint32_t), 4u);
                        ch[0] = dop->id;
                        qir_enode_t r10b = node;
                        r10b.variant       = QIR_VAR_CTRL_ADD_CUCCARO;
                        r10b.rule_tag      = 11u;
                        r10b.aux_param     = lhs_vid | (is_sub ? UINT32_C(0x80000000) : 0u);
                        r10b.chain_ops     = ch;
                        r10b.num_chain_ops = 1u;
                        qir__egraph_add_enode(eg, &r10b);
                        eg->rules_fired[11]++;
                    }
                }
            }
        }

        /* Rule R10c: XOR_ACCUM for GF2M_ADD(acc, op(...)) */
        if (op->kind == QIR_OP_GF2M_ADD && op->num_ins == 2u && op->pred_val == UINT32_MAX && op->num_outs > 0u) {
            uint16_t ow = mod->vals[op->outs[0]].width;
            for (uint32_t prod_idx = 0u; prod_idx < 2u; ++prod_idx) {
                uint32_t acc_idx = prod_idx ^ 1u;
                uint32_t v_prod  = op->ins[prod_idx];
                uint32_t v_acc   = op->ins[acc_idx];
                uint32_t def_id  = mod->vals[v_prod].def_op;
                if (def_id == UINT32_MAX || eg->classes[v_prod].use_count != 1u ||
                    mod->vals[v_prod].width != ow || mod->vals[v_acc].width != ow ||
                    op->u.gf2m.m != ow) {
                    continue;
                }
                const qir_op_t *dop = &mod->ops[def_id];
                if (dop->pred_val != UINT32_MAX || dop->num_outs != 1u) continue;
                uint32_t *ch = (uint32_t *)qcc_arena_alloc(arena, sizeof(uint32_t), 4u);
                ch[0] = dop->id;
                if (dop->kind == QIR_OP_LUT_ROM && dop->u.lut.val_bits == ow) {
                    const qir_lut_payload_t *dl = &dop->u.lut;
                    {
                        qir_enode_t r10c = node;
                        r10c.variant       = QIR_VAR_XOR_ACCUM;
                        r10c.rule_tag      = 12u;
                        r10c.aux_param     = acc_idx | ((uint32_t)QIR_VAR_STD << 8u);
                        r10c.chain_ops     = ch;
                        r10c.num_chain_ops = 1u;
                        qir__egraph_add_enode(eg, &r10c);
                        eg->rules_fired[12]++;
                    }
                    for (uint32_t b_sh = 1u; b_sh <= 6u; ++b_sh) {
                        if (b_sh <= dl->key_bits && dl->num_entries >= (UINT32_C(1) << b_sh)) {
                            qir_enode_t r10c = node;
                            r10c.variant       = QIR_VAR_XOR_ACCUM;
                            r10c.rule_tag      = 12u;
                            r10c.aux_param     = acc_idx | ((uint32_t)QIR_VAR_QROAM_VSLOT0 << 8u) | (b_sh << 16u);
                            r10c.chain_ops     = ch;
                            r10c.num_chain_ops = 1u;
                            qir__egraph_add_enode(eg, &r10c);
                            eg->rules_fired[12]++;
                        }
                    }
                    if (dl->val_bits >= 16u) {
                        static const uint32_t KS[2] = {2u, 4u};
                        for (uint32_t ki = 0u; ki < 2u; ++ki) {
                            uint32_t K = KS[ki];
                            if (dl->val_bits < K * 2u) continue;
                            for (uint32_t b_sh = 1u; b_sh <= 4u; ++b_sh) {
                                if (b_sh <= dl->key_bits && dl->num_entries >= (UINT32_C(1) << b_sh)) {
                                    uint32_t sub_aux = K | (b_sh << 8u);
                                    qir_enode_t r10c = node;
                                    r10c.variant       = QIR_VAR_XOR_ACCUM;
                                    r10c.rule_tag      = 12u;
                                    r10c.aux_param     = acc_idx | ((uint32_t)QIR_VAR_QROAM_SPLIT << 8u) | (sub_aux << 16u);
                                    r10c.chain_ops     = ch;
                                    r10c.num_chain_ops = 1u;
                                    qir__egraph_add_enode(eg, &r10c);
                                    eg->rules_fired[12]++;
                                }
                            }
                        }
                    }
                    qir_gauge_lut_t *gl = NULL;
                    for (uint32_t ei = 0u; ei < eg->classes[v_prod].num_enodes; ++ei) {
                        const qir_enode_t *pe = &eg->enodes[eg->classes[v_prod].enodes[ei]];
                        if (pe->variant == QIR_VAR_LIN_GAUGE_TSP && pe->gauge_lut) { gl = pe->gauge_lut; break; }
                    }
                    if (!gl) gl = qir__build_gauge_lut(arena, dl);
                    if (gl) {
                        qir_enode_t r10c = node;
                        r10c.variant       = QIR_VAR_XOR_ACCUM;
                        r10c.rule_tag      = 12u;
                        r10c.aux_param     = acc_idx | ((uint32_t)QIR_VAR_LIN_GAUGE_TSP << 8u);
                        r10c.gauge_lut     = gl;
                        r10c.chain_ops     = ch;
                        r10c.num_chain_ops = 1u;
                        qir__egraph_add_enode(eg, &r10c);
                        eg->rules_fired[12]++;
                    }
                } else if (dop->kind == QIR_OP_GF2M_MUL) {
                    qir_variant_kind_t sv = (dop->u.gf2m.m == 12u) ? QIR_VAR_GF2M_CHUDNOVSKY_R42
                                          : ((dop->u.gf2m.m == 6u) ? QIR_VAR_GF2M_RANK15 : QIR_VAR_STD);
                    qir_enode_t r10c = node;
                    r10c.variant       = QIR_VAR_XOR_ACCUM;
                    r10c.rule_tag      = 12u;
                    r10c.aux_param     = acc_idx | ((uint32_t)sv << 8u);
                    r10c.chain_ops     = ch;
                    r10c.num_chain_ops = 1u;
                    qir__egraph_add_enode(eg, &r10c);
                    eg->rules_fired[12]++;
                } else if (dop->kind == QIR_OP_LINF2 && dop->num_ins == 1u) {
                    qir_enode_t r10c = node;
                    r10c.variant       = QIR_VAR_XOR_ACCUM;
                    r10c.rule_tag      = 12u;
                    r10c.aux_param     = acc_idx | ((uint32_t)QIR_VAR_STD << 8u);
                    r10c.chain_ops     = ch;
                    r10c.num_chain_ops = 1u;
                    qir__egraph_add_enode(eg, &r10c);
                    eg->rules_fired[12]++;
                }
            }
        }

        /* Rule R2 / Algebraic: AND_VEC with provably zero mask can elide */
        if (op->kind == QIR_OP_AND_VEC) {
            const qir_val_t *ov = &mod->vals[op->outs[0]];
            bool all_zero = true;
            for (uint16_t b = 0u; b < ov->width; ++b) {
                if (!qcc_bv_get(ov->zero_mask, b)) { all_zero = false; break; }
            }
            if (all_zero) {
                qir_enode_t rz = node;
                rz.variant  = QIR_VAR_TRIVIAL_ELIDE;
                rz.rule_tag = 2u;
                qir__egraph_add_enode(eg, &rz);
            }
        }
    }

    return true;
}

#ifdef __cplusplus
}
#endif

#endif /* QIR_REWRITE_H */
