/*
 * qir_opt.h — Symbiotic SEARU x FLOW 2-Jet QIR Optimizer & Bit-Level Clean/Dirty Borrow Allocator
 * Part of QCC (Quantum Circuit Compiler) for FLOW v7.0.1.
 * Pure C17, zero external dependencies, zero indirect function calls.
 */
#ifndef QIR_OPT_H
#define QIR_OPT_H

#include "qir_codegen.h"
#include "searu.h"
#include "flow_opt_params.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t        max_qubits;       /* Hard qubit budget (0 = unconstrained) */
    qcc_cost_mode_t cost_mode;        /* Cost metric (DEFAULT, QC, Q2T, QT) */
    uint32_t        searu_iters;      /* Max SEARU iterations (0 => default 256) */
    uint32_t        flow_iters;       /* Max FLOW 2-Jet relaxation steps (0 => default 16) */
    uint64_t        seed;
    uint32_t        enabled_rules;    /* Bitmask bit r (1..15) => Rule Rr enabled (0xFFFF = all) */
    bool            aggressive_uncompute;
    bool            require_clean_scratch;
    bool            emit_clean_reg;   /* Declare non-IO clean qubits under "@clean" register */
} qcc_opt_cfg_t;

static inline qcc_opt_cfg_t qcc_opt_cfg_default(void) {
    qcc_opt_cfg_t cfg;
    cfg.max_qubits            = 0u;
    cfg.cost_mode             = QCC_COST_DEFAULT;
    cfg.searu_iters           = 256u;
    cfg.flow_iters            = 16u;
    cfg.seed                  = UINT64_C(0x51C0FFEE77771234);
    cfg.enabled_rules         = 0xFFFFu;
    cfg.aggressive_uncompute  = true;
    cfg.require_clean_scratch = false;
    cfg.emit_clean_reg        = false;
    return cfg;
}

static inline bool qir__op_supports_uncompute(const qir_op_t *op, bool inplace) {
    if (!op) return false;
    switch (op->kind) {
        case QIR_OP_CONST:
        case QIR_OP_SLICE:
        case QIR_OP_CONCAT:
        case QIR_OP_AND_VEC:
        case QIR_OP_ADD:
        case QIR_OP_SUB:
        case QIR_OP_CS_ADD3:
        case QIR_OP_MUL_FXP:
        case QIR_OP_SQR_FXP:
        case QIR_OP_GF2M_ADD:
        case QIR_OP_GF2M_MUL:
        case QIR_OP_GF2M_INV:
        case QIR_OP_LUT_ROM:
        case QIR_OP_CLZ_NORM:
        case QIR_OP_CTRL_SEL:
        case QIR_OP_SHL_VAR:
        case QIR_OP_SHR_VAR:
        case QIR_OP_EQ_ZERO:
            return true;
        case QIR_OP_LINF2:
            return !inplace;
        default:
            return false;
    }
}

static inline void qir__query_lut_res(
    const qir_lut_payload_t *lut,
    uint16_t vb,
    qir_variant_kind_t var,
    uint32_t aux_param,
    const qir_gauge_lut_t *gauge_lut,
    bool uncomp,
    uint32_t *out_temp_anc,
    uint32_t *out_cbits
) {
    if (uncomp && var == QIR_VAR_HMR_PARITY_LUT && lut->key_bits >= 2u && lut->key_bits <= 13u) {
        uint32_t nm = lut->key_bits >> 1u;
        uint32_t nl = lut->key_bits - nm;
        uint32_t leaves = (UINT32_C(1) << nm) + (UINT32_C(1) << nl);
        uint32_t borrow = vb < leaves ? vb : leaves;
        *out_temp_anc = leaves - borrow;
        *out_cbits    = vb + leaves + 2u;
    } else if (var == QIR_VAR_QROAM_VSLOT0) {
        qk_qroam_vslot0_resources(lut->table, lut->num_entries, lut->val_words, 0u, vb,
                                  lut->key_bits, aux_param ? aux_param : 1u,
                                  out_temp_anc, out_cbits);
    } else if (var == QIR_VAR_QROAM_SPLIT) {
        uint32_t K = aux_param & 0xFFu;
        uint32_t b_sh = (aux_param >> 8u) & 0xFFu;
        if (K < 2u) K = 2u;
        if (b_sh == 0u) b_sh = 1u;
        uint32_t chunk = (vb + K - 1u) / K;
        if (chunk == 0u) chunk = 1u;
        uint32_t mx_anc = 0u, sum_cb = 0u;
        for (uint32_t lo = 0u; lo < vb; lo += chunk) {
            uint32_t hi = (lo + chunk < vb) ? (lo + chunk) : vb;
            uint32_t ca = 0u, cc = 0u;
            qk_qroam_vslot0_resources(lut->table, lut->num_entries, lut->val_words, lo, hi,
                                      lut->key_bits, b_sh, &ca, &cc);
            if (ca > mx_anc) mx_anc = ca;
            sum_cb += cc;
        }
        *out_temp_anc = mx_anc;
        *out_cbits    = sum_cb;
    } else if (var == QIR_VAR_LIN_GAUGE_TSP && gauge_lut) {
        uint32_t rb = gauge_lut->rank_bits;
        uint32_t extra = lut->key_bits;
        if (lut->key_bits >= 2u && lut->key_bits <= 13u) {
            uint32_t nm = lut->key_bits >> 1u;
            uint32_t nl = lut->key_bits - nm;
            uint32_t leaves = (UINT32_C(1) << nm) + (UINT32_C(1) << nl);
            uint32_t borrow = rb < leaves ? rb : leaves;
            if (leaves - borrow > extra) extra = leaves - borrow;
        }
        uint32_t cb1 = (lut->key_bits < 16u) ? (UINT32_C(1) << lut->key_bits) : 4096u;
        *out_temp_anc = rb + extra;
        *out_cbits    = cb1 + 64u;
    } else if (var == QIR_VAR_LUT_COL_SPLIT && aux_param >= 2u) {
        uint32_t cb1 = (lut->key_bits < 16u) ? (UINT32_C(1) << lut->key_bits) : 4096u;
        *out_temp_anc = lut->key_bits;
        *out_cbits    = aux_param * cb1;
    } else {
        *out_temp_anc = lut->key_bits;
        *out_cbits    = (lut->key_bits < 16u) ? (UINT32_C(1) << lut->key_bits) : 4096u;
    }
}

/* Query exact step temp ancilla requirement cleanly (matches qir_codegen.h & qir_kernels.h) */
static inline void qir__query_step_temp(
    const qir_module_t *m,
    const qir_phys_plan_t *plan_plans_only,
    const qir_sched_step_t *st,
    uint32_t *out_temp_anc,
    uint32_t *out_cbits,
    uint32_t *out_aux_qubits
) {
    *out_temp_anc   = 0u;
    *out_cbits      = 0u;
    *out_aux_qubits = 0u;
    const qir_op_t *op = &m->ops[st->op_id];
    uint16_t out_w = (op->num_outs > 0u) ? m->vals[op->outs[0]].width : 0u;
    bool uncomp = (st->kind == QIR_STEP_UNCOMPUTE);
    (void)plan_plans_only;

    switch (op->kind) {
        case QIR_OP_AND_VEC:
            if (uncomp) *out_cbits = out_w;
            break;
        case QIR_OP_ADD:
        case QIR_OP_SUB:
            if (st->variant == QIR_VAR_CSA_TREE && st->num_chain_ops > 1u) {
                uint32_t n_sum = st->aux_param >= 3u ? st->aux_param : 3u;
                uint32_t stages = n_sum - 2u;
                *out_temp_anc = 2u * stages * out_w + (out_w > 0u ? out_w - 1u : 0u);
                *out_cbits    = stages * out_w + 2u;
            } else {
                *out_temp_anc = (out_w > 1u) ? (out_w - 1u) : 0u;
                *out_cbits    = 1u;
            }
            break;
        case QIR_OP_CS_ADD3:
            if (uncomp) *out_cbits = out_w + 2u;
            break;
        case QIR_OP_MUL_FXP: {
            int32_t sh_s = (int32_t)op->u.mul.frac_a + (int32_t)op->u.mul.frac_b - (int32_t)op->u.mul.frac_out;
            uint32_t shift = sh_s > 0 ? (uint32_t)sh_s : 0u;
            uint32_t full  = shift + out_w;
            uint32_t keep  = (op->flags & QIR_FLAG_APPROX) ? st->approx_keep_bits : 0u;
            uint32_t minp  = (keep > 0u && keep < full) ? (full - keep) : 0u;
            if (minp > shift) minp = shift;
            uint32_t acc_bits = full - minp;
            uint32_t low_junk = shift - minp;
            *out_temp_anc   = acc_bits;
            *out_aux_qubits = uncomp ? 0u : low_junk;
            if (st->variant == QIR_VAR_MUL_CONST_ADD && op->num_ins >= 2u) {
                *out_cbits = m->vals[op->ins[st->aux_param & 1u]].width + 4u;
            } else {
                *out_cbits = m->vals[op->ins[0]].width + 4u;
            }
            break;
        }
        case QIR_OP_SQR_FXP: {
            uint16_t wa = m->vals[op->ins[0]].width;
            int32_t sh_s = 2 * (int32_t)op->u.mul.frac_a - (int32_t)op->u.mul.frac_out;
            uint32_t shift = sh_s > 0 ? (uint32_t)sh_s : 0u;
            uint32_t full  = shift + out_w;
            uint32_t keep  = (op->flags & QIR_FLAG_APPROX) ? st->approx_keep_bits : 0u;
            uint32_t minp  = (keep > 0u && keep < full) ? (full - keep) : 0u;
            if (minp > shift) minp = shift;
            uint32_t acc_bits = full - minp;
            uint32_t low_junk = shift - minp;
            *out_temp_anc   = (st->variant == QIR_VAR_SQR_TRIANGLE) ? acc_bits : (wa + acc_bits);
            *out_aux_qubits = uncomp ? 0u : low_junk;
            *out_cbits      = wa + 4u;
            break;
        }
        case QIR_OP_GF2M_ADD:
            if (st->variant == QIR_VAR_XOR_ACCUM && st->num_chain_ops > 0u) {
                const qir_op_t *dop = &m->ops[st->chain_ops[0]];
                qir_variant_kind_t sub_var = (qir_variant_kind_t)((st->aux_param >> 8u) & 0xFFu);
                uint32_t sub_aux = st->aux_param >> 16u;
                if (dop->kind == QIR_OP_LUT_ROM) {
                    qir__query_lut_res(&dop->u.lut, dop->u.lut.val_bits, sub_var, sub_aux, st->gauge_lut, false,
                                       out_temp_anc, out_cbits);
                } else if (dop->kind == QIR_OP_GF2M_MUL) {
                    *out_temp_anc = 1u;
                    *out_cbits    = (dop->u.gf2m.m == 12u || dop->u.gf2m.m == 8u || dop->u.gf2m.m == 6u)
                        ? 4u : (uint32_t)dop->u.gf2m.m * dop->u.gf2m.m;
                }
            }
            break;
        case QIR_OP_GF2M_MUL:
            if (st->variant == QIR_VAR_GF2M_CONST_LINF2) {
                break;
            }
            *out_temp_anc = 1u;
            *out_cbits    = (uncomp && (op->u.gf2m.m == 12u || op->u.gf2m.m == 8u || op->u.gf2m.m == 6u))
                ? (op->u.gf2m.m + 1u)
                : ((op->u.gf2m.m == 12u || op->u.gf2m.m == 8u || op->u.gf2m.m == 6u) ? 4u : (uint32_t)op->u.gf2m.m * op->u.gf2m.m);
            break;
        case QIR_OP_GF2M_INV:
            if (op->u.gf2m.m == 12u) {
                *out_temp_anc = (op->num_ins >= 2u) ? 22u : 19u;
                *out_cbits    = 48u;
            } else if (op->u.gf2m.m == 8u && st->variant == QIR_VAR_GF2M_TOWER_INV &&
                       plan_plans_only && plan_plans_only->has_gf256_tower &&
                       plan_plans_only->gf256_tower_poly == op->u.gf2m.irred_poly) {
                *out_temp_anc = (op->num_ins >= 2u) ? 13u : 11u;
                *out_cbits    = 32u;
            } else if (op->num_ins >= 2u) {
                *out_temp_anc = 2u * op->u.gf2m.m;
                *out_cbits    = (UINT32_C(1) << op->u.gf2m.m) * 2u;
            } else {
                *out_temp_anc = op->u.gf2m.m;
                *out_cbits    = UINT32_C(1) << op->u.gf2m.m;
            }
            break;
        case QIR_OP_LUT_ROM: {
            const qir_lut_payload_t *lut = &op->u.lut;
            uint16_t vb = lut->val_bits;
            if ((op->flags & QIR_FLAG_APPROX) && st->approx_keep_bits > 0u && st->approx_keep_bits < vb) {
                vb = st->approx_keep_bits;
            }
            qir__query_lut_res(lut, vb, st->variant, st->aux_param, st->gauge_lut, uncomp, out_temp_anc, out_cbits);
            break;
        }
        case QIR_OP_CTRL_SEL:
            if (st->variant == QIR_VAR_CTRL_ADD_CUCCARO && st->num_chain_ops > 0u) {
                uint32_t lhs_vid = st->aux_param & UINT32_C(0x7FFFFFFF);
                uint16_t k_in = m->vals[lhs_vid].width;
                uint32_t k = (k_in < out_w) ? k_in : out_w;
                *out_temp_anc = (k == 1u) ? ((out_w == 1u) ? 0u : out_w - 1u) : ((k == out_w) ? 1u : (out_w - k));
                *out_cbits    = 2u;
            }
            break;
        case QIR_OP_CLZ_NORM: {
            uint16_t iw = m->vals[op->ins[0]].width;
            uint16_t cw = m->vals[op->outs[1]].width;
            *out_temp_anc = 2u * iw + cw + 4u;
            *out_cbits    = 6u * iw;
            break;
        }
        case QIR_OP_SHL_VAR: {
            uint16_t w_in = (op->num_ins > 0u) ? m->vals[op->ins[0]].width : 0u;
            uint32_t W = ((uint32_t)w_in + op->u.shift.max_shift > out_w)
                ? ((uint32_t)w_in + op->u.shift.max_shift) : out_w;
            *out_aux_qubits = (!uncomp && W > out_w) ? (W - out_w) : 0u;
            break;
        }
        case QIR_OP_SHR_VAR: {
            uint16_t w_in = (op->num_ins > 0u) ? m->vals[op->ins[0]].width : 0u;
            uint32_t base_w = (w_in > out_w) ? w_in : out_w;
            uint32_t W = base_w + op->u.shift.max_shift;
            *out_aux_qubits = (!uncomp && W > out_w) ? (W - out_w) : 0u;
            break;
        }
        case QIR_OP_EQ_ZERO: {
            uint16_t iw = (op->num_ins > 0u) ? m->vals[op->ins[0]].width : 0u;
            *out_temp_anc = (iw >= 3u) ? (uint32_t)(iw - 2u) : 0u;
            *out_cbits    = 2u;
            break;
        }
        case QIR_OP_PHASE_DOT:
            *out_cbits = (op->num_ins >= 3u) ? 2u : 0u;
            break;
        case QIR_OP_PHASE_POLY: {
            uint16_t max_deg = 0u;
            for (uint32_t t = 0u; t < op->u.phase_poly.num_terms; ++t) {
                if (!op->u.phase_poly.terms[t].is_xor) {
                    uint16_t pc = (uint16_t)qcc_bv_popcount(op->u.phase_poly.terms[t].mask);
                    if (pc > max_deg) max_deg = pc;
                }
            }
            if (max_deg >= 2u) {
                *out_temp_anc = max_deg - 1u;
                *out_cbits    = op->u.phase_poly.num_terms;
            }
            break;
        }
        case QIR_OP_OPAQUE_KMB:
            if (op->u.opaque.kmb) *out_cbits = op->u.opaque.kmb->num_bits + 16u;
            break;
        case QIR_OP_HMR_MEAS:
        case QIR_OP_BIT_OP:
        default:
            break;
    }
}

/* Dynamic Qubit Pool Bitset (up to 65536 qubits cleanly) */
#define QIR_POOL_WORDS 1024u
#define QIR_POOL_MAX_Q (QIR_POOL_WORDS * 64u)
typedef struct {
    uint64_t dirty[QIR_POOL_WORDS];     /* Bit q=1 => currently allocated live / non-zero */
    uint64_t dirty_perm[QIR_POOL_WORDS];/* Bit q=1 => permanently dirty tail (never borrowable) */
} qir__qpool_t;

static inline uint32_t qir__pool_alloc_clean(qir__qpool_t *p, uint32_t *out_peak) {
    for (uint32_t w = 0u; w < QIR_POOL_WORDS; ++w) {
        uint64_t busy = p->dirty[w] | p->dirty_perm[w];
        if (busy != UINT64_MAX) {
            uint32_t b = (uint32_t)__builtin_ctzll((unsigned long long)(~busy));
            uint32_t q = 64u * w + b;
            p->dirty[w] |= (UINT64_C(1) << b);
            if (out_peak && q + 1u > *out_peak) *out_peak = q + 1u;
            return q;
        }
    }
    return QIR_POOL_MAX_Q - 1u;
}

static inline void qir__pool_free_clean(qir__qpool_t *p, uint32_t q) {
    if (q < QIR_POOL_MAX_Q) p->dirty[q >> 6u] &= ~(UINT64_C(1) << (q & 63u));
}

static inline uint16_t qir__collect_op_inputs(
    const qir_module_t *m,
    const qir_enode_t *en,
    const qir_op_t *op,
    const uint8_t *op_skip,
    uint32_t out_ins[32]
) {
    (void)op_skip;
    uint16_t n = 0u;
    if (en && en->variant == QIR_VAR_LINF2_GAUSS && en->fused_linf2 && en->num_chain_ops > 0u) {
        const qir_op_t *prev = &m->ops[en->aux_param];
        for (uint16_t i = 0u; i < prev->num_ins && n < 32u; ++i) out_ins[n++] = prev->ins[i];
    } else if (en && en->variant == QIR_VAR_CSA_TREE && en->num_chain_ops > 1u) {
        for (uint16_t ci = 0u; ci < en->num_chain_ops; ++ci) {
            const qir_op_t *cop = &m->ops[en->chain_ops[ci]];
            for (uint16_t i = 0u; i < cop->num_ins; ++i) {
                uint32_t iv = cop->ins[i];
                bool inside = false;
                if (m->vals[iv].def_op != UINT32_MAX) {
                    for (uint16_t cj = 0u; cj < en->num_chain_ops; ++cj) {
                        if (en->chain_ops[cj] == m->vals[iv].def_op) { inside = true; break; }
                    }
                }
                if (!inside && n < 32u) out_ins[n++] = iv;
            }
        }
    } else if (en && en->variant == QIR_VAR_PHASE_CSA_DIRECT) {
        if (en->phase_split) {
            const qir_phase_split_t *ps = en->phase_split;
            for (uint16_t li = 0u; li < ps->num_leaves && n < 32u; ++li) {
                out_ins[n++] = ps->leaf_vals[li];
            }
            for (uint16_t i = 0u; i < op->num_ins && n < 32u; ++i) {
                if (op->ins[i] != ps->target_val) out_ins[n++] = op->ins[i];
            }
        } else {
            const qir_op_t *dead_add = &m->ops[en->aux_param];
            out_ins[n++] = dead_add->ins[0];
            out_ins[n++] = dead_add->ins[1];
            for (uint16_t i = 0u; i < op->num_ins && n < 32u; ++i) {
                if (op->ins[i] != dead_add->outs[0]) out_ins[n++] = op->ins[i];
            }
        }
    } else if (en && en->variant == QIR_VAR_MUL_CONST_ADD && op->num_ins >= 2u) {
        out_ins[n++] = op->ins[(en->aux_param & 1u) ^ 1u];
    } else if (en && en->variant == QIR_VAR_GF2M_CONST_LINF2 && op->num_ins >= 2u) {
        out_ins[n++] = op->ins[en->aux_param & 1u];
    } else if (en && en->variant == QIR_VAR_CTRL_ADD_CUCCARO && op->num_ins >= 3u) {
        out_ins[n++] = op->ins[0];
        out_ins[n++] = op->ins[1];
        out_ins[n++] = en->aux_param & UINT32_C(0x7FFFFFFF);
    } else if (en && en->variant == QIR_VAR_XOR_ACCUM && en->num_chain_ops > 0u && op->num_ins >= 2u) {
        out_ins[n++] = op->ins[en->aux_param & 3u];
        const qir_op_t *dop = &m->ops[en->chain_ops[0]];
        for (uint16_t i = 0u; i < dop->num_ins && n < 32u; ++i) out_ins[n++] = dop->ins[i];
    } else {
        for (uint16_t i = 0u; i < op->num_ins && n < 32u; ++i) out_ins[n++] = op->ins[i];
    }
    if (op->pred_val != UINT32_MAX && n < 32u) {
        out_ins[n++] = op->pred_val;
    }
    return n;
}

static inline bool qir__pool_is_busy(const qir__qpool_t *p, uint32_t q) {
    if (q >= QIR_POOL_MAX_Q) return true;
    return (((p->dirty[q >> 6u] | p->dirty_perm[q >> 6u]) >> (q & 63u)) & 1u) != 0u;
}

static inline bool qir__can_exec_inplace(
    const qir_module_t *m,
    const qir_phys_plan_t *plan,
    const uint8_t *op_skip,
    const uint8_t *will_alias,
    const uint8_t *is_in,
    const uint8_t *is_pinned_out,
    const uint8_t *is_direct_out,
    const uint32_t *fwd_refs,
    const uint32_t *uncomp_pin,
    const uint32_t *back_vals,
    const uint16_t *back_cnt,
    const qir_op_t *op,
    uint32_t src_v,
    const uint32_t *other_vals,
    uint32_t n_others,
    bool eff_clean
) {
    if (op->num_outs == 0u || src_v >= m->num_vals) return false;
    uint32_t dst = op->outs[0];
    if (m->vals[src_v].rep != QIR_REP_STD || m->vals[dst].rep != QIR_REP_STD) return false;
    uint16_t w = m->vals[dst].width;
    if (w == 0u || m->vals[src_v].width != w) return false;

    const qir_val_alloc_t *sva = &plan->val_allocs[src_v];
    if (!sva->qubits) return false;
    for (uint16_t a = 0u; a < w; ++a) {
        for (uint16_t b = a + 1u; b < w; ++b) {
            if (sva->qubits[a] == sva->qubits[b]) return false;
        }
    }
    for (uint32_t oi = 0u; oi < n_others; ++oi) {
        uint32_t ov = other_vals[oi];
        if (ov == src_v || ov >= m->num_vals) continue;
        const qir_val_alloc_t *ova = &plan->val_allocs[ov];
        if (!ova->qubits) continue;
        for (uint16_t a = 0u; a < w; ++a) {
            for (uint16_t b = 0u; b < ova->width; ++b) {
                if (sva->qubits[a] == ova->qubits[b]) return false;
            }
        }
    }
    if (back_cnt[src_v] == 0u) return false;
    for (uint16_t bi = 0u; bi < back_cnt[src_v]; ++bi) {
        uint32_t bv = back_vals[(size_t)src_v * 64u + bi];
        if (fwd_refs[bv] != 1u || uncomp_pin[bv] != 0u || is_pinned_out[bv]) return false;
    }

    if (!eff_clean) {
        if (sva->is_alias || will_alias[src_v]) return false;
        uint32_t bv = src_v;
        return (is_in[bv] || (m->vals[bv].def_op != UINT32_MAX &&
                              !qir__op_supports_uncompute(&m->ops[m->vals[bv].def_op], false)));
    }

    for (uint16_t bi = 0u; bi < back_cnt[src_v]; ++bi) {
        uint32_t bv = back_vals[(size_t)src_v * 64u + bi];
        if (!is_in[bv]) return false;
    }

    /* Map every bit of src_v to its module input register (vid, bit) */
    uint32_t *bit_in_vid = (uint32_t *)qcc_xmalloc((size_t)w * sizeof(uint32_t));
    uint16_t *bit_in_pos = (uint16_t *)qcc_xmalloc((size_t)w * sizeof(uint16_t));
    for (uint16_t b = 0u; b < w; ++b) {
        bit_in_vid[b] = UINT32_MAX;
        for (uint32_t ii = 0u; ii < m->num_inputs; ++ii) {
            uint32_t iv = m->inputs[ii];
            const qir_val_alloc_t *iva = &plan->val_allocs[iv];
            if (!iva->qubits) continue;
            for (uint16_t k = 0u; k < iva->width; ++k) {
                if (iva->qubits[k] == sva->qubits[b]) {
                    bit_in_vid[b] = iv;
                    bit_in_pos[b] = k;
                    break;
                }
            }
            if (bit_in_vid[b] != UINT32_MAX) break;
        }
        if (bit_in_vid[b] == UINT32_MAX) {
            free(bit_in_vid); free(bit_in_pos);
            return false;
        }
    }

    bool ok = false;
    if (is_direct_out[dst] && back_cnt[src_v] == 1u && is_in[src_v]) {
        const char *dn = m->vals[dst].name;
        const char *sn = m->vals[src_v].name;
        bool has_cons = false;
        for (uint32_t ci = op->id + 1u; ci < m->num_ops && !has_cons; ++ci) {
            if (op_skip[ci]) continue;
            for (uint16_t ii = 0u; ii < m->ops[ci].num_ins; ++ii) {
                if (m->ops[ci].ins[ii] == dst) { has_cons = true; break; }
            }
        }
        if (!has_cons && dn && sn && strcmp(dn, sn) == 0) ok = true;
    } else if (!is_direct_out[dst]) {
        uint8_t *cov = (uint8_t *)qcc_xcalloc(w, 1u);
        ok = true;
        for (uint32_t ci = op->id + 1u; ci < m->num_ops && ok; ++ci) {
            if (op_skip[ci]) continue;
            const qir_op_t *cop = &m->ops[ci];
            bool uses_dst = false;
            for (uint16_t ii = 0u; ii < cop->num_ins; ++ii) {
                if (cop->ins[ii] == dst) { uses_dst = true; break; }
            }
            if (!uses_dst) continue;
            if (cop->kind != QIR_OP_SLICE || cop->num_outs == 0u ||
                !will_alias[cop->outs[0]] || !is_direct_out[cop->outs[0]]) {
                ok = false; break;
            }
            uint32_t so = cop->outs[0];
            uint16_t lo = cop->u.slice.lo, sw = m->vals[so].width;
            const char *onm = m->vals[so].name;
            if (lo + sw > w || !onm) { ok = false; break; }
            for (uint16_t sb = 0u; sb < sw; ++sb) {
                uint32_t iv = bit_in_vid[lo + sb];
                const char *inm = m->vals[iv].name;
                if (cov[lo + sb] != 0u || bit_in_pos[lo + sb] != sb ||
                    m->vals[iv].width != sw || !inm || strcmp(onm, inm) != 0) {
                    ok = false; break;
                }
                cov[lo + sb] = 1u;
            }
        }
        for (uint16_t b = 0u; b < w && ok; ++b) {
            if (cov[b] == 0u) ok = false;
        }
        free(cov);
    }
    free(bit_in_vid);
    free(bit_in_pos);
    return ok;
}

static inline void qir__heap_push(uint32_t *heap, uint32_t *len, uint8_t *in_ready, uint32_t val) {
    if (in_ready[val]) return;
    in_ready[val] = 1u;
    uint32_t i = (*len)++;
    heap[i] = val;
    while (i > 0u) {
        uint32_t p = (i - 1u) >> 1u;
        if (heap[p] >= heap[i]) break;
        uint32_t t = heap[p]; heap[p] = heap[i]; heap[i] = t;
        i = p;
    }
}

static inline uint32_t qir__heap_pop(uint32_t *heap, uint32_t *len, uint8_t *in_ready) {
    uint32_t top = heap[0];
    in_ready[top] = 0u;
    uint32_t n = --(*len);
    if (n > 0u) {
        heap[0] = heap[n];
        uint32_t i = 0u;
        for (;;) {
            uint32_t l = 2u * i + 1u, r = l + 1u, best = i;
            if (l < n && heap[l] > heap[best]) best = l;
            if (r < n && heap[r] > heap[best]) best = r;
            if (best == i) break;
            uint32_t t = heap[i]; heap[i] = heap[best]; heap[best] = t;
            i = best;
        }
    }
    return top;
}

static inline void qir__maybe_push_uncomp(
    const qir_module_t *m,
    const uint8_t *computed,
    const uint8_t *uncomputed,
    const uint8_t *can_uncomp,
    const uint32_t *fwd_refs,
    const uint32_t *uncomp_pin,
    const uint8_t *is_pinned_out,
    uint32_t uo,
    uint32_t *heap,
    uint32_t *len,
    uint8_t *in_ready
) {
    if (uo == UINT32_MAX || uo >= m->num_ops) return;
    if (!computed[uo] || uncomputed[uo] || !can_uncomp[uo] || in_ready[uo]) return;
    const qir_op_t *dop = &m->ops[uo];
    for (uint8_t ko = 0u; ko < dop->num_outs; ++ko) {
        uint32_t ov = dop->outs[ko];
        if (fwd_refs[ov] > 0u || uncomp_pin[ov] > 0u || is_pinned_out[ov]) return;
    }
    qir__heap_push(heap, len, in_ready, uo);
}

static inline void qir__alloc_step_temp(
    const qir_module_t *m,
    const qir_phys_plan_t *plan,
    qir_sched_step_t *st,
    qcc_arena_t *arena,
    qir__qpool_t *pool,
    uint32_t *peak_q,
    bool *certified
) {
    uint32_t need_temp = st->num_temp_anc;
    if (need_temp == 0u) return;
    const qir_op_t *op = &m->ops[st->op_id];
    bool uncomp = (st->kind == QIR_STEP_UNCOMPUTE);
    bool can_borrow_zero = false;
    if (op->num_outs > 0u) {
        if (op->kind == QIR_OP_LUT_ROM && st->variant != QIR_VAR_LIN_GAUGE_TSP &&
            !(uncomp && st->variant == QIR_VAR_HMR_PARITY_LUT)) {
            can_borrow_zero = true;
        } else if (op->kind == QIR_OP_GF2M_ADD && st->variant == QIR_VAR_XOR_ACCUM && st->num_chain_ops > 0u) {
            const qir_op_t *dop = &m->ops[st->chain_ops[0]];
            uint32_t sub_var = (st->aux_param >> 8u) & 0xFFu;
            if (dop->kind == QIR_OP_LUT_ROM && sub_var != QIR_VAR_LIN_GAUGE_TSP) {
                can_borrow_zero = true;
            }
        }
    }

    st->temp_anc = (uint32_t *)qcc_arena_xalloc(arena, (size_t)need_temp * sizeof(uint32_t), 4u);
    uint8_t *borrowed = (uint8_t *)qcc_xcalloc(need_temp, 1u);
    uint32_t filled = 0u;
    if (can_borrow_zero) {
        uint32_t vid0 = op->outs[0];
        const qir_val_alloc_t *va0 = &plan->val_allocs[vid0];
        const qir_val_t *ov = &m->vals[vid0];
        if (va0->qubits) {
            for (uint16_t b = 0u; b < va0->width && filled < need_temp; ++b) {
                if (qcc_bv_get(ov->zero_mask, b)) {
                    borrowed[filled] = 1u;
                    st->temp_anc[filled++] = va0->qubits[b];
                }
            }
        }
    }
    while (filled < need_temp) {
        st->temp_anc[filled++] = qir__pool_alloc_clean(pool, peak_q);
    }
    for (uint32_t a = 0u; a < need_temp; ++a) {
        for (uint32_t b = a + 1u; b < need_temp; ++b) {
            if (st->temp_anc[a] == st->temp_anc[b]) *certified = false;
        }
        if (!borrowed[a] && op->num_outs > 0u && !st->inplace_on_in0) {
            const qir_val_alloc_t *va0 = &plan->val_allocs[op->outs[0]];
            if (va0->qubits) {
                for (uint16_t b = 0u; b < va0->width; ++b) {
                    if (st->temp_anc[a] == va0->qubits[b]) *certified = false;
                }
            }
            if (va0->aux_qubits) {
                for (uint16_t b = 0u; b < va0->aux_width; ++b) {
                    if (st->temp_anc[a] == va0->aux_qubits[b]) *certified = false;
                }
            }
        }
    }
    for (uint32_t a = 0u; a < need_temp; ++a) {
        if (!borrowed[a]) qir__pool_free_clean(pool, st->temp_anc[a]);
    }
    free(borrowed);
}

/* Construct physical schedule + execute Bit-Level Borrow Allocator + Disjointness Certificate */
static inline bool qir_regalloc_pack(
    const qir_egraph_t *eg,
    const uint32_t *chosen_enodes_in,
    const uint8_t *uncomp_flags,
    const uint16_t *approx_widths,
    bool do_uncompute,
    bool force_clean,
    qir_phys_plan_t *out_plan,
    uint32_t *out_peak_op_at_max
) {
    const qir_module_t *m = eg->mod;
    memset(out_plan, 0, sizeof(*out_plan));
    out_plan->arena = eg->arena;
    out_plan->mod   = m;
    out_plan->eg    = eg;

    size_t nv = m->num_vals ? m->num_vals : 1u;
    size_t no = m->num_ops ? m->num_ops : 1u;
    bool eff_clean = force_clean || m->require_clean_scratch;

    out_plan->chosen_enodes = (uint32_t *)qcc_xmalloc(nv * sizeof(uint32_t));
    out_plan->val_allocs    = (qir_val_alloc_t *)qcc_xcalloc(nv, sizeof(qir_val_alloc_t));
    if (m->num_vals > 0u) memcpy(out_plan->chosen_enodes, chosen_enodes_in, nv * sizeof(uint32_t));

    /* Pre-build cached GF plans when needed */
    for (uint32_t oi = 0u; oi < m->num_ops; ++oi) {
        const qir_op_t *op = &m->ops[oi];
        if ((op->kind == QIR_OP_GF2M_MUL || (op->kind == QIR_OP_PHASE_DOT && op->num_ins >= 3u)) &&
            op->u.gf2m.m == 12u && !out_plan->has_gf12_r42) {
            out_plan->has_gf12_r42 = qir_build_gf12_rank42(op->u.gf2m.irred_poly, &out_plan->gf12_r42);
            if (out_plan->has_gf12_r42) out_plan->gf12_r42_poly = op->u.gf2m.irred_poly;
        } else if ((op->kind == QIR_OP_GF2M_MUL || (op->kind == QIR_OP_PHASE_DOT && op->num_ins >= 3u)) &&
                   op->u.gf2m.m == 8u && !out_plan->has_gf256_r27) {
            out_plan->has_gf256_r27 = qir_build_gf256_rank27(op->u.gf2m.irred_poly, &out_plan->gf256_r27);
            if (out_plan->has_gf256_r27) out_plan->gf256_r27_poly = op->u.gf2m.irred_poly;
        } else if ((op->kind == QIR_OP_GF2M_MUL || (op->kind == QIR_OP_PHASE_DOT && op->num_ins >= 3u)) &&
                   op->u.gf2m.m == 6u && !out_plan->has_gf64_r15) {
            out_plan->has_gf64_r15 = qir_build_gf64_rank15(op->u.gf2m.irred_poly, &out_plan->gf64_r15);
            if (out_plan->has_gf64_r15) out_plan->gf64_r15_poly = op->u.gf2m.irred_poly;
        } else if (op->kind == QIR_OP_GF2M_INV && op->u.gf2m.m == 12u && !out_plan->has_gf12_tower) {
            out_plan->has_gf12_tower = qir_build_gf12_tower_inv(op->u.gf2m.irred_poly, &out_plan->gf12_tower);
            if (out_plan->has_gf12_tower) out_plan->gf12_tower_poly = op->u.gf2m.irred_poly;
        } else if (op->kind == QIR_OP_GF2M_INV && op->u.gf2m.m == 8u && !out_plan->has_gf256_tower) {
            out_plan->has_gf256_tower = qir_build_gf256_tower_inv(op->u.gf2m.irred_poly, &out_plan->gf256_tower);
            if (out_plan->has_gf256_tower) out_plan->gf256_tower_poly = op->u.gf2m.irred_poly;
        }
    }

    uint8_t  *op_skip       = (uint8_t *)qcc_xcalloc(no, 1u);
    uint8_t  *will_alias    = (uint8_t *)qcc_xcalloc(nv, 1u);
    uint8_t  *is_direct_out = (uint8_t *)qcc_xcalloc(nv, 1u);
    uint8_t  *is_in         = (uint8_t *)qcc_xcalloc(nv, 1u);
    uint8_t  *is_pinned_out = (uint8_t *)qcc_xcalloc(nv, 1u);
    uint8_t  *hmr_measured  = (uint8_t *)qcc_xcalloc(nv, 1u);
    uint8_t  *can_uncomp    = (uint8_t *)qcc_xcalloc(no, 1u);
    uint8_t  *computed      = (uint8_t *)qcc_xcalloc(no, 1u);
    uint8_t  *uncomputed    = (uint8_t *)qcc_xcalloc(no, 1u);
    uint32_t *fwd_refs      = (uint32_t *)qcc_xcalloc(nv, sizeof(uint32_t));
    uint32_t *uncomp_pin    = (uint32_t *)qcc_xcalloc(nv, sizeof(uint32_t));
    uint32_t *back_vals     = (uint32_t *)qcc_xcalloc(nv * 64u, sizeof(uint32_t));
    uint16_t *back_cnt      = (uint16_t *)qcc_xcalloc(nv, sizeof(uint16_t));
    uint32_t *step_back_ins = (uint32_t *)qcc_xcalloc(no * 64u, sizeof(uint32_t));
    uint16_t *step_back_cnt = (uint16_t *)qcc_xcalloc(no, sizeof(uint16_t));

    /* Enforce coupled rewrite invariants in reverse topological order so downstream folds take precedence */
    for (int32_t oi_s = (int32_t)m->num_ops - 1; oi_s >= 0; --oi_s) {
        uint32_t oi = (uint32_t)oi_s;
        if (op_skip[oi]) continue;
        const qir_op_t *op = &m->ops[oi];
        if (op->num_outs == 0u) continue;
        const qir_enode_t *en = &eg->enodes[out_plan->chosen_enodes[op->outs[0]]];
        if (en->variant == QIR_VAR_TRIVIAL_ELIDE) {
            if (en->rule_tag == 4u) {
                uint32_t cons_vid = m->ops[en->aux_param].outs[0];
                if (eg->enodes[out_plan->chosen_enodes[cons_vid]].variant == QIR_VAR_PHASE_CSA_DIRECT) {
                    op_skip[oi] = 1u;
                } else {
                    out_plan->chosen_enodes[op->outs[0]] = eg->classes[op->outs[0]].enodes[0];
                }
            } else {
                op_skip[oi] = 1u;
            }
        }
        if (op_skip[oi]) continue;
        if (en->variant == QIR_VAR_PHASE_CSA_DIRECT) {
            if (en->num_chain_ops > 0u) {
                for (uint16_t ci = 0u; ci < en->num_chain_ops; ++ci) op_skip[en->chain_ops[ci]] = 1u;
            } else {
                op_skip[en->aux_param] = 1u;
            }
        }
        if (en->variant == QIR_VAR_LINF2_GAUSS && en->fused_linf2 && en->num_chain_ops > 0u) {
            for (uint16_t ci = 0u; ci < en->num_chain_ops; ++ci) op_skip[en->chain_ops[ci]] = 1u;
        }
        if (en->variant == QIR_VAR_CSA_TREE && en->num_chain_ops > 1u) {
            for (uint16_t ci = 0u; ci + 1u < en->num_chain_ops; ++ci) op_skip[en->chain_ops[ci]] = 1u;
        }
        if (en->variant == QIR_VAR_MUL_CONST_ADD && en->num_chain_ops > 0u) {
            for (uint16_t ci = 0u; ci < en->num_chain_ops; ++ci) op_skip[en->chain_ops[ci]] = 1u;
        }
        if ((en->variant == QIR_VAR_CTRL_ADD_CUCCARO ||
             en->variant == QIR_VAR_XOR_ACCUM ||
             en->variant == QIR_VAR_GF2M_CONST_LINF2) && en->num_chain_ops > 0u) {
            for (uint16_t ci = 0u; ci < en->num_chain_ops; ++ci) op_skip[en->chain_ops[ci]] = 1u;
        }
    }

    for (uint32_t i = 0u; i < m->num_inputs; ++i) is_in[m->inputs[i]] = 1u;
    for (uint32_t i = 0u; i < m->num_outputs; ++i) is_direct_out[m->outputs[i]] = 1u;
    for (uint32_t i = 0u; i < m->num_classical_outputs; ++i) {
        if (m->classical_outputs[i] < m->num_vals) is_direct_out[m->classical_outputs[i]] = 1u;
    }

    /* Backing-Set Alias Analysis through SLICE / CONCAT */
    for (uint32_t v = 0u; v < m->num_vals; ++v) {
        uint32_t def_oi = m->vals[v].def_op;
        bool try_alias = (def_oi != UINT32_MAX && !op_skip[def_oi] &&
                          m->ops[def_oi].pred_val == UINT32_MAX &&
                          (m->ops[def_oi].kind == QIR_OP_SLICE || m->ops[def_oi].kind == QIR_OP_CONCAT));
        if (try_alias) {
            const qir_op_t *dop = &m->ops[def_oi];
            uint16_t cnt = 0u;
            bool fits = true;
            for (uint16_t ii = 0u; ii < dop->num_ins && fits; ++ii) {
                uint32_t iv = dop->ins[ii];
                for (uint16_t bi = 0u; bi < back_cnt[iv]; ++bi) {
                    uint32_t bv = back_vals[(size_t)iv * 64u + bi];
                    bool dup = false;
                    for (uint16_t k = 0u; k < cnt; ++k) {
                        if (back_vals[(size_t)v * 64u + k] == bv) { dup = true; break; }
                    }
                    if (!dup) {
                        if (cnt >= 64u) { fits = false; break; }
                        back_vals[(size_t)v * 64u + cnt++] = bv;
                    }
                }
            }
            if (fits && cnt > 0u) {
                will_alias[v] = 1u;
                back_cnt[v]   = cnt;
                continue;
            }
        }
        will_alias[v]              = 0u;
        back_vals[(size_t)v * 64u] = v;
        back_cnt[v]                = 1u;
    }

    for (uint32_t i = 0u; i < m->num_outputs; ++i) {
        uint32_t ov = m->outputs[i];
        for (uint16_t bi = 0u; bi < back_cnt[ov]; ++bi) {
            is_pinned_out[back_vals[(size_t)ov * 64u + bi]] = 1u;
        }
    }
    for (uint32_t i = 0u; i < m->num_classical_outputs; ++i) {
        uint32_t ov = m->classical_outputs[i];
        if (ov < m->num_vals) {
            for (uint16_t bi = 0u; bi < back_cnt[ov]; ++bi) {
                is_pinned_out[back_vals[(size_t)ov * 64u + bi]] = 1u;
            }
        }
    }

    /* Precompute deduplicated backing inputs and refcounts for every non-skipped non-alias op */
    for (uint32_t oi = 0u; oi < m->num_ops; ++oi) {
        if (op_skip[oi]) continue;
        const qir_op_t *op = &m->ops[oi];
        if ((op->kind == QIR_OP_SLICE || op->kind == QIR_OP_CONCAT) &&
            op->num_outs > 0u && will_alias[op->outs[0]]) {
            continue;
        }
        const qir_enode_t *en = (op->num_outs > 0u) ? &eg->enodes[out_plan->chosen_enodes[op->outs[0]]] : NULL;
        uint32_t op_ins[32];
        uint16_t n_ins = qir__collect_op_inputs(m, en, op, op_skip, op_ins);
        uint16_t cnt = 0u;
        for (uint16_t ii = 0u; ii < n_ins; ++ii) {
            uint32_t iv = op_ins[ii];
            for (uint16_t bi = 0u; bi < back_cnt[iv]; ++bi) {
                uint32_t bv = back_vals[(size_t)iv * 64u + bi];
                bool dup = false;
                for (uint16_t k = 0u; k < cnt; ++k) {
                    if (step_back_ins[(size_t)oi * 64u + k] == bv) { dup = true; break; }
                }
                if (!dup && cnt < 64u) {
                    step_back_ins[(size_t)oi * 64u + cnt++] = bv;
                    fwd_refs[bv]++;
                }
            }
        }
        step_back_cnt[oi] = cnt;
    }

    /* Elide dead unconditional QIR_OP_CONST, QIR_OP_SLICE, QIR_OP_CONCAT in reverse topological order */
    for (int32_t oi_s = (int32_t)m->num_ops - 1; oi_s >= 0; --oi_s) {
        uint32_t oi = (uint32_t)oi_s;
        if (op_skip[oi]) continue;
        const qir_op_t *op = &m->ops[oi];
        if ((op->kind == QIR_OP_CONST || op->kind == QIR_OP_SLICE || op->kind == QIR_OP_CONCAT) &&
            op->pred_val == UINT32_MAX && op->num_outs > 0u) {
            uint32_t ov = op->outs[0];
            if (fwd_refs[ov] == 0u && !is_pinned_out[ov] && !is_direct_out[ov]) {
                bool has_consumer = false;
                for (uint32_t ci = oi + 1u; ci < m->num_ops && !has_consumer; ++ci) {
                    if (op_skip[ci]) continue;
                    const qir_op_t *cop = &m->ops[ci];
                    const qir_enode_t *cen = (cop->num_outs > 0u)
                        ? &eg->enodes[out_plan->chosen_enodes[cop->outs[0]]] : NULL;
                    uint32_t cins[32];
                    uint16_t ncins = qir__collect_op_inputs(m, cen, cop, op_skip, cins);
                    for (uint16_t ii = 0u; ii < ncins; ++ii) {
                        if (cins[ii] == ov) { has_consumer = true; break; }
                    }
                    for (uint16_t ii = 0u; ii < cop->num_ins && !has_consumer; ++ii) {
                        if (cop->ins[ii] == ov) { has_consumer = true; break; }
                    }
                }
                if (!has_consumer) op_skip[oi] = 1u;
            }
        }
    }

    /* Pre-detect unconditional HMR_MEAS targets (HMR already resets non-alias measured values to |0>) */
    for (uint32_t oi = 0u; oi < m->num_ops; ++oi) {
        if (op_skip[oi]) continue;
        const qir_op_t *op = &m->ops[oi];
        if (op->kind == QIR_OP_HMR_MEAS && op->pred_val == UINT32_MAX && op->num_ins > 0u) {
            uint32_t iv = op->ins[0];
            if (!will_alias[iv]) hmr_measured[iv] = 1u;
        }
    }

    qir__qpool_t *pool = (qir__qpool_t *)calloc(1u, sizeof(qir__qpool_t));
    if (!pool) {
        free(op_skip); free(will_alias); free(is_direct_out); free(is_in); free(is_pinned_out);
        free(hmr_measured); free(can_uncomp); free(computed); free(uncomputed);
        free(fwd_refs); free(uncomp_pin); free(back_vals); free(back_cnt);
        free(step_back_ins); free(step_back_cnt);
        return false;
    }
    uint32_t peak_q = 0u;
    uint32_t live_peak_q = 0u;
    uint32_t peak_op = 0u;
    uint32_t cbit = 0u;

    /* Pin all module inputs cleanly at [0 .. total_in_bits) */
    for (uint32_t i = 0u; i < m->num_inputs; ++i) {
        uint32_t vid = m->inputs[i];
        uint16_t w = m->vals[vid].width;
        qir_val_alloc_t *va = &out_plan->val_allocs[vid];
        va->is_pinned_io = 1u;
        if (m->vals[vid].rep == QIR_REP_CLASSICAL) {
            va->cwidth = w;
            va->cbits  = w ? (uint32_t *)qcc_arena_xalloc(eg->arena, (size_t)w * sizeof(uint32_t), 4u) : NULL;
            for (uint16_t b = 0u; b < w; ++b) va->cbits[b] = cbit++;
        } else {
            va->width  = w;
            va->qubits = w ? (uint32_t *)qcc_arena_xalloc(eg->arena, (size_t)w * sizeof(uint32_t), 4u) : NULL;
            for (uint16_t b = 0u; b < w; ++b) {
                va->qubits[b] = qir__pool_alloc_clean(pool, &peak_q);
            }
        }
    }

    size_t max_steps = (size_t)m->num_ops * 2u + 16u;
    out_plan->steps  = (qir_sched_step_t *)qcc_xcalloc(max_steps, sizeof(qir_sched_step_t));
    uint32_t *ready_heap = (uint32_t *)qcc_xmalloc(no * sizeof(uint32_t));
    uint8_t  *in_ready   = (uint8_t *)qcc_xcalloc(no, 1u);
    uint32_t  rh_len     = 0u;
    uint32_t ns      = 0u;
    bool certified   = true;

    for (uint32_t oi = 0u; oi < m->num_ops; ++oi) {
        if (op_skip[oi]) continue;
        const qir_op_t *op = &m->ops[oi];
        const qir_enode_t *en = (op->num_outs > 0u) ? &eg->enodes[out_plan->chosen_enodes[op->outs[0]]] : NULL;

        /* Handle zero-cost structural wire aliases (SLICE / CONCAT) */
        if (op->kind == QIR_OP_SLICE && op->num_outs > 0u && will_alias[op->outs[0]]) {
            uint32_t src = op->ins[0];
            uint32_t dst = op->outs[0];
            uint16_t w = m->vals[dst].width;
            const qir_val_alloc_t *sva = &out_plan->val_allocs[src];
            qir_val_alloc_t *dva = &out_plan->val_allocs[dst];
            dva->is_alias = 1u;
            if (m->vals[dst].rep == QIR_REP_CLASSICAL) {
                if (!sva->cbits) continue;
                dva->cwidth = w;
                dva->cbits  = w ? (uint32_t *)qcc_arena_xalloc(eg->arena, (size_t)w * sizeof(uint32_t), 4u) : NULL;
                for (uint16_t b = 0u; b < w; ++b) dva->cbits[b] = sva->cbits[op->u.slice.lo + b];
            } else {
                if (!sva->qubits) continue;
                dva->width  = w;
                dva->qubits = w ? (uint32_t *)qcc_arena_xalloc(eg->arena, (size_t)w * sizeof(uint32_t), 4u) : NULL;
                for (uint16_t b = 0u; b < w; ++b) dva->qubits[b] = sva->qubits[op->u.slice.lo + b];
            }
            continue;
        }
        if (op->kind == QIR_OP_CONCAT && op->num_outs > 0u && will_alias[op->outs[0]]) {
            uint32_t dst = op->outs[0];
            uint16_t w = m->vals[dst].width;
            qir_val_alloc_t *dva = &out_plan->val_allocs[dst];
            dva->is_alias = 1u;
            if (m->vals[dst].rep == QIR_REP_CLASSICAL) {
                dva->cwidth = w;
                dva->cbits  = w ? (uint32_t *)qcc_arena_xalloc(eg->arena, (size_t)w * sizeof(uint32_t), 4u) : NULL;
                uint16_t pos = 0u;
                for (uint16_t i = 0u; i < op->num_ins; ++i) {
                    const qir_val_alloc_t *iva = &out_plan->val_allocs[op->ins[i]];
                    if (!iva->cbits) { pos += m->vals[op->ins[i]].width; continue; }
                    for (uint16_t b = 0u; b < iva->cwidth && pos < w; ++b, ++pos) dva->cbits[pos] = iva->cbits[b];
                }
            } else {
                dva->width  = w;
                dva->qubits = w ? (uint32_t *)qcc_arena_xalloc(eg->arena, (size_t)w * sizeof(uint32_t), 4u) : NULL;
                uint16_t pos = 0u;
                for (uint16_t i = 0u; i < op->num_ins; ++i) {
                    const qir_val_alloc_t *iva = &out_plan->val_allocs[op->ins[i]];
                    if (!iva->qubits) { pos += m->vals[op->ins[i]].width; continue; }
                    for (uint16_t b = 0u; b < iva->width && pos < w; ++b, ++pos) dva->qubits[pos] = iva->qubits[b];
                }
            }
            continue;
        }

        qir_sched_step_t st;
        memset(&st, 0, sizeof(st));
        st.kind             = QIR_STEP_COMPUTE;
        st.op_id            = oi;
        st.variant          = en ? en->variant : QIR_VAR_STD;
        st.aux_param        = en ? en->aux_param : 0u;
        st.approx_keep_bits = (approx_widths && (op->flags & QIR_FLAG_APPROX)) ? approx_widths[oi] : (op->num_outs ? m->vals[op->outs[0]].width : 0u);
        st.fused_linf2      = en ? en->fused_linf2 : NULL;
        st.gauge_lut        = en ? en->gauge_lut : NULL;
        st.phase_split      = en ? en->phase_split : NULL;
        st.chain_ops        = en ? en->chain_ops : NULL;
        st.num_chain_ops    = en ? en->num_chain_ops : 0u;
        st.pred_cbit        = UINT32_MAX;
        st.inplace_in_idx   = 0xFFu;
        if (op->pred_val != UINT32_MAX && op->pred_val < m->num_vals && out_plan->val_allocs[op->pred_val].cbits) {
            st.pred_cbit = out_plan->val_allocs[op->pred_val].cbits[0];
        }

        /* Determine if op acts cleanly in-place on an operand (dying wire or same-name IO partner) */
        if (op->pred_val == UINT32_MAX && op->num_outs == 1u) {
            uint8_t cand_idx = 0xFFu;
            uint32_t cand_v = UINT32_MAX;
            uint32_t others[32];
            uint32_t n_others = 0u;
            bool allow = false;

            if (op->kind == QIR_OP_OPAQUE_KMB && op->num_ins == 1u) {
                cand_idx = 0u; cand_v = op->ins[0]; allow = true;
            } else if (op->kind == QIR_OP_LINF2 && st.variant == QIR_VAR_LINF2_GAUSS) {
                const qir_op_t *src_op = (st.fused_linf2 && st.num_chain_ops > 0u) ? &m->ops[st.aux_param] : op;
                if (src_op->num_ins == 1u) {
                    cand_idx = 0u; cand_v = src_op->ins[0];
                    uint16_t w = m->vals[op->outs[0]].width;
                    const qir_linf2_payload_t *lf = st.fused_linf2 ? st.fused_linf2 : &op->u.linf2;
                    qcc_bv256_t *M = (qcc_bv256_t *)qcc_xmalloc((size_t)(w ? w : 1u) * sizeof(qcc_bv256_t));
                    memcpy(M, lf->rows, (size_t)w * sizeof(qcc_bv256_t));
                    bool inv_ok = true;
                    for (uint32_t col = 0u; col < w; ++col) {
                        uint32_t piv = col;
                        while (piv < w && !qcc_bv_get(M[piv].w, col)) ++piv;
                        if (piv == w) { inv_ok = false; break; }
                        if (piv != col) { qcc_bv256_t t = M[col]; M[col] = M[piv]; M[piv] = t; }
                        for (uint32_t r = 0u; r < w; ++r) {
                            if (r != col && qcc_bv_get(M[r].w, col)) {
                                for (uint32_t ww = 0u; ww < QIR_WORDS; ++ww) M[r].w[ww] ^= M[col].w[ww];
                            }
                        }
                    }
                    free(M);
                    allow = inv_ok;
                }
            } else if (op->kind == QIR_OP_GF2M_MUL && st.variant == QIR_VAR_GF2M_CONST_LINF2 && op->num_ins == 2u) {
                uint32_t var_idx = st.aux_param & 1u;
                uint32_t cvid = op->ins[var_idx ^ 1u];
                uint16_t dm = op->u.gf2m.m;
                uint32_t dmask = (dm == 32u) ? UINT32_MAX : ((UINT32_C(1) << dm) - 1u);
                if ((m->vals[cvid].const_val[0] & dmask) != 0u) {
                    cand_idx = (uint8_t)var_idx;
                    cand_v   = op->ins[var_idx];
                    allow    = true;
                }
            } else if ((op->kind == QIR_OP_ADD || op->kind == QIR_OP_SUB || op->kind == QIR_OP_GF2M_ADD) &&
                       st.variant == QIR_VAR_INPLACE_ARITH && op->num_ins == 2u) {
                uint32_t tidx = st.aux_param & 1u;
                cand_idx = (uint8_t)tidx;
                cand_v   = op->ins[tidx];
                others[n_others++] = op->ins[tidx ^ 1u];
                allow = (op->ins[0] != op->ins[1]);
            } else if (op->kind == QIR_OP_CTRL_SEL && st.variant == QIR_VAR_CTRL_ADD_CUCCARO && st.num_chain_ops > 0u && op->num_ins == 3u) {
                uint32_t lhs_vid = st.aux_param & UINT32_C(0x7FFFFFFF);
                cand_idx = 1u;
                cand_v   = op->ins[1];
                others[n_others++] = op->ins[0];
                others[n_others++] = lhs_vid;
                allow = (cand_v != op->ins[0] && cand_v != lhs_vid);
            } else if (op->kind == QIR_OP_GF2M_ADD && st.variant == QIR_VAR_XOR_ACCUM && st.num_chain_ops > 0u && op->num_ins == 2u) {
                uint32_t acc_idx = st.aux_param & 3u;
                const qir_op_t *dop = &m->ops[st.chain_ops[0]];
                cand_idx = (uint8_t)acc_idx;
                cand_v   = op->ins[acc_idx];
                allow    = true;
                for (uint16_t i = 0u; i < dop->num_ins && n_others < 32u; ++i) {
                    if (dop->ins[i] == cand_v) allow = false;
                    others[n_others++] = dop->ins[i];
                }
            }

            if (allow && qir__can_exec_inplace(m, out_plan, op_skip, will_alias, is_in, is_pinned_out,
                                               is_direct_out, fwd_refs, uncomp_pin, back_vals, back_cnt,
                                               op, cand_v, others, n_others, eff_clean)) {
                st.inplace_on_in0 = true;
                st.inplace_in_idx = cand_idx;
                uint32_t dst = op->outs[0];
                out_plan->val_allocs[dst] = out_plan->val_allocs[cand_v];
                out_plan->val_allocs[dst].is_pinned_io = is_pinned_out[dst];
                out_plan->val_allocs[cand_v].is_alias = 1u;
            }
        }

        /* Allocate outputs (classical bits or clean quantum bits) */
        if (!st.inplace_on_in0) {
            for (uint8_t ko = 0u; ko < op->num_outs; ++ko) {
                uint32_t vid = op->outs[ko];
                uint16_t w   = m->vals[vid].width;
                qir_val_alloc_t *va = &out_plan->val_allocs[vid];
                va->is_pinned_io = is_pinned_out[vid];
                if (m->vals[vid].rep == QIR_REP_CLASSICAL) {
                    va->cwidth = w;
                    va->cbits  = w ? (uint32_t *)qcc_arena_xalloc(eg->arena, (size_t)w * sizeof(uint32_t), 4u) : NULL;
                    for (uint16_t b = 0u; b < w; ++b) va->cbits[b] = cbit++;
                } else {
                    va->width  = w;
                    va->qubits = w ? (uint32_t *)qcc_arena_xalloc(eg->arena, (size_t)w * sizeof(uint32_t), 4u) : NULL;
                    for (uint16_t b = 0u; b < w; ++b) va->qubits[b] = qir__pool_alloc_clean(pool, &peak_q);
                }
            }
        }

        /* Query step temp_anc and persistent clean low-bit aux_qubits */
        uint32_t need_temp = 0u, need_cb = 0u, need_aux = 0u;
        qir__query_step_temp(m, out_plan, &st, &need_temp, &need_cb, &need_aux);
        if (need_aux > 0u && op->num_outs > 0u && !st.inplace_on_in0) {
            qir_val_alloc_t *va0 = &out_plan->val_allocs[op->outs[0]];
            va0->aux_width  = (uint16_t)need_aux;
            va0->aux_qubits = (uint32_t *)qcc_arena_xalloc(eg->arena, (size_t)need_aux * sizeof(uint32_t), 4u);
            for (uint32_t a = 0u; a < need_aux; ++a) {
                va0->aux_qubits[a] = qir__pool_alloc_clean(pool, &peak_q);
            }
        }
        st.num_temp_anc = need_temp;
        st.cbit_base    = cbit;
        cbit           += need_cb;
        qir__alloc_step_temp(m, out_plan, &st, eg->arena, pool, &peak_q, &certified);

        if (peak_q > live_peak_q) {
            live_peak_q = peak_q;
            peak_op = oi;
        }
        if (ns >= max_steps) {
            max_steps *= 2u;
            out_plan->steps = (qir_sched_step_t *)realloc(out_plan->steps, max_steps * sizeof(qir_sched_step_t));
        }
        out_plan->steps[ns++] = st;
        computed[oi] = 1u;

        bool want_u = do_uncompute && (eff_clean || (uncomp_flags ? (uncomp_flags[oi] != 0u) : true));
        bool outs_uncomp_ok = (op->num_outs > 0u && m->vals[op->outs[0]].rep != QIR_REP_CLASSICAL &&
                               !st.inplace_on_in0 && qir__op_supports_uncompute(op, false));
        if (outs_uncomp_ok) {
            for (uint8_t ko = 0u; ko < op->num_outs; ++ko) {
                uint32_t ov = op->outs[ko];
                if (is_pinned_out[ov] || hmr_measured[ov]) { outs_uncomp_ok = false; break; }
            }
        }
        if (want_u && outs_uncomp_ok) {
            can_uncomp[oi] = 1u;
            for (uint16_t bi = 0u; bi < step_back_cnt[oi]; ++bi) {
                uncomp_pin[step_back_ins[(size_t)oi * 64u + bi]]++;
            }
        }

        /* Decrement forward refcounts on consumed backing operands & enqueue ready uncomputes */
        for (uint16_t bi = 0u; bi < step_back_cnt[oi]; ++bi) {
            uint32_t bv = step_back_ins[(size_t)oi * 64u + bi];
            if (fwd_refs[bv] > 0u) fwd_refs[bv]--;
            if (do_uncompute && fwd_refs[bv] == 0u && uncomp_pin[bv] == 0u) {
                qir__maybe_push_uncomp(m, computed, uncomputed, can_uncomp, fwd_refs, uncomp_pin,
                                       is_pinned_out, m->vals[bv].def_op, ready_heap, &rh_len, in_ready);
            }
        }
        if (do_uncompute) {
            qir__maybe_push_uncomp(m, computed, uncomputed, can_uncomp, fwd_refs, uncomp_pin,
                                   is_pinned_out, oi, ready_heap, &rh_len, in_ready);
        }

        /* HMR_MEAS resets non-alias non-pinned input qubits cleanly back to |0> */
        if (op->kind == QIR_OP_HMR_MEAS && op->pred_val == UINT32_MAX && op->num_ins > 0u) {
            uint32_t iv = op->ins[0];
            if (!will_alias[iv] && !is_in[iv] && !is_pinned_out[iv] &&
                out_plan->val_allocs[iv].aux_width == 0u &&
                fwd_refs[iv] == 0u && uncomp_pin[iv] == 0u) {
                const qir_val_alloc_t *iva = &out_plan->val_allocs[iv];
                for (uint16_t b = 0u; b < iva->width; ++b) {
                    qir__pool_free_clean(pool, iva->qubits[b]);
                }
            }
        }

        /* O(N log N) Cascading Reverse-Topological Clean Uncomputation via Ready Max-Heap */
        if (do_uncompute) {
            while (rh_len > 0u) {
                uint32_t uo = qir__heap_pop(ready_heap, &rh_len, in_ready);
                if (!computed[uo] || uncomputed[uo] || !can_uncomp[uo]) continue;
                const qir_op_t *dop = &m->ops[uo];
                bool all_outs_free = true;
                for (uint8_t ko = 0u; ko < dop->num_outs; ++ko) {
                    uint32_t ov = dop->outs[ko];
                    if (fwd_refs[ov] > 0u || uncomp_pin[ov] > 0u || is_pinned_out[ov]) {
                        all_outs_free = false;
                        break;
                    }
                }
                if (!all_outs_free) continue;

                const qir_enode_t *den = &eg->enodes[out_plan->chosen_enodes[dop->outs[0]]];
                qir_sched_step_t ust;
                memset(&ust, 0, sizeof(ust));
                ust.kind             = QIR_STEP_UNCOMPUTE;
                ust.op_id            = uo;
                ust.variant          = den->variant;
                ust.aux_param        = den->aux_param;
                ust.approx_keep_bits = (approx_widths && (dop->flags & QIR_FLAG_APPROX)) ? approx_widths[uo] : m->vals[dop->outs[0]].width;
                ust.fused_linf2      = den->fused_linf2;
                ust.gauge_lut        = den->gauge_lut;
                ust.phase_split      = den->phase_split;
                ust.chain_ops        = den->chain_ops;
                ust.num_chain_ops    = den->num_chain_ops;
                ust.pred_cbit        = UINT32_MAX;
                ust.inplace_in_idx   = 0xFFu;
                if (dop->pred_val != UINT32_MAX && dop->pred_val < m->num_vals && out_plan->val_allocs[dop->pred_val].cbits) {
                    ust.pred_cbit = out_plan->val_allocs[dop->pred_val].cbits[0];
                }

                uint32_t u_temp = 0u, u_cb = 0u, u_aux = 0u;
                qir__query_step_temp(m, out_plan, &ust, &u_temp, &u_cb, &u_aux);
                ust.num_temp_anc = u_temp;
                ust.cbit_base    = cbit;
                cbit            += u_cb;
                qir__alloc_step_temp(m, out_plan, &ust, eg->arena, pool, &peak_q, &certified);
                if (ns >= max_steps) {
                    max_steps *= 2u;
                    out_plan->steps = (qir_sched_step_t *)realloc(out_plan->steps, max_steps * sizeof(qir_sched_step_t));
                }
                out_plan->steps[ns++] = ust;

                for (uint8_t ko = 0u; ko < dop->num_outs; ++ko) {
                    const qir_val_alloc_t *uva = &out_plan->val_allocs[dop->outs[ko]];
                    for (uint16_t b = 0u; b < uva->width; ++b) qir__pool_free_clean(pool, uva->qubits[b]);
                    for (uint16_t b = 0u; b < uva->aux_width; ++b) qir__pool_free_clean(pool, uva->aux_qubits[b]);
                }
                uncomputed[uo] = 1u;
                for (uint16_t bi = 0u; bi < step_back_cnt[uo]; ++bi) {
                    uint32_t bv = step_back_ins[(size_t)uo * 64u + bi];
                    if (uncomp_pin[bv] > 0u) uncomp_pin[bv]--;
                    if (fwd_refs[bv] == 0u && uncomp_pin[bv] == 0u) {
                        qir__maybe_push_uncomp(m, computed, uncomputed, can_uncomp, fwd_refs, uncomp_pin,
                                               is_pinned_out, m->vals[bv].def_op, ready_heap, &rh_len, in_ready);
                    }
                }
            }
        }
    }

    /* Certify Clean Scratch: check whether any non-I/O qubit remains busy */
    uint64_t *io_busy = (uint64_t *)calloc(QIR_POOL_WORDS, sizeof(uint64_t));
    if (io_busy) {
        for (uint32_t i = 0u; i < m->num_inputs; ++i) {
            const qir_val_alloc_t *va = &out_plan->val_allocs[m->inputs[i]];
            for (uint16_t b = 0u; b < va->width; ++b) {
                uint32_t q = va->qubits[b];
                if (q < QIR_POOL_MAX_Q) io_busy[q >> 6u] |= (UINT64_C(1) << (q & 63u));
            }
        }
        for (uint32_t i = 0u; i < m->num_outputs; ++i) {
            const qir_val_alloc_t *va = &out_plan->val_allocs[m->outputs[i]];
            for (uint16_t b = 0u; b < va->width; ++b) {
                uint32_t q = va->qubits[b];
                if (q < QIR_POOL_MAX_Q) io_busy[q >> 6u] |= (UINT64_C(1) << (q & 63u));
            }
        }
    }
    bool clean_ok = (io_busy != NULL);
    for (uint32_t q = 0u; q < peak_q && clean_ok; ++q) {
        if (qir__pool_is_busy(pool, q) && !((io_busy[q >> 6u] >> (q & 63u)) & 1u)) {
            clean_ok = false;
            break;
        }
    }
    free(io_busy);
    free(pool);

    out_plan->num_steps               = ns;
    out_plan->disjointness_certified  = certified;
    out_plan->clean_scratch_certified = clean_ok;
    out_plan->max_live_frustration_op = peak_op;
    if (out_peak_op_at_max) *out_peak_op_at_max = peak_op;

    free(ready_heap); free(in_ready);
    free(op_skip); free(will_alias); free(is_direct_out); free(is_in); free(is_pinned_out);
    free(hmr_measured); free(can_uncomp); free(computed); free(uncomputed);
    free(fwd_refs); free(uncomp_pin); free(back_vals); free(back_cnt);
    free(step_back_ins); free(step_back_cnt);

    qcc_kmb_stats_t est;
    if (!qir_codegen_or_count(out_plan, NULL, &est)) {
        return false;
    }
    out_plan->est_qubits   = est.num_qubits > peak_q ? est.num_qubits : peak_q;
    out_plan->est_bits     = est.num_bits;
    out_plan->est_toffoli  = est.toffoli;
    out_plan->est_z_pow    = est.z_pow;
    out_plan->est_t_like   = est.t_like;
    out_plan->est_hmr      = est.hmr;
    out_plan->est_cnot     = est.cnot;
    return true;
}

static inline uint64_t qir__score_plan(const qir_phys_plan_t *p, const qcc_opt_cfg_t *cfg) {
    uint64_t q   = p->est_qubits;
    uint64_t ccx = p->est_toffoli;
    uint64_t zp  = p->est_z_pow;
    uint64_t tl  = p->est_t_like;
    uint64_t pen = 0u;
    if ((cfg->require_clean_scratch || p->mod->require_clean_scratch) && !p->clean_scratch_certified) {
        pen += UINT64_C(1000000000000);
    }
    if (cfg->max_qubits > 0u && q > cfg->max_qubits) {
        pen += UINT64_C(1000000000) * (uint64_t)(q - cfg->max_qubits);
    }
    uint64_t core = 0u;
    switch (cfg->cost_mode) {
        case QCC_COST_QC:  core = q * (ccx + zp + 1u); break;
        case QCC_COST_Q2T: core = q * q * (4u * ccx + tl + 1u); break;
        case QCC_COST_QT:  core = q * (4u * ccx + zp + 1u); break;
        default:           core = 1000u * ccx + 250u * zp + 50u * q + 10u * p->est_hmr + p->est_cnot; break;
    }
    return pen + core;
}

static inline bool qir__eval_approx_candidate(
    const qir_egraph_t *eg,
    const uint32_t *base_choices,
    const uint8_t *uncomp_flags,
    uint16_t *approx_widths,
    uint32_t target_oi,
    uint16_t cand_w,
    const qcc_opt_cfg_t *cfg,
    bool eff_clean,
    qir_phys_plan_t *best_plan,
    uint64_t *best_score,
    uint16_t *best_w,
    uint64_t *out_sc
) {
    approx_widths[target_oi] = cand_w;
    qir_phys_plan_t tp;
    uint32_t pk = 0u;
    if (!qir_regalloc_pack(eg, base_choices, uncomp_flags, approx_widths, true, eff_clean, &tp, &pk)) {
        *out_sc = UINT64_MAX;
        return false;
    }
    uint64_t sc = qir__score_plan(&tp, cfg);
    *out_sc = sc;
    if (sc <= *best_score) {
        qir_phys_plan_free(best_plan);
        *best_plan  = tp;
        *best_score = sc;
        *best_w     = cand_w;
    } else {
        qir_phys_plan_free(&tp);
    }
    return true;
}

/* Main QIR Optimizer: E-Graph Saturation + Impact-Ranked SEARU Genome Search + FLOW 2-Jet Relaxation */
static inline bool qir_optimize_module(
    const qir_module_t *mod,
    qcc_arena_t *arena,
    const qcc_opt_cfg_t *cfg_in,
    qir_egraph_t *out_eg_opt,
    qir_phys_plan_t *out_plan
) {
    if (!mod || !arena || !out_plan) return false;
    qcc_opt_cfg_t cfg = cfg_in ? *cfg_in : qcc_opt_cfg_default();
    bool eff_clean = cfg.require_clean_scratch || mod->require_clean_scratch;

    qir_egraph_t eg_local;
    qir_egraph_t *eg = out_eg_opt ? out_eg_opt : &eg_local;
    if (!qir_egraph_build_and_saturate(eg, arena, mod)) return false;

    size_t nc = eg->num_classes ? eg->num_classes : 1u;
    size_t no = mod->num_ops ? mod->num_ops : 1u;

    /* Impact-ranked selection of multi-choice E-Classes mapped to 1..3 bits in SEARU's 256-bit canvas */
    uint32_t *all_multi   = (uint32_t *)qcc_xmalloc(nc * sizeof(uint32_t));
    uint64_t *multi_score = (uint64_t *)qcc_xmalloc(nc * sizeof(uint64_t));
    uint32_t n_all_multi  = 0u;
    for (uint32_t c = 0u; c < eg->num_classes; ++c) {
        if (eg->classes[c].num_enodes <= 1u) continue;
        uint64_t sc = (uint64_t)(eg->classes[c].width ? eg->classes[c].width : 1u) * (eg->classes[c].use_count + 1u) + 10u;
        uint32_t dop = mod->vals[c].def_op;
        if (dop != UINT32_MAX) {
            qir_op_kind_t ok = mod->ops[dop].kind;
            if (ok == QIR_OP_LUT_ROM) sc += (uint64_t)mod->ops[dop].u.lut.num_entries * 20u + 2000u;
            else if (ok == QIR_OP_PHASE_POLY) sc += 5000u;
            else if (ok == QIR_OP_GF2M_ADD) sc += 2500u;
            else if (ok == QIR_OP_MUL_FXP || ok == QIR_OP_SQR_FXP) sc += 1500u;
            else if (ok == QIR_OP_GF2M_INV || ok == QIR_OP_GF2M_MUL || ok == QIR_OP_CTRL_SEL) sc += 1200u;
            else if (ok == QIR_OP_ADD || ok == QIR_OP_SUB || ok == QIR_OP_LINF2) sc += 400u;
        }
        all_multi[n_all_multi]   = c;
        multi_score[n_all_multi] = sc;
        n_all_multi++;
    }
    for (uint32_t i = 0u; i < n_all_multi; ++i) {
        for (uint32_t j = i + 1u; j < n_all_multi; ++j) {
            if (multi_score[j] > multi_score[i]) {
                uint64_t ts = multi_score[i]; multi_score[i] = multi_score[j]; multi_score[j] = ts;
                uint32_t tc = all_multi[i];   all_multi[i]   = all_multi[j];   all_multi[j]   = tc;
            }
        }
    }

    uint32_t multi_cls[192];
    uint8_t  cls_bits[192];
    uint16_t cls_bit_off[192];
    uint32_t n_multi = 0u;
    uint32_t n_egraph_vars = 0u;
    for (uint32_t i = 0u; i < n_all_multi && n_multi < 192u; ++i) {
        uint32_t cid = all_multi[i];
        uint32_t ne = eg->classes[cid].num_enodes;
        uint8_t b = (ne > 16u) ? 5u : ((ne > 8u) ? 4u : ((ne > 4u) ? 3u : ((ne > 2u) ? 2u : 1u)));
        if (n_egraph_vars + b > 192u) break;
        multi_cls[n_multi]   = cid;
        cls_bits[n_multi]    = b;
        cls_bit_off[n_multi] = (uint16_t)n_egraph_vars;
        n_egraph_vars       += b;
        n_multi++;
    }
    free(multi_score);

    /* Rank uncompute-capable ops by output width & cost so SEARU targets high-value uncomputation first */
    uint32_t *all_uncomp = (uint32_t *)qcc_xmalloc(no * sizeof(uint32_t));
    uint64_t *uncomp_sc  = (uint64_t *)qcc_xmalloc(no * sizeof(uint64_t));
    uint32_t n_all_unc   = 0u;
    for (uint32_t oi = 0u; oi < mod->num_ops; ++oi) {
        const qir_op_t *op = &mod->ops[oi];
        if (op->num_outs > 0u && mod->vals[op->outs[0]].rep != QIR_REP_CLASSICAL &&
            qir__op_supports_uncompute(op, false)) {
            uint16_t w = mod->vals[op->outs[0]].width;
            uint64_t sc = (uint64_t)w * 10u + (op->kind == QIR_OP_LUT_ROM ? 1000u : 100u);
            all_uncomp[n_all_unc] = oi;
            uncomp_sc[n_all_unc]  = sc;
            n_all_unc++;
        }
    }
    for (uint32_t i = 0u; i < n_all_unc; ++i) {
        for (uint32_t j = i + 1u; j < n_all_unc; ++j) {
            if (uncomp_sc[j] > uncomp_sc[i]) {
                uint64_t ts = uncomp_sc[i]; uncomp_sc[i]  = uncomp_sc[j];  uncomp_sc[j]  = ts;
                uint32_t to = all_uncomp[i]; all_uncomp[i] = all_uncomp[j]; all_uncomp[j] = to;
            }
        }
    }
    uint32_t max_unc_slots = (256u > n_egraph_vars) ? (256u - n_egraph_vars) : 0u;
    uint32_t n_uncomp_vars = (n_all_unc > max_unc_slots) ? max_unc_slots : n_all_unc;
    uint32_t uncomp_ops[256];
    for (uint32_t i = 0u; i < n_uncomp_vars; ++i) uncomp_ops[i] = all_uncomp[i];
    free(all_uncomp);
    free(uncomp_sc);

    /* Backward reachability from PHASE_POLY / PHASE_DOT so phase-feeding approx ops remain bit-exact */
    uint8_t *feeds_phase_val = (uint8_t *)qcc_xcalloc(eg->num_classes ? eg->num_classes : 1u, 1u);
    uint8_t *feeds_phase_op  = (uint8_t *)qcc_xcalloc(no, 1u);
    for (uint32_t oi = 0u; oi < mod->num_ops; ++oi) {
        const qir_op_t *op = &mod->ops[oi];
        if (op->kind == QIR_OP_PHASE_POLY || op->kind == QIR_OP_PHASE_DOT) {
            feeds_phase_op[oi] = 1u;
            for (uint16_t ii = 0u; ii < op->num_ins; ++ii) feeds_phase_val[op->ins[ii]] = 1u;
        }
    }
    for (int32_t oi_s = (int32_t)mod->num_ops - 1; oi_s >= 0; --oi_s) {
        uint32_t oi = (uint32_t)oi_s;
        const qir_op_t *op = &mod->ops[oi];
        bool out_feeds = false;
        for (uint8_t ko = 0u; ko < op->num_outs; ++ko) {
            if (feeds_phase_val[op->outs[ko]]) { out_feeds = true; break; }
        }
        if (out_feeds) {
            feeds_phase_op[oi] = 1u;
            for (uint16_t ii = 0u; ii < op->num_ins; ++ii) feeds_phase_val[op->ins[ii]] = 1u;
        }
    }

    uint16_t *approx_widths = (uint16_t *)qcc_xcalloc(no, sizeof(uint16_t));
    uint8_t  *uncomp_flags  = (uint8_t *)qcc_xmalloc(no);
    uint32_t *base_choices  = (uint32_t *)qcc_xcalloc(nc, sizeof(uint32_t));
    for (uint32_t oi = 0u; oi < mod->num_ops; ++oi) {
        const qir_op_t *op = &mod->ops[oi];
        uint16_t out_w = (op->num_outs > 0u) ? mod->vals[op->outs[0]].width : 0u;
        if (op->kind == QIR_OP_MUL_FXP) {
            int32_t sh = (int32_t)op->u.mul.frac_a + (int32_t)op->u.mul.frac_b - (int32_t)op->u.mul.frac_out;
            approx_widths[oi] = (uint16_t)(out_w + (sh > 0 ? (uint32_t)sh : 0u));
        } else if (op->kind == QIR_OP_SQR_FXP) {
            int32_t sh = 2 * (int32_t)op->u.mul.frac_a - (int32_t)op->u.mul.frac_out;
            approx_widths[oi] = (uint16_t)(out_w + (sh > 0 ? (uint32_t)sh : 0u));
        } else {
            approx_widths[oi] = out_w;
        }
        uncomp_flags[oi] = (eff_clean || cfg.aggressive_uncompute) ? 1u : 0u;
    }
    for (uint32_t c = 0u; c < eg->num_classes; ++c) {
        base_choices[c] = eg->classes[c].num_enodes ? eg->classes[c].enodes[0] : 0u;
    }
    if (n_all_multi > n_multi || mod->num_ops > 200u) {
        static const qir_variant_kind_t tail_pref[] = {
            QIR_VAR_GF2M_TOWER_INV,
            QIR_VAR_GF2M_RANK27,
            QIR_VAR_GF2M_CHUDNOVSKY_R42,
            QIR_VAR_GF2M_RANK15,
            QIR_VAR_GF2M_CONST_LINF2,
            QIR_VAR_LINF2_GAUSS,
            QIR_VAR_CTRL_ADD_CUCCARO,
            QIR_VAR_XOR_ACCUM,
            QIR_VAR_CSA_TREE,
            QIR_VAR_INPLACE_ARITH,
            QIR_VAR_HMR_PARITY_LUT
        };
        uint32_t seed_start = (mod->num_ops > 200u) ? 0u : n_multi;
        for (uint32_t i = seed_start; i < n_all_multi; ++i) {
            uint32_t cid = all_multi[i];
            const qir_eclass_t *ec = &eg->classes[cid];
            uint32_t dop = mod->vals[cid].def_op;
            qir_op_kind_t ok = (dop != UINT32_MAX) ? mod->ops[dop].kind : QIR_OP_COUNT;
            bool picked = false;
            for (size_t pi = 0u; pi < sizeof(tail_pref) / sizeof(tail_pref[0]) && !picked; ++pi) {
                if (tail_pref[pi] == QIR_VAR_HMR_PARITY_LUT && ok != QIR_OP_GF2M_MUL) continue;
                for (uint32_t ei = 0u; ei < ec->num_enodes; ++ei) {
                    uint32_t eid = ec->enodes[ei];
                    const qir_enode_t *en = &eg->enodes[eid];
                    if (en->variant != tail_pref[pi]) continue;
                    if (en->rule_tag != 0u && !((cfg.enabled_rules >> en->rule_tag) & 1u)) continue;
                    base_choices[cid] = eid;
                    picked = true;
                    break;
                }
            }
        }
    }
    free(all_multi);

    /* Phase 1: Baseline SEED plan evaluation */
    qir_phys_plan_t best_plan;
    uint32_t peak_op = 0u;
    if (!qir_regalloc_pack(eg, base_choices, uncomp_flags, approx_widths, true, eff_clean, &best_plan, &peak_op)) {
        free(feeds_phase_val); free(feeds_phase_op);
        free(approx_widths); free(uncomp_flags); free(base_choices);
        if (!out_eg_opt) qir_egraph_free(eg);
        return false;
    }
    uint64_t best_score = qir__score_plan(&best_plan, &cfg);

    /* Phase 1b: Single-class coordinate sweep so SEARU starts from exact rule-enabled best */
    uint32_t sweep_lim = (mod->num_ops > 200u && n_multi > 24u) ? 24u : n_multi;
    for (uint32_t mi = 0u; mi < sweep_lim; ++mi) {
        uint32_t cid = multi_cls[mi];
        const qir_eclass_t *ec = &eg->classes[cid];
        uint32_t best_eid = base_choices[cid];
        for (uint32_t ei = 0u; ei < ec->num_enodes; ++ei) {
            uint32_t eid = ec->enodes[ei];
            uint8_t rt = eg->enodes[eid].rule_tag;
            if (rt != 0u && !((cfg.enabled_rules >> rt) & 1u)) continue;
            base_choices[cid] = eid;
            qir_phys_plan_t cplan;
            uint32_t pk_op = 0u;
            if (qir_regalloc_pack(eg, base_choices, uncomp_flags, approx_widths, true, eff_clean, &cplan, &pk_op)) {
                uint64_t sc = qir__score_plan(&cplan, &cfg);
                if (sc < best_score) {
                    qir_phys_plan_free(&best_plan);
                    best_plan  = cplan;
                    best_score = sc;
                    best_eid   = eid;
                } else {
                    qir_phys_plan_free(&cplan);
                }
            }
        }
        base_choices[cid] = best_eid;
    }

    /* Phase 2: Genuine FLOW 2-Jet Continuous Precision Relaxation for QIR_FLAG_APPROX ops */
    uint32_t flow_steps = cfg.flow_iters ? cfg.flow_iters : 16u;
    bool has_approx = false;
    for (uint32_t oi = 0u; oi < mod->num_ops; ++oi) {
        const qir_op_t *op = &mod->ops[oi];
        bool is_cg = ((op->kind == QIR_OP_MUL_FXP || op->kind == QIR_OP_SQR_FXP) &&
                      (op->flags & QIR_FLAG_APPROX) && op->u.mul.keep_exact == 2u);
        if ((op->flags & QIR_FLAG_APPROX) && (!feeds_phase_op[oi] || is_cg) &&
            op->approx_tol_q16 >= (UINT32_C(1) << 16u) && op->num_outs > 0u) {
            has_approx = true;
            break;
        }
    }
    if (has_approx) {
        flow_cell *flow_arena = (flow_cell *)calloc(32u, sizeof(flow_cell));
        for (uint32_t oi = 0u; oi < mod->num_ops; ++oi) {
            const qir_op_t *op = &mod->ops[oi];
            bool is_cg = ((op->kind == QIR_OP_MUL_FXP || op->kind == QIR_OP_SQR_FXP) &&
                          (op->flags & QIR_FLAG_APPROX) && op->u.mul.keep_exact == 2u);
            if (!(op->flags & QIR_FLAG_APPROX) || (feeds_phase_op[oi] && !is_cg) ||
                op->approx_tol_q16 < (UINT32_C(1) << 16u) || op->num_outs == 0u) {
                continue;
            }
            uint16_t full_w  = approx_widths[oi];
            uint16_t out_w   = mod->vals[op->outs[0]].width;
            uint32_t tol_lsb = op->approx_tol_q16 >> 16u;
            uint16_t min_keep = full_w;
            if (is_cg) {
                min_keep = (out_w + tol_lsb < full_w) ? (uint16_t)(out_w + tol_lsb) : full_w;
            } else if (op->kind == QIR_OP_MUL_FXP || op->kind == QIR_OP_SQR_FXP) {
                min_keep = (tol_lsb >= 2u) ? out_w : full_w;
            } else {
                uint32_t drop_bits = 0u;
                while ((UINT32_C(1) << (drop_bits + 1u)) <= tol_lsb && drop_bits < 30u) drop_bits++;
                min_keep = (out_w > drop_bits) ? (uint16_t)(out_w - drop_bits) : 2u;
                if (min_keep < 2u && out_w >= 2u) min_keep = 2u;
            }
            if (min_keep >= full_w) continue;

            uint64_t bound[8] = {0};
            bound[3] = ((uint64_t)(uint32_t)((int32_t)min_keep * FLOW_Q16_ONE) << 32U) |
                       (uint64_t)(uint32_t)((int32_t)full_w * FLOW_Q16_ONE);
            bound[4] = ((uint64_t)(uint32_t)(2 * FLOW_Q16_ONE) << 32U) |
                       (uint64_t)(uint32_t)FLOW_Q16_ONE;
            flow_arena_publish_bound(flow_arena, 0u, 31u, 1u, bound);

            uint16_t w_cur  = full_w;
            uint16_t best_w = full_w;
            for (uint32_t it = 0u; it < flow_steps && w_cur > min_keep; ++it) {
                uint16_t w1 = (w_cur > min_keep) ? (uint16_t)(w_cur - 1u) : min_keep;
                uint16_t w2 = (w_cur >= (uint16_t)(min_keep + 2u)) ? (uint16_t)(w_cur - 2u) : min_keep;
                uint64_t s0 = 0u, s1 = 0u, s2 = 0u;
                qir__eval_approx_candidate(eg, base_choices, uncomp_flags, approx_widths, oi, w_cur,
                                           &cfg, eff_clean, &best_plan, &best_score, &best_w, &s0);
                qir__eval_approx_candidate(eg, base_choices, uncomp_flags, approx_widths, oi, w1,
                                           &cfg, eff_clean, &best_plan, &best_score, &best_w, &s1);
                qir__eval_approx_candidate(eg, base_choices, uncomp_flags, approx_widths, oi, w2,
                                           &cfg, eff_clean, &best_plan, &best_score, &best_w, &s2);

                int64_t d1 = (s0 >= s1) ? (int64_t)(s0 - s1) : -(int64_t)(s1 - s0);
                int64_t d2 = (s1 >= s2) ? (int64_t)(s1 - s2) : -(int64_t)(s2 - s1);
                int64_t d2_fd = d1 - d2;
                int64_t max_d = (d1 < 0 ? -d1 : d1);
                if ((d2_fd < 0 ? -d2_fd : d2_fd) > max_d) max_d = (d2_fd < 0 ? -d2_fd : d2_fd);
                uint32_t sh = 0u;
                while ((max_d >> sh) > (FLOW_MAX_Q16 >> 4) && sh < 48u) sh++;

                int64_t j1 = d1 >> sh;
                int64_t j2 = d2_fd >> sh;
                int64_t j0 = (j1 > 0 ? j1 : 0) + (int64_t)(w_cur - min_keep);
                if (j1 < -FLOW_MAX_J1_Q16) j1 = -FLOW_MAX_J1_Q16;
                if (j1 >  FLOW_MAX_J1_Q16) j1 =  FLOW_MAX_J1_Q16;
                if (j2 < -FLOW_MAX_Q16)    j2 = -FLOW_MAX_Q16;
                if (j2 >  FLOW_MAX_Q16)    j2 =  FLOW_MAX_Q16;

                int64_t pressure = (cfg.max_qubits > 0u && best_plan.est_qubits > cfg.max_qubits)
                    ? ((int64_t)(best_plan.est_qubits - cfg.max_qubits) * FLOW_Q16_ONE)
                    : 0;
                uint64_t st_in[8] = {
                    (uint64_t)j0, (uint64_t)j1, (uint64_t)j2,
                    (uint64_t)(uint32_t)((int32_t)w_cur * FLOW_Q16_ONE),
                    0u, 0u, 0u, 0u
                };
                uint64_t st_out[8] = {0};
                flow_step_atlas(st_out, st_in, (const flow_cell *)flow_arena, 31u, pressure, 0u);
                int64_t x_next   = (int64_t)(int32_t)(uint32_t)st_out[3];
                int64_t w_next_s = (x_next + (FLOW_Q16_ONE / 2)) >> FLOW_Q16_SHIFT;
                if (w_next_s < (int64_t)min_keep) w_next_s = (int64_t)min_keep;
                if (w_next_s > (int64_t)w_cur)    w_next_s = (int64_t)w_cur;

                if ((uint16_t)w_next_s < w_cur) {
                    w_cur = (uint16_t)w_next_s;
                } else if (s1 <= best_score && w_cur > min_keep) {
                    w_cur = w1;
                } else {
                    break;
                }
            }
            approx_widths[oi] = best_w;
        }
        free(flow_arena);
    }

    /* Phase 3: Symbiotic SEARU Inline Discrete-Structural Optimizer */
    uint32_t num_vars = n_egraph_vars + n_uncomp_vars;
    if (num_vars == 0u) num_vars = 1u;
    if (num_vars > 256u) num_vars = 256u;

    uint64_t init_c[4] = {0}, pin_c[4] = {0};
    for (uint32_t ui = 0u; ui < n_uncomp_vars; ++ui) {
        uint32_t oi = uncomp_ops[ui];
        if (eff_clean) {
            qcc_bv_set(init_c, n_egraph_vars + ui);
            qcc_bv_set(pin_c, n_egraph_vars + ui);
        } else if (uncomp_flags[oi]) {
            qcc_bv_set(init_c, n_egraph_vars + ui);
        }
    }
    for (uint32_t mi = 0u; mi < n_multi; ++mi) {
        uint32_t off = cls_bit_off[mi];
        const qir_eclass_t *ec = &eg->classes[multi_cls[mi]];
        uint32_t code = 0u;
        for (uint32_t ei = 0u; ei < ec->num_enodes; ++ei) {
            if (ec->enodes[ei] == base_choices[multi_cls[mi]]) { code = ei; break; }
        }
        for (uint8_t bi = 0u; bi < cls_bits[mi]; ++bi) {
            if ((code >> bi) & 1u) qcc_bv_set(init_c, off + bi);
        }
    }

    searu_hyperparams hp = SEARU_ACTIVE_OPT_HP;
    searu_result sres;
    uint32_t *cand_choices = (uint32_t *)qcc_xmalloc(nc * sizeof(uint32_t));
    uint8_t  *cand_uncomp  = (uint8_t *)qcc_xmalloc(no);

    SEARU_SOLVE_INLINE_HP(
        &sres, num_vars, init_c, pin_c, cfg.seed,
        (cfg.searu_iters ? cfg.searu_iters : 128u), &hp,
        cand_v, frust_v, energy_v,
        {
            memcpy(cand_choices, base_choices, nc * sizeof(uint32_t));
            memcpy(cand_uncomp, uncomp_flags, no);
            for (uint32_t mi = 0u; mi < n_multi; ++mi) {
                uint32_t cid = multi_cls[mi];
                const qir_eclass_t *ec = &eg->classes[cid];
                uint32_t off = cls_bit_off[mi];
                uint32_t code = 0u;
                for (uint8_t bi = 0u; bi < cls_bits[mi]; ++bi) {
                    if (qcc_bv_get(cand_v, off + bi)) code |= (UINT32_C(1) << bi);
                }
                uint32_t idx = code % ec->num_enodes;
                uint32_t eid = ec->enodes[idx];
                uint8_t rt = eg->enodes[eid].rule_tag;
                if (rt == 0u || ((cfg.enabled_rules >> rt) & 1u)) cand_choices[cid] = eid;
            }
            for (uint32_t ui = 0u; ui < n_uncomp_vars; ++ui) {
                cand_uncomp[uncomp_ops[ui]] = (eff_clean || qcc_bv_get(cand_v, n_egraph_vars + ui)) ? 1u : 0u;
            }

            qir_phys_plan_t cplan;
            uint32_t pk_op = 0u;
            if (!qir_regalloc_pack(eg, cand_choices, cand_uncomp, approx_widths, true, eff_clean, &cplan, &pk_op)) {
                energy_v = 1e15;
            } else {
                uint64_t sc = qir__score_plan(&cplan, &cfg);
                energy_v = (double)sc;
                /* Dynamic Peak-Live Frustration Mask guides SEARU to uncompute / ColSplit overflow ops! */
                if (cfg.max_qubits > 0u && cplan.est_qubits > cfg.max_qubits) {
                    for (uint32_t ui = 0u; ui < n_uncomp_vars; ++ui) {
                        if (uncomp_ops[ui] == pk_op) {
                            qcc_bv_set(frust_v, n_egraph_vars + ui);
                            break;
                        }
                    }
                    const qir_op_t *pop = &mod->ops[pk_op];
                    if (pop->num_outs > 0u) {
                        for (uint32_t mi = 0u; mi < n_multi; ++mi) {
                            if (multi_cls[mi] == pop->outs[0]) {
                                for (uint8_t bi = 0u; bi < cls_bits[mi]; ++bi) {
                                    qcc_bv_set(frust_v, cls_bit_off[mi] + bi);
                                }
                            }
                        }
                    }
                }
                if (sc < best_score) {
                    qir_phys_plan_free(&best_plan);
                    best_plan  = cplan;
                    best_score = sc;
                } else {
                    qir_phys_plan_free(&cplan);
                }
            }
        }
    );

    best_plan.emit_clean_reg = cfg.emit_clean_reg;
    *out_plan = best_plan;
    free(feeds_phase_val);
    free(feeds_phase_op);
    free(cand_choices);
    free(cand_uncomp);
    free(approx_widths);
    free(uncomp_flags);
    free(base_choices);
    if (!out_eg_opt) qir_egraph_free(eg);
    return true;
}

#ifdef __cplusplus
}
#endif

#endif /* QIR_OPT_H */
