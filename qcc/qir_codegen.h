/*
 * qir_codegen.h — Zero-Discrepancy Physical Lowering (QIR Plan -> .kmb v2)
 * Part of QCC (Quantum Circuit Compiler) for FLOW v7.0.1.
 * Pure C17, zero external dependencies, zero indirect function calls.
 */
#ifndef QIR_CODEGEN_H
#define QIR_CODEGEN_H

#include "qir_kernels.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    QCC_COST_DEFAULT = 0, /* Balanced Toffoli + Z_POW + Qubit Tie-Break */
    QCC_COST_QC      = 1, /* Minimize Q * (CCX + Z_POW) */
    QCC_COST_Q2T     = 2, /* Minimize Q^2 * T-depth / Toffoli */
    QCC_COST_QT      = 3  /* Minimize Q * Toffoli */
} qcc_cost_mode_t;

typedef enum {
    QIR_STEP_COMPUTE   = 0,
    QIR_STEP_UNCOMPUTE = 1
} qir_step_kind_t;

typedef struct {
    qir_step_kind_t    kind;
    uint32_t           op_id;
    qir_variant_kind_t variant;
    uint32_t           aux_param;
    uint16_t           approx_keep_bits;
    const qir_linf2_payload_t *fused_linf2;
    const qir_gauge_lut_t     *gauge_lut;
    const qir_phase_split_t   *phase_split;
    const uint32_t            *chain_ops;
    uint16_t                   num_chain_ops;
    /* Regalloc assignments (populated by qir_regalloc_pack) */
    uint32_t          *temp_anc;
    uint32_t           num_temp_anc;
    uint32_t           cbit_base;
    uint32_t           pred_cbit; /* UINT32_MAX if unconditional */
    bool               inplace_on_in0;
    uint8_t            inplace_in_idx; /* 0xFF when out-of-place, else input index */
} qir_sched_step_t;

typedef struct {
    uint32_t *qubits;       /* [width] qubit indices assigned to SSA value */
    uint16_t  width;
    uint32_t *aux_qubits;   /* [aux_width] clean low-bit scratch for reversible MUL/SQR */
    uint16_t  aux_width;
    uint32_t *cbits;        /* [cwidth] classical bit indices when QIR_REP_CLASSICAL */
    uint16_t  cwidth;
    uint8_t   is_alias;     /* 1 for zero-cost SLICE / CONCAT / pass-through */
    uint8_t   is_pinned_io; /* 1 for module input or output */
    uint32_t  dirty_tail_start;
    uint32_t  dirty_tail_len;
} qir_val_alloc_t;

typedef struct {
    qcc_arena_t          *arena;
    const qir_module_t   *mod;
    const qir_egraph_t   *eg;
    qir_sched_step_t     *steps;
    uint32_t              num_steps;
    qir_val_alloc_t      *val_allocs; /* [mod->num_vals] */
    uint32_t             *chosen_enodes; /* [mod->num_vals] */

    /* Cached algebraic tower plans (shared across GF ops) */
    bool                       has_gf12_r42;
    uint32_t                   gf12_r42_poly;
    qir_gf12_rank42_t          gf12_r42;
    bool                       has_gf12_tower;
    uint32_t                   gf12_tower_poly;
    qir_gf12_tower_inv_plan_t  gf12_tower;
    bool                       has_gf256_r27;
    uint32_t                   gf256_r27_poly;
    qir_gf256_rank27_t         gf256_r27;
    bool                       has_gf256_tower;
    uint32_t                   gf256_tower_poly;
    qir_gf256_tower_inv_plan_t gf256_tower;
    bool                       has_gf64_r15;
    uint32_t                   gf64_r15_poly;
    qir_gf64_rank15_t          gf64_r15;

    /* Hard-verified physical resource certificates (Zero-Discrepancy Invariant) */
    uint32_t              est_qubits;
    uint32_t              est_bits;
    uint64_t              est_toffoli;
    uint64_t              est_z_pow;
    uint64_t              est_t_like;
    uint64_t              est_hmr;
    uint64_t              est_cnot;
    uint64_t              est_reaction_depth;
    bool                  disjointness_certified;
    bool                  clean_scratch_certified;
    bool                  emit_clean_reg;
    uint32_t              max_live_frustration_op;
} qir_phys_plan_t;

static inline void qir_phys_plan_free(qir_phys_plan_t *p) {
    if (!p) return;
    free(p->steps);
    free(p->val_allocs);
    free(p->chosen_enodes);
    memset(p, 0, sizeof(*p));
}

/* Helper: CLZ_NORM kernel with optional KB emission */
static inline void qk__toggle_all_zero(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *bits,
    uint32_t len,
    uint32_t dst,
    const uint32_t *cand,
    uint32_t *sb
) {
    for (uint32_t i = 0u; i < len; ++i) qk_x(kb, cost, bits[i]);
    if (len == 1u) {
        qk_cx(kb, cost, bits[0], dst);
    } else {
        qk_and_compute(kb, cost, bits[0], bits[1], cand[0]);
        for (uint32_t i = 2u; i < len; ++i) {
            qk_and_compute(kb, cost, cand[i - 2u], bits[i], cand[i - 1u]);
        }
        qk_cx(kb, cost, cand[len - 2u], dst);
        for (uint32_t i = len - 1u; i >= 2u; --i) {
            qk_and_uncompute(kb, cost, cand[i - 2u], bits[i], cand[i - 1u], (*sb)++);
        }
        qk_and_uncompute(kb, cost, bits[0], bits[1], cand[0], (*sb)++);
    }
    for (uint32_t i = 0u; i < len; ++i) qk_x(kb, cost, bits[i]);
}

static inline void qk_clz_norm(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *x_q,
    uint16_t iw,
    const uint32_t *mant_q,
    uint16_t mw,
    const uint32_t *clz_q,
    uint16_t clz_w,
    const uint32_t *anc, /* length >= iw + iw + clz_w + 4 */
    uint32_t *next_cbit
) {
    if (iw == 0u) return;
    /* anc layout: work[iw], cond_anc[iw], raw_clz[clz_w], is_zero_flag[1], tmp[1] */
    uint32_t *work = (uint32_t *)qcc_xmalloc((size_t)iw * sizeof(uint32_t));
    uint32_t *cand = (uint32_t *)qcc_xmalloc((size_t)iw * sizeof(uint32_t));
    uint32_t *rclz = (uint32_t *)qcc_xmalloc((size_t)(clz_w ? clz_w : 1u) * sizeof(uint32_t));
    for (uint16_t i = 0u; i < iw; ++i) {
        work[i] = qk__q(anc, i);
        cand[i] = qk__q(anc, iw + i);
        qk_r(kb, cost, work[i]);
        qk_cx(kb, cost, qk__q(x_q, i), work[i]);
    }
    for (uint16_t k = 0u; k < clz_w; ++k) {
        rclz[k] = qk__q(anc, 2u * iw + k);
        qk_r(kb, cost, rclz[k]);
    }
    uint32_t sb = next_cbit ? *next_cbit : 0u;
    uint32_t is_z = qk__q(anc, 2u * iw + clz_w);
    if (cost && (uint32_t)(2u * iw + clz_w + 1u) > cost->temp_anc_needed) {
        cost->temp_anc_needed = 2u * iw + clz_w + 1u;
    }

    /* Compute is_z = all_zero(work[0..iw-1]) */
    qk_r(kb, cost, is_z);
    qk__toggle_all_zero(kb, cost, work, iw, is_z, cand, &sb);

    /* Barrel left-shift from MSB stage s = 1 << (clz_w-1) down to 1 */
    for (int32_t k = (int32_t)clz_w - 1; k >= 0; --k) {
        uint32_t s = UINT32_C(1) << (uint32_t)k;
        if (s >= iw) continue;
        qk__toggle_all_zero(kb, cost, work + (iw - s), s, rclz[k], cand, &sb);

        /* Conditionally shift left by s when rclz[k] == 1 (when x==0 work is 0 so CSWAP on 0s does nothing) */
        for (int32_t i = (int32_t)iw - 1; i >= (int32_t)s; --i) {
            qk_cswap(kb, cost, rclz[k], work[i], work[i - (int32_t)s]);
        }
    }

    /* Copy mantissa (top mw bits of work) into clean mant_q */
    uint16_t drop = (iw > mw) ? (uint16_t)(iw - mw) : 0u;
    for (uint16_t b = 0u; b < mw; ++b) {
        qk_cx(kb, cost, work[drop + b], qk__q(mant_q, b));
    }
    /* Output clz = is_z ? iw : rclz */
    for (uint16_t k = 0u; k < clz_w; ++k) {
        uint32_t wbit = (iw >> k) & 1u;
        uint32_t qk = qk__q(clz_q, k);
        /* qk = (is_z & wbit) ^ ((1 ^ is_z) & rclz[k]) = rclz[k] ^ (is_z & (rclz[k] ^ wbit)) */
        qk_cx(kb, cost, rclz[k], qk);
        if (wbit) qk_x(kb, cost, rclz[k]);
        qk_ccx(kb, cost, is_z, rclz[k], qk);
        if (wbit) qk_x(kb, cost, rclz[k]);
    }

    /* Cleanly uncompute barrel shift, rclz, is_z, and work back to |0> */
    for (uint16_t k = 0u; k < clz_w; ++k) {
        uint32_t s = UINT32_C(1) << k;
        if (s >= iw) continue;
        for (uint32_t i = s; i < iw; ++i) {
            qk_cswap(kb, cost, rclz[k], work[i], work[i - s]);
        }
        qk__toggle_all_zero(kb, cost, work + (iw - s), s, rclz[k], cand, &sb);
    }
    qk__toggle_all_zero(kb, cost, work, iw, is_z, cand, &sb);
    for (uint16_t i = 0u; i < iw; ++i) {
        qk_cx(kb, cost, qk__q(x_q, i), work[i]);
    }
    if (next_cbit) *next_cbit = sb;
    free(work);
    free(cand);
    free(rclz);
}

/* ============================================================================
 * Single Unified Plan Emitter / Counter (Zero-Discrepancy Core)
 * ============================================================================ */

static inline bool qir_codegen_or_count(
    const qir_phys_plan_t *plan,
    qcc_kmb_circuit_t *kb_or_null,
    qcc_kmb_stats_t *out_stats
) {
    if (!plan || !plan->mod || !plan->val_allocs) return false;
    const qir_module_t *m = plan->mod;
    qir_kernel_cost_t kc;
    memset(&kc, 0, sizeof(kc));

    uint32_t max_q = 0u, max_b = 0u;
    for (uint32_t v = 0u; v < m->num_vals; ++v) {
        const qir_val_alloc_t *va = &plan->val_allocs[v];
        for (uint16_t i = 0u; i < va->width; ++i) {
            if (va->qubits && va->qubits[i] + 1u > max_q) max_q = va->qubits[i] + 1u;
        }
        for (uint16_t i = 0u; i < va->aux_width; ++i) {
            if (va->aux_qubits && va->aux_qubits[i] + 1u > max_q) max_q = va->aux_qubits[i] + 1u;
        }
        for (uint16_t i = 0u; i < va->cwidth; ++i) {
            if (va->cbits && va->cbits[i] + 1u > max_b) max_b = va->cbits[i] + 1u;
        }
        if (va->dirty_tail_len > 0u && va->dirty_tail_start + va->dirty_tail_len > max_q) {
            max_q = va->dirty_tail_start + va->dirty_tail_len;
        }
    }

    /* O(num_steps) precompute folded ops from R1/R3/R4/R9/R10b/R10c/R11 super-nodes */
    uint8_t *folded_op = m->num_ops ? (uint8_t *)calloc(m->num_ops, 1u) : NULL;
    for (uint32_t sj = 0u; sj < plan->num_steps; ++sj) {
        const qir_sched_step_t *s2 = &plan->steps[sj];
        if (s2->num_chain_ops == 0u || !s2->chain_ops) continue;
        if (s2->variant == QIR_VAR_CSA_TREE) {
            for (uint16_t ci = 0u; ci + 1u < s2->num_chain_ops; ++ci) {
                uint32_t co = s2->chain_ops[ci];
                if (co < m->num_ops && co != s2->op_id) folded_op[co] = 1u;
            }
        } else if (s2->variant == QIR_VAR_LINF2_GAUSS ||
                   s2->variant == QIR_VAR_PHASE_CSA_DIRECT ||
                   s2->variant == QIR_VAR_MUL_CONST_ADD ||
                   s2->variant == QIR_VAR_CTRL_ADD_CUCCARO ||
                   s2->variant == QIR_VAR_XOR_ACCUM ||
                   s2->variant == QIR_VAR_GF2M_CONST_LINF2) {
            for (uint16_t ci = 0u; ci < s2->num_chain_ops; ++ci) {
                uint32_t co = s2->chain_ops[ci];
                if (co < m->num_ops && co != s2->op_id) folded_op[co] = 1u;
            }
        }
    }

    for (uint32_t si = 0u; si < plan->num_steps; ++si) {
        const qir_sched_step_t *st = &plan->steps[si];
        const qir_op_t *op = &m->ops[st->op_id];
        uint32_t cbit = st->cbit_base;
        for (uint32_t a = 0u; a < st->num_temp_anc; ++a) {
            if (st->temp_anc && st->temp_anc[a] + 1u > max_q) max_q = st->temp_anc[a] + 1u;
        }

        if (st->variant == QIR_VAR_TRIVIAL_ELIDE) continue;
        if (folded_op && folded_op[st->op_id]) continue;

        bool has_pred = (st->pred_cbit != UINT32_MAX);
        if (has_pred) {
            if (st->pred_cbit + 1u > max_b) max_b = st->pred_cbit + 1u;
            qk_push_cond(kb_or_null, &kc, st->pred_cbit);
        }

        const uint32_t *out_q = (op->num_outs > 0u) ? plan->val_allocs[op->outs[0]].qubits : NULL;
        const uint32_t *in0_q = (op->num_ins > 0u)  ? plan->val_allocs[op->ins[0]].qubits : NULL;
        const uint32_t *in1_q = (op->num_ins > 1u)  ? plan->val_allocs[op->ins[1]].qubits : NULL;
        const uint32_t *in2_q = (op->num_ins > 2u)  ? plan->val_allocs[op->ins[2]].qubits : NULL;
        uint16_t out_w = (op->num_outs > 0u) ? m->vals[op->outs[0]].width : 0u;
        bool uncomp = (st->kind == QIR_STEP_UNCOMPUTE);

        switch (op->kind) {
            case QIR_OP_CONST: {
                if (op->num_outs > 0u && m->vals[op->outs[0]].rep == QIR_REP_CLASSICAL) {
                    const uint32_t *out_b = plan->val_allocs[op->outs[0]].cbits;
                    for (uint16_t b = 0u; b < out_w; ++b) {
                        uint32_t cb = qk__q(out_b, b);
                        if (qcc_bv_get(op->u.cst.val, b)) qk_bit_store1(kb_or_null, &kc, cb);
                        else                              qk_bit_store0(kb_or_null, &kc, cb);
                    }
                } else {
                    for (uint16_t b = 0u; b < out_w; ++b) {
                        if (qcc_bv_get(op->u.cst.val, b)) qk_x(kb_or_null, &kc, qk__q(out_q, b));
                    }
                }
                break;
            }
            case QIR_OP_SLICE:
            case QIR_OP_CONCAT: {
                /* Zero-cost wire aliasing when marked alias */
                if (op->num_outs > 0u && !plan->val_allocs[op->outs[0]].is_alias) {
                    if (m->vals[op->outs[0]].rep == QIR_REP_CLASSICAL) {
                        const uint32_t *out_b = plan->val_allocs[op->outs[0]].cbits;
                        if (op->kind == QIR_OP_SLICE) {
                            const uint32_t *in0_b = plan->val_allocs[op->ins[0]].cbits;
                            for (uint16_t b = 0u; b < out_w; ++b) {
                                uint32_t dst = qk__q(out_b, b);
                                qk_bit_store0(kb_or_null, &kc, dst);
                                qk_bit_invert_if(kb_or_null, &kc, dst, qk__q(in0_b, op->u.slice.lo + b));
                            }
                        } else {
                            uint16_t pos = 0u;
                            for (uint16_t i = 0u; i < op->num_ins; ++i) {
                                const uint32_t *ib = plan->val_allocs[op->ins[i]].cbits;
                                uint16_t iw = m->vals[op->ins[i]].width;
                                for (uint16_t b = 0u; b < iw && pos < out_w; ++b, ++pos) {
                                    uint32_t dst = qk__q(out_b, pos);
                                    qk_bit_store0(kb_or_null, &kc, dst);
                                    qk_bit_invert_if(kb_or_null, &kc, dst, qk__q(ib, b));
                                }
                            }
                        }
                    } else if (op->kind == QIR_OP_SLICE) {
                        for (uint16_t b = 0u; b < out_w; ++b) {
                            qk_cx(kb_or_null, &kc, qk__q(in0_q, op->u.slice.lo + b), qk__q(out_q, b));
                        }
                    } else {
                        uint16_t pos = 0u;
                        for (uint16_t i = 0u; i < op->num_ins; ++i) {
                            const uint32_t *iq = plan->val_allocs[op->ins[i]].qubits;
                            uint16_t iw = m->vals[op->ins[i]].width;
                            for (uint16_t b = 0u; b < iw && pos < out_w; ++b, ++pos) {
                                qk_cx(kb_or_null, &kc, qk__q(iq, b), qk__q(out_q, pos));
                            }
                        }
                    }
                }
                break;
            }
            case QIR_OP_LINF2: {
                const qir_linf2_payload_t *lf = st->fused_linf2 ? st->fused_linf2 : &op->u.linf2;
                const qir_op_t *src_op = (st->fused_linf2 && st->num_chain_ops > 0u)
                                             ? &m->ops[st->aux_param]
                                             : op;
                const uint32_t *src_in0_q = (src_op->num_ins > 0u) ? plan->val_allocs[src_op->ins[0]].qubits : in0_q;
                if (st->inplace_on_in0) {
                    if (!uncomp) {
                        qk_linf2_inplace(kb_or_null, &kc, src_in0_q, lf->out_bits, lf->rows, &lf->const_vec);
                    } else {
                        /* Invert LINF2 matrix over F_2 to restore in-place if needed */
                        uint32_t w = lf->out_bits;
                        qcc_bv256_t *M = (qcc_bv256_t *)qcc_xmalloc((size_t)(w ? w : 1u) * sizeof(qcc_bv256_t));
                        qcc_bv256_t *I = (qcc_bv256_t *)qcc_xcalloc(w ? w : 1u, sizeof(qcc_bv256_t));
                        memcpy(M, lf->rows, (size_t)w * sizeof(qcc_bv256_t));
                        for (uint32_t r = 0u; r < w; ++r) qcc_bv_set(I[r].w, r);
                        for (uint32_t col = 0u; col < w; ++col) {
                            uint32_t piv = col;
                            while (piv < w && !qcc_bv_get(M[piv].w, col)) ++piv;
                            if (piv == w) continue;
                            if (piv != col) {
                                qcc_bv256_t tm = M[col]; M[col] = M[piv]; M[piv] = tm;
                                qcc_bv256_t ti = I[col]; I[col] = I[piv]; I[piv] = ti;
                            }
                            for (uint32_t r = 0u; r < w; ++r) {
                                if (r != col && qcc_bv_get(M[r].w, col)) {
                                    for (uint32_t ww = 0u; ww < QIR_WORDS; ++ww) {
                                        M[r].w[ww] ^= M[col].w[ww];
                                        I[r].w[ww] ^= I[col].w[ww];
                                    }
                                }
                            }
                        }
                        for (uint32_t r = 0u; r < w; ++r) {
                            if (qcc_bv_get(lf->const_vec.w, r)) qk_x(kb_or_null, &kc, qk__q(src_in0_q, r));
                        }
                        qk_linf2_inplace(kb_or_null, &kc, src_in0_q, w, I, NULL);
                        free(M);
                        free(I);
                    }
                } else {
                    uint32_t *flat_in = (uint32_t *)qcc_xmalloc((size_t)(lf->in_bits ? lf->in_bits : 1u) * sizeof(uint32_t));
                    uint32_t p = 0u;
                    for (uint16_t i = 0u; i < src_op->num_ins; ++i) {
                        uint32_t vid = src_op->ins[i];
                        const uint32_t *vq = plan->val_allocs[vid].qubits;
                        for (uint16_t b = 0u; b < m->vals[vid].width; ++b) {
                            flat_in[p++] = qk__q(vq, b);
                        }
                    }
                    qk_linf2_out_of_place(kb_or_null, &kc, flat_in, lf->in_bits, out_q, lf->out_bits, lf->rows, &lf->const_vec);
                    free(flat_in);
                }
                break;
            }
            case QIR_OP_AND_VEC: {
                const uint64_t *zm = m->vals[op->outs[0]].zero_mask;
                for (uint16_t b = 0u; b < out_w; ++b) {
                    if (qcc_bv_get(zm, b)) continue;
                    if (!uncomp) qk_and_compute(kb_or_null, &kc, in0_q[b], in1_q[b], out_q[b]);
                    else         qk_and_uncompute(kb_or_null, &kc, in0_q[b], in1_q[b], out_q[b], cbit++);
                }
                break;
            }
            case QIR_OP_ADD: {
                if (st->variant == QIR_VAR_CSA_TREE && st->num_chain_ops > 1u) {
                    /* Collect summand registers from chained binary ADDs */
                    uint32_t sum_vals[16];
                    uint32_t n_sum = 0u;
                    for (uint16_t ci = 0u; ci < st->num_chain_ops; ++ci) {
                        const qir_op_t *cop = &m->ops[st->chain_ops[ci]];
                        for (uint16_t ii = 0u; ii < 2u; ++ii) {
                            uint32_t iv = cop->ins[ii];
                            bool produced_in_chain = false;
                            if (m->vals[iv].def_op != UINT32_MAX) {
                                for (uint16_t cj = 0u; cj < st->num_chain_ops; ++cj) {
                                    if (st->chain_ops[cj] == m->vals[iv].def_op) { produced_in_chain = true; break; }
                                }
                            }
                            if (!produced_in_chain && n_sum < 16u) sum_vals[n_sum++] = iv;
                        }
                    }
                    if (n_sum >= 3u) {
                        /* Stage reductions into temp_anc pairs `(S_k, C_k)` of width out_w */
                        const uint32_t *cur_s = plan->val_allocs[sum_vals[0]].qubits;
                        const uint32_t *cur_c = plan->val_allocs[sum_vals[1]].qubits;
                        uint32_t stages = n_sum - 2u;
                        for (uint32_t s = 0u; s < stages; ++s) {
                            const uint32_t *nxt_s = st->temp_anc ? (st->temp_anc + (size_t)(2u * s) * out_w) : NULL;
                            const uint32_t *nxt_c = nxt_s ? (nxt_s + out_w) : NULL;
                            qk_cs_add3(kb_or_null, &kc, cur_s, cur_c, plan->val_allocs[sum_vals[2u + s]].qubits, nxt_s, nxt_c, out_w);
                            cur_s = nxt_s;
                            cur_c = nxt_c;
                        }
                        const uint32_t *add_anc = st->temp_anc ? (st->temp_anc + (size_t)(2u * stages) * out_w) : NULL;
                        if (!uncomp) {
                            for (uint16_t b = 0u; b < out_w; ++b) qk_cx(kb_or_null, &kc, qk__q(cur_c, b), qk__q(out_q, b));
                            qk_gidney_iadd(kb_or_null, &kc, cur_s, out_q, add_anc, out_w, &cbit);
                        } else {
                            qk_gidney_isub_asym(kb_or_null, &kc, cur_s, out_w, out_q, out_w, add_anc, &cbit);
                            for (uint16_t b = 0u; b < out_w; ++b) qk_cx(kb_or_null, &kc, qk__q(cur_c, b), qk__q(out_q, b));
                        }
                        /* Uncompute intermediate CSA pairs cleanly in reverse using 0-CCX HMR dual! */
                        for (int32_t s = (int32_t)stages - 1; s >= 0; --s) {
                            const uint32_t *ps = (s == 0) ? plan->val_allocs[sum_vals[0]].qubits
                                                          : (st->temp_anc ? (st->temp_anc + (size_t)(2u * ((uint32_t)s - 1u)) * out_w) : NULL);
                            const uint32_t *pc = (s == 0) ? plan->val_allocs[sum_vals[1]].qubits
                                                          : (ps ? (ps + out_w) : NULL);
                            const uint32_t *ns = st->temp_anc ? (st->temp_anc + (size_t)(2u * (uint32_t)s) * out_w) : NULL;
                            const uint32_t *nc = ns ? (ns + out_w) : NULL;
                            qk_cs_add3_uncompute(kb_or_null, &kc, ps, pc, plan->val_allocs[sum_vals[2u + (uint32_t)s]].qubits, ns, nc, out_w, &cbit);
                        }
                        break;
                    }
                }
                if (st->variant == QIR_VAR_INPLACE_ARITH && st->inplace_on_in0 && st->inplace_in_idx < 2u) {
                    uint8_t tidx = st->inplace_in_idx;
                    uint32_t lhs_vid = op->ins[tidx ^ 1u];
                    const uint32_t *lhs_q = plan->val_allocs[lhs_vid].qubits;
                    uint16_t wl = m->vals[lhs_vid].width;
                    if (!uncomp) qk_gidney_iadd_asym(kb_or_null, &kc, lhs_q, wl, out_q, out_w, st->temp_anc, &cbit);
                    else         qk_gidney_isub_asym(kb_or_null, &kc, lhs_q, wl, out_q, out_w, st->temp_anc, &cbit);
                    break;
                }
                uint16_t w0 = m->vals[op->ins[0]].width;
                uint16_t w1 = m->vals[op->ins[1]].width;
                uint16_t c0 = (w0 < out_w) ? w0 : out_w;
                uint16_t c1 = (w1 < out_w) ? w1 : out_w;
                if (!uncomp) {
                    for (uint16_t b = 0u; b < c1; ++b) qk_cx(kb_or_null, &kc, qk__q(in1_q, b), qk__q(out_q, b));
                    qk_gidney_iadd_asym(kb_or_null, &kc, in0_q, c0, out_q, out_w, st->temp_anc, &cbit);
                } else {
                    qk_gidney_isub_asym(kb_or_null, &kc, in0_q, c0, out_q, out_w, st->temp_anc, &cbit);
                    for (uint16_t b = 0u; b < c1; ++b) qk_cx(kb_or_null, &kc, qk__q(in1_q, b), qk__q(out_q, b));
                }
                break;
            }
            case QIR_OP_SUB: {
                if (st->variant == QIR_VAR_INPLACE_ARITH && st->inplace_on_in0 && st->inplace_in_idx < 2u) {
                    uint8_t tidx = st->inplace_in_idx;
                    uint32_t lhs_vid = op->ins[tidx ^ 1u];
                    const uint32_t *lhs_q = plan->val_allocs[lhs_vid].qubits;
                    uint16_t wl = m->vals[lhs_vid].width;
                    if (tidx == 0u) {
                        if (!uncomp) qk_gidney_isub_asym(kb_or_null, &kc, lhs_q, wl, out_q, out_w, st->temp_anc, &cbit);
                        else         qk_gidney_iadd_asym(kb_or_null, &kc, lhs_q, wl, out_q, out_w, st->temp_anc, &cbit);
                    } else {
                        for (uint16_t b = 0u; b < out_w; ++b) qk_x(kb_or_null, &kc, qk__q(out_q, b));
                        qk_gidney_iadd_p1(kb_or_null, &kc, lhs_q, out_q, st->temp_anc, out_w, &cbit);
                    }
                    break;
                }
                uint16_t w0 = m->vals[op->ins[0]].width;
                uint16_t w1 = m->vals[op->ins[1]].width;
                uint16_t c0 = (w0 < out_w) ? w0 : out_w;
                uint16_t c1 = (w1 < out_w) ? w1 : out_w;
                if (!uncomp) {
                    for (uint16_t b = 0u; b < c0; ++b) qk_cx(kb_or_null, &kc, qk__q(in0_q, b), qk__q(out_q, b));
                    qk_gidney_isub_asym(kb_or_null, &kc, in1_q, c1, out_q, out_w, st->temp_anc, &cbit);
                } else {
                    qk_gidney_iadd_asym(kb_or_null, &kc, in1_q, c1, out_q, out_w, st->temp_anc, &cbit);
                    for (uint16_t b = 0u; b < c0; ++b) qk_cx(kb_or_null, &kc, qk__q(in0_q, b), qk__q(out_q, b));
                }
                break;
            }
            case QIR_OP_CS_ADD3: {
                const uint32_t *s_q = plan->val_allocs[op->outs[0]].qubits;
                const uint32_t *c_q = plan->val_allocs[op->outs[1]].qubits;
                if (!uncomp) qk_cs_add3(kb_or_null, &kc, in0_q, in1_q, in2_q, s_q, c_q, out_w);
                else         qk_cs_add3_uncompute(kb_or_null, &kc, in0_q, in1_q, in2_q, s_q, c_q, out_w, &cbit);
                break;
            }
            case QIR_OP_MUL_FXP: {
                int32_t sh_s = (int32_t)op->u.mul.frac_a + (int32_t)op->u.mul.frac_b - (int32_t)op->u.mul.frac_out;
                uint32_t shift = sh_s > 0 ? (uint32_t)sh_s : 0u;
                uint32_t keep = (op->flags & QIR_FLAG_APPROX) ? st->approx_keep_bits : 0u;
                const uint32_t *low_q = plan->val_allocs[op->outs[0]].aux_qubits;
                if (st->variant == QIR_VAR_MUL_CONST_ADD && op->num_ins >= 2u) {
                    uint32_t c_idx = (st->aux_param & 1u);
                    uint32_t cvid  = op->ins[c_idx];
                    uint32_t vvid  = op->ins[c_idx ^ 1u];
                    const uint32_t *vq = plan->val_allocs[vvid].qubits;
                    qk_mul_const_fxp(kb_or_null, &kc,
                                     vq, m->vals[vvid].width, m->vals[vvid].zero_mask,
                                     m->vals[cvid].const_val, m->vals[cvid].width,
                                     low_q, out_q, out_w, shift, keep,
                                     st->temp_anc, &cbit, uncomp);
                } else {
                    qk_mul_fxp_ex(kb_or_null, &kc,
                                  in0_q, m->vals[op->ins[0]].width, m->vals[op->ins[0]].zero_mask,
                                  in1_q, m->vals[op->ins[1]].width, m->vals[op->ins[1]].zero_mask,
                                  low_q, out_q, out_w, shift, keep, st->temp_anc, &cbit, uncomp);
                }
                break;
            }
            case QIR_OP_SQR_FXP: {
                uint16_t wa = m->vals[op->ins[0]].width;
                int32_t sh_s = 2 * (int32_t)op->u.mul.frac_a - (int32_t)op->u.mul.frac_out;
                uint32_t shift = sh_s > 0 ? (uint32_t)sh_s : 0u;
                uint32_t keep = (op->flags & QIR_FLAG_APPROX) ? st->approx_keep_bits : 0u;
                const uint32_t *low_q = plan->val_allocs[op->outs[0]].aux_qubits;
                if (st->variant == QIR_VAR_SQR_TRIANGLE) {
                    qk_sqr_fxp_tri(kb_or_null, &kc,
                                   in0_q, wa, m->vals[op->ins[0]].zero_mask,
                                   low_q, out_q, out_w, shift, keep, st->temp_anc, &cbit, uncomp);
                } else {
                    uint32_t full = shift + out_w;
                    uint32_t minp = (keep > 0u && keep < full) ? (full - keep) : 0u;
                    if (minp > shift) minp = shift;
                    uint32_t acc_bits = full - minp;
                    const uint32_t *mul_anc = st->temp_anc;
                    const uint32_t *a_copy  = st->temp_anc ? (st->temp_anc + acc_bits) : NULL;
                    for (uint16_t b = 0u; b < wa; ++b) {
                        uint32_t acq = qk__q(a_copy, b);
                        qk_r(kb_or_null, &kc, acq);
                        qk_cx(kb_or_null, &kc, qk__q(in0_q, b), acq);
                    }
                    qk_mul_fxp_ex(kb_or_null, &kc,
                                  a_copy, wa, m->vals[op->ins[0]].zero_mask,
                                  in0_q, wa, m->vals[op->ins[0]].zero_mask,
                                  low_q, out_q, out_w, shift, keep, mul_anc, &cbit, uncomp);
                    for (uint16_t b = 0u; b < wa; ++b) {
                        qk_cx(kb_or_null, &kc, qk__q(in0_q, b), qk__q(a_copy, b));
                    }
                }
                break;
            }
            case QIR_OP_GF2M_ADD: {
                if (st->variant == QIR_VAR_XOR_ACCUM && st->num_chain_ops > 0u) {
                    const qir_op_t *dop = &m->ops[st->chain_ops[0]];
                    uint32_t acc_idx = st->aux_param & 3u;
                    qir_variant_kind_t sub_var = (qir_variant_kind_t)((st->aux_param >> 8u) & 0xFFu);
                    uint32_t sub_aux = st->aux_param >> 16u;
                    const uint32_t *acc_q = plan->val_allocs[op->ins[acc_idx]].qubits;
                    bool ip = st->inplace_on_in0 && (st->inplace_in_idx == acc_idx);
                    if (!uncomp && !ip) {
                        for (uint16_t b = 0u; b < out_w; ++b) qk_cx(kb_or_null, &kc, qk__q(acc_q, b), qk__q(out_q, b));
                    }
                    if (dop->kind == QIR_OP_LUT_ROM) {
                        const qir_lut_payload_t *lut = &dop->u.lut;
                        const uint32_t *k_q = plan->val_allocs[dop->ins[0]].qubits;
                        uint16_t vb = lut->val_bits;
                        if (sub_var == QIR_VAR_QROAM_VSLOT0) {
                            qk_qroam_vslot0(kb_or_null, &kc, lut->table, lut->num_entries, lut->val_words,
                                            0u, vb, k_q, lut->key_bits, out_q, sub_aux ? sub_aux : 1u,
                                            st->temp_anc, st->num_temp_anc, &cbit);
                        } else if (sub_var == QIR_VAR_QROAM_SPLIT) {
                            uint32_t K = sub_aux & 0xFFu;
                            uint32_t b_sh = (sub_aux >> 8u) & 0xFFu;
                            if (K < 2u) K = 2u;
                            if (b_sh == 0u) b_sh = 1u;
                            uint32_t chunk = (vb + K - 1u) / K;
                            if (chunk == 0u) chunk = 1u;
                            for (uint32_t lo = 0u; lo < vb; lo += chunk) {
                                uint32_t hi = (lo + chunk < vb) ? (lo + chunk) : vb;
                                qk_qroam_vslot0(kb_or_null, &kc, lut->table, lut->num_entries, lut->val_words,
                                                lo, hi, k_q, lut->key_bits, out_q ? (out_q + lo) : NULL,
                                                b_sh, st->temp_anc, st->num_temp_anc, &cbit);
                            }
                        } else if (sub_var == QIR_VAR_LIN_GAUGE_TSP && st->gauge_lut) {
                            const qir_gauge_lut_t *gl = st->gauge_lut;
                            uint32_t rb = gl->rank_bits;
                            const uint32_t *comp_q   = st->temp_anc;
                            const uint32_t *qrom_anc = st->temp_anc ? (st->temp_anc + rb) : NULL;
                            for (uint32_t r = 0u; r < rb; ++r) qk_r(kb_or_null, &kc, qk__q(comp_q, r));
                            qk_lut_rom_qrom(kb_or_null, &kc, gl->comp_table, gl->num_entries, gl->comp_words,
                                            0u, rb, k_q, lut->key_bits, comp_q, qrom_anc, &cbit);
                            qk_linf2_out_of_place(kb_or_null, &kc, comp_q, rb, out_q, gl->val_bits, gl->recon_rows, &gl->const_shift);
                            if (!qk_hmr_parity_lut_uncompute(kb_or_null, &kc, gl->comp_table, gl->num_entries, gl->comp_words,
                                                             0u, rb, k_q, lut->key_bits, comp_q, (uint32_t *)qrom_anc, &cbit)) {
                                qk_lut_rom_qrom(kb_or_null, &kc, gl->comp_table, gl->num_entries, gl->comp_words,
                                                0u, rb, k_q, lut->key_bits, comp_q, qrom_anc, &cbit);
                            }
                        } else {
                            qk_lut_rom_qrom(kb_or_null, &kc, lut->table, lut->num_entries, lut->val_words,
                                            0u, vb, k_q, lut->key_bits, out_q, st->temp_anc, &cbit);
                        }
                    } else if (dop->kind == QIR_OP_GF2M_MUL) {
                        const uint32_t *da_q = plan->val_allocs[dop->ins[0]].qubits;
                        const uint32_t *db_q = plan->val_allocs[dop->ins[1]].qubits;
                        if (dop->u.gf2m.m == 12u && plan->has_gf12_r42 && plan->gf12_r42_poly == dop->u.gf2m.irred_poly) {
                            qk_gf12_mul_r42_compute(kb_or_null, &kc, da_q, db_q, out_q, st->temp_anc[0], &cbit, &plan->gf12_r42);
                        } else if (dop->u.gf2m.m == 8u && plan->has_gf256_r27 && plan->gf256_r27_poly == dop->u.gf2m.irred_poly) {
                            qk_gf256_mul_r27_compute(kb_or_null, &kc, da_q, db_q, out_q, st->temp_anc[0], &cbit, &plan->gf256_r27);
                        } else if (dop->u.gf2m.m == 6u && plan->has_gf64_r15 && plan->gf64_r15_poly == dop->u.gf2m.irred_poly) {
                            qk_gf64_mul_r15_compute(kb_or_null, &kc, da_q, db_q, out_q, st->temp_anc[0], &cbit, &plan->gf64_r15);
                        } else {
                            uint16_t mdeg = dop->u.gf2m.m;
                            uint32_t poly = dop->u.gf2m.irred_poly;
                            for (uint16_t j = 0u; j < mdeg; ++j) {
                                for (uint16_t i = 0u; i < mdeg; ++i) {
                                    uint32_t shift_elem = UINT32_C(1) << i;
                                    for (uint16_t s = 0u; s < j; ++s) {
                                        shift_elem <<= 1u;
                                        if (shift_elem & (UINT32_C(1) << mdeg)) shift_elem ^= poly;
                                    }
                                    uint32_t p = st->temp_anc[0];
                                    qk_and_compute(kb_or_null, &kc, da_q[i], db_q[j], p);
                                    for (uint16_t k = 0u; k < mdeg; ++k) {
                                        if ((shift_elem >> k) & 1u) qk_cx(kb_or_null, &kc, p, out_q[k]);
                                    }
                                    qk_and_uncompute(kb_or_null, &kc, da_q[i], db_q[j], p, cbit++);
                                }
                            }
                        }
                    } else if (dop->kind == QIR_OP_LINF2) {
                        const qir_linf2_payload_t *lf = &dop->u.linf2;
                        const uint32_t *dq = plan->val_allocs[dop->ins[0]].qubits;
                        qk_linf2_out_of_place(kb_or_null, &kc, dq, lf->in_bits, out_q, lf->out_bits, lf->rows, &lf->const_vec);
                    }
                    if (uncomp && !ip) {
                        for (uint16_t b = 0u; b < out_w; ++b) qk_cx(kb_or_null, &kc, qk__q(acc_q, b), qk__q(out_q, b));
                    }
                    break;
                }
                if (st->variant == QIR_VAR_INPLACE_ARITH && st->inplace_on_in0 && st->inplace_in_idx < 2u) {
                    const uint32_t *other_q = plan->val_allocs[op->ins[st->inplace_in_idx ^ 1u]].qubits;
                    uint16_t deg = op->u.gf2m.m;
                    for (uint16_t b = 0u; b < deg; ++b) qk_cx(kb_or_null, &kc, qk__q(other_q, b), qk__q(out_q, b));
                    break;
                }
                uint16_t deg = op->u.gf2m.m;
                for (uint16_t b = 0u; b < deg; ++b) {
                    qk_cx(kb_or_null, &kc, in0_q[b], out_q[b]);
                    qk_cx(kb_or_null, &kc, in1_q[b], out_q[b]);
                }
                break;
            }
            case QIR_OP_GF2M_MUL: {
                if (st->variant == QIR_VAR_GF2M_CONST_LINF2 && st->fused_linf2) {
                    const qir_linf2_payload_t *lf = st->fused_linf2;
                    uint32_t var_idx = st->aux_param & 1u;
                    const uint32_t *vq = plan->val_allocs[op->ins[var_idx]].qubits;
                    uint32_t dm = lf->out_bits;
                    uint32_t c_idx = var_idx ^ 1u;
                    uint32_t dmask = (dm == 32u) ? UINT32_MAX : ((UINT32_C(1) << dm) - 1u);
                    uint32_t celem = (uint32_t)(m->vals[op->ins[c_idx]].const_val[0] & dmask);
                    if (celem == 0u) break;
                    if (st->inplace_on_in0 && st->inplace_in_idx == var_idx) {
                        if (!uncomp) {
                            qk_linf2_inplace(kb_or_null, &kc, vq, dm, lf->rows, NULL);
                        } else {
                            uint32_t cinv = qir__gf2m_inv_u32(celem, (uint16_t)dm, op->u.gf2m.irred_poly);
                            qcc_bv256_t *inv_rows = (qcc_bv256_t *)qcc_xcalloc(dm ? dm : 1u, sizeof(qcc_bv256_t));
                            for (uint32_t c = 0u; c < dm; ++c) {
                                uint32_t col = qir__gf2m_mul_u32(cinv, UINT32_C(1) << c, (uint16_t)dm, op->u.gf2m.irred_poly);
                                for (uint32_t r = 0u; r < dm; ++r) {
                                    if ((col >> r) & 1u) qcc_bv_set(inv_rows[r].w, c);
                                }
                            }
                            qk_linf2_inplace(kb_or_null, &kc, vq, dm, inv_rows, NULL);
                            free(inv_rows);
                        }
                    } else {
                        qk_linf2_out_of_place(kb_or_null, &kc, vq, lf->in_bits, out_q, dm, lf->rows, NULL);
                    }
                    break;
                }
                if (op->u.gf2m.m == 12u && plan->has_gf12_r42 && plan->gf12_r42_poly == op->u.gf2m.irred_poly) {
                    if (!uncomp) {
                        qk_gf12_mul_r42_compute(kb_or_null, &kc, in0_q, in1_q, out_q, st->temp_anc[0], &cbit, &plan->gf12_r42);
                    } else if (st->variant == QIR_VAR_HMR_PARITY_LUT) {
                        qk_gf12_mul_r42_uncompute_hmr(kb_or_null, &kc, in0_q, in1_q, out_q, &cbit, &plan->gf12_r42);
                    } else {
                        qk_gf12_mul_r42_compute(kb_or_null, &kc, in0_q, in1_q, out_q, st->temp_anc[0], &cbit, &plan->gf12_r42);
                    }
                } else if (op->u.gf2m.m == 8u && plan->has_gf256_r27 && plan->gf256_r27_poly == op->u.gf2m.irred_poly) {
                    if (!uncomp) {
                        qk_gf256_mul_r27_compute(kb_or_null, &kc, in0_q, in1_q, out_q, st->temp_anc[0], &cbit, &plan->gf256_r27);
                    } else if (st->variant == QIR_VAR_HMR_PARITY_LUT) {
                        qk_gf256_mul_r27_uncompute_hmr(kb_or_null, &kc, in0_q, in1_q, out_q, &cbit, &plan->gf256_r27);
                    } else {
                        qk_gf256_mul_r27_compute(kb_or_null, &kc, in0_q, in1_q, out_q, st->temp_anc[0], &cbit, &plan->gf256_r27);
                    }
                } else if (op->u.gf2m.m == 6u && plan->has_gf64_r15 && plan->gf64_r15_poly == op->u.gf2m.irred_poly) {
                    if (!uncomp) qk_gf64_mul_r15_compute(kb_or_null, &kc, in0_q, in1_q, out_q, st->temp_anc[0], &cbit, &plan->gf64_r15);
                    else         qk_gf64_mul_r15_uncompute_hmr(kb_or_null, &kc, in0_q, in1_q, out_q, &cbit, &plan->gf64_r15);
                } else {
                    /* Generic GF(2^m) schoolbook shift-and-reduce */
                    uint16_t mdeg = op->u.gf2m.m;
                    uint32_t poly = op->u.gf2m.irred_poly;
                    for (uint16_t j = 0u; j < mdeg; ++j) {
                        for (uint16_t i = 0u; i < mdeg; ++i) {
                            uint32_t shift_elem = UINT32_C(1) << i;
                            for (uint16_t s = 0u; s < j; ++s) {
                                shift_elem <<= 1u;
                                if (shift_elem & (UINT32_C(1) << mdeg)) shift_elem ^= poly;
                            }
                            uint32_t p = st->temp_anc[0];
                            qk_and_compute(kb_or_null, &kc, in0_q[i], in1_q[j], p);
                            for (uint16_t k = 0u; k < mdeg; ++k) {
                                if ((shift_elem >> k) & 1u) qk_cx(kb_or_null, &kc, p, out_q[k]);
                            }
                            qk_and_uncompute(kb_or_null, &kc, in0_q[i], in1_q[j], p, cbit++);
                        }
                    }
                }
                break;
            }
            case QIR_OP_GF2M_INV: {
                uint32_t ctrl_q = (op->num_ins >= 2u && in1_q) ? qk__q(in1_q, 0u) : UINT32_MAX;
                if (op->u.gf2m.m == 12u && plan->has_gf12_tower && plan->gf12_tower_poly == op->u.gf2m.irred_poly) {
                    if (!uncomp) {
                        qk_gf12_inv_tower(kb_or_null, &kc, in0_q, out_q, ctrl_q, st->temp_anc, &cbit, &plan->gf12_tower);
                    } else {
                        qk_gf12_inv_tower_uncompute(kb_or_null, &kc, in0_q, out_q, ctrl_q, st->temp_anc, &cbit, &plan->gf12_tower);
                    }
                } else if (op->u.gf2m.m == 8u && st->variant == QIR_VAR_GF2M_TOWER_INV &&
                           plan->has_gf256_tower && plan->gf256_tower_poly == op->u.gf2m.irred_poly) {
                    if (!uncomp) {
                        qk_gf256_inv_tower(kb_or_null, &kc, in0_q, out_q, ctrl_q, st->temp_anc, &cbit, &plan->gf256_tower);
                    } else {
                        qk_gf256_inv_tower_uncompute(kb_or_null, &kc, in0_q, out_q, ctrl_q, st->temp_anc, &cbit, &plan->gf256_tower);
                    }
                } else {
                    /* Direct LUT for small field (e.g. deg <= 8) */
                    uint16_t deg = op->u.gf2m.m;
                    uint32_t nent = UINT32_C(1) << deg;
                    uint64_t *tbl = (uint64_t *)qcc_xcalloc(nent, sizeof(uint64_t));
                    for (uint32_t a = 0u; a < nent; ++a) tbl[a] = qir__gf2m_inv_u32(a, deg, op->u.gf2m.irred_poly);
                    if (ctrl_q == UINT32_MAX) {
                        qk_lut_rom_qrom(kb_or_null, &kc, tbl, nent, 1u, 0u, deg, in0_q, deg, out_q, st->temp_anc, &cbit);
                    } else {
                        const uint32_t *inv_tmp  = st->temp_anc;
                        const uint32_t *qrom_anc = st->temp_anc ? (st->temp_anc + deg) : NULL;
                        for (uint16_t b = 0u; b < deg; ++b) qk_r(kb_or_null, &kc, qk__q(inv_tmp, b));
                        qk_lut_rom_qrom(kb_or_null, &kc, tbl, nent, 1u, 0u, deg, in0_q, deg, inv_tmp, qrom_anc, &cbit);
                        for (uint16_t b = 0u; b < deg; ++b) qk_ccx(kb_or_null, &kc, ctrl_q, qk__q(inv_tmp, b), qk__q(out_q, b));
                        qk_lut_rom_qrom(kb_or_null, &kc, tbl, nent, 1u, 0u, deg, in0_q, deg, inv_tmp, qrom_anc, &cbit);
                    }
                    free(tbl);
                }
                break;
            }
            case QIR_OP_LUT_ROM: {
                const qir_lut_payload_t *lut = &op->u.lut;
                uint16_t vb = lut->val_bits;
                if ((op->flags & QIR_FLAG_APPROX) && st->approx_keep_bits > 0u && st->approx_keep_bits < vb) {
                    vb = st->approx_keep_bits;
                }
                if (uncomp && st->variant == QIR_VAR_HMR_PARITY_LUT &&
                    qk_hmr_parity_lut_uncompute(kb_or_null, &kc, lut->table, lut->num_entries, lut->val_words,
                                                0u, vb, in0_q, lut->key_bits, out_q, st->temp_anc, &cbit)) {
                    break;
                }
                if (st->variant == QIR_VAR_QROAM_VSLOT0) {
                    qk_qroam_vslot0(kb_or_null, &kc, lut->table, lut->num_entries, lut->val_words,
                                    0u, vb, in0_q, lut->key_bits, out_q,
                                    st->aux_param ? st->aux_param : 1u,
                                    st->temp_anc, st->num_temp_anc, &cbit);
                    break;
                }
                if (st->variant == QIR_VAR_QROAM_SPLIT) {
                    uint32_t K = st->aux_param & 0xFFu;
                    uint32_t b_sh = (st->aux_param >> 8u) & 0xFFu;
                    if (K < 2u) K = 2u;
                    if (b_sh == 0u) b_sh = 1u;
                    uint32_t chunk = (vb + K - 1u) / K;
                    if (chunk == 0u) chunk = 1u;
                    for (uint32_t lo = 0u; lo < vb; lo += chunk) {
                        uint32_t hi = (lo + chunk < vb) ? (lo + chunk) : vb;
                        qk_qroam_vslot0(kb_or_null, &kc, lut->table, lut->num_entries, lut->val_words,
                                        lo, hi, in0_q, lut->key_bits, out_q ? (out_q + lo) : NULL,
                                        b_sh, st->temp_anc, st->num_temp_anc, &cbit);
                    }
                    break;
                }
                if (st->variant == QIR_VAR_LIN_GAUGE_TSP && st->gauge_lut) {
                    const qir_gauge_lut_t *gl = st->gauge_lut;
                    uint32_t rb = gl->rank_bits;
                    const uint32_t *comp_q   = st->temp_anc;
                    const uint32_t *qrom_anc = st->temp_anc ? (st->temp_anc + rb) : NULL;
                    for (uint32_t r = 0u; r < rb; ++r) qk_r(kb_or_null, &kc, qk__q(comp_q, r));
                    qk_lut_rom_qrom(kb_or_null, &kc, gl->comp_table, gl->num_entries, gl->comp_words,
                                    0u, rb, in0_q, lut->key_bits, comp_q, qrom_anc, &cbit);
                    qk_linf2_out_of_place(kb_or_null, &kc, comp_q, rb, out_q, gl->val_bits, gl->recon_rows, &gl->const_shift);
                    if (!qk_hmr_parity_lut_uncompute(kb_or_null, &kc, gl->comp_table, gl->num_entries, gl->comp_words,
                                                     0u, rb, in0_q, lut->key_bits, comp_q, (uint32_t *)qrom_anc, &cbit)) {
                        qk_lut_rom_qrom(kb_or_null, &kc, gl->comp_table, gl->num_entries, gl->comp_words,
                                        0u, rb, in0_q, lut->key_bits, comp_q, qrom_anc, &cbit);
                    }
                    break;
                }
                if (st->variant == QIR_VAR_LUT_COL_SPLIT && st->aux_param >= 2u) {
                    uint32_t K = st->aux_param;
                    uint32_t chunk = (vb + K - 1u) / K;
                    if (chunk == 0u) chunk = 1u;
                    for (uint32_t lo = 0u; lo < vb; lo += chunk) {
                        uint32_t hi = (lo + chunk < vb) ? (lo + chunk) : vb;
                        qk_lut_rom_qrom(kb_or_null, &kc, lut->table, lut->num_entries, lut->val_words,
                                        lo, hi, in0_q, lut->key_bits, out_q ? (out_q + lo) : NULL, st->temp_anc, &cbit);
                    }
                    break;
                }
                qk_lut_rom_qrom(kb_or_null, &kc, lut->table, lut->num_entries, lut->val_words,
                                0u, vb, in0_q, lut->key_bits, out_q, st->temp_anc, &cbit);
                break;
            }
            case QIR_OP_CLZ_NORM: {
                const uint32_t *mant_q = plan->val_allocs[op->outs[0]].qubits;
                const uint32_t *clz_q  = plan->val_allocs[op->outs[1]].qubits;
                qk_clz_norm(kb_or_null, &kc, in0_q, m->vals[op->ins[0]].width,
                            mant_q, m->vals[op->outs[0]].width,
                            clz_q, m->vals[op->outs[1]].width,
                            st->temp_anc, &cbit);
                break;
            }
            case QIR_OP_CTRL_SEL: {
                uint32_t sel = qk__q(in0_q, 0u);
                if (st->variant == QIR_VAR_CTRL_ADD_CUCCARO && st->num_chain_ops > 0u) {
                    uint32_t lhs_vid = st->aux_param & UINT32_C(0x7FFFFFFF);
                    bool is_sub = (st->aux_param >> 31u) != 0u;
                    const uint32_t *lhs_q = plan->val_allocs[lhs_vid].qubits;
                    uint16_t k_in = m->vals[lhs_vid].width;
                    bool ip = st->inplace_on_in0 && (st->inplace_in_idx == 1u);
                    if (!uncomp) {
                        if (!ip) {
                            for (uint16_t b = 0u; b < out_w; ++b) qk_cx(kb_or_null, &kc, qk__q(in1_q, b), qk__q(out_q, b));
                        }
                        if (!is_sub) qk_ctrl_iadd_cuccaro(kb_or_null, &kc, sel, lhs_q, k_in, out_q, out_w, st->temp_anc, &cbit);
                        else         qk_ctrl_isub_cuccaro(kb_or_null, &kc, sel, lhs_q, k_in, out_q, out_w, st->temp_anc, &cbit);
                    } else {
                        if (!is_sub) qk_ctrl_isub_cuccaro(kb_or_null, &kc, sel, lhs_q, k_in, out_q, out_w, st->temp_anc, &cbit);
                        else         qk_ctrl_iadd_cuccaro(kb_or_null, &kc, sel, lhs_q, k_in, out_q, out_w, st->temp_anc, &cbit);
                        if (!ip) {
                            for (uint16_t b = 0u; b < out_w; ++b) qk_cx(kb_or_null, &kc, qk__q(in1_q, b), qk__q(out_q, b));
                        }
                    }
                    break;
                }
                for (uint16_t b = 0u; b < out_w; ++b) {
                    uint32_t q1 = qk__q(in1_q, b);
                    uint32_t q2 = qk__q(in2_q, b);
                    uint32_t qo = qk__q(out_q, b);
                    if (q1 == q2) {
                        qk_cx(kb_or_null, &kc, q1, qo);
                        continue;
                    }
                    qk_cx(kb_or_null, &kc, q1, qo);
                    qk_cx(kb_or_null, &kc, q1, q2);
                    qk_ccx(kb_or_null, &kc, sel, q2, qo);
                    qk_cx(kb_or_null, &kc, q1, q2);
                }
                break;
            }
            case QIR_OP_PHASE_POLY: {
                uint32_t flat_q[QIR_MAX_VAL_BITS];
                uint32_t pos = 0u;
                for (uint16_t i = 0u; i < op->num_ins; ++i) {
                    const uint32_t *iq = plan->val_allocs[op->ins[i]].qubits;
                    uint16_t iw = m->vals[op->ins[i]].width;
                    for (uint16_t b = 0u; b < iw && pos < QIR_MAX_VAL_BITS; ++b, ++pos) {
                        flat_q[pos] = qk__q(iq, b);
                    }
                }
                if (st->variant == QIR_VAR_PHASE_CSA_DIRECT) {
                    const qir_phase_split_t *ps = st->phase_split;
                    uint16_t bit_lo = 0u, sum_w = 0u;
                    uint32_t leaves[16];
                    int8_t   lsigns[16];
                    uint16_t n_leaves = 0u;
                    if (ps) {
                        bit_lo   = ps->bit_lo;
                        sum_w    = ps->width;
                        n_leaves = ps->num_leaves;
                        for (uint16_t li = 0u; li < n_leaves; ++li) {
                            leaves[li] = ps->leaf_vals[li];
                            lsigns[li] = ps->leaf_signs[li];
                        }
                    } else {
                        const qir_op_t *dead_op = &m->ops[st->aux_param];
                        uint32_t sum_vid = dead_op->outs[0];
                        sum_w = m->vals[sum_vid].width;
                        for (uint16_t i = 0u; i < op->num_ins; ++i) {
                            if (op->ins[i] == sum_vid) break;
                            bit_lo = (uint16_t)(bit_lo + m->vals[op->ins[i]].width);
                        }
                        leaves[0] = dead_op->ins[0]; lsigns[0] = 1;
                        leaves[1] = dead_op->ins[1]; lsigns[1] = (dead_op->kind == QIR_OP_SUB) ? -1 : 1;
                        n_leaves = 2u;
                    }
                    for (uint32_t t = 0u; t < op->u.phase_poly.num_terms; ++t) {
                        const qir_phase_term_t *pt = &op->u.phase_poly.terms[t];
                        int32_t hit_b = -1;
                        for (uint16_t b = 0u; b < sum_w; ++b) {
                            if (qcc_bv_get(pt->mask, (uint32_t)(bit_lo + b))) { hit_b = b; break; }
                        }
                        if (hit_b >= 0) {
                            flow_u128 ang = uncomp ? (flow_u128)(0u - pt->angle) : pt->angle;
                            for (uint16_t li = 0u; li < n_leaves; ++li) {
                                flow_u128 lang = (lsigns[li] < 0) ? (flow_u128)(0u - ang) : ang;
                                const uint32_t *lq = plan->val_allocs[leaves[li]].qubits;
                                if (lq) qk_z_pow_u128(kb_or_null, &kc, lq[hit_b], lang);
                            }
                        } else {
                            qir_phase_term_t pt_eff = *pt;
                            if (uncomp) pt_eff.angle = (flow_u128)(0u - pt_eff.angle);
                            qk_emit_phase_terms(kb_or_null, &kc, flat_q, op->u.phase_poly.total_in_bits,
                                                &pt_eff, 1u, st->temp_anc, &cbit);
                        }
                    }
                    break;
                }
                if (!uncomp) {
                    qk_emit_phase_terms(kb_or_null, &kc, flat_q, op->u.phase_poly.total_in_bits,
                                        op->u.phase_poly.terms, op->u.phase_poly.num_terms, st->temp_anc, &cbit);
                } else {
                    for (uint32_t t = 0u; t < op->u.phase_poly.num_terms; ++t) {
                        qir_phase_term_t pt_eff = op->u.phase_poly.terms[t];
                        pt_eff.angle = (flow_u128)(0u - pt_eff.angle);
                        qk_emit_phase_terms(kb_or_null, &kc, flat_q, op->u.phase_poly.total_in_bits,
                                            &pt_eff, 1u, st->temp_anc, &cbit);
                    }
                }
                break;
            }
            case QIR_OP_HMR_MEAS: {
                if (op->num_outs > 0u) {
                    const uint32_t *out_b = plan->val_allocs[op->outs[0]].cbits;
                    uint16_t iw = (op->num_ins > 0u) ? m->vals[op->ins[0]].width : out_w;
                    qk_hmr_meas(kb_or_null, &kc, in0_q, out_b, iw);
                }
                break;
            }
            case QIR_OP_BIT_OP: {
                if (op->num_outs > 0u) {
                    const uint32_t *a_b   = (op->num_ins > 0u) ? plan->val_allocs[op->ins[0]].cbits : NULL;
                    const uint32_t *b_b   = (op->num_ins > 1u) ? plan->val_allocs[op->ins[1]].cbits : NULL;
                    const uint32_t *out_b = plan->val_allocs[op->outs[0]].cbits;
                    uint32_t bw = (op->u.bit_op.mode == 3u && op->num_ins > 0u) ? m->vals[op->ins[0]].width : out_w;
                    qk_bit_op(kb_or_null, &kc, op->u.bit_op.mode, a_b, b_b, out_b, bw);
                }
                break;
            }
            case QIR_OP_SHL_VAR: {
                if (op->num_outs > 0u && op->num_ins >= 2u) {
                    const uint32_t *aux_q = plan->val_allocs[op->outs[0]].aux_qubits;
                    qk_shl_var(kb_or_null, &kc,
                               in0_q, m->vals[op->ins[0]].width,
                               in1_q, m->vals[op->ins[1]].width,
                               out_q, out_w, aux_q, op->u.shift.max_shift, uncomp);
                }
                break;
            }
            case QIR_OP_SHR_VAR: {
                if (op->num_outs > 0u && op->num_ins >= 2u) {
                    const uint32_t *aux_q = plan->val_allocs[op->outs[0]].aux_qubits;
                    qk_shr_var(kb_or_null, &kc,
                               in0_q, m->vals[op->ins[0]].width,
                               in1_q, m->vals[op->ins[1]].width,
                               out_q, out_w, aux_q, op->u.shift.max_shift, uncomp);
                }
                break;
            }
            case QIR_OP_EQ_ZERO: {
                if (op->num_outs > 0u && op->num_ins >= 1u) {
                    uint16_t iw = m->vals[op->ins[0]].width;
                    if (!uncomp) qk_eq_zero_compute(kb_or_null, &kc, in0_q, iw, out_q, st->temp_anc, &cbit);
                    else         qk_eq_zero_uncompute(kb_or_null, &kc, in0_q, iw, out_q, st->temp_anc, &cbit);
                }
                break;
            }
            case QIR_OP_PHASE_DOT: {
                if (op->num_ins == 2u) {
                    const uint32_t *cb = plan->val_allocs[op->ins[1]].cbits;
                    qk_phase_dot(kb_or_null, &kc, in0_q, cb, m->vals[op->ins[0]].width);
                } else if (op->num_ins >= 3u) {
                    const uint32_t *cb = plan->val_allocs[op->ins[2]].cbits;
                    qk_phase_dot_gf2m(kb_or_null, &kc, in0_q, in1_q, cb,
                                      op->u.gf2m.m, op->u.gf2m.irred_poly,
                                      (plan->has_gf12_r42 && plan->gf12_r42_poly == op->u.gf2m.irred_poly) ? &plan->gf12_r42 : NULL,
                                      (plan->has_gf256_r27 && plan->gf256_r27_poly == op->u.gf2m.irred_poly) ? &plan->gf256_r27 : NULL,
                                      (plan->has_gf64_r15 && plan->gf64_r15_poly == op->u.gf2m.irred_poly) ? &plan->gf64_r15 : NULL,
                                      &cbit);
                }
                break;
            }
            case QIR_OP_OPAQUE_KMB: {
                const qcc_kmb_circuit_t *sk = op->u.opaque.kmb;
                if (!sk) break;
                uint32_t sub_nq = op->u.opaque.num_q ? op->u.opaque.num_q : sk->num_qubits;
                uint32_t *qmap = (uint32_t *)qcc_xmalloc((size_t)(sub_nq ? sub_nq : 1u) * sizeof(uint32_t));
                for (uint32_t q = 0u; q < sub_nq; ++q) qmap[q] = q;
                if (st->inplace_on_in0) {
                    uint32_t pos = 0u;
                    for (uint16_t i = 0u; i < op->num_ins; ++i) {
                        const uint32_t *iq = plan->val_allocs[op->ins[i]].qubits;
                        uint16_t iw = m->vals[op->ins[i]].width;
                        for (uint16_t b = 0u; b < iw && pos < op->u.opaque.in_bits; ++b, ++pos) {
                            uint32_t sq = op->u.opaque.in_q_idx[pos];
                            if (sq < sub_nq) qmap[sq] = iq[b];
                        }
                    }
                    for (uint16_t b = 0u; b < op->u.opaque.out_bits && out_q; ++b) {
                        uint32_t sq = op->u.opaque.out_q_idx[b];
                        if (sq < sub_nq) qmap[sq] = out_q[b];
                    }
                } else {
                    /* Non-inplace: bind output slots (copied clean from inputs first) */
                    for (uint16_t b = 0u; b < op->u.opaque.out_bits && out_q; ++b) {
                        uint32_t sq = op->u.opaque.out_q_idx[b];
                        if (sq < sub_nq) qmap[sq] = out_q[b];
                    }
                    uint32_t pos = 0u;
                    for (uint16_t i = 0u; i < op->num_ins; ++i) {
                        const uint32_t *iq = plan->val_allocs[op->ins[i]].qubits;
                        uint16_t iw = m->vals[op->ins[i]].width;
                        for (uint16_t b = 0u; b < iw && pos < op->u.opaque.in_bits; ++b, ++pos) {
                            uint32_t sq = op->u.opaque.in_q_idx[pos];
                            if (sq < sub_nq) {
                                /* If output also uses sq, copy iq[b] -> qmap[sq] rather than aliasing! */
                                bool shared_with_out = false;
                                for (uint16_t ob = 0u; ob < op->u.opaque.out_bits; ++ob) {
                                    if (op->u.opaque.out_q_idx[ob] == sq) { shared_with_out = true; break; }
                                }
                                if (shared_with_out) qk_cx(kb_or_null, &kc, iq[b], qmap[sq]);
                                else                 qmap[sq] = iq[b];
                            }
                        }
                    }
                }
                qcc_kmb_cursor_t cur;
                qcc_kmb_cursor_init(&cur);
                qcc_kmb_op_view_t v;
                while (qcc_kmb_cursor_next(sk, &cur, &v)) {
                    for (uint8_t k = 0u; k < v.n_q; ++k) {
                        if (v.q[k] < sub_nq) v.q[k] = qmap[v.q[k]];
                        if (v.q[k] + 1u > max_q) max_q = v.q[k] + 1u;
                    }
                    if (v.has_b) {
                        v.b += st->cbit_base;
                        if (v.b + 1u > max_b) max_b = v.b + 1u;
                    }
                    if (v.has_cond) {
                        v.cond += st->cbit_base;
                        if (v.cond + 1u > max_b) max_b = v.cond + 1u;
                    }
                    if (kb_or_null) {
                        qcc_kmb_emit_view(kb_or_null, &v);
                    } else {
                        uint8_t base = v.op & (uint8_t)(~QCC_OP_FLAG_COND);
                        if (base == QCC_OP_CCX || base == QCC_OP_CCZ || base == QCC_OP_CSWAP) kc.toffoli++;
                        else if (base == QCC_OP_Z_POW) {
                            kc.z_pow++;
                            if (qcc_u128_clifford_level(qcc_u128_from_angle(v.angle)) == 3u) kc.t_like++;
                        } else if (base == QCC_OP_HMR) kc.hmr++;
                        else if (base == QCC_OP_CX) kc.cnot++;
                        else if (base == QCC_OP_CZ) kc.cz++;
                        else if (base == QCC_OP_R) kc.r_count++;
                    }
                }
                free(qmap);
                break;
            }
            default:
                break;
        }
        if (has_pred) {
            qk_pop_cond(kb_or_null, &kc);
        }
        if (cbit > max_b) max_b = cbit;
    }
    free(folded_op);

    if (out_stats) {
        if (kb_or_null) {
            if (!qcc_kmb_recompute_stats(kb_or_null, out_stats)) return false;
        } else {
            memset(out_stats, 0, sizeof(*out_stats));
            out_stats->num_qubits = max_q;
            out_stats->num_bits   = max_b;
            out_stats->toffoli    = kc.toffoli;
            out_stats->z_pow      = kc.z_pow;
            out_stats->t_like     = kc.t_like;
            out_stats->hmr        = kc.hmr;
            out_stats->cnot       = kc.cnot;
            out_stats->cz         = kc.cz;
            out_stats->r_count    = kc.r_count;
        }
    }
    return true;
}

static inline bool qir_codegen_emit(const qir_phys_plan_t *plan, qcc_kmb_circuit_t *out_kmb) {
    if (!plan || !out_kmb) return false;
    qcc_kmb_init(out_kmb);
    const qir_module_t *m = plan->mod;

    /* Declare quantum and classical I/O registers cleanly */
    for (uint32_t i = 0u; i < m->num_inputs; ++i) {
        uint32_t vid = m->inputs[i];
        const qir_val_alloc_t *va = &plan->val_allocs[vid];
        if (va->width > 0u && va->qubits) {
            char rname[64];
            snprintf(rname, sizeof(rname), "%s", m->vals[vid].name ? m->vals[vid].name : "in");
            qcc_kmb_add_quantum_reg(out_kmb, rname, va->qubits, va->width);
        } else if (va->cwidth > 0u && va->cbits) {
            char rname[64];
            snprintf(rname, sizeof(rname), "%s", m->vals[vid].name ? m->vals[vid].name : "cin");
            qcc_kmb_add_classical_reg(out_kmb, rname, va->cbits, va->cwidth);
        }
    }
    for (uint32_t i = 0u; i < m->num_outputs; ++i) {
        uint32_t vid = m->outputs[i];
        const qir_val_alloc_t *va = &plan->val_allocs[vid];
        if (va->width > 0u && va->qubits) {
            const char *oname = m->vals[vid].name ? m->vals[vid].name : "out";
            bool same_as_in = false;
            for (uint32_t j = 0u; j < m->num_inputs; ++j) {
                uint32_t iv = m->inputs[j];
                if (iv == vid) { same_as_in = true; break; }
                const qir_val_alloc_t *iva = &plan->val_allocs[iv];
                const char *iname = m->vals[iv].name ? m->vals[iv].name : "in";
                if (iva->width == va->width && iva->qubits &&
                    memcmp(iva->qubits, va->qubits, (size_t)va->width * sizeof(uint32_t)) == 0 &&
                    strcmp(oname, iname) == 0) {
                    same_as_in = true;
                    break;
                }
            }
            if (!same_as_in) {
                char rname[64];
                snprintf(rname, sizeof(rname), "%s", oname);
                qcc_kmb_add_quantum_reg(out_kmb, rname, va->qubits, va->width);
            }
        }
    }
    for (uint32_t i = 0u; i < m->num_classical_outputs; ++i) {
        uint32_t vid = m->classical_outputs[i];
        if (vid < m->num_vals) {
            const qir_val_alloc_t *va = &plan->val_allocs[vid];
            if (va->cwidth > 0u && va->cbits) {
                char rname[64];
                snprintf(rname, sizeof(rname), "%s", m->vals[vid].name ? m->vals[vid].name : "cout");
                qcc_kmb_add_classical_reg(out_kmb, rname, va->cbits, va->cwidth);
            }
        }
    }

    qcc_kmb_stats_t actual;
    if (!qir_codegen_or_count(plan, out_kmb, &actual)) return false;

    /* Hard Zero-Discrepancy Compiler Invariant Assertions */
    if (actual.toffoli != plan->est_toffoli ||
        actual.z_pow != plan->est_z_pow ||
        actual.t_like != plan->est_t_like ||
        actual.hmr != plan->est_hmr ||
        actual.num_qubits > plan->est_qubits) {
        fprintf(stderr,
                "[QCC CODEGEN INVARIANT VIOLATION] est(q=%u,ccx=%llu,zpow=%llu,hmr=%llu) != actual(q=%u,ccx=%llu,zpow=%llu,hmr=%llu)\n",
                (unsigned)plan->est_qubits, (unsigned long long)plan->est_toffoli,
                (unsigned long long)plan->est_z_pow, (unsigned long long)plan->est_hmr,
                (unsigned)actual.num_qubits, (unsigned long long)actual.toffoli,
                (unsigned long long)actual.z_pow, (unsigned long long)actual.hmr);
        return false;
    }

    if (plan->emit_clean_reg && plan->clean_scratch_certified && actual.num_qubits > 0u) {
        uint8_t *io_q = (uint8_t *)qcc_xcalloc(actual.num_qubits, 1u);
        for (uint32_t r = 0u; r < out_kmb->num_regs; ++r) {
            const qcc_kmb_reg_t *reg = &out_kmb->regs[r];
            for (uint32_t k = 0u; k < reg->len; ++k) {
                if ((reg->targets[k] & UINT32_C(0x80000000)) != 0u) {
                    uint32_t qid = reg->targets[k] & UINT32_C(0x7FFFFFFF);
                    if (qid < actual.num_qubits) io_q[qid] = 1u;
                }
            }
        }
        uint32_t *clean_qs = (uint32_t *)qcc_xmalloc(actual.num_qubits * sizeof(uint32_t));
        uint32_t n_clean = 0u;
        for (uint32_t q = 0u; q < actual.num_qubits; ++q) {
            if (!io_q[q]) clean_qs[n_clean++] = q;
        }
        if (n_clean > 0u) qcc_kmb_add_quantum_reg(out_kmb, "@clean", clean_qs, n_clean);
        free(io_q);
        free(clean_qs);
    }
    return true;
}

#ifdef __cplusplus
}
#endif

#endif /* QIR_CODEGEN_H */
