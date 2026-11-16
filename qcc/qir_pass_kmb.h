/*
 * qir_pass_kmb.h — Post-Codegen Phase-Poly Fusion, Commutation Cancellation & Live-Range Compaction
 * Part of QCC (Quantum Circuit Compiler) for FLOW v7.0.1.
 * Pure C17, zero external dependencies.
 */
#ifndef QIR_PASS_KMB_H
#define QIR_PASS_KMB_H

#include "qir_codegen.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t   ops_before;
    size_t   ops_after;
    uint64_t z_pow_fused_or_cancelled;
    uint64_t self_inv_cancelled;
    uint32_t qubits_before;
    uint32_t qubits_after;
} qir_kmb_pass_stats_t;

static inline bool qir__kmb_views_equal(const qcc_kmb_op_view_t *a, const qcc_kmb_op_view_t *b) {
    if (a->op != b->op || a->n_q != b->n_q || a->has_b != b->has_b || a->has_cond != b->has_cond) return false;
    for (uint8_t i = 0u; i < a->n_q; ++i) if (a->q[i] != b->q[i]) return false;
    if (a->has_b && a->b != b->b) return false;
    if (a->has_cond && a->cond != b->cond) return false;
    if (a->has_angle) {
        if (a->angle.lo != b->angle.lo || a->angle.hi != b->angle.hi) return false;
    }
    return true;
}

static inline void qir__kmb_canonicalize_view(qcc_kmb_op_view_t *v) {
    uint8_t base = v->op & (uint8_t)(~QCC_OP_FLAG_COND);
    if (base == QCC_OP_CZ || base == QCC_OP_SWAP || base == QCC_OP_CCX) {
        if (v->q[0] > v->q[1]) {
            uint32_t t = v->q[0]; v->q[0] = v->q[1]; v->q[1] = t;
        }
    } else if (base == QCC_OP_CCZ) {
        if (v->q[0] > v->q[1]) { uint32_t t = v->q[0]; v->q[0] = v->q[1]; v->q[1] = t; }
        if (v->q[1] > v->q[2]) { uint32_t t = v->q[1]; v->q[1] = v->q[2]; v->q[2] = t; }
        if (v->q[0] > v->q[1]) { uint32_t t = v->q[0]; v->q[0] = v->q[1]; v->q[1] = t; }
    } else if (base == QCC_OP_CSWAP) {
        if (v->q[1] > v->q[2]) {
            uint32_t t = v->q[1]; v->q[1] = v->q[2]; v->q[2] = t;
        }
    }
}

static inline bool qir__kmb_is_self_inverse(const qcc_kmb_op_view_t *v) {
    if (v->has_cond || v->has_b || v->has_angle || v->has_info_tag) return false;
    switch (v->op) {
        case QCC_OP_X:
        case QCC_OP_Z:
        case QCC_OP_CX:
        case QCC_OP_CZ:
        case QCC_OP_SWAP:
        case QCC_OP_CCX:
        case QCC_OP_CCZ:
        case QCC_OP_CSWAP:
            return true;
        default:
            return false;
    }
}

static inline bool qir__kmb_is_single_z(const qcc_kmb_op_view_t *v, uint32_t *out_q, flow_u128 *out_ang) {
    if (v->has_cond || v->has_b || v->has_info_tag) return false;
    if (v->op == QCC_OP_Z) {
        *out_q = v->q[0];
        *out_ang = (flow_u128)1u << 127u;
        return true;
    }
    if (v->op == QCC_OP_Z_POW) {
        *out_q = v->q[0];
        *out_ang = qcc_u128_from_angle(v->angle);
        return true;
    }
    return false;
}

static inline bool qir__kmb_commutes_with_zq(const qcc_kmb_op_view_t *v, uint32_t q) {
    if (v->op == QCC_OP_PUSH_CONDITION || v->op == QCC_OP_POP_CONDITION) return false;
    if (v->op == QCC_OP_INFO_TAG || v->op == QCC_OP_DEBUG_PRINT_EMPTY) return true;
    uint8_t base = v->op & (uint8_t)(~QCC_OP_FLAG_COND);
    switch (base) {
        case QCC_OP_Z:
        case QCC_OP_Z_POW:
        case QCC_OP_CZ:
        case QCC_OP_CCZ:
        case QCC_OP_NEG:
        case QCC_OP_BIT_STORE0:
        case QCC_OP_BIT_STORE1:
        case QCC_OP_BIT_INVERT:
            return true;
        case QCC_OP_X:
        case QCC_OP_R:
        case QCC_OP_HMR:
            return v->q[0] != q;
        case QCC_OP_CX:
            return v->q[1] != q;
        case QCC_OP_SWAP:
            return (v->q[0] != q && v->q[1] != q);
        case QCC_OP_CCX:
            return v->q[2] != q;
        case QCC_OP_CSWAP:
            return (v->q[1] != q && v->q[2] != q);
        default:
            return false;
    }
}

typedef struct {
    uint32_t gen[8];
    uint8_t  len;
    uint8_t  c;
} qir__kmb_parity_t;

typedef struct {
    qir__kmb_parity_t key;
    size_t            op_idx;
    flow_u128         angle;
    uint8_t           valid;
} qir__kmb_zslot_t;

static inline bool qir__kmb_parity_eq(const qir__kmb_parity_t *a, const qir__kmb_parity_t *b) {
    if (a->len != b->len || a->c != b->c) return false;
    for (uint8_t i = 0u; i < a->len; ++i) {
        if (a->gen[i] != b->gen[i]) return false;
    }
    return true;
}

static inline void qir__kmb_set_z_view(qcc_kmb_op_view_t *v, uint32_t q, flow_u128 ang) {
    if (ang == ((flow_u128)1u << 127u)) {
        v->op        = QCC_OP_Z;
        v->q[0]      = q;
        v->has_angle = 0u;
        v->angle     = qcc_angle_from_u128(0u);
    } else {
        v->op        = QCC_OP_Z_POW;
        v->q[0]      = q;
        v->has_angle = 1u;
        v->angle     = qcc_angle_from_u128(ang);
    }
}

/* Peephole Canonicalization, Parity-Tracked Z/Z_POW Fusion, O(n) Horizon Cancellation & Qubit Compaction */
static inline bool qir_pass_kmb_optimize_with_plan(
    qcc_kmb_circuit_t *c,
    qir_phys_plan_t *plan_or_null,
    qir_kmb_pass_stats_t *out_stats
) {
    if (!c) return false;
    qir_kmb_pass_stats_t st;
    memset(&st, 0, sizeof(st));
    st.ops_before = c->num_ops;

    qcc_kmb_stats_t ks0;
    if (!qcc_kmb_recompute_stats(c, &ks0)) return false;
    st.qubits_before = ks0.num_qubits;

    size_t n = c->num_ops;
    if (n == 0u) {
        st.qubits_after = ks0.num_qubits;
        if (out_stats) *out_stats = st;
        return true;
    }

    qcc_kmb_op_view_t *views = (qcc_kmb_op_view_t *)calloc(n, sizeof(qcc_kmb_op_view_t));
    uint8_t           *alive = (uint8_t *)malloc(n);
    if (!views || !alive) { free(views); free(alive); return false; }
    memset(alive, 1, n);

    qcc_kmb_cursor_t cur;
    qcc_kmb_cursor_init(&cur);
    uint32_t max_q_seen = ks0.num_qubits;
    for (size_t i = 0u; i < n; ++i) {
        if (!qcc_kmb_cursor_next(c, &cur, &views[i])) {
            free(views); free(alive); return false;
        }
        /* Pass 0: Symmetric Gate Canonicalization */
        qir__kmb_canonicalize_view(&views[i]);
        for (uint8_t k = 0u; k < views[i].n_q; ++k) {
            if (views[i].q[k] + 1u > max_q_seen) max_q_seen = views[i].q[k] + 1u;
        }
    }
    uint32_t old_nq = max_q_seen ? max_q_seen : 1u;

    /* Pass A1: Parity-Tracked Z/Z_POW Fusion across linear Clifford/X/CX/SWAP networks */
    qir__kmb_parity_t *basis = (qir__kmb_parity_t *)calloc(old_nq, sizeof(qir__kmb_parity_t));
    qir__kmb_zslot_t   zslots[256];
    memset(zslots, 0, sizeof(zslots));
    uint32_t next_gen = 1u;
    for (uint32_t q = 0u; q < old_nq; ++q) {
        basis[q].gen[0] = next_gen++;
        basis[q].len    = 1u;
        basis[q].c      = 0u;
    }
    uint32_t z_cursor = 0u;
    size_t cond_depth = 0u;

    for (size_t i = 0u; i < n; ++i) {
        if (!alive[i]) continue;
        const qcc_kmb_op_view_t *v = &views[i];
        if (v->op == QCC_OP_PUSH_CONDITION) { cond_depth++; continue; }
        if (v->op == QCC_OP_POP_CONDITION)  { if (cond_depth > 0u) cond_depth--; continue; }

        uint32_t zq = 0u;
        flow_u128 zang = 0u;
        if (cond_depth == 0u && qir__kmb_is_single_z(v, &zq, &zang)) {
            if (basis[zq].len == 0u && basis[zq].c == 0u) {
                /* Provably |0>: diagonal Z_POW is global phase 1 */
                alive[i] = 0u;
                st.z_pow_fused_or_cancelled++;
                continue;
            }
            int32_t hit = -1;
            for (int32_t s = 0; s < 256; ++s) {
                if (zslots[s].valid && alive[zslots[s].op_idx] &&
                    qir__kmb_parity_eq(&zslots[s].key, &basis[zq])) {
                    hit = s;
                    break;
                }
            }
            if (hit >= 0) {
                size_t first = zslots[hit].op_idx;
                flow_u128 sum = zslots[hit].angle + zang;
                alive[i] = 0u;
                st.z_pow_fused_or_cancelled++;
                if (sum == 0u) {
                    alive[first] = 0u;
                    zslots[hit].valid = 0u;
                    st.z_pow_fused_or_cancelled++;
                } else {
                    zslots[hit].angle = sum;
                    qir__kmb_set_z_view(&views[first], views[first].q[0], sum);
                }
            } else {
                int32_t dst_slot = -1;
                for (int32_t s = 0; s < 256; ++s) {
                    if (!zslots[s].valid || !alive[zslots[s].op_idx]) { dst_slot = s; break; }
                }
                if (dst_slot < 0) {
                    dst_slot = (int32_t)(z_cursor & 255u);
                    z_cursor++;
                }
                zslots[dst_slot].valid  = 1u;
                zslots[dst_slot].key    = basis[zq];
                zslots[dst_slot].op_idx = i;
                zslots[dst_slot].angle  = zang;
            }
            continue;
        }

        uint8_t base = v->op & (uint8_t)(~QCC_OP_FLAG_COND);
        bool uncond = (cond_depth == 0u && !v->has_cond);
        if (uncond && base == QCC_OP_X) {
            basis[v->q[0]].c ^= 1u;
        } else if (uncond && base == QCC_OP_SWAP) {
            qir__kmb_parity_t tmp = basis[v->q[0]];
            basis[v->q[0]] = basis[v->q[1]];
            basis[v->q[1]] = tmp;
        } else if (uncond && base == QCC_OP_CX) {
            const qir__kmb_parity_t *src = &basis[v->q[0]];
            qir__kmb_parity_t       *dst = &basis[v->q[1]];
            uint32_t merged[16];
            uint8_t ml = 0u, ia = 0u, ib = 0u;
            while (ia < dst->len && ib < src->len) {
                if (dst->gen[ia] < src->gen[ib])      merged[ml++] = dst->gen[ia++];
                else if (src->gen[ib] < dst->gen[ia]) merged[ml++] = src->gen[ib++];
                else { ia++; ib++; }
            }
            while (ia < dst->len) merged[ml++] = dst->gen[ia++];
            while (ib < src->len) merged[ml++] = src->gen[ib++];
            if (ml <= 8u) {
                dst->len = ml;
                for (uint8_t k = 0u; k < ml; ++k) dst->gen[k] = merged[k];
                dst->c ^= src->c;
            } else {
                dst->len = 1u; dst->gen[0] = next_gen++; dst->c = 0u;
            }
        } else if (uncond && (base == QCC_OP_R || base == QCC_OP_HMR)) {
            basis[v->q[0]].len = 0u;
            basis[v->q[0]].c   = 0u;
        } else if (base == QCC_OP_CCX) {
            basis[v->q[2]].len = 1u; basis[v->q[2]].gen[0] = next_gen++; basis[v->q[2]].c = 0u;
        } else if (base == QCC_OP_CSWAP) {
            basis[v->q[1]].len = 1u; basis[v->q[1]].gen[0] = next_gen++; basis[v->q[1]].c = 0u;
            basis[v->q[2]].len = 1u; basis[v->q[2]].gen[0] = next_gen++; basis[v->q[2]].c = 0u;
        } else if (!uncond) {
            if (base == QCC_OP_X || base == QCC_OP_R || base == QCC_OP_HMR) {
                basis[v->q[0]].len = 1u; basis[v->q[0]].gen[0] = next_gen++; basis[v->q[0]].c = 0u;
            } else if (base == QCC_OP_CX) {
                basis[v->q[1]].len = 1u; basis[v->q[1]].gen[0] = next_gen++; basis[v->q[1]].c = 0u;
            } else if (base == QCC_OP_SWAP) {
                basis[v->q[0]].len = 1u; basis[v->q[0]].gen[0] = next_gen++; basis[v->q[0]].c = 0u;
                basis[v->q[1]].len = 1u; basis[v->q[1]].gen[0] = next_gen++; basis[v->q[1]].c = 0u;
            }
        }
    }
    free(basis);

    /* Pass A2: Local commuting-control Z/Z_POW cleanup across 64-op lookahead */
    for (size_t i = 0u; i < n; ++i) {
        if (!alive[i]) continue;
        uint32_t q = 0u;
        flow_u128 ang1 = 0u;
        if (!qir__kmb_is_single_z(&views[i], &q, &ang1)) continue;
        size_t steps = 0u;
        for (size_t j = i + 1u; j < n && steps < 64u; ++j) {
            if (!alive[j]) continue;
            uint32_t q2 = 0u;
            flow_u128 ang2 = 0u;
            if (qir__kmb_is_single_z(&views[j], &q2, &ang2) && q2 == q) {
                flow_u128 sum = ang1 + ang2;
                alive[j] = 0u;
                st.z_pow_fused_or_cancelled++;
                ang1 = sum;
                if (sum == 0u) {
                    alive[i] = 0u;
                    st.z_pow_fused_or_cancelled++;
                    break;
                }
                qir__kmb_set_z_view(&views[i], q, sum);
                continue;
            }
            if (!qir__kmb_commutes_with_zq(&views[j], q)) break;
            steps++;
        }
    }

    /* Pass B: O(n) Per-Qubit Horizon Self-Inverse Cancellation */
    size_t *last_q = (size_t *)malloc((size_t)old_nq * sizeof(size_t));
    size_t (*prev_q)[3] = (size_t (*)[3])malloc(n * sizeof(*prev_q));
    for (uint32_t q = 0u; q < old_nq; ++q) last_q[q] = SIZE_MAX;
    for (size_t i = 0u; i < n; ++i) prev_q[i][0] = prev_q[i][1] = prev_q[i][2] = SIZE_MAX;

    cond_depth = 0u;
    for (size_t i = 0u; i < n; ++i) {
        if (!alive[i]) continue;
        const qcc_kmb_op_view_t *v = &views[i];
        if (v->op == QCC_OP_PUSH_CONDITION) { cond_depth++; continue; }
        if (v->op == QCC_OP_POP_CONDITION)  { if (cond_depth > 0u) cond_depth--; continue; }
        if (v->op == QCC_OP_INFO_TAG || v->n_q == 0u) continue;

        if (cond_depth == 0u && qir__kmb_is_self_inverse(v)) {
            size_t p = last_q[v->q[0]];
            bool same_horizon = (p != SIZE_MAX && alive[p]);
            for (uint8_t k = 1u; k < v->n_q && same_horizon; ++k) {
                if (last_q[v->q[k]] != p) same_horizon = false;
            }
            if (same_horizon && qir__kmb_views_equal(&views[p], v)) {
                alive[p] = 0u;
                alive[i] = 0u;
                st.self_inv_cancelled += 2u;
                for (uint8_t k = 0u; k < views[p].n_q; ++k) {
                    last_q[views[p].q[k]] = prev_q[p][k];
                }
                continue;
            }
        }
        for (uint8_t k = 0u; k < v->n_q; ++k) {
            prev_q[i][k] = last_q[v->q[k]];
            last_q[v->q[k]] = i;
        }
    }
    free(last_q);
    free(prev_q);

    /* Pass C: Physical Qubit Hole Compaction (syncing c->regs and plan_or_null) */
    uint8_t  *q_used  = (uint8_t *)calloc(old_nq, 1u);
    uint32_t *q_remap = (uint32_t *)calloc(old_nq, sizeof(uint32_t));
    for (uint32_t r = 0u; r < c->num_regs; ++r) {
        const qcc_kmb_reg_t *reg = &c->regs[r];
        if (reg->len > 0u && ((reg->targets[0] & UINT32_C(0x80000000)) != 0u)) {
            for (uint32_t k = 0u; k < reg->len; ++k) {
                uint32_t qid = reg->targets[k] & UINT32_C(0x7FFFFFFF);
                if (qid < old_nq) q_used[qid] = 1u;
            }
        }
    }
    for (size_t i = 0u; i < n; ++i) {
        if (!alive[i]) continue;
        for (uint8_t k = 0u; k < views[i].n_q; ++k) {
            if (views[i].q[k] < old_nq) q_used[views[i].q[k]] = 1u;
        }
    }
    uint32_t new_nq = 0u;
    for (uint32_t q = 0u; q < old_nq; ++q) {
        if (q_used[q]) q_remap[q] = new_nq++;
    }
    for (size_t i = 0u; i < n; ++i) {
        if (!alive[i]) continue;
        for (uint8_t k = 0u; k < views[i].n_q; ++k) {
            if (views[i].q[k] < old_nq) views[i].q[k] = q_remap[views[i].q[k]];
        }
    }
    if (plan_or_null && plan_or_null->val_allocs && plan_or_null->mod) {
        for (uint32_t v = 0u; v < plan_or_null->mod->num_vals; ++v) {
            qir_val_alloc_t *va = &plan_or_null->val_allocs[v];
            if (va->qubits) {
                for (uint16_t b = 0u; b < va->width; ++b) {
                    if (va->qubits[b] < old_nq) va->qubits[b] = q_remap[va->qubits[b]];
                }
            }
            if (va->aux_qubits) {
                for (uint16_t b = 0u; b < va->aux_width; ++b) {
                    if (va->aux_qubits[b] < old_nq) va->aux_qubits[b] = q_remap[va->aux_qubits[b]];
                }
            }
        }
    }

    /* Rebuild compacted circuit */
    qcc_kmb_circuit_t next_c;
    qcc_kmb_init(&next_c);
    for (uint32_t r = 0u; r < c->num_regs; ++r) {
        const qcc_kmb_reg_t *reg = &c->regs[r];
        bool is_q = (reg->len > 0u) && ((reg->targets[0] & UINT32_C(0x80000000)) != 0u);
        uint32_t *ids = reg->len ? (uint32_t *)malloc(reg->len * sizeof(uint32_t)) : NULL;
        for (uint32_t k = 0u; k < reg->len; ++k) {
            uint32_t tid = reg->targets[k] & UINT32_C(0x7FFFFFFF);
            ids[k] = (is_q && tid < old_nq) ? q_remap[tid] : tid;
        }
        if (is_q) qcc_kmb_add_quantum_reg(&next_c, reg->name, ids, reg->len);
        else      qcc_kmb_add_classical_reg(&next_c, reg->name, ids, reg->len);
        free(ids);
    }
    for (size_t i = 0u; i < n; ++i) {
        if (alive[i]) qcc_kmb_emit_view(&next_c, &views[i]);
    }
    free(q_used);
    free(q_remap);
    free(views);
    free(alive);

    qcc_kmb_stats_t ks1;
    if (!qcc_kmb_recompute_stats(&next_c, &ks1)) {
        qcc_kmb_free(&next_c);
        return false;
    }
    st.ops_after    = next_c.num_ops;
    st.qubits_after = ks1.num_qubits;

    qcc_kmb_free(c);
    *c = next_c;
    if (out_stats) *out_stats = st;
    return true;
}

static inline bool qir_pass_kmb_optimize(qcc_kmb_circuit_t *c, qir_kmb_pass_stats_t *out_stats) {
    return qir_pass_kmb_optimize_with_plan(c, NULL, out_stats);
}

#ifdef __cplusplus
}
#endif

#endif /* QIR_PASS_KMB_H */
