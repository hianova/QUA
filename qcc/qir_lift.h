/*
 * qir_lift.h — Scale-Free .kmb v2 Structural Pattern Lifter to Executable Semantic QIR
 * Part of QCC (Quantum Circuit Compiler) for FLOW v7.0.1.
 * Pure C17, zero external dependencies.
 */
#ifndef QIR_LIFT_H
#define QIR_LIFT_H

#include "qir.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t lifted_linf2;
    uint32_t lifted_phase;
    uint32_t lifted_and_vec;
    uint32_t lifted_opaque;
    uint32_t lifted_ops_total;
} qir_lift_stats_t;

static inline bool qir__is_linf2_op(uint8_t op) {
    return (op == QCC_OP_X || op == QCC_OP_CX || op == QCC_OP_SWAP);
}

static inline bool qir__is_phase_op(uint8_t op) {
    return (op == QCC_OP_Z || op == QCC_OP_CZ || op == QCC_OP_CCZ || op == QCC_OP_Z_POW);
}

static inline bool qir__lift_is_and_compute(const qcc_kmb_op_view_t *views, size_t nops, size_t pos) {
    if (pos + 1u >= nops) return false;
    if (views[pos].op != QCC_OP_R || views[pos + 1u].op != QCC_OP_CCX) return false;
    uint32_t anc = views[pos].q[0];
    uint32_t c0  = views[pos + 1u].q[0];
    uint32_t c1  = views[pos + 1u].q[1];
    return (views[pos + 1u].q[2] == anc && c0 != anc && c1 != anc);
}

static inline bool qir__lift_is_and_uncompute(
    const qir_module_t *m,
    const uint32_t *q_wire,
    const qcc_kmb_op_view_t *views,
    size_t nops,
    size_t pos
) {
    if (pos + 1u >= nops) return false;
    if (views[pos].op != QCC_OP_HMR || views[pos + 1u].op != QCC_OP_CZ_IF) return false;
    uint32_t anc = views[pos].q[0];
    uint32_t b   = views[pos].b;
    uint32_t q0  = views[pos + 1u].q[0];
    uint32_t q1  = views[pos + 1u].q[1];
    if (views[pos + 1u].cond != b || q0 == anc || q1 == anc || q0 == q1) return false;

    /* Ensure classical scratch bit `b` is never read again before being overwritten */
    for (size_t p = pos + 2u; p < nops; ++p) {
        if (views[p].has_cond && views[p].cond == b) return false;
        if (views[p].has_b && views[p].b == b) {
            uint8_t base = views[p].op & (uint8_t)(~QCC_OP_FLAG_COND);
            if (base == QCC_OP_BIT_STORE0 || base == QCC_OP_BIT_STORE1 || base == QCC_OP_HMR) break;
            return false;
        }
    }

    uint32_t w_anc = q_wire[anc];
    uint32_t def   = m->vals[w_anc].def_op;
    if (def == UINT32_MAX) return false;
    const qir_op_t *op = &m->ops[def];
    if (op->kind != QIR_OP_AND_VEC || op->num_ins != 2u) return false;
    return ((op->ins[0] == q_wire[q0] && op->ins[1] == q_wire[q1]) ||
            (op->ins[0] == q_wire[q1] && op->ins[1] == q_wire[q0]));
}

static inline uint32_t qir__lift_find_q(const uint32_t *arr, uint32_t n, uint32_t q) {
    uint32_t lo = 0u, hi = n;
    while (lo < hi) {
        uint32_t mid = lo + ((hi - lo) >> 1u);
        if (arr[mid] < q) lo = mid + 1u;
        else hi = mid;
    }
    return (lo < n && arr[lo] == q) ? lo : UINT32_MAX;
}

static inline bool qir__lift_try_add_op_qubits(
    const qcc_kmb_op_view_t *v,
    uint32_t act_q[QIR_MAX_VAL_BITS],
    uint32_t *k_len
) {
    uint32_t tmp[3], nt = 0u;
    for (uint8_t i = 0u; i < v->n_q; ++i) {
        uint32_t q = v->q[i];
        if (qir__lift_find_q(act_q, *k_len, q) == UINT32_MAX) {
            bool dup = false;
            for (uint32_t j = 0u; j < nt; ++j) if (tmp[j] == q) dup = true;
            if (!dup) tmp[nt++] = q;
        }
    }
    if (*k_len + nt > QIR_MAX_VAL_BITS) return false;
    for (uint32_t t = 0u; t < nt; ++t) {
        uint32_t q = tmp[t];
        uint32_t pos = *k_len;
        while (pos > 0u && act_q[pos - 1u] > q) {
            act_q[pos] = act_q[pos - 1u];
            pos--;
        }
        act_q[pos] = q;
        (*k_len)++;
    }
    return true;
}

static inline uint32_t qir__lift_gather_bus(
    qir_module_t *m, const uint32_t *q_wire, const uint32_t *act_q, uint32_t k
) {
    if (k == 1u) return q_wire[act_q[0]];
    uint32_t parts[QIR_MAX_VAL_BITS];
    for (uint32_t i = 0u; i < k; ++i) parts[i] = q_wire[act_q[i]];
    return qir_emit_concat(m, parts, (uint16_t)k, NULL);
}

/* Scale-free lossless .kmb -> SSA QIR lifter handling arbitrary qubit counts (> 256 supported via chunked buses). */
static inline bool qir_lift_from_kmb(
    qir_module_t *m,
    qcc_arena_t *arena,
    const qcc_kmb_circuit_t *k,
    const char *mod_name,
    qir_lift_stats_t *out_stats
) {
    if (!m || !arena || !k) return false;
    qir_module_init(m, arena, mod_name ? mod_name : "lifted_kmb");
    if (out_stats) memset(out_stats, 0, sizeof(*out_stats));

    qcc_kmb_stats_t ks;
    if (!qcc_kmb_recompute_stats(k, &ks)) return false;

    uint32_t nq = ks.num_qubits > k->num_qubits ? ks.num_qubits : k->num_qubits;
    if (nq == 0u) nq = 1u;

    uint32_t *q_wire = (uint32_t *)malloc((size_t)nq * sizeof(uint32_t));
    if (!q_wire) return false;

    if (nq <= QIR_MAX_VAL_BITS) {
        uint32_t bus = qir_add_input(m, (uint16_t)nq, "qbus");
        for (uint32_t q = 0u; q < nq; ++q) {
            q_wire[q] = (nq == 1u) ? bus : qir_emit_slice(m, bus, (uint16_t)q, (uint16_t)(q + 1u), NULL);
        }
    } else {
        for (uint32_t base = 0u, ci = 0u; base < nq; base += 128u, ++ci) {
            uint16_t w = (uint16_t)((nq - base > 128u) ? 128u : (nq - base));
            char rname[32];
            snprintf(rname, sizeof(rname), "qbus_%u", (unsigned)ci);
            uint32_t bus = qir_add_input(m, w, rname);
            for (uint16_t b = 0u; b < w; ++b) {
                q_wire[base + b] = (w == 1u) ? bus : qir_emit_slice(m, bus, b, (uint16_t)(b + 1u), NULL);
            }
        }
    }

    size_t nops = k->num_ops;
    qcc_kmb_op_view_t *views = nops ? (qcc_kmb_op_view_t *)calloc(nops, sizeof(qcc_kmb_op_view_t)) : NULL;
    if (nops && !views) { free(q_wire); return false; }
    qcc_kmb_cursor_t cur;
    qcc_kmb_cursor_init(&cur);
    for (size_t i = 0u; i < nops; ++i) {
        if (!qcc_kmb_cursor_next(k, &cur, &views[i])) {
            free(views);
            free(q_wire);
            return false;
        }
    }

    size_t pos = 0u;
    int32_t cond_depth = 0;
    while (pos < nops) {
        uint8_t op0 = views[pos].op;

        if (cond_depth == 0 &&
            (op0 == QCC_OP_INFO_TAG || op0 == QCC_OP_DEBUG_PRINT_EMPTY ||
             op0 == QCC_OP_DEBUG_PRINT_Q || op0 == QCC_OP_DEBUG_PRINT_C)) {
            pos++;
            continue;
        }

        /* 1A. Gidney Temporary-AND Compute: R anc; CCX c0 c1 anc */
        if (cond_depth == 0 && qir__lift_is_and_compute(views, nops, pos)) {
            uint32_t c0  = views[pos + 1u].q[0];
            uint32_t c1  = views[pos + 1u].q[1];
            uint32_t anc = views[pos + 1u].q[2];
            q_wire[anc] = qir_emit_and_vec(m, q_wire[c0], q_wire[c1], QIR_FLAG_CLEAN_ANC, NULL);
            if (out_stats) {
                out_stats->lifted_and_vec++;
                out_stats->lifted_ops_total += 2u;
            }
            pos += 2u;
            continue;
        }

        /* 1B. HMR-Dual Temporary-AND Uncompute: HMR anc b; CZ_IF q0 q1 b */
        if (cond_depth == 0 && qir__lift_is_and_uncompute(m, q_wire, views, nops, pos)) {
            uint32_t anc = views[pos].q[0];
            q_wire[anc] = qir_emit_const_u64(m, 1u, 0u, NULL);
            if (out_stats) {
                out_stats->lifted_and_vec++;
                out_stats->lifted_ops_total += 2u;
            }
            pos += 2u;
            continue;
        }

        /* 2. Linear F2 Island {X, CX, SWAP} up to 256 active qubits */
        if (cond_depth == 0 && qir__is_linf2_op(op0)) {
            uint32_t act_q[QIR_MAX_VAL_BITS];
            uint32_t k_act = 0u;
            size_t end = pos;
            while (end < nops && qir__is_linf2_op(views[end].op)) {
                if (!qir__lift_try_add_op_qubits(&views[end], act_q, &k_act)) break;
                end++;
            }
            if (k_act == 0u) { pos = end; continue; }

            qcc_bv256_t *rows = (qcc_bv256_t *)qcc_arena_xalloc_zero(arena, (size_t)k_act * sizeof(qcc_bv256_t), 8u);
            qcc_bv256_t cvec = qcc_bv_zero();
            for (uint32_t i = 0u; i < k_act; ++i) qcc_bv_set(rows[i].w, i);

            for (size_t p = pos; p < end; ++p) {
                if (views[p].op == QCC_OP_X) {
                    uint32_t lt = qir__lift_find_q(act_q, k_act, views[p].q[0]);
                    if (qcc_bv_get(cvec.w, lt)) qcc_bv_clr(cvec.w, lt);
                    else qcc_bv_set(cvec.w, lt);
                } else if (views[p].op == QCC_OP_CX) {
                    uint32_t lc = qir__lift_find_q(act_q, k_act, views[p].q[0]);
                    uint32_t lt = qir__lift_find_q(act_q, k_act, views[p].q[1]);
                    for (uint32_t w = 0u; w < QIR_WORDS; ++w) rows[lt].w[w] ^= rows[lc].w[w];
                    if (qcc_bv_get(cvec.w, lc)) {
                        if (qcc_bv_get(cvec.w, lt)) qcc_bv_clr(cvec.w, lt);
                        else qcc_bv_set(cvec.w, lt);
                    }
                } else if (views[p].op == QCC_OP_SWAP) {
                    uint32_t la = qir__lift_find_q(act_q, k_act, views[p].q[0]);
                    uint32_t lb = qir__lift_find_q(act_q, k_act, views[p].q[1]);
                    if (la != lb) {
                        qcc_bv256_t tmp = rows[la]; rows[la] = rows[lb]; rows[lb] = tmp;
                        bool ba = qcc_bv_get(cvec.w, la), bb = qcc_bv_get(cvec.w, lb);
                        if (ba != bb) {
                            if (bb) qcc_bv_set(cvec.w, la); else qcc_bv_clr(cvec.w, la);
                            if (ba) qcc_bv_set(cvec.w, lb); else qcc_bv_clr(cvec.w, lb);
                        }
                    }
                }
            }

            uint32_t bus = qir__lift_gather_bus(m, q_wire, act_q, k_act);
            uint32_t res = qir_emit_linf2(m, &bus, 1u, (uint16_t)k_act, rows, &cvec, NULL);
            for (uint32_t i = 0u; i < k_act; ++i) {
                q_wire[act_q[i]] = (k_act == 1u) ? res : qir_emit_slice(m, res, (uint16_t)i, (uint16_t)(i + 1u), NULL);
            }
            if (out_stats) {
                out_stats->lifted_linf2++;
                out_stats->lifted_ops_total += (uint32_t)(end - pos);
            }
            pos = end;
            continue;
        }

        /* 3. Diagonal Phase-Polynomial Island {Z, CZ, CCZ, Z_POW} up to 256 active qubits */
        if (cond_depth == 0 && qir__is_phase_op(op0)) {
            uint32_t act_q[QIR_MAX_VAL_BITS];
            uint32_t k_act = 0u;
            size_t end = pos;
            while (end < nops && qir__is_phase_op(views[end].op)) {
                if (!qir__lift_try_add_op_qubits(&views[end], act_q, &k_act)) break;
                end++;
            }
            size_t nterms = end - pos;
            qir_phase_term_t *terms = nterms
                ? (qir_phase_term_t *)qcc_arena_xalloc_zero(arena, nterms * sizeof(qir_phase_term_t), 16u)
                : NULL;
            for (size_t p = pos; p < end; ++p) {
                qir_phase_term_t *pt = &terms[p - pos];
                const qcc_kmb_op_view_t *v = &views[p];
                pt->is_xor = 0u;
                for (uint8_t qi = 0u; qi < v->n_q; ++qi) {
                    qcc_bv_set(pt->mask, qir__lift_find_q(act_q, k_act, v->q[qi]));
                }
                pt->angle = (v->op == QCC_OP_Z_POW) ? qcc_u128_from_angle(v->angle) : ((flow_u128)1u << 127u);
            }
            uint32_t bus = qir__lift_gather_bus(m, q_wire, act_q, k_act);
            (void)qir_emit_phase_poly(m, &bus, 1u, terms, (uint32_t)nterms, NULL);
            if (out_stats) {
                out_stats->lifted_phase++;
                out_stats->lifted_ops_total += (uint32_t)nterms;
            }
            pos = end;
            continue;
        }

        /* 4. Opaque Balanced KMB Island bounded to <= QIR_MAX_VAL_BITS touched qubits */
        uint32_t act_q[QIR_MAX_VAL_BITS];
        uint32_t k_act = 0u;
        size_t start = pos;
        while (pos < nops) {
            if (cond_depth == 0 && pos > start) {
                if (views[pos].op == QCC_OP_PUSH_CONDITION) {
                    /* Allow balanced region unless already large */
                }
            }
            if (!qir__lift_try_add_op_qubits(&views[pos], act_q, &k_act)) {
                if (pos == start) {
                    free(views);
                    free(q_wire);
                    return false;
                }
                if (cond_depth == 0) break;
                free(views);
                free(q_wire);
                return false;
            }
            uint8_t op = views[pos].op;
            if (op == QCC_OP_PUSH_CONDITION) cond_depth++;
            else if (op == QCC_OP_POP_CONDITION) {
                if (cond_depth > 0) cond_depth--;
            }
            pos++;
            if (cond_depth == 0 && pos < nops) {
                if (qir__is_linf2_op(views[pos].op) || qir__is_phase_op(views[pos].op) ||
                    qir__lift_is_and_compute(views, nops, pos) ||
                    qir__lift_is_and_uncompute(m, q_wire, views, nops, pos)) {
                    break;
                }
                /* Also check if next op won't fit our 256-qubit island cap */
                uint32_t probe_q[QIR_MAX_VAL_BITS];
                memcpy(probe_q, act_q, (size_t)k_act * sizeof(uint32_t));
                uint32_t probe_k = k_act;
                if (!qir__lift_try_add_op_qubits(&views[pos], probe_q, &probe_k)) break;
            }
        }

        qcc_kmb_circuit_t sub;
        qcc_kmb_init(&sub);
        sub.num_qubits = k_act;
        uint32_t max_src_b = ks.num_bits > k->num_bits ? ks.num_bits : k->num_bits;
        uint32_t *b_map = max_src_b ? (uint32_t *)qcc_xmalloc((size_t)max_src_b * sizeof(uint32_t)) : NULL;
        if (b_map) {
            for (uint32_t bi = 0u; bi < max_src_b; ++bi) b_map[bi] = UINT32_MAX;
        }
        uint32_t nb_act = 0u;
        for (size_t p = start; p < pos; ++p) {
            qcc_kmb_op_view_t lv = views[p];
            for (uint8_t qi = 0u; qi < lv.n_q; ++qi) {
                lv.q[qi] = qir__lift_find_q(act_q, k_act, lv.q[qi]);
            }
            if (b_map && lv.has_b && lv.b < max_src_b) {
                if (b_map[lv.b] == UINT32_MAX) b_map[lv.b] = nb_act++;
                lv.b = b_map[lv.b];
            }
            if (b_map && lv.has_cond && lv.cond < max_src_b) {
                if (b_map[lv.cond] == UINT32_MAX) b_map[lv.cond] = nb_act++;
                lv.cond = b_map[lv.cond];
            }
            qcc_kmb_emit_view(&sub, &lv);
        }
        free(b_map);
        sub.num_bits = nb_act;

        if (k_act > 0u) {
            uint32_t id_map[QIR_MAX_VAL_BITS];
            for (uint32_t i = 0u; i < k_act; ++i) id_map[i] = i;
            uint32_t bus = qir__lift_gather_bus(m, q_wire, act_q, k_act);
            uint32_t res = qir_emit_opaque_kmb(m, &bus, 1u, (uint16_t)k_act, &sub, id_map, (uint16_t)k_act, id_map, NULL);
            for (uint32_t i = 0u; i < k_act; ++i) {
                q_wire[act_q[i]] = (k_act == 1u) ? res : qir_emit_slice(m, res, (uint16_t)i, (uint16_t)(i + 1u), NULL);
            }
        } else {
            (void)qir_emit_opaque_kmb(m, NULL, 0u, 0u, &sub, NULL, 0u, NULL, NULL);
        }
        qcc_kmb_free(&sub);

        if (out_stats) {
            out_stats->lifted_opaque++;
            out_stats->lifted_ops_total += (uint32_t)(pos - start);
        }
    }

    free(views);
    if (nq <= QIR_MAX_VAL_BITS) {
        qir_add_output(m, (nq == 1u) ? q_wire[0] : qir_emit_concat(m, q_wire, (uint16_t)nq, NULL));
    } else {
        for (uint32_t base = 0u; base < nq; base += 128u) {
            uint16_t w = (uint16_t)((nq - base > 128u) ? 128u : (nq - base));
            uint32_t out_v = (w == 1u) ? q_wire[base] : qir_emit_concat(m, q_wire + base, w, NULL);
            qir_add_output(m, out_v);
        }
    }
    free(q_wire);
    qir_recompute_zero_masks(m);
    return true;
}

#ifdef __cplusplus
}
#endif

#endif /* QIR_LIFT_H */
