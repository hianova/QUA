/*
 * qir_verify.h — Pure C17 64-Lane Bit-Sliced .kmb Simulator & Dual-Seed Formal Verifier
 * Part of QCC (Quantum Circuit Compiler) for FLOW v7.0.1.
 * Pure C17, zero external dependencies.
 */
#ifndef QIR_VERIFY_H
#define QIR_VERIFY_H

#include "qir_opt.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    bool      passed;
    uint64_t  states_tested;
    uint64_t  value_mismatches;
    uint64_t  seed_drift_mismatches;
    uint64_t  phase_mismatches;
    uint64_t  dirty_ancilla_mismatches;
    uint64_t  max_abs_err;
    uint64_t  first_fail_state;
} qir_verify_report_t;

static inline bool qir_sim_kmb_64_ex(
    const qcc_kmb_circuit_t *k,
    const uint64_t *init_qstate, /* [k->num_qubits], NULL => all 0 */
    const uint64_t *init_bstate, /* [k->num_bits], NULL => all 0 */
    uint64_t seed,
    uint64_t *out_qstate,        /* [k->num_qubits] */
    uint64_t *out_bstate,        /* Optional [k->num_bits] */
    flow_u128 out_phase[64]
) {
    if (!k) return false;
    qcc_kmb_stats_t ks;
    if (!qcc_kmb_recompute_stats(k, &ks)) return false;

    uint32_t nq = ks.num_qubits > k->num_qubits ? ks.num_qubits : k->num_qubits;
    uint32_t nb = ks.num_bits > k->num_bits ? ks.num_bits : k->num_bits;

    uint64_t *qstate = nq ? (uint64_t *)calloc(nq, sizeof(uint64_t)) : NULL;
    uint64_t *bstate = nb ? (uint64_t *)calloc(nb, sizeof(uint64_t)) : NULL;
    size_t cond_cap = 64u, cond_len = 0u;
    uint64_t *cond_stack = (uint64_t *)calloc(cond_cap, sizeof(uint64_t));
    if ((nq && !qstate) || (nb && !bstate) || !cond_stack) {
        free(qstate); free(bstate); free(cond_stack);
        return false;
    }

    if (init_qstate && k->num_qubits > 0u) memcpy(qstate, init_qstate, k->num_qubits * sizeof(uint64_t));
    if (init_bstate && k->num_bits > 0u)   memcpy(bstate, init_bstate, k->num_bits * sizeof(uint64_t));
    flow_u128 phase[64];
    memset(phase, 0, sizeof(phase));

    uint64_t rng = seed ? seed : UINT64_C(0x987654321FEDCBA0);
    uint64_t cond_all = UINT64_MAX;
    const flow_u128 pi_turn = (flow_u128)1u << 127u;

    qcc_kmb_cursor_t cur;
    qcc_kmb_cursor_init(&cur);
    qcc_kmb_op_view_t v;
    while (qcc_kmb_cursor_next(k, &cur, &v)) {
        if (v.op == QCC_OP_PUSH_CONDITION) {
            if (cond_len >= cond_cap) {
                cond_cap *= 2u;
                cond_stack = (uint64_t *)realloc(cond_stack, cond_cap * sizeof(uint64_t));
            }
            cond_stack[cond_len++] = bstate[v.cond];
            cond_all &= bstate[v.cond];
            continue;
        }
        if (v.op == QCC_OP_POP_CONDITION) {
            if (cond_len > 0u) cond_len--;
            cond_all = UINT64_MAX;
            for (size_t s = 0u; s < cond_len; ++s) cond_all &= cond_stack[s];
            continue;
        }

        uint64_t active = cond_all & (v.has_cond ? bstate[v.cond] : UINT64_MAX);
        uint8_t base = v.op & (uint8_t)(~QCC_OP_FLAG_COND);
        switch (base) {
            case QCC_OP_NEG:
                for (uint32_t l = 0u; l < 64u; ++l) if ((active >> l) & 1u) phase[l] += pi_turn;
                break;
            case QCC_OP_BIT_STORE0:
                bstate[v.b] &= ~active;
                break;
            case QCC_OP_BIT_STORE1:
                bstate[v.b] |= active;
                break;
            case QCC_OP_BIT_INVERT:
                bstate[v.b] ^= active;
                break;
            case QCC_OP_X:
                qstate[v.q[0]] ^= active;
                break;
            case QCC_OP_Z: {
                uint64_t hit = active & qstate[v.q[0]];
                for (uint32_t l = 0u; l < 64u; ++l) if ((hit >> l) & 1u) phase[l] += pi_turn;
                break;
            }
            case QCC_OP_R: {
                uint64_t rbits = qcc_splitmix64(&rng);
                uint64_t hit = active & qstate[v.q[0]] & rbits;
                for (uint32_t l = 0u; l < 64u; ++l) if ((hit >> l) & 1u) phase[l] += pi_turn;
                qstate[v.q[0]] &= ~active;
                break;
            }
            case QCC_OP_Z_POW: {
                flow_u128 ang = qcc_u128_from_angle(v.angle);
                uint64_t hit = active & qstate[v.q[0]];
                for (uint32_t l = 0u; l < 64u; ++l) if ((hit >> l) & 1u) phase[l] += ang;
                break;
            }
            case QCC_OP_HMR: {
                uint64_t rbits = qcc_splitmix64(&rng);
                uint64_t flip = active & rbits & qstate[v.q[0]];
                for (uint32_t l = 0u; l < 64u; ++l) if ((flip >> l) & 1u) phase[l] += pi_turn;
                bstate[v.b] = (bstate[v.b] & ~active) | (rbits & active);
                qstate[v.q[0]] &= ~active;
                break;
            }
            case QCC_OP_CX:
                qstate[v.q[1]] ^= (active & qstate[v.q[0]]);
                break;
            case QCC_OP_CZ: {
                uint64_t hit = active & qstate[v.q[0]] & qstate[v.q[1]];
                for (uint32_t l = 0u; l < 64u; ++l) if ((hit >> l) & 1u) phase[l] += pi_turn;
                break;
            }
            case QCC_OP_SWAP: {
                uint64_t diff = (qstate[v.q[0]] ^ qstate[v.q[1]]) & active;
                qstate[v.q[0]] ^= diff;
                qstate[v.q[1]] ^= diff;
                break;
            }
            case QCC_OP_CCX:
                qstate[v.q[2]] ^= (active & qstate[v.q[0]] & qstate[v.q[1]]);
                break;
            case QCC_OP_CCZ: {
                uint64_t hit = active & qstate[v.q[0]] & qstate[v.q[1]] & qstate[v.q[2]];
                for (uint32_t l = 0u; l < 64u; ++l) if ((hit >> l) & 1u) phase[l] += pi_turn;
                break;
            }
            case QCC_OP_CSWAP: {
                uint64_t act2 = active & qstate[v.q[0]];
                uint64_t diff = (qstate[v.q[1]] ^ qstate[v.q[2]]) & act2;
                qstate[v.q[1]] ^= diff;
                qstate[v.q[2]] ^= diff;
                break;
            }
            default:
                break;
        }
    }

    if (out_qstate && k->num_qubits > 0u) memcpy(out_qstate, qstate, k->num_qubits * sizeof(uint64_t));
    if (out_bstate && k->num_bits > 0u) memcpy(out_bstate, bstate, k->num_bits * sizeof(uint64_t));
    if (out_phase) memcpy(out_phase, phase, sizeof(phase));

    free(qstate);
    free(bstate);
    free(cond_stack);
    return true;
}

static inline bool qir_sim_kmb_64(
    const qcc_kmb_circuit_t *k,
    const uint64_t *init_qstate, /* [k->num_qubits], NULL => all 0 */
    uint64_t seed,
    uint64_t *out_qstate,        /* [k->num_qubits] */
    uint64_t *out_bstate,        /* Optional [k->num_bits] */
    flow_u128 out_phase[64]
) {
    return qir_sim_kmb_64_ex(k, init_qstate, NULL, seed, out_qstate, out_bstate, out_phase);
}

static inline bool qir_verify_kmb_vs_qir(
    const qcc_kmb_circuit_t *kmb,
    const qir_module_t *mod,
    const qir_phys_plan_t *plan,
    uint32_t num_batches,
    uint64_t seed,
    bool allow_approx,
    uint64_t max_approx_lsb_err,
    qir_verify_report_t *out_rep
) {
    qir_verify_report_t rep;
    memset(&rep, 0, sizeof(rep));
    if (!kmb || !mod || !plan) return false;

    qcc_arena_t scratch;
    if (!qcc_arena_init(&scratch, 512u * 1024u)) return false;

    uint32_t k_nq = kmb->num_qubits ? kmb->num_qubits : 1u;
    uint32_t k_nb = kmb->num_bits ? kmb->num_bits : 1u;
    uint64_t *init_q = (uint64_t *)calloc(k_nq, sizeof(uint64_t));
    uint64_t *init_b = (uint64_t *)calloc(k_nb, sizeof(uint64_t));
    uint64_t *out_q1 = (uint64_t *)calloc(k_nq, sizeof(uint64_t));
    uint64_t *out_q2 = (uint64_t *)calloc(k_nq, sizeof(uint64_t));
    uint64_t *out_b1 = (uint64_t *)calloc(k_nb, sizeof(uint64_t));
    uint64_t *out_b2 = (uint64_t *)calloc(k_nb, sizeof(uint64_t));
    uint8_t  *is_io_q = (uint8_t *)calloc(k_nq, 1u);
    qir_val_vec64_t *in_vecs = mod->num_inputs ? (qir_val_vec64_t *)calloc(mod->num_inputs, sizeof(qir_val_vec64_t)) : NULL;

    for (uint32_t ii = 0u; ii < mod->num_inputs; ++ii) {
        uint32_t vid = mod->inputs[ii];
        if (plan->val_allocs[vid].qubits) {
            for (uint16_t bit = 0u; bit < mod->vals[vid].width; ++bit) {
                uint32_t qid = plan->val_allocs[vid].qubits[bit];
                if (qid < k_nq) is_io_q[qid] = 1u;
            }
        }
    }
    for (uint32_t oi = 0u; oi < mod->num_outputs; ++oi) {
        uint32_t vid = mod->outputs[oi];
        if (plan->val_allocs[vid].qubits) {
            for (uint16_t bit = 0u; bit < mod->vals[vid].width; ++bit) {
                uint32_t qid = plan->val_allocs[vid].qubits[bit];
                if (qid < k_nq) is_io_q[qid] = 1u;
            }
        }
    }

    uint64_t rng = seed ? seed : UINT64_C(0xDEADBEEF12345678);
    flow_u128 global_phase_offset = 0u;
    bool has_offset = false;

    for (uint32_t b = 0u; b < num_batches; ++b) {
        memset(init_q, 0, k_nq * sizeof(uint64_t));
        memset(init_b, 0, k_nb * sizeof(uint64_t));
        memset(out_q1, 0, k_nq * sizeof(uint64_t));
        memset(out_q2, 0, k_nq * sizeof(uint64_t));
        memset(out_b1, 0, k_nb * sizeof(uint64_t));
        memset(out_b2, 0, k_nb * sizeof(uint64_t));
        for (uint32_t ii = 0u; ii < mod->num_inputs; ++ii) {
            uint32_t vid = mod->inputs[ii];
            uint16_t w = mod->vals[vid].width;
            qcc_bv256_t mask = qcc_bv_mask_n(w);
            for (uint32_t l = 0u; l < 64u; ++l) {
                for (uint32_t ww = 0u; ww < QIR_WORDS; ++ww) {
                    uint64_t rv = (b == 0u && l == 0u) ? 0u : qcc_splitmix64(&rng);
                    if (b == 0u && l == 1u) rv = UINT64_MAX;
                    in_vecs[ii].lanes[l][ww] = rv & mask.w[ww];
                }
            }
            if (mod->vals[vid].rep == QIR_REP_CLASSICAL && plan->val_allocs[vid].cbits) {
                for (uint16_t bit = 0u; bit < w; ++bit) {
                    uint32_t bid = plan->val_allocs[vid].cbits[bit];
                    uint64_t lane_word = 0u;
                    for (uint32_t l = 0u; l < 64u; ++l) {
                        if (qcc_bv_get(in_vecs[ii].lanes[l], bit)) lane_word |= (UINT64_C(1) << l);
                    }
                    if (bid < k_nb) init_b[bid] = lane_word;
                }
            } else if (plan->val_allocs[vid].qubits) {
                for (uint16_t bit = 0u; bit < w; ++bit) {
                    uint32_t qid = plan->val_allocs[vid].qubits[bit];
                    uint64_t lane_word = 0u;
                    for (uint32_t l = 0u; l < 64u; ++l) {
                        if (qcc_bv_get(in_vecs[ii].lanes[l], bit)) lane_word |= (UINT64_C(1) << l);
                    }
                    if (qid < k_nq) init_q[qid] = lane_word;
                }
            }
        }

        qcc_arena_reset(&scratch);
        qir_val_vec64_t *ref_vals = NULL;
        flow_u128 ref_phase[64], p1[64], p2[64];
        if (!qir_eval_64(mod, in_vecs, seed ^ (uint64_t)b, &scratch, &ref_vals, ref_phase) ||
            !qir_sim_kmb_64_ex(kmb, init_q, init_b, (seed + 111u) ^ (uint64_t)b, out_q1, out_b1, p1) ||
            !qir_sim_kmb_64_ex(kmb, init_q, init_b, (seed + 99991u) ^ (uint64_t)b, out_q2, out_b2, p2)) {
            rep.value_mismatches++;
            break;
        }

        for (uint32_t l = 0u; l < 64u; ++l) {
            rep.states_tested++;
            if (p1[l] != p2[l]) {
                rep.seed_drift_mismatches++;
                rep.phase_mismatches++;
            }
            flow_u128 diff = p1[l] - ref_phase[l];
            if (!has_offset) {
                global_phase_offset = diff;
                has_offset = true;
            } else if (diff != global_phase_offset) {
                rep.phase_mismatches++;
            }

            for (uint32_t oi = 0u; oi < mod->num_outputs; ++oi) {
                uint32_t vid = mod->outputs[oi];
                uint16_t w = mod->vals[vid].width;
                if (w == 0u || !plan->val_allocs[vid].qubits) continue;
                uint64_t got_w0 = 0u, got2_w0 = 0u, ref_w0 = ref_vals[vid].lanes[l][0];
                bool bit_mismatch = false;
                for (uint16_t bit = 0u; bit < w; ++bit) {
                    uint32_t qid = plan->val_allocs[vid].qubits[bit];
                    uint64_t b1 = (qid < k_nq) ? ((out_q1[qid] >> l) & 1u) : 0u;
                    uint64_t b2 = (qid < k_nq) ? ((out_q2[qid] >> l) & 1u) : 0u;
                    uint64_t br = qcc_bv_get(ref_vals[vid].lanes[l], bit) ? 1u : 0u;
                    if (b1 != b2) rep.seed_drift_mismatches++;
                    if (bit < 64u) {
                        got_w0  |= (b1 << bit);
                        got2_w0 |= (b2 << bit);
                    }
                    if (b1 != br) bit_mismatch = true;
                }
                if (bit_mismatch) {
                    if (allow_approx && w <= 64u) {
                        uint64_t err = (got_w0 > ref_w0) ? (got_w0 - ref_w0) : (ref_w0 - got_w0);
                        if (err > rep.max_abs_err) rep.max_abs_err = err;
                        if (err > max_approx_lsb_err) {
                            rep.value_mismatches++;
                            if (rep.value_mismatches == 1u) rep.first_fail_state = rep.states_tested - 1u;
                        }
                    } else {
                        rep.value_mismatches++;
                        if (rep.value_mismatches == 1u) rep.first_fail_state = rep.states_tested - 1u;
                    }
                }
                (void)got2_w0;
            }

            for (uint32_t ci = 0u; ci < mod->num_classical_outputs; ++ci) {
                uint32_t vid = mod->classical_outputs[ci];
                if (vid >= mod->num_vals) continue;
                uint16_t w = mod->vals[vid].width;
                const uint32_t *cb = plan->val_allocs[vid].cbits;
                for (uint16_t bit = 0u; bit < w; ++bit) {
                    uint32_t bid = cb ? cb[bit] : UINT32_MAX;
                    uint64_t b1 = (bid < k_nb) ? ((out_b1[bid] >> l) & 1u) : 0u;
                    uint64_t b2 = (bid < k_nb) ? ((out_b2[bid] >> l) & 1u) : 0u;
                    uint64_t br = qcc_bv_get(ref_vals[vid].lanes[l], bit) ? 1u : 0u;
                    if (b1 != b2) rep.seed_drift_mismatches++;
                    if (b1 != br) {
                        rep.value_mismatches++;
                        if (rep.value_mismatches == 1u) rep.first_fail_state = rep.states_tested - 1u;
                    }
                }
            }

            if (plan->clean_scratch_certified) {
                for (uint32_t q = 0u; q < k_nq; ++q) {
                    if (!is_io_q[q]) {
                        if (((out_q1[q] >> l) & 1u) != 0u || ((out_q2[q] >> l) & 1u) != 0u) {
                            rep.dirty_ancilla_mismatches++;
                        }
                    }
                }
            }
        }
    }

    rep.passed = (rep.value_mismatches == 0u &&
                  rep.seed_drift_mismatches == 0u &&
                  rep.phase_mismatches == 0u &&
                  rep.dirty_ancilla_mismatches == 0u);
    if (out_rep) *out_rep = rep;

    free(init_q);
    free(init_b);
    free(out_q1);
    free(out_q2);
    free(out_b1);
    free(out_b2);
    free(is_io_q);
    free(in_vecs);
    qcc_arena_free(&scratch);
    return rep.passed;
}

static inline bool qir_verify_kmb_vs_kmb(
    const qcc_kmb_circuit_t *cand,
    const qcc_kmb_circuit_t *ref,
    uint32_t num_io_qubits,
    uint32_t num_batches,
    uint64_t seed,
    qir_verify_report_t *out_rep
) {
    qir_verify_report_t rep;
    memset(&rep, 0, sizeof(rep));
    if (!cand || !ref) return false;

    uint32_t nq = cand->num_qubits > ref->num_qubits ? cand->num_qubits : ref->num_qubits;
    if (num_io_qubits > nq) nq = num_io_qubits;
    if (nq == 0u) nq = 1u;
    uint32_t io = num_io_qubits ? num_io_qubits : nq;

    uint64_t *init_q = (uint64_t *)calloc(nq, sizeof(uint64_t));
    uint64_t *qc1 = (uint64_t *)calloc(nq, sizeof(uint64_t));
    uint64_t *qc2 = (uint64_t *)calloc(nq, sizeof(uint64_t));
    uint64_t *qr1 = (uint64_t *)calloc(nq, sizeof(uint64_t));

    uint64_t rng = seed ? seed : UINT64_C(0x13579BDF2468ACE0);
    flow_u128 global_off = 0u;
    bool has_off = false;

    for (uint32_t b = 0u; b < num_batches; ++b) {
        memset(init_q, 0, nq * sizeof(uint64_t));
        for (uint32_t q = 0u; q < io; ++q) {
            init_q[q] = (b == 0u) ? ((q & 1u) ? UINT64_C(0xAAAAAAAAAAAAAAAA) : UINT64_C(0x5555555555555555)) : qcc_splitmix64(&rng);
        }
        flow_u128 pc1[64], pc2[64], pr1[64];
        memset(qc1, 0, nq * sizeof(uint64_t));
        memset(qc2, 0, nq * sizeof(uint64_t));
        memset(qr1, 0, nq * sizeof(uint64_t));
        if (!qir_sim_kmb_64(cand, init_q, (seed + 17u) ^ b, qc1, NULL, pc1) ||
            !qir_sim_kmb_64(cand, init_q, (seed + 7919u) ^ b, qc2, NULL, pc2) ||
            !qir_sim_kmb_64(ref,  init_q, (seed + 3571u) ^ b, qr1, NULL, pr1)) {
            rep.value_mismatches++;
            break;
        }
        for (uint32_t l = 0u; l < 64u; ++l) {
            rep.states_tested++;
            if (pc1[l] != pc2[l]) { rep.seed_drift_mismatches++; rep.phase_mismatches++; }
            flow_u128 d = pc1[l] - pr1[l];
            if (!has_off) { global_off = d; has_off = true; }
            else if (d != global_off) rep.phase_mismatches++;

            for (uint32_t q = 0u; q < io; ++q) {
                uint64_t b1 = (qc1[q] >> l) & 1u;
                uint64_t b2 = (qc2[q] >> l) & 1u;
                uint64_t br = (qr1[q] >> l) & 1u;
                if (b1 != b2) rep.seed_drift_mismatches++;
                if (b1 != br) {
                    rep.value_mismatches++;
                    if (rep.value_mismatches == 1u) rep.first_fail_state = rep.states_tested - 1u;
                }
            }
        }
    }

    rep.passed = (rep.value_mismatches == 0u && rep.seed_drift_mismatches == 0u && rep.phase_mismatches == 0u);
    if (out_rep) *out_rep = rep;
    free(init_q); free(qc1); free(qc2); free(qr1);
    return rep.passed;
}

#ifdef __cplusplus
}
#endif

#endif /* QIR_VERIFY_H */
