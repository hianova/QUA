/*
 * pqc_main.c — PQC / NIST Security Baseline Quantum Cryptanalysis Compiler & Resource Estimator
 * Uses ONLY declarative QIR builders (pqc_oracles.h) compiled via QCC (qir_pass_kmb.h).
 * Pure C17, zero <math.h>, zero raw qcc_kmb_* gate emission calls.
 */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#include "qir_pass_kmb.h"
#include "qir_verify.h"
#include "pqc_oracles.h"
#include "pqc_nist_est.h"

#include <sys/stat.h>
#include <time.h>

static double pqc__now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

typedef enum {
    PQC_KAT_NONE = 0,
    PQC_KAT_AES,
    PQC_KAT_SHA256,
    PQC_KAT_ASCON,
    PQC_KAT_KECCAK400,
    PQC_KAT_MCELIECE_GOPPA,
    PQC_KAT_MCELIECE_ISD,
    PQC_KAT_MLKEM_LWE,
    PQC_KAT_MQ24
} pqc_kat_kind_t;

typedef struct {
    pqc_kat_kind_t kind;
    uint32_t       rounds;
    uint32_t       key_bits;
    uint8_t        pt[16];
    uint8_t        target_bytes[32];
    uint64_t       ascon_iv_n[3];
    uint64_t       ascon_tag[2];
    uint16_t       keccak_rate[8];
    uint16_t       goppa_coeffs[9];
    uint32_t       goppa_poly12;
    uint16_t       goppa_target_inv;
    uint64_t       isd_h[24];
    uint32_t       isd_syn;
    uint8_t        isd_wt;
    uint16_t       lwe_a[4][8];
    uint16_t       lwe_b[4];
    uint16_t       lwe_eta;
    pqc_mq24_sys_t mq_sys;
    uint32_t       mq_target;
    uint8_t        planted_sol[32];
    uint32_t       planted_sol_bits;
} pqc_instance_t;

static void pqc__unpack_lanes_bytes(
    const uint64_t *lanes_q,
    const uint32_t *q_map,
    uint32_t n_bits,
    uint32_t lane,
    uint8_t *out_bytes
) {
    uint32_t nb = (n_bits + 7u) / 8u;
    memset(out_bytes, 0, nb);
    for (uint32_t b = 0u; b < n_bits; ++b) {
        if ((lanes_q[q_map[b]] >> lane) & 1u) pqc__set_bit(out_bytes, b, 1u);
    }
}

static bool pqc__eval_planted_hit(const pqc_instance_t *inst, const uint8_t *in_bytes) {
    switch (inst->kind) {
        case PQC_KAT_AES: {
            uint8_t ct[16];
            pqc_ref_aes_encrypt(inst->pt, in_bytes, inst->key_bits, inst->rounds, ct);
            return memcmp(ct, inst->target_bytes, 16u) == 0;
        }
        case PQC_KAT_SHA256: {
            uint8_t dig[32];
            pqc_ref_sha256_msg256(in_bytes, inst->rounds, dig);
            return memcmp(dig, inst->target_bytes, 32u) == 0;
        }
        case PQC_KAT_ASCON: {
            uint64_t s[5];
            s[0] = inst->ascon_iv_n[0];
            s[1] = pqc__pack_u64_le(in_bytes);
            s[2] = pqc__pack_u64_le(in_bytes + 8);
            s[3] = inst->ascon_iv_n[1];
            s[4] = inst->ascon_iv_n[2];
            pqc_ref_ascon_p(s, inst->rounds);
            return s[0] == inst->ascon_tag[0] && s[1] == inst->ascon_tag[1];
        }
        case PQC_KAT_KECCAK400: {
            uint16_t st[25];
            memset(st, 0, sizeof(st));
            for (uint32_t i = 0u; i < 8u; ++i) {
                st[i] = (uint16_t)in_bytes[2u * i] | ((uint16_t)in_bytes[2u * i + 1u] << 8u);
            }
            pqc_ref_keccak_f400(st, inst->rounds);
            return memcmp(st, inst->keccak_rate, sizeof(inst->keccak_rate)) == 0;
        }
        case PQC_KAT_MCELIECE_GOPPA: {
            uint16_t alpha = (uint16_t)(in_bytes[0] | ((uint16_t)(in_bytes[1] & 0x0Fu) << 8u));
            return pqc_ref_mceliece_goppa_inv(alpha, inst->goppa_coeffs, inst->goppa_poly12) == inst->goppa_target_inv;
        }
        case PQC_KAT_MCELIECE_ISD: {
            uint64_t evec = pqc__pack_u64_le(in_bytes);
            uint32_t syn = 0u;
            uint8_t  wt  = 0u;
            pqc_ref_mceliece_isd(evec, inst->isd_h, &syn, &wt);
            return syn == inst->isd_syn && wt == inst->isd_wt;
        }
        case PQC_KAT_MLKEM_LWE: {
            uint32_t sp = (uint32_t)in_bytes[0] | ((uint32_t)in_bytes[1] << 8u) | ((uint32_t)in_bytes[2] << 16u);
            uint16_t res[4];
            uint32_t worst = pqc_ref_mlkem_lwe(sp, inst->lwe_a, inst->lwe_b, res);
            return worst <= inst->lwe_eta;
        }
        case PQC_KAT_MQ24: {
            uint32_t x24 = (uint32_t)in_bytes[0] | ((uint32_t)in_bytes[1] << 8u) | ((uint32_t)in_bytes[2] << 16u);
            return pqc_ref_mq24(x24, &inst->mq_sys) == inst->mq_target;
        }
        default:
            return false;
    }
}

static bool pqc_compile_and_verify(
    const char *label,
    const char *family,
    uint32_t search_bits,
    qir_module_t *mod,
    qcc_arena_t *arena,
    const qcc_opt_cfg_t *cfg_in,
    const pqc_instance_t *inst,
    bool phase_oracle,
    const char *kmb_out_path_opt,
    uint32_t verify_batches,
    pqc_nist_report_t *out_rep,
    uint64_t *out_verified_states,
    double *out_compile_sec
) {
    double t0 = pqc__now_sec();
    qir_recompute_zero_masks(mod);
    qcc_opt_cfg_t cfg = cfg_in ? *cfg_in : qcc_opt_cfg_default();
    cfg.require_clean_scratch = true;
    if (mod->num_ops > 150u) {
        if (cfg.searu_iters > 16u) cfg.searu_iters = 16u;
        if (cfg.flow_iters > 4u)   cfg.flow_iters  = 4u;
    }

    qir_egraph_t eg;
    qir_phys_plan_t plan;
    if (!qir_optimize_module(mod, arena, &cfg, &eg, &plan)) {
        fprintf(stderr, "[FAIL] %s: qir_optimize_module failed\n", label);
        return false;
    }
    if (!plan.disjointness_certified || !plan.clean_scratch_certified) {
        fprintf(stderr, "[FAIL] %s: disjoint=%d clean=%d\n",
                label, (int)plan.disjointness_certified, (int)plan.clean_scratch_certified);
        qir_phys_plan_free(&plan);
        qir_egraph_free(&eg);
        return false;
    }

    qcc_kmb_circuit_t kmb;
    if (!qir_codegen_emit(&plan, &kmb)) {
        fprintf(stderr, "[FAIL] %s: qir_codegen_emit failed\n", label);
        qir_phys_plan_free(&plan);
        qir_egraph_free(&eg);
        return false;
    }
    qir_kmb_pass_stats_t pstats;
    if (!qir_pass_kmb_optimize_with_plan(&kmb, &plan, &pstats)) {
        fprintf(stderr, "[FAIL] %s: qir_pass_kmb failed\n", label);
        qcc_kmb_free(&kmb);
        qir_phys_plan_free(&plan);
        qir_egraph_free(&eg);
        return false;
    }
    qcc_kmb_stats_t ks;
    if (!qcc_kmb_recompute_stats(&kmb, &ks)) {
        qcc_kmb_free(&kmb);
        qir_phys_plan_free(&plan);
        qir_egraph_free(&eg);
        return false;
    }

    qir_verify_report_t vrep;
    bool qir_ok = qir_verify_kmb_vs_qir(&kmb, mod, &plan, verify_batches ? verify_batches : 2u,
                                        UINT64_C(0x50C12345), false, 0u, &vrep);
    if (!qir_ok) {
        fprintf(stderr, "[FAIL] %s: QIR dual-seed verify failed (val=%llu seed=%llu phase=%llu dirty=%llu)\n",
                label,
                (unsigned long long)vrep.value_mismatches,
                (unsigned long long)vrep.seed_drift_mismatches,
                (unsigned long long)vrep.phase_mismatches,
                (unsigned long long)vrep.dirty_ancilla_mismatches);
        qcc_kmb_free(&kmb);
        qir_phys_plan_free(&plan);
        qir_egraph_free(&eg);
        return false;
    }

    /* Independent planted-solution / classical-KAT verification over 64-lane dual-seed KMB simulation */
    uint64_t extra_states = 0u;
    if (inst && mod->num_inputs == 1u) {
        uint32_t in_vid = mod->inputs[0];
        uint16_t in_w   = mod->vals[in_vid].width;
        uint32_t k_nq   = ks.num_qubits ? ks.num_qubits : 1u;
        uint32_t k_nb   = ks.num_bits ? ks.num_bits : 1u;
        uint64_t *init_q = (uint64_t *)calloc(k_nq, sizeof(uint64_t));
        uint64_t *out_q1 = (uint64_t *)calloc(k_nq, sizeof(uint64_t));
        uint64_t *out_q2 = (uint64_t *)calloc(k_nq, sizeof(uint64_t));
        uint64_t *init_b = (uint64_t *)calloc(k_nb, sizeof(uint64_t));
        uint64_t *out_b1 = (uint64_t *)calloc(k_nb, sizeof(uint64_t));
        uint64_t *out_b2 = (uint64_t *)calloc(k_nb, sizeof(uint64_t));
        uint8_t  *io_q   = (uint8_t *)calloc(k_nq, 1u);
        if (!init_q || !out_q1 || !out_q2 || !init_b || !out_b1 || !out_b2 || !io_q) {
            free(init_q); free(out_q1); free(out_q2);
            free(init_b); free(out_b1); free(out_b2); free(io_q);
            qcc_kmb_free(&kmb);
            qir_phys_plan_free(&plan);
            qir_egraph_free(&eg);
            return false;
        }

        for (uint16_t b = 0u; b < in_w; ++b) io_q[plan.val_allocs[in_vid].qubits[b]] = 1u;
        for (uint32_t oi = 0u; oi < mod->num_outputs; ++oi) {
            uint32_t ov = mod->outputs[oi];
            for (uint16_t b = 0u; b < mod->vals[ov].width; ++b) io_q[plan.val_allocs[ov].qubits[b]] = 1u;
        }

        uint64_t rng = UINT64_C(0x9E3779B97F4A7C15);
        for (uint32_t l = 0u; l < 64u; ++l) {
            for (uint16_t b = 0u; b < in_w; ++b) {
                uint8_t bit = 0u;
                if (l == 0u && inst->planted_sol_bits == in_w) {
                    bit = pqc__get_bit(inst->planted_sol, b);
                } else {
                    bit = (uint8_t)((qcc_splitmix64(&rng) >> (b & 63u)) & 1u);
                }
                if (bit) init_q[plan.val_allocs[in_vid].qubits[b]] |= (UINT64_C(1) << l);
            }
        }
        /* For lane 1 ensure a guaranteed non-solution when phase_oracle by flipping bit 0 of planted_sol */
        if (phase_oracle && inst->planted_sol_bits == in_w && in_w > 0u) {
            uint32_t q0 = plan.val_allocs[in_vid].qubits[0];
            uint64_t b0 = (init_q[q0] & 1u) ^ 1u;
            init_q[q0] = (init_q[q0] & ~UINT64_C(2)) | (b0 << 1u);
            for (uint16_t b = 1u; b < in_w; ++b) {
                uint32_t qb = plan.val_allocs[in_vid].qubits[b];
                uint64_t bb = (init_q[qb] & 1u);
                init_q[qb] = (init_q[qb] & ~UINT64_C(2)) | (bb << 1u);
            }
        }

        flow_u128 p1[64], p2[64];
        bool ok_sim = qir_sim_kmb_64_ex(&kmb, init_q, init_b, UINT64_C(0x123456789ABCDEF0), out_q1, out_b1, p1) &&
                      qir_sim_kmb_64_ex(&kmb, init_q, init_b, UINT64_C(0x0FEDCBA987654321), out_q2, out_b2, p2);
        bool kat_ok = ok_sim;
        const flow_u128 pi_turn = (flow_u128)1u << 127u;
        flow_u128 base_nonhit = 0u;
        bool has_nonhit = false;

        for (uint32_t q = 0u; q < k_nq && kat_ok; ++q) {
            if (out_q1[q] != out_q2[q]) kat_ok = false;
            if (!io_q[q] && out_q1[q] != 0u) kat_ok = false;
            if (phase_oracle && io_q[q] && out_q1[q] != init_q[q]) kat_ok = false;
        }
        for (uint32_t l = 0u; l < 64u && kat_ok; ++l) {
            extra_states++;
            if (p1[l] != p2[l]) { kat_ok = false; break; }
            uint8_t ibuf[32];
            pqc__unpack_lanes_bytes(init_q, plan.val_allocs[in_vid].qubits, in_w, l, ibuf);
            if (phase_oracle) {
                bool is_hit = pqc__eval_planted_hit(inst, ibuf);
                if (!is_hit && !has_nonhit) {
                    base_nonhit = p1[l];
                    has_nonhit  = true;
                }
            }
        }
        if (phase_oracle && kat_ok) {
            bool saw_hit = false;
            for (uint32_t l = 0u; l < 64u; ++l) {
                uint8_t ibuf[32];
                pqc__unpack_lanes_bytes(init_q, plan.val_allocs[in_vid].qubits, in_w, l, ibuf);
                bool is_hit = pqc__eval_planted_hit(inst, ibuf);
                if (is_hit) saw_hit = true;
                flow_u128 expect_diff = is_hit ? pi_turn : 0u;
                if ((p1[l] - base_nonhit) != expect_diff) {
                    fprintf(stderr, "[FAIL] %s: lane %u phase mismatch (hit=%d, p1_hi=0x%016llx, base_hi=0x%016llx)\n",
                            label, l, (int)is_hit,
                            (unsigned long long)(p1[l] >> 64u),
                            (unsigned long long)(base_nonhit >> 64u));
                    kat_ok = false;
                    break;
                }
            }
            if (!saw_hit || !has_nonhit) {
                fprintf(stderr, "[FAIL] %s: planted test did not exercise both hit & non-hit\n", label);
                kat_ok = false;
            }
        } else if (!phase_oracle && kat_ok && mod->num_outputs == 1u) {
            uint32_t out_vid = mod->outputs[0];
            uint16_t out_w   = mod->vals[out_vid].width;
            for (uint32_t l = 0u; l < 64u && kat_ok; ++l) {
                uint8_t ibuf[32], obuf[32];
                pqc__unpack_lanes_bytes(init_q, plan.val_allocs[in_vid].qubits, in_w, l, ibuf);
                memset(obuf, 0, sizeof(obuf));
                for (uint16_t b = 0u; b < out_w; ++b) {
                    if ((out_q1[plan.val_allocs[out_vid].qubits[b]] >> l) & 1u) pqc__set_bit(obuf, b, 1u);
                }
                if (inst->kind == PQC_KAT_AES) {
                    uint8_t ref_ct[16];
                    pqc_ref_aes_encrypt(inst->pt, ibuf, inst->key_bits, inst->rounds, ref_ct);
                    if (memcmp(obuf, ref_ct, 16u) != 0) kat_ok = false;
                } else if (inst->kind == PQC_KAT_SHA256) {
                    uint8_t ref_dig[32];
                    pqc_ref_sha256_msg256(ibuf, inst->rounds, ref_dig);
                    if (memcmp(obuf, ref_dig, 32u) != 0) kat_ok = false;
                } else if (inst->kind == PQC_KAT_ASCON) {
                    uint64_t s[5] = {
                        inst->ascon_iv_n[0], pqc__pack_u64_le(ibuf), pqc__pack_u64_le(ibuf + 8),
                        inst->ascon_iv_n[1], inst->ascon_iv_n[2]
                    };
                    pqc_ref_ascon_p(s, inst->rounds);
                    if (pqc__pack_u64_le(obuf) != s[0] || pqc__pack_u64_le(obuf + 8) != s[1]) kat_ok = false;
                } else if (inst->kind == PQC_KAT_KECCAK400) {
                    uint16_t st[25];
                    memset(st, 0, sizeof(st));
                    for (uint32_t i = 0u; i < 8u; ++i) {
                        st[i] = (uint16_t)ibuf[2u * i] | ((uint16_t)ibuf[2u * i + 1u] << 8u);
                    }
                    pqc_ref_keccak_f400(st, inst->rounds);
                    if (memcmp(obuf, st, 16u) != 0) kat_ok = false;
                } else if (inst->kind == PQC_KAT_MCELIECE_GOPPA) {
                    uint16_t alpha = (uint16_t)ibuf[0] | ((uint16_t)(ibuf[1] & 0x0Fu) << 8u);
                    uint16_t got = (uint16_t)obuf[0] | ((uint16_t)(obuf[1] & 0x0Fu) << 8u);
                    if (got != pqc_ref_mceliece_goppa_inv(alpha, inst->goppa_coeffs, inst->goppa_poly12)) kat_ok = false;
                } else if (inst->kind == PQC_KAT_MLKEM_LWE) {
                    uint32_t sp = (uint32_t)ibuf[0] | ((uint32_t)ibuf[1] << 8u) | ((uint32_t)ibuf[2] << 16u);
                    uint16_t ref_res[4];
                    (void)pqc_ref_mlkem_lwe(sp, inst->lwe_a, inst->lwe_b, ref_res);
                    for (uint32_t i = 0u; i < 4u; ++i) {
                        uint16_t got_r = 0u;
                        for (uint16_t b = 0u; b < 12u; ++b) {
                            if (pqc__get_bit(obuf, 12u * i + b)) got_r |= (uint16_t)(1u << b);
                        }
                        if (got_r != ref_res[i]) { kat_ok = false; break; }
                    }
                } else if (inst->kind == PQC_KAT_MQ24) {
                    uint32_t x24 = (uint32_t)ibuf[0] | ((uint32_t)ibuf[1] << 8u) | ((uint32_t)ibuf[2] << 16u);
                    uint32_t got = (uint32_t)obuf[0] | ((uint32_t)obuf[1] << 8u) | ((uint32_t)obuf[2] << 16u);
                    if (got != pqc_ref_mq24(x24, &inst->mq_sys)) kat_ok = false;
                }
            }
        }

        free(init_q); free(out_q1); free(out_q2);
        free(init_b); free(out_b1); free(out_b2); free(io_q);
        if (!kat_ok) {
            fprintf(stderr, "[FAIL] %s: classical KAT/phase verification failed\n", label);
            qcc_kmb_free(&kmb);
            qir_phys_plan_free(&plan);
            qir_egraph_free(&eg);
            return false;
        }
    }

    if (kmb_out_path_opt) {
        qcc_kmb_fiber_t wf = qcc_kmb_write(&kmb, kmb_out_path_opt);
        if (wf != QCC_KMB_FIBER_COMMITTED) {
            fprintf(stderr, "[FAIL] %s: write %s failed (%s, code=%d)\n",
                    label, kmb_out_path_opt, qcc_kmb_fiber_name(wf), (int)wf);
            qcc_kmb_free(&kmb);
            qir_phys_plan_free(&plan);
            qir_egraph_free(&eg);
            return false;
        }
    }

    pqc_nist_report_t rep;
    pqc_nist_estimate(label, family, search_bits, &ks, &rep);
    if (out_rep) *out_rep = rep;
    if (out_verified_states) *out_verified_states = vrep.states_tested + extra_states;
    if (out_compile_sec) *out_compile_sec = pqc__now_sec() - t0;

    qcc_kmb_free(&kmb);
    qir_phys_plan_free(&plan);
    qir_egraph_free(&eg);
    return true;
}

static void pqc_init_sample_instances(
    uint8_t aes_pt[16], uint8_t aes_key[32],
    uint8_t sha_msg[32],
    uint64_t ascon_iv_n[3], uint8_t ascon_key[16],
    uint8_t keccak_seed[16],
    uint16_t goppa_coeffs[9],
    uint64_t isd_h[24],
    uint16_t lwe_a[4][8], uint16_t lwe_b[4], uint8_t lwe_sol[3],
    pqc_mq24_sys_t *mq_sys
) {
    /* FIPS-197 Appendix C keys & plaintext */
    for (uint32_t i = 0u; i < 16u; ++i) aes_pt[i] = (uint8_t)(i * 0x11u);
    for (uint32_t i = 0u; i < 32u; ++i) aes_key[i] = (uint8_t)i;

    for (uint32_t i = 0u; i < 32u; ++i) sha_msg[i] = (uint8_t)(0xA5u ^ (i * 7u));

    ascon_iv_n[0] = UINT64_C(0x00001000808C0001);
    ascon_iv_n[1] = UINT64_C(0x0123456789ABCDEF);
    ascon_iv_n[2] = UINT64_C(0xFEDCBA9876543210);
    for (uint32_t i = 0u; i < 16u; ++i) ascon_key[i] = (uint8_t)(0x10u + i);

    for (uint32_t i = 0u; i < 16u; ++i) keccak_seed[i] = (uint8_t)(0x3Cu + 5u * i);

    /* Degree-8 monic Goppa polynomial over GF(2^12) = 0x1009 */
    static const uint16_t g_init[9] = {0x317u, 0x9A1u, 0x05Eu, 0x7D4u, 0x208u, 0xB43u, 0x411u, 0x8F2u, 0x001u};
    memcpy(goppa_coeffs, g_init, sizeof(g_init));

    uint64_t rng = UINT64_C(0x51C0FFEE20260001);
    for (uint32_t r = 0u; r < 24u; ++r) isd_h[r] = qcc_splitmix64(&rng);

    /* ML-KEM q=3329 LWE instance with planted secret s_j in {-2,-1,0,1,2} encoded as (sv + 2) */
    uint32_t s_pack = 0u;
    for (uint32_t j = 0u; j < 8u; ++j) {
        uint32_t code = (j * 3u + 1u) % 5u; /* strictly in 0..4 */
        s_pack |= (code << (3u * j));
    }
    lwe_sol[0] = (uint8_t)(s_pack & 0xFFu);
    lwe_sol[1] = (uint8_t)((s_pack >> 8u) & 0xFFu);
    lwe_sol[2] = (uint8_t)((s_pack >> 16u) & 0xFFu);
    for (uint32_t i = 0u; i < 4u; ++i) {
        int32_t acc = 0;
        for (uint32_t j = 0u; j < 8u; ++j) {
            lwe_a[i][j] = (uint16_t)(qcc_splitmix64(&rng) % PQC_KYBER_Q);
            int32_t sv = (int32_t)((s_pack >> (3u * j)) & 7u) - 2;
            acc += (int32_t)lwe_a[i][j] * sv;
        }
        int32_t noise = (int32_t)(i % 5u) - 2; /* bounded in [-2, +2] */
        int32_t bv = (acc + noise) % (int32_t)PQC_KYBER_Q;
        if (bv < 0) bv += (int32_t)PQC_KYBER_Q;
        lwe_b[i] = (uint16_t)bv;
    }

    memset(mq_sys, 0, sizeof(*mq_sys));
    mq_sys->const_bits = (uint32_t)(qcc_splitmix64(&rng) & UINT32_C(0xFFFFFF));
    for (uint32_t k = 0u; k < 24u; ++k) {
        mq_sys->lin_coeff[k] = (uint32_t)(qcc_splitmix64(&rng) & UINT32_C(0xFFFFFF));
        for (uint32_t d = 0u; d < 12u; ++d) {
            uint32_t mask = (d == 11u) ? UINT32_C(0xFFF) : UINT32_C(0xFFFFFF);
            mq_sys->quad_coeff[k][d] = (uint32_t)(qcc_splitmix64(&rng) & mask);
        }
    }
}

static bool pqc_run_classical_kats(void) {
    /* 1. FIPS-197 AES-128 / 192 / 256 Appendix C vectors */
    uint8_t pt[16], key[32], ct[16];
    for (uint32_t i = 0u; i < 16u; ++i) pt[i] = (uint8_t)(i * 0x11u);
    for (uint32_t i = 0u; i < 32u; ++i) key[i] = (uint8_t)i;

    static const uint8_t KAT_AES128[16] = {
        0x69,0xc4,0xe0,0xd8,0x6a,0x7b,0x04,0x30,0xd8,0xcd,0xb7,0x80,0x70,0xb4,0xc5,0x5a
    };
    static const uint8_t KAT_AES192[16] = {
        0xdd,0xa9,0x7c,0xa4,0x86,0x4c,0xdf,0xe0,0x6e,0xaf,0x70,0xa0,0xec,0x0d,0x71,0x91
    };
    static const uint8_t KAT_AES256[16] = {
        0x8e,0xa2,0xb7,0xca,0x51,0x67,0x45,0xbf,0xea,0xfc,0x49,0x90,0x4b,0x49,0x60,0x89
    };
    pqc_ref_aes_encrypt(pt, key, 128u, 10u, ct);
    if (memcmp(ct, KAT_AES128, 16u) != 0) { fprintf(stderr, "KAT AES-128 failed\n"); return false; }
    pqc_ref_aes_encrypt(pt, key, 192u, 12u, ct);
    if (memcmp(ct, KAT_AES192, 16u) != 0) { fprintf(stderr, "KAT AES-192 failed\n"); return false; }
    pqc_ref_aes_encrypt(pt, key, 256u, 14u, ct);
    if (memcmp(ct, KAT_AES256, 16u) != 0) { fprintf(stderr, "KAT AES-256 failed\n"); return false; }

    /* 2. FIPS 180-4 SHA-256("abc") */
    static const uint8_t KAT_SHA_ABC[32] = {
        0xba,0x78,0x16,0xbf,0x8f,0x01,0xcf,0xea,0x41,0x41,0x40,0xde,0x5d,0xae,0x22,0x23,
        0xb0,0x03,0x61,0xa3,0x96,0x17,0x7a,0x9c,0xb4,0x10,0xff,0x61,0xf2,0x00,0x15,0xad
    };
    uint8_t dig[32];
    pqc_ref_sha256_bytes((const uint8_t *)"abc", 3u, dig);
    if (memcmp(dig, KAT_SHA_ABC, 32u) != 0) { fprintf(stderr, "KAT SHA-256('abc') failed\n"); return false; }

    /* 3. Ascon-p[12] KAT (v1.2 Ascon-Hash IV = 0x00400c0000000100 || 0^256 -> S0 = 0xee9398aadb67f03d) */
    uint64_t as[5] = {UINT64_C(0x00400c0000000100), 0, 0, 0, 0};
    pqc_ref_ascon_p(as, 12u);
    if (as[0] != UINT64_C(0xee9398aadb67f03d)) {
        fprintf(stderr, "KAT Ascon-p[12] failed: got 0x%016llx\n", (unsigned long long)as[0]);
        return false;
    }

    /* 4. Fixed-point log2 & NIST MAXDEPTH=2^64 category sanity */
    uint32_t l8 = pqc_log2_q16_u64(8u);
    if (l8 != (3u << 16u)) { fprintf(stderr, "log2(8) check failed\n"); return false; }
    if (pqc_classify_nist_cat(91u << 16u) != 1u ||
        pqc_classify_nist_cat(156u << 16u) != 3u ||
        pqc_classify_nist_cat(220u << 16u) != 5u) {
        fprintf(stderr, "NIST category classification sanity check failed\n");
        return false;
    }
    return true;
}

static int cmd_selftest(int argc, char **argv) {
    const char *sample_path = NULL;
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: pqc selftest [--emit-sample <path.kmb>]\n");
            return 0;
        } else if (strcmp(argv[i], "--emit-sample") == 0 && i + 1 < argc) {
            sample_path = argv[++i];
        } else {
            fprintf(stderr, "Error: unknown or incomplete option '%s'\n", argv[i]);
            return 2;
        }
    }
    printf("=== PQC Quantum Cryptanalysis Compiler Selftest ===\n");
    if (!pqc_run_classical_kats()) return 1;
    printf("  KAT_FIPS197_FIPS180_SP800232       | AES-128/192/256 + SHA-256 + Ascon-p[12] | PASS\n");

    uint8_t aes_pt[16], aes_key[32], sha_msg[32], ascon_key[16], keccak_seed[16], lwe_sol[3];
    uint64_t ascon_iv_n[3], isd_h[24];
    uint16_t goppa_coeffs[9], lwe_a[4][8], lwe_b[4];
    pqc_mq24_sys_t mq_sys;
    pqc_init_sample_instances(aes_pt, aes_key, sha_msg, ascon_iv_n, ascon_key, keccak_seed,
                              goppa_coeffs, isd_h, lwe_a, lwe_b, lwe_sol, &mq_sys);

    int fails = 0;
    qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
    cfg.require_clean_scratch = true;
    cfg.cost_mode = QCC_COST_QT;

    /* Helper macro-like loop over compact verified forward & Grover oracle tests */
    {
        /* 1. AES-128 Forward (10R, emits build/selftest_aes128_fwd.kmb for qcc verify) & Grover (2R) */
        uint8_t aes_ct10[16], aes_ct2[16];
        pqc_ref_aes_encrypt(aes_pt, aes_key, 128u, 10u, aes_ct10);
        pqc_ref_aes_encrypt(aes_pt, aes_key, 128u, 2u, aes_ct2);
        for (int pass = 0; pass < 2; ++pass) {
            bool is_ph = (pass == 1);
            pqc_instance_t inst;
            memset(&inst, 0, sizeof(inst));
            inst.kind             = PQC_KAT_AES;
            inst.rounds           = is_ph ? 2u : 10u;
            inst.key_bits         = 128u;
            inst.planted_sol_bits = 128u;
            memcpy(inst.pt, aes_pt, 16u);
            memcpy(inst.planted_sol, aes_key, 16u);
            memcpy(inst.target_bytes, is_ph ? aes_ct2 : aes_ct10, 16u);

            qcc_arena_t ar; qcc_arena_init(&ar, 32u * 1024u * 1024u);
            qir_module_t m; qir_module_init(&m, &ar, is_ph ? "aes128_2r_ph" : "aes128_10r_fwd");
            m.require_clean_scratch = true;
            pqc_build_aes_qir(&m, 128u, inst.rounds, aes_pt, inst.target_bytes, is_ph);
            pqc_nist_report_t rep;
            uint64_t states = 0u;
            double sec = 0.0;
            const char *path = is_ph ? NULL : sample_path;
            if (!pqc_compile_and_verify(is_ph ? "AES128_2R_Grover" : "AES128_10R_Fwd", "AES", 128u,
                                        &m, &ar, &cfg, &inst, is_ph, path, 2u, &rep, &states, &sec)) {
                fails++;
            } else {
                printf("  %-34s | Q=%5u CCX=%6llu D_t=%5llu | PASS (%llu states, %.2fs)\n",
                       rep.name, rep.oracle_stats.num_qubits,
                       (unsigned long long)rep.oracle_stats.toffoli,
                       (unsigned long long)rep.oracle_stats.toffoli_depth,
                       (unsigned long long)states, sec);
            }
            qir_module_free(&m); qcc_arena_free(&ar);
        }
    }

    {
        /* 2. AES-192 (2R Fwd) & AES-256 (2R Fwd) verifying KeySchedule variants */
        const uint32_t kbits[2] = {192u, 256u};
        for (int ki = 0; ki < 2; ++ki) {
            pqc_instance_t inst;
            memset(&inst, 0, sizeof(inst));
            inst.kind             = PQC_KAT_AES;
            inst.rounds           = 2u;
            inst.key_bits         = kbits[ki];
            inst.planted_sol_bits = kbits[ki];
            memcpy(inst.pt, aes_pt, 16u);
            memcpy(inst.planted_sol, aes_key, kbits[ki] / 8u);
            pqc_ref_aes_encrypt(aes_pt, aes_key, kbits[ki], 2u, inst.target_bytes);

            qcc_arena_t ar; qcc_arena_init(&ar, 16u * 1024u * 1024u);
            qir_module_t m; qir_module_init(&m, &ar, ki == 0 ? "aes192_2r" : "aes256_2r");
            m.require_clean_scratch = true;
            pqc_build_aes_qir(&m, kbits[ki], 2u, aes_pt, inst.target_bytes, true);
            pqc_nist_report_t rep;
            uint64_t states = 0u;
            double sec = 0.0;
            const char *lbl = (ki == 0) ? "AES192_2R_Grover" : "AES256_2R_Grover";
            if (!pqc_compile_and_verify(lbl, "AES", kbits[ki], &m, &ar, &cfg, &inst, true, NULL, 2u, &rep, &states, &sec)) {
                fails++;
            } else {
                printf("  %-34s | Q=%5u CCX=%6llu D_t=%5llu | PASS (%llu states, %.2fs)\n",
                       rep.name, rep.oracle_stats.num_qubits,
                       (unsigned long long)rep.oracle_stats.toffoli,
                       (unsigned long long)rep.oracle_stats.toffoli_depth,
                       (unsigned long long)states, sec);
            }
            qir_module_free(&m); qcc_arena_free(&ar);
        }
    }

    {
        /* 3. SHA-256 (64R full Grover) */
        pqc_instance_t inst;
        memset(&inst, 0, sizeof(inst));
        inst.kind             = PQC_KAT_SHA256;
        inst.rounds           = 64u;
        inst.planted_sol_bits = 256u;
        memcpy(inst.planted_sol, sha_msg, 32u);
        pqc_ref_sha256_msg256(sha_msg, 64u, inst.target_bytes);

        qcc_arena_t ar; qcc_arena_init(&ar, 64u * 1024u * 1024u);
        qir_module_t m; qir_module_init(&m, &ar, "sha256_64r");
        m.require_clean_scratch = true;
        pqc_build_sha256_qir(&m, 64u, inst.target_bytes, true);
        pqc_nist_report_t rep;
        uint64_t states = 0u;
        double sec = 0.0;
        if (!pqc_compile_and_verify("SHA256_64R_Grover", "Hash", 256u, &m, &ar, &cfg, &inst, true, NULL, 2u, &rep, &states, &sec)) {
            fails++;
        } else {
            printf("  %-34s | Q=%5u CCX=%6llu D_t=%5llu | PASS (%llu states, %.2fs)\n",
                   rep.name, rep.oracle_stats.num_qubits,
                   (unsigned long long)rep.oracle_stats.toffoli,
                   (unsigned long long)rep.oracle_stats.toffoli_depth,
                   (unsigned long long)states, sec);
        }
        qir_module_free(&m); qcc_arena_free(&ar);
    }

    {
        /* 4. Ascon-p[12] (12R Grover) */
        pqc_instance_t inst;
        memset(&inst, 0, sizeof(inst));
        inst.kind             = PQC_KAT_ASCON;
        inst.rounds           = 12u;
        inst.planted_sol_bits = 128u;
        memcpy(inst.ascon_iv_n, ascon_iv_n, sizeof(ascon_iv_n));
        memcpy(inst.planted_sol, ascon_key, 16u);
        uint64_t s[5] = {ascon_iv_n[0], pqc__pack_u64_le(ascon_key), pqc__pack_u64_le(ascon_key + 8), ascon_iv_n[1], ascon_iv_n[2]};
        pqc_ref_ascon_p(s, 12u);
        inst.ascon_tag[0] = s[0]; inst.ascon_tag[1] = s[1];

        qcc_arena_t ar; qcc_arena_init(&ar, 32u * 1024u * 1024u);
        qir_module_t m; qir_module_init(&m, &ar, "ascon12");
        m.require_clean_scratch = true;
        pqc_build_ascon_qir(&m, 12u, ascon_iv_n, inst.ascon_tag, true);
        pqc_nist_report_t rep;
        uint64_t states = 0u;
        double sec = 0.0;
        if (!pqc_compile_and_verify("Ascon128_12R_Grover", "Sponge", 128u, &m, &ar, &cfg, &inst, true, NULL, 2u, &rep, &states, &sec)) {
            fails++;
        } else {
            printf("  %-34s | Q=%5u CCX=%6llu D_t=%5llu | PASS (%llu states, %.2fs)\n",
                   rep.name, rep.oracle_stats.num_qubits,
                   (unsigned long long)rep.oracle_stats.toffoli,
                   (unsigned long long)rep.oracle_stats.toffoli_depth,
                   (unsigned long long)states, sec);
        }
        qir_module_free(&m); qcc_arena_free(&ar);
    }

    {
        /* 5. Keccak-f[400] (20R Grover) */
        pqc_instance_t inst;
        memset(&inst, 0, sizeof(inst));
        inst.kind             = PQC_KAT_KECCAK400;
        inst.rounds           = 20u;
        inst.planted_sol_bits = 128u;
        memcpy(inst.planted_sol, keccak_seed, 16u);
        uint16_t kst[25];
        memset(kst, 0, sizeof(kst));
        for (uint32_t i = 0u; i < 8u; ++i) {
            kst[i] = (uint16_t)keccak_seed[2u * i] | ((uint16_t)keccak_seed[2u * i + 1u] << 8u);
        }
        pqc_ref_keccak_f400(kst, 20u);
        memcpy(inst.keccak_rate, kst, sizeof(inst.keccak_rate));

        qcc_arena_t ar; qcc_arena_init(&ar, 64u * 1024u * 1024u);
        qir_module_t m; qir_module_init(&m, &ar, "keccak400");
        m.require_clean_scratch = true;
        pqc_build_keccak_f400_qir(&m, 20u, inst.keccak_rate, true);
        pqc_nist_report_t rep;
        uint64_t states = 0u;
        double sec = 0.0;
        if (!pqc_compile_and_verify("Keccak_f400_20R_Grover", "Sponge", 128u, &m, &ar, &cfg, &inst, true, NULL, 2u, &rep, &states, &sec)) {
            fails++;
        } else {
            printf("  %-34s | Q=%5u CCX=%6llu D_t=%5llu | PASS (%llu states, %.2fs)\n",
                   rep.name, rep.oracle_stats.num_qubits,
                   (unsigned long long)rep.oracle_stats.toffoli,
                   (unsigned long long)rep.oracle_stats.toffoli_depth,
                   (unsigned long long)states, sec);
        }
        qir_module_free(&m); qcc_arena_free(&ar);
    }

    {
        /* 6. Classic McEliece Goppa Inv (GF(2^12) 0x1009) + ISD */
        pqc_instance_t inst;
        memset(&inst, 0, sizeof(inst));
        inst.kind             = PQC_KAT_MCELIECE_GOPPA;
        inst.goppa_poly12     = 0x1009u;
        memcpy(inst.goppa_coeffs, goppa_coeffs, sizeof(goppa_coeffs));
        uint16_t alpha_sol    = 0x7A3u;
        inst.goppa_target_inv = pqc_ref_mceliece_goppa_inv(alpha_sol, goppa_coeffs, 0x1009u);
        inst.planted_sol[0]   = (uint8_t)(alpha_sol & 0xFFu);
        inst.planted_sol[1]   = (uint8_t)((alpha_sol >> 8u) & 0x0Fu);
        inst.planted_sol_bits = 12u;

        qcc_arena_t ar; qcc_arena_init(&ar, 4u * 1024u * 1024u);
        qir_module_t m; qir_module_init(&m, &ar, "mceliece_goppa");
        m.require_clean_scratch = true;
        pqc_build_mceliece_goppa_qir(&m, goppa_coeffs, 0x1009u, inst.goppa_target_inv, true);
        pqc_nist_report_t rep;
        uint64_t states = 0u;
        double sec = 0.0;
        if (!pqc_compile_and_verify("McEliece_Goppa_GF12_0x1009", "Code", 12u, &m, &ar, &cfg, &inst, true, NULL, 4u, &rep, &states, &sec)) {
            fails++;
        } else {
            printf("  %-34s | Q=%5u CCX=%6llu D_t=%5llu | PASS (%llu states, %.2fs)\n",
                   rep.name, rep.oracle_stats.num_qubits,
                   (unsigned long long)rep.oracle_stats.toffoli,
                   (unsigned long long)rep.oracle_stats.toffoli_depth,
                   (unsigned long long)states, sec);
        }
        qir_module_free(&m); qcc_arena_free(&ar);
    }

    {
        pqc_instance_t inst;
        memset(&inst, 0, sizeof(inst));
        inst.kind             = PQC_KAT_MCELIECE_ISD;
        memcpy(inst.isd_h, isd_h, sizeof(isd_h));
        uint64_t err_sol      = UINT64_C(0x0010080204008011);
        pqc_ref_mceliece_isd(err_sol, isd_h, &inst.isd_syn, &inst.isd_wt);
        pqc__unpack_u64_le(err_sol, inst.planted_sol);
        inst.planted_sol_bits = 64u;

        qcc_arena_t ar; qcc_arena_init(&ar, 8u * 1024u * 1024u);
        qir_module_t m; qir_module_init(&m, &ar, "mceliece_isd");
        m.require_clean_scratch = true;
        pqc_build_mceliece_isd_qir(&m, isd_h, inst.isd_syn, inst.isd_wt, true);
        pqc_nist_report_t rep;
        uint64_t states = 0u;
        double sec = 0.0;
        if (!pqc_compile_and_verify("McEliece_ISD_64b_SynWt", "Code", 64u, &m, &ar, &cfg, &inst, true, NULL, 2u, &rep, &states, &sec)) {
            fails++;
        } else {
            printf("  %-34s | Q=%5u CCX=%6llu D_t=%5llu | PASS (%llu states, %.2fs)\n",
                   rep.name, rep.oracle_stats.num_qubits,
                   (unsigned long long)rep.oracle_stats.toffoli,
                   (unsigned long long)rep.oracle_stats.toffoli_depth,
                   (unsigned long long)states, sec);
        }
        qir_module_free(&m); qcc_arena_free(&ar);
    }

    {
        /* 7. ML-KEM (q=3329 Module-LWE Bounded Search) */
        pqc_instance_t inst;
        memset(&inst, 0, sizeof(inst));
        inst.kind             = PQC_KAT_MLKEM_LWE;
        memcpy(inst.lwe_a, lwe_a, sizeof(lwe_a));
        memcpy(inst.lwe_b, lwe_b, sizeof(lwe_b));
        inst.lwe_eta          = 2u;
        memcpy(inst.planted_sol, lwe_sol, 3u);
        inst.planted_sol_bits = 24u;

        qcc_arena_t ar; qcc_arena_init(&ar, 8u * 1024u * 1024u);
        qir_module_t m; qir_module_init(&m, &ar, "mlkem_lwe");
        m.require_clean_scratch = true;
        pqc_build_mlkem_lwe_qir(&m, lwe_a, lwe_b, 2u, true);
        pqc_nist_report_t rep;
        uint64_t states = 0u;
        double sec = 0.0;
        if (!pqc_compile_and_verify("MLKEM_q3329_LWE_Grover", "Lattice", 24u, &m, &ar, &cfg, &inst, true, NULL, 2u, &rep, &states, &sec)) {
            fails++;
        } else {
            printf("  %-34s | Q=%5u CCX=%6llu D_t=%5llu | PASS (%llu states, %.2fs)\n",
                   rep.name, rep.oracle_stats.num_qubits,
                   (unsigned long long)rep.oracle_stats.toffoli,
                   (unsigned long long)rep.oracle_stats.toffoli_depth,
                   (unsigned long long)states, sec);
        }
        qir_module_free(&m); qcc_arena_free(&ar);
    }

    {
        /* 8. Multivariate Quadratic (MQ n=24, m=24) */
        pqc_instance_t inst;
        memset(&inst, 0, sizeof(inst));
        inst.kind             = PQC_KAT_MQ24;
        inst.mq_sys           = mq_sys;
        uint32_t x_sol        = UINT32_C(0x5A3C9F);
        inst.mq_target        = pqc_ref_mq24(x_sol, &mq_sys);
        inst.planted_sol[0]   = (uint8_t)(x_sol & 0xFFu);
        inst.planted_sol[1]   = (uint8_t)((x_sol >> 8u) & 0xFFu);
        inst.planted_sol[2]   = (uint8_t)((x_sol >> 16u) & 0xFFu);
        inst.planted_sol_bits = 24u;

        qcc_arena_t ar; qcc_arena_init(&ar, 4u * 1024u * 1024u);
        qir_module_t m; qir_module_init(&m, &ar, "mq24");
        m.require_clean_scratch = true;
        pqc_build_mq24_qir(&m, &mq_sys, inst.mq_target, true);
        pqc_nist_report_t rep;
        uint64_t states = 0u;
        double sec = 0.0;
        if (!pqc_compile_and_verify("MQ_n24_m24_Grover", "MQ", 24u, &m, &ar, &cfg, &inst, true, NULL, 2u, &rep, &states, &sec)) {
            fails++;
        } else {
            printf("  %-34s | Q=%5u CCX=%6llu D_t=%5llu | PASS (%llu states, %.2fs)\n",
                   rep.name, rep.oracle_stats.num_qubits,
                   (unsigned long long)rep.oracle_stats.toffoli,
                   (unsigned long long)rep.oracle_stats.toffoli_depth,
                   (unsigned long long)states, sec);
        }
        qir_module_free(&m); qcc_arena_free(&ar);
    }

    printf("=== PQC Selftest Summary: %s (fails=%d) ===\n", fails == 0 ? "ALL PASS" : "FAILED", fails);
    return fails == 0 ? 0 : 1;
}

static void pqc_write_json_report(const char *path, const pqc_nist_report_t *reps, uint32_t n) {
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "Error: cannot write JSON to %s\n", path); return; }
    fprintf(f, "{\n  \"benchmarks\": [\n");
    for (uint32_t i = 0u; i < n; ++i) {
        const pqc_nist_report_t *r = &reps[i];
        char sb[32], md40t[32], md64t[32], md96t[32], md64q[32], md64gd[32];
        pqc_format_pow2_q16(r->log2_serial_t_q16, sb, sizeof(sb));
        pqc_format_pow2_q16(r->maxdepth[0].log2_total_t_q16, md40t, sizeof(md40t));
        pqc_format_pow2_q16(r->maxdepth[1].log2_total_t_q16, md64t, sizeof(md64t));
        pqc_format_pow2_q16(r->maxdepth[2].log2_total_t_q16, md96t, sizeof(md96t));
        pqc_format_pow2_q16(r->maxdepth[1].log2_total_qubits_q16, md64q, sizeof(md64q));
        pqc_format_pow2_q16(r->maxdepth[1].log2_g_times_d_q16, md64gd, sizeof(md64gd));
        fprintf(f,
                "    {\n"
                "      \"name\": \"%s\",\n"
                "      \"family\": \"%s\",\n"
                "      \"search_bits\": %u,\n"
                "      \"qubits\": %u,\n"
                "      \"cbits\": %u,\n"
                "      \"toffoli\": %llu,\n"
                "      \"hmr\": %llu,\n"
                "      \"cnot\": %llu,\n"
                "      \"t_eq_gidney\": %llu,\n"
                "      \"t_eq_unitary\": %llu,\n"
                "      \"toffoli_depth\": %llu,\n"
                "      \"circuit_depth\": %llu,\n"
                "      \"reaction_depth\": %llu,\n"
                "      \"serial_t\": \"%s\",\n"
                "      \"md40_t\": \"%s\",\n"
                "      \"md64_t\": \"%s\",\n"
                "      \"md64_qubits\": \"%s\",\n"
                "      \"md64_gd\": \"%s\",\n"
                "      \"md96_t\": \"%s\",\n"
                "      \"nist_category_equiv\": %u\n"
                "    }%s\n",
                r->name, r->family, (unsigned)r->search_bits,
                (unsigned)r->oracle_stats.num_qubits, (unsigned)r->oracle_stats.num_bits,
                (unsigned long long)r->oracle_stats.toffoli,
                (unsigned long long)r->oracle_stats.hmr,
                (unsigned long long)r->oracle_stats.cnot,
                (unsigned long long)r->t_eq_gidney,
                (unsigned long long)r->t_eq_unitary,
                (unsigned long long)r->oracle_stats.toffoli_depth,
                (unsigned long long)r->oracle_stats.circuit_depth,
                (unsigned long long)r->oracle_stats.reaction_depth,
                sb, md40t, md64t, md64q, md64gd, md96t,
                (unsigned)r->nist_category_equiv,
                (i + 1u < n) ? "," : "");
    }
    fprintf(f, "  ]\n}\n");
    fclose(f);
}

static void pqc_write_md_report(const char *path, const pqc_nist_report_t *reps, uint32_t n) {
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "Error: cannot write Markdown to %s\n", path); return; }
    fprintf(f, "# PQC / NIST Post-Quantum Cryptanalytic Resource Report\n\n");
    fprintf(f, "| Oracle / Benchmark | Family | k | Q | CCX | T-eq (Gidney) | T-depth | Circ-depth | D_react | Serial T | MD 2^40 T | MD 2^64 T | MD 2^64 Q | MD 2^64 G*D | MD 2^96 T | NIST Cat |\n");
    fprintf(f, "| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n");
    for (uint32_t i = 0u; i < n; ++i) {
        const pqc_nist_report_t *r = &reps[i];
        char sb[32], md40t[32], md64t[32], md96t[32], md64q[32], md64gd[32];
        pqc_format_pow2_q16(r->log2_serial_t_q16, sb, sizeof(sb));
        pqc_format_pow2_q16(r->maxdepth[0].log2_total_t_q16, md40t, sizeof(md40t));
        pqc_format_pow2_q16(r->maxdepth[1].log2_total_t_q16, md64t, sizeof(md64t));
        pqc_format_pow2_q16(r->maxdepth[2].log2_total_t_q16, md96t, sizeof(md96t));
        pqc_format_pow2_q16(r->maxdepth[1].log2_total_qubits_q16, md64q, sizeof(md64q));
        pqc_format_pow2_q16(r->maxdepth[1].log2_g_times_d_q16, md64gd, sizeof(md64gd));
        fprintf(f, "| %s | %s | %u | %u | %llu | %llu | %llu | %llu | %llu | %s | %s | %s | %s | %s | %s | Cat %u |\n",
                r->name, r->family, (unsigned)r->search_bits,
                (unsigned)r->oracle_stats.num_qubits,
                (unsigned long long)r->oracle_stats.toffoli,
                (unsigned long long)r->t_eq_gidney,
                (unsigned long long)r->oracle_stats.toffoli_depth,
                (unsigned long long)r->oracle_stats.circuit_depth,
                (unsigned long long)r->oracle_stats.reaction_depth,
                sb, md40t, md64t, md64q, md64gd, md96t,
                (unsigned)r->nist_category_equiv);
    }
    fclose(f);
}

static int cmd_bench(int argc, char **argv) {
    const char *out_dir   = "build";
    const char *json_path = NULL;
    const char *md_path   = NULL;
    bool have_out = false;
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: pqc bench [out_dir] [--json <path>] [--md <path>]\n");
            return 0;
        } else if (strcmp(argv[i], "--json") == 0 && i + 1 < argc) {
            json_path = argv[++i];
        } else if (strcmp(argv[i], "--md") == 0 && i + 1 < argc) {
            md_path = argv[++i];
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "Error: unknown or incomplete option '%s'\n", argv[i]);
            return 2;
        } else if (!have_out) {
            out_dir = argv[i];
            have_out = true;
        }
    }
    (void)mkdir(out_dir, 0755);

    uint8_t aes_pt[16], aes_key[32], sha_msg[32], ascon_key[16], keccak_seed[16], lwe_sol[3];
    uint64_t ascon_iv_n[3], isd_h[24];
    uint16_t goppa_coeffs[9], lwe_a[4][8], lwe_b[4];
    pqc_mq24_sys_t mq_sys;
    pqc_init_sample_instances(aes_pt, aes_key, sha_msg, ascon_iv_n, ascon_key, keccak_seed,
                              goppa_coeffs, isd_h, lwe_a, lwe_b, lwe_sol, &mq_sys);

    qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
    cfg.require_clean_scratch = true;
    cfg.cost_mode = QCC_COST_QT;

    pqc_nist_report_t reps[16];
    uint32_t nr = 0u;
    int fails = 0;

    printf("=== PQC / NIST Quantum Cryptanalysis Resource Benchmark Suite ===\n");
    printf("  %-24s | %6s | %8s | %9s | %7s | %8s | %8s | %8s | %5s\n",
           "Benchmark", "Q", "CCX", "T-eq", "T-depth", "Serial T", "MD64 T", "MD64 Q", "Time");
    printf("  -------------------------+--------+----------+-----------+---------+----------+----------+----------+------\n");

    typedef struct {
        const char *name;
        const char *family;
        uint32_t    search_bits;
        bool        phase_oracle;
        pqc_kat_kind_t kind;
        uint32_t    rounds;
        uint32_t    key_bits;
    } bench_spec_t;

    static const bench_spec_t SPECS[11] = {
        {"aes128_fwd",      "AES-Cat1", 128u, false, PQC_KAT_AES,            10u, 128u},
        {"aes128_grover",   "AES-Cat1", 128u, true,  PQC_KAT_AES,            10u, 128u},
        {"aes192_grover",   "AES-Cat3", 192u, true,  PQC_KAT_AES,            12u, 192u},
        {"aes256_grover",   "AES-Cat5", 256u, true,  PQC_KAT_AES,            14u, 256u},
        {"sha256_grover",   "Hash",     256u, true,  PQC_KAT_SHA256,         64u, 256u},
        {"ascon128_grover", "Sponge",   128u, true,  PQC_KAT_ASCON,          12u, 128u},
        {"keccak400_grover","Sponge",   128u, true,  PQC_KAT_KECCAK400,      20u, 128u},
        {"mceliece_goppa12","Code-GF12", 12u, true,  PQC_KAT_MCELIECE_GOPPA,  8u,  12u},
        {"mceliece_isd",    "Code-ISD",  64u, true,  PQC_KAT_MCELIECE_ISD,    1u,  64u},
        {"mlkem_lwe3329",   "Lattice",   24u, true,  PQC_KAT_MLKEM_LWE,       1u,  24u},
        {"mq24_grover",     "MQ",        24u, true,  PQC_KAT_MQ24,            1u,  24u}
    };

    for (uint32_t si = 0u; si < 11u; ++si) {
        const bench_spec_t *sp = &SPECS[si];
        pqc_instance_t inst;
        memset(&inst, 0, sizeof(inst));
        inst.kind = sp->kind;
        inst.rounds = sp->rounds;
        inst.key_bits = sp->key_bits;

        qcc_arena_t ar; qcc_arena_init(&ar, 64u * 1024u * 1024u);
        qir_module_t m; qir_module_init(&m, &ar, sp->name);
        m.require_clean_scratch = true;

        bool built = false;
        if (sp->kind == PQC_KAT_AES) {
            inst.planted_sol_bits = sp->key_bits;
            memcpy(inst.pt, aes_pt, 16u);
            memcpy(inst.planted_sol, aes_key, sp->key_bits / 8u);
            pqc_ref_aes_encrypt(aes_pt, aes_key, sp->key_bits, sp->rounds, inst.target_bytes);
            built = pqc_build_aes_qir(&m, sp->key_bits, sp->rounds, aes_pt, inst.target_bytes, sp->phase_oracle);
        } else if (sp->kind == PQC_KAT_SHA256) {
            inst.planted_sol_bits = 256u;
            memcpy(inst.planted_sol, sha_msg, 32u);
            pqc_ref_sha256_msg256(sha_msg, sp->rounds, inst.target_bytes);
            built = pqc_build_sha256_qir(&m, sp->rounds, inst.target_bytes, sp->phase_oracle);
        } else if (sp->kind == PQC_KAT_ASCON) {
            inst.planted_sol_bits = 128u;
            memcpy(inst.ascon_iv_n, ascon_iv_n, sizeof(ascon_iv_n));
            memcpy(inst.planted_sol, ascon_key, 16u);
            uint64_t s[5] = {ascon_iv_n[0], pqc__pack_u64_le(ascon_key), pqc__pack_u64_le(ascon_key + 8), ascon_iv_n[1], ascon_iv_n[2]};
            pqc_ref_ascon_p(s, sp->rounds);
            inst.ascon_tag[0] = s[0]; inst.ascon_tag[1] = s[1];
            built = pqc_build_ascon_qir(&m, sp->rounds, ascon_iv_n, inst.ascon_tag, sp->phase_oracle);
        } else if (sp->kind == PQC_KAT_KECCAK400) {
            inst.planted_sol_bits = 128u;
            memcpy(inst.planted_sol, keccak_seed, 16u);
            uint16_t kst[25]; memset(kst, 0, sizeof(kst));
            for (uint32_t i = 0u; i < 8u; ++i) {
                kst[i] = (uint16_t)keccak_seed[2u * i] | ((uint16_t)keccak_seed[2u * i + 1u] << 8u);
            }
            pqc_ref_keccak_f400(kst, sp->rounds);
            memcpy(inst.keccak_rate, kst, sizeof(inst.keccak_rate));
            built = pqc_build_keccak_f400_qir(&m, sp->rounds, inst.keccak_rate, sp->phase_oracle);
        } else if (sp->kind == PQC_KAT_MCELIECE_GOPPA) {
            inst.goppa_poly12 = 0x1009u;
            memcpy(inst.goppa_coeffs, goppa_coeffs, sizeof(goppa_coeffs));
            uint16_t alpha_sol = 0x7A3u;
            inst.goppa_target_inv = pqc_ref_mceliece_goppa_inv(alpha_sol, goppa_coeffs, 0x1009u);
            inst.planted_sol[0] = (uint8_t)(alpha_sol & 0xFFu);
            inst.planted_sol[1] = (uint8_t)((alpha_sol >> 8u) & 0x0Fu);
            inst.planted_sol_bits = 12u;
            built = pqc_build_mceliece_goppa_qir(&m, goppa_coeffs, 0x1009u, inst.goppa_target_inv, sp->phase_oracle);
        } else if (sp->kind == PQC_KAT_MCELIECE_ISD) {
            memcpy(inst.isd_h, isd_h, sizeof(isd_h));
            uint64_t err_sol = UINT64_C(0x0010080204008011);
            pqc_ref_mceliece_isd(err_sol, isd_h, &inst.isd_syn, &inst.isd_wt);
            pqc__unpack_u64_le(err_sol, inst.planted_sol);
            inst.planted_sol_bits = 64u;
            built = pqc_build_mceliece_isd_qir(&m, isd_h, inst.isd_syn, inst.isd_wt, sp->phase_oracle);
        } else if (sp->kind == PQC_KAT_MLKEM_LWE) {
            memcpy(inst.lwe_a, lwe_a, sizeof(lwe_a));
            memcpy(inst.lwe_b, lwe_b, sizeof(lwe_b));
            inst.lwe_eta = 2u;
            memcpy(inst.planted_sol, lwe_sol, 3u);
            inst.planted_sol_bits = 24u;
            built = pqc_build_mlkem_lwe_qir(&m, lwe_a, lwe_b, 2u, sp->phase_oracle);
        } else if (sp->kind == PQC_KAT_MQ24) {
            inst.mq_sys = mq_sys;
            uint32_t x_sol = UINT32_C(0x5A3C9F);
            inst.mq_target = pqc_ref_mq24(x_sol, &mq_sys);
            inst.planted_sol[0] = (uint8_t)(x_sol & 0xFFu);
            inst.planted_sol[1] = (uint8_t)((x_sol >> 8u) & 0xFFu);
            inst.planted_sol[2] = (uint8_t)((x_sol >> 16u) & 0xFFu);
            inst.planted_sol_bits = 24u;
            built = pqc_build_mq24_qir(&m, &mq_sys, inst.mq_target, sp->phase_oracle);
        }

        char kmb_path[512];
        snprintf(kmb_path, sizeof(kmb_path), "%s/%s.kmb", out_dir, sp->name);
        pqc_nist_report_t rep;
        uint64_t states = 0u;
        double sec = 0.0;
        if (!built || !pqc_compile_and_verify(sp->name, sp->family, sp->search_bits, &m, &ar, &cfg, &inst,
                                              sp->phase_oracle, kmb_path, 2u, &rep, &states, &sec)) {
            fails++;
        } else {
            char sb[24], md64t[24], md64q[24];
            pqc_format_pow2_q16(rep.log2_serial_t_q16, sb, sizeof(sb));
            pqc_format_pow2_q16(rep.maxdepth[1].log2_total_t_q16, md64t, sizeof(md64t));
            pqc_format_pow2_q16(rep.maxdepth[1].log2_total_qubits_q16, md64q, sizeof(md64q));
            printf("  %-24s | %6u | %8llu | %9llu | %7llu | %8s | %8s | %8s | %4.2fs\n",
                   rep.name, (unsigned)rep.oracle_stats.num_qubits,
                   (unsigned long long)rep.oracle_stats.toffoli,
                   (unsigned long long)rep.t_eq_gidney,
                   (unsigned long long)rep.oracle_stats.toffoli_depth,
                   sb, md64t, md64q, sec);
            reps[nr++] = rep;
        }
        qir_module_free(&m); qcc_arena_free(&ar);
    }

    if (json_path && nr > 0u) pqc_write_json_report(json_path, reps, nr);
    if (md_path && nr > 0u)   pqc_write_md_report(md_path, reps, nr);
    return fails == 0 ? 0 : 1;
}

static int cmd_analyze(int argc, char **argv) {
    const char *in_path = NULL;
    const char *label   = "custom_oracle";
    uint32_t bits = 0u;
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: pqc analyze <circuit.kmb> --bits <k> [--name <label>]\n");
            return 0;
        } else if (strcmp(argv[i], "--bits") == 0 && i + 1 < argc) {
            bits = (uint32_t)strtoul(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--name") == 0 && i + 1 < argc) {
            label = argv[++i];
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "Error: unknown or incomplete option '%s'\n", argv[i]);
            return 2;
        } else if (!in_path) {
            in_path = argv[i];
        }
    }
    if (!in_path || bits == 0u) {
        fprintf(stderr, "Usage: pqc analyze <circuit.kmb> --bits <k> [--name <label>]\n");
        return 2;
    }

    qcc_kmb_circuit_t kmb;
    qcc_kmb_fiber_t rf = qcc_kmb_read(&kmb, in_path);
    if (rf != QCC_KMB_FIBER_COMMITTED) {
        fprintf(stderr, "Error: failed to read %s (%s, code=%d)\n",
                in_path, qcc_kmb_fiber_name(rf), (int)rf);
        return 1;
    }
    qcc_kmb_stats_t ks;
    qcc_kmb_recompute_stats(&kmb, &ks);
    pqc_nist_report_t rep;
    pqc_nist_estimate(label, "Custom", bits, &ks, &rep);

    char sb[24], sd[24];
    pqc_format_pow2_q16(rep.log2_serial_t_q16, sb, sizeof(sb));
    pqc_format_pow2_q16(rep.log2_serial_depth_q16, sd, sizeof(sd));
    printf("[pqc analyze] %s (%s) | k=%u | Q=%u CCX=%llu HMR=%llu CNOT=%llu\n",
           rep.name, in_path, (unsigned)rep.search_bits,
           (unsigned)ks.num_qubits, (unsigned long long)ks.toffoli,
           (unsigned long long)ks.hmr, (unsigned long long)ks.cnot);
    printf("  Oracle cost: T-eq(Gidney)=%llu, T-eq(Unitary)=%llu | T-depth=%llu, Circ-depth=%llu, D_react=%llu\n",
           (unsigned long long)rep.t_eq_gidney, (unsigned long long)rep.t_eq_unitary,
           (unsigned long long)ks.toffoli_depth, (unsigned long long)ks.circuit_depth,
           (unsigned long long)ks.reaction_depth);
    printf("  Serial Grover: Total_T=%s, Total_Depth=%s | NIST Equiv: Cat %u\n",
           sb, sd, (unsigned)rep.nist_category_equiv);
    for (uint32_t i = 0u; i < 3u; ++i) {
        char p[24], q[24], tt[24], gd[24];
        pqc_format_pow2_q16(rep.maxdepth[i].log2_instances_q16, p, sizeof(p));
        pqc_format_pow2_q16(rep.maxdepth[i].log2_total_qubits_q16, q, sizeof(q));
        pqc_format_pow2_q16(rep.maxdepth[i].log2_total_t_q16, tt, sizeof(tt));
        pqc_format_pow2_q16(rep.maxdepth[i].log2_g_times_d_q16, gd, sizeof(gd));
        printf("  MAXDEPTH=2^%u: P=%s Total_Q=%s Total_T=%s G*D=%s\n",
               (unsigned)rep.maxdepth[i].maxdepth_log2, p, q, tt, gd);
    }
    qcc_kmb_free(&kmb);
    return 0;
}

static void print_usage(FILE *out) {
    fprintf(out,
            "PQC — NIST Post-Quantum Cryptanalysis Resource Compiler (QCC v3)\n"
            "Subcommands:\n"
            "  pqc selftest\n"
            "  pqc bench [out_dir] [--json <path>] [--md <path>]\n"
            "  pqc analyze <circuit.kmb> --bits <k> [--name <label>]\n");
}

int main(int argc, char **argv) {
    if (argc < 2) {
        print_usage(stderr);
        return 2;
    }
    if (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "help") == 0) {
        print_usage(stdout);
        return 0;
    }
    if (strcmp(argv[1], "selftest") == 0) return cmd_selftest(argc, argv);
    if (strcmp(argv[1], "bench") == 0)    return cmd_bench(argc, argv);
    if (strcmp(argv[1], "analyze") == 0)  return cmd_analyze(argc, argv);
    fprintf(stderr, "Unknown subcommand: %s\n", argv[1]);
    print_usage(stderr);
    return 2;
}
