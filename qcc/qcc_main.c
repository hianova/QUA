/*
 * qcc_main.c — Universal Fault-Tolerant Quantum Circuit Compiler CLI & Selftest Suite
 * Part of QCC (Quantum Circuit Compiler) for FLOW v7.0.1.
 * Pure C17, zero external dependencies.
 */
#include "qir_verify.h"
#include "qir_pass_kmb.h"

static bool qcc__ends_with(const char *s, const char *suf) {
    size_t ls = strlen(s), lsf = strlen(suf);
    return (ls >= lsf) && (strcmp(s + (ls - lsf), suf) == 0);
}

static qcc_cost_mode_t qcc__parse_metric(const char *s) {
    if (!s) return QCC_COST_DEFAULT;
    if (strcmp(s, "qc") == 0)  return QCC_COST_QC;
    if (strcmp(s, "q2t") == 0) return QCC_COST_Q2T;
    if (strcmp(s, "qt") == 0)  return QCC_COST_QT;
    return QCC_COST_DEFAULT;
}

static int cmd_lift(int argc, char **argv) {
    const char *in_path = NULL, *out_path = NULL;
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) out_path = argv[++i];
        else if (!in_path) in_path = argv[i];
    }
    if (!in_path || !out_path) {
        fprintf(stderr, "Usage: qcc lift <in.kmb> -o <out.qir>\n");
        return 2;
    }
    qcc_kmb_circuit_t kmb;
    if (qcc_kmb_read(&kmb, in_path) != 0) {
        fprintf(stderr, "Error: failed to read %s\n", in_path);
        return 1;
    }
    qcc_arena_t arena;
    qcc_arena_init(&arena, 512u * 1024u);
    qir_module_t mod;
    qir_lift_stats_t lst;
    if (!qir_lift_from_kmb(&mod, &arena, &kmb, "lifted", &lst)) {
        fprintf(stderr, "Error: failed to lift %s\n", in_path);
        qcc_kmb_free(&kmb);
        qcc_arena_free(&arena);
        return 1;
    }
    if (!qir_write_file(&mod, out_path)) {
        fprintf(stderr, "Error: failed to write %s\n", out_path);
        qir_module_free(&mod);
        qcc_kmb_free(&kmb);
        qcc_arena_free(&arena);
        return 1;
    }
    printf("[qcc lift] %s -> %s (ops=%zu, linf2=%u, phase=%u, and_vec=%u, opaque=%u)\n",
           in_path, out_path, kmb.num_ops,
           (unsigned)lst.lifted_linf2, (unsigned)lst.lifted_phase,
           (unsigned)lst.lifted_and_vec, (unsigned)lst.lifted_opaque);
    qir_module_free(&mod);
    qcc_kmb_free(&kmb);
    qcc_arena_free(&arena);
    return 0;
}

static int cmd_opt(int argc, char **argv) {
    const char *in_path = NULL, *out_path = NULL, *emit_qir = NULL;
    qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
    uint32_t verify_batches = 64u;
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) out_path = argv[++i];
        else if (strcmp(argv[i], "--max-q") == 0 && i + 1 < argc) cfg.max_qubits = (uint32_t)strtoul(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--metric") == 0 && i + 1 < argc) cfg.cost_mode = qcc__parse_metric(argv[++i]);
        else if (strcmp(argv[i], "--emit-qir") == 0 && i + 1 < argc) emit_qir = argv[++i];
        else if (strcmp(argv[i], "--iters") == 0 && i + 1 < argc) cfg.searu_iters = (uint32_t)strtoul(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--verify-batches") == 0 && i + 1 < argc) verify_batches = (uint32_t)strtoul(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--clean-anc") == 0) cfg.require_clean_scratch = true;
        else if (!in_path) in_path = argv[i];
    }
    if (!in_path || !out_path) {
        fprintf(stderr, "Usage: qcc opt <in.qir|in.kmb> -o <out.kmb> [--max-q N] [--metric default|qc|q2t|qt] [--emit-qir post.qir] [--clean-anc]\n");
        return 2;
    }

    qcc_arena_t arena;
    qcc_arena_init(&arena, 1024u * 1024u);
    qir_module_t mod;
    bool ok = false;

    if (qcc__ends_with(in_path, ".kmb")) {
        qcc_kmb_circuit_t kmb;
        if (qcc_kmb_read(&kmb, in_path) == 0) {
            ok = qir_lift_from_kmb(&mod, &arena, &kmb, "lifted_opt", NULL);
            qcc_kmb_free(&kmb);
        }
    } else {
        ok = qir_parse_file(&mod, &arena, in_path);
    }
    if (!ok) {
        fprintf(stderr, "Error: failed to load %s\n", in_path);
        qcc_arena_free(&arena);
        return 1;
    }

    if (emit_qir) qir_write_file(&mod, emit_qir);

    qir_egraph_t eg;
    qir_phys_plan_t plan;
    if (!qir_optimize_module(&mod, &arena, &cfg, &eg, &plan)) {
        fprintf(stderr, "Error: optimization failed\n");
        qir_module_free(&mod);
        qcc_arena_free(&arena);
        return 1;
    }

    qcc_kmb_circuit_t out_kmb;
    if (!qir_codegen_emit(&plan, &out_kmb)) {
        fprintf(stderr, "Error: codegen invariant check failed\n");
        qir_phys_plan_free(&plan);
        qir_egraph_free(&eg);
        qir_module_free(&mod);
        qcc_arena_free(&arena);
        return 1;
    }

    qir_kmb_pass_stats_t pstats;
    qir_pass_kmb_optimize_with_plan(&out_kmb, &plan, &pstats);

    qcc_kmb_stats_t final_st;
    qcc_kmb_recompute_stats(&out_kmb, &final_st);

    if (qcc_kmb_write(&out_kmb, out_path) != 0) {
        fprintf(stderr, "Error: failed to write %s\n", out_path);
        qcc_kmb_free(&out_kmb);
        qir_phys_plan_free(&plan);
        qir_egraph_free(&eg);
        qir_module_free(&mod);
        qcc_arena_free(&arena);
        return 1;
    }

    bool any_approx = false;
    uint64_t approx_max = 0u;
    for (uint32_t oi = 0u; oi < mod.num_ops; ++oi) {
        if (mod.ops[oi].flags & QIR_FLAG_APPROX) {
            any_approx = true;
            uint64_t tol = (mod.ops[oi].approx_tol_q16 >> 16u) + 4u;
            if (tol > approx_max) approx_max = tol;
        }
    }

    qir_verify_report_t vrep;
    bool v_ok = qir_verify_kmb_vs_qir(&out_kmb, &mod, &plan, verify_batches, cfg.seed, any_approx, approx_max, &vrep);

    printf("[qcc opt] %s -> %s | Q=%u (cap=%u) CCX=%llu Z_POW=%llu T_LIKE=%llu HMR=%llu D_react=%llu | pass_kmb(-%llu z,-%llu inv) | disjoint=%d clean_anc=%d | verify=%s (%llu states, max_err=%llu)\n",
           in_path, out_path,
           (unsigned)final_st.num_qubits, (unsigned)cfg.max_qubits,
           (unsigned long long)final_st.toffoli,
           (unsigned long long)final_st.z_pow,
           (unsigned long long)final_st.t_like,
           (unsigned long long)final_st.hmr,
           (unsigned long long)final_st.reaction_depth,
           (unsigned long long)pstats.z_pow_fused_or_cancelled,
           (unsigned long long)pstats.self_inv_cancelled,
           plan.disjointness_certified ? 1 : 0,
           plan.clean_scratch_certified ? 1 : 0,
           v_ok ? "PASS" : "FAIL",
           (unsigned long long)vrep.states_tested,
           (unsigned long long)vrep.max_abs_err);

    qcc_kmb_free(&out_kmb);
    qir_phys_plan_free(&plan);
    qir_egraph_free(&eg);
    qir_module_free(&mod);
    qcc_arena_free(&arena);
    return v_ok ? 0 : 1;
}

static int cmd_verify(int argc, char **argv) {
    const char *cand_path = NULL, *ref_path = NULL;
    uint32_t num_states = 4096u;
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--ref") == 0 && i + 1 < argc) ref_path = argv[++i];
        else if (strcmp(argv[i], "--states") == 0 && i + 1 < argc) num_states = (uint32_t)strtoul(argv[++i], NULL, 10);
        else if (!cand_path) cand_path = argv[i];
    }
    if (!cand_path) {
        fprintf(stderr, "Usage: qcc verify <circuit.kmb> [--ref <ref.kmb|ref.qir>] [--states N]\n");
        return 2;
    }

    qcc_kmb_circuit_t cand;
    if (qcc_kmb_read(&cand, cand_path) != 0) {
        fprintf(stderr, "Error: failed to read %s\n", cand_path);
        return 1;
    }
    qcc_kmb_stats_t cst;
    qcc_kmb_recompute_stats(&cand, &cst);
    uint32_t batches = (num_states + 63u) / 64u;
    if (batches == 0u) batches = 1u;

    if (!ref_path) {
        /* Self-consistency dual-seed check on declared input/IO registers (keeping @clean/work ancillas |0>) */
        uint32_t io_q = 0u;
        if (cand.num_regs > 0u) {
            for (uint32_t r = 0u; r < cand.num_regs; ++r) {
                if (cand.regs[r].name[0] == '@') continue;
                for (uint32_t j = 0u; j < cand.regs[r].len; ++j) {
                    uint32_t tgt = cand.regs[r].targets[j];
                    if (tgt & UINT32_C(0x80000000)) {
                        uint32_t q = tgt & UINT32_C(0x7FFFFFFF);
                        if (q + 1u > io_q) io_q = q + 1u;
                    }
                }
            }
        }
        if (io_q == 0u) io_q = cst.num_qubits;
        qir_verify_report_t rep;
        bool ok = qir_verify_kmb_vs_kmb(&cand, &cand, io_q, batches, UINT64_C(0xCAFEBABE), &rep);
        printf("[qcc verify] %s | Q=%u B=%u CCX=%llu Z_POW=%llu HMR=%llu D_react=%llu | dual-seed=%s (%llu states)\n",
               cand_path, (unsigned)cst.num_qubits, (unsigned)cst.num_bits,
               (unsigned long long)cst.toffoli, (unsigned long long)cst.z_pow,
               (unsigned long long)cst.hmr, (unsigned long long)cst.reaction_depth,
               ok ? "PASS" : "FAIL", (unsigned long long)rep.states_tested);
        qcc_kmb_free(&cand);
        return ok ? 0 : 1;
    }

    if (qcc__ends_with(ref_path, ".kmb")) {
        qcc_kmb_circuit_t ref;
        if (qcc_kmb_read(&ref, ref_path) != 0) {
            fprintf(stderr, "Error: failed to read ref %s\n", ref_path);
            qcc_kmb_free(&cand);
            return 1;
        }
        qcc_kmb_stats_t rst;
        qcc_kmb_recompute_stats(&ref, &rst);
        uint32_t io_q = rst.num_qubits < cst.num_qubits ? rst.num_qubits : cst.num_qubits;
        if (cand.num_regs > 0u && cand.regs[0].len > 0u) io_q = cand.regs[0].len;
        qir_verify_report_t rep;
        bool ok = qir_verify_kmb_vs_kmb(&cand, &ref, io_q, batches, UINT64_C(0x1122334455), &rep);
        printf("[qcc verify] %s vs %s | %s (%llu states, val_err=%llu, seed_drift=%llu, phase_err=%llu)\n",
               cand_path, ref_path, ok ? "PASS" : "FAIL",
               (unsigned long long)rep.states_tested,
               (unsigned long long)rep.value_mismatches,
               (unsigned long long)rep.seed_drift_mismatches,
               (unsigned long long)rep.phase_mismatches);
        qcc_kmb_free(&ref);
        qcc_kmb_free(&cand);
        return ok ? 0 : 1;
    } else {
        qcc_arena_t arena;
        qcc_arena_init(&arena, 512u * 1024u);
        qir_module_t ref_mod;
        if (!qir_parse_file(&ref_mod, &arena, ref_path)) {
            fprintf(stderr, "Error: failed to parse %s\n", ref_path);
            qcc_kmb_free(&cand);
            qcc_arena_free(&arena);
            return 1;
        }
        qir_phys_plan_t dummy_plan;
        memset(&dummy_plan, 0, sizeof(dummy_plan));
        dummy_plan.mod = &ref_mod;
        dummy_plan.val_allocs = (qir_val_alloc_t *)calloc(ref_mod.num_vals, sizeof(qir_val_alloc_t));
        uint8_t reg_used[256] = {0};
        uint32_t qoff = 0u;
        for (uint32_t i = 0u; i < ref_mod.num_inputs; ++i) {
            uint32_t vid = ref_mod.inputs[i];
            uint16_t w = ref_mod.vals[vid].width;
            dummy_plan.val_allocs[vid].width  = w;
            dummy_plan.val_allocs[vid].qubits = w ? (uint32_t *)qcc_arena_alloc(&arena, (size_t)w * sizeof(uint32_t), 4u) : NULL;
            const char *vname = ref_mod.vals[vid].name;
            int32_t hit = -1;
            for (uint32_t r = 0u; r < cand.num_regs && r < 256u; ++r) {
                if (reg_used[r] || cand.regs[r].len == 0u || ((cand.regs[r].targets[0] & UINT32_C(0x80000000)) == 0u)) continue;
                if (vname && strcmp(cand.regs[r].name, vname) == 0) { hit = (int32_t)r; break; }
            }
            if (hit < 0) {
                for (uint32_t r = 0u; r < cand.num_regs && r < 256u; ++r) {
                    if (!reg_used[r] && cand.regs[r].len >= w && ((cand.regs[r].targets[0] & UINT32_C(0x80000000)) != 0u)) { hit = (int32_t)r; break; }
                }
            }
            if (hit >= 0) {
                reg_used[hit] = 1u;
                for (uint16_t b = 0u; b < w; ++b) {
                    dummy_plan.val_allocs[vid].qubits[b] = (b < cand.regs[hit].len)
                        ? (cand.regs[hit].targets[b] & UINT32_C(0x7FFFFFFF)) : qoff++;
                }
            } else {
                for (uint16_t b = 0u; b < w; ++b) dummy_plan.val_allocs[vid].qubits[b] = qoff++;
            }
        }
        for (uint32_t i = 0u; i < ref_mod.num_outputs; ++i) {
            uint32_t vid = ref_mod.outputs[i];
            if (dummy_plan.val_allocs[vid].qubits) continue;
            uint16_t w = ref_mod.vals[vid].width;
            const char *vname = ref_mod.vals[vid].name;
            int32_t hit = -1;
            for (uint32_t r = 0u; r < cand.num_regs && r < 256u; ++r) {
                if (reg_used[r] || cand.regs[r].len == 0u || ((cand.regs[r].targets[0] & UINT32_C(0x80000000)) == 0u)) continue;
                if (vname && strcmp(cand.regs[r].name, vname) == 0) { hit = (int32_t)r; break; }
            }
            if (hit < 0) {
                for (uint32_t r = 0u; r < cand.num_regs && r < 256u; ++r) {
                    if (!reg_used[r] && cand.regs[r].len >= w && ((cand.regs[r].targets[0] & UINT32_C(0x80000000)) != 0u)) { hit = (int32_t)r; break; }
                }
            }
            if (hit >= 0) {
                reg_used[hit] = 1u;
                dummy_plan.val_allocs[vid].width  = w;
                dummy_plan.val_allocs[vid].qubits = w ? (uint32_t *)qcc_arena_alloc(&arena, (size_t)w * sizeof(uint32_t), 4u) : NULL;
                for (uint16_t b = 0u; b < w && b < cand.regs[hit].len; ++b) {
                    dummy_plan.val_allocs[vid].qubits[b] = cand.regs[hit].targets[b] & UINT32_C(0x7FFFFFFF);
                }
            }
        }
        for (uint32_t i = 0u; i < ref_mod.num_classical_outputs; ++i) {
            uint32_t vid = ref_mod.classical_outputs[i];
            if (vid >= ref_mod.num_vals || dummy_plan.val_allocs[vid].cbits) continue;
            uint16_t w = ref_mod.vals[vid].width;
            const char *vname = ref_mod.vals[vid].name;
            int32_t hit = -1;
            for (uint32_t r = 0u; r < cand.num_regs && r < 256u; ++r) {
                if (reg_used[r] || cand.regs[r].len == 0u || ((cand.regs[r].targets[0] & UINT32_C(0x80000000)) != 0u)) continue;
                if (vname && strcmp(cand.regs[r].name, vname) == 0) { hit = (int32_t)r; break; }
            }
            if (hit < 0) {
                for (uint32_t r = 0u; r < cand.num_regs && r < 256u; ++r) {
                    if (!reg_used[r] && cand.regs[r].len >= w && ((cand.regs[r].targets[0] & UINT32_C(0x80000000)) == 0u)) { hit = (int32_t)r; break; }
                }
            }
            if (hit >= 0) {
                reg_used[hit] = 1u;
                dummy_plan.val_allocs[vid].cwidth = w;
                dummy_plan.val_allocs[vid].cbits  = w ? (uint32_t *)qcc_arena_alloc(&arena, (size_t)w * sizeof(uint32_t), 4u) : NULL;
                for (uint16_t b = 0u; b < w && b < cand.regs[hit].len; ++b) {
                    dummy_plan.val_allocs[vid].cbits[b] = cand.regs[hit].targets[b] & UINT32_C(0x7FFFFFFF);
                }
            }
        }
        bool any_approx = false;
        uint64_t approx_max = 0u;
        for (uint32_t oi = 0u; oi < ref_mod.num_ops; ++oi) {
            if (ref_mod.ops[oi].flags & QIR_FLAG_APPROX) {
                any_approx = true;
                uint64_t tol = (ref_mod.ops[oi].approx_tol_q16 >> 16u) + 4u;
                if (tol > approx_max) approx_max = tol;
            }
        }
        dummy_plan.clean_scratch_certified = ref_mod.require_clean_scratch;
        qir_verify_report_t rep;
        bool ok = qir_verify_kmb_vs_qir(&cand, &ref_mod, &dummy_plan, batches, UINT64_C(0xCAFE1234), any_approx, approx_max, &rep);
        printf("[qcc verify] %s vs %s | %s (%llu states, val_err=%llu, seed_drift=%llu, phase_err=%llu, dirty_anc=%llu)\n",
               cand_path, ref_path, ok ? "PASS" : "FAIL",
               (unsigned long long)rep.states_tested,
               (unsigned long long)rep.value_mismatches,
               (unsigned long long)rep.seed_drift_mismatches,
               (unsigned long long)rep.phase_mismatches,
               (unsigned long long)rep.dirty_ancilla_mismatches);
        free(dummy_plan.val_allocs);
        qir_module_free(&ref_mod);
        qcc_kmb_free(&cand);
        qcc_arena_free(&arena);
        return ok ? 0 : 1;
    }
}

/* ============================================================================
 * Automated Selftest Suite (Rules R1..R7, Zero-Discrepancy, Constraints, KMB)
 * ============================================================================ */

static bool test_compile_and_verify(
    const char *label,
    qir_module_t *mod,
    qcc_arena_t *arena,
    const qcc_opt_cfg_t *cfg,
    uint32_t batches,
    bool expect_max_q,
    bool expect_zero_ccx,
    bool allow_approx,
    uint64_t max_approx_err,
    qcc_kmb_stats_t *out_final_st
) {
    qir_egraph_t eg;
    qir_phys_plan_t plan;
    if (!qir_optimize_module(mod, arena, cfg, &eg, &plan)) {
        printf("FAIL [%s]: qir_optimize_module failed\n", label);
        return false;
    }
    if (!plan.disjointness_certified) {
        printf("FAIL [%s]: disjointness certificate failed\n", label);
        qir_phys_plan_free(&plan); qir_egraph_free(&eg);
        return false;
    }
    if ((cfg->require_clean_scratch || mod->require_clean_scratch) && !plan.clean_scratch_certified) {
        printf("FAIL [%s]: clean_scratch_certified failed\n", label);
        qir_phys_plan_free(&plan); qir_egraph_free(&eg);
        return false;
    }
    qcc_kmb_circuit_t kmb;
    if (!qir_codegen_emit(&plan, &kmb)) {
        printf("FAIL [%s]: codegen zero-discrepancy assertion failed\n", label);
        qcc_kmb_free(&kmb); qir_phys_plan_free(&plan); qir_egraph_free(&eg);
        return false;
    }
    qir_kmb_pass_stats_t ps;
    qir_pass_kmb_optimize_with_plan(&kmb, &plan, &ps);

    qcc_kmb_stats_t st;
    qcc_kmb_recompute_stats(&kmb, &st);
    if (out_final_st) *out_final_st = st;

    if (expect_max_q && cfg->max_qubits > 0u && st.num_qubits > cfg->max_qubits) {
        printf("FAIL [%s]: Q=%u exceeds max_qubits=%u\n", label, (unsigned)st.num_qubits, (unsigned)cfg->max_qubits);
        qcc_kmb_free(&kmb); qir_phys_plan_free(&plan); qir_egraph_free(&eg);
        return false;
    }
    if (expect_zero_ccx && st.toffoli != 0u) {
        printf("FAIL [%s]: expected 0 CCX, got %llu\n", label, (unsigned long long)st.toffoli);
        qcc_kmb_free(&kmb); qir_phys_plan_free(&plan); qir_egraph_free(&eg);
        return false;
    }

    qir_verify_report_t vrep;
    bool ok = qir_verify_kmb_vs_qir(&kmb, mod, &plan, batches, UINT64_C(0xABCD1234), allow_approx, max_approx_err, &vrep);
    printf("  %-34s | Q=%3u CCX=%4llu ZP=%3llu HMR=%3llu | %s (%llu states",
           label, (unsigned)st.num_qubits, (unsigned long long)st.toffoli,
           (unsigned long long)st.z_pow, (unsigned long long)st.hmr,
           ok ? "PASS" : "FAIL", (unsigned long long)vrep.states_tested);
    if (!ok) {
        printf(", val_err=%llu, seed_drift=%llu, phase_err=%llu, dirty_anc=%llu, max_abs=%llu",
               (unsigned long long)vrep.value_mismatches,
               (unsigned long long)vrep.seed_drift_mismatches,
               (unsigned long long)vrep.phase_mismatches,
               (unsigned long long)vrep.dirty_ancilla_mismatches,
               (unsigned long long)vrep.max_abs_err);
    }
    printf(")\n");

    qcc_kmb_free(&kmb);
    qir_phys_plan_free(&plan);
    qir_egraph_free(&eg);
    return ok;
}

static int cmd_selftest(void) {
    printf("=== QCC Universal Compiler Selftest Suite ===\n");
    int fails = 0;

    /* T1: .kmb binary round-trip + lift -> codegen equivalence */
    {
        qcc_kmb_circuit_t k1, k2;
        qcc_kmb_init(&k1);
        uint32_t qa[4] = {0, 1, 2, 3}, qb[1] = {0};
        qcc_kmb_add_quantum_reg(&k1, "q", qa, 4u);
        qcc_kmb_add_classical_reg(&k1, "c", qb, 1u);
        qcc_kmb_x(&k1, 0u);
        qcc_kmb_cx(&k1, 0u, 1u);
        qcc_kmb_info_tag(&k1, "@enter:test");
        qcc_kmb_and_compute(&k1, 0u, 1u, 2u);
        qcc_kmb_z_pow_f64(&k1, 2u, 0.125);
        qcc_kmb_and_uncompute(&k1, 0u, 1u, 2u, 0u);
        qcc_kmb_info_tag(&k1, "@exit:test");

        const char *tmp = "/tmp/qcc_selftest_rt.kmb";
        if (qcc_kmb_write(&k1, tmp) != 0 || qcc_kmb_read(&k2, tmp) != 0 ||
            k1.num_ops != k2.num_ops || k1.num_info_tags != k2.num_info_tags ||
            memcmp(k1.op_types, k2.op_types, k1.num_ops) != 0) {
            printf("  FAIL [T1 .kmb round-trip]\n");
            fails++;
        } else {
            qcc_arena_t ar; qcc_arena_init(&ar, 256u * 1024u);
            qir_module_t m;
            qir_lift_stats_t ls;
            qir_lift_from_kmb(&m, &ar, &k2, "rt_lift", &ls);
            if (ls.lifted_opaque != 0u) {
                printf("  FAIL [T1 expected 0 opaque islands, got %u]\n", ls.lifted_opaque);
                fails++;
            }
            qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
            qcc_kmb_stats_t st;
            if (!test_compile_and_verify("T1_KMB_Roundtrip_Lift", &m, &ar, &cfg, 128u, false, false, false, 0u, &st)) fails++;
            qir_module_free(&m);
            qcc_arena_free(&ar);
        }
        qcc_kmb_free(&k1);
        qcc_kmb_free(&k2);
        remove(tmp);
    }

    /* T2: Rule R1 (LinF2_Fuse) + DSL parse round-trip */
    {
        qcc_arena_t ar; qcc_arena_init(&ar, 256u * 1024u);
        qir_module_t m, m_parsed;
        qir_module_init(&m, &ar, "r1_linf2");
        uint32_t x = qir_add_input(&m, 8u, "x");
        qcc_bv256_t r1[8], r2[8], c1 = qcc_bv_zero(), c2 = qcc_bv_zero();
        for (uint32_t i = 0u; i < 8u; ++i) {
            r1[i] = qcc_bv_zero(); qcc_bv_set(r1[i].w, i); qcc_bv_set(r1[i].w, (i + 1u) % 8u);
            r2[i] = qcc_bv_zero(); qcc_bv_set(r2[i].w, i); qcc_bv_set(r2[i].w, (i + 3u) % 8u);
        }
        qcc_bv_set(c1.w, 2u); qcc_bv_set(c2.w, 5u);
        uint32_t t1 = qir_emit_linf2(&m, &x, 1u, 8u, r1, &c1, "t1");
        uint32_t t2 = qir_emit_linf2(&m, &t1, 1u, 8u, r2, &c2, "t2");
        qir_add_output(&m, t2);
        qir_recompute_zero_masks(&m);

        const char *dsl_tmp = "/tmp/qcc_selftest_dsl.qir";
        qir_write_file(&m, dsl_tmp);
        if (!qir_parse_file(&m_parsed, &ar, dsl_tmp)) {
            printf("  FAIL [T2 DSL round-trip parse]\n");
            fails++;
        } else {
            qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
            qcc_kmb_stats_t st;
            if (!test_compile_and_verify("T2_R1_LinF2_Fuse_DSL", &m_parsed, &ar, &cfg, 1024u, false, true, false, 0u, &st)) fails++;
            qir_module_free(&m_parsed);
        }
        remove(dsl_tmp);
        qir_module_free(&m);
        qcc_arena_free(&ar);
    }

    /* T3: Rule R2 (HMR_Dual) LUT Compute + Eager Uncompute */
    {
        qcc_arena_t ar; qcc_arena_init(&ar, 256u * 1024u);
        qir_module_t m;
        qir_module_init(&m, &ar, "r2_hmr_lut");
        uint32_t idx = qir_add_input(&m, 6u, "idx");
        uint32_t y_in = qir_add_input(&m, 8u, "y");
        uint64_t tbl1[64], tbl2[64];
        for (uint32_t i = 0u; i < 64u; ++i) {
            tbl1[i] = ((i * i * 3u + 7u) & 0xFFu);
            tbl2[i] = ((i * 5u + 11u) & 0xFFu);
        }
        uint32_t lut1 = qir_emit_lut_rom(&m, idx, 8u, 64u, tbl1, 0u, "lut1");
        uint32_t sum1 = qir_emit_add(&m, lut1, y_in, 0u, "sum1");
        uint32_t idx2 = qir_emit_slice(&m, sum1, 0u, 6u, "idx2");
        uint32_t lut2 = qir_emit_lut_rom(&m, idx2, 8u, 64u, tbl2, 0u, "lut2");
        uint32_t out  = qir_emit_add(&m, lut2, sum1, 0u, "out");
        qir_add_output(&m, out);
        qir_recompute_zero_masks(&m);

        qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
        cfg.max_qubits = 46u;
        qcc_kmb_stats_t st;
        if (!test_compile_and_verify("T3_R2_HMR_ParityLUT_Uncomp", &m, &ar, &cfg, 512u, true, false, false, 0u, &st)) fails++;
        qir_module_free(&m);
        qcc_arena_free(&ar);
    }

    /* T4: Rule R3 (CSA_Tree) Multi-Operand Add Tree */
    {
        qcc_arena_t ar; qcc_arena_init(&ar, 256u * 1024u);
        qir_module_t m;
        qir_module_init(&m, &ar, "r3_csa_tree");
        uint32_t a = qir_add_input(&m, 8u, "a");
        uint32_t b = qir_add_input(&m, 8u, "b");
        uint32_t c = qir_add_input(&m, 8u, "c");
        uint32_t s1 = qir_emit_add(&m, a, b, 0u, "s1");
        uint32_t s2 = qir_emit_add(&m, s1, c, 0u, "s2");
        qir_add_output(&m, s2);
        qir_recompute_zero_masks(&m);

        qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
        qcc_kmb_stats_t st;
        if (!test_compile_and_verify("T4_R3_CSA_Tree_Add3", &m, &ar, &cfg, 512u, false, false, false, 0u, &st)) fails++;
        qir_module_free(&m);
        qcc_arena_free(&ar);
    }

    /* T5: Rule R4 (Phase_CSA_Split) Exact 0-CCX Elimination */
    {
        qcc_arena_t ar; qcc_arena_init(&ar, 256u * 1024u);
        qir_module_t m;
        qir_module_init(&m, &ar, "r4_phase_csa");
        uint32_t a = qir_add_input(&m, 8u, "a");
        uint32_t b = qir_add_input(&m, 8u, "b");
        uint32_t sum = qir_emit_add(&m, a, b, 0u, "sum");
        qir_phase_term_t terms[8];
        memset(terms, 0, sizeof(terms));
        for (uint32_t j = 0u; j < 8u; ++j) {
            qcc_bv_set(terms[j].mask, j);
            terms[j].angle = ((flow_u128)1u << (120u + j)); /* k=1/256 * 2^j */
        }
        (void)qir_emit_phase_poly(&m, &sum, 1u, terms, 8u, "ph");
        qir_add_output(&m, a);
        qir_recompute_zero_masks(&m);

        qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
        qcc_kmb_stats_t st;
        if (!test_compile_and_verify("T5_R4_Phase_CSA_Split_0CCX", &m, &ar, &cfg, 512u, false, true, false, 0u, &st)) fails++;
        qir_module_free(&m);
        qcc_arena_free(&ar);
    }

    /* T6: Rule R5 (LUT_ColSplit) + R7 (BasisIso_Gauge) Low-Rank Wide LUT */
    {
        qcc_arena_t ar; qcc_arena_init(&ar, 512u * 1024u);
        qir_module_t m;
        qir_module_init(&m, &ar, "r5_r7_wide_lut");
        uint32_t idx = qir_add_input(&m, 6u, "idx");
        uint64_t tbl[64];
        /* 24-bit LUT whose image lives inside a 6-dimensional affine subspace! */
        for (uint32_t i = 0u; i < 64u; ++i) {
            uint64_t v = ((uint64_t)i) | (((uint64_t)i ^ 0x2Au) << 8u) | (((uint64_t)i ^ 0x15u) << 16u);
            tbl[i] = v;
        }
        uint32_t val = qir_emit_lut_rom(&m, idx, 24u, 64u, tbl, 0u, "wide");
        qir_add_output(&m, val);
        qir_recompute_zero_masks(&m);

        qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
        cfg.max_qubits = 46u;
        qcc_kmb_stats_t st;
        if (!test_compile_and_verify("T6_R5_R7_WideLUT_Gauge", &m, &ar, &cfg, 512u, true, false, false, 0u, &st)) fails++;
        qir_module_free(&m);
        qcc_arena_free(&ar);
    }

    /* T7: Rule R6 (GF2m_Tower_Rank): GF(2^12) R42 + GF(64) R15 + GF(2^12) Tower Inv */
    {
        qcc_arena_t ar; qcc_arena_init(&ar, 512u * 1024u);
        qir_module_t m;
        qir_module_init(&m, &ar, "r6_gf12_tower");
        uint32_t a = qir_add_input(&m, 12u, "a");
        uint32_t b = qir_add_input(&m, 12u, "b");
        uint32_t inv_a = qir_emit_gf2m_inv(&m, a, 12u, 0x1053u, "inv_a");
        uint32_t prod  = qir_emit_gf2m_mul(&m, inv_a, b, 12u, 0x1053u, "prod");
        qir_add_output(&m, prod);
        qir_recompute_zero_masks(&m);

        qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
        qcc_kmb_stats_t st;
        if (!test_compile_and_verify("T7_R6_GF12_TowerInv_Mul42", &m, &ar, &cfg, 512u, false, false, false, 0u, &st)) fails++;
        qir_module_free(&m);
        qcc_arena_free(&ar);
    }

    /* T8: CLZ_NORM + SQR_FXP + MUL_FXP + FLOW 2-Jet Bounded Approx */
    {
        qcc_arena_t ar; qcc_arena_init(&ar, 512u * 1024u);
        qir_module_t m;
        qir_module_init(&m, &ar, "t8_clz_fxp_flow");
        uint32_t x = qir_add_input(&m, 8u, "x");
        uint32_t y = qir_add_input(&m, 6u, "y");
        uint32_t mant = 0u, clz = 0u;
        qir_emit_clz_norm(&m, x, 6u, &mant, &clz, "mant", "clz");
        uint32_t sq = qir_emit_sqr_fxp(&m, mant, 7u, 5u, 5u, QIR_FLAG_APPROX, "sq");
        m.ops[m.num_ops - 1u].approx_tol_q16 = (2u << 16u);
        uint32_t p  = qir_emit_mul_fxp(&m, sq, y, 8u, 5u, 4u, 4u, QIR_FLAG_APPROX, "p");
        m.ops[m.num_ops - 1u].approx_tol_q16 = (4u << 16u);
        qir_add_output(&m, p);
        qir_add_output(&m, clz);
        qir_recompute_zero_masks(&m);

        qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
        qcc_kmb_stats_t st;
        if (!test_compile_and_verify("T8_CLZ_SQR_MUL_FLOW2Jet", &m, &ar, &cfg, 256u, false, false, true, 16u, &st)) fails++;
        qir_module_free(&m);
        qcc_arena_free(&ar);
    }

    /* T9: Post-Codegen Peephole Phase Fusion & Self-Inverse Cancellation */
    {
        qcc_kmb_circuit_t k;
        qcc_kmb_init(&k);
        uint32_t qs[2] = {0, 1};
        qcc_kmb_add_quantum_reg(&k, "q", qs, 2u);
        qcc_kmb_z_pow_f64(&k, 0u, 0.125);
        qcc_kmb_cx(&k, 0u, 1u);
        qcc_kmb_z_pow_f64(&k, 0u, 0.125); /* commutes past ctrl q[0] -> fuses to 0.25 (S)! */
        qcc_kmb_cx(&k, 0u, 1u);           /* CX;CX -> cancels! */
        qcc_kmb_x(&k, 1u);
        qcc_kmb_x(&k, 1u);                /* X;X -> cancels! */

        qir_kmb_pass_stats_t ps;
        qir_pass_kmb_optimize(&k, &ps);
        qcc_kmb_stats_t st;
        qcc_kmb_recompute_stats(&k, &st);
        bool ok = (st.z_pow == 1u && st.cnot == 0u && k.num_ops == 1u);
        printf("  %-34s | ops: 6->%zu ZP=%llu CNOT=%llu | %s\n",
               "T9_PassKMB_PhaseFuse_Cancel", k.num_ops,
               (unsigned long long)st.z_pow, (unsigned long long)st.cnot,
               ok ? "PASS" : "FAIL");
        if (!ok) fails++;
        qcc_kmb_free(&k);
    }

    /* T10: Rule R8 (QROAM_VSlot0) + LUT zero_mask column stripping */
    {
        qcc_arena_t ar; qcc_arena_init(&ar, 512u * 1024u);
        qir_module_t m;
        qir_module_init(&m, &ar, "t10_r8_qroam");
        uint32_t idx = qir_add_input(&m, 6u, "idx");
        uint64_t tbl[64];
        for (uint32_t i = 0u; i < 64u; ++i) {
            /* Top 4 MSB bits [12..15] strictly 0 -> stripped cleanly by zero_mask + R8 */
            tbl[i] = ((i * i * 7u + 13u * i + 3u) & UINT64_C(0x0FFF));
        }
        uint32_t val = qir_emit_lut_rom(&m, idx, 16u, 64u, tbl, 0u, "val");
        qir_add_output(&m, val);
        qir_recompute_zero_masks(&m);

        bool zm_ok = ((m.vals[val].zero_mask[0] & UINT64_C(0xF000)) == UINT64_C(0xF000));
        qcc_opt_cfg_t cfg_std = qcc_opt_cfg_default();
        cfg_std.enabled_rules = 0u;
        qir_egraph_t eg0; qir_phys_plan_t p0;
        qir_optimize_module(&m, &ar, &cfg_std, &eg0, &p0);
        uint64_t std_ccx = p0.est_toffoli;
        qir_phys_plan_free(&p0); qir_egraph_free(&eg0);

        qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
        qir_egraph_t eg1; qir_phys_plan_t p1;
        qir_optimize_module(&m, &ar, &cfg, &eg1, &p1);
        bool used_r8 = false;
        for (uint32_t si = 0u; si < p1.num_steps; ++si) {
            if (p1.steps[si].variant == QIR_VAR_QROAM_VSLOT0) { used_r8 = true; break; }
        }
        qir_phys_plan_free(&p1); qir_egraph_free(&eg1);

        qcc_kmb_stats_t st;
        bool v_ok = test_compile_and_verify("T10_R8_QROAM_VSlot0_ZeroMSB", &m, &ar, &cfg, 512u, false, false, false, 0u, &st);
        if (!v_ok || !zm_ok || !used_r8 || st.toffoli >= std_ccx) {
            printf("  FAIL [T10 R8 check: zm=%d r8=%d ccx=%llu vs std=%llu]\n",
                   zm_ok ? 1 : 0, used_r8 ? 1 : 0, (unsigned long long)st.toffoli, (unsigned long long)std_ccx);
            fails++;
        }
        qir_module_free(&m);
        qcc_arena_free(&ar);
    }

    /* T11: Rule R9 (MUL_CONST_ADD & SQR_TRIANGLE) Algebraic Strength Reduction */
    {
        qcc_arena_t ar; qcc_arena_init(&ar, 512u * 1024u);
        qir_module_t m;
        qir_module_init(&m, &ar, "t11_r9_constmul_trisqr");
        uint32_t x  = qir_add_input(&m, 6u, "x");
        uint32_t sq = qir_emit_sqr_fxp(&m, x, 6u, 5u, 5u, 0u, "sq");
        uint32_t c  = qir_emit_const_u64(&m, 6u, 19u, "c");
        uint32_t p  = qir_emit_mul_fxp(&m, sq, c, 6u, 5u, 5u, 5u, 0u, "p");
        qir_add_output(&m, p);
        qir_recompute_zero_masks(&m);

        qcc_opt_cfg_t cfg_std = qcc_opt_cfg_default();
        cfg_std.enabled_rules = 0u;
        qir_egraph_t eg0; qir_phys_plan_t p0;
        qir_optimize_module(&m, &ar, &cfg_std, &eg0, &p0);
        uint32_t std_q = p0.est_qubits;
        uint64_t std_ccx = p0.est_toffoli;
        qir_phys_plan_free(&p0); qir_egraph_free(&eg0);

        qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
        qir_egraph_t eg1; qir_phys_plan_t p1;
        qir_optimize_module(&m, &ar, &cfg, &eg1, &p1);
        bool used_tri = false, used_cmul = false;
        for (uint32_t si = 0u; si < p1.num_steps; ++si) {
            if (p1.steps[si].variant == QIR_VAR_SQR_TRIANGLE)  used_tri  = true;
            if (p1.steps[si].variant == QIR_VAR_MUL_CONST_ADD) used_cmul = true;
        }
        qir_phys_plan_free(&p1); qir_egraph_free(&eg1);

        qcc_kmb_stats_t st;
        bool v_ok = test_compile_and_verify("T11_R9_ConstMul_TriSqr", &m, &ar, &cfg, 512u, false, false, false, 0u, &st);
        if (!v_ok || !used_tri || !used_cmul || st.toffoli >= std_ccx || st.num_qubits >= std_q) {
            printf("  FAIL [T11 R9 check: tri=%d cmul=%d Q=%u vs %u, CCX=%llu vs %llu]\n",
                   used_tri ? 1 : 0, used_cmul ? 1 : 0, (unsigned)st.num_qubits, (unsigned)std_q,
                   (unsigned long long)st.toffoli, (unsigned long long)std_ccx);
            fails++;
        }
        qir_module_free(&m);
        qcc_arena_free(&ar);
    }

    /* T12: Cascading Zero-Garbage Uncomputation Oracle (clean_anc + backing-set aliases) */
    {
        const char *dsl =
            "(qir \"t12_clean_oracle\"\n"
            "  (clean_anc)\n"
            "  (input %x u4)\n"
            "  (let %lut_v u6 (lut_rom %x (0x01 0x03 0x07 0x0c 0x12 0x19 0x1f 0x28 0x2d 0x33 0x37 0x3b 0x3d 0x3e 0x3f 0x00)))\n"
            "  (let %s4 u4 (slice %lut_v 1 5))\n"
            "  (let %sq u4 (sqr_fxp %s4 :frac 3 3))\n"
            "  (let %k u4 (const 5))\n"
            "  (let %prod u4 (mul_fxp %sq %k :frac 3 3 3))\n"
            "  (let %subv u4 (sub %prod %s4))\n"
            "  (let %ph u0 (phase_poly (%subv) (and 0x3 :turns 0.125) (and 0x5 :turns 0.25)))\n"
            "  (output %x)\n"
            ")\n";
        qcc_arena_t ar; qcc_arena_init(&ar, 512u * 1024u);
        qir_module_t m;
        if (!qir_parse_str(&m, &ar, dsl, strlen(dsl))) {
            printf("  FAIL [T12 DSL parse]\n");
            fails++;
        } else {
            qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
            qcc_kmb_stats_t st;
            if (!test_compile_and_verify("T12_Cascading_ZeroGarbage_Oracle", &m, &ar, &cfg, 256u, false, false, false, 0u, &st)) fails++;
            qir_module_free(&m);
        }
        qcc_arena_free(&ar);
    }

    /* T13: Scale-Free .kmb Lifter (nq = 320 > 256, SWAP, HMR+CZ_IF AND-uncompute) */
    {
        qcc_kmb_circuit_t k;
        qcc_kmb_init(&k);
        uint32_t qs[320];
        for (uint32_t i = 0u; i < 320u; ++i) qs[i] = i;
        qcc_kmb_add_quantum_reg(&k, "q", qs, 320u);
        uint32_t cbits[1] = {0u};
        qcc_kmb_add_classical_reg(&k, "c", cbits, 1u);

        qcc_kmb_x(&k, 0u);
        qcc_kmb_x(&k, 65u);
        qcc_kmb_x(&k, 260u);
        qcc_kmb_x(&k, 319u);
        qcc_kmb_cx(&k, 0u, 65u);
        qcc_kmb_cx(&k, 130u, 200u);
        qcc_kmb_cx(&k, 255u, 257u);
        qcc_kmb_cx(&k, 310u, 319u);
        qcc_kmb_swap(&k, 3u, 290u);
        qcc_kmb_and_compute(&k, 10u, 270u, 280u);
        qcc_kmb_z_pow_f64(&k, 280u, 0.125);
        qcc_kmb_and_uncompute(&k, 10u, 270u, 280u, 0u);

        qcc_arena_t ar; qcc_arena_init(&ar, 2048u * 1024u);
        qir_module_t m;
        qir_lift_stats_t ls;
        bool lifted = qir_lift_from_kmb(&m, &ar, &k, "lift320", &ls);
        if (!lifted || ls.lifted_opaque != 0u) {
            printf("  FAIL [T13 lift 320Q failed: lifted=%d opaque=%u]\n", lifted ? 1 : 0, ls.lifted_opaque);
            fails++;
        } else {
            qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
            qir_egraph_t eg; qir_phys_plan_t plan;
            qcc_kmb_circuit_t opt_k;
            bool ok = qir_optimize_module(&m, &ar, &cfg, &eg, &plan) && qir_codegen_emit(&plan, &opt_k);
            qir_kmb_pass_stats_t ps;
            if (ok) qir_pass_kmb_optimize_with_plan(&opt_k, &plan, &ps);

            qir_verify_report_t qrep;
            if (ok) ok = qir_verify_kmb_vs_qir(&opt_k, &m, &plan, 32u, UINT64_C(0x320320), false, 0u, &qrep);

            /* Also cross-check opt_k directly against original 320Q circuit k */
            uint32_t k_nq = 320u;
            uint32_t opt_nq = opt_k.num_qubits ? opt_k.num_qubits : 320u;
            uint64_t *init_orig = (uint64_t *)calloc(k_nq, sizeof(uint64_t));
            uint64_t *init_opt  = (uint64_t *)calloc(opt_nq, sizeof(uint64_t));
            uint64_t *out_orig  = (uint64_t *)calloc(k_nq, sizeof(uint64_t));
            uint64_t *out_opt   = (uint64_t *)calloc(opt_nq, sizeof(uint64_t));
            uint64_t rng = UINT64_C(0x320320AABBCCDD);
            uint64_t cross_err = 0u;
            flow_u128 goff = 0u;
            bool has_goff = false;

            for (uint32_t b = 0u; b < 32u && ok; ++b) {
                for (uint32_t q = 0u; q < 320u; ++q) {
                    init_orig[q] = (b == 0u && q == 0u) ? 0u : qcc_splitmix64(&rng);
                }
                init_orig[280] = 0u;
                memset(init_opt, 0, opt_nq * sizeof(uint64_t));
                uint32_t pos = 0u;
                for (uint32_t ii = 0u; ii < m.num_inputs; ++ii) {
                    uint32_t vid = m.inputs[ii];
                    for (uint16_t bit = 0u; bit < m.vals[vid].width && pos < 320u; ++bit, ++pos) {
                        uint32_t pq = plan.val_allocs[vid].qubits[bit];
                        if (pq < opt_nq) init_opt[pq] = init_orig[pos];
                    }
                }
                flow_u128 po[64], pk[64];
                if (!qir_sim_kmb_64(&k, init_orig, UINT64_C(0x1111) ^ b, out_orig, NULL, po) ||
                    !qir_sim_kmb_64(&opt_k, init_opt, UINT64_C(0x9999) ^ b, out_opt, NULL, pk)) {
                    cross_err++;
                    break;
                }
                for (uint32_t l = 0u; l < 64u; ++l) {
                    flow_u128 d = pk[l] - po[l];
                    if (!has_goff) { goff = d; has_goff = true; }
                    else if (d != goff) cross_err++;
                }
                pos = 0u;
                for (uint32_t oi = 0u; oi < m.num_outputs; ++oi) {
                    uint32_t vid = m.outputs[oi];
                    for (uint16_t bit = 0u; bit < m.vals[vid].width && pos < 320u; ++bit, ++pos) {
                        uint32_t pq = plan.val_allocs[vid].qubits[bit];
                        uint64_t got = (pq < opt_nq) ? out_opt[pq] : 0u;
                        if (got != out_orig[pos]) cross_err++;
                    }
                }
            }
            bool t13_ok = (ok && cross_err == 0u);
            printf("  %-34s | Q=%3u->%3u ops=%zu->%zu opaque=0 | %s (2048 states)\n",
                   "T13_ScaleFree_Lift_320Q", 320u, (unsigned)opt_k.num_qubits,
                   k.num_ops, opt_k.num_ops, t13_ok ? "PASS" : "FAIL");
            if (!t13_ok) fails++;

            free(init_orig); free(init_opt); free(out_orig); free(out_opt);
            if (ok) qcc_kmb_free(&opt_k);
            qir_phys_plan_free(&plan);
            qir_egraph_free(&eg);
            qir_module_free(&m);
        }
        qcc_arena_free(&ar);
        qcc_kmb_free(&k);
    }

    /* T14: Classical HMR Measurement, Kickback Cancellation & Classical Control Predicates */
    {
        const char *dsl =
            "(qir \"t14_classical\"\n"
            "  (input %a u1)\n"
            "  (input %b u1)\n"
            "  (input %x u2)\n"
            "  (let %r u1 (and_vec %a %b))\n"
            "  (let %m c1 (hmr_meas %r))\n"
            "  (let %ph u0 (phase_poly (%a %b) (and 0x3 :turns 0.5) :if %m))\n"
            "  (let %nm c1 (bit_not %m))\n"
            "  (let %t1 c1 (bit_xor %m %nm))\n"
            "  (let %c0 c1 (bit_and %m %nm))\n"
            "  (let %k u2 (const 1))\n"
            "  (let %out u2 (add %x %k :if %t1))\n"
            "  (let %out2 u2 (add %out %k :if %c0))\n"
            "  (output %a)\n"
            "  (output %b)\n"
            "  (output %out2)\n"
            "  (output_bits %t1)\n"
            "  (output_bits %c0)\n"
            ")\n";
        qcc_arena_t ar; qcc_arena_init(&ar, 512u * 1024u);
        qir_module_t m, m2;
        const char *tmp_qir = "/tmp/qcc_selftest_t14.qir";
        if (!qir_parse_str(&m, &ar, dsl, strlen(dsl)) || !qir_write_file(&m, tmp_qir) || !qir_parse_file(&m2, &ar, tmp_qir)) {
            printf("  FAIL [T14 DSL round-trip]\n");
            fails++;
        } else {
            qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
            qcc_kmb_stats_t st;
            if (!test_compile_and_verify("T14_Classical_Meas_Cond_Pred", &m2, &ar, &cfg, 256u, false, false, false, 0u, &st)) fails++;
            qir_module_free(&m2);
            qir_module_free(&m);
        }
        remove(tmp_qir);
        qcc_arena_free(&ar);
    }

    /* T15: Declarative approx_fn (invsqrt & recip) + FLOW 2-Jet Numerical Accuracy */
    {
        const char *dsl =
            "(qir \"t15_approx_fn\"\n"
            "  (input %x u8)\n"
            "  (let %y u8 (approx_fn %x :fn invsqrt :seg 4 :deg 2 :frac 5))\n"
            "  (let %r u8 (approx_fn %x :fn recip :seg 4 :deg 2 :frac 7))\n"
            "  (output %y)\n"
            "  (output %r)\n"
            ")\n";
        qcc_arena_t ar; qcc_arena_init(&ar, 2048u * 1024u);
        qir_module_t m;
        if (!qir_parse_str(&m, &ar, dsl, strlen(dsl))) {
            printf("  FAIL [T15 DSL parse]\n");
            fails++;
        } else {
            qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
            qcc_kmb_stats_t st;
            bool v_ok = test_compile_and_verify("T15_Declarative_ApproxFn_2Jet", &m, &ar, &cfg, 256u, false, false, true, 8u, &st);

            /* Numerically verify QIR approx_fn outputs against exact 1/sqrt(x) and 1/x for all x=1..255 */
            qir_val_vec64_t *in_v = (qir_val_vec64_t *)calloc(1u, sizeof(qir_val_vec64_t));
            qcc_arena_t ev_ar; qcc_arena_init(&ev_ar, 512u * 1024u);
            uint64_t max_y_err = 0u, max_r_err = 0u;
            for (uint32_t base = 0u; base < 256u; base += 64u) {
                memset(in_v, 0, sizeof(*in_v));
                for (uint32_t l = 0u; l < 64u; ++l) in_v->lanes[l][0] = base + l;
                qir_val_vec64_t *out_vals = NULL;
                flow_u128 ph[64];
                qcc_arena_reset(&ev_ar);
                if (!qir_eval_64(&m, in_v, UINT64_C(1), &ev_ar, &out_vals, ph)) {
                    max_y_err = 999u;
                    break;
                }
                for (uint32_t l = 0u; l < 64u; ++l) {
                    uint32_t x = base + l;
                    if (x == 0u) continue;
                    uint64_t got_y = out_vals[m.outputs[0]].lanes[l][0];
                    uint64_t got_r = out_vals[m.outputs[1]].lanes[l][0];
                    uint64_t ref_y = (flow_isqrt_round_u128(((flow_u128)1024u << 64u) / (flow_u128)x) + (UINT64_C(1) << 31u)) >> 32u;
                    uint64_t ref_r = (128u + x / 2u) / x;
                    uint64_t ey = (got_y > ref_y) ? (got_y - ref_y) : (ref_y - got_y);
                    uint64_t er = (got_r > ref_r) ? (got_r - ref_r) : (ref_r - got_r);
                    if (ey > max_y_err) max_y_err = ey;
                    if (er > max_r_err) max_r_err = er;
                }
            }
            if (!v_ok || max_y_err > 2u || max_r_err > 2u) {
                printf("  FAIL [T15 math accuracy: max_y_err=%llu max_r_err=%llu]\n",
                       (unsigned long long)max_y_err, (unsigned long long)max_r_err);
                fails++;
            }
            free(in_v);
            qcc_arena_free(&ev_ar);
            qir_module_free(&m);
        }
        qcc_arena_free(&ar);
    }

    printf("=== Selftest Summary: %s (fails=%d) ===\n", fails == 0 ? "ALL PASS" : "FAILED", fails);
    return fails == 0 ? 0 : 1;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "QCC — Universal Fault-Tolerant Quantum Circuit Compiler\n"
                        "Subcommands:\n"
                        "  qcc lift <in.kmb> -o <out.qir>\n"
                        "  qcc opt <in.qir|in.kmb> -o <out.kmb> [--max-q N] [--metric default|qc|q2t|qt] [--emit-qir post.qir] [--clean-anc]\n"
                        "  qcc verify <circuit.kmb> [--ref <ref.kmb|ref.qir>] [--states N]\n"
                        "  qcc selftest\n");
        return 2;
    }
    if (strcmp(argv[1], "lift") == 0)     return cmd_lift(argc, argv);
    if (strcmp(argv[1], "opt") == 0)      return cmd_opt(argc, argv);
    if (strcmp(argv[1], "verify") == 0)   return cmd_verify(argc, argv);
    if (strcmp(argv[1], "selftest") == 0) return cmd_selftest();
    fprintf(stderr, "Unknown subcommand: %s\n", argv[1]);
    return 2;
}
