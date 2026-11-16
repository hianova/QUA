/*
 * qua_synth.c — Thin Benchmark Driver for QCC (Quantum Circuit Compiler)
 * Pure C17, zero raw gate emission, zero bespoke optimization heuristics.
 */
#include "qir_pass_kmb.h"

static bool compile_and_log(
    const char *tag,
    qir_module_t *m,
    qcc_arena_t *ar,
    const qcc_opt_cfg_t *cfg,
    const char *path,
    qcc_kmb_stats_t *out_st
) {
    qcc_kmb_stats_t st;
    if (!qcc_compile_module_to_file(m, ar, cfg, path, &st)) {
        fprintf(stderr, "[ERROR] QCC compilation failed for %s -> %s\n", tag, path);
        return false;
    }
    if (out_st) *out_st = st;
    printf("%s -> %s | Q=%u CCX=%llu Z_POW=%llu T_LIKE=%llu HMR=%llu\n",
           tag, path, st.num_qubits,
           (unsigned long long)st.toffoli,
           (unsigned long long)st.z_pow,
           (unsigned long long)st.t_like,
           (unsigned long long)st.hmr);
    return true;
}

static void mul_f2_mats_12(const uint16_t a_rows[12], const uint16_t b_rows[12], qcc_bv256_t out[12]) {
    for (uint32_t r = 0u; r < 12u; ++r) {
        uint16_t acc = 0u;
        for (uint32_t k = 0u; k < 12u; ++k) {
            if ((a_rows[r] >> k) & 1u) acc ^= b_rows[k];
        }
        memset(&out[r], 0, sizeof(out[r]));
        out[r].w[0] = (uint64_t)acc;
    }
}

static int run_selftests(const char *out_dir) {
    char path[512];
    qcc_opt_cfg_t base_cfg = qcc_opt_cfg_default();
    base_cfg.require_clean_scratch = true;
    base_cfg.emit_clean_reg        = true;
    qcc_kmb_stats_t st;

    /* 1. selftest_adder16.kmb */
    {
        qcc_arena_t ar; qcc_arena_init(&ar, 512u * 1024u);
        qir_module_t m; qir_module_init(&m, &ar, "selftest_adder16");
        m.require_clean_scratch = true;
        uint32_t lhs = qir_add_input(&m, 16u, "lhs");
        uint32_t tgt = qir_add_input(&m, 16u, "target");
        uint32_t out = qir_emit_add(&m, lhs, tgt, 0u, "target");
        qir_add_output(&m, out);
        snprintf(path, sizeof(path), "%s/selftest_adder16.kmb", out_dir);
        bool ok = compile_and_log("[1/8]", &m, &ar, &base_cfg, path, &st);
        qir_module_free(&m); qcc_arena_free(&ar);
        if (!ok || st.num_qubits != 47u || st.toffoli != 15u) return 1;
    }

    /* 2. selftest_gf2_12_mul.kmb */
    {
        qcc_arena_t ar; qcc_arena_init(&ar, 512u * 1024u);
        qir_module_t m; qir_module_init(&m, &ar, "selftest_gf2_12_mul");
        m.require_clean_scratch = true;
        uint32_t a   = qir_add_input(&m, 12u, "a");
        uint32_t b   = qir_add_input(&m, 12u, "b");
        uint32_t out = qir_add_input(&m, 12u, "out");
        uint32_t p   = qir_emit_gf2m_mul(&m, a, b, 12u, 0x1053u, "p");
        uint32_t o2  = qir_emit_gf2m_add(&m, out, p, 12u, "out");
        qir_add_output(&m, o2);
        snprintf(path, sizeof(path), "%s/selftest_gf2_12_mul.kmb", out_dir);
        bool ok = compile_and_log("[2/8]", &m, &ar, &base_cfg, path, &st);
        qir_module_free(&m); qcc_arena_free(&ar);
        if (!ok || st.num_qubits != 37u || st.toffoli != 42u) return 1;
    }

    /* 3. selftest_gf64_mul.kmb */
    {
        qcc_arena_t ar; qcc_arena_init(&ar, 512u * 1024u);
        qir_module_t m; qir_module_init(&m, &ar, "selftest_gf64_mul");
        m.require_clean_scratch = true;
        uint32_t a   = qir_add_input(&m, 6u, "a");
        uint32_t b   = qir_add_input(&m, 6u, "b");
        uint32_t out = qir_add_input(&m, 6u, "out");
        uint32_t p   = qir_emit_gf2m_mul(&m, a, b, 6u, 0u, "p");
        uint32_t o2  = qir_emit_gf2m_add(&m, out, p, 6u, "out");
        qir_add_output(&m, o2);
        snprintf(path, sizeof(path), "%s/selftest_gf64_mul.kmb", out_dir);
        bool ok = compile_and_log("[3/8]", &m, &ar, &base_cfg, path, &st);
        qir_module_free(&m); qcc_arena_free(&ar);
        if (!ok || st.num_qubits != 19u || st.toffoli != 15u) return 1;
    }

    /* 4. selftest_qroam_vslot0.kmb */
    {
        uint64_t tab37[37];
        for (uint32_t r = 0u; r < 37u; ++r) tab37[r] = (uint64_t)((r * 0x5BDu + 0x137u) & 0x7FFu);
        qcc_arena_t ar; qcc_arena_init(&ar, 512u * 1024u);
        qir_module_t m; qir_module_init(&m, &ar, "selftest_qroam_vslot0");
        m.require_clean_scratch = true;
        uint32_t sel = qir_add_input(&m, 6u, "selection");
        uint32_t tgt = qir_add_input(&m, 12u, "target");
        uint32_t lv  = qir_emit_lut_rom(&m, sel, 12u, 37u, tab37, 0u, "v");
        uint32_t o2  = qir_emit_gf2m_add(&m, tgt, lv, 12u, "target");
        qir_add_output(&m, o2);
        snprintf(path, sizeof(path), "%s/selftest_qroam_vslot0.kmb", out_dir);
        bool ok = compile_and_log("[4/8]", &m, &ar, &base_cfg, path, &st);
        qir_module_free(&m); qcc_arena_free(&ar);
        if (!ok) return 1;
    }

    /* 5. selftest_zpow128.kmb */
    {
        qir_phase_term_t terms[4];
        memset(terms, 0, sizeof(terms));
        for (uint32_t i = 0u; i < 4u; ++i) qcc_bv_set(terms[i].mask, i);
        terms[0].angle = ((flow_u128)UINT64_C(0x2000000000000000) << 64u);
        terms[1].angle = ((flow_u128)UINT64_C(0x6000000000000000) << 64u);
        terms[2].angle = ((flow_u128)UINT64_C(0x4000000000000000) << 64u);
        terms[3].angle = ((flow_u128)UINT64_C(0x123456789ABCDEF0) << 64u) | UINT64_C(0xFEDCBA9876543210);
        qcc_arena_t ar; qcc_arena_init(&ar, 512u * 1024u);
        qir_module_t m; qir_module_init(&m, &ar, "selftest_zpow128");
        uint32_t x = qir_add_input(&m, 4u, "x");
        (void)qir_emit_phase_poly(&m, &x, 1u, terms, 4u, "ph");
        qir_add_output(&m, x);
        qcc_opt_cfg_t zcfg = base_cfg;
        zcfg.emit_clean_reg = false;
        snprintf(path, sizeof(path), "%s/selftest_zpow128.kmb", out_dir);
        bool ok = compile_and_log("[5/8]", &m, &ar, &zcfg, path, &st);
        qir_module_free(&m); qcc_arena_free(&ar);
        if (!ok || st.num_qubits != 4u || st.z_pow != 4u || st.t_like != 2u) return 1;
    }

    /* 6. selftest_hmr_lut.kmb */
    {
        uint64_t tab37[37];
        for (uint32_t r = 0u; r < 37u; ++r) tab37[r] = (uint64_t)((r * 0x5BDu + 0x137u) & 0xFFFu);
        qcc_arena_t ar; qcc_arena_init(&ar, 512u * 1024u);
        qir_module_t m; qir_module_init(&m, &ar, "selftest_hmr_lut");
        m.require_clean_scratch = true;
        uint32_t sel = qir_add_input(&m, 6u, "selection");
        uint32_t out = qir_add_input(&m, 12u, "out");
        uint32_t lv  = qir_emit_lut_rom(&m, sel, 12u, 37u, tab37, 0u, "v");
        uint32_t o2  = qir_emit_gf2m_add(&m, out, lv, 12u, "out");
        qir_add_output(&m, o2);
        qcc_opt_cfg_t hcfg = base_cfg;
        hcfg.enabled_rules &= ~((1u << 8u) | (1u << 7u) | (1u << 12u));
        snprintf(path, sizeof(path), "%s/selftest_hmr_lut.kmb", out_dir);
        bool ok = compile_and_log("[6/8]", &m, &ar, &hcfg, path, &st);
        qir_module_free(&m); qcc_arena_free(&ar);
        if (!ok || st.hmr == 0u) return 1;
    }

    /* 7. selftest_linf2.kmb */
    {
        uint16_t l1[12], u1[12], l2[12], u2[12];
        for (uint32_t r = 0u; r < 12u; ++r) {
            uint16_t lo_m = (uint16_t)((1u << r) - 1u);
            uint16_t hi_m = (uint16_t)(0xFFFu ^ ((1u << (r + 1u)) - 1u));
            l1[r] = (uint16_t)((1u << r) | (((r + 1u) * 0x5ADu) & lo_m));
            u1[r] = (uint16_t)((1u << r) | (((r + 3u) * 0x36Bu) & hi_m));
            l2[r] = (uint16_t)((1u << r) | (((r + 5u) * 0x6D3u) & lo_m));
            u2[r] = (uint16_t)((1u << r) | (((r + 2u) * 0x59Du) & hi_m));
        }
        qcc_bv256_t m1_rows[12], m2_rows[12], zero_c;
        memset(&zero_c, 0, sizeof(zero_c));
        mul_f2_mats_12(l1, u1, m1_rows);
        mul_f2_mats_12(u2, l2, m2_rows);

        qcc_arena_t ar; qcc_arena_init(&ar, 512u * 1024u);
        qir_module_t m; qir_module_init(&m, &ar, "selftest_linf2");
        m.require_clean_scratch = true;
        uint32_t x = qir_add_input(&m, 12u, "x");
        uint32_t t = qir_emit_linf2(&m, &x, 1u, 12u, m1_rows, &zero_c, "t");
        uint32_t y = qir_emit_linf2(&m, &t, 1u, 12u, m2_rows, &zero_c, "x");
        qir_add_output(&m, y);
        snprintf(path, sizeof(path), "%s/selftest_linf2.kmb", out_dir);
        bool ok = compile_and_log("[7/8]", &m, &ar, &base_cfg, path, &st);
        qir_module_free(&m); qcc_arena_free(&ar);
        if (!ok || st.num_qubits != 12u || st.toffoli != 0u) return 1;
    }

    /* 8. selftest_cuccaro.kmb */
    {
        qcc_arena_t ar; qcc_arena_init(&ar, 512u * 1024u);
        qir_module_t m; qir_module_init(&m, &ar, "selftest_cuccaro");
        m.require_clean_scratch = true;
        uint32_t ctrl = qir_add_input(&m, 1u, "ctrl");
        uint32_t lhs  = qir_add_input(&m, 16u, "lhs");
        uint32_t tgt  = qir_add_input(&m, 16u, "target");
        uint32_t sum  = qir_emit_add(&m, tgt, lhs, 0u, "s");
        uint32_t out  = qir_emit_ctrl_sel(&m, ctrl, tgt, sum, 0u, "target");
        qir_add_output(&m, out);
        snprintf(path, sizeof(path), "%s/selftest_cuccaro.kmb", out_dir);
        bool ok = compile_and_log("[8/8]", &m, &ar, &base_cfg, path, &st);
        qir_module_free(&m); qcc_arena_free(&ar);
        if (!ok || st.num_qubits != 34u || st.toffoli != 46u) return 1;
    }
    return 0;
}

static int run_invsqrt(const char *out_dir) {
    char path[512];
    snprintf(path, sizeof(path), "%s/invsqrt.kmb", out_dir);
    qcc_arena_t ar;
    if (!qcc_arena_init(&ar, 16u * 1024u * 1024u)) return 1;
    qir_module_t m;
    qir_module_init(&m, &ar, "invsqrt");
    m.require_clean_scratch = true;
    uint32_t x = qir_add_input(&m, 22u, "x");
    if (!qir_emit_dyadic_phase_fn(&m, x, 1u, UINT64_C(3133452), 8u, 37u, 28u, qir__builtin_invsqrt_ref, NULL)) {
        qir_module_free(&m); qcc_arena_free(&ar);
        return 1;
    }
    qir_add_output(&m, x);

    qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
    cfg.require_clean_scratch = true;
    cfg.emit_clean_reg        = false;
    cfg.cost_mode             = QCC_COST_QC;
    bool ok = compile_and_log("[invsqrt]", &m, &ar, &cfg, path, NULL);
    qir_module_free(&m);
    qcc_arena_free(&ar);
    return ok ? 0 : 1;
}

static int run_dfthcsos(const char *out_dir) {
    char words_path[512], thc_dir[512];
    snprintf(thc_dir, sizeof(thc_dir), "%s/dfthcsos", out_dir);
    snprintf(words_path, sizeof(words_path), "%s/words.bin", thc_dir);
    FILE *f = fopen(words_path, "rb");
    if (!f) f = fopen("dfthcsos_opt/words.bin", "rb");
    if (!f) {
        fprintf(stderr, "[ERROR] words.bin not found at %s or dfthcsos_opt/words.bin\n", words_path);
        return 1;
    }
    int32_t hdr[4] = {0};
    if (fread(hdr, sizeof(int32_t), 4, f) != 4) { fclose(f); return 1; }
    uint32_t R = (uint32_t)hdr[0], B = (uint32_t)hdr[1], brot = (uint32_t)hdr[3];
    uint32_t num_rows = 76u + R * B;
    uint32_t addr_bits = 1u;
    while ((UINT32_C(1) << addr_bits) < num_rows) ++addr_bits;

    const uint32_t num_angles = 75u;
    size_t n_words = (size_t)num_rows * num_angles;
    uint16_t *words = (uint16_t *)malloc(n_words * sizeof(uint16_t));
    if (!words || fread(words, sizeof(uint16_t), n_words, f) != n_words) {
        free(words); fclose(f); return 1;
    }
    fclose(f);

    for (uint32_t k = 0u; k < num_angles; ++k) {
        char stale[512];
        snprintf(stale, sizeof(stale), "%s/angle_k_%u.kmb", thc_dir, k);
        remove(stale);
    }

    const uint32_t J = 5u;
    uint32_t K = (num_angles + J - 1u) / J;
    uint32_t w_bits = J * brot;
    uint16_t stride_u64 = (uint16_t)((w_bits + 63u) / 64u);
    if (stride_u64 == 0u) stride_u64 = 1u;
    uint64_t *tab = (uint64_t *)calloc((size_t)num_rows * stride_u64, sizeof(uint64_t));
    if (!tab) { free(words); return 1; }

    qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
    cfg.max_qubits            = 128u;
    cfg.cost_mode             = QCC_COST_QT;
    cfg.require_clean_scratch = true;
    cfg.emit_clean_reg        = false;

    uint32_t max_q = 0u;
    uint64_t tot_ccx = 0u;
    for (uint32_t k = 0u; k < K; ++k) {
        memset(tab, 0, (size_t)num_rows * stride_u64 * sizeof(uint64_t));
        for (uint32_t r = 0u; r < num_rows; ++r) {
            for (uint32_t j = 0u; j < J; ++j) {
                uint32_t ic = k * J + j;
                uint16_t wc = (ic < num_angles) ? words[(size_t)r * num_angles + ic] : 0u;
                uint16_t wp = (k > 0u && (k - 1u) * J + j < num_angles)
                    ? words[(size_t)r * num_angles + (k - 1u) * J + j] : 0u;
                uint16_t val = (k == 0u) ? wc : (uint16_t)(wc ^ wp);
                for (uint32_t b = 0u; b < brot; ++b) {
                    if ((val >> b) & 1u) {
                        uint32_t bit = j * brot + b;
                        tab[(size_t)r * stride_u64 + (bit >> 6u)] |= (UINT64_C(1) << (bit & 63u));
                    }
                }
            }
        }

        qcc_arena_t ar;
        if (!qcc_arena_init(&ar, 4u * 1024u * 1024u)) { free(tab); free(words); return 1; }
        qir_module_t m;
        qir_module_init(&m, &ar, "dfthcsos_batch");
        m.require_clean_scratch = true;
        uint32_t sel = qir_add_input(&m, (uint16_t)addr_bits, "selection");
        uint32_t w_in[16];
        for (uint32_t j = 0u; j < J; ++j) {
            char wn[32];
            snprintf(wn, sizeof(wn), "word_%u", j);
            w_in[j] = qir_add_input(&m, (uint16_t)brot, wn);
        }
        uint32_t y_bus  = qir_emit_concat(&m, w_in, (uint16_t)J, "y_bus");
        uint32_t lut    = qir_emit_lut_rom(&m, sel, (uint16_t)w_bits, num_rows, tab, 0u, "lut");
        uint32_t y_next = qir_emit_gf2m_add(&m, y_bus, lut, (uint16_t)w_bits, "y_next");
        for (uint32_t j = 0u; j < J; ++j) {
            char wn[32];
            snprintf(wn, sizeof(wn), "word_%u", j);
            uint32_t wo = qir_emit_slice(&m, y_next, (uint16_t)(j * brot), (uint16_t)((j + 1u) * brot), wn);
            qir_add_output(&m, wo);
        }

        char kmb_path[512];
        snprintf(kmb_path, sizeof(kmb_path), "%s/angle_k_%u.kmb", thc_dir, k);
        qcc_kmb_stats_t st;
        bool ok = qcc_compile_module_to_file(&m, &ar, &cfg, kmb_path, &st);
        qir_module_free(&m);
        qcc_arena_free(&ar);
        if (!ok || st.num_qubits > 128u || st.z_pow != 0u) {
            fprintf(stderr, "[ERROR] dfthcsos batch %u failed (ok=%d Q=%u ZP=%llu)\n",
                    k, (int)ok, st.num_qubits, (unsigned long long)st.z_pow);
            free(tab); free(words);
            return 1;
        }
        if (st.num_qubits > max_q) max_q = st.num_qubits;
        tot_ccx += st.toffoli;
    }

    printf("[dfthcsos] Emitted %u circuits (J=%u, brot=%u, rows=%u) into %s | max_Q=%u<=128, 1pass_CCX=%llu\n",
           K, J, brot, num_rows, thc_dir, max_q, (unsigned long long)tot_ccx);
    free(tab);
    free(words);
    return 0;
}

int main(int argc, char **argv) {
    const char *out_dir = (argc > 1) ? argv[1] : "build";
    const char *mode    = (argc > 2) ? argv[2] : "all";
    bool do_all      = (strcmp(mode, "all") == 0);
    bool do_selftest = do_all || (strcmp(mode, "selftest") == 0);
    bool do_invsqrt  = do_all || (strcmp(mode, "invsqrt") == 0);
    bool do_dfthcsos = do_all || (strcmp(mode, "dfthcsos") == 0);
    if (!do_selftest && !do_invsqrt && !do_dfthcsos) {
        fprintf(stderr, "Usage: %s [out_dir] [all|selftest|invsqrt|dfthcsos]\n", argv[0]);
        return 1;
    }
    if (do_selftest && run_selftests(out_dir) != 0) return 1;
    if (do_invsqrt  && run_invsqrt(out_dir) != 0) return 1;
    if (do_dfthcsos && run_dfthcsos(out_dir) != 0) return 1;
    return 0;
}
