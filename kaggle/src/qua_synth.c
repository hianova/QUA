/*
 * qua_synth.c — ELE-Style FLOW v7.0.1 x SEARU Reversible Quantum Physical-Structure Compiler
 * Competition: Google Quantum AI — Verifiable Quantum Circuit Compilation
 *
 * Bridges FLOW's 128-bit 0-ULP 2-Jet arithmetic (flow.h + flow_opt_params.h) and
 * SEARU's CLZ Zipf Sparse-k + Graph-Coloring + Bilevel solvers (searu.h) directly
 * to Kickmix v2 (.kmb) via qua_ele.h:
 *   [Pass 1] Reversible SSA Pebbling & Bilevel Cost-Law Co-Optimization (Q*C, Q^2*T, Q*T)
 *   [Pass 2] SEARU Live-Range Interference Graph Coloring & Disjoint Qubit Packing
 *   [Pass 3] Rank-42 Chudnovsky P^1(GF(2^3)) GF(2^12) Bilinear Multiplier (42 CCX, Q=37)
 *            + Gidney HMR Ripple-Carry Adder + Virtual-Slot-0 Select-Swap QROM Planner
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>

#include "flow_opt_params.h"
#include "flow.h"
#include "searu.h"
#include "qua_kmb.h"
#include "qua_ele.h"
#include "qua_invsqrt.h"
#include "qua_dqi.h"
#include "qua_dfthcsos.h"

/* ============================================================================
 * 1. Benchmark I (invsqrt): SEARU x 2-Jet Architecture Co-Optimizer (Min Q * C)
 * ============================================================================ */
typedef struct {
    uint32_t k_lut;
    uint32_t poly_deg;
    uint32_t b_swap;
    uint32_t newton_steps;
    uint32_t guard_bits;
    uint64_t peak_qubits;
    uint64_t toffolis;
    uint64_t cost_qc;
    uint32_t eff_prec_bits;
} qua_invsqrt_arch_t;

static void qua_eval_invsqrt_arch(uint64_t word0, uint32_t target_prec_bits, qua_invsqrt_arch_t *out) {
    uint32_t k_lut = 5u + (uint32_t)(word0 & 7u); /* 5..12 QROM mantissa bits */
    uint32_t d_raw = (uint32_t)((word0 >> 3u) & 3u);
    uint32_t deg   = (d_raw == 0u) ? 1u : ((d_raw == 1u) ? 2u : ((d_raw == 2u) ? 3u : 2u));
    uint32_t b_sh  = (uint32_t)((word0 >> 5u) & 7u);
    if (b_sh > k_lut) b_sh = k_lut;
    uint32_t b_sw  = 1u << b_sh;
    uint32_t m_raw = (uint32_t)((word0 >> 8u) & 3u);
    uint32_t m_new = (m_raw > 2u) ? 0u : m_raw;
    static const uint32_t GUARD_TAB[4] = {2u, 3u, 4u, 5u};
    uint32_t g_bits = GUARD_TAB[(word0 >> 10u) & 3u];

    uint32_t seed_bits = (deg + 1u) * (k_lut - 1u) + 2u;
    uint32_t eff_bits  = seed_bits << m_new;
    uint32_t w         = target_prec_bits + g_bits;

    /* Split-stage QROM (ELE Pass 2 live-range coloring: load c2 -> eval -> del_lookup -> load c1):
     * Only ONE coefficient register is live at a time alongside the accumulator! */
    uint64_t qrom_entries = UINT64_C(1) << k_lut;
    uint64_t qrom_sel_tof = (qrom_entries + b_sw - 1u) / b_sw;
    uint64_t qrom_swp_tof = (uint64_t)(b_sw - 1u) * (uint64_t)w;
    uint64_t qrom_tof     = (uint64_t)(deg + 1u) * (qrom_sel_tof + qrom_swp_tof);

    /* Reduced-width delta Taylor terms: wd1 ~ w - k_lut, wd2 ~ w - 2*k_lut */
    uint32_t wd1 = (w > k_lut) ? (w - k_lut) : 4u;
    uint64_t mul_tof    = (uint64_t)wd1 * (uint64_t)(wd1 + 1u) / 2ull;
    uint64_t horner_tof = (uint64_t)deg * mul_tof;
    uint64_t newton_tof = (uint64_t)m_new * (3ull * (uint64_t)w * (uint64_t)(w - 1u) / 2ull);

    uint64_t total_tof  = qrom_tof + horner_tof + newton_tof;
    uint64_t qrom_bus_q = (uint64_t)b_sw * (uint64_t)w;
    uint64_t arith_q    = (uint64_t)(2u + 2u * m_new) * (uint64_t)w + (uint64_t)k_lut;
    uint64_t peak_q     = 2ull * target_prec_bits + ((qrom_bus_q > arith_q) ? qrom_bus_q : arith_q);

    out->k_lut         = k_lut;
    out->poly_deg      = deg;
    out->b_swap        = b_sw;
    out->newton_steps  = m_new;
    out->guard_bits    = g_bits;
    out->peak_qubits   = peak_q;
    out->toffolis      = total_tof;
    out->eff_prec_bits = eff_bits;
    out->cost_qc       = peak_q * total_tof;
}

static uint64_t qua_eval_invsqrt_frust(uint64_t word0, uint32_t target_prec_bits,
                                       const qua_invsqrt_arch_t *plan) {
    if (plan->eff_prec_bits < target_prec_bits) {
        /* Precision deficit: steer k_lut [0..2], poly_deg [3..4], newton [8..9] */
        return UINT64_C(0x31F);
    }
    if (plan->newton_steps > 0u) {
        /* Newton iterations bloat both Q and Toffoli quadratically: steer away from newton + adjust k_lut/deg */
        return UINT64_C(0x31F);
    }
    /* Precision satisfied with Newton=0: steer k_lut [0..2], poly_deg [3..4], b_swap [5..7], guard [10..11] */
    (void)word0;
    return UINT64_C(0xCFF);
}

static qua_invsqrt_arch_t qua_optimize_invsqrt_plan_hp(uint32_t target_prec_bits,
                                                       const searu_hyperparams *hp,
                                                       uint64_t seed) {
    searu_result res;
    uint64_t init_c[4] = {0u, 0u, 0u, 0u};
    SEARU_SOLVE_INLINE_HP(&res, 12u, init_c, NULL, seed, 4096u, hp, cand, frust, energy, {
        qua_invsqrt_arch_t plan;
        qua_eval_invsqrt_arch(cand[0], target_prec_bits, &plan);
        if (plan.eff_prec_bits < target_prec_bits) {
            uint64_t deficit = (uint64_t)(target_prec_bits - plan.eff_prec_bits);
            energy = UINT64_C(1000000000000) + deficit * UINT64_C(10000000000) + plan.cost_qc;
        } else {
            energy = plan.cost_qc;
        }
        frust[0] = qua_eval_invsqrt_frust(cand[0], target_prec_bits, &plan);
    });
    qua_invsqrt_arch_t best;
    qua_eval_invsqrt_arch(res.canvas[0], target_prec_bits, &best);
    return best;
}

typedef struct {
    const qua_ele_dag_t *dag;
    uint32_t target_prec_bits;
} qua_bilevel_ctx_t;

static uint64_t qua_bilevel_inner_eval(const searu_hyperparams *hp, uint64_t seed, void *ctx) {
    const qua_bilevel_ctx_t *bctx = (const qua_bilevel_ctx_t *)ctx;
    qua_ele_pebble_plan_t p = qua_ele_solve_pebble_dag_hp(bctx->dag, hp, seed, 512u);
    qua_ele_reg_alloc_t ra  = qua_ele_alloc_qubits_hp(bctx->dag, &p, hp, seed ^ UINT64_C(0xA5A5));
    qua_invsqrt_arch_t inv  = qua_optimize_invsqrt_plan_hp(bctx->target_prec_bits, hp, seed ^ UINT64_C(0x51A7));
    uint64_t inv_cost = (inv.eff_prec_bits < bctx->target_prec_bits)
                            ? UINT64_C(1000000000000)
                            : inv.cost_qc;
    uint64_t bank_cost = ra.verified_disjoint
                             ? ((uint64_t)ra.phys_qubits_used * p.total_toffolis)
                             : (UINT64_C(1) << 40u);
    return p.objective_cost + bank_cost + inv_cost;
}

/* ============================================================================
 * 2. Benchmark III (dfthcsos): Virtual-Slot-0 <= 128-Qubit Select-Swap Planner
 * ============================================================================ */
typedef struct {
    uint32_t addr_bits;
    uint32_t angle_bits;
    uint32_t best_b;
    uint32_t total_qubits;
    uint64_t est_toffolis;
} qua_dfthc_plan_t;

static qua_dfthc_plan_t qua_plan_dfthc_angle_qrom(uint32_t num_entries, uint32_t angle_bits) {
    uint32_t addr_bits = 1u;
    while ((1u << addr_bits) < num_entries) ++addr_bits;

    qua_dfthc_plan_t best = {addr_bits, angle_bits, 1u, 0u, UINT64_MAX};
    for (uint32_t b = 1u; b <= (1u << addr_bits); b <<= 1u) {
        uint32_t log_b = 0u;
        while ((1u << log_b) < b) ++log_b;
        uint32_t unary_anc = (addr_bits > log_b) ? (addr_bits - log_b) : 1u;
        /* Virtual-Slot-0 trick: slot 0 lives directly on the phase-kickback target register,
         * so clean swap bus only needs (b - 1) * w_p qubits instead of b * angle_bits! */
        uint32_t q_total = addr_bits + unary_anc + (b - 1u) * angle_bits;
        if (q_total > 128u) break;

        uint64_t n_blocks = ((uint64_t)num_entries + b - 1u) / b;
        /* Pebbled unary iteration + Fredkin routing + HMR phase erase (1x CSWAP, not 2x!) */
        uint64_t tof = (n_blocks > 0u ? (n_blocks - 1u) : 0u) + (uint64_t)(b - 1u) * angle_bits;
        if (tof < best.est_toffolis) {
            best.best_b       = b;
            best.total_qubits = q_total;
            best.est_toffolis = tof;
        }
    }
    return best;
}

/* ============================================================================
 * Self-Test & Emit Verified Physical-Structure .kmb Circuits for Kickmix Check
 * ============================================================================ */
int main(int argc, char **argv) {
    const char *out_dir = (argc > 1) ? argv[1] : "build";
    const char *mode    = (argc > 2) ? argv[2] : "all";
    bool do_all      = (strcmp(mode, "all") == 0);
    bool do_selftest = do_all || (strcmp(mode, "selftest") == 0);
    bool do_invsqrt  = do_all || (strcmp(mode, "invsqrt") == 0);
    bool do_dqi      = do_all || (strcmp(mode, "dqi") == 0);
    bool do_dfthcsos = do_all || (strcmp(mode, "dfthcsos") == 0);
    if (!do_selftest && !do_invsqrt && !do_dqi && !do_dfthcsos) {
        fprintf(stderr, "Usage: %s [out_dir] [all|selftest|invsqrt|dqi|dfthcsos]\n", argv[0]);
        return 1;
    }

    if (do_selftest) {
    /* 0. Verify ELE Pass 1 (SSA Pebbling) + Pass 2 (SEARU Live-Range Coloring) on a
     *    10-node Split-2-Jet + Half-EEA prototype DAG */
    qua_ele_dag_t dag;
    memset(&dag, 0, sizeof(dag));
    dag.num_nodes      = 8u;
    dag.base_io_qubits = 64u;
    dag.cost_law       = QUA_COST_QC;
    /* Node 0: LZC shift_x (25q, persistent until end) */
    dag.nodes[0] = (qua_ele_ssa_node_t){25u, 110u, 110u, 0u, 0u, 0x00u};
    /* Node 1: delta^2 square (9q, HMR-erasable via lookup/Gidney) */
    dag.nodes[1] = (qua_ele_ssa_node_t){ 9u,  36u,  36u, 1u, 0u, 0x01u};
    /* Node 2: QROM c2 (9q, HMR-erasable via del_lookup) */
    dag.nodes[2] = (qua_ele_ssa_node_t){ 9u,  44u,  32u, 1u, 0u, 0x01u};
    /* Node 3: quad_term = c2 * sq (14q, accumulated into acc) */
    dag.nodes[3] = (qua_ele_ssa_node_t){14u,  96u,  96u, 1u, 0u, 0x06u};
    /* Node 4: QROM c1 (17q, allocated AFTER c2 & sq are released -> reuses same qubits!) */
    dag.nodes[4] = (qua_ele_ssa_node_t){17u,  68u,  48u, 1u, 0u, 0x01u};
    /* Node 5: lin_term = c1 * delta (20q) */
    dag.nodes[5] = (qua_ele_ssa_node_t){20u, 180u, 180u, 1u, 0u, 0x11u};
    /* Node 6: QROM c0 + sum (26q) */
    dag.nodes[6] = (qua_ele_ssa_node_t){26u,  92u,  64u, 1u, 0u, 0x29u};
    /* Node 7: Output denormalization (persistent IO) */
    dag.nodes[7] = (qua_ele_ssa_node_t){ 0u, 140u,   0u, 0u, 1u, 0x41u};

    qua_bilevel_ctx_t bctx = {&dag, 24u};
    qua_ele_self_opt_stats_t sopt_stats;
    searu_hyperparams opt_hp =
        qua_ele_self_opt_hp(qua_bilevel_inner_eval, &bctx, UINT64_C(0xB11E2026), 64u, &sopt_stats);

    qua_ele_pebble_plan_t naive_plan, opt_plan;
    qua_ele_eval_pebble_plan(&dag, 0u, &naive_plan);
    opt_plan = qua_ele_solve_pebble_dag_hp(&dag, &opt_hp, UINT64_C(0xE1E2026), 2048u);
    qua_ele_reg_alloc_t reg_alloc = qua_ele_alloc_qubits_hp(&dag, &opt_plan, &opt_hp, UINT64_C(0xC010A));
    if (!reg_alloc.verified_disjoint) {
        fprintf(stderr, "[ERROR] ELE Pass 2 qubit live-range overlap detected!\n");
        return 1;
    }
    printf("[Bilevel Meta-SEARU x FLOW 2-Jet] train=%llu->%llu | val=%llu->%llu (accepted=%d, cool_q16=%llu)\n",
           (unsigned long long)sopt_stats.base_train_cost, (unsigned long long)sopt_stats.opt_train_cost,
           (unsigned long long)sopt_stats.base_val_cost, (unsigned long long)sopt_stats.opt_val_cost,
           (int)sopt_stats.accepted, (unsigned long long)opt_hp.inl_cool_q16);
    printf("[ELE Pass 1+2] SSA Pebbling + Hetero-Width SEARU Packing: Q %u -> %u (banked=%u, banks=%u, disjoint=%d), Q*C %llu -> %llu\n",
           naive_plan.peak_qubits, opt_plan.peak_qubits, reg_alloc.phys_qubits_used,
           reg_alloc.num_colors_used, (int)reg_alloc.verified_disjoint,
           (unsigned long long)naive_plan.objective_cost, (unsigned long long)opt_plan.objective_cost);

    /* 1. Verify FLOW v7.0.1 0-ULP isqrt + optimize invsqrt architectural Pareto point */
    flow_u128 test_x = ((flow_u128)1 << 64);
    uint64_t root_chk = flow_isqrt_round_u128(test_x);
    qua_invsqrt_arch_t inv_plan = qua_optimize_invsqrt_plan_hp(24u, &opt_hp, UINT64_C(0x51A7E2026));

    /* 2. Emit 16-bit Gidney HMR Ripple-Carry Adder (.kmb) for Kickmix Fuzz Test */
    char path_adder[512];
    snprintf(path_adder, sizeof(path_adder), "%s/selftest_adder16.kmb", out_dir);
    {
        qua_builder_t kb;
        qua_builder_init(&kb);
        uint32_t lhs[16], tgt[16], anc[15];
        for (uint32_t i = 0; i < 16u; ++i) { lhs[i] = i; tgt[i] = 16u + i; }
        for (uint32_t i = 0; i < 15u; ++i) { anc[i] = 32u + i; }
        qua_add_quantum_reg(&kb, "lhs", lhs, 16u);
        qua_add_quantum_reg(&kb, "target", tgt, 16u);
        qua_add_quantum_reg(&kb, "@clean", anc, 15u);
        qua_emit_gidney_iadd(&kb, lhs, tgt, anc, 16u, 0u);
        if (qua_builder_write_kmb(&kb, path_adder) != 0) {
            fprintf(stderr, "Failed to write %s\n", path_adder);
            return 1;
        }
        printf("[1/8] Emitted %s | Q=%u, CCX=%llu, HMR=%llu (Gidney 16-bit iadd: exact n-1=15 Toffolis)\n",
               path_adder, kb.num_qubits, (unsigned long long)kb.count_toffoli, (unsigned long long)kb.count_hmr);
        qua_builder_free(&kb);
    }

    /* 3. Emit Rank-42 Chudnovsky P^1(GF(2^3)) GF(2^12) Multiplier (.kmb) */
    char path_gf12[512];
    snprintf(path_gf12, sizeof(path_gf12), "%s/selftest_gf2_12_mul.kmb", out_dir);
    {
        qua_ele_gf12_rank42_t gf12_plan;
        if (!qua_ele_build_gf12_rank42(0x1053u, &gf12_plan)) {
            fprintf(stderr, "Failed to synthesize Rank-42 GF(2^12) bilinear decomposition!\n");
            return 1;
        }
        qua_builder_t kb;
        qua_builder_init(&kb);
        uint32_t qa[12], qc[12], qout[12], qanc = 36u;
        for (uint32_t i = 0; i < 12u; ++i) {
            qa[i]   = i;
            qc[i]   = 12u + i;
            qout[i] = 24u + i;
        }
        qua_add_quantum_reg(&kb, "a", qa, 12u);
        qua_add_quantum_reg(&kb, "b", qc, 12u);
        qua_add_quantum_reg(&kb, "out", qout, 12u);
        qua_add_quantum_reg(&kb, "@clean", &qanc, 1u);
        qua_ele_emit_gf2_12_rank42_mul(&kb, qa, qc, qout, qanc, 0u, &gf12_plan);
        if (qua_builder_write_kmb(&kb, path_gf12) != 0) {
            fprintf(stderr, "Failed to write %s\n", path_gf12);
            return 1;
        }
        printf("[2/8] Emitted %s | Q=%u (1 ancilla!), CCX=%llu vs 54 Karatsuba / 144 naive (%.2fx Toffoli compression)\n",
               path_gf12, kb.num_qubits, (unsigned long long)kb.count_toffoli, 144.0 / (double)kb.count_toffoli);
        qua_builder_free(&kb);
    }

    /* 3b. Emit Rank-15 Chudnovsky P^1(GF(4)) GF(64) Multiplier (.kmb) */
    char path_gf64[512];
    snprintf(path_gf64, sizeof(path_gf64), "%s/selftest_gf64_mul.kmb", out_dir);
    {
        qua_ele_gf64_rank15_t gf64_plan;
        if (!qua_ele_build_gf64_rank15(qua_dqi_gf64_mul, &gf64_plan)) {
            fprintf(stderr, "Failed to synthesize Rank-15 GF(64) bilinear decomposition!\n");
            return 1;
        }
        qua_builder_t kb;
        qua_builder_init(&kb);
        uint32_t qa[6], qb6[6], qout[6], qanc = 18u;
        for (uint32_t i = 0u; i < 6u; ++i) {
            qa[i]   = i;
            qb6[i]  = 6u + i;
            qout[i] = 12u + i;
        }
        qua_add_quantum_reg(&kb, "a", qa, 6u);
        qua_add_quantum_reg(&kb, "b", qb6, 6u);
        qua_add_quantum_reg(&kb, "out", qout, 6u);
        qua_add_quantum_reg(&kb, "@clean", &qanc, 1u);
        qua_ele_emit_gf64_rank15_mul(&kb, qa, qb6, qout, qanc, 0u, &gf64_plan);
        if (qua_builder_write_kmb(&kb, path_gf64) != 0) {
            fprintf(stderr, "Failed to write %s\n", path_gf64);
            return 1;
        }
        printf("[3/8] Emitted %s | Q=%u (1 ancilla!), CCX=%llu vs 18 Karatsuba / 36 naive\n",
               path_gf64, kb.num_qubits, (unsigned long long)kb.count_toffoli);
        qua_builder_free(&kb);
    }

    /* 3c. Emit Virtual-Slot-0 QROAM + 0-CCX MBUC Phase-Fixup Self-Test (.kmb) */
    char path_qroam[512];
    snprintf(path_qroam, sizeof(path_qroam), "%s/selftest_qroam_vslot0.kmb", out_dir);
    {
        uint64_t tab37[37];
        for (uint32_t r = 0u; r < 37u; ++r) {
            /* Deliberately keep bit 11 identically 0 to test constant-0 column pruning + non-power-of-2 rows */
            tab37[r] = ((uint64_t)((r * 0x5BDu + 0x137u) & 0x7FFu));
        }
        qua_bitop_qroam_res_t qres;
        qua_bitop_qroam_vslot0_resources(tab37, 37u, 1u, 6u, 12u, 2u, 2u, 0u, &qres);
        uint32_t anc_need = qres.anc_needed;
        uint64_t est_tof = qres.total_ccx;
        qua_builder_t kb;
        qua_builder_init(&kb);
        uint32_t sel_q[6], tgt_q[12], anc_q[128];
        for (uint32_t i = 0u; i < 6u; ++i) sel_q[i] = i;
        for (uint32_t i = 0u; i < 12u; ++i) tgt_q[i] = 6u + i;
        for (uint32_t i = 0u; i < anc_need; ++i) anc_q[i] = 18u + i;
        qua_add_quantum_reg(&kb, "selection", sel_q, 6u);
        qua_add_quantum_reg(&kb, "target", tgt_q, 12u);
        qua_add_quantum_reg(&kb, "@clean", anc_q, anc_need);
        uint32_t next_cbit = 0u;
        qua_bitop_emit_qroam_vslot0(&kb, tab37, 37u, 1u, sel_q, 6u, tgt_q, 12u, 2u, 2u, anc_q, anc_need, &next_cbit);
        if (qua_builder_write_kmb(&kb, path_qroam) != 0) {
            fprintf(stderr, "Failed to write %s\n", path_qroam);
            return 1;
        }
        printf("[4/8] Emitted %s | Q=%u (anc=%u), CCX=%llu (est=%llu), HMR=%llu\n",
               path_qroam, kb.num_qubits, anc_need,
               (unsigned long long)kb.count_toffoli, (unsigned long long)est_tof,
               (unsigned long long)kb.count_hmr);
        qua_builder_free(&kb);
    }

    /* 4. Emit 128-bit FixedPrecisionAngle128 Z_POW Phase Circuit (.kmb) & Plan DF-THC QROM */
    char path_zpow[512];
    snprintf(path_zpow, sizeof(path_zpow), "%s/selftest_zpow128.kmb", out_dir);
    qua_dfthc_plan_t thc_plan = qua_plan_dfthc_angle_qrom(1024u, 16u);
    {
        qua_builder_t kb;
        qua_builder_init(&kb);
        uint32_t q[4] = {0u, 1u, 2u, 3u};
        qua_add_quantum_reg(&kb, "x", q, 4u);
        qua_z_pow_turns_u128(&kb, 0u, ((flow_u128)0x2000000000000000ull << 64));
        qua_z_pow_turns_u128(&kb, 1u, ((flow_u128)0x6000000000000000ull << 64));
        qua_z_pow_turns_u128(&kb, 2u, ((flow_u128)0x4000000000000000ull << 64));
        qua_z_pow_turns_u128(&kb, 3u, ((flow_u128)0x123456789ABCDEF0ull << 64) | 0xFEDCBA9876543210ull);
        if (qua_builder_write_kmb(&kb, path_zpow) != 0) {
            fprintf(stderr, "Failed to write %s\n", path_zpow);
            return 1;
        }
        printf("[5/8] Emitted %s | Q=%u, Z_POW=%llu, T-like=%llu | Virtual-Slot-0 DF-THC(1024x16b) B*=%u, Q=%u<=128, EstTof=%llu\n",
               path_zpow, kb.num_qubits, (unsigned long long)kb.count_z_pow, (unsigned long long)kb.count_t_gate,
               thc_plan.best_b, thc_plan.total_qubits, (unsigned long long)thc_plan.est_toffolis);
        qua_builder_free(&kb);
    }

    /* 4b. Emit HMR_Dual Parity-LUT Uncompute Self-Test (.kmb) */
    char path_hmr_lut[512];
    snprintf(path_hmr_lut, sizeof(path_hmr_lut), "%s/selftest_hmr_lut.kmb", out_dir);
    {
        uint64_t tab37[37];
        for (uint32_t r = 0u; r < 37u; ++r) {
            tab37[r] = (uint64_t)((r * 0x5BDu + 0x137u) & 0xFFFu);
        }
        qua_bitop_qroam_res_t qres;
        qua_bitop_qroam_vslot0_resources(tab37, 37u, 1u, 6u, 12u, 2u, 2u, 0u, &qres);
        qua_bitop_hmr_lut_res_t hres;
        qua_bitop_hmr_parity_lut_resources(37u, 6u, 12u, 0u, &hres);
        if (hres.p_peel != 2u || hres.pool_anc_needed != 0u) {
            fprintf(stderr, "[ERROR] hmr_parity_lut expected p_peel=2, pool_anc=0, got p=%u, pool=%u\n",
                    hres.p_peel, hres.pool_anc_needed);
            return 1;
        }
        if (12u + qres.anc_needed > 256u) {
            fprintf(stderr, "[ERROR] hmr_lut ancilla overflow: %u\n", qres.anc_needed);
            return 1;
        }
        uint64_t est_ccx = qres.total_ccx + hres.total_ccx;

        qua_builder_t kb;
        qua_builder_init(&kb);
        uint32_t sel_q[6], out_q[12], clean_q[256];
        for (uint32_t i = 0u; i < 6u; ++i) sel_q[i] = i;
        for (uint32_t i = 0u; i < 12u; ++i) out_q[i] = 6u + i;
        for (uint32_t i = 0u; i < 12u + qres.anc_needed; ++i) clean_q[i] = 18u + i;
        const uint32_t *tgt_q = clean_q;
        const uint32_t *anc_q = clean_q + 12u;
        qua_add_quantum_reg(&kb, "selection", sel_q, 6u);
        qua_add_quantum_reg(&kb, "out", out_q, 12u);
        qua_add_quantum_reg(&kb, "@clean", clean_q, 12u + qres.anc_needed);

        uint32_t next_cbit = 0u;
        qua_bitop_emit_qroam_vslot0(&kb, tab37, 37u, 1u, sel_q, 6u, tgt_q, 12u, 2u, 2u, anc_q, qres.anc_needed, &next_cbit);
        for (uint32_t i = 0u; i < 12u; ++i) qua_cx(&kb, tgt_q[i], out_q[i]);
        qua_bitop_emit_hmr_parity_lut(&kb, tab37, 37u, sel_q, 6u, tgt_q, 12u, NULL, 0u, &next_cbit);

        if (kb.count_toffoli != est_ccx) {
            fprintf(stderr, "[ERROR] %s CCX mismatch: actual=%llu vs est=%llu (qroam=%llu + hmr=%llu)\n",
                    path_hmr_lut, (unsigned long long)kb.count_toffoli, (unsigned long long)est_ccx,
                    (unsigned long long)qres.total_ccx, (unsigned long long)hres.total_ccx);
            qua_builder_free(&kb);
            return 1;
        }
        if (qua_builder_write_kmb(&kb, path_hmr_lut) != 0) {
            fprintf(stderr, "Failed to write %s\n", path_hmr_lut);
            qua_builder_free(&kb);
            return 1;
        }
        printf("[6/8] Emitted %s | Q=%u (tgt=12, anc=%u, p_peel=%u), CCX=%llu (est=%llu = %llu+%llu), HMR=%llu\n",
               path_hmr_lut, kb.num_qubits, qres.anc_needed, hres.p_peel,
               (unsigned long long)kb.count_toffoli, (unsigned long long)est_ccx,
               (unsigned long long)qres.total_ccx, (unsigned long long)hres.total_ccx,
               (unsigned long long)kb.count_hmr);
        qua_builder_free(&kb);
    }

    /* 4c. Emit LinF2_Fuse GL_12(F_2) Composition + In-Place Synthesis Self-Test (.kmb) */
    char path_linf2[512];
    snprintf(path_linf2, sizeof(path_linf2), "%s/selftest_linf2.kmb", out_dir);
    {
        uint64_t L1[12], U1[12], L2[12], U2[12], M1[12], M2[12], M21[12];
        for (uint32_t r = 0u; r < 12u; ++r) {
            uint64_t lo_mask = (UINT64_C(1) << r) - 1ull;
            uint64_t hi_mask = 0xFFFull ^ ((UINT64_C(1) << (r + 1u)) - 1ull);
            L1[r] = (UINT64_C(1) << r) | (((uint64_t)(r + 1u) * 0x5ADull) & lo_mask);
            U1[r] = (UINT64_C(1) << r) | (((uint64_t)(r + 3u) * 0x36Bull) & hi_mask);
            L2[r] = (UINT64_C(1) << r) | (((uint64_t)(r + 5u) * 0x6D3ull) & lo_mask);
            U2[r] = (UINT64_C(1) << r) | (((uint64_t)(r + 2u) * 0x59Dull) & hi_mask);
        }
        qua_bitop_compose_linf2_u64(L1, U1, M1, 12u);
        qua_bitop_compose_linf2_u64(U2, L2, M2, 12u);
        qua_bitop_compose_linf2_u64(M2, M1, M21, 12u);

        qua_builder_t kb;
        qua_builder_init(&kb);
        uint32_t x_q[12];
        for (uint32_t i = 0u; i < 12u; ++i) x_q[i] = i;
        qua_add_quantum_reg(&kb, "x", x_q, 12u);
        qua_bitop_emit_linf2_inplace_u64(&kb, x_q, 12u, M21);
        if (kb.num_qubits != 12u || kb.count_toffoli != 0ull) {
            fprintf(stderr, "[ERROR] %s expected Q=12, CCX=0, got Q=%u, CCX=%llu\n",
                    path_linf2, kb.num_qubits, (unsigned long long)kb.count_toffoli);
            qua_builder_free(&kb);
            return 1;
        }
        if (qua_builder_write_kmb(&kb, path_linf2) != 0) {
            fprintf(stderr, "Failed to write %s\n", path_linf2);
            qua_builder_free(&kb);
            return 1;
        }
        printf("[7/8] Emitted %s | Q=%u, CCX=%llu, CX=%zu (GL_12(F_2) M2*M1 in-place)\n",
               path_linf2, kb.num_qubits, (unsigned long long)kb.count_toffoli, kb.num_2q);
        qua_builder_free(&kb);
    }

    /* 4d. Emit Cuccaro 1-Ancilla Controlled 16-bit Adder Self-Test (.kmb) */
    char path_cuccaro[512];
    snprintf(path_cuccaro, sizeof(path_cuccaro), "%s/selftest_cuccaro.kmb", out_dir);
    {
        uint32_t cuc_anc = 0u;
        uint64_t cuc_ccx = 0ull;
        qua_bitop_ctrl_iadd_cuccaro_resources(16u, 16u, &cuc_anc, &cuc_ccx);
        if (cuc_anc != 1u || cuc_ccx != 46ull) {
            fprintf(stderr, "[ERROR] ctrl_iadd_cuccaro(16,16) expected anc=1, ccx=46, got anc=%u, ccx=%llu\n",
                    cuc_anc, (unsigned long long)cuc_ccx);
            return 1;
        }
        qua_builder_t kb;
        qua_builder_init(&kb);
        uint32_t ctrl_q = 0u, lhs_q[16], tgt_q[16], anc_q = 33u;
        for (uint32_t i = 0u; i < 16u; ++i) { lhs_q[i] = 1u + i; tgt_q[i] = 17u + i; }
        qua_add_quantum_reg(&kb, "ctrl", &ctrl_q, 1u);
        qua_add_quantum_reg(&kb, "lhs", lhs_q, 16u);
        qua_add_quantum_reg(&kb, "target", tgt_q, 16u);
        qua_add_quantum_reg(&kb, "@clean", &anc_q, 1u);
        uint32_t next_cbit = 0u;
        qua_bitop_emit_ctrl_iadd_cuccaro(&kb, ctrl_q, lhs_q, 16u, tgt_q, 16u, false, &anc_q, &next_cbit);
        if (kb.num_qubits != 34u || kb.count_toffoli != cuc_ccx) {
            fprintf(stderr, "[ERROR] %s expected Q=34, CCX=%llu, got Q=%u, CCX=%llu\n",
                    path_cuccaro, (unsigned long long)cuc_ccx, kb.num_qubits, (unsigned long long)kb.count_toffoli);
            qua_builder_free(&kb);
            return 1;
        }
        if (qua_builder_write_kmb(&kb, path_cuccaro) != 0) {
            fprintf(stderr, "Failed to write %s\n", path_cuccaro);
            qua_builder_free(&kb);
            return 1;
        }
        printf("[8/8] Emitted %s | Q=%u (1 ancilla vs 15 ripple!), CCX=%llu (est=%llu)\n",
               path_cuccaro, kb.num_qubits, (unsigned long long)kb.count_toffoli, (unsigned long long)cuc_ccx);
        qua_builder_free(&kb);
    }

    printf("[SEARU x FLOW v7.0.1] invsqrt(24b) Split-Stage Plan: k_lut=%u, deg=%u, B=%u, Newton=%u -> eff=%ub, Q=%llu, Tof=%llu, Q*C=%llu (isqrt_chk=0x%llx)\n",
           inv_plan.k_lut, inv_plan.poly_deg, inv_plan.b_swap, inv_plan.newton_steps,
           inv_plan.eff_prec_bits, (unsigned long long)inv_plan.peak_qubits,
           (unsigned long long)inv_plan.toffolis, (unsigned long long)inv_plan.cost_qc,
           (unsigned long long)root_chk);
    }

    /* 5. Pure-Search Bilevel SEARU x FLOW 2-Jet Benchmark I (`invsqrt.kmb`) */
    if (do_invsqrt) {
        char path_inv[512];
        snprintf(path_inv, sizeof(path_inv), "%s/invsqrt.kmb", out_dir);
        if (qua_invsqrt_run_pure_search(path_inv) != 0) {
            return 1;
        }
    }

    /* 6. Pure-Search Bilevel SEARU x FLOW 2-Jet Benchmark II (`dqi.kmb`) */
    if (do_dqi) {
        char path_dqi[512];
        snprintf(path_dqi, sizeof(path_dqi), "%s/dqi.kmb", out_dir);
        if (qua_dqi_run_pure_search(path_dqi) != 0) {
            return 1;
        }
    }

    /* 7. Pure-Search Bilevel SEARU x FLOW 2-Jet Benchmark III (`dfthcsos/angle_k_*.kmb`) */
    if (do_dfthcsos) {
        char path_thc_words[512], path_thc_dir[512];
        snprintf(path_thc_words, sizeof(path_thc_words), "%s/dfthcsos/words.bin", out_dir);
        snprintf(path_thc_dir, sizeof(path_thc_dir), "%s/dfthcsos", out_dir);
        if (qua_dfthc_run_pure_search(path_thc_words, path_thc_dir) != 0) {
            return 1;
        }
    }
    return 0;
}
