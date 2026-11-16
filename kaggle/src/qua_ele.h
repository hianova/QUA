/*
 * qua_ele.h — ELE-Style Quantum Physical-Structure Compiler (SEARU x FLOW 2-Jet)
 * Target: Google Quantum AI — Verifiable Quantum Circuit Compilation (Kickmix v2)
 *
 * Mirrors ELE's 3-Pass Physical-Structure Compilation Architecture for Fault-Tolerant
 * Reversible Quantum Circuits:
 *   - Pass 1 (Reversible SSA Pebbling & Bilevel Retiming):
 *       Co-optimizes eager HMR measurement-based erasure vs unitary recomputation
 *       under exact benchmark cost laws (Q*C, R*Q^2*T, lambda*Q*T) via SEARU + FLOW 2-Jet.
 *   - Pass 2 (Live-Range Interference Graph Coloring):
 *       Packs logical SSA quantum registers into minimal physical qubits Q using
 *       searu_solve_graph_color + formal zero-overlap live-range verification.
 *   - Pass 3 (Algebraic Multiplicative-Complexity & Virtual-Slot-0 QROM Lowering):
 *       Synthesizes exact Rank-42 Chudnovsky P^1(GF(2^3)) bilinear GF(2^12) multipliers
 *       (42 CCX, 1 clean ancilla, 0 HMR phase kickback) with SEARU CNOT ordering,
 *       and Virtual-Slot-0 Pebbled Select-Swap QROM circuits directly into Kickmix .kmb.
 */
#ifndef QUA_ELE_H
#define QUA_ELE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "flow_opt_params.h"
#include "flow.h"
#include "searu.h"
#include "qua_kmb.h"

#define QUA_ELE_MAX_NODES 64u
#define QUA_ELE_GF12_RANK 42u
#define QUA_ELE_GF64_RANK 15u

/* ============================================================================
 * Pass 1: Reversible SSA DAG Pebbling & Cost-Law Co-Optimizer
 * ============================================================================ */
typedef enum {
    QUA_COST_QC   = 0, /* Benchmark I   (invsqrt)  : Q * C */
    QUA_COST_Q2T  = 1, /* Benchmark II  (dqi)      : Q^2 * T */
    QUA_COST_QT   = 2  /* Benchmark III (dfthcsos) : Q * T */
} qua_ele_cost_law_t;

typedef struct {
    uint16_t width_qubits;       /* Logical register width in qubits */
    uint32_t compute_toffolis;   /* Forward synthesis Toffoli cost */
    uint32_t uncompute_toffolis; /* Unitary inverse Toffoli cost (if not HMR-erasable) */
    uint8_t  can_hmr_erase;      /* 1 if phase-corrected 0-Toffoli HMR erasure is valid */
    uint8_t  is_persistent_io;   /* 1 if input/output register must stay live to end */
    uint64_t fanin_mask;         /* Bitmask of predecessor SSA nodes */
} qua_ele_ssa_node_t;

typedef struct {
    uint32_t num_nodes;
    uint32_t base_io_qubits;
    qua_ele_cost_law_t cost_law;
    qua_ele_ssa_node_t nodes[QUA_ELE_MAX_NODES];
} qua_ele_dag_t;

typedef struct {
    uint64_t eager_release_mask; /* Bit i = 1 => uncompute/erase node i at earliest last-use */
    uint32_t birth_step[QUA_ELE_MAX_NODES];
    uint32_t release_step[QUA_ELE_MAX_NODES];
    uint32_t peak_qubits;
    uint64_t total_toffolis;
    uint64_t objective_cost;
} qua_ele_pebble_plan_t;

static inline void qua_ele_eval_pebble_plan(
    const qua_ele_dag_t *dag,
    uint64_t eager_mask,
    qua_ele_pebble_plan_t *plan
) {
    const uint32_t n = dag->num_nodes;
    memset(plan, 0, sizeof(*plan));
    plan->eager_release_mask = eager_mask;

    /* Forward topological steps: 0 .. n-1; reverse cleanup steps: n .. 2n-1 */
    uint32_t last_use[QUA_ELE_MAX_NODES];
    for (uint32_t i = 0u; i < n; ++i) {
        plan->birth_step[i] = i;
        last_use[i] = i;
    }
    for (uint32_t j = 0u; j < n; ++j) {
        uint64_t m = dag->nodes[j].fanin_mask;
        for (uint32_t i = 0u; i < j; ++i) {
            if ((m >> i) & 1ull) {
                if (j > last_use[i]) last_use[i] = j;
            }
        }
    }

    uint64_t toff = 0u;
    for (uint32_t i = 0u; i < n; ++i) {
        const qua_ele_ssa_node_t *nd = &dag->nodes[i];
        toff += nd->compute_toffolis;
        if (nd->is_persistent_io) {
            plan->release_step[i] = 2u * n;
            continue;
        }
        bool eager = ((eager_mask >> i) & 1ull) != 0ull;
        if (eager) {
            /* Eagerly erase/uncompute right after last consumer finishes.
             * If unitary uncompute is required while fanin predecessors were already eagerly
             * released, account for reversible rematerialization overhead of missing fanins. */
            plan->release_step[i] = last_use[i];
            if (!nd->can_hmr_erase) {
                toff += nd->uncompute_toffolis;
                uint64_t fm = nd->fanin_mask;
                for (uint32_t p = 0u; p < i; ++p) {
                    if (((fm >> p) & 1ull) && plan->release_step[p] < last_use[i]) {
                        toff += (uint64_t)dag->nodes[p].compute_toffolis +
                                (dag->nodes[p].can_hmr_erase ? 0u : dag->nodes[p].uncompute_toffolis);
                    }
                }
            }
        } else {
            /* Bennett pebble: hold live until reverse sweep step (2n - 1 - i) */
            plan->release_step[i] = 2u * n - 1u - i;
            if (!nd->can_hmr_erase) {
                toff += nd->uncompute_toffolis;
            }
        }
    }

    uint32_t peak_q = dag->base_io_qubits;
    for (uint32_t t = 0u; t <= 2u * n; ++t) {
        uint32_t live_q = dag->base_io_qubits;
        for (uint32_t i = 0u; i < n; ++i) {
            if (plan->birth_step[i] <= t && t <= plan->release_step[i]) {
                live_q += dag->nodes[i].width_qubits;
            }
        }
        if (live_q > peak_q) peak_q = live_q;
    }

    plan->peak_qubits = peak_q;
    plan->total_toffolis = toff;
    if (dag->cost_law == QUA_COST_Q2T) {
        plan->objective_cost = (uint64_t)peak_q * (uint64_t)peak_q * toff;
    } else {
        plan->objective_cost = (uint64_t)peak_q * toff;
    }
}

static inline uint64_t qua_ele_pebble_dynamic_frust(
    const qua_ele_dag_t *dag,
    const qua_ele_pebble_plan_t *plan
) {
    const uint32_t n = dag->num_nodes;
    uint32_t t_peak = 0u;
    uint32_t peak_live = 0u;
    for (uint32_t t = 0u; t <= 2u * n; ++t) {
        uint32_t live_q = dag->base_io_qubits;
        for (uint32_t i = 0u; i < n; ++i) {
            if (plan->birth_step[i] <= t && t <= plan->release_step[i]) {
                live_q += dag->nodes[i].width_qubits;
            }
        }
        if (live_q > peak_live) {
            peak_live = live_q;
            t_peak = t;
        }
    }

    uint64_t f = 0u;
    for (uint32_t i = 0u; i < n; ++i) {
        const qua_ele_ssa_node_t *nd = &dag->nodes[i];
        if (nd->is_persistent_io || nd->width_qubits == 0u) continue;
        bool live_at_peak = (plan->birth_step[i] <= t_peak) && (t_peak <= plan->release_step[i]);
        bool is_eager = ((plan->eager_release_mask >> i) & 1ull) != 0ull;
        /* If live at peak and can be eagerly erased/uncomputed, or if eager caused remat bloat */
        if (live_at_peak && !is_eager) {
            f |= (UINT64_C(1) << i);
        } else if (is_eager && !nd->can_hmr_erase) {
            uint64_t fm = nd->fanin_mask;
            for (uint32_t p = 0u; p < i; ++p) {
                if (((fm >> p) & 1ull) && plan->release_step[p] < plan->release_step[i]) {
                    f |= (UINT64_C(1) << i) | (UINT64_C(1) << p);
                }
            }
        }
    }
    if (f == 0u) {
        for (uint32_t i = 0u; i < n; ++i) {
            if (!dag->nodes[i].is_persistent_io && dag->nodes[i].width_qubits > 0u) {
                f |= (UINT64_C(1) << i);
            }
        }
    }
    return f;
}

static inline qua_ele_pebble_plan_t qua_ele_solve_pebble_dag_hp(
    const qua_ele_dag_t *dag,
    const searu_hyperparams *hp,
    uint64_t seed,
    uint32_t max_iters
) {
    uint64_t pin_mask = 0u;
    for (uint32_t i = 0u; i < dag->num_nodes; ++i) {
        if (dag->nodes[i].is_persistent_io) pin_mask |= (UINT64_C(1) << i);
    }
    uint64_t init_c[4] = {0u, 0u, 0u, 0u};
    uint64_t pin_c[4]  = {pin_mask, 0u, 0u, 0u};
    searu_result res;
    SEARU_SOLVE_INLINE_HP(&res, dag->num_nodes, init_c, pin_c, seed, max_iters, hp,
                          cand, frust, energy, {
        qua_ele_pebble_plan_t p;
        qua_ele_eval_pebble_plan(dag, cand[0], &p);
        energy = p.objective_cost;
        frust[0] = qua_ele_pebble_dynamic_frust(dag, &p);
    });
    qua_ele_pebble_plan_t best;
    qua_ele_eval_pebble_plan(dag, res.canvas[0], &best);
    return best;
}

static inline qua_ele_pebble_plan_t qua_ele_solve_pebble_dag(
    const qua_ele_dag_t *dag,
    uint64_t seed,
    uint32_t max_iters
) {
    const searu_hyperparams def_hp = SEARU_ACTIVE_OPT_HP;
    return qua_ele_solve_pebble_dag_hp(dag, &def_hp, seed, max_iters);
}

/* ============================================================================
 * Pass 2: SEARU Heterogeneous-Width Live-Range Packing & Disjointness Proof
 * ============================================================================ */
typedef struct {
    uint32_t num_colors_used;
    uint32_t phys_qubits_used;
    uint32_t slot_color[QUA_ELE_MAX_NODES];
    uint32_t phys_offset[QUA_ELE_MAX_NODES];
    bool     verified_disjoint;
} qua_ele_reg_alloc_t;

static inline bool qua_ele_intervals_overlap(
    uint32_t b0, uint32_t r0,
    uint32_t b1, uint32_t r1
) {
    return (b0 <= r1) && (b1 <= r0);
}

static inline uint32_t qua_ele_extract_bank_code(
    const uint64_t canvas[4],
    uint32_t node_idx,
    uint32_t bits_per_node
) {
    uint32_t code = 0u;
    for (uint32_t b = 0u; b < bits_per_node; ++b) {
        uint32_t v = node_idx * bits_per_node + b;
        if (((canvas[v >> 6u] >> (v & 63u)) & 1ull) != 0ull) {
            code |= (1u << b);
        }
    }
    return code;
}

static inline uint64_t qua_ele_eval_hetero_pack(
    const qua_ele_dag_t *dag,
    const qua_ele_pebble_plan_t *plan,
    uint32_t pool_banks,
    uint32_t bits_per_node,
    const uint64_t cand[4],
    uint64_t frust[4]
) {
    const uint32_t n = dag->num_nodes;
    uint32_t bank_of[QUA_ELE_MAX_NODES];
    uint32_t bank_max_w[QUA_ELE_MAX_NODES];
    memset(bank_max_w, 0, sizeof(bank_max_w));

    for (uint32_t u = 0u; u < n; ++u) {
        uint32_t b = qua_ele_extract_bank_code(cand, u, bits_per_node) % pool_banks;
        bank_of[u] = b;
        if (dag->nodes[u].width_qubits > bank_max_w[b]) {
            bank_max_w[b] = dag->nodes[u].width_qubits;
        }
    }

    uint64_t sum_w = 0u;
    for (uint32_t b = 0u; b < pool_banks; ++b) {
        sum_w += (uint64_t)bank_max_w[b];
    }

    uint64_t overlap_cnt = 0u;
    uint64_t conflict_mask = 0u;
    for (uint32_t u = 0u; u < n; ++u) {
        if (dag->nodes[u].width_qubits == 0u) continue;
        for (uint32_t v = u + 1u; v < n; ++v) {
            if (dag->nodes[v].width_qubits == 0u) continue;
            if (bank_of[u] == bank_of[v] &&
                qua_ele_intervals_overlap(plan->birth_step[u], plan->release_step[u],
                                          plan->birth_step[v], plan->release_step[v])) {
                ++overlap_cnt;
                conflict_mask |= (UINT64_C(1) << u) | (UINT64_C(1) << v);
            }
        }
    }

    frust[0] = 0u;
    frust[1] = 0u;
    frust[2] = 0u;
    frust[3] = 0u;
    for (uint32_t u = 0u; u < n; ++u) {
        if (dag->nodes[u].width_qubits == 0u) continue;
        bool bad = ((conflict_mask >> u) & 1ull) != 0ull;
        bool padded = (bank_max_w[bank_of[u]] > dag->nodes[u].width_qubits);
        if (bad || (overlap_cnt == 0u && padded)) {
            for (uint32_t b = 0u; b < bits_per_node; ++b) {
                uint32_t v = u * bits_per_node + b;
                frust[v >> 6u] |= (UINT64_C(1) << (v & 63u));
            }
        }
    }

    return (overlap_cnt << 24u) + sum_w;
}

static inline qua_ele_reg_alloc_t qua_ele_alloc_qubits_hp(
    const qua_ele_dag_t *dag,
    const qua_ele_pebble_plan_t *plan,
    const searu_hyperparams *hp,
    uint64_t seed
) {
    qua_ele_reg_alloc_t out;
    memset(&out, 0, sizeof(out));
    const uint32_t n = dag->num_nodes;
    if (n == 0u) {
        out.verified_disjoint = true;
        return out;
    }

    /* 1. Greedy chordal interval left-edge coloring as valid initial canvas */
    uint32_t max_c = 0u;
    for (uint32_t u = 0u; u < n; ++u) {
        if (dag->nodes[u].width_qubits == 0u) {
            out.slot_color[u] = 0u;
            continue;
        }
        uint64_t used = 0u;
        for (uint32_t v = 0u; v < u; ++v) {
            if (dag->nodes[v].width_qubits == 0u) continue;
            if (qua_ele_intervals_overlap(plan->birth_step[u], plan->release_step[u],
                                          plan->birth_step[v], plan->release_step[v])) {
                used |= (UINT64_C(1) << out.slot_color[v]);
            }
        }
        uint32_t c = 0u;
        while ((used >> c) & 1ull) ++c;
        out.slot_color[u] = c;
        if (c + 1u > max_c) max_c = c + 1u;
    }
    if (max_c == 0u) max_c = 1u;
    out.num_colors_used = max_c;

    /* 2. SEARU heterogeneous-width bank packing (mirrors ele_bootstrap.c Pass 2B) */
    if (max_c >= 2u && max_c <= 16u) {
        uint32_t bits_per_node = 1u;
        while ((1u << bits_per_node) < max_c) ++bits_per_node;
        if (n * bits_per_node <= SEARU_CANVAS_BITS) {
            uint64_t init_c[4] = {0u, 0u, 0u, 0u};
            uint64_t pin_c[4]  = {0u, 0u, 0u, 0u};
            for (uint32_t u = 0u; u < n; ++u) {
                uint32_t code = out.slot_color[u];
                for (uint32_t b = 0u; b < bits_per_node; ++b) {
                    uint32_t v = u * bits_per_node + b;
                    if ((code >> b) & 1u) {
                        init_c[v >> 6u] |= (UINT64_C(1) << (v & 63u));
                    }
                    if (dag->nodes[u].width_qubits == 0u) {
                        pin_c[v >> 6u] |= (UINT64_C(1) << (v & 63u));
                    }
                }
            }
            searu_result res;
            SEARU_SOLVE_INLINE_HP(&res, n * bits_per_node, init_c, pin_c, seed, 4096u, hp,
                                  cand, frust, energy, {
                energy = qua_ele_eval_hetero_pack(dag, plan, max_c, bits_per_node, cand, frust);
            });
            if (res.best_energy < (UINT64_C(1) << 24u)) {
                uint32_t remap[QUA_ELE_MAX_NODES];
                uint32_t compact = 0u;
                for (uint32_t b = 0u; b < max_c; ++b) remap[b] = UINT32_MAX;
                for (uint32_t u = 0u; u < n; ++u) {
                    if (dag->nodes[u].width_qubits == 0u) continue;
                    uint32_t b = qua_ele_extract_bank_code(res.canvas, u, bits_per_node) % max_c;
                    if (remap[b] == UINT32_MAX) remap[b] = compact++;
                    out.slot_color[u] = remap[b];
                }
                out.num_colors_used = (compact > 0u) ? compact : 1u;
            }
        }
    }

    /* 3. Size each physical color bank to max width of registers assigned to that color */
    uint32_t bank_width[QUA_ELE_MAX_NODES];
    uint32_t bank_base[QUA_ELE_MAX_NODES];
    memset(bank_width, 0, sizeof(bank_width));
    for (uint32_t u = 0u; u < n; ++u) {
        uint32_t c = out.slot_color[u];
        if (dag->nodes[u].width_qubits > bank_width[c]) {
            bank_width[c] = dag->nodes[u].width_qubits;
        }
    }
    uint32_t cursor = dag->base_io_qubits;
    for (uint32_t c = 0u; c < out.num_colors_used; ++c) {
        bank_base[c] = cursor;
        cursor += bank_width[c];
    }
    for (uint32_t u = 0u; u < n; ++u) {
        out.phys_offset[u] = bank_base[out.slot_color[u]];
    }
    out.phys_qubits_used = cursor;

    /* 4. Formal Disjointness Verification (mirrors ele_bootstrap.c Pass 2) */
    out.verified_disjoint = true;
    for (uint32_t u = 0u; u < n; ++u) {
        uint32_t u_lo = out.phys_offset[u];
        uint32_t u_hi = u_lo + dag->nodes[u].width_qubits;
        for (uint32_t v = u + 1u; v < n; ++v) {
            if (!qua_ele_intervals_overlap(plan->birth_step[u], plan->release_step[u],
                                           plan->birth_step[v], plan->release_step[v])) {
                continue;
            }
            uint32_t v_lo = out.phys_offset[v];
            uint32_t v_hi = v_lo + dag->nodes[v].width_qubits;
            if (u_lo < v_hi && v_lo < u_hi) {
                out.verified_disjoint = false;
            }
        }
    }
    return out;
}

static inline qua_ele_reg_alloc_t qua_ele_alloc_qubits(
    const qua_ele_dag_t *dag,
    const qua_ele_pebble_plan_t *plan,
    uint64_t seed
) {
    const searu_hyperparams def_hp = SEARU_ACTIVE_OPT_HP;
    return qua_ele_alloc_qubits_hp(dag, plan, &def_hp, seed);
}

/* ============================================================================
 * Pass 0 / Outer Meta-SEARU: Bilevel Self-Optimization + FLOW 2-Jet Polishing
 * (Mirrors flow.c cmd_self_opt + polish_cool_2jet)
 * ============================================================================ */
typedef uint64_t (*qua_ele_inner_eval_fn)(const searu_hyperparams *hp, uint64_t seed, void *ctx);

typedef struct {
    uint64_t base_train_cost;
    uint64_t opt_train_cost;
    uint64_t base_val_cost;
    uint64_t opt_val_cost;
    bool     accepted;
} qua_ele_self_opt_stats_t;

static inline void qua_ele_decode_inl_genome(uint64_t word0, searu_hyperparams *hp) {
    static const uint32_t WALK_TAB[4] = {1u, 3u, 7u, 15u};
    uint32_t cmin_g = (uint32_t)(word0 & 3u);
    uint32_t cmax_g = (uint32_t)((word0 >> 2u) & 3u);
    uint32_t stag_g = (uint32_t)((word0 >> 4u) & 3u);
    uint32_t bst_g  = (uint32_t)((word0 >> 6u) & 3u);
    uint32_t wlk_g  = (uint32_t)((word0 >> 8u) & 3u);
    uint32_t itmp_g = (uint32_t)((word0 >> 10u) & 7u);
    uint32_t cool_g = (uint32_t)((word0 >> 13u) & 15u);
    uint32_t rstag_g = (uint32_t)((word0 >> 17u) & 3u);
    uint32_t rtmp_g = (uint32_t)((word0 >> 19u) & 3u);

    hp->chaos_min_q8      = 96u + cmin_g * 32u;
    hp->chaos_max_q8      = 256u + cmax_g * 64u;
    if (hp->chaos_max_q8 <= hp->chaos_min_q8) hp->chaos_max_q8 = hp->chaos_min_q8 + 64u;
    hp->stag_shift        = stag_g;
    hp->barrier_boost_q8  = 16u + bst_g * 16u;
    hp->inl_walk_mask     = WALK_TAB[wlk_g];
    /* Exact Q12 fixed-point scaling matching searu_crucible_accept (delta << 12) */
    hp->inl_init_temp_q12 = (UINT64_C(64) << itmp_g) << 12u;
    hp->inl_cool_q16      = UINT64_C(63600) + (uint64_t)cool_g * UINT64_C(120);
    hp->inl_reheat_stag   = 16u << rstag_g;
    hp->inl_reheat_q12    = (UINT64_C(8) << rtmp_g) << 12u;
}

static inline uint64_t qua_ele_eval_hp_suite(
    qua_ele_inner_eval_fn fn,
    void *ctx,
    const searu_hyperparams *hp,
    const uint64_t *seeds,
    uint32_t n_seeds
) {
    uint64_t total = 0u;
    for (uint32_t i = 0u; i < n_seeds; ++i) {
        total += fn(hp, seeds[i], ctx);
    }
    return total;
}

static inline void qua_ele_polish_cool_2jet(
    qua_ele_inner_eval_fn fn,
    void *ctx,
    searu_hyperparams *hp,
    const uint64_t *seeds,
    uint32_t n_seeds
) {
    uint64_t cur_cost = qua_ele_eval_hp_suite(fn, ctx, hp, seeds, n_seeds);
    for (int64_t h = INT64_C(256); h >= INT64_C(32); h >>= 1u) {
        uint64_t orig = hp->inl_cool_q16;
        searu_hyperparams hp_p = *hp;
        searu_hyperparams hp_m = *hp;
        hp_p.inl_cool_q16 = (orig + (uint64_t)h <= UINT64_C(65480)) ? (orig + (uint64_t)h) : UINT64_C(65480);
        hp_m.inl_cool_q16 = (orig > (uint64_t)h + UINT64_C(60000)) ? (orig - (uint64_t)h) : UINT64_C(60000);
        uint64_t ep = qua_ele_eval_hp_suite(fn, ctx, &hp_p, seeds, n_seeds);
        uint64_t em = qua_ele_eval_hp_suite(fn, ctx, &hp_m, seeds, n_seeds);
        int64_t g1 = ((int64_t)ep - (int64_t)em) >> 1u;
        int64_t g2 = (int64_t)ep - 2 * (int64_t)cur_cost + (int64_t)em;
        if (g2 > 0 && g1 != 0) {
            uint64_t st[8] = {
                (uint64_t)(g1 * FLOW_Q16_ONE), (uint64_t)(g2 * FLOW_Q16_ONE),
                0ULL, 0ULL, 0ULL, 0ULL, 0ULL, 0ULL
            };
            uint64_t bd[8] = {
                0ULL, 0ULL, 0ULL,
                ((uint64_t)(uint32_t)(int32_t)(-4 * FLOW_Q16_ONE) << 32u) | (uint64_t)(uint32_t)(int32_t)(4 * FLOW_Q16_ONE),
                ((uint64_t)(uint32_t)(int32_t)(2 * FLOW_Q16_ONE) << 32u) | 1ULL,
                0ULL, 0ULL, 0ULL
            };
            flow_step(st, st, bd);
            int64_t step_c = (int64_t)orig + (((int64_t)st[3] * h) >> FLOW_Q16_SHIFT);
            if (step_c >= INT64_C(60000) && step_c <= INT64_C(65480)) {
                hp->inl_cool_q16 = (uint64_t)step_c;
                uint64_t et = qua_ele_eval_hp_suite(fn, ctx, hp, seeds, n_seeds);
                if (et < cur_cost) {
                    cur_cost = et;
                } else {
                    hp->inl_cool_q16 = orig;
                }
            }
        }
        if (ep < cur_cost) {
            hp->inl_cool_q16 = hp_p.inl_cool_q16;
            cur_cost = ep;
        } else if (em < cur_cost) {
            hp->inl_cool_q16 = hp_m.inl_cool_q16;
            cur_cost = em;
        }
    }
}

static inline searu_hyperparams qua_ele_self_opt_hp(
    qua_ele_inner_eval_fn fn,
    void *ctx,
    uint64_t meta_seed,
    uint32_t meta_iters,
    qua_ele_self_opt_stats_t *out_stats
) {
    static const uint64_t train_seeds[4] = {
        UINT64_C(0x13579BDF2468ACE1), UINT64_C(0x9E3779B97F4A7C15),
        UINT64_C(0xBF58476D1CE4E5B9), UINT64_C(0x94D049BB133111EB),
    };
    static const uint64_t val_seeds[4] = {
        UINT64_C(0x5EED5A71A55A3C11), UINT64_C(0xCAFEBABE12345678),
        UINT64_C(0x6A09E667F3BCC908), UINT64_C(0xBB67AE8584CAA73B),
    };

    const searu_hyperparams base_hp = SEARU_ACTIVE_OPT_HP;
    uint64_t base_train = qua_ele_eval_hp_suite(fn, ctx, &base_hp, train_seeds, 4u);
    uint64_t base_val   = qua_ele_eval_hp_suite(fn, ctx, &base_hp, val_seeds, 4u);

    /* Genome init matching SEARU_ACTIVE_OPT_HP Q12 encoding */
    const uint64_t init_g[4] = {
        (1ULL) | (1ULL << 2u) | (0ULL << 4u) | (0ULL << 6u) |
        (0ULL << 8u) | (3ULL << 10u) | (1ULL << 13u) | (2ULL << 17u) | (2ULL << 19u),
        0ULL, 0ULL, 0ULL
    };
    searu_result res_g;
    SEARU_SOLVE_INLINE_HP(&res_g, 21u, init_g, NULL, meta_seed, meta_iters, &base_hp,
                          cand, frust, energy, {
        searu_hyperparams chp = base_hp;
        qua_ele_decode_inl_genome(cand[0], &chp);
        energy = qua_ele_eval_hp_suite(fn, ctx, &chp, train_seeds, 4u);
        frust[0] = UINT64_C(0x1FFFFF);
    });

    searu_hyperparams pre_polish_hp = base_hp;
    qua_ele_decode_inl_genome(res_g.canvas[0], &pre_polish_hp);
    searu_hyperparams opt_hp = pre_polish_hp;
    qua_ele_polish_cool_2jet(fn, ctx, &opt_hp, train_seeds, 4u);

    uint64_t opt_train = qua_ele_eval_hp_suite(fn, ctx, &opt_hp, train_seeds, 4u);
    uint64_t opt_val   = qua_ele_eval_hp_suite(fn, ctx, &opt_hp, val_seeds, 4u);
    if (opt_val > base_val) {
        uint64_t pre_train = qua_ele_eval_hp_suite(fn, ctx, &pre_polish_hp, train_seeds, 4u);
        uint64_t pre_val   = qua_ele_eval_hp_suite(fn, ctx, &pre_polish_hp, val_seeds, 4u);
        if (pre_train <= base_train && pre_val <= base_val) {
            opt_hp    = pre_polish_hp;
            opt_train = pre_train;
            opt_val   = pre_val;
        }
    }
    bool accepted = (opt_train <= base_train) && (opt_val <= base_val);

    if (out_stats != NULL) {
        out_stats->base_train_cost = base_train;
        out_stats->opt_train_cost  = accepted ? opt_train : base_train;
        out_stats->base_val_cost   = base_val;
        out_stats->opt_val_cost    = accepted ? opt_val : base_val;
        out_stats->accepted        = accepted;
    }
    return accepted ? opt_hp : base_hp;
}

/* ============================================================================
 * Pass 3A: Exact Rank-42 Chudnovsky P^1(GF(2^3)) GF(2^12) Bilinear Synthesizer
 * ============================================================================
 * Decomposes multiplication in GF(2)[x]/(prim_poly) (deg = 12) via the unique
 * degree-4 tower extension over GF(2^3):
 *   - |P^1(GF(2^3))| = 8 + 1 = 9 >= 7 rational evaluation points!
 *   - Interpolating degree-3 x degree-3 -> degree-6 polynomials over GF(2^3) at
 *     7 points {0, 1, t, t^2, t+1, t^2+t, inf} uses 7 multiplications in GF(2^3).
 *   - Each GF(2^3) = GF(2)[t]/(t^3 + t + 1) product uses 6-AND Karatsuba.
 *   - Total multiplicative complexity: 7 * 6 = 42 Toffoli (CCX) gates!
 *   - Executed in-place on `a[12], b[12]` with 1 clean ancilla + 0-CCX HMR + CZ_IF
 *     phase cleanup => Q = 37 qubits, CCX = 42 (3.43x compression vs 144 naive,
 *     22.2% fewer CCX than 54-CCX 2-level Karatsuba at identical Q = 37!).
 */
typedef struct {
    uint16_t l_mask[QUA_ELE_GF12_RANK]; /* 12-bit linear form on `a` */
    uint16_t r_mask[QUA_ELE_GF12_RANK]; /* 12-bit linear form on `b` */
    uint16_t w_mask[QUA_ELE_GF12_RANK]; /* 12-bit output fanout mask on `out` */
    uint32_t order[QUA_ELE_GF12_RANK];  /* SEARU CNOT-transition optimized schedule */
} qua_ele_gf12_rank42_t;

static inline uint16_t qua_gf12_mul_ref(uint16_t a, uint16_t b, uint32_t prim_poly) {
    uint32_t res = 0u;
    uint32_t x = a & 0xFFFu;
    for (uint32_t i = 0u; i < 12u; ++i) {
        if ((b >> i) & 1u) res ^= x;
        x <<= 1u;
        if (x & (1u << 12u)) x ^= prim_poly;
    }
    return (uint16_t)(res & 0xFFFu);
}

static inline uint8_t qua_gf8_mul(uint8_t a, uint8_t b) {
    /* GF(2^3) = GF(2)[t] / (t^3 + t + 1), modulus = 0xB (binary 1011) */
    uint8_t r = 0u;
    for (uint32_t i = 0u; i < 3u; ++i) {
        if ((b >> i) & 1u) r ^= a;
        a = (uint8_t)((a << 1u) ^ ((a & 4u) ? 0xBu : 0u));
    }
    return (uint8_t)(r & 7u);
}

static inline uint8_t qua_gf8_inv(uint8_t a) {
    /* a^(2^3 - 2) = a^6 in GF(2^3)* */
    uint8_t a2 = qua_gf8_mul(a, a);
    uint8_t a4 = qua_gf8_mul(a2, a2);
    return qua_gf8_mul(a2, a4);
}

static inline bool qua_ele_build_gf12_rank42(uint32_t prim_poly, qua_ele_gf12_rank42_t *out) {
    memset(out, 0, sizeof(*out));

    /* 1. Locate t in GF(2^12) satisfying t^3 + t + 1 = 0 */
    uint16_t t_elem = 0u;
    for (uint16_t cand = 2u; cand < 4096u; ++cand) {
        uint16_t c2 = qua_gf12_mul_ref(cand, cand, prim_poly);
        uint16_t c3 = qua_gf12_mul_ref(c2, cand, prim_poly);
        if ((c3 ^ cand ^ 1u) == 0u) {
            t_elem = cand;
            break;
        }
    }
    if (t_elem == 0u) return false;

    /* 2. Find generator W such that {t^i * W^j : 0<=i<3, 0<=j<4} spans GF(2^12) */
    uint16_t t_pow[3] = {1u, t_elem, qua_gf12_mul_ref(t_elem, t_elem, prim_poly)};
    uint16_t basis[12];
    uint16_t inv_cols[12]; /* Maps standard basis bit s -> 12-bit tower coordinates */
    bool found_basis = false;

    for (uint16_t w_cand = 2u; w_cand < 4096u && !found_basis; ++w_cand) {
        uint16_t w_pow = 1u;
        for (uint32_t j = 0u; j < 4u; ++j) {
            for (uint32_t i = 0u; i < 3u; ++i) {
                basis[3u * j + i] = qua_gf12_mul_ref(t_pow[i], w_pow, prim_poly);
            }
            w_pow = qua_gf12_mul_ref(w_pow, w_cand, prim_poly);
        }
        /* Invert 12x12 binary matrix M where column col = basis[col] */
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
            /* inv[coord] gives the linear form on standard input bits! */
            for (uint32_t coord = 0u; coord < 12u; ++coord) inv_cols[coord] = inv[coord];
            /* Also extend basis powers W^0 .. W^6 for degree-6 product reconstruction! */
            found_basis = true;
            /* Store w_cand in t_pow[0]UNUSED? No, compute w_pows7 directly right here: */
            uint16_t w7[7];
            w7[0] = 1u;
            for (uint32_t d = 1u; d < 7u; ++d) w7[d] = qua_gf12_mul_ref(w7[d - 1u], w_cand, prim_poly);

            /* 3. Build 7x7 Vandermonde evaluation & inverse over GF(2^3) at
             *    pts = {0, 1, 2(t), 4(t^2), 3(t+1), 6(t^2+t), inf} */
            static const uint8_t pts[6] = {0u, 1u, 2u, 4u, 3u, 6u};
            uint8_t vmat[7][7], vinv[7][7];
            memset(vmat, 0, sizeof(vmat));
            memset(vinv, 0, sizeof(vinv));
            for (uint32_t r = 0u; r < 6u; ++r) {
                uint8_t p_pow = 1u;
                for (uint32_t c = 0u; c < 7u; ++c) {
                    vmat[r][c] = p_pow;
                    p_pow = qua_gf8_mul(p_pow, pts[r]);
                }
                vinv[r][r] = 1u;
            }
            vmat[6][6] = 1u; /* Point at infinity extracts degree-6 leading coeff a3*b3 */
            vinv[6][6] = 1u;

            /* Gauss-Jordan over GF(2^3) */
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
                uint8_t scale = qua_gf8_inv(vmat[col][col]);
                for (uint32_t k = 0u; k < 7u; ++k) {
                    vmat[col][k] = qua_gf8_mul(vmat[col][k], scale);
                    vinv[col][k] = qua_gf8_mul(vinv[col][k], scale);
                }
                for (uint32_t r = 0u; r < 7u; ++r) {
                    if (r != col && vmat[r][col] != 0u) {
                        uint8_t f = vmat[r][col];
                        for (uint32_t k = 0u; k < 7u; ++k) {
                            vmat[r][k] ^= qua_gf8_mul(f, vmat[col][k]);
                            vinv[r][k] ^= qua_gf8_mul(f, vinv[col][k]);
                        }
                    }
                }
            }

            /* 4. Emit 7 * 6 = 42 bilinear triples!
             *    6-AND Karatsuba for u, v in GF(2^3) = GF(2)[t]/(t^3 + t + 1):
             *      p0 = u0 * v0                 -> 1 ^ t        (mask 0b011 = 3? Wait: check exact!)
             *    Let's compute the 3-bit output of each of the 6 Karatsuba products algebraically!
             *      K0 = u0 * v0                 contributes to deg0 (1), and cancels from deg1 (t) & deg2 (t^2) => 1 ^ t ^ t^2
             *      K1 = u1 * v1                 contributes to deg2 (t^2), cancels from deg1 (t) & deg3 (t+1)   => t^2 ^ t ^ (t^1 ^ 1) = t^2 ^ 1
             *      K2 = u2 * v2                 contributes to deg4 (t^2+t), cancels from deg2 (t^2) & deg3 (t+1)=> (t^2^t) ^ t^2 ^ (t^1) = 1
             *      K3 = (u0^u1)*(v0^v1)         contributes to deg1 (t)
             *      K4 = (u0^u2)*(v0^v2)         contributes to deg2 (t^2)
             *      K5 = (u1^u2)*(v1^v2)         contributes to deg3 (t+1) */
            static const uint8_t k6_u[6]   = {0x1u, 0x2u, 0x4u, 0x3u, 0x5u, 0x6u};
            static const uint8_t k6_v[6]   = {0x1u, 0x2u, 0x4u, 0x3u, 0x5u, 0x6u};
            static const uint8_t k6_out[6] = {0x7u, 0x5u, 0x1u, 0x2u, 0x4u, 0x3u};

            for (uint32_t m = 0u; m < 7u; ++m) {
                /* Linear forms u_bit[0..2], v_bit[0..2] in standard 12-bit input space */
                uint16_t u_lin[3] = {0u, 0u, 0u};
                uint16_t v_lin[3] = {0u, 0u, 0u};
                if (m < 6u) {
                    uint8_t p_pow = 1u;
                    for (uint32_t j = 0u; j < 4u; ++j) {
                        /* Multiply tower coefficient a_j = (inv_cols[3j+0], inv_cols[3j+1], inv_cols[3j+2]) by p_pow in GF(2^3) */
                        for (uint32_t i = 0u; i < 3u; ++i) {
                            uint8_t prod = qua_gf8_mul((uint8_t)(1u << i), p_pow);
                            for (uint32_t b = 0u; b < 3u; ++b) {
                                if ((prod >> b) & 1u) {
                                    u_lin[b] ^= inv_cols[3u * j + i];
                                    v_lin[b] ^= inv_cols[3u * j + i];
                                }
                            }
                        }
                        p_pow = qua_gf8_mul(p_pow, pts[m]);
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
                    /* k6_out[k] in GF(2^3) at evaluation point m contributes vinv[d][m] * k6_out[k] to coeff c_d of W^d */
                    for (uint32_t d = 0u; d < 7u; ++d) {
                        uint8_t coeff_gf8 = qua_gf8_mul(vinv[d][m], k6_out[k]);
                        for (uint32_t i = 0u; i < 3u; ++i) {
                            if ((coeff_gf8 >> i) & 1u) {
                                wm ^= qua_gf12_mul_ref(t_pow[i], w7[d], prim_poly);
                            }
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

    /* 5. Algebraic Completeness Proof on all 12 x 12 = 144 basis pairs */
    for (uint32_t i = 0u; i < 12u; ++i) {
        for (uint32_t j = 0u; j < 12u; ++j) {
            uint16_t a = (uint16_t)(1u << i);
            uint16_t b = (uint16_t)(1u << j);
            uint16_t got = 0u;
            for (uint32_t p = 0u; p < QUA_ELE_GF12_RANK; ++p) {
                uint32_t la = (uint32_t)__builtin_parity((unsigned)(out->l_mask[p] & a));
                uint32_t rb = (uint32_t)__builtin_parity((unsigned)(out->r_mask[p] & b));
                if (la & rb) got ^= out->w_mask[p];
            }
            if (got != qua_gf12_mul_ref(a, b, prim_poly)) return false;
        }
    }

    /* 6. Greedy nearest-neighbor CNOT ordering on (l_mask, r_mask) pairs */
    bool used[QUA_ELE_GF12_RANK];
    memset(used, 0, sizeof(used));
    uint32_t cur = 0u;
    for (uint32_t step = 0u; step < QUA_ELE_GF12_RANK; ++step) {
        uint32_t best_p = 0u;
        uint32_t best_cost = 9999u;
        for (uint32_t p = 0u; p < QUA_ELE_GF12_RANK; ++p) {
            if (used[p]) continue;
            uint32_t cost = (step == 0u)
                ? (uint32_t)(__builtin_popcount(out->l_mask[p]) + __builtin_popcount(out->r_mask[p]))
                : (uint32_t)(__builtin_popcount((unsigned)(out->l_mask[p] ^ out->l_mask[cur])) +
                             __builtin_popcount((unsigned)(out->r_mask[p] ^ out->r_mask[cur])));
            if (cost < best_cost) {
                best_cost = cost;
                best_p = p;
            }
        }
        used[best_p] = true;
        out->order[step] = best_p;
        cur = best_p;
    }
    return true;
}

static inline void qua_ele_emit_gf2_12_rank42_mul(
    qua_builder_t *b,
    const uint32_t qa[12],
    const uint32_t qb[12],
    const uint32_t qout[12],
    uint32_t anc_prod,
    uint32_t hmr_bit,
    const qua_ele_gf12_rank42_t *plan
) {
    for (uint32_t step = 0u; step < QUA_ELE_GF12_RANK; ++step) {
        uint32_t p = plan->order[step];
        uint16_t lm = plan->l_mask[p];
        uint16_t rm = plan->r_mask[p];
        uint16_t wm = plan->w_mask[p];
        if (lm == 0u || rm == 0u || wm == 0u) continue;

        uint32_t piv_a = (uint32_t)__builtin_ctz((unsigned)lm);
        uint32_t piv_b = (uint32_t)__builtin_ctz((unsigned)rm);

        /* Prepare linear combinations <l_p, a> on qa[piv_a] and <r_p, b> on qb[piv_b] */
        for (uint32_t i = piv_a + 1u; i < 12u; ++i) {
            if ((lm >> i) & 1u) qua_cx(b, qa[i], qa[piv_a]);
        }
        for (uint32_t j = piv_b + 1u; j < 12u; ++j) {
            if ((rm >> j) & 1u) qua_cx(b, qb[j], qb[piv_b]);
        }

        /* 1-CCX temporary AND + CNOT fanout + 0-CCX HMR phase-exact erasure */
        qua_and_compute(b, qa[piv_a], qb[piv_b], anc_prod);
        for (uint32_t k = 0u; k < 12u; ++k) {
            if ((wm >> k) & 1u) qua_cx(b, anc_prod, qout[k]);
        }
        qua_and_uncompute(b, qa[piv_a], qb[piv_b], anc_prod, hmr_bit);

        /* Restore qa[piv_a] and qb[piv_b] cleanly */
        for (uint32_t j = 11u; j > piv_b; --j) {
            if ((rm >> j) & 1u) qua_cx(b, qb[j], qb[piv_b]);
        }
        for (uint32_t i = 11u; i > piv_a; --i) {
            if ((lm >> i) & 1u) qua_cx(b, qa[i], qa[piv_a]);
        }
    }
}

/* ============================================================================
 * Pass 3B: Exact Rank-15 Chudnovsky P^1(GF(4)) GF(64) Bilinear Synthesizer
 * ============================================================================
 * Decomposes any GF(64) multiplication (`mul_fn`) via the unique degree-3
 * extension GF(4)[z]/(z^3 + z + 1):
 *   - |P^1(GF(4))| = 4 + 1 = 5 rational evaluation points {0, 1, w, w^2, inf}.
 *   - Interpolates degree-2 x degree-2 -> degree-4 polynomials over GF(4) using
 *     5 multiplications in GF(4).
 *   - Each GF(4) = GF(2)[w]/(w^2 + w + 1) product uses 3-AND Karatsuba.
 *   - Total multiplicative complexity: 5 * 3 = 15 Toffoli (CCX) gates!
 *     (Saves 3 CCX per GF(64) product vs 18-CCX GF(8) Karatsuba; ~38,100 CCX across DQI.)
 */
typedef struct {
    uint8_t  l_mask[QUA_ELE_GF64_RANK]; /* 6-bit linear form on `a` */
    uint8_t  r_mask[QUA_ELE_GF64_RANK]; /* 6-bit linear form on `b` */
    uint8_t  w_mask[QUA_ELE_GF64_RANK]; /* 6-bit output fanout mask on `out` */
    uint32_t order[QUA_ELE_GF64_RANK];  /* CNOT-transition schedule */
} qua_ele_gf64_rank15_t;

static inline uint8_t qua_gf4_mul(uint8_t a, uint8_t b) {
    /* GF(4) = GF(2)[w] / (w^2 + w + 1), modulus = 0x7 (binary 111) */
    uint8_t r = 0u;
    for (uint32_t i = 0u; i < 2u; ++i) {
        if ((b >> i) & 1u) r ^= a;
        a = (uint8_t)((a << 1u) ^ ((a & 2u) ? 0x7u : 0u));
    }
    return (uint8_t)(r & 3u);
}

static inline uint8_t qua_gf4_inv(uint8_t a) {
    return qua_gf4_mul(a, a);
}

static inline bool qua_ele_build_gf64_rank15(
    uint8_t (*mul_fn)(uint8_t, uint8_t),
    qua_ele_gf64_rank15_t *out
) {
    memset(out, 0, sizeof(*out));

    /* 1. Locate w in GF(64) satisfying w^2 + w + 1 = 0 */
    uint8_t w_elem = 0u;
    for (uint8_t cand = 2u; cand < 64u; ++cand) {
        uint8_t c2 = mul_fn(cand, cand);
        if ((c2 ^ cand ^ 1u) == 0u) {
            w_elem = cand;
            break;
        }
    }
    if (w_elem == 0u) return false;

    /* 2. Locate z in GF(64) satisfying z^3 + z + 1 = 0 such that {w^i * z^j} spans GF(64) */
    uint8_t w_pow[2] = {1u, w_elem};
    uint8_t inv_cols[6];
    uint8_t z5[5];
    bool found_basis = false;

    for (uint8_t z_cand = 2u; z_cand < 64u && !found_basis; ++z_cand) {
        uint8_t z2 = mul_fn(z_cand, z_cand);
        uint8_t z3 = mul_fn(z2, z_cand);
        if ((z3 ^ z_cand ^ 1u) != 0u) continue;

        uint8_t basis[6];
        uint8_t zp = 1u;
        for (uint32_t j = 0u; j < 3u; ++j) {
            for (uint32_t i = 0u; i < 2u; ++i) {
                basis[2u * j + i] = mul_fn(w_pow[i], zp);
            }
            zp = mul_fn(zp, z_cand);
        }

        uint8_t mat[6], inv[6];
        for (uint32_t r = 0u; r < 6u; ++r) {
            uint8_t row = 0u;
            for (uint32_t c = 0u; c < 6u; ++c) {
                if ((basis[c] >> r) & 1u) row |= (uint8_t)(1u << c);
            }
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
            for (uint32_t d = 1u; d < 5u; ++d) z5[d] = mul_fn(z5[d - 1u], z_cand);
            found_basis = true;
        }
    }
    if (!found_basis) return false;

    /* 3. Invert 5x5 Vandermonde over P^1(GF(4)) = {0, 1, 2(w), 3(w^2), inf} */
    static const uint8_t pts[4] = {0u, 1u, 2u, 3u};
    uint8_t vmat[5][5], vinv[5][5];
    memset(vmat, 0, sizeof(vmat));
    memset(vinv, 0, sizeof(vinv));
    for (uint32_t r = 0u; r < 4u; ++r) {
        uint8_t p_pow = 1u;
        for (uint32_t c = 0u; c < 5u; ++c) {
            vmat[r][c] = p_pow;
            p_pow = qua_gf4_mul(p_pow, pts[r]);
        }
        vinv[r][r] = 1u;
    }
    vmat[4][4] = 1u;
    vinv[4][4] = 1u;

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
        uint8_t scale = qua_gf4_inv(vmat[col][col]);
        for (uint32_t k = 0u; k < 5u; ++k) {
            vmat[col][k] = qua_gf4_mul(vmat[col][k], scale);
            vinv[col][k] = qua_gf4_mul(vinv[col][k], scale);
        }
        for (uint32_t r = 0u; r < 5u; ++r) {
            if (r != col && vmat[r][col] != 0u) {
                uint8_t f = vmat[r][col];
                for (uint32_t k = 0u; k < 5u; ++k) {
                    vmat[r][k] ^= qua_gf4_mul(f, vmat[col][k]);
                    vinv[r][k] ^= qua_gf4_mul(f, vinv[col][k]);
                }
            }
        }
    }

    /* 4. Emit 5 * 3 = 15 bilinear triples via 3-AND GF(4) Karatsuba */
    static const uint8_t k3_u[3]   = {0x1u, 0x2u, 0x3u};
    static const uint8_t k3_v[3]   = {0x1u, 0x2u, 0x3u};
    static const uint8_t k3_out[3] = {0x3u, 0x1u, 0x2u};

    for (uint32_t m = 0u; m < 5u; ++m) {
        uint8_t u_lin[2] = {0u, 0u};
        uint8_t v_lin[2] = {0u, 0u};
        if (m < 4u) {
            uint8_t p_pow = 1u;
            for (uint32_t j = 0u; j < 3u; ++j) {
                for (uint32_t i = 0u; i < 2u; ++i) {
                    uint8_t prod = qua_gf4_mul((uint8_t)(1u << i), p_pow);
                    for (uint32_t b = 0u; b < 2u; ++b) {
                        if ((prod >> b) & 1u) {
                            u_lin[b] ^= inv_cols[2u * j + i];
                            v_lin[b] ^= inv_cols[2u * j + i];
                        }
                    }
                }
                p_pow = qua_gf4_mul(p_pow, pts[m]);
            }
        } else {
            u_lin[0] = inv_cols[4];
            u_lin[1] = inv_cols[5];
            v_lin[0] = inv_cols[4];
            v_lin[1] = inv_cols[5];
        }

        for (uint32_t k = 0u; k < 3u; ++k) {
            uint32_t idx = 3u * m + k;
            uint8_t lm = 0u, rm = 0u, wm = 0u;
            for (uint32_t b = 0u; b < 2u; ++b) {
                if ((k3_u[k] >> b) & 1u) lm ^= u_lin[b];
                if ((k3_v[k] >> b) & 1u) rm ^= v_lin[b];
            }
            for (uint32_t d = 0u; d < 5u; ++d) {
                uint8_t coeff_gf4 = qua_gf4_mul(vinv[d][m], k3_out[k]);
                for (uint32_t i = 0u; i < 2u; ++i) {
                    if ((coeff_gf4 >> i) & 1u) {
                        wm ^= mul_fn(w_pow[i], z5[d]);
                    }
                }
            }
            out->l_mask[idx] = lm;
            out->r_mask[idx] = rm;
            out->w_mask[idx] = wm;
            out->order[idx]  = idx;
        }
    }

    /* 5. Exhaustive verification over all 64 x 64 = 4,096 pairs in GF(64) */
    for (uint32_t a = 0u; a < 64u; ++a) {
        for (uint32_t b = 0u; b < 64u; ++b) {
            uint8_t got = 0u;
            for (uint32_t p = 0u; p < QUA_ELE_GF64_RANK; ++p) {
                uint32_t la = (uint32_t)__builtin_parity((unsigned)(out->l_mask[p] & a));
                uint32_t rb = (uint32_t)__builtin_parity((unsigned)(out->r_mask[p] & b));
                if (la & rb) got ^= out->w_mask[p];
            }
            if (got != mul_fn((uint8_t)a, (uint8_t)b)) return false;
        }
    }

    /* 6. Greedy CNOT-transition ordering */
    bool used[QUA_ELE_GF64_RANK];
    memset(used, 0, sizeof(used));
    uint32_t cur = 0u;
    for (uint32_t step = 0u; step < QUA_ELE_GF64_RANK; ++step) {
        uint32_t best_p = 0u, best_cost = 9999u;
        for (uint32_t p = 0u; p < QUA_ELE_GF64_RANK; ++p) {
            if (used[p]) continue;
            uint32_t cost = (step == 0u)
                ? (uint32_t)(__builtin_popcount(out->l_mask[p]) + __builtin_popcount(out->r_mask[p]))
                : (uint32_t)(__builtin_popcount((unsigned)(out->l_mask[p] ^ out->l_mask[cur])) +
                             __builtin_popcount((unsigned)(out->r_mask[p] ^ out->r_mask[cur])));
            if (cost < best_cost) { best_cost = cost; best_p = p; }
        }
        used[best_p] = true;
        out->order[step] = best_p;
        cur = best_p;
    }
    return true;
}

static inline void qua_ele_emit_gf64_rank15_mul(
    qua_builder_t *b,
    const uint32_t qa[6],
    const uint32_t qb[6],
    const uint32_t qout[6],
    uint32_t anc_prod,
    uint32_t hmr_bit,
    const qua_ele_gf64_rank15_t *plan
) {
    for (uint32_t step = 0u; step < QUA_ELE_GF64_RANK; ++step) {
        uint32_t p = plan->order[step];
        uint8_t lm = plan->l_mask[p];
        uint8_t rm = plan->r_mask[p];
        uint8_t wm = plan->w_mask[p];
        if (lm == 0u || rm == 0u || wm == 0u) continue;

        uint32_t piv_a = (uint32_t)__builtin_ctz((unsigned)lm);
        uint32_t piv_b = (uint32_t)__builtin_ctz((unsigned)rm);

        for (uint32_t i = piv_a + 1u; i < 6u; ++i) {
            if ((lm >> i) & 1u) qua_cx(b, qa[i], qa[piv_a]);
        }
        for (uint32_t j = piv_b + 1u; j < 6u; ++j) {
            if ((rm >> j) & 1u) qua_cx(b, qb[j], qb[piv_b]);
        }

        qua_and_compute(b, qa[piv_a], qb[piv_b], anc_prod);
        for (uint32_t k = 0u; k < 6u; ++k) {
            if ((wm >> k) & 1u) qua_cx(b, anc_prod, qout[k]);
        }
        qua_and_uncompute(b, qa[piv_a], qb[piv_b], anc_prod, hmr_bit);

        for (uint32_t j = 5u; j > piv_b; --j) {
            if ((rm >> j) & 1u) qua_cx(b, qb[j], qb[piv_b]);
        }
        for (uint32_t i = 5u; i > piv_a; --i) {
            if ((lm >> i) & 1u) qua_cx(b, qa[i], qa[piv_a]);
        }
    }
}

/* ============================================================================
 * Pass 4: Unified .kmb Peephole Self-Inverse Gate Cancellation (All Benchmarks)
 * ============================================================================
 * Cancels adjacent identical self-inverse gates (X, CX, CZ, CCX, CSWAP) on
 * the same qubits without crossing conditions, measurements (HMR), or Z_POW.
 * Runs in O(num_ops) in-place on any `qua_builder_t`.
 */
static inline uint64_t qua_ele_optimize_builder(qua_builder_t *b) {
    if (!b || b->num_ops < 2u) return 0u;
    size_t n = b->num_ops;
    uint8_t *dead = (uint8_t *)calloc(n, 1u);
    if (!dead) return 0u;

    /* Map each op index to its payload index in 1q / 2q / 3q arrays */
    uint32_t *op_idx = (uint32_t *)malloc(n * sizeof(uint32_t));
    if (!op_idx) { free(dead); return 0u; }

    size_t i1 = 0u, i2 = 0u, i3 = 0u;
    for (size_t i = 0u; i < n; ++i) {
        uint8_t op = b->op_types[i];
        if ((op & QUA_OP_HAS_3Q) == QUA_OP_HAS_3Q) op_idx[i] = (uint32_t)(i3++);
        else if ((op & QUA_OP_HAS_2Q) == QUA_OP_HAS_2Q) op_idx[i] = (uint32_t)(i2++);
        else if ((op & QUA_OP_HAS_1Q) == QUA_OP_HAS_1Q) op_idx[i] = (uint32_t)(i1++);
        else op_idx[i] = 0u;
    }

    /* Track last active op index per qubit (-1 if none or invalidated) */
    uint32_t nq = b->num_qubits;
    int64_t *last_op = (int64_t *)malloc((size_t)(nq ? nq : 1u) * sizeof(int64_t));
    if (!last_op) { free(op_idx); free(dead); return 0u; }
    for (uint32_t q = 0u; q < nq; ++q) last_op[q] = -1;

    uint64_t cancelled = 0u;
    for (size_t i = 0u; i < n; ++i) {
        uint8_t op = b->op_types[i];
        /* Global condition/tag barriers invalidate all qubit history */
        if (op == QUA_OP_PUSH_CONDITION || op == QUA_OP_POP_CONDITION ||
            op == QUA_OP_INFO_TAG || op == QCC_OP_DEBUG_PRINT_EMPTY ||
            op == QCC_OP_DEBUG_PRINT_EMPTY_IF) {
            for (uint32_t q = 0u; q < nq; ++q) last_op[q] = -1;
            continue;
        }
        uint8_t qmask = op & QUA_OP_HAS_3Q;
        if (qmask == 0u) {
            /* Pure classical or global-phase 0-qubit op (NEG, BIT_STORE*, BIT_INVERT*, etc.) */
            continue;
        }
        if ((op & QUA_OP_HAS_COND) || op == QUA_OP_HMR || op == QUA_OP_Z_POW ||
            op == QUA_OP_R || op == QCC_OP_DEBUG_PRINT_Q) {
            uint32_t idx = op_idx[i];
            if (qmask == QUA_OP_HAS_3Q) {
                if (b->q3_0[idx] < nq) last_op[b->q3_0[idx]] = -1;
                if (b->q3_1[idx] < nq) last_op[b->q3_1[idx]] = -1;
                if (b->q3_2[idx] < nq) last_op[b->q3_2[idx]] = -1;
            } else if (qmask == QUA_OP_HAS_2Q) {
                if (b->q2_0[idx] < nq) last_op[b->q2_0[idx]] = -1;
                if (b->q2_1[idx] < nq) last_op[b->q2_1[idx]] = -1;
            } else if (qmask == QUA_OP_HAS_1Q) {
                if (b->q1_0[idx] < nq) last_op[b->q1_0[idx]] = -1;
            }
            continue;
        }

        if (op == QUA_OP_X || op == QUA_OP_Z) {
            uint32_t q = b->q1_0[op_idx[i]];
            int64_t prev = (q < nq) ? last_op[q] : -1;
            if (prev >= 0 && b->op_types[prev] == op && !dead[prev]) {
                dead[prev] = 1u;
                dead[i] = 1u;
                last_op[q] = -1;
                cancelled += 2u;
            } else if (q < nq) {
                last_op[q] = (int64_t)i;
            }
        } else if (op == QUA_OP_CX || op == QUA_OP_CZ || op == QUA_OP_SWAP) {
            uint32_t idx = op_idx[i];
            uint32_t q0 = b->q2_0[idx], q1 = b->q2_1[idx];
            int64_t p0 = (q0 < nq) ? last_op[q0] : -1;
            int64_t p1 = (q1 < nq) ? last_op[q1] : -1;
            if (p0 >= 0 && p0 == p1 && b->op_types[p0] == op && !dead[p0]) {
                uint32_t pidx = op_idx[p0];
                bool same = (b->q2_0[pidx] == q0 && b->q2_1[pidx] == q1) ||
                            ((op == QUA_OP_CZ || op == QUA_OP_SWAP) &&
                             b->q2_0[pidx] == q1 && b->q2_1[pidx] == q0);
                if (same) {
                    dead[p0] = 1u;
                    dead[i] = 1u;
                    last_op[q0] = -1;
                    last_op[q1] = -1;
                    cancelled += 2u;
                    continue;
                }
            }
            if (q0 < nq) last_op[q0] = (int64_t)i;
            if (q1 < nq) last_op[q1] = (int64_t)i;
        } else if (op == QUA_OP_CCX || op == QUA_OP_CCZ || op == QUA_OP_CSWAP) {
            uint32_t idx = op_idx[i];
            uint32_t q0 = b->q3_0[idx], q1 = b->q3_1[idx], q2 = b->q3_2[idx];
            int64_t p0 = (q0 < nq) ? last_op[q0] : -1;
            int64_t p1 = (q1 < nq) ? last_op[q1] : -1;
            int64_t p2 = (q2 < nq) ? last_op[q2] : -1;
            if (p0 >= 0 && p0 == p1 && p0 == p2 && b->op_types[p0] == op && !dead[p0]) {
                uint32_t pidx = op_idx[p0];
                bool same = false;
                if (op == QUA_OP_CCX) {
                    same = (b->q3_2[pidx] == q2) &&
                           ((b->q3_0[pidx] == q0 && b->q3_1[pidx] == q1) ||
                            (b->q3_0[pidx] == q1 && b->q3_1[pidx] == q0));
                } else if (op == QUA_OP_CSWAP) {
                    same = (b->q3_0[pidx] == q0) &&
                           ((b->q3_1[pidx] == q1 && b->q3_2[pidx] == q2) ||
                            (b->q3_1[pidx] == q2 && b->q3_2[pidx] == q1));
                } else { /* QUA_OP_CCZ: symmetric in all 3 qubits */
                    uint32_t a0 = b->q3_0[pidx], a1 = b->q3_1[pidx], a2 = b->q3_2[pidx];
                    same = (a0 == q0 && ((a1 == q1 && a2 == q2) || (a1 == q2 && a2 == q1))) ||
                           (a0 == q1 && ((a1 == q0 && a2 == q2) || (a1 == q2 && a2 == q0))) ||
                           (a0 == q2 && ((a1 == q0 && a2 == q1) || (a1 == q1 && a2 == q0)));
                }
                if (same) {
                    dead[p0] = 1u;
                    dead[i] = 1u;
                    last_op[q0] = -1;
                    last_op[q1] = -1;
                    last_op[q2] = -1;
                    cancelled += 2u;
                    if (b->count_toffoli >= 2u) b->count_toffoli -= 2u;
                    continue;
                }
            }
            if (q0 < nq) last_op[q0] = (int64_t)i;
            if (q1 < nq) last_op[q1] = (int64_t)i;
            if (q2 < nq) last_op[q2] = (int64_t)i;
        } else {
            uint32_t idx = op_idx[i];
            if (qmask == QUA_OP_HAS_3Q) {
                if (b->q3_0[idx] < nq) last_op[b->q3_0[idx]] = -1;
                if (b->q3_1[idx] < nq) last_op[b->q3_1[idx]] = -1;
                if (b->q3_2[idx] < nq) last_op[b->q3_2[idx]] = -1;
            } else if (qmask == QUA_OP_HAS_2Q) {
                if (b->q2_0[idx] < nq) last_op[b->q2_0[idx]] = -1;
                if (b->q2_1[idx] < nq) last_op[b->q2_1[idx]] = -1;
            } else if (qmask == QUA_OP_HAS_1Q) {
                if (b->q1_0[idx] < nq) last_op[b->q1_0[idx]] = -1;
            }
        }
    }

    if (cancelled > 0u) {
        size_t w_op = 0u, w1 = 0u, w2 = 0u, w3 = 0u, w1b = 0u, wc = 0u, wa = 0u;
        size_t r1 = 0u, r2 = 0u, r3 = 0u, r1b = 0u, rc = 0u, ra = 0u;
        for (size_t i = 0u; i < n; ++i) {
            uint8_t op = b->op_types[i];
            bool has_3q = (op & QUA_OP_HAS_3Q) == QUA_OP_HAS_3Q;
            bool has_2q = !has_3q && ((op & QUA_OP_HAS_2Q) == QUA_OP_HAS_2Q);
            bool has_1q = !has_3q && !has_2q && ((op & QUA_OP_HAS_1Q) == QUA_OP_HAS_1Q);
            bool has_1b = (op & QUA_OP_HAS_1B) == QUA_OP_HAS_1B;
            bool has_cd = (op & QUA_OP_HAS_COND) == QUA_OP_HAS_COND;
            bool has_ag = (op == QUA_OP_Z_POW || op == QUA_OP_Z_POW_IF);

            if (!dead[i]) {
                b->op_types[w_op++] = op;
                if (has_3q) { b->q3_0[w3] = b->q3_0[r3]; b->q3_1[w3] = b->q3_1[r3]; b->q3_2[w3] = b->q3_2[r3]; w3++; }
                if (has_2q) { b->q2_0[w2] = b->q2_0[r2]; b->q2_1[w2] = b->q2_1[r2]; w2++; }
                if (has_1q) { b->q1_0[w1] = b->q1_0[r1]; w1++; }
                if (has_1b) { b->b1_0[w1b] = b->b1_0[r1b]; w1b++; }
                if (has_cd) { b->bcond[wc] = b->bcond[rc]; wc++; }
                if (has_ag) { b->angles[wa] = b->angles[ra]; wa++; }
            }
            if (has_3q) r3++;
            if (has_2q) r2++;
            if (has_1q) r1++;
            if (has_1b) r1b++;
            if (has_cd) rc++;
            if (has_ag) ra++;
        }
        b->num_ops = w_op;
        b->num_3q = w3;
        b->num_2q = w2;
        b->num_1q = w1;
        b->num_1b = w1b;
        b->num_cond = wc;
        b->num_angles = wa;
    }

    free(last_op);
    free(op_idx);
    free(dead);
    return cancelled;
}

/* ============================================================================
 * Pass 5: Unified .kmb Ancilla Live-Range Compactor (All Benchmarks)
 * ============================================================================
 * Computes exact [first_op, last_op] intervals for every non-IO ancilla qubit
 * in `b` and packs disjoint lifetimes onto the minimum physical qubit bank.
 */
static inline uint32_t qua_ele_compact_qubits(qua_builder_t *b) {
    if (!b || b->num_qubits == 0u || b->num_ops == 0u) return 0u;
    uint32_t nq = b->num_qubits;
    uint8_t *pinned = (uint8_t *)calloc(nq, 1u);
    int64_t *first_op = (int64_t *)malloc((size_t)nq * sizeof(int64_t));
    int64_t *last_op  = (int64_t *)malloc((size_t)nq * sizeof(int64_t));
    uint32_t *remap   = (uint32_t *)malloc((size_t)nq * sizeof(uint32_t));
    if (!pinned || !first_op || !last_op || !remap) {
        free(pinned); free(first_op); free(last_op); free(remap);
        return 0u;
    }
    for (uint32_t q = 0u; q < nq; ++q) {
        first_op[q] = -1;
        last_op[q]  = -1;
        remap[q]    = q;
    }
    uint32_t max_pinned_plus1 = 0u;
    for (uint32_t r = 0u; r < b->num_regs; ++r) {
        for (uint32_t i = 0u; i < b->regs[r].len; ++i) {
            uint32_t t = b->regs[r].targets[i];
            if (t & 0x80000000u) {
                uint32_t q = t & 0x7FFFFFFFu;
                if (q < nq) {
                    pinned[q] = 1u;
                    if (q + 1u > max_pinned_plus1) max_pinned_plus1 = q + 1u;
                }
            }
        }
    }

    size_t i1 = 0u, i2 = 0u, i3 = 0u;
    for (size_t i = 0u; i < b->num_ops; ++i) {
        uint8_t op = b->op_types[i];
        bool has_3q = (op & QUA_OP_HAS_3Q) == QUA_OP_HAS_3Q;
        bool has_2q = !has_3q && ((op & QUA_OP_HAS_2Q) == QUA_OP_HAS_2Q);
        bool has_1q = !has_3q && !has_2q && ((op & QUA_OP_HAS_1Q) == QUA_OP_HAS_1Q);
        if (has_3q) {
            uint32_t qs[3] = {b->q3_0[i3], b->q3_1[i3], b->q3_2[i3]};
            for (int k = 0; k < 3; ++k) {
                uint32_t q = qs[k];
                if (q < nq && !pinned[q]) {
                    if (first_op[q] < 0) first_op[q] = (int64_t)i;
                    last_op[q] = (int64_t)i;
                }
            }
            i3++;
        } else if (has_2q) {
            uint32_t qs[2] = {b->q2_0[i2], b->q2_1[i2]};
            for (int k = 0; k < 2; ++k) {
                uint32_t q = qs[k];
                if (q < nq && !pinned[q]) {
                    if (first_op[q] < 0) first_op[q] = (int64_t)i;
                    last_op[q] = (int64_t)i;
                }
            }
            i2++;
        } else if (has_1q) {
            uint32_t q = b->q1_0[i1];
            if (q < nq && !pinned[q]) {
                if (first_op[q] < 0) first_op[q] = (int64_t)i;
                last_op[q] = (int64_t)i;
            }
            i1++;
        }
    }

    /* Collect active ancillas sorted by first_op (naturally ordered by q or birth) */
    uint32_t *anc_list = (uint32_t *)malloc((size_t)nq * sizeof(uint32_t));
    int64_t *slot_end  = (int64_t *)malloc((size_t)nq * sizeof(int64_t));
    if (!anc_list || !slot_end) {
        free(anc_list); free(slot_end);
        free(pinned); free(first_op); free(last_op); free(remap);
        return 0u;
    }
    uint32_t n_anc = 0u;
    for (uint32_t q = 0u; q < nq; ++q) {
        if (!pinned[q] && first_op[q] >= 0) anc_list[n_anc++] = q;
    }
    /* Insertion/quick sort by first_op ascending */
    for (uint32_t i = 1u; i < n_anc; ++i) {
        uint32_t key = anc_list[i];
        int64_t kf = first_op[key];
        uint32_t j = i;
        while (j > 0u && first_op[anc_list[j - 1u]] > kf) {
            anc_list[j] = anc_list[j - 1u];
            --j;
        }
        anc_list[j] = key;
    }

    uint32_t n_slots = 0u;
    for (uint32_t i = 0u; i < n_anc; ++i) {
        uint32_t q = anc_list[i];
        int64_t f = first_op[q], l = last_op[q];
        uint32_t chosen = UINT32_MAX;
        for (uint32_t s = 0u; s < n_slots; ++s) {
            if (slot_end[s] < f) {
                chosen = s;
                break;
            }
        }
        if (chosen == UINT32_MAX) {
            chosen = n_slots++;
        }
        slot_end[chosen] = l;
        remap[q] = max_pinned_plus1 + chosen;
    }

    uint32_t new_nq = max_pinned_plus1 + n_slots;
    if (new_nq < nq) {
        for (size_t i = 0u; i < b->num_1q; ++i) b->q1_0[i] = remap[b->q1_0[i]];
        for (size_t i = 0u; i < b->num_2q; ++i) {
            b->q2_0[i] = remap[b->q2_0[i]];
            b->q2_1[i] = remap[b->q2_1[i]];
        }
        for (size_t i = 0u; i < b->num_3q; ++i) {
            b->q3_0[i] = remap[b->q3_0[i]];
            b->q3_1[i] = remap[b->q3_1[i]];
            b->q3_2[i] = remap[b->q3_2[i]];
        }
        b->num_qubits = new_nq;
    }

    uint32_t saved = (new_nq < nq) ? (nq - new_nq) : 0u;
    free(anc_list); free(slot_end);
    free(pinned); free(first_op); free(last_op); free(remap);
    return saved;
}

#endif /* QUA_ELE_H */
