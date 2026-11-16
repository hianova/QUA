/*
 * qua_bitop.h — Grand Unified Reversible Bit-Operator DAG & Bilevel SEARU x FLOW 2-Jet Search Engine
 *
 * Bottom Layer (Problem-Agnostic):
 *   All quantum circuits decompose into a unified reversible bit-operator DAG
 *   (IO, XOR/CX, AND/CCX, CSWAP/MUX, QROM, MUL, PHASE) solved by ONE bilevel
 *   Meta-SEARU + FLOW 2-Jet engine optimizing:
 *     1. Operator topology & decomposition
 *     2. Intermediate register bit-widths & precision/error budget
 *     3. Pebbling / uncomputation strategy (HMR vs Adjoint) + heterogeneous-width bank packing
 *     4. Post-emission .kmb peephole cancellation + global live-range qubit compaction
 *
 * Top Layer (Spec Only):
 *   Each competition task provides ONLY:
 *     - Mathematical truth / error evaluator (`eval_spec`)
 *     - Cost Law (`Q*C`, `R(eps)*Q^2*T`, `lambda*Q*T`) + contract limits (`allow_z_pow`, `max_qubits`)
 */
#ifndef QUA_BITOP_H
#define QUA_BITOP_H

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "flow_opt_params.h"
#include "flow.h"
#include "searu.h"
#include "qua_kmb.h"
#include "qua_ele.h"

#define QUA_BITOP_MAX_NODES QUA_ELE_MAX_NODES

typedef enum {
    QUA_BOP_IO    = 0,
    QUA_BOP_XOR   = 1,
    QUA_BOP_AND   = 2,
    QUA_BOP_CSWAP = 3,
    QUA_BOP_QROM  = 4,
    QUA_BOP_MUL   = 5,
    QUA_BOP_PHASE = 6
} qua_bop_kind_t;

/* Explicit Quantum-Algebraic Representation Type for Composable Shortcut Routing:
 *   - KET:         Materialized computational-basis register |x>|f(x)>
 *   - LINF2:       In-place / zero-Toffoli GL_w(F_2) linear basis bijection (0 CCX)
 *   - PHASE:       Diagonal phase kickback exp(2*pi*i*phi(x))|x> (0 workspace qubits, additive)
 *   - PARITY_DUAL: Post-HMR classical mask m carrying residual sign (-1)^{<m, f(x)>}|x> */
typedef enum {
    QUA_REP_KET         = 0,
    QUA_REP_LINF2       = 1,
    QUA_REP_PHASE       = 2,
    QUA_REP_PARITY_DUAL = 3
} qua_rep_kind_t;

typedef struct {
    qua_bop_kind_t kind;
    qua_rep_kind_t rep;
    uint32_t       width_bits;
    uint32_t       tof_compute;
    uint32_t       tof_uncompute; /* 0 if HMR+CZ_IF */
    uint32_t       z_pow_count;
    uint8_t        is_hmr_cleanable;
    uint8_t        is_output;
    uint64_t       deps_mask;
} qua_bitop_node_t;

typedef struct {
    uint32_t           num_nodes;
    uint32_t           base_io_qubits;
    qua_ele_cost_law_t cost_law;
    double             q_violation;   /* >= 0: soft excess qubits over contract/budget */
    double             err_violation; /* >= 0: soft excess normalized error max(0, err/tol - 1) */
    qua_bitop_node_t   nodes[QUA_BITOP_MAX_NODES];
} qua_bitop_dag_t;

typedef struct {
    const char        *name;
    qua_ele_cost_law_t cost_law;
    bool               allow_z_pow;
    uint32_t       max_qubits;   /* 0 => unlimited */
    uint32_t       genome_bits;  /* <= 256 */
    double         energy_scale; /* scales raw_cost into ~2048..8192 SA temperature regime */
    uint32_t       meta_iters;
    uint32_t       inner_iters;
    uint32_t       restarts;
    bool (*eval_spec)(
        const uint64_t *genome,
        qua_bitop_dag_t *dag,
        double *out_raw_cost,
        uint64_t *out_frust,
        void *user_ctx
    );
    void *user_ctx;
    bool (*relax_2jet)(
        const uint64_t *genome,
        qua_bitop_dag_t *dag,
        double *inout_raw_cost,
        uint64_t *inout_frust,
        void *user_ctx
    );
} qua_bitop_spec_t;

static inline qua_rep_kind_t qua_bitop_default_rep(qua_bop_kind_t kind) {
    if (kind == QUA_BOP_XOR)   return QUA_REP_LINF2;
    if (kind == QUA_BOP_PHASE) return QUA_REP_PHASE;
    return QUA_REP_KET;
}

static inline uint32_t qua_bitop_dag_push_rep(
    qua_bitop_dag_t *dag,
    qua_bop_kind_t kind,
    qua_rep_kind_t rep,
    uint32_t width_bits,
    uint32_t tof_compute,
    uint32_t tof_uncompute,
    uint32_t z_pow_count,
    bool is_hmr_cleanable,
    bool is_output,
    uint64_t deps_mask
) {
    if (dag->num_nodes >= QUA_BITOP_MAX_NODES) return dag->num_nodes;
    uint32_t idx = dag->num_nodes++;
    dag->nodes[idx].kind             = kind;
    dag->nodes[idx].rep              = rep;
    dag->nodes[idx].width_bits       = width_bits;
    dag->nodes[idx].tof_compute      = tof_compute;
    dag->nodes[idx].tof_uncompute    = tof_uncompute;
    dag->nodes[idx].z_pow_count      = z_pow_count;
    dag->nodes[idx].is_hmr_cleanable = is_hmr_cleanable ? 1u : 0u;
    dag->nodes[idx].is_output        = is_output ? 1u : 0u;
    dag->nodes[idx].deps_mask        = deps_mask;
    return idx;
}

static inline uint32_t qua_bitop_dag_push(
    qua_bitop_dag_t *dag,
    qua_bop_kind_t kind,
    uint32_t width_bits,
    uint32_t tof_compute,
    uint32_t tof_uncompute,
    uint32_t z_pow_count,
    bool is_hmr_cleanable,
    bool is_output,
    uint64_t deps_mask
) {
    return qua_bitop_dag_push_rep(dag, kind, qua_bitop_default_rep(kind), width_bits,
                                  tof_compute, tof_uncompute, z_pow_count,
                                  is_hmr_cleanable, is_output, deps_mask);
}

static inline void qua_bitop_to_ele_dag(const qua_bitop_dag_t *bdag, qua_ele_dag_t *edag) {
    memset(edag, 0, sizeof(*edag));
    edag->num_nodes      = bdag->num_nodes;
    edag->base_io_qubits = bdag->base_io_qubits;
    edag->cost_law       = bdag->cost_law;
    for (uint32_t i = 0u; i < bdag->num_nodes && i < QUA_ELE_MAX_NODES; ++i) {
        edag->nodes[i].width_qubits       = (uint16_t)bdag->nodes[i].width_bits;
        edag->nodes[i].compute_toffolis   = bdag->nodes[i].tof_compute;
        edag->nodes[i].uncompute_toffolis = bdag->nodes[i].tof_uncompute;
        edag->nodes[i].can_hmr_erase      = bdag->nodes[i].is_hmr_cleanable;
        edag->nodes[i].is_persistent_io   = bdag->nodes[i].is_output;
        edag->nodes[i].fanin_mask         = bdag->nodes[i].deps_mask;
    }
}

static inline void qua_bitop_eval_dag_pebble(
    const qua_bitop_dag_t *bdag,
    uint64_t eager_mask,
    qua_ele_pebble_plan_t *out_plan,
    qua_ele_reg_alloc_t *out_alloc,
    uint64_t seed
) {
    qua_ele_dag_t edag;
    qua_bitop_to_ele_dag(bdag, &edag);
    qua_ele_pebble_plan_t plan;
    qua_ele_eval_pebble_plan(&edag, eager_mask, &plan);
    if (out_plan) *out_plan = plan;
    if (out_alloc) *out_alloc = qua_ele_alloc_qubits(&edag, &plan, seed);
}

typedef struct {
    const qua_bitop_spec_t *spec;
    uint32_t                eval_iters;
} qua_bitop_bilevel_ctx_t;

/*
 * 4-Stage Homotopy Penalty Annealing + Symbiotic Bilevel (SEARU x FLOW 2-Jet)
 *
 * Stages beta_s in {0.25, 1.0, 4.0, 32.0}:
 *   - Early stages allow discrete SEARU moves to tunnel across soft qubit-budget
 *     and precision-boundary barriers instead of bouncing off flat 10^12 walls.
 *   - Promising discrete topologies invoke inner `relax_2jet` fiber optimization
 *     during search (SymbioticDuality.lean Theorem 2).
 *   - Strictly feasible solutions (`q_violation <= 0 && err_violation <= 0`) are
 *     tracked separately in `best_feas_genome` so penalized infeasible states can
 *     never overwrite a valid circuit.
 */
static inline uint64_t qua_bitop_splitmix64(uint64_t *state) {
    uint64_t z = (*state += UINT64_C(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30u)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27u)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31u);
}

static inline uint64_t qua_bitop_run_inner(
    const qua_bitop_spec_t *spec,
    const searu_hyperparams *hp,
    uint64_t seed,
    uint32_t iters,
    uint64_t out_genome[4],
    double *out_best_raw
) {
    static const double BETA_SCHED[4] = {0.25, 1.0, 4.0, 32.0};
    const double escale = (spec->energy_scale > 0.0) ? spec->energy_scale : 1.0;
    const uint64_t full_mask = (spec->genome_bits >= 64u)
        ? UINT64_MAX
        : ((UINT64_C(1) << spec->genome_bits) - 1u);

    uint64_t s_init = seed ^ UINT64_C(0x6A09E667F3BCC909);
    uint64_t cur_init[4] = {
        qua_bitop_splitmix64(&s_init) & full_mask,
        qua_bitop_splitmix64(&s_init),
        qua_bitop_splitmix64(&s_init),
        qua_bitop_splitmix64(&s_init)
    };
    uint64_t best_feas_g[4] = {0u, 0u, 0u, 0u};
    double   best_feas_raw = 1e30;
    bool     have_feas = false;
    uint64_t last_best_energy = UINT64_C(1000000000000);

    uint32_t n_stages = (iters >= 64u) ? 4u : 1u;
    uint32_t stage_iters = (iters + n_stages - 1u) / n_stages;
    if (stage_iters < 16u) stage_iters = 16u;

    for (uint32_t stg = 0u; stg < n_stages; ++stg) {
        const double beta = (n_stages == 1u) ? 32.0 : BETA_SCHED[stg];
        searu_hyperparams shp = *hp;
        if (stg > 0u) {
            /* Gentle stage-coupled initial temperature so tunneling cools smoothly */
            shp.inl_init_temp_q12 = (hp->inl_init_temp_q12 >> stg);
            if (shp.inl_init_temp_q12 < (UINT64_C(16) << 12u)) {
                shp.inl_init_temp_q12 = (UINT64_C(16) << 12u);
            }
        }
        searu_result res;
        uint64_t stg_seed = seed ^ (UINT64_C(0xD1B54A32D192ED03) * (stg + 1u));
        SEARU_SOLVE_INLINE_HP(&res, spec->genome_bits, cur_init, NULL, stg_seed, stage_iters, &shp,
                              cand, frust, energy, {
            qua_bitop_dag_t dag;
            memset(&dag, 0, sizeof(dag));
            dag.cost_law = spec->cost_law;
            double raw_cost = 1e30;
            uint64_t f_words[4] = {full_mask, UINT64_MAX, UINT64_MAX, UINT64_MAX};
            bool ok = spec->eval_spec(cand, &dag, &raw_cost, f_words, spec->user_ctx);

            /* Inner FLOW 2-Jet continuous/quantized fiber relaxation on promising candidates */
            if (spec->relax_2jet != NULL && isfinite(raw_cost) && raw_cost > 0.0 &&
                dag.q_violation <= 2.0 && dag.err_violation <= 0.35) {
                bool ok2 = spec->relax_2jet(cand, &dag, &raw_cost, f_words, spec->user_ctx);
                if (ok2) ok = true;
            }

            if (f_words[0] == 0u && spec->genome_bits <= 64u) f_words[0] = full_mask;
            frust[0] = f_words[0];
            frust[1] = f_words[1];
            frust[2] = f_words[2];
            frust[3] = f_words[3];

            bool contract_ok = true;
            if (dag.num_nodes > 0u && !spec->allow_z_pow) {
                uint64_t z_total = 0u;
                for (uint32_t i = 0u; i < dag.num_nodes; ++i) z_total += dag.nodes[i].z_pow_count;
                if (z_total > 0u) contract_ok = false;
            }
            if (spec->max_qubits > 0u && dag.num_nodes > 0u && dag.q_violation == 0.0) {
                qua_ele_pebble_plan_t pp;
                qua_bitop_eval_dag_pebble(&dag, UINT64_MAX, &pp, NULL, 0u);
                if (pp.peak_qubits > spec->max_qubits) {
                    dag.q_violation = (double)(pp.peak_qubits - spec->max_qubits);
                }
            }

            double qv = (dag.q_violation > 0.0) ? dag.q_violation : 0.0;
            double ev = (dag.err_violation > 0.0) ? dag.err_violation : 0.0;
            bool strict_feas = ok && contract_ok && isfinite(raw_cost) && (raw_cost > 0.0) &&
                               (qv <= 0.0) && (ev <= 0.0);
            if (strict_feas && raw_cost < best_feas_raw) {
                best_feas_raw = raw_cost;
                for (int k = 0; k < 4; ++k) best_feas_g[k] = cand[k];
                have_feas = true;
            }

            if (!contract_ok || !isfinite(raw_cost) || raw_cost <= 0.0) {
                energy = UINT64_C(1000000000);
            } else {
                double base_e = raw_cost / escale;
                double ev_log = log2(1.0 + ev);
                double ev_lin = (ev < 2.0) ? ev : (2.0 + 0.25 * ev_log);
                double pen_e  = beta * (64.0 * qv + 8.0 * qv * qv +
                                        512.0 * ev_lin + 384.0 * ev_log * ev_log);
                double tot_e  = base_e + pen_e;
                if (tot_e < 1.0) tot_e = 1.0;
                if (tot_e > 1e11) tot_e = 1e11;
                energy = (uint64_t)llround(tot_e);
            }
        });
        for (int k = 0; k < 4; ++k) cur_init[k] = res.canvas[k];
        last_best_energy = res.best_energy;
    }

    if (have_feas) {
        if (out_genome) {
            for (int i = 0; i < 4; ++i) out_genome[i] = best_feas_g[i];
        }
        if (out_best_raw) *out_best_raw = best_feas_raw;
        double e_feas = best_feas_raw / escale;
        if (e_feas < 1.0) e_feas = 1.0;
        if (e_feas > 1e11) e_feas = 1e11;
        return (uint64_t)llround(e_feas);
    }

    if (out_genome) {
        for (int i = 0; i < 4; ++i) out_genome[i] = cur_init[i];
    }
    if (out_best_raw) *out_best_raw = 1e30;
    return last_best_energy + UINT64_C(1000000);
}

static inline uint64_t qua_bitop_bilevel_eval(const searu_hyperparams *hp, uint64_t seed, void *ctx_p) {
    const qua_bitop_bilevel_ctx_t *bctx = (const qua_bitop_bilevel_ctx_t *)ctx_p;
    return qua_bitop_run_inner(bctx->spec, hp, seed, bctx->eval_iters, NULL, NULL);
}

static inline void qua_bitop_solve(
    const qua_bitop_spec_t *spec,
    uint64_t base_seed,
    uint64_t best_genome[4],
    double *out_best_raw
) {
    qua_bitop_bilevel_ctx_t bctx = {spec, spec->inner_iters > 0u ? (spec->inner_iters / 4u) : 256u};
    if (bctx.eval_iters < 128u) bctx.eval_iters = 128u;
    qua_ele_self_opt_stats_t sopt_stats;
    searu_hyperparams opt_hp = qua_ele_self_opt_hp(
        qua_bitop_bilevel_eval, &bctx, base_seed,
        spec->meta_iters > 0u ? spec->meta_iters : 64u, &sopt_stats
    );
    printf("[%s Unified Bilevel SEARU x 2-Jet] train=%llu->%llu | val=%llu->%llu (accepted=%d, cool_q16=%llu)\n",
           spec->name,
           (unsigned long long)sopt_stats.base_train_cost, (unsigned long long)sopt_stats.opt_train_cost,
           (unsigned long long)sopt_stats.base_val_cost, (unsigned long long)sopt_stats.opt_val_cost,
           (int)sopt_stats.accepted, (unsigned long long)opt_hp.inl_cool_q16);

    double best_raw = 1e30;
    uint64_t win_g[4] = {0u, 0u, 0u, 0u};
    uint32_t restarts = spec->restarts > 0u ? spec->restarts : 4u;
    for (uint32_t r = 0u; r < restarts; ++r) {
        uint64_t g[4] = {0u, 0u, 0u, 0u};
        double raw = 1e30;
        qua_bitop_run_inner(spec, &opt_hp, base_seed ^ (UINT64_C(0x9E3779B97F4A7C15) * (r + 1u)),
                            spec->inner_iters > 0u ? spec->inner_iters : 2048u, g, &raw);
        if (raw < best_raw) {
            best_raw = raw;
            memcpy(win_g, g, sizeof(win_g));
        }
    }
    memcpy(best_genome, win_g, sizeof(win_g));
    if (out_best_raw) *out_best_raw = best_raw;
}

/* Unified post-emission .kmb optimization + contract enforcement + write */
static inline int qua_bitop_finalize_kmb(
    qua_builder_t *kb,
    const qua_bitop_spec_t *spec,
    const char *out_path
) {
    uint64_t cancelled = qua_ele_optimize_builder(kb);
    uint32_t q_before  = kb->num_qubits;
    uint32_t q_saved   = qua_ele_compact_qubits(kb);

    if (!spec->allow_z_pow && kb->count_z_pow > 0u) {
        fprintf(stderr, "[CONTRACT VIOLATION] %s forbids Z_POW, got %llu\n",
                spec->name, (unsigned long long)kb->count_z_pow);
        return 1;
    }
    if (spec->max_qubits > 0u && kb->num_qubits > spec->max_qubits) {
        fprintf(stderr, "[CONTRACT VIOLATION] %s requires Q <= %u, got %u\n",
                spec->name, spec->max_qubits, kb->num_qubits);
        return 1;
    }
    int rc = qua_builder_write_kmb(kb, out_path);
    printf("[%s Unified Finalize] %s | Q=%u (compacted %u->%u, saved=%u), CCX=%llu, HMR=%llu, Z_POW=%llu, peephole=-%llu\n",
           spec->name, out_path, kb->num_qubits, q_before, kb->num_qubits, q_saved,
           (unsigned long long)kb->count_toffoli, (unsigned long long)kb->count_hmr,
           (unsigned long long)kb->count_z_pow, (unsigned long long)cancelled);
    return rc;
}

/* ============================================================================
 * Unified Select-Swap QROM + HMR Phase-Erase (b_sh=0: Segment-Tree, b_sh>=1: Select-Swap)
 * ============================================================================ */

typedef struct {
    qua_builder_t  *kb;
    const uint64_t *table;
    uint32_t        num_rows;
    uint32_t        stride_u64;
    uint32_t        w_bits;
    const uint32_t *addr_msb;
    uint32_t        hi_bits;
    const uint32_t *ancillas;
    const uint32_t *slot_qs;
    uint32_t        b_ss;
    uint32_t       *next_cbit;
    bool            phase_fixup_mode;
    const uint32_t *m_bits;
} qua_bitop_qrom_ctx_t;

static inline bool qua_bitop_qrom_bit(
    const qua_bitop_qrom_ctx_t *ctx,
    uint32_t row,
    uint32_t bit_idx
) {
    if (row >= ctx->num_rows) return false;
    uint64_t w = ctx->table[(size_t)row * ctx->stride_u64 + (bit_idx >> 6u)];
    return ((w >> (bit_idx & 63u)) & 1u) != 0u;
}

static inline bool qua_bitop_qrom_block_const(
    const qua_bitop_qrom_ctx_t *ctx,
    uint32_t lo_blk,
    uint32_t hi_blk
) {
    if (lo_blk + 1u >= hi_blk) return true;
    uint32_t end_row = hi_blk * ctx->b_ss;
    if (end_row > ctx->num_rows) end_row = ctx->num_rows;
    uint32_t base_row = lo_blk * ctx->b_ss;
    if (base_row >= end_row) return true;
    for (uint32_t b = lo_blk + 1u; b < hi_blk; ++b) {
        uint32_t r0 = b * ctx->b_ss;
        if (r0 >= end_row) break;
        for (uint32_t s = 0u; s < ctx->b_ss; ++s) {
            uint32_t r_cur = r0 + s;
            uint32_t r_ref = base_row + s;
            if (r_cur >= end_row) break;
            for (uint32_t w = 0u; w < ctx->stride_u64; ++w) {
                uint64_t v_cur = ctx->table[(size_t)r_cur * ctx->stride_u64 + w];
                uint64_t v_ref = ctx->table[(size_t)r_ref * ctx->stride_u64 + w];
                if (v_cur != v_ref) return false;
            }
        }
    }
    return true;
}

static inline void qua_bitop_qrom_leaf(
    const qua_bitop_qrom_ctx_t *ctx,
    uint32_t blk,
    uint32_t ctrl
) {
    for (uint32_t s = 0u; s < ctx->b_ss; ++s) {
        uint32_t row = blk * ctx->b_ss + s;
        if (row >= ctx->num_rows) break;
        for (uint32_t j = 0u; j < ctx->w_bits; ++j) {
            if (qua_bitop_qrom_bit(ctx, row, j)) {
                uint32_t idx = s * ctx->w_bits + j;
                if (ctx->phase_fixup_mode) {
                    if (ctrl != UINT32_MAX) {
                        qua_z_if(ctx->kb, ctrl, ctx->m_bits[idx]);
                    }
                } else {
                    uint32_t q = ctx->slot_qs[idx];
                    if (ctrl == UINT32_MAX) qua_x(ctx->kb, q);
                    else qua_cx(ctx->kb, ctrl, q);
                }
            }
        }
    }
}

static inline void qua_bitop_qrom_rec(
    const qua_bitop_qrom_ctx_t *ctx,
    uint32_t ctrl,
    uint32_t level,
    uint32_t lo_blk,
    uint32_t hi_blk
) {
    uint32_t n_blocks = (ctx->num_rows + ctx->b_ss - 1u) / ctx->b_ss;
    if (lo_blk >= n_blocks) return;
    if (hi_blk <= n_blocks && (hi_blk - lo_blk == 1u || qua_bitop_qrom_block_const(ctx, lo_blk, hi_blk))) {
        qua_bitop_qrom_leaf(ctx, lo_blk, ctrl);
        return;
    }
    if (level >= ctx->hi_bits) {
        qua_bitop_qrom_leaf(ctx, lo_blk, ctrl);
        return;
    }
    uint32_t mid = (lo_blk + hi_blk) / 2u;
    if (n_blocks <= mid) {
        qua_bitop_qrom_rec(ctx, ctrl, level + 1u, lo_blk, mid);
        return;
    }
    uint32_t bit_q = ctx->addr_msb[level];
    if (ctrl == UINT32_MAX) {
        qua_x(ctx->kb, bit_q);
        qua_bitop_qrom_rec(ctx, bit_q, level + 1u, lo_blk, mid);
        qua_x(ctx->kb, bit_q);
        qua_bitop_qrom_rec(ctx, bit_q, level + 1u, mid, hi_blk);
    } else {
        uint32_t child = ctx->ancillas[level - 1u];
        qua_x(ctx->kb, bit_q);
        qua_and_compute(ctx->kb, ctrl, bit_q, child);
        qua_x(ctx->kb, bit_q);
        qua_bitop_qrom_rec(ctx, child, level + 1u, lo_blk, mid);
        qua_cx(ctx->kb, ctrl, child);
        qua_bitop_qrom_rec(ctx, child, level + 1u, mid, hi_blk);
        qua_and_uncompute(ctx->kb, ctrl, bit_q, child, (*ctx->next_cbit)++);
    }
}

static inline void qua_bitop_emit_qrom_ss(
    qua_builder_t *kb,
    const uint64_t *table,
    uint32_t num_rows,
    uint32_t stride_u64,
    const uint32_t *addr_lsb,
    uint32_t addr_bits,
    const uint32_t *target_qs,
    uint32_t w_bits,
    uint32_t b_sh,
    const uint32_t *anc_pool,
    uint32_t *next_cbit,
    bool inverse_hmr_erase
) {
    if (b_sh > addr_bits) b_sh = addr_bits;
    uint32_t b_ss = 1u << b_sh;
    uint32_t hi_bits = addr_bits - b_sh;
    uint32_t addr_msb[32];
    for (uint32_t i = 0u; i < hi_bits; ++i) {
        addr_msb[i] = addr_lsb[addr_bits - 1u - i];
    }

    if (b_sh == 0u) {
        if (!inverse_hmr_erase) {
            qua_bitop_qrom_ctx_t ctx = {
                kb, table, num_rows, stride_u64, w_bits, addr_msb, hi_bits,
                anc_pool, target_qs, 1u, next_cbit, false, NULL
            };
            qua_bitop_qrom_rec(&ctx, UINT32_MAX, 0u, 0u, 1u << hi_bits);
        } else {
            uint32_t *m_bits = (uint32_t *)malloc((size_t)w_bits * sizeof(uint32_t));
            if (!m_bits) abort();
            for (uint32_t j = 0u; j < w_bits; ++j) {
                m_bits[j] = (*next_cbit)++;
                qua_hmr(kb, target_qs[j], m_bits[j]);
            }
            qua_bitop_qrom_ctx_t ctx = {
                kb, table, num_rows, stride_u64, w_bits, addr_msb, hi_bits,
                anc_pool, target_qs, 1u, next_cbit, true, m_bits
            };
            qua_bitop_qrom_rec(&ctx, UINT32_MAX, 0u, 0u, 1u << hi_bits);
            free(m_bits);
        }
        return;
    }

    uint32_t slot_total = b_ss * w_bits;
    const uint32_t *slots = anc_pool;
    const uint32_t *unary_anc = anc_pool + slot_total;
    uint32_t *m_tgt = NULL;
    if (inverse_hmr_erase) {
        m_tgt = (uint32_t *)malloc((size_t)w_bits * sizeof(uint32_t));
        if (!m_tgt) abort();
        for (uint32_t j = 0u; j < w_bits; ++j) {
            m_tgt[j] = (*next_cbit)++;
            qua_hmr(kb, target_qs[j], m_tgt[j]);
        }
    }

    /* 1. Forward Unary SELECT over hi_bits into slots[0..b_ss-1] */
    qua_bitop_qrom_ctx_t fwd_ctx = {
        kb, table, num_rows, stride_u64, w_bits, addr_msb, hi_bits,
        unary_anc, slots, b_ss, next_cbit, false, NULL
    };
    qua_bitop_qrom_rec(&fwd_ctx, UINT32_MAX, 0u, 0u, 1u << hi_bits);

    /* 2. Forward SWAP tree over addr_lsb[0..b_sh-1] */
    for (uint32_t p = 0u; p < b_sh; ++p) {
        uint32_t stride = 1u << p;
        for (uint32_t s = 0u; s < b_ss; s += (2u << p)) {
            for (uint32_t j = 0u; j < w_bits; ++j) {
                qua_cswap(kb, addr_lsb[p], slots[s * w_bits + j], slots[(s + stride) * w_bits + j]);
            }
        }
    }

    /* 3. Extract or cancel target_qs */
    if (!inverse_hmr_erase) {
        for (uint32_t j = 0u; j < w_bits; ++j) qua_cx(kb, slots[j], target_qs[j]);
    } else {
        for (uint32_t j = 0u; j < w_bits; ++j) qua_z_if(kb, slots[j], m_tgt[j]);
        free(m_tgt);
    }

    /* 4. Reverse SWAP tree */
    for (int32_t p = (int32_t)b_sh - 1; p >= 0; --p) {
        uint32_t stride = 1u << (uint32_t)p;
        for (uint32_t s = 0u; s < b_ss; s += (2u << (uint32_t)p)) {
            for (uint32_t j = 0u; j < w_bits; ++j) {
                qua_cswap(kb, addr_lsb[p], slots[s * w_bits + j], slots[(s + stride) * w_bits + j]);
            }
        }
    }

    /* 5. HMR phase-erase on slots */
    uint32_t *m_slots = (uint32_t *)malloc((size_t)slot_total * sizeof(uint32_t));
    if (!m_slots) abort();
    for (uint32_t idx = 0u; idx < slot_total; ++idx) {
        m_slots[idx] = (*next_cbit)++;
        qua_hmr(kb, slots[idx], m_slots[idx]);
    }
    qua_bitop_qrom_ctx_t fix_ctx = {
        kb, table, num_rows, stride_u64, w_bits, addr_msb, hi_bits,
        unary_anc, slots, b_ss, next_cbit, true, m_slots
    };
    qua_bitop_qrom_rec(&fix_ctx, UINT32_MAX, 0u, 0u, 1u << hi_bits);
    free(m_slots);
}

/* ============================================================================
 * Stride-Indexed In-Place One-Hot Decoder + 0-Toffoli MBUC Uncomputation
 *   - Compute:   max(0, 2^n - 2) CCX, 0 extra ancillas beyond `leaves[0..2^n-1]`.
 *   - Uncompute: 0 CCX via `CX + HMR + CZ_IF`!
 * ============================================================================ */
static inline void qua_bitop_onehot_compute(
    qua_builder_t *kb,
    const uint32_t *addr_msb,
    uint32_t n_bits,
    const uint32_t *leaves
) {
    if (n_bits == 0u) {
        qua_x(kb, leaves[0]);
        return;
    }
    uint32_t step0 = 1u << (n_bits - 1u);
    qua_x(kb, addr_msb[0]);
    qua_cx(kb, addr_msb[0], leaves[0]);
    qua_x(kb, addr_msb[0]);
    qua_cx(kb, addr_msb[0], leaves[step0]);
    for (uint32_t lvl = 1u; lvl < n_bits; ++lvl) {
        uint32_t step = 1u << (n_bits - 1u - lvl);
        uint32_t cnt  = 1u << lvl;
        uint32_t bit  = addr_msb[lvl];
        for (uint32_t k = 0u; k < cnt; ++k) {
            uint32_t d0 = leaves[2u * k * step];
            uint32_t d1 = leaves[(2u * k + 1u) * step];
            qua_and_compute(kb, d0, bit, d1);
            qua_cx(kb, d1, d0);
        }
    }
}

static inline void qua_bitop_onehot_uncompute_mbuc(
    qua_builder_t *kb,
    const uint32_t *addr_msb,
    uint32_t n_bits,
    const uint32_t *leaves,
    uint32_t *next_cbit
) {
    if (n_bits == 0u) {
        qua_x(kb, leaves[0]);
        return;
    }
    for (uint32_t l_rev = n_bits - 1u; l_rev >= 1u; --l_rev) {
        uint32_t step = 1u << (n_bits - 1u - l_rev);
        uint32_t cnt  = 1u << l_rev;
        uint32_t bit  = addr_msb[l_rev];
        for (uint32_t k = 0u; k < cnt; ++k) {
            uint32_t d0 = leaves[2u * k * step];
            uint32_t d1 = leaves[(2u * k + 1u) * step];
            qua_cx(kb, d1, d0);
            qua_and_uncompute(kb, d0, bit, d1, (*next_cbit)++);
        }
    }
    uint32_t step0 = 1u << (n_bits - 1u);
    qua_cx(kb, addr_msb[0], leaves[step0]);
    qua_x(kb, addr_msb[0]);
    qua_cx(kb, addr_msb[0], leaves[0]);
    qua_x(kb, addr_msb[0]);
}

/* ============================================================================
 * Virtual-Slot-0 One-Hot Gather + Berry-Gidney 1-Bit Phase-Fixup QROAM
 *   1. Constant-0 column pruning (`w_act <= w_bits`).
 *   2. Virtual Slot 0: loads `T[q, 0]` directly into `target_qs` (`0` slot-0 ancillas!)
 *      and `D[q, d] = T[q, d] ^ T[q, 0]` into `(B_ss - 1) * w_act` clean slot ancillas.
 *   3. Optional Top-Spine Collapse (`spine_collapse = 3` saves 1 spine qubit;
 *      `spine_collapse = 4` saves 2 spine qubits via immediate MBUC pair release).
 *   4. Non-Zero One-Hot Gather (`h_1..h_{B_ss-1}` in `B_ss - b_sh - 1` CCX + 0-CCX MBUC,
 *      `1x` CCX gather instead of `2x` CSWAP forward+reverse).
 *   5. 0-CCX HMR slot reset + 1-Bit Split-Unary Phase Fixup (`BIT_INVERT_IF` + `CZ_IF`).
 * ============================================================================ */

static inline void qua_bitop_qrom_sub_rec(
    const qua_bitop_qrom_ctx_t *ctx,
    uint32_t ctrl,
    uint32_t sub_lvl,
    uint32_t sub_bits,
    const uint32_t *sub_addr_msb,
    const uint32_t *sub_anc,
    uint32_t lo_blk,
    uint32_t hi_blk
) {
    uint32_t n_blocks = (ctx->num_rows + ctx->b_ss - 1u) / ctx->b_ss;
    if (lo_blk >= n_blocks) return;
    if (hi_blk <= n_blocks && (hi_blk - lo_blk == 1u || qua_bitop_qrom_block_const(ctx, lo_blk, hi_blk))) {
        qua_bitop_qrom_leaf(ctx, lo_blk, ctrl);
        return;
    }
    if (sub_lvl >= sub_bits) {
        qua_bitop_qrom_leaf(ctx, lo_blk, ctrl);
        return;
    }
    uint32_t mid = (lo_blk + hi_blk) / 2u;
    if (n_blocks <= mid) {
        qua_bitop_qrom_sub_rec(ctx, ctrl, sub_lvl + 1u, sub_bits, sub_addr_msb, sub_anc, lo_blk, mid);
        return;
    }
    uint32_t bit_q = sub_addr_msb[sub_lvl];
    uint32_t child = sub_anc[sub_lvl];
    qua_x(ctx->kb, bit_q);
    qua_and_compute(ctx->kb, ctrl, bit_q, child);
    qua_x(ctx->kb, bit_q);
    qua_bitop_qrom_sub_rec(ctx, child, sub_lvl + 1u, sub_bits, sub_addr_msb, sub_anc, lo_blk, mid);
    qua_cx(ctx->kb, ctrl, child);
    qua_bitop_qrom_sub_rec(ctx, child, sub_lvl + 1u, sub_bits, sub_addr_msb, sub_anc, mid, hi_blk);
    qua_and_uncompute(ctx->kb, ctrl, bit_q, child, (*ctx->next_cbit)++);
}

static inline void qua_bitop_qrom_walk_collapsed(
    const qua_bitop_qrom_ctx_t *ctx,
    uint32_t spine_collapse
) {
    uint32_t hi_bits = ctx->hi_bits;
    uint32_t n_blocks = (ctx->num_rows + ctx->b_ss - 1u) / ctx->b_ss;
    if (spine_collapse == 5u && hi_bits >= 8u) {
        uint32_t sub_bits = hi_bits - 5u;
        uint32_t sub_span = 1u << sub_bits;
        uint32_t root_ctrl = ctx->ancillas[0];
        const uint32_t *sub_anc = ctx->ancillas + 1u; /* length >= hi_bits - 5 >= 3 */
        uint32_t p01 = sub_anc[0], p23 = sub_anc[1], p0123 = sub_anc[2];
        for (uint32_t pref = 0u; pref < 32u; ++pref) {
            uint32_t lo_blk = pref * sub_span;
            if (lo_blk >= n_blocks) break;
            uint32_t hi_blk = lo_blk + sub_span;
            for (uint32_t b = 0u; b < 5u; ++b) {
                if (((pref >> (4u - b)) & 1u) == 0u) qua_x(ctx->kb, ctx->addr_msb[b]);
            }
            qua_and_compute(ctx->kb, ctx->addr_msb[0], ctx->addr_msb[1], p01);
            qua_and_compute(ctx->kb, ctx->addr_msb[2], ctx->addr_msb[3], p23);
            qua_and_compute(ctx->kb, p01, p23, p0123);
            qua_and_compute(ctx->kb, p0123, ctx->addr_msb[4], root_ctrl);
            qua_and_uncompute(ctx->kb, p01, p23, p0123, (*ctx->next_cbit)++);
            qua_and_uncompute(ctx->kb, ctx->addr_msb[2], ctx->addr_msb[3], p23, (*ctx->next_cbit)++);
            qua_and_uncompute(ctx->kb, ctx->addr_msb[0], ctx->addr_msb[1], p01, (*ctx->next_cbit)++);

            qua_bitop_qrom_sub_rec(ctx, root_ctrl, 0u, sub_bits, ctx->addr_msb + 5u, sub_anc, lo_blk, hi_blk);

            qua_and_compute(ctx->kb, ctx->addr_msb[0], ctx->addr_msb[1], p01);
            qua_and_compute(ctx->kb, ctx->addr_msb[2], ctx->addr_msb[3], p23);
            qua_and_compute(ctx->kb, p01, p23, p0123);
            qua_and_uncompute(ctx->kb, p0123, ctx->addr_msb[4], root_ctrl, (*ctx->next_cbit)++);
            qua_and_uncompute(ctx->kb, p01, p23, p0123, (*ctx->next_cbit)++);
            qua_and_uncompute(ctx->kb, ctx->addr_msb[2], ctx->addr_msb[3], p23, (*ctx->next_cbit)++);
            qua_and_uncompute(ctx->kb, ctx->addr_msb[0], ctx->addr_msb[1], p01, (*ctx->next_cbit)++);
            for (uint32_t b = 0u; b < 5u; ++b) {
                if (((pref >> (4u - b)) & 1u) == 0u) qua_x(ctx->kb, ctx->addr_msb[b]);
            }
        }
        return;
    }
    if (spine_collapse == 4u && hi_bits >= 6u) {
        uint32_t sub_bits = hi_bits - 4u;
        uint32_t sub_span = 1u << sub_bits;
        uint32_t root_ctrl = ctx->ancillas[0];
        const uint32_t *sub_anc = ctx->ancillas + 1u; /* length >= hi_bits - 4 >= 2 */
        uint32_t p01 = sub_anc[0], p23 = sub_anc[1];
        for (uint32_t pref = 0u; pref < 16u; ++pref) {
            uint32_t lo_blk = pref * sub_span;
            if (lo_blk >= n_blocks) break;
            uint32_t hi_blk = lo_blk + sub_span;
            for (uint32_t b = 0u; b < 4u; ++b) {
                if (((pref >> (3u - b)) & 1u) == 0u) qua_x(ctx->kb, ctx->addr_msb[b]);
            }
            qua_and_compute(ctx->kb, ctx->addr_msb[0], ctx->addr_msb[1], p01);
            qua_and_compute(ctx->kb, ctx->addr_msb[2], ctx->addr_msb[3], p23);
            qua_and_compute(ctx->kb, p01, p23, root_ctrl);
            qua_and_uncompute(ctx->kb, ctx->addr_msb[2], ctx->addr_msb[3], p23, (*ctx->next_cbit)++);
            qua_and_uncompute(ctx->kb, ctx->addr_msb[0], ctx->addr_msb[1], p01, (*ctx->next_cbit)++);

            qua_bitop_qrom_sub_rec(ctx, root_ctrl, 0u, sub_bits, ctx->addr_msb + 4u, sub_anc, lo_blk, hi_blk);

            qua_and_compute(ctx->kb, ctx->addr_msb[0], ctx->addr_msb[1], p01);
            qua_and_compute(ctx->kb, ctx->addr_msb[2], ctx->addr_msb[3], p23);
            qua_and_uncompute(ctx->kb, p01, p23, root_ctrl, (*ctx->next_cbit)++);
            qua_and_uncompute(ctx->kb, ctx->addr_msb[2], ctx->addr_msb[3], p23, (*ctx->next_cbit)++);
            qua_and_uncompute(ctx->kb, ctx->addr_msb[0], ctx->addr_msb[1], p01, (*ctx->next_cbit)++);
            for (uint32_t b = 0u; b < 4u; ++b) {
                if (((pref >> (3u - b)) & 1u) == 0u) qua_x(ctx->kb, ctx->addr_msb[b]);
            }
        }
        return;
    }
    if (spine_collapse == 3u && hi_bits >= 4u) {
        uint32_t sub_bits = hi_bits - 3u;
        uint32_t sub_span = 1u << sub_bits;
        uint32_t root_ctrl = ctx->ancillas[0];
        const uint32_t *sub_anc = ctx->ancillas + 1u; /* length >= hi_bits - 3 >= 1 */
        uint32_t p01 = sub_anc[0];
        for (uint32_t pref = 0u; pref < 8u; ++pref) {
            uint32_t lo_blk = pref * sub_span;
            if (lo_blk >= n_blocks) break;
            uint32_t hi_blk = lo_blk + sub_span;
            for (uint32_t b = 0u; b < 3u; ++b) {
                if (((pref >> (2u - b)) & 1u) == 0u) qua_x(ctx->kb, ctx->addr_msb[b]);
            }
            qua_and_compute(ctx->kb, ctx->addr_msb[0], ctx->addr_msb[1], p01);
            qua_and_compute(ctx->kb, p01, ctx->addr_msb[2], root_ctrl);
            qua_and_uncompute(ctx->kb, ctx->addr_msb[0], ctx->addr_msb[1], p01, (*ctx->next_cbit)++);

            qua_bitop_qrom_sub_rec(ctx, root_ctrl, 0u, sub_bits, ctx->addr_msb + 3u, sub_anc, lo_blk, hi_blk);

            qua_and_compute(ctx->kb, ctx->addr_msb[0], ctx->addr_msb[1], p01);
            qua_and_uncompute(ctx->kb, p01, ctx->addr_msb[2], root_ctrl, (*ctx->next_cbit)++);
            qua_and_uncompute(ctx->kb, ctx->addr_msb[0], ctx->addr_msb[1], p01, (*ctx->next_cbit)++);
            for (uint32_t b = 0u; b < 3u; ++b) {
                if (((pref >> (2u - b)) & 1u) == 0u) qua_x(ctx->kb, ctx->addr_msb[b]);
            }
        }
        return;
    }
    qua_bitop_qrom_rec(ctx, UINT32_MAX, 0u, 0u, 1u << hi_bits);
}

static inline uint32_t qua_bitop_spine_qubits(uint32_t hi_bits, uint32_t spine_collapse) {
    if (hi_bits <= 1u) return 0u;
    if (spine_collapse == 5u && hi_bits >= 8u) return hi_bits - 4u;
    if (spine_collapse == 4u && hi_bits >= 6u) return hi_bits - 3u;
    if (spine_collapse == 3u && hi_bits >= 4u) return hi_bits - 2u;
    return hi_bits - 1u;
}

typedef struct {
    uint32_t w_act;
    uint32_t anc_needed;
    uint64_t total_ccx;
} qua_bitop_qroam_res_t;

static inline void qua_bitop_qroam_vslot0_resources(
    const uint64_t *table,
    uint32_t num_rows,
    uint32_t stride_u64,
    uint32_t addr_bits,
    uint32_t w_bits,
    uint32_t b_sh,
    uint32_t spine_collapse,
    uint32_t free_anc_cap,
    qua_bitop_qroam_res_t *out_res
) {
    if (out_res) memset(out_res, 0, sizeof(*out_res));
    if (!table || w_bits == 0u || num_rows == 0u) return;

    qua_bitop_qrom_ctx_t probe = {
        NULL, table, num_rows, stride_u64, w_bits, NULL, 0u, NULL, NULL, 1u, NULL, false, NULL
    };
    uint32_t w_act = 0u;
    uint32_t act_cols_stack[256];
    uint32_t *act_cols = (w_bits <= 256u) ? act_cols_stack : (uint32_t *)malloc((size_t)w_bits * sizeof(uint32_t));
    if (!act_cols) abort();
    for (uint32_t j = 0u; j < w_bits; ++j) {
        for (uint32_t r = 0u; r < num_rows; ++r) {
            if (qua_bitop_qrom_bit(&probe, r, j)) {
                act_cols[w_act++] = j;
                break;
            }
        }
    }
    if (w_act == 0u) {
        if (act_cols != act_cols_stack) free(act_cols);
        return;
    }

    if (b_sh > addr_bits) b_sh = addr_bits;
    uint32_t b_ss = 1u << b_sh;
    uint32_t hi_bits = addr_bits - b_sh;
    uint32_t n_blocks = (num_rows + b_ss - 1u) / b_ss;
    uint32_t spine_q = qua_bitop_spine_qubits(hi_bits, spine_collapse);

    uint64_t walk_ccx = 0u;
    if (spine_collapse == 5u && hi_bits >= 8u) {
        uint32_t sub_span = 1u << (hi_bits - 5u);
        uint32_t n_pref = (n_blocks + sub_span - 1u) / sub_span;
        walk_ccx = (uint64_t)n_pref * 7ull + ((n_blocks > n_pref) ? (uint64_t)(n_blocks - n_pref) : 0ull);
    } else if (spine_collapse == 4u && hi_bits >= 6u) {
        uint32_t sub_span = 1u << (hi_bits - 4u);
        uint32_t n_pref = (n_blocks + sub_span - 1u) / sub_span;
        walk_ccx = (uint64_t)n_pref * 5ull + ((n_blocks > n_pref) ? (uint64_t)(n_blocks - n_pref) : 0ull);
    } else if (spine_collapse == 3u && hi_bits >= 4u) {
        uint32_t sub_span = 1u << (hi_bits - 3u);
        uint32_t n_pref = (n_blocks + sub_span - 1u) / sub_span;
        walk_ccx = (uint64_t)n_pref * 3ull + ((n_blocks > n_pref) ? (uint64_t)(n_blocks - n_pref) : 0ull);
    } else {
        walk_ccx = (n_blocks >= 2u) ? (uint64_t)(n_blocks - 2u) : 0ull;
    }

    if (b_sh == 0u) {
        if (out_res) {
            out_res->w_act      = w_act;
            out_res->anc_needed = spine_q;
            out_res->total_ccx  = walk_ccx;
        }
        if (act_cols != act_cols_stack) free(act_cols);
        return;
    }

    uint64_t dec_ccx = (b_sh >= 2u) ? (uint64_t)(b_ss - b_sh - 1u) : 0ull;
    uint64_t gath_ccx = 0u;
    for (uint32_t d = 1u; d < b_ss; ++d) {
        for (uint32_t w = 0u; w < w_act; ++w) {
            uint32_t col = act_cols[w];
            for (uint32_t q = 0u; q < n_blocks; ++q) {
                uint32_t r0 = q * b_ss, rd = r0 + d;
                bool b0 = qua_bitop_qrom_bit(&probe, r0, col);
                bool bd = (rd < num_rows) ? qua_bitop_qrom_bit(&probe, rd, col) : false;
                if (b0 ^ bd) { ++gath_ccx; break; }
            }
        }
    }
    if (act_cols != act_cols_stack) free(act_cols);

    /* Zero-differential slots need 0 physical qubits! */
    uint32_t slot_q = (uint32_t)gath_ccx;
    uint32_t gather_q = (b_sh >= 2u) ? (b_ss - 1u) : 0u;
    uint32_t peak_q = slot_q + ((spine_q > gather_q) ? spine_q : gather_q);
    uint32_t avail_phase = (peak_q > free_anc_cap) ? peak_q : free_anc_cap;

    uint64_t phase_ccx = 0u;
    if (hi_bits >= 2u) {
        uint32_t nh = hi_bits / 2u, nl = hi_bits - nh;
        uint32_t need0 = (1u << nh) + (1u << nl);
        if (avail_phase >= need0) {
            phase_ccx = (uint64_t)(((1u << nh) >= 2u ? (1u << nh) - 2u : 0u) +
                                   ((1u << nl) >= 2u ? (1u << nl) - 2u : 0u));
        } else if (hi_bits >= 3u) {
            uint32_t rem1 = hi_bits - 1u;
            uint32_t nm1 = rem1 / 2u, nl1 = rem1 - nm1;
            uint32_t need1 = 1u + (1u << nm1) + (1u << nl1);
            if (avail_phase >= need1) {
                uint32_t act_m = (n_blocks > (1u << rem1)) ? (1u << nm1) : ((n_blocks + (1u << nl1) - 1u) >> nl1);
                phase_ccx = (uint64_t)(((1u << nm1) >= 2u ? (1u << nm1) - 2u : 0u) +
                                       ((1u << nl1) >= 2u ? (1u << nl1) - 2u : 0u) + act_m);
            } else if (hi_bits >= 4u) {
                uint32_t rem2 = hi_bits - 2u;
                uint32_t nm2 = rem2 / 2u, nl2 = rem2 - nm2;
                uint32_t need2 = 2u + (1u << nm2) + (1u << nl2);
                if (peak_q < need2 && free_anc_cap < need2) peak_q = need2;
                uint32_t sub2 = 1u << rem2;
                uint32_t n_top2 = (n_blocks + sub2 - 1u) / sub2;
                phase_ccx = (uint64_t)(((1u << nm2) >= 2u ? (1u << nm2) - 2u : 0u) +
                                       ((1u << nl2) >= 2u ? (1u << nl2) - 2u : 0u) +
                                       n_top2 * (1u + (1u << nm2)));
            } else {
                if (peak_q < need1) peak_q = need1;
                phase_ccx = 4u;
            }
        } else {
            if (peak_q < need0) peak_q = need0;
        }
    }

    if (out_res) {
        out_res->w_act      = w_act;
        out_res->anc_needed = peak_q;
        out_res->total_ccx  = walk_ccx + dec_ccx + gath_ccx + phase_ccx;
    }
}

static inline void qua_bitop_emit_qroam_vslot0(
    qua_builder_t *kb,
    const uint64_t *table,
    uint32_t num_rows,
    uint32_t stride_u64,
    const uint32_t *addr_lsb,
    uint32_t addr_bits,
    const uint32_t *target_qs,
    uint32_t w_bits,
    uint32_t b_sh,
    uint32_t spine_collapse,
    const uint32_t *anc_pool,
    uint32_t anc_pool_len,
    uint32_t *next_cbit
) {
    if (num_rows == 0u || w_bits == 0u) return;
    if (b_sh > addr_bits) b_sh = addr_bits;

    /* 1. Prune identically-zero columns across [0, num_rows) */
    uint32_t *act_cols = (uint32_t *)malloc((size_t)w_bits * sizeof(uint32_t));
    uint32_t *act_tgts = (uint32_t *)malloc((size_t)w_bits * sizeof(uint32_t));
    if (!act_cols || !act_tgts) abort();
    uint32_t w_act = 0u;
    qua_bitop_qrom_ctx_t probe = {
        kb, table, num_rows, stride_u64, w_bits, NULL, 0u, NULL, NULL, 1u, next_cbit, false, NULL
    };
    for (uint32_t j = 0u; j < w_bits; ++j) {
        bool nz = false;
        for (uint32_t r = 0u; r < num_rows; ++r) {
            if (qua_bitop_qrom_bit(&probe, r, j)) { nz = true; break; }
        }
        if (nz) {
            act_cols[w_act] = j;
            act_tgts[w_act] = target_qs[j];
            ++w_act;
        }
    }
    if (w_act == 0u) {
        free(act_cols);
        free(act_tgts);
        return;
    }

    uint32_t b_ss = 1u << b_sh;
    uint32_t hi_bits = addr_bits - b_sh;
    uint32_t n_blocks = (num_rows + b_ss - 1u) / b_ss;
    uint32_t addr_msb[32];
    for (uint32_t i = 0u; i < hi_bits; ++i) {
        addr_msb[i] = addr_lsb[addr_bits - 1u - i];
    }

    uint32_t pack_stride = (w_act + 63u) / 64u;
    uint32_t total_rows_pack = n_blocks * b_ss;
    uint64_t *pack_tab = (uint64_t *)calloc((size_t)total_rows_pack * (size_t)pack_stride, sizeof(uint64_t));
    if (!pack_tab) abort();

    if (b_sh == 0u) {
        for (uint32_t r = 0u; r < num_rows; ++r) {
            for (uint32_t w = 0u; w < w_act; ++w) {
                if (qua_bitop_qrom_bit(&probe, r, act_cols[w])) {
                    pack_tab[(size_t)r * pack_stride + (w >> 6u)] |= (UINT64_C(1) << (w & 63u));
                }
            }
        }
        qua_bitop_qrom_ctx_t ctx0 = {
            kb, pack_tab, num_rows, pack_stride, w_act, addr_msb, hi_bits,
            anc_pool, act_tgts, 1u, next_cbit, false, NULL
        };
        qua_bitop_qrom_walk_collapsed(&ctx0, spine_collapse);
        free(pack_tab);
        free(act_cols);
        free(act_tgts);
        return;
    }

    /* 2. Build Virtual-Slot-0 table:
     *    row (q*b_ss + 0) holds T[q*b_ss + 0] (XORed directly into act_tgts!)
     *    row (q*b_ss + d) holds D[q, d] = T[q*b_ss + d] ^ T[q*b_ss + 0] (d = 1..b_ss-1) */
    uint32_t slot_total = (b_ss - 1u) * w_act;
    uint8_t *slot_ever_nz = (uint8_t *)calloc(slot_total, 1u);
    if (!slot_ever_nz) abort();
    for (uint32_t q = 0u; q < n_blocks; ++q) {
        uint32_t r0 = q * b_ss;
        for (uint32_t w = 0u; w < w_act; ++w) {
            bool b0 = qua_bitop_qrom_bit(&probe, r0, act_cols[w]);
            if (b0) {
                pack_tab[(size_t)r0 * pack_stride + (w >> 6u)] |= (UINT64_C(1) << (w & 63u));
            }
            for (uint32_t d = 1u; d < b_ss; ++d) {
                uint32_t rd = r0 + d;
                bool bd = (rd < num_rows) ? qua_bitop_qrom_bit(&probe, rd, act_cols[w]) : false;
                if (bd ^ b0) {
                    pack_tab[(size_t)rd * pack_stride + (w >> 6u)] |= (UINT64_C(1) << (w & 63u));
                    slot_ever_nz[(d - 1u) * w_act + w] = 1u;
                }
            }
        }
    }

    /* Allocate physical slot qubits ONLY for non-zero differential columns */
    uint32_t *slot_qs = (uint32_t *)malloc((size_t)(slot_total ? slot_total : 1u) * sizeof(uint32_t));
    if (!slot_qs) abort();
    uint32_t nz_slots = 0u;
    for (uint32_t idx = 0u; idx < slot_total; ++idx) {
        slot_qs[idx] = slot_ever_nz[idx] ? anc_pool[nz_slots++] : UINT32_MAX;
    }
    const uint32_t *work_anc = anc_pool + nz_slots;
    uint32_t *all_slots = (uint32_t *)malloc((size_t)(b_ss * w_act) * sizeof(uint32_t));
    if (!all_slots) abort();
    for (uint32_t w = 0u; w < w_act; ++w) all_slots[w] = act_tgts[w];
    for (uint32_t idx = 0u; idx < slot_total; ++idx) all_slots[w_act + idx] = slot_qs[idx];

    /* Step A: Forward Unary Walk over hi_bits (total_rows_pack covers trailing D = 0^b0 entries!) */
    qua_bitop_qrom_ctx_t fwd_ctx = {
        kb, pack_tab, total_rows_pack, pack_stride, w_act, addr_msb, hi_bits,
        work_anc, all_slots, b_ss, next_cbit, false, NULL
    };
    qua_bitop_qrom_walk_collapsed(&fwd_ctx, spine_collapse);

    /* Step B: Non-Zero One-Hot Gather over addr_lsb[0 .. b_sh - 1] */
    if (b_sh == 1u) {
        for (uint32_t w = 0u; w < w_act; ++w) {
            if (slot_ever_nz[w]) qua_ccx(kb, addr_lsb[0], slot_qs[w], act_tgts[w]);
        }
    } else {
        uint32_t r_msb[16];
        for (uint32_t i = 0u; i < b_sh; ++i) r_msb[i] = addr_lsb[b_sh - 1u - i];
        /* h_qs[1 .. b_ss - 1] mapped onto work_anc[0 .. b_ss - 2] without pointer underflow UB */
        uint32_t step0 = 1u << (b_sh - 2u);
        qua_and_compute(kb, r_msb[0], r_msb[1], work_anc[(3u * step0) - 1u]);
        qua_cx(kb, r_msb[0], work_anc[(2u * step0) - 1u]);
        qua_cx(kb, work_anc[(3u * step0) - 1u], work_anc[(2u * step0) - 1u]);
        qua_cx(kb, r_msb[1], work_anc[(1u * step0) - 1u]);
        qua_cx(kb, work_anc[(3u * step0) - 1u], work_anc[(1u * step0) - 1u]);

        for (uint32_t lvl = 2u; lvl < b_sh; ++lvl) {
            uint32_t step = 1u << (b_sh - 1u - lvl);
            uint32_t m_prev = 1u << lvl;
            uint32_t bit = r_msb[lvl];
            for (uint32_t k = 1u; k < m_prev; ++k) {
                uint32_t d0 = 2u * k * step, d1 = (2u * k + 1u) * step;
                qua_and_compute(kb, work_anc[d0 - 1u], bit, work_anc[d1 - 1u]);
                qua_cx(kb, work_anc[d1 - 1u], work_anc[d0 - 1u]);
            }
            qua_cx(kb, bit, work_anc[step - 1u]);
            for (uint32_t k = 1u; k < m_prev; ++k) {
                qua_cx(kb, work_anc[((2u * k + 1u) * step) - 1u], work_anc[step - 1u]);
            }
        }

        for (uint32_t d = 1u; d < b_ss; ++d) {
            for (uint32_t w = 0u; w < w_act; ++w) {
                uint32_t idx = (d - 1u) * w_act + w;
                if (slot_ever_nz[idx]) qua_ccx(kb, work_anc[d - 1u], slot_qs[idx], act_tgts[w]);
            }
        }

        /* 0-Toffoli MBUC uncompute of h_qs[1 .. b_ss - 1] */
        for (uint32_t l_rev = b_sh - 1u; l_rev >= 2u; --l_rev) {
            uint32_t step = 1u << (b_sh - 1u - l_rev);
            uint32_t m_prev = 1u << l_rev;
            uint32_t bit = r_msb[l_rev];
            for (uint32_t k = 1u; k < m_prev; ++k) {
                qua_cx(kb, work_anc[((2u * k + 1u) * step) - 1u], work_anc[step - 1u]);
            }
            qua_cx(kb, bit, work_anc[step - 1u]);
            for (uint32_t k = 1u; k < m_prev; ++k) {
                uint32_t d0 = 2u * k * step, d1 = (2u * k + 1u) * step;
                qua_cx(kb, work_anc[d1 - 1u], work_anc[d0 - 1u]);
                qua_and_uncompute(kb, work_anc[d0 - 1u], bit, work_anc[d1 - 1u], (*next_cbit)++);
            }
        }
        qua_cx(kb, work_anc[(3u * step0) - 1u], work_anc[(1u * step0) - 1u]);
        qua_cx(kb, r_msb[1], work_anc[(1u * step0) - 1u]);
        qua_cx(kb, work_anc[(3u * step0) - 1u], work_anc[(2u * step0) - 1u]);
        qua_cx(kb, r_msb[0], work_anc[(2u * step0) - 1u]);
        qua_and_uncompute(kb, r_msb[0], r_msb[1], work_anc[(3u * step0) - 1u], (*next_cbit)++);
    }

    /* Step C: HMR all active slot qubits + 1-Bit Split-Unary Phase Fixup */
    uint32_t *m_slots = (uint32_t *)malloc((size_t)slot_total * sizeof(uint32_t));
    if (!m_slots) abort();
    for (uint32_t idx = 0u; idx < slot_total; ++idx) {
        if (slot_ever_nz[idx]) {
            m_slots[idx] = (*next_cbit)++;
            qua_hmr(kb, slot_qs[idx], m_slots[idx]);
        } else {
            m_slots[idx] = UINT32_MAX;
        }
    }

    uint32_t pbit = (*next_cbit)++;
    if (hi_bits == 0u) {
        bool any = false;
        qua_bit_store0(kb, pbit);
        for (uint32_t d = 1u; d < b_ss; ++d) {
            for (uint32_t w = 0u; w < w_act; ++w) {
                uint64_t wd = pack_tab[(size_t)d * pack_stride + (w >> 6u)];
                if ((wd >> (w & 63u)) & 1ull) {
                    qua_bit_invert_if(kb, pbit, m_slots[(d - 1u) * w_act + w]);
                    any = true;
                }
            }
        }
        if (any) {
            qcc_kmb_neg_if(kb, pbit);
        }
    } else if (hi_bits == 1u) {
        for (uint32_t q = 0u; q < n_blocks; ++q) {
            bool any = false;
            qua_bit_store0(kb, pbit);
            for (uint32_t d = 1u; d < b_ss; ++d) {
                uint32_t rd = q * b_ss + d;
                for (uint32_t w = 0u; w < w_act; ++w) {
                    uint64_t wd = pack_tab[(size_t)rd * pack_stride + (w >> 6u)];
                    if ((wd >> (w & 63u)) & 1ull) {
                        qua_bit_invert_if(kb, pbit, m_slots[(d - 1u) * w_act + w]);
                        any = true;
                    }
                }
            }
            if (any) {
                if (q == 0u) qua_x(kb, addr_msb[0]);
                qua_z_if(kb, addr_msb[0], pbit);
                if (q == 0u) qua_x(kb, addr_msb[0]);
            }
        }
    } else {
        /* All anc_pool[0 .. anc_pool_len - 1] are now clean |0>! */
        uint32_t nh0 = hi_bits / 2u, nl0 = hi_bits - nh0;
        uint32_t need0 = (1u << nh0) + (1u << nl0);
        if (anc_pool_len >= need0) {
            const uint32_t *hi_lv = anc_pool;
            const uint32_t *lo_lv = anc_pool + (1u << nh0);
            qua_bitop_onehot_compute(kb, addr_msb, nh0, hi_lv);
            qua_bitop_onehot_compute(kb, addr_msb + nh0, nl0, lo_lv);
            uint32_t lo_mask = (1u << nl0) - 1u;
            for (uint32_t q = 0u; q < n_blocks; ++q) {
                bool any = false;
                qua_bit_store0(kb, pbit);
                for (uint32_t d = 1u; d < b_ss; ++d) {
                    uint32_t rd = q * b_ss + d;
                    for (uint32_t w = 0u; w < w_act; ++w) {
                        uint64_t wd = pack_tab[(size_t)rd * pack_stride + (w >> 6u)];
                        if ((wd >> (w & 63u)) & 1ull) {
                            qua_bit_invert_if(kb, pbit, m_slots[(d - 1u) * w_act + w]);
                            any = true;
                        }
                    }
                }
                if (any) {
                    qua_cz_if(kb, hi_lv[q >> nl0], lo_lv[q & lo_mask], pbit);
                }
            }
            qua_bitop_onehot_uncompute_mbuc(kb, addr_msb + nh0, nl0, lo_lv, next_cbit);
            qua_bitop_onehot_uncompute_mbuc(kb, addr_msb, nh0, hi_lv, next_cbit);
        } else {
            uint32_t rem1 = hi_bits - 1u;
            uint32_t nm1 = rem1 / 2u, nl1 = rem1 - nm1;
            uint32_t need1 = 1u + (1u << nm1) + (1u << nl1);
            if (anc_pool_len >= need1) {
                uint32_t anc_m = anc_pool[0];
                const uint32_t *mid_lv = anc_pool + 1u;
                const uint32_t *lo_lv  = mid_lv + (1u << nm1);
                qua_bitop_onehot_compute(kb, addr_msb + 1u, nm1, mid_lv);
                qua_bitop_onehot_compute(kb, addr_msb + 1u + nm1, nl1, lo_lv);
                uint32_t m_cnt = 1u << nm1, l_cnt = 1u << nl1;
                for (uint32_t m = 0u; m < m_cnt; ++m) {
                    uint32_t q0_base = (0u << rem1) | (m << nl1);
                    uint32_t q1_base = (1u << rem1) | (m << nl1);
                    if (q0_base >= n_blocks) break;
                    qua_x(kb, addr_msb[0]);
                    qua_and_compute(kb, addr_msb[0], mid_lv[m], anc_m);
                    qua_x(kb, addr_msb[0]);
                    for (uint32_t l = 0u; l < l_cnt && (q0_base + l) < n_blocks; ++l) {
                        uint32_t q = q0_base + l;
                        bool any = false;
                        qua_bit_store0(kb, pbit);
                        for (uint32_t d = 1u; d < b_ss; ++d) {
                            uint32_t rd = q * b_ss + d;
                            for (uint32_t w = 0u; w < w_act; ++w) {
                                uint64_t wd = pack_tab[(size_t)rd * pack_stride + (w >> 6u)];
                                if ((wd >> (w & 63u)) & 1ull) {
                                    qua_bit_invert_if(kb, pbit, m_slots[(d - 1u) * w_act + w]);
                                    any = true;
                                }
                            }
                        }
                        if (any) qua_cz_if(kb, anc_m, lo_lv[l], pbit);
                    }
                    if (q1_base < n_blocks) {
                        qua_cx(kb, mid_lv[m], anc_m);
                        for (uint32_t l = 0u; l < l_cnt && (q1_base + l) < n_blocks; ++l) {
                            uint32_t q = q1_base + l;
                            bool any = false;
                            qua_bit_store0(kb, pbit);
                            for (uint32_t d = 1u; d < b_ss; ++d) {
                                uint32_t rd = q * b_ss + d;
                                for (uint32_t w = 0u; w < w_act; ++w) {
                                    uint64_t wd = pack_tab[(size_t)rd * pack_stride + (w >> 6u)];
                                    if ((wd >> (w & 63u)) & 1ull) {
                                        qua_bit_invert_if(kb, pbit, m_slots[(d - 1u) * w_act + w]);
                                        any = true;
                                    }
                                }
                            }
                            if (any) qua_cz_if(kb, anc_m, lo_lv[l], pbit);
                        }
                        qua_and_uncompute(kb, addr_msb[0], mid_lv[m], anc_m, (*next_cbit)++);
                    } else {
                        qua_x(kb, addr_msb[0]);
                        qua_and_uncompute(kb, addr_msb[0], mid_lv[m], anc_m, (*next_cbit)++);
                        qua_x(kb, addr_msb[0]);
                    }
                }
                qua_bitop_onehot_uncompute_mbuc(kb, addr_msb + 1u + nm1, nl1, lo_lv, next_cbit);
                qua_bitop_onehot_uncompute_mbuc(kb, addr_msb + 1u, nm1, mid_lv, next_cbit);
            } else {
                uint32_t rem2 = hi_bits - 2u;
                uint32_t nm2 = rem2 / 2u, nl2 = rem2 - nm2;
                uint32_t p01 = anc_pool[0], anc_m = anc_pool[1];
                const uint32_t *mid_lv = anc_pool + 2u;
                const uint32_t *lo_lv  = mid_lv + (1u << nm2);
                qua_bitop_onehot_compute(kb, addr_msb + 2u, nm2, mid_lv);
                qua_bitop_onehot_compute(kb, addr_msb + 2u + nm2, nl2, lo_lv);
                uint32_t m_cnt = 1u << nm2, l_cnt = 1u << nl2;
                for (uint32_t tv = 0u; tv < 4u; ++tv) {
                    uint32_t q_top = tv << rem2;
                    if (q_top >= n_blocks) break;
                    if ((tv & 2u) == 0u) qua_x(kb, addr_msb[0]);
                    if ((tv & 1u) == 0u) qua_x(kb, addr_msb[1]);
                    qua_and_compute(kb, addr_msb[0], addr_msb[1], p01);
                    for (uint32_t m = 0u; m < m_cnt; ++m) {
                        uint32_t q_base = q_top | (m << nl2);
                        if (q_base >= n_blocks) break;
                        qua_and_compute(kb, p01, mid_lv[m], anc_m);
                        for (uint32_t l = 0u; l < l_cnt && (q_base + l) < n_blocks; ++l) {
                            uint32_t q = q_base + l;
                            bool any = false;
                            qua_bit_store0(kb, pbit);
                            for (uint32_t d = 1u; d < b_ss; ++d) {
                                uint32_t rd = q * b_ss + d;
                                for (uint32_t w = 0u; w < w_act; ++w) {
                                    uint64_t wd = pack_tab[(size_t)rd * pack_stride + (w >> 6u)];
                                    if ((wd >> (w & 63u)) & 1ull) {
                                        qua_bit_invert_if(kb, pbit, m_slots[(d - 1u) * w_act + w]);
                                        any = true;
                                    }
                                }
                            }
                            if (any) qua_cz_if(kb, anc_m, lo_lv[l], pbit);
                        }
                        qua_and_uncompute(kb, p01, mid_lv[m], anc_m, (*next_cbit)++);
                    }
                    qua_and_uncompute(kb, addr_msb[0], addr_msb[1], p01, (*next_cbit)++);
                    if ((tv & 1u) == 0u) qua_x(kb, addr_msb[1]);
                    if ((tv & 2u) == 0u) qua_x(kb, addr_msb[0]);
                }
                qua_bitop_onehot_uncompute_mbuc(kb, addr_msb + 2u + nm2, nl2, lo_lv, next_cbit);
                qua_bitop_onehot_uncompute_mbuc(kb, addr_msb + 2u, nm2, mid_lv, next_cbit);
            }
        }
    }

    free(m_slots);
    free(all_slots);
    free(slot_qs);
    free(slot_ever_nz);
    free(pack_tab);
    free(act_cols);
    free(act_tgts);
}

/* ============================================================================
 * 5. Composable Quantum-Algebraic Shortcut Combinators
 *    (HMR_Dual, LinF2_Fuse, PhaseSplit + Cuccaro, HGCD_Band, BasisIso_TSP)
 * ============================================================================ */

/* --- Combinator 1: HMR_Dual (Adaptive-Peel Ket[w] -> ParityDual -> |0>^w) ---
 * Measures `tgt[0..n_tgt-1]` in X-basis (0 CCX), immediately reuses those
 * freshly zeroed `tgt` qubits as one-hot leaves, and selects `p_peel in [0..3]`
 * so total extra pool ancillas fit within `free_anc_cap` whenever possible. */
typedef struct {
    uint32_t p_peel;
    uint32_t leaves_needed;
    uint32_t pool_anc_needed;
    uint64_t total_ccx;
} qua_bitop_hmr_lut_res_t;

static inline void qua_bitop_hmr_parity_lut_resources(
    uint32_t n_entries,
    uint32_t n_addr,
    uint32_t n_tgt,
    uint32_t free_anc_cap,
    qua_bitop_hmr_lut_res_t *out
) {
    memset(out, 0, sizeof(*out));
    if (n_tgt == 0u || n_addr <= 1u || n_entries <= 1u) return;

    uint32_t chosen_p = UINT32_MAX;
    for (uint32_t p = 0u; p <= 3u; ++p) {
        if (n_addr < p + 2u) break;
        if (p == 0u && n_addr > 8u) continue;
        uint32_t rem = n_addr - p;
        if (rem > 8u) continue;
        uint32_t nm = rem >> 1u;
        uint32_t nl = rem - nm;
        uint32_t need = p + (1u << nm) + (1u << nl);
        uint32_t extra = (need > n_tgt) ? (need - n_tgt) : 0u;
        chosen_p = p;
        if (extra <= free_anc_cap) break;
    }
    if (chosen_p == UINT32_MAX) {
        chosen_p = (n_addr > 2u) ? (n_addr - 2u) : 0u;
        if (chosen_p > 3u) chosen_p = 3u;
    }

    uint32_t rem = n_addr - chosen_p;
    uint32_t nm = rem >> 1u;
    uint32_t nl = rem - nm;
    uint32_t need = chosen_p + (1u << nm) + (1u << nl);
    uint64_t ccx = (uint64_t)((nm >= 2u) ? ((1u << nm) - 2u) : 0u) +
                   (uint64_t)((nl >= 2u) ? ((1u << nl) - 2u) : 0u);
    if (chosen_p == 1u) {
        uint32_t m_cnt = 1u << nm;
        for (uint32_t m = 0u; m < m_cnt; ++m) {
            if ((m << nl) >= n_entries) break;
            ccx += 1ull;
        }
    } else if (chosen_p >= 2u) {
        uint32_t top_cnt = 1u << chosen_p;
        uint32_t m_cnt   = 1u << nm;
        for (uint32_t tv = 0u; tv < top_cnt; ++tv) {
            uint32_t q_top = tv << rem;
            if (q_top >= n_entries) break;
            ccx += (uint64_t)(chosen_p - 1u);
            for (uint32_t m = 0u; m < m_cnt; ++m) {
                if ((q_top | (m << nl)) >= n_entries) break;
                ccx += 1ull;
            }
        }
    }
    out->p_peel          = chosen_p;
    out->leaves_needed   = need;
    out->pool_anc_needed = (need > n_tgt) ? (need - n_tgt) : 0u;
    out->total_ccx       = ccx;
}

static inline void qua_bitop_emit_hmr_parity_lut(
    qua_builder_t *kb,
    const uint64_t *table,
    uint32_t n_entries,
    const uint32_t *addr_lsb,
    uint32_t n_addr,
    const uint32_t *tgt,
    uint32_t n_tgt,
    const uint32_t *extra_anc,
    uint32_t extra_anc_len,
    uint32_t *next_cbit
) {
    if (n_tgt == 0u) return;
    if (n_tgt > 128u || n_addr > 16u) abort();
    uint32_t hmr_bits[128];
    for (uint32_t j = 0u; j < n_tgt; ++j) {
        hmr_bits[j] = (*next_cbit)++;
        qua_hmr(kb, tgt[j], hmr_bits[j]);
    }
    if (n_addr == 0u || n_entries == 0u) return;

    uint32_t addr_msb[16];
    for (uint32_t i = 0u; i < n_addr; ++i) addr_msb[i] = addr_lsb[n_addr - 1u - i];
    uint32_t pbit = (*next_cbit)++;

    if (n_addr == 1u || n_entries == 1u) {
        uint32_t lim = (n_addr == 1u && n_entries >= 2u) ? 2u : 1u;
        for (uint32_t a = 0u; a < lim; ++a) {
            bool any = false;
            qua_bit_store0(kb, pbit);
            uint64_t val = table[a];
            for (uint32_t j = 0u; j < n_tgt; ++j) {
                if ((val >> j) & 1ull) {
                    qua_bit_invert_if(kb, pbit, hmr_bits[j]);
                    any = true;
                }
            }
            if (any) {
                if (n_addr == 0u) continue;
                if (a == 0u) qua_x(kb, addr_msb[0]);
                qua_z_if(kb, addr_msb[0], pbit);
                if (a == 0u) qua_x(kb, addr_msb[0]);
            }
        }
        return;
    }

    qua_bitop_hmr_lut_res_t res;
    qua_bitop_hmr_parity_lut_resources(n_entries, n_addr, n_tgt, extra_anc_len, &res);
    if (res.leaves_needed > 128u) abort();

    uint32_t leaves_buf[128];
    uint32_t borrow = (n_tgt < res.leaves_needed) ? n_tgt : res.leaves_needed;
    for (uint32_t i = 0u; i < borrow; ++i) leaves_buf[i] = tgt[i];
    for (uint32_t i = 0u; i < res.pool_anc_needed; ++i) leaves_buf[borrow + i] = extra_anc[i];

    uint32_t p = res.p_peel;
    uint32_t rem = n_addr - p;
    uint32_t nm = rem >> 1u;
    uint32_t nl = rem - nm;
    uint32_t m_cnt = 1u << nm;
    uint32_t l_cnt = 1u << nl;

    if (p == 0u) {
        const uint32_t *hi_lv = leaves_buf;
        const uint32_t *lo_lv = leaves_buf + m_cnt;
        qua_bitop_onehot_compute(kb, addr_msb, nm, hi_lv);
        qua_bitop_onehot_compute(kb, addr_msb + nm, nl, lo_lv);
        for (uint32_t a = 0u; a < n_entries; ++a) {
            bool any = false;
            qua_bit_store0(kb, pbit);
            uint64_t val = table[a];
            for (uint32_t j = 0u; j < n_tgt; ++j) {
                if ((val >> j) & 1ull) {
                    qua_bit_invert_if(kb, pbit, hmr_bits[j]);
                    any = true;
                }
            }
            if (any) qua_cz_if(kb, hi_lv[a >> nl], lo_lv[a & (l_cnt - 1u)], pbit);
        }
        qua_bitop_onehot_uncompute_mbuc(kb, addr_msb + nm, nl, lo_lv, next_cbit);
        qua_bitop_onehot_uncompute_mbuc(kb, addr_msb, nm, hi_lv, next_cbit);
        return;
    }

    if (p == 1u) {
        uint32_t anc_m = leaves_buf[0];
        const uint32_t *mid_lv = leaves_buf + 1u;
        const uint32_t *lo_lv  = mid_lv + m_cnt;
        qua_bitop_onehot_compute(kb, addr_msb + 1u, nm, mid_lv);
        qua_bitop_onehot_compute(kb, addr_msb + 1u + nm, nl, lo_lv);
        for (uint32_t m = 0u; m < m_cnt; ++m) {
            uint32_t a0_base = m << nl;
            uint32_t a1_base = (1u << rem) | (m << nl);
            if (a0_base >= n_entries) break;
            qua_x(kb, addr_msb[0]);
            qua_and_compute(kb, addr_msb[0], mid_lv[m], anc_m);
            qua_x(kb, addr_msb[0]);
            for (uint32_t l = 0u; l < l_cnt && (a0_base + l) < n_entries; ++l) {
                bool any = false;
                qua_bit_store0(kb, pbit);
                uint64_t val = table[a0_base + l];
                for (uint32_t j = 0u; j < n_tgt; ++j) {
                    if ((val >> j) & 1ull) {
                        qua_bit_invert_if(kb, pbit, hmr_bits[j]);
                        any = true;
                    }
                }
                if (any) qua_cz_if(kb, anc_m, lo_lv[l], pbit);
            }
            if (a1_base < n_entries) {
                qua_cx(kb, mid_lv[m], anc_m);
                for (uint32_t l = 0u; l < l_cnt && (a1_base + l) < n_entries; ++l) {
                    bool any = false;
                    qua_bit_store0(kb, pbit);
                    uint64_t val = table[a1_base + l];
                    for (uint32_t j = 0u; j < n_tgt; ++j) {
                        if ((val >> j) & 1ull) {
                            qua_bit_invert_if(kb, pbit, hmr_bits[j]);
                            any = true;
                        }
                    }
                    if (any) qua_cz_if(kb, anc_m, lo_lv[l], pbit);
                }
                qua_and_uncompute(kb, addr_msb[0], mid_lv[m], anc_m, (*next_cbit)++);
            } else {
                qua_x(kb, addr_msb[0]);
                qua_and_uncompute(kb, addr_msb[0], mid_lv[m], anc_m, (*next_cbit)++);
                qua_x(kb, addr_msb[0]);
            }
        }
        qua_bitop_onehot_uncompute_mbuc(kb, addr_msb + 1u + nm, nl, lo_lv, next_cbit);
        qua_bitop_onehot_uncompute_mbuc(kb, addr_msb + 1u, nm, mid_lv, next_cbit);
        return;
    }

    /* p in {2, 3}: prefix AND chain on addr_msb[0..p-1] + anc_m + (mid_lv, lo_lv) */
    const uint32_t *pref_anc = leaves_buf;          /* length p - 1 */
    uint32_t        anc_m    = leaves_buf[p - 1u];
    const uint32_t *mid_lv   = leaves_buf + p;
    const uint32_t *lo_lv    = mid_lv + m_cnt;
    qua_bitop_onehot_compute(kb, addr_msb + p, nm, mid_lv);
    qua_bitop_onehot_compute(kb, addr_msb + p + nm, nl, lo_lv);
    uint32_t top_cnt = 1u << p;
    for (uint32_t tv = 0u; tv < top_cnt; ++tv) {
        uint32_t q_top = tv << rem;
        if (q_top >= n_entries) break;
        for (uint32_t bi = 0u; bi < p; ++bi) {
            if (((tv >> (p - 1u - bi)) & 1u) == 0u) qua_x(kb, addr_msb[bi]);
        }
        qua_and_compute(kb, addr_msb[0], addr_msb[1], pref_anc[0]);
        for (uint32_t bi = 2u; bi < p; ++bi) {
            qua_and_compute(kb, pref_anc[bi - 2u], addr_msb[bi], pref_anc[bi - 1u]);
        }
        uint32_t p_last = pref_anc[p - 2u];
        for (uint32_t m = 0u; m < m_cnt; ++m) {
            uint32_t q_base = q_top | (m << nl);
            if (q_base >= n_entries) break;
            qua_and_compute(kb, p_last, mid_lv[m], anc_m);
            for (uint32_t l = 0u; l < l_cnt && (q_base + l) < n_entries; ++l) {
                bool any = false;
                qua_bit_store0(kb, pbit);
                uint64_t val = table[q_base + l];
                for (uint32_t j = 0u; j < n_tgt; ++j) {
                    if ((val >> j) & 1ull) {
                        qua_bit_invert_if(kb, pbit, hmr_bits[j]);
                        any = true;
                    }
                }
                if (any) qua_cz_if(kb, anc_m, lo_lv[l], pbit);
            }
            qua_and_uncompute(kb, p_last, mid_lv[m], anc_m, (*next_cbit)++);
        }
        for (uint32_t bi = p - 1u; bi >= 2u; --bi) {
            qua_and_uncompute(kb, pref_anc[bi - 2u], addr_msb[bi], pref_anc[bi - 1u], (*next_cbit)++);
        }
        qua_and_uncompute(kb, addr_msb[0], addr_msb[1], pref_anc[0], (*next_cbit)++);
        for (uint32_t bi = p; bi > 0u; --bi) {
            if (((tv >> (p - bi)) & 1u) == 0u) qua_x(kb, addr_msb[bi - 1u]);
        }
    }
    qua_bitop_onehot_uncompute_mbuc(kb, addr_msb + p + nm, nl, lo_lv, next_cbit);
    qua_bitop_onehot_uncompute_mbuc(kb, addr_msb + p, nm, mid_lv, next_cbit);
}

/* --- Combinator 2: LinF2_Fuse (GL_w(F_2) Composition + Min-Fanout In-Place Synthesis) --- */
static inline void qua_bitop_compose_linf2_u64(
    const uint64_t *B_rows,
    const uint64_t *A_rows,
    uint64_t *out_BA_rows,
    uint32_t w
) {
    uint64_t tmp[64];
    for (uint32_t r = 0u; r < w; ++r) {
        uint64_t brow = B_rows[r];
        uint64_t acc  = 0ull;
        for (uint32_t k = 0u; k < w; ++k) {
            if ((brow >> k) & 1ull) {
                for (uint32_t c = 0u; c < w; ++c) {
                    if ((A_rows[k] >> c) & 1ull) acc ^= (UINT64_C(1) << c);
                }
            }
        }
        tmp[r] = acc;
    }
    for (uint32_t r = 0u; r < w; ++r) out_BA_rows[r] = tmp[r];
}

static inline void qua_bitop_emit_linf2_inplace_u64(
    qua_builder_t *kb,
    const uint32_t *q,
    uint32_t w,
    const uint64_t *mat_rows
) {
    if (w == 0u || w > 64u) return;
    uint64_t A[64];
    for (uint32_t r = 0u; r < w; ++r) A[r] = mat_rows[r];

    typedef struct { uint8_t kind, a, b; } linf2_op_t;
    linf2_op_t ops[4096];
    uint32_t nops = 0u;

    for (uint32_t col = 0u; col < w; ++col) {
        uint32_t best_piv = col;
        uint32_t best_cost = UINT32_MAX;
        for (uint32_t r = col; r < w; ++r) {
            if (((A[r] >> col) & 1ull) == 0ull) continue;
            uint32_t cost = (uint32_t)__builtin_popcountll((unsigned long long)A[r]) + ((r != col) ? 2u : 0u);
            if (cost < best_cost) {
                best_cost = cost;
                best_piv  = r;
            }
        }
        if (((A[best_piv] >> col) & 1ull) == 0ull) continue;
        if (best_piv != col) {
            if (nops >= 4096u) abort();
            uint64_t t = A[col]; A[col] = A[best_piv]; A[best_piv] = t;
            ops[nops++] = (linf2_op_t){1u, (uint8_t)col, (uint8_t)best_piv};
        }
        for (uint32_t r = 0u; r < w; ++r) {
            if (r != col && ((A[r] >> col) & 1ull)) {
                if (nops >= 4096u) abort();
                A[r] ^= A[col];
                ops[nops++] = (linf2_op_t){0u, (uint8_t)r, (uint8_t)col};
            }
        }
    }

    for (uint32_t k = nops; k > 0u; --k) {
        linf2_op_t op = ops[k - 1u];
        if (op.kind == 1u) {
            qua_cx(kb, q[op.a], q[op.b]);
            qua_cx(kb, q[op.b], q[op.a]);
            qua_cx(kb, q[op.a], q[op.b]);
        } else {
            qua_cx(kb, q[op.b], q[op.a]);
        }
    }
}

/* --- Combinator 3: Cuccaro 1-Ancilla Controlled Adder + PhaseSplit DAG Router ---
 * Computes `if (ctrl) tgt[0..n-1] +=/-= lhs[0..k-1] (mod 2^n)` using in-place
 * MAJ + CTRL_UMA on the bottom `k` bits (1 carry-in ancilla instead of `k`!). */
static inline void qua_bitop_ctrl_iadd_cuccaro_resources(
    uint32_t k_in,
    uint32_t n_tgt,
    uint32_t *out_anc,
    uint64_t *out_ccx
) {
    if (k_in == 0u || n_tgt == 0u) {
        *out_anc = 0u;
        *out_ccx = 0u;
        return;
    }
    uint32_t k = (k_in < n_tgt) ? k_in : n_tgt;
    if (k == 1u) {
        if (n_tgt == 1u) { *out_anc = 0u; *out_ccx = 1ull; }
        else             { *out_anc = n_tgt - 1u; *out_ccx = (uint64_t)n_tgt; }
        return;
    }
    if (k == n_tgt) {
        *out_anc = 1u;
        *out_ccx = 3ull * (uint64_t)n_tgt - 2ull;
        return;
    }
    uint32_t tail = n_tgt - k;
    *out_anc = (tail == 1u) ? 1u : tail;
    *out_ccx = 2ull * (uint64_t)k + (uint64_t)n_tgt;
}

static inline void qua_bitop_emit_ctrl_iadd_cuccaro(
    qua_builder_t *kb,
    uint32_t ctrl,
    const uint32_t *lhs,
    uint32_t k_in,
    const uint32_t *tgt,
    uint32_t n_tgt,
    bool subtract,
    const uint32_t *anc,
    uint32_t *next_cbit
) {
    if (k_in == 0u || n_tgt == 0u) return;
    uint32_t k = (k_in < n_tgt) ? k_in : n_tgt;
    if (subtract) {
        for (uint32_t i = 0u; i < n_tgt; ++i) qua_x(kb, tgt[i]);
    }
    if (k == 1u) {
        if (n_tgt == 1u) {
            qua_ccx(kb, ctrl, lhs[0], tgt[0]);
        } else {
            uint32_t p0 = anc[0];
            qua_and_compute(kb, ctrl, lhs[0], p0);
            if (n_tgt == 2u) {
                qua_ccx(kb, p0, tgt[0], tgt[1]);
                qua_cx(kb, p0, tgt[0]);
            } else {
                const uint32_t *carries = anc + 1u; /* length n_tgt - 2 */
                uint32_t m = n_tgt;
                qua_and_compute(kb, p0, tgt[0], carries[0]);
                for (uint32_t i = 1u; i < m - 2u; ++i) {
                    qua_and_compute(kb, carries[i - 1u], tgt[i], carries[i]);
                }
                qua_ccx(kb, carries[m - 3u], tgt[m - 2u], tgt[m - 1u]);
                for (uint32_t i = m - 3u; i >= 1u; --i) {
                    qua_and_uncompute(kb, carries[i - 1u], tgt[i], carries[i], (*next_cbit)++);
                    qua_cx(kb, carries[i - 1u], tgt[i]);
                }
                qua_and_uncompute(kb, p0, tgt[0], carries[0], (*next_cbit)++);
                qua_cx(kb, p0, tgt[0]);
            }
            qua_and_uncompute(kb, ctrl, lhs[0], p0, (*next_cbit)++);
        }
        if (subtract) {
            for (uint32_t i = 0u; i < n_tgt; ++i) qua_x(kb, tgt[i]);
        }
        return;
    }

    uint32_t c0 = anc[0];
    uint32_t maj_hi = (k == n_tgt) ? (n_tgt - 1u) : k;
    for (uint32_t i = 0u; i < maj_hi; ++i) {
        uint32_t cin = (i == 0u) ? c0 : lhs[i - 1u];
        qua_cx(kb, lhs[i], tgt[i]);
        qua_cx(kb, lhs[i], cin);
        qua_ccx(kb, cin, tgt[i], lhs[i]);
    }

    if (k == n_tgt) {
        uint32_t c_last = lhs[n_tgt - 2u];
        qua_cx(kb, c_last, lhs[n_tgt - 1u]);
        qua_ccx(kb, ctrl, lhs[n_tgt - 1u], tgt[n_tgt - 1u]);
        qua_cx(kb, c_last, lhs[n_tgt - 1u]);
    } else {
        uint32_t tail = n_tgt - k;
        uint32_t c_out = lhs[k - 1u];
        if (tail == 1u) {
            qua_ccx(kb, ctrl, c_out, tgt[k]);
        } else {
            uint32_t p_top = anc[1];
            qua_and_compute(kb, ctrl, c_out, p_top);
            if (tail == 2u) {
                qua_ccx(kb, p_top, tgt[k], tgt[k + 1u]);
                qua_cx(kb, p_top, tgt[k]);
            } else {
                const uint32_t *carries = anc + 2u; /* length tail - 2 */
                qua_and_compute(kb, p_top, tgt[k], carries[0]);
                for (uint32_t i = 1u; i < tail - 2u; ++i) {
                    qua_and_compute(kb, carries[i - 1u], tgt[k + i], carries[i]);
                }
                qua_ccx(kb, carries[tail - 3u], tgt[k + tail - 2u], tgt[k + tail - 1u]);
                for (uint32_t i = tail - 3u; i >= 1u; --i) {
                    qua_and_uncompute(kb, carries[i - 1u], tgt[k + i], carries[i], (*next_cbit)++);
                    qua_cx(kb, carries[i - 1u], tgt[k + i]);
                }
                qua_and_uncompute(kb, p_top, tgt[k], carries[0], (*next_cbit)++);
                qua_cx(kb, p_top, tgt[k]);
            }
            qua_and_uncompute(kb, ctrl, c_out, p_top, (*next_cbit)++);
        }
    }

    for (uint32_t ii = maj_hi; ii > 0u; --ii) {
        uint32_t i = ii - 1u;
        uint32_t cin = (i == 0u) ? c0 : lhs[i - 1u];
        qua_ccx(kb, cin, tgt[i], lhs[i]);
        qua_cx(kb, lhs[i], cin);
        qua_cx(kb, cin, lhs[i]);
        qua_ccx(kb, ctrl, lhs[i], tgt[i]);
        qua_cx(kb, cin, lhs[i]);
        qua_cx(kb, lhs[i], tgt[i]);
    }

    if (subtract) {
        for (uint32_t i = 0u; i < n_tgt; ++i) qua_x(kb, tgt[i]);
    }
}

/* --- Combinator 4: Telescoping Multi-Band HGCD Window Verifier --- */
static inline bool qua_bitop_verify_hgcd_bands(
    uint32_t n_syn,
    uint32_t n_disjoint_lo,
    uint32_t n_stages,
    const uint32_t *r_end,
    const uint32_t *s_lo
) {
    if (n_stages == 0u || n_stages > 8u) return false;
    if (s_lo[0] > n_disjoint_lo && n_stages > 1u) return false;
    for (uint32_t p = 0u; p + 1u < n_stages; ++p) {
        if (r_end[p] == 0u || r_end[p] >= r_end[p + 1u]) return false;
        if (2u * r_end[p] > n_syn) return false;
        if (s_lo[p] > n_syn - 2u * r_end[p]) return false;
        if (p + 2u < n_stages && s_lo[p + 1u] >= s_lo[p]) return false;
    }
    return true;
}

/* --- Combinator 5: BasisIso Hamming-TSP Row Co-Permuter ---
 * Finds a permutation `perm[0..n_rows-1]` minimizing total adjacent Hamming
 * distance across rows of `words_per_row` uint64_t words via nearest-neighbor
 * + 2-opt local search (100% deterministic). */
static inline uint32_t qua_bitop_row_hamming_u64(
    const uint64_t *a,
    const uint64_t *b,
    uint32_t words_per_row
) {
    uint32_t d = 0u;
    for (uint32_t w = 0u; w < words_per_row; ++w) {
        d += (uint32_t)__builtin_popcountll((unsigned long long)(a[w] ^ b[w]));
    }
    return d;
}

static inline void qua_bitop_tsp_gray_order_u64(
    const uint64_t *rows,
    uint32_t n_rows,
    uint32_t words_per_row,
    uint32_t *out_perm
) {
    if (n_rows == 0u) return;
    for (uint32_t i = 0u; i < n_rows; ++i) out_perm[i] = i;
    if (n_rows <= 2u) return;

    uint8_t *used = (uint8_t *)calloc(n_rows, 1u);
    if (!used) return;

    /* Start from minimum Hamming-weight row to anchor Virtual-Slot-0 near 0 */
    uint32_t start = 0u, min_w = UINT32_MAX;
    for (uint32_t i = 0u; i < n_rows; ++i) {
        uint32_t w = 0u;
        for (uint32_t k = 0u; k < words_per_row; ++k) {
            w += (uint32_t)__builtin_popcountll((unsigned long long)rows[(size_t)i * words_per_row + k]);
        }
        if (w < min_w) { min_w = w; start = i; }
    }

    out_perm[0] = start;
    used[start] = 1u;
    for (uint32_t step = 1u; step < n_rows; ++step) {
        const uint64_t *prev = rows + (size_t)out_perm[step - 1u] * words_per_row;
        uint32_t best_j = 0u, best_d = UINT32_MAX;
        for (uint32_t j = 0u; j < n_rows; ++j) {
            if (used[j]) continue;
            uint32_t d = qua_bitop_row_hamming_u64(prev, rows + (size_t)j * words_per_row, words_per_row);
            if (d < best_d) { best_d = d; best_j = j; }
        }
        used[best_j] = 1u;
        out_perm[step] = best_j;
    }
    free(used);

    /* Open-path 2-opt refinement */
    uint32_t max_passes = (n_rows <= 256u) ? 6u : 2u;
    for (uint32_t pass = 0u; pass < max_passes; ++pass) {
        bool improved = false;
        for (uint32_t i = 1u; i + 1u < n_rows; ++i) {
            const uint64_t *ri_m1 = rows + (size_t)out_perm[i - 1u] * words_per_row;
            const uint64_t *ri    = rows + (size_t)out_perm[i] * words_per_row;
            uint32_t d_left = qua_bitop_row_hamming_u64(ri_m1, ri, words_per_row);
            for (uint32_t k = i + 1u; k < n_rows; ++k) {
                const uint64_t *rk = rows + (size_t)out_perm[k] * words_per_row;
                uint32_t d_right = (k + 1u < n_rows)
                    ? qua_bitop_row_hamming_u64(rk, rows + (size_t)out_perm[k + 1u] * words_per_row, words_per_row)
                    : 0u;
                uint32_t new_left  = qua_bitop_row_hamming_u64(ri_m1, rk, words_per_row);
                uint32_t new_right = (k + 1u < n_rows)
                    ? qua_bitop_row_hamming_u64(ri, rows + (size_t)out_perm[k + 1u] * words_per_row, words_per_row)
                    : 0u;
                if (new_left + new_right < d_left + d_right) {
                    uint32_t lo = i, hi = k;
                    while (lo < hi) {
                        uint32_t tmp = out_perm[lo];
                        out_perm[lo] = out_perm[hi];
                        out_perm[hi] = tmp;
                        ++lo; --hi;
                    }
                    improved = true;
                    ri     = rows + (size_t)out_perm[i] * words_per_row;
                    d_left = qua_bitop_row_hamming_u64(ri_m1, ri, words_per_row);
                }
            }
        }
        if (!improved) break;
    }
}

#endif /* QUA_BITOP_H */
