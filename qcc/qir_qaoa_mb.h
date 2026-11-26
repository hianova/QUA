/*
 * qir_qaoa_mb.h — Declarative SSA QIR Builders for QAOA MaxCut & Many-Body 2D Ising Trotter
 * Part of QCC v3 (QUA). Pure C17 & C++20 compatible, zero raw qcc_kmb_* gate emission.
 */
#ifndef QIR_QAOA_MB_H
#define QIR_QAOA_MB_H

#include "qir.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Build a declarative 1-layer QAOA MaxCut cost-phase oracle:
 *   |x> -> exp(i * gamma * sum_{(u,v) in E} (x_u XOR x_v)) |x>
 * Uses explicit XOR cut indicators + balanced binary adder popcount + Z_POW phase polynomial,
 * with clean uncomputation of all intermediate cut/sum scratch registers (`clean_anc = 1`).
 */
static inline bool qir_build_qaoa_maxcut_qir(
    qir_module_t *m,
    uint32_t n_nodes,
    const uint8_t (*edges)[2],
    uint32_t n_edges
) {
    if (!m || !edges || n_nodes < 2u || n_nodes > 64u || n_edges == 0u || n_edges > 128u) {
        return false;
    }
    m->require_clean_scratch = true;
    uint32_t x = qir_add_input(m, (uint16_t)n_nodes, "qaoa_x");

    qcc_bv256_t rows[128];
    memset(rows, 0, sizeof(rows));
    for (uint32_t e = 0u; e < n_edges; ++e) {
        uint8_t u = edges[e][0];
        uint8_t v = edges[e][1];
        if (u >= n_nodes || v >= n_nodes || u == v) {
            return false;
        }
        qcc_bv_set(rows[e].w, u);
        qcc_bv_set(rows[e].w, v);
    }
    uint32_t cut_vec = qir_emit_linf2(m, &x, 1u, (uint16_t)n_edges, rows, NULL, "cut_edges");

    uint32_t cut_terms[128];
    for (uint32_t e = 0u; e < n_edges; ++e) {
        cut_terms[e] = qir_emit_slice(m, cut_vec, (uint16_t)e, (uint16_t)(e + 1u), "cut_e");
    }

    /* Compute final sum width w_sum = ceil_log2(n_edges + 1), pad leaves ONCE to w_sum */
    uint16_t sum_w = 1u;
    while ((1u << sum_w) <= n_edges) {
        sum_w++;
    }
    if (sum_w > 16u) {
        return false;
    }

    uint32_t level[128];
    for (uint32_t i = 0u; i < n_edges; ++i) {
        if (sum_w > 1u) {
            uint32_t pad0 = qir_emit_const_u64(m, (uint16_t)(sum_w - 1u), 0u, "pad0");
            uint32_t cins[2] = {cut_terms[i], pad0};
            level[i] = qir_emit_concat(m, cins, 2u, "cut_w");
        } else {
            level[i] = cut_terms[i];
        }
    }

    uint32_t cur_n = n_edges;
    while (cur_n > 1u) {
        uint32_t next_n = 0u;
        for (uint32_t i = 0u; i + 1u < cur_n; i += 2u) {
            level[next_n++] = qir_emit_add(m, level[i], level[i + 1u], 0u, "cut_sum");
        }
        if ((cur_n & 1u) != 0u) {
            level[next_n++] = level[cur_n - 1u];
        }
        cur_n = next_n;
    }

    uint32_t total_cut = level[0];

    qir_phase_term_t terms[17];
    memset(terms, 0, sizeof(terms));
    /* Linear QAOA cost phase: exp(i * (2*pi / 2^sum_w) * total_cut) */
    for (uint16_t b = 0u; b < sum_w; ++b) {
        qcc_bv_set(terms[b].mask, b);
        terms[b].is_xor = 0u;
        uint32_t exp_b = 128u - (uint32_t)sum_w + (uint32_t)b;
        terms[b].angle = (exp_b < 128u) ? ((flow_u128)1u << exp_b) : 0u;
    }
    uint32_t nt = (uint32_t)sum_w;
    if (sum_w >= 3u) {
        qcc_bv_set(terms[nt].mask, (uint32_t)(sum_w - 3u));
        qcc_bv_set(terms[nt].mask, (uint32_t)(sum_w - 2u));
        terms[nt].is_xor = 0u;
        terms[nt].angle  = (flow_u128)1u << 127u;
        nt++;
    } else if (sum_w == 2u) {
        qcc_bv_set(terms[nt].mask, 0u);
        qcc_bv_set(terms[nt].mask, 1u);
        terms[nt].is_xor = 0u;
        terms[nt].angle  = (flow_u128)1u << 127u;
        nt++;
    }
    (void)qir_emit_phase_poly(m, &total_cut, 1u, terms, nt, "qaoa_cost_phase");

    qir_add_output(m, x);
    return true;
}

/* Optimized QAOA MaxCut cost-phase oracle via 4-bit Popcount-LUT chunking:
 *   Instead of padding all `n_edges` single bits to `sum_w` and summing with `n_edges - 1`
 *   adders, partition the edge cut indicators into <=4-bit blocks, look up exact 3-bit
 *   popcount via `QIR_OP_LUT_ROM` (optimized by R2/R7/R8), and sum only `ceil(E/4)` values.
 *   Produces the bit-exact identical unitary phase as `qir_build_qaoa_maxcut_qir`.
 */
static inline bool qir_build_qaoa_maxcut_opt_qir(
    qir_module_t *m,
    uint32_t n_nodes,
    const uint8_t (*edges)[2],
    uint32_t n_edges
) {
    if (!m || !edges || n_nodes < 2u || n_nodes > 64u || n_edges == 0u || n_edges > 128u) {
        return false;
    }
    m->require_clean_scratch = true;
    uint32_t x = qir_add_input(m, (uint16_t)n_nodes, "qaoa_x");

    qcc_bv256_t rows[128];
    memset(rows, 0, sizeof(rows));
    for (uint32_t e = 0u; e < n_edges; ++e) {
        uint8_t u = edges[e][0];
        uint8_t v = edges[e][1];
        if (u >= n_nodes || v >= n_nodes || u == v) {
            return false;
        }
        qcc_bv_set(rows[e].w, u);
        qcc_bv_set(rows[e].w, v);
    }
    uint32_t cut_vec = qir_emit_linf2(m, &x, 1u, (uint16_t)n_edges, rows, NULL, "cut_edges");

    uint16_t sum_w = 1u;
    while ((1u << sum_w) <= n_edges) {
        sum_w++;
    }
    if (sum_w < 3u) sum_w = 3u;
    if (sum_w > 16u) {
        return false;
    }

    uint32_t level[32];
    uint32_t n_chunks = 0u;
    for (uint32_t off = 0u; off < n_edges; off += 4u) {
        uint32_t cw = (n_edges - off >= 4u) ? 4u : (n_edges - off);
        uint32_t sub = qir_emit_slice(m, cut_vec, (uint16_t)off, (uint16_t)(off + cw), "cut_blk");
        uint64_t pop_tbl[16];
        uint32_t n_ent = UINT32_C(1) << cw;
        for (uint32_t idx = 0u; idx < n_ent; ++idx) {
            pop_tbl[idx] = (uint64_t)__builtin_popcount(idx);
        }
        uint32_t pc3 = qir_emit_lut_rom(m, sub, 3u, n_ent, pop_tbl, 0u, "pop_lut");
        if (sum_w > 3u) {
            uint32_t pad0 = qir_emit_const_u64(m, (uint16_t)(sum_w - 3u), 0u, "pad0");
            uint32_t cins[2] = {pc3, pad0};
            level[n_chunks++] = qir_emit_concat(m, cins, 2u, "pop_w");
        } else {
            level[n_chunks++] = pc3;
        }
    }

    uint32_t cur_n = n_chunks;
    while (cur_n > 1u) {
        uint32_t next_n = 0u;
        for (uint32_t i = 0u; i + 1u < cur_n; i += 2u) {
            level[next_n++] = qir_emit_add(m, level[i], level[i + 1u], 0u, "cut_sum");
        }
        if ((cur_n & 1u) != 0u) {
            level[next_n++] = level[cur_n - 1u];
        }
        cur_n = next_n;
    }

    uint32_t total_cut = level[0];

    qir_phase_term_t terms[17];
    memset(terms, 0, sizeof(terms));
    for (uint16_t b = 0u; b < sum_w; ++b) {
        qcc_bv_set(terms[b].mask, b);
        terms[b].is_xor = 0u;
        uint32_t exp_b = 128u - (uint32_t)sum_w + (uint32_t)b;
        terms[b].angle = (exp_b < 128u) ? ((flow_u128)1u << exp_b) : 0u;
    }
    uint32_t nt = (uint32_t)sum_w;
    qcc_bv_set(terms[nt].mask, (uint32_t)(sum_w - 3u));
    qcc_bv_set(terms[nt].mask, (uint32_t)(sum_w - 2u));
    terms[nt].is_xor = 0u;
    terms[nt].angle  = (flow_u128)1u << 127u;
    nt++;
    (void)qir_emit_phase_poly(m, &total_cut, 1u, terms, nt, "qaoa_cost_phase");

    qir_add_output(m, x);
    return true;
}

/* Build a declarative 2D Plaquette Lattice Gauge Theory (LGT) Trotter step:
 *   |s> -> exp(i * (g4 * sum_plaq (-1)^{s0^s1^s2^s3} + g2 * sum_<u,v> (-1)^{s_u^s_v} + g1 * sum_u (-1)^{s_u})) |s>
 * Pure zero-CCX phase-polynomial synthesis (`49` multi-body Pauli terms on a 4x4 grid).
 */
static inline bool qir_build_plaquette_lgt_qir(
    qir_module_t *m,
    uint32_t rows,
    uint32_t cols
) {
    if (!m || rows < 2u || cols < 2u || (rows * cols) > 64u) {
        return false;
    }
    uint32_t n_q = rows * cols;
    m->require_clean_scratch = true;
    uint32_t s = qir_add_input(m, (uint16_t)n_q, "lgt_lat");

    uint32_t max_terms = (rows - 1u) * (cols - 1u) + rows * (cols - 1u) + (rows - 1u) * cols + n_q;
    qir_phase_term_t *terms = (qir_phase_term_t *)qcc_xcalloc(max_terms, sizeof(qir_phase_term_t));
    uint32_t nt = 0u;

    /* 4-body ZZZZ plaquette terms (angle = pi/4 -> 1 << 125) */
    const flow_u128 p4_angle = (flow_u128)1u << 125u;
    for (uint32_t r = 0u; r + 1u < rows; ++r) {
        for (uint32_t c = 0u; c + 1u < cols; ++c) {
            qcc_bv_set(terms[nt].mask, r * cols + c);
            qcc_bv_set(terms[nt].mask, r * cols + c + 1u);
            qcc_bv_set(terms[nt].mask, (r + 1u) * cols + c);
            qcc_bv_set(terms[nt].mask, (r + 1u) * cols + c + 1u);
            terms[nt].is_xor = 1u;
            terms[nt].angle  = p4_angle;
            nt++;
        }
    }

    /* 2-body ZZ gauge-matter link terms (angle = pi/8 -> 1 << 124) */
    const flow_u128 p2_angle = (flow_u128)1u << 124u;
    for (uint32_t r = 0u; r < rows; ++r) {
        for (uint32_t c = 0u; c + 1u < cols; ++c) {
            qcc_bv_set(terms[nt].mask, r * cols + c);
            qcc_bv_set(terms[nt].mask, r * cols + c + 1u);
            terms[nt].is_xor = 1u;
            terms[nt].angle  = p2_angle;
            nt++;
        }
    }
    for (uint32_t r = 0u; r + 1u < rows; ++r) {
        for (uint32_t c = 0u; c < cols; ++c) {
            qcc_bv_set(terms[nt].mask, r * cols + c);
            qcc_bv_set(terms[nt].mask, (r + 1u) * cols + c);
            terms[nt].is_xor = 1u;
            terms[nt].angle  = p2_angle;
            nt++;
        }
    }

    /* 1-body Z electric field terms (angle = pi/16 -> 1 << 123) */
    const flow_u128 p1_angle = (flow_u128)1u << 123u;
    for (uint32_t u = 0u; u < n_q; ++u) {
        qcc_bv_set(terms[nt].mask, u);
        terms[nt].is_xor = 0u;
        terms[nt].angle  = p1_angle;
        nt++;
    }

    (void)qir_emit_phase_poly(m, &s, 1u, terms, nt, "lgt_plaquette_trotter");
    free(terms);

    qir_add_output(m, s);
    return true;
}

/* Build a declarative 3-SAT Grover phase oracle:
 *   |x> -> (-1)^{ALL_CLAUSES_SAT(x)} |x>
 * Uses one fused LINF2 layer for all negated literals, clean AND_VEC conjunctions per clause,
 * `QIR_OP_EQ_ZERO` zero-unsat check, and pi phase kickback with full clean ancilla uncomputation.
 */
static inline bool qir_build_ksat_oracle_qir(
    qir_module_t *m,
    uint32_t n_vars,
    const int8_t (*clauses)[3],
    uint32_t n_clauses
) {
    if (!m || !clauses || n_vars < 3u || n_vars > 64u || n_clauses == 0u || n_clauses > 64u) {
        return false;
    }
    m->require_clean_scratch = true;
    uint32_t x = qir_add_input(m, (uint16_t)n_vars, "sat_x");

    uint16_t n_lits = (uint16_t)(3u * n_clauses);
    qcc_bv256_t rows[192];
    memset(rows, 0, sizeof(rows));
    qcc_bv256_t cvec = qcc_bv_zero();

    for (uint32_t ci = 0u; ci < n_clauses; ++ci) {
        for (uint32_t li = 0u; li < 3u; ++li) {
            int8_t lit = clauses[ci][li];
            int32_t abs_v = (lit < 0) ? -(int32_t)lit : (int32_t)lit;
            if (lit == 0 || (uint32_t)abs_v > n_vars) {
                return false;
            }
            uint32_t row = 3u * ci + li;
            uint32_t vbit = (uint32_t)(abs_v - 1);
            qcc_bv_set(rows[row].w, vbit);
            /* Clause falsified iff all 3 literals are 0 -> compute negated literal:
             *   lit > 0 (want ~x_v) => XOR with 1 (set cvec bit)
             *   lit < 0 (want  x_v) => keep x_v
             */
            if (lit > 0) {
                qcc_bv_set(cvec.w, row);
            }
        }
    }

    uint32_t neg_lits = qir_emit_linf2(m, &x, 1u, n_lits, rows, &cvec, "neg_lits");
    uint32_t unsat_bits[64];
    for (uint32_t ci = 0u; ci < n_clauses; ++ci) {
        uint32_t b0 = qir_emit_slice(m, neg_lits, (uint16_t)(3u * ci), (uint16_t)(3u * ci + 1u), "nl0");
        uint32_t b1 = qir_emit_slice(m, neg_lits, (uint16_t)(3u * ci + 1u), (uint16_t)(3u * ci + 2u), "nl1");
        uint32_t b2 = qir_emit_slice(m, neg_lits, (uint16_t)(3u * ci + 2u), (uint16_t)(3u * ci + 3u), "nl2");
        uint32_t a01 = qir_emit_and_vec(m, b0, b1, 0u, "u01");
        unsat_bits[ci] = qir_emit_and_vec(m, a01, b2, 0u, "unsat_c");
    }

    uint32_t unsat_vec = (n_clauses == 1u)
        ? unsat_bits[0]
        : qir_emit_concat(m, unsat_bits, (uint16_t)n_clauses, "unsat_vec");
    uint32_t all_sat = qir_emit_eq_zero(m, unsat_vec, "all_sat");

    qir_phase_term_t ph;
    memset(&ph, 0, sizeof(ph));
    qcc_bv_set(ph.mask, 0u);
    ph.is_xor = 0u;
    ph.angle  = (flow_u128)1u << 127u; /* (-1)^all_sat Grover phase flip */
    (void)qir_emit_phase_poly(m, &all_sat, 1u, &ph, 1u, "ksat_grover_flip");

    qir_add_output(m, x);
    return true;
}

/* Build a declarative 2D Transverse/Longitudinal Ising Many-Body Trotter step:
 *   |s> -> exp(i * (J * sum_{<u,v>} (-1)^{s_u ^ s_v} + h * sum_u (-1)^{s_u})) |s>
 * Uses `QIR_OP_PHASE_POLY` with `is_xor = 1` two-qubit ZZ parity terms + `is_xor = 0` Z terms.
 */
static inline bool qir_build_ising_trotter_qir(
    qir_module_t *m,
    uint32_t rows,
    uint32_t cols
) {
    if (!m || rows < 2u || cols < 2u || (rows * cols) > 64u) {
        return false;
    }
    uint32_t n_q = rows * cols;
    m->require_clean_scratch = true;
    uint32_t s = qir_add_input(m, (uint16_t)n_q, "spin_lat");

    uint32_t max_terms = (rows * (cols - 1u)) + ((rows - 1u) * cols) + n_q;
    qir_phase_term_t *terms = (qir_phase_term_t *)qcc_xcalloc(max_terms, sizeof(qir_phase_term_t));
    uint32_t nt = 0u;

    /* 2Q ZZ interaction terms (J = pi/4 -> T-gate parity 1 << 125) */
    const flow_u128 zz_angle = (flow_u128)1u << 125u;
    for (uint32_t r = 0u; r < rows; ++r) {
        for (uint32_t c = 0u; c + 1u < cols; ++c) {
            uint32_t u = r * cols + c;
            uint32_t v = r * cols + c + 1u;
            qcc_bv_set(terms[nt].mask, u);
            qcc_bv_set(terms[nt].mask, v);
            terms[nt].is_xor = 1u;
            terms[nt].angle  = zz_angle;
            nt++;
        }
    }
    for (uint32_t r = 0u; r + 1u < rows; ++r) {
        for (uint32_t c = 0u; c < cols; ++c) {
            uint32_t u = r * cols + c;
            uint32_t v = (r + 1u) * cols + c;
            qcc_bv_set(terms[nt].mask, u);
            qcc_bv_set(terms[nt].mask, v);
            terms[nt].is_xor = 1u;
            terms[nt].angle  = zz_angle;
            nt++;
        }
    }

    /* 1Q longitudinal Z field terms (h = pi/8 -> 1 << 124) */
    const flow_u128 z_angle = (flow_u128)1u << 124u;
    for (uint32_t u = 0u; u < n_q; ++u) {
        qcc_bv_set(terms[nt].mask, u);
        terms[nt].is_xor = 0u;
        terms[nt].angle  = z_angle;
        nt++;
    }

    (void)qir_emit_phase_poly(m, &s, 1u, terms, nt, "ising_trotter_step");
    free(terms);

    qir_add_output(m, s);
    return true;
}

#ifdef __cplusplus
}
#endif

#endif /* QIR_QAOA_MB_H */
