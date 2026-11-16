/*
 * qir_kernels.h — Unified Fault-Tolerant Kernels (Emit-or-Count Zero-Discrepancy Core)
 * Part of QCC (Quantum Circuit Compiler) for FLOW v7.0.1.
 * Pure C17, zero external dependencies, zero indirect function calls.
 *
 * Every kernel takes (qcc_kmb_circuit_t *kb_or_null, qir_kernel_cost_t *cost, ...).
 * Passing kb_or_null == NULL performs exact control-flow cost enumeration so
 * estimated Toffoli/Z_POW/HMR/Ancilla counts equal emitted counts identically.
 */
#ifndef QIR_KERNELS_H
#define QIR_KERNELS_H

#include "qir_rewrite.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint64_t toffoli;
    uint64_t z_pow;
    uint64_t t_like;
    uint64_t hmr;
    uint64_t cnot;
    uint64_t cz;
    uint64_t r_count;
    uint32_t temp_anc_needed;
    uint32_t scratch_bits_needed;
    uint32_t max_live_temp;
} qir_kernel_cost_t;

/* Low-level Emit-or-Count Primitives */
static inline void qk_x(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t q) {
    if (kb) qcc_kmb_x(kb, q);
    (void)c;
}
static inline void qk_z(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t q) {
    if (kb) qcc_kmb_z(kb, q);
    (void)c;
}
static inline void qk_r(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t q) {
    if (kb) qcc_kmb_r(kb, q);
    if (c) c->r_count++;
}
static inline void qk_cx(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t ctrl, uint32_t tgt) {
    if (ctrl == tgt) return;
    if (kb) qcc_kmb_cx(kb, ctrl, tgt);
    if (c) c->cnot++;
}
static inline void qk_cz(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t q0, uint32_t q1) {
    if (kb) qcc_kmb_cz(kb, q0, q1);
    if (c) c->cz++;
}
static inline void qk_ccx(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t c0, uint32_t c1, uint32_t tgt) {
    if (kb) qcc_kmb_ccx(kb, c0, c1, tgt);
    if (c) c->toffoli++;
}
static inline void qk_ccz(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t q0, uint32_t q1, uint32_t q2) {
    if (kb) qcc_kmb_ccz(kb, q0, q1, q2);
    if (c) c->toffoli++;
}
static inline void qk_cswap(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t ctrl, uint32_t q0, uint32_t q1) {
    if (kb) qcc_kmb_cswap(kb, ctrl, q0, q1);
    if (c) c->toffoli++;
}
static inline void qk_hmr(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t q, uint32_t bit) {
    if (kb) qcc_kmb_hmr(kb, q, bit);
    if (c) c->hmr++;
}
static inline void qk_bit_store0(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t bit) {
    if (kb) qcc_kmb_bit_store0(kb, bit);
    (void)c;
}
static inline void qk_bit_store1(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t bit) {
    if (kb) qcc_kmb_bit_store1(kb, bit);
    (void)c;
}
static inline void qk_bit_invert(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t bit) {
    if (kb) qcc_kmb_bit_invert(kb, bit);
    (void)c;
}
static inline void qk_bit_invert_if(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t bit, uint32_t cond) {
    if (kb) {
        qcc_kmb_push_op_raw(kb, QCC_OP_BIT_INVERT_IF);
        qcc_kmb_push_1b(kb, bit);
        qcc_kmb_push_cond(kb, cond);
    }
    (void)c;
}
static inline void qk_neg_if(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t cond) {
    if (kb) qcc_kmb_neg_if(kb, cond);
    (void)c;
}
static inline void qk_z_if(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t q, uint32_t cond) {
    if (kb) qcc_kmb_z_if(kb, q, cond);
    (void)c;
}
static inline void qk_push_cond(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t cond) {
    if (kb) qcc_kmb_push_condition(kb, cond);
    (void)c;
}
static inline void qk_pop_cond(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c) {
    if (kb) qcc_kmb_pop_condition(kb);
    (void)c;
}
static inline void qk_cz_if(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t q0, uint32_t q1, uint32_t cond) {
    if (kb) qcc_kmb_cz_if(kb, q0, q1, cond);
    if (c) c->cz++;
}
static inline void qk_z_pow_u128(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t q, flow_u128 ang) {
    if (ang == 0u) return;
    uint64_t lo = (uint64_t)ang, hi = (uint64_t)(ang >> 64u);
    if (lo == 0u && hi == UINT64_C(0x8000000000000000)) {
        qk_z(kb, c, q);
        return;
    }
    if (kb) qcc_kmb_z_pow(kb, q, lo, hi);
    if (c) {
        c->z_pow++;
        if (qcc_u128_clifford_level(ang) == 3u) c->t_like++;
    }
}
static inline void qk_and_compute(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t c0, uint32_t c1, uint32_t anc) {
    qk_r(kb, c, anc);
    qk_ccx(kb, c, c0, c1, anc);
}
static inline void qk_and_uncompute(qcc_kmb_circuit_t *kb, qir_kernel_cost_t *c, uint32_t c0, uint32_t c1, uint32_t anc, uint32_t sbit) {
    qk_hmr(kb, c, anc, sbit);
    qk_cz_if(kb, c, c0, c1, sbit);
}

/* Scratch index helper when kb == NULL */
static inline uint32_t qk__q(const uint32_t *arr, uint32_t idx) {
    return arr ? arr[idx] : idx;
}

/* ============================================================================
 * Kernel 1: Gidney n-bit In-Place Adder / Subtractor (n-1 ancillas, n-1 CCX)
 * ============================================================================ */
static inline void qk_gidney_iadd(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *lhs,
    const uint32_t *target,
    const uint32_t *anc,
    uint32_t n,
    uint32_t *next_cbit
) {
    if (n == 0u) return;
    if (n == 1u) {
        qk_cx(kb, cost, qk__q(lhs, 0u), qk__q(target, 0u));
        return;
    }
    if (cost && (n - 1u) > cost->temp_anc_needed) cost->temp_anc_needed = n - 1u;
    if (cost && cost->scratch_bits_needed < 1u) cost->scratch_bits_needed = 1u;
    uint32_t sb = next_cbit ? *next_cbit : 0u;

    qk_and_compute(kb, cost, qk__q(lhs, 0u), qk__q(target, 0u), qk__q(anc, 0u));
    for (uint32_t i = 1u; i + 1u < n; ++i) {
        qk_cx(kb, cost, qk__q(anc, i - 1u), qk__q(lhs, i));
        qk_cx(kb, cost, qk__q(anc, i - 1u), qk__q(target, i));
        qk_and_compute(kb, cost, qk__q(lhs, i), qk__q(target, i), qk__q(anc, i));
        qk_cx(kb, cost, qk__q(anc, i - 1u), qk__q(anc, i));
    }
    qk_cx(kb, cost, qk__q(anc, n - 2u), qk__q(target, n - 1u));
    qk_cx(kb, cost, qk__q(lhs, n - 1u), qk__q(target, n - 1u));
    for (uint32_t k = n - 2u; k >= 1u; --k) {
        qk_cx(kb, cost, qk__q(anc, k - 1u), qk__q(anc, k));
        qk_and_uncompute(kb, cost, qk__q(lhs, k), qk__q(target, k), qk__q(anc, k), sb);
        qk_cx(kb, cost, qk__q(anc, k - 1u), qk__q(lhs, k));
        qk_cx(kb, cost, qk__q(lhs, k), qk__q(target, k));
    }
    qk_and_uncompute(kb, cost, qk__q(lhs, 0u), qk__q(target, 0u), qk__q(anc, 0u), sb);
    qk_cx(kb, cost, qk__q(lhs, 0u), qk__q(target, 0u));
    if (next_cbit) (*next_cbit)++;
}

static inline void qk_gidney_iadd_asym(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *lhs,
    uint32_t k_in,
    const uint32_t *target,
    uint32_t n_tgt,
    const uint32_t *anc,
    uint32_t *next_cbit
) {
    if (k_in == 0u || n_tgt == 0u) return;
    uint32_t k = (k_in < n_tgt) ? k_in : n_tgt;
    if (k == n_tgt) {
        qk_gidney_iadd(kb, cost, lhs, target, anc, n_tgt, next_cbit);
        return;
    }
    if (cost && (n_tgt - 1u) > cost->temp_anc_needed) cost->temp_anc_needed = n_tgt - 1u;
    if (cost && cost->scratch_bits_needed < 1u) cost->scratch_bits_needed = 1u;
    uint32_t sb = next_cbit ? *next_cbit : 0u;

    qk_and_compute(kb, cost, qk__q(lhs, 0u), qk__q(target, 0u), qk__q(anc, 0u));
    for (uint32_t i = 1u; i < k && i + 1u < n_tgt; ++i) {
        qk_cx(kb, cost, qk__q(anc, i - 1u), qk__q(lhs, i));
        qk_cx(kb, cost, qk__q(anc, i - 1u), qk__q(target, i));
        qk_and_compute(kb, cost, qk__q(lhs, i), qk__q(target, i), qk__q(anc, i));
        qk_cx(kb, cost, qk__q(anc, i - 1u), qk__q(anc, i));
    }
    for (uint32_t i = k; i + 1u < n_tgt; ++i) {
        qk_and_compute(kb, cost, qk__q(anc, i - 1u), qk__q(target, i), qk__q(anc, i));
    }
    qk_cx(kb, cost, qk__q(anc, n_tgt - 2u), qk__q(target, n_tgt - 1u));
    for (int32_t ii = (int32_t)n_tgt - 2; ii >= (int32_t)k; --ii) {
        uint32_t i = (uint32_t)ii;
        qk_and_uncompute(kb, cost, qk__q(anc, i - 1u), qk__q(target, i), qk__q(anc, i), sb);
        qk_cx(kb, cost, qk__q(anc, i - 1u), qk__q(target, i));
    }
    for (int32_t ii = (int32_t)k - 1; ii >= 1; --ii) {
        uint32_t i = (uint32_t)ii;
        qk_cx(kb, cost, qk__q(anc, i - 1u), qk__q(anc, i));
        qk_and_uncompute(kb, cost, qk__q(lhs, i), qk__q(target, i), qk__q(anc, i), sb);
        qk_cx(kb, cost, qk__q(anc, i - 1u), qk__q(lhs, i));
        qk_cx(kb, cost, qk__q(lhs, i), qk__q(target, i));
    }
    qk_and_uncompute(kb, cost, qk__q(lhs, 0u), qk__q(target, 0u), qk__q(anc, 0u), sb);
    qk_cx(kb, cost, qk__q(lhs, 0u), qk__q(target, 0u));
    if (next_cbit) (*next_cbit)++;
}

static inline void qk_gidney_isub_asym(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *lhs,
    uint32_t k_in,
    const uint32_t *target,
    uint32_t n_tgt,
    const uint32_t *anc,
    uint32_t *next_cbit
) {
    if (k_in == 0u || n_tgt == 0u) return;
    for (uint32_t i = 0u; i < n_tgt; ++i) qk_x(kb, cost, qk__q(target, i));
    qk_gidney_iadd_asym(kb, cost, lhs, k_in, target, n_tgt, anc, next_cbit);
    for (uint32_t i = 0u; i < n_tgt; ++i) qk_x(kb, cost, qk__q(target, i));
}

/* Computes target = (target + lhs + 1) mod 2^n -> when target starts as ~b, target becomes a - b mod 2^n */
static inline void qk_gidney_iadd_p1(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *lhs,
    const uint32_t *target,
    const uint32_t *anc,
    uint32_t n,
    uint32_t *next_cbit
) {
    if (n == 0u) return;
    if (n == 1u) {
        qk_cx(kb, cost, qk__q(lhs, 0u), qk__q(target, 0u));
        qk_x(kb, cost, qk__q(target, 0u));
        return;
    }
    if (cost && (n - 1u) > cost->temp_anc_needed) cost->temp_anc_needed = n - 1u;
    if (cost && cost->scratch_bits_needed < 1u) cost->scratch_bits_needed = 1u;
    uint32_t sb = next_cbit ? *next_cbit : 0u;

    uint32_t l0 = qk__q(lhs, 0u), t0 = qk__q(target, 0u), a0 = qk__q(anc, 0u);
    qk_x(kb, cost, l0);
    qk_x(kb, cost, t0);
    qk_and_compute(kb, cost, l0, t0, a0);
    qk_x(kb, cost, a0);
    qk_x(kb, cost, l0);
    qk_x(kb, cost, t0);

    for (uint32_t i = 1u; i + 1u < n; ++i) {
        qk_cx(kb, cost, qk__q(anc, i - 1u), qk__q(lhs, i));
        qk_cx(kb, cost, qk__q(anc, i - 1u), qk__q(target, i));
        qk_and_compute(kb, cost, qk__q(lhs, i), qk__q(target, i), qk__q(anc, i));
        qk_cx(kb, cost, qk__q(anc, i - 1u), qk__q(anc, i));
    }
    qk_cx(kb, cost, qk__q(anc, n - 2u), qk__q(target, n - 1u));
    qk_cx(kb, cost, qk__q(lhs, n - 1u), qk__q(target, n - 1u));
    for (uint32_t k = n - 2u; k >= 1u; --k) {
        qk_cx(kb, cost, qk__q(anc, k - 1u), qk__q(anc, k));
        qk_and_uncompute(kb, cost, qk__q(lhs, k), qk__q(target, k), qk__q(anc, k), sb);
        qk_cx(kb, cost, qk__q(anc, k - 1u), qk__q(lhs, k));
        qk_cx(kb, cost, qk__q(lhs, k), qk__q(target, k));
    }
    qk_x(kb, cost, a0);
    qk_x(kb, cost, l0);
    qk_x(kb, cost, t0);
    qk_and_uncompute(kb, cost, l0, t0, a0, sb);
    qk_x(kb, cost, l0);
    qk_x(kb, cost, t0);
    qk_cx(kb, cost, l0, t0);
    qk_x(kb, cost, t0);
    if (next_cbit) (*next_cbit)++;
}

/* Out-of-place 3:2 Carry-Save Adder: S = a^b^c, C = (MAJ(a,b,c)<<1) mod 2^w (w-1 CCX, 0 temp anc) */
static inline void qk_cs_add3(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *a,
    const uint32_t *b,
    const uint32_t *c,
    const uint32_t *out_s,
    const uint32_t *out_c,
    uint32_t w
) {
    for (uint32_t i = 0u; i + 1u < w; ++i) {
        uint32_t ai = qk__q(a, i), bi = qk__q(b, i), ci = qk__q(c, i), co = qk__q(out_c, i + 1u);
        qk_r(kb, cost, co);
        qk_cx(kb, cost, ci, ai);
        qk_cx(kb, cost, ci, bi);
        qk_ccx(kb, cost, ai, bi, co);
        qk_cx(kb, cost, ci, bi);
        qk_cx(kb, cost, ci, ai);
        qk_cx(kb, cost, ci, co);
    }
    for (uint32_t i = 0u; i < w; ++i) {
        uint32_t so = qk__q(out_s, i);
        qk_r(kb, cost, so);
        qk_cx(kb, cost, qk__q(a, i), so);
        qk_cx(kb, cost, qk__q(b, i), so);
        qk_cx(kb, cost, qk__q(c, i), so);
    }
    if (w > 0u) qk_r(kb, cost, qk__q(out_c, 0u));
}

/* Uncompute out-of-place 3:2 Carry-Save (reversible) */
static inline void qk_cs_add3_uncompute(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *a,
    const uint32_t *b,
    const uint32_t *c,
    const uint32_t *out_s,
    const uint32_t *out_c,
    uint32_t w,
    uint32_t *next_cbit
) {
    /* Clear out_s directly via CNOTs or HMR */
    for (uint32_t i = 0u; i < w; ++i) {
        qk_cx(kb, cost, qk__q(c, i), qk__q(out_s, i));
        qk_cx(kb, cost, qk__q(b, i), qk__q(out_s, i));
        qk_cx(kb, cost, qk__q(a, i), qk__q(out_s, i));
    }
    /* Clear out_c[i+1] via 0-Toffoli HMR dual of MAJ(a, b, c) = c ^ ((a^c)&(b^c)) */
    uint32_t sb = next_cbit ? *next_cbit : 0u;
    if (cost && w > 1u && cost->scratch_bits_needed < 1u) cost->scratch_bits_needed = 1u;
    for (uint32_t i = 0u; i + 1u < w; ++i) {
        uint32_t ai = qk__q(a, i), bi = qk__q(b, i), ci = qk__q(c, i), co = qk__q(out_c, i + 1u);
        qk_cx(kb, cost, ci, ai);
        qk_cx(kb, cost, ci, bi);
        qk_hmr(kb, cost, co, sb);
        qk_cz_if(kb, cost, ai, bi, sb);
        qk_z_pow_u128(kb, cost, ci, 0u);
        if (kb) qcc_kmb_z_if(kb, ci, sb);
        qk_cx(kb, cost, ci, bi);
        qk_cx(kb, cost, ci, ai);
    }
    if (next_cbit && w > 1u) (*next_cbit)++;
}

/* ============================================================================
 * Kernel 2: Cuccaro 1-Ancilla Controlled Adder (if (ctrl) tgt += lhs mod 2^n)
 * ============================================================================ */
static inline void qk_ctrl_iadd_cuccaro(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    uint32_t ctrl,
    const uint32_t *lhs,
    uint32_t k_in,
    const uint32_t *tgt,
    uint32_t n_tgt,
    const uint32_t *anc,
    uint32_t *next_cbit
) {
    if (k_in == 0u || n_tgt == 0u) return;
    uint32_t k = (k_in < n_tgt) ? k_in : n_tgt;
    uint32_t anc_need = (k == 1u) ? ((n_tgt == 1u) ? 0u : n_tgt - 1u) : ((k == n_tgt) ? 1u : (n_tgt - k));
    if (cost && anc_need > cost->temp_anc_needed) cost->temp_anc_needed = anc_need;
    if (cost && anc_need > 0u && cost->scratch_bits_needed < 1u) cost->scratch_bits_needed = 1u;
    uint32_t sb = next_cbit ? *next_cbit : 0u;

    if (k == 1u) {
        if (n_tgt == 1u) {
            qk_ccx(kb, cost, ctrl, qk__q(lhs, 0u), qk__q(tgt, 0u));
        } else {
            uint32_t p0 = qk__q(anc, 0u);
            qk_and_compute(kb, cost, ctrl, qk__q(lhs, 0u), p0);
            if (n_tgt == 2u) {
                qk_ccx(kb, cost, p0, qk__q(tgt, 0u), qk__q(tgt, 1u));
                qk_cx(kb, cost, p0, qk__q(tgt, 0u));
            } else {
                for (uint32_t i = 0u; i < n_tgt - 2u; ++i) {
                    uint32_t prev = (i == 0u) ? p0 : qk__q(anc, i);
                    qk_and_compute(kb, cost, prev, qk__q(tgt, i), qk__q(anc, 1u + i));
                }
                qk_ccx(kb, cost, qk__q(anc, n_tgt - 2u), qk__q(tgt, n_tgt - 2u), qk__q(tgt, n_tgt - 1u));
                qk_cx(kb, cost, qk__q(anc, n_tgt - 2u), qk__q(tgt, n_tgt - 2u));
                for (uint32_t i = n_tgt - 3u; i >= 1u; --i) {
                    qk_and_uncompute(kb, cost, qk__q(anc, i), qk__q(tgt, i), qk__q(anc, 1u + i), sb);
                    qk_cx(kb, cost, qk__q(anc, i), qk__q(tgt, i));
                }
                qk_and_uncompute(kb, cost, p0, qk__q(tgt, 0u), qk__q(anc, 1u), sb);
                qk_cx(kb, cost, p0, qk__q(tgt, 0u));
            }
            qk_and_uncompute(kb, cost, ctrl, qk__q(lhs, 0u), p0, sb);
        }
        if (next_cbit) (*next_cbit)++;
        return;
    }

    uint32_t c0 = qk__q(anc, 0u);
    qk_r(kb, cost, c0);
    uint32_t maj_hi = (k == n_tgt) ? (n_tgt - 1u) : k;
    for (uint32_t i = 0u; i < maj_hi; ++i) {
        uint32_t cin = (i == 0u) ? c0 : qk__q(lhs, i - 1u);
        qk_cx(kb, cost, qk__q(lhs, i), qk__q(tgt, i));
        qk_cx(kb, cost, qk__q(lhs, i), cin);
        qk_ccx(kb, cost, cin, qk__q(tgt, i), qk__q(lhs, i));
    }

    if (k == n_tgt) {
        uint32_t c_last = qk__q(lhs, n_tgt - 2u);
        qk_cx(kb, cost, c_last, qk__q(lhs, n_tgt - 1u));
        qk_ccx(kb, cost, ctrl, qk__q(lhs, n_tgt - 1u), qk__q(tgt, n_tgt - 1u));
        qk_cx(kb, cost, c_last, qk__q(lhs, n_tgt - 1u));
    } else {
        uint32_t tail = n_tgt - k;
        uint32_t c_out = qk__q(lhs, k - 1u);
        if (tail == 1u) {
            qk_ccx(kb, cost, ctrl, c_out, qk__q(tgt, k));
        } else {
            uint32_t p_top = qk__q(anc, 1u);
            qk_and_compute(kb, cost, ctrl, c_out, p_top);
            if (tail == 2u) {
                qk_ccx(kb, cost, p_top, qk__q(tgt, k), qk__q(tgt, k + 1u));
                qk_cx(kb, cost, p_top, qk__q(tgt, k));
            } else {
                for (uint32_t i = 0u; i < tail - 2u; ++i) {
                    uint32_t prev = (i == 0u) ? p_top : qk__q(anc, 1u + i);
                    qk_and_compute(kb, cost, prev, qk__q(tgt, k + i), qk__q(anc, 2u + i));
                }
                qk_ccx(kb, cost, qk__q(anc, tail - 1u), qk__q(tgt, k + tail - 2u), qk__q(tgt, k + tail - 1u));
                qk_cx(kb, cost, qk__q(anc, tail - 1u), qk__q(tgt, k + tail - 2u));
                for (uint32_t i = tail - 3u; i >= 1u; --i) {
                    qk_and_uncompute(kb, cost, qk__q(anc, 1u + i), qk__q(tgt, k + i), qk__q(anc, 2u + i), sb);
                    qk_cx(kb, cost, qk__q(anc, 1u + i), qk__q(tgt, k + i));
                }
                qk_and_uncompute(kb, cost, p_top, qk__q(tgt, k), qk__q(anc, 2u), sb);
                qk_cx(kb, cost, p_top, qk__q(tgt, k));
            }
            qk_and_uncompute(kb, cost, ctrl, c_out, p_top, sb);
        }
    }

    for (uint32_t ii = maj_hi; ii > 0u; --ii) {
        uint32_t i = ii - 1u;
        uint32_t cin = (i == 0u) ? c0 : qk__q(lhs, i - 1u);
        qk_ccx(kb, cost, cin, qk__q(tgt, i), qk__q(lhs, i));
        qk_cx(kb, cost, qk__q(lhs, i), cin);
        qk_cx(kb, cost, cin, qk__q(lhs, i));
        qk_ccx(kb, cost, ctrl, qk__q(lhs, i), qk__q(tgt, i));
        qk_cx(kb, cost, cin, qk__q(lhs, i));
        qk_cx(kb, cost, qk__q(lhs, i), qk__q(tgt, i));
    }
    if (next_cbit) (*next_cbit)++;
}

/* Exact Controlled Subtractor: if (ctrl) tgt -= lhs mod 2^n_tgt via bitwise complement */
static inline void qk_ctrl_isub_cuccaro(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    uint32_t ctrl,
    const uint32_t *lhs,
    uint32_t k_in,
    const uint32_t *tgt,
    uint32_t n_tgt,
    const uint32_t *anc,
    uint32_t *next_cbit
) {
    if (k_in == 0u || n_tgt == 0u) return;
    for (uint32_t i = 0u; i < n_tgt; ++i) qk_x(kb, cost, qk__q(tgt, i));
    qk_ctrl_iadd_cuccaro(kb, cost, ctrl, lhs, k_in, tgt, n_tgt, anc, next_cbit);
    for (uint32_t i = 0u; i < n_tgt; ++i) qk_x(kb, cost, qk__q(tgt, i));
}

/* ============================================================================
 * Kernel 3: Gauss-Jordan Min-Fanout In-Place GL_w(F_2) Synthesis (up to 256b)
 * ============================================================================ */
static inline void qk_linf2_inplace(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *q,
    uint32_t w,
    const qcc_bv256_t *rows,
    const qcc_bv256_t *const_vec
) {
    if (w == 0u || w > QIR_MAX_VAL_BITS) return;
    qcc_bv256_t *A = (qcc_bv256_t *)malloc((size_t)w * sizeof(qcc_bv256_t));
    memcpy(A, rows, (size_t)w * sizeof(qcc_bv256_t));

    typedef struct { uint8_t kind; uint16_t a, b; } step_t;
    size_t cap = 4096u, nops = 0u;
    step_t *ops = (step_t *)malloc(cap * sizeof(step_t));

    for (uint32_t col = 0u; col < w; ++col) {
        uint32_t best_piv = col, best_cost = UINT32_MAX;
        for (uint32_t r = col; r < w; ++r) {
            if (!qcc_bv_get(A[r].w, col)) continue;
            uint32_t pc = qcc_bv_popcount(A[r].w) + ((r != col) ? 2u : 0u);
            if (pc < best_cost) { best_cost = pc; best_piv = r; }
        }
        if (!qcc_bv_get(A[best_piv].w, col)) continue;
        if (best_piv != col) {
            qcc_bv256_t tmp = A[col]; A[col] = A[best_piv]; A[best_piv] = tmp;
            if (nops >= cap) { cap *= 2u; ops = (step_t *)realloc(ops, cap * sizeof(step_t)); }
            ops[nops++] = (step_t){1u, (uint16_t)col, (uint16_t)best_piv};
        }
        for (uint32_t r = 0u; r < w; ++r) {
            if (r != col && qcc_bv_get(A[r].w, col)) {
                for (uint32_t ww = 0u; ww < QIR_WORDS; ++ww) A[r].w[ww] ^= A[col].w[ww];
                if (nops >= cap) { cap *= 2u; ops = (step_t *)realloc(ops, cap * sizeof(step_t)); }
                ops[nops++] = (step_t){0u, (uint16_t)r, (uint16_t)col};
            }
        }
    }

    for (size_t k = nops; k > 0u; --k) {
        step_t st = ops[k - 1u];
        uint32_t qa = qk__q(q, st.a), qb = qk__q(q, st.b);
        if (st.kind == 1u) {
            qk_cx(kb, cost, qa, qb);
            qk_cx(kb, cost, qb, qa);
            qk_cx(kb, cost, qa, qb);
        } else {
            qk_cx(kb, cost, qb, qa);
        }
    }
    if (const_vec) {
        for (uint32_t r = 0u; r < w; ++r) {
            if (qcc_bv_get(const_vec->w, r)) qk_x(kb, cost, qk__q(q, r));
        }
    }
    free(A);
    free(ops);
}

/* Out-of-place Affine F_2 Synthesis: tgt[0..out_w-1] ^= M*in ^ c into clean tgt */
static inline void qk_linf2_out_of_place(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *in_q,
    uint32_t in_w,
    const uint32_t *tgt_q,
    uint32_t out_w,
    const qcc_bv256_t *rows,
    const qcc_bv256_t *const_vec
) {
    for (uint32_t r = 0u; r < out_w; ++r) {
        uint32_t qr = qk__q(tgt_q, r);
        for (uint32_t c = 0u; c < in_w; ++c) {
            if (qcc_bv_get(rows[r].w, c)) qk_cx(kb, cost, qk__q(in_q, c), qr);
        }
        if (const_vec && qcc_bv_get(const_vec->w, r)) qk_x(kb, cost, qr);
    }
}

/* ============================================================================
 * Kernel 4: Unified Out-of-Place QROAM Binary-Tree Unary ROM (Column Sliceable)
 * Guarantees EXACT est == actual because cost & codegen execute identical code!
 * ============================================================================ */

static inline bool qk__qrom_range_equal(
    const uint64_t *table,
    uint32_t num_rows,
    uint16_t val_words,
    uint32_t col_lo,
    uint32_t col_hi,
    uint32_t lo_a,
    uint32_t hi_a
) {
    if (lo_a >= num_rows) return true;
    uint32_t last = (hi_a < num_rows) ? hi_a : num_rows;
    if (last <= lo_a + 1u) return true;
    for (uint32_t c = col_lo; c < col_hi; ++c) {
        uint32_t cw = c >> 6u, cb = c & 63u;
        uint64_t b0 = (table[(size_t)lo_a * val_words + cw] >> cb) & 1u;
        for (uint32_t a = lo_a + 1u; a < last; ++a) {
            uint64_t ba = (table[(size_t)a * val_words + cw] >> cb) & 1u;
            if (ba != b0) return false;
        }
    }
    return true;
}

static inline void qk__qrom_rec(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint64_t *table,
    uint32_t num_rows,
    uint16_t val_words,
    uint32_t col_lo,
    uint32_t col_hi,
    const uint32_t *addr_lsb,
    uint32_t addr_bits,
    const uint32_t *tgt,       /* length >= (col_hi - col_lo), indexed by (c - col_lo) */
    const uint32_t *anc,       /* length >= addr_bits */
    uint32_t *next_cbit,
    uint32_t ctrl,
    uint32_t lvl_msb,          /* 0 = MSB .. addr_bits */
    uint32_t lo_a,
    uint32_t hi_a
) {
    if (lo_a >= num_rows) return;
    if (lvl_msb >= addr_bits || qk__qrom_range_equal(table, num_rows, val_words, col_lo, col_hi, lo_a, hi_a)) {
        for (uint32_t c = col_lo; c < col_hi; ++c) {
            uint32_t cw = c >> 6u, cb = c & 63u;
            if ((table[(size_t)lo_a * val_words + cw] >> cb) & 1u) {
                uint32_t qt = qk__q(tgt, c - col_lo);
                if (ctrl == UINT32_MAX) qk_x(kb, cost, qt);
                else                    qk_cx(kb, cost, ctrl, qt);
            }
        }
        return;
    }

    uint32_t mid = lo_a + ((hi_a - lo_a) >> 1u);
    uint32_t qbit = qk__q(addr_lsb, addr_bits - 1u - lvl_msb);

    if (mid >= num_rows) {
        /* Only left child [lo_a, mid) active -> condition on !qbit if right side could exist? No, num_rows < mid means valid addresses < num_rows < mid have qbit==0! */
        qk__qrom_rec(kb, cost, table, num_rows, val_words, col_lo, col_hi, addr_lsb, addr_bits, tgt, anc, next_cbit, ctrl, lvl_msb + 1u, lo_a, mid);
        return;
    }

    if (ctrl == UINT32_MAX) {
        qk_x(kb, cost, qbit);
        qk__qrom_rec(kb, cost, table, num_rows, val_words, col_lo, col_hi, addr_lsb, addr_bits, tgt, anc, next_cbit, qbit, lvl_msb + 1u, lo_a, mid);
        qk_x(kb, cost, qbit);
        qk__qrom_rec(kb, cost, table, num_rows, val_words, col_lo, col_hi, addr_lsb, addr_bits, tgt, anc, next_cbit, qbit, lvl_msb + 1u, mid, hi_a);
        return;
    }

    uint32_t c_child = qk__q(anc, lvl_msb);
    uint32_t sb = next_cbit ? *next_cbit : 0u;
    if (next_cbit) (*next_cbit)++;
    if (cost && lvl_msb + 1u > cost->temp_anc_needed) cost->temp_anc_needed = lvl_msb + 1u;
    if (cost && cost->scratch_bits_needed < 1u) cost->scratch_bits_needed = 1u;

    qk_x(kb, cost, qbit);
    qk_and_compute(kb, cost, ctrl, qbit, c_child);
    qk_x(kb, cost, qbit);

    qk__qrom_rec(kb, cost, table, num_rows, val_words, col_lo, col_hi, addr_lsb, addr_bits, tgt, anc, next_cbit, c_child, lvl_msb + 1u, lo_a, mid);

    /* Switch active branch from (ctrl & !qbit) to (ctrl & qbit) via CNOT(ctrl, c_child)! */
    qk_cx(kb, cost, ctrl, c_child);
    qk__qrom_rec(kb, cost, table, num_rows, val_words, col_lo, col_hi, addr_lsb, addr_bits, tgt, anc, next_cbit, c_child, lvl_msb + 1u, mid, hi_a);

    qk_and_uncompute(kb, cost, ctrl, qbit, c_child, sb);
}

static inline void qk_lut_rom_qrom(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint64_t *table,
    uint32_t num_rows,
    uint16_t val_words,
    uint32_t col_lo,
    uint32_t col_hi,
    const uint32_t *addr_lsb,
    uint32_t addr_bits,
    const uint32_t *tgt,
    const uint32_t *anc,
    uint32_t *next_cbit
) {
    uint32_t span = (addr_bits < 31u) ? (UINT32_C(1) << addr_bits) : num_rows;
    qk__qrom_rec(kb, cost, table, num_rows, val_words, col_lo, col_hi, addr_lsb, addr_bits,
                 tgt, anc, next_cbit, UINT32_MAX, 0u, 0u, span);
}

/* Kernel 4B: Adaptive-Peel / Binary-Tree HMR Uncomputation of any materialized LUT register `tgt` (0-CCX if linear or unary tree with HMR!) */
static inline void qk_onehot_compute(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *addr_msb,
    uint32_t n_bits,
    const uint32_t *leaves
) {
    if (n_bits == 0u) {
        qk_x(kb, cost, qk__q(leaves, 0u));
        return;
    }
    uint32_t step0 = UINT32_C(1) << (n_bits - 1u);
    uint32_t q0 = qk__q(addr_msb, 0u);
    qk_x(kb, cost, q0);
    qk_cx(kb, cost, q0, qk__q(leaves, 0u));
    qk_x(kb, cost, q0);
    qk_cx(kb, cost, q0, qk__q(leaves, step0));
    for (uint32_t lvl = 1u; lvl < n_bits; ++lvl) {
        uint32_t step = UINT32_C(1) << (n_bits - 1u - lvl);
        uint32_t cnt  = UINT32_C(1) << lvl;
        uint32_t bit  = qk__q(addr_msb, lvl);
        for (uint32_t k = 0u; k < cnt; ++k) {
            uint32_t d0 = qk__q(leaves, 2u * k * step);
            uint32_t d1 = qk__q(leaves, (2u * k + 1u) * step);
            qk_and_compute(kb, cost, d0, bit, d1);
            qk_cx(kb, cost, d1, d0);
        }
    }
}

static inline void qk_onehot_uncompute_mbuc(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *addr_msb,
    uint32_t n_bits,
    const uint32_t *leaves,
    uint32_t *next_cbit
) {
    if (n_bits == 0u) {
        qk_x(kb, cost, qk__q(leaves, 0u));
        return;
    }
    uint32_t sb = next_cbit ? *next_cbit : 0u;
    for (uint32_t l_rev = n_bits - 1u; l_rev >= 1u; --l_rev) {
        uint32_t step = UINT32_C(1) << (n_bits - 1u - l_rev);
        uint32_t cnt  = UINT32_C(1) << l_rev;
        uint32_t bit  = qk__q(addr_msb, l_rev);
        for (uint32_t k = 0u; k < cnt; ++k) {
            uint32_t d0 = qk__q(leaves, 2u * k * step);
            uint32_t d1 = qk__q(leaves, (2u * k + 1u) * step);
            qk_cx(kb, cost, d1, d0);
            qk_and_uncompute(kb, cost, d0, bit, d1, sb);
            sb++;
        }
    }
    uint32_t step0 = UINT32_C(1) << (n_bits - 1u);
    uint32_t q0 = qk__q(addr_msb, 0u);
    qk_cx(kb, cost, q0, qk__q(leaves, step0));
    qk_x(kb, cost, q0);
    qk_cx(kb, cost, q0, qk__q(leaves, 0u));
    qk_x(kb, cost, q0);
    if (next_cbit) *next_cbit = sb;
}

/* HMR Parity-LUT uncompute for n_addr <= 10 (returns 1 if used, 0 to fallback to reverse QROM) */
static inline bool qk_hmr_parity_lut_uncompute(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint64_t *table,
    uint32_t num_rows,
    uint16_t val_words,
    uint32_t col_lo,
    uint32_t col_hi,
    const uint32_t *addr_lsb,
    uint32_t n_addr,
    const uint32_t *tgt,
    const uint32_t *extra_anc, /* Note: after HMR, all (col_hi-col_lo) tgt qubits are clean & reusable! */
    uint32_t *next_cbit
) {
    uint32_t n_tgt = col_hi - col_lo;
    if (n_addr < 2u || n_addr > 13u || n_tgt == 0u) return false;
    uint32_t nm = n_addr >> 1u;
    uint32_t nl = n_addr - nm;
    uint32_t m_cnt = UINT32_C(1) << nm;
    uint32_t l_cnt = UINT32_C(1) << nl;
    uint32_t total_leaves = m_cnt + l_cnt;
    uint32_t borrow = (n_tgt < total_leaves) ? n_tgt : total_leaves;
    uint32_t extra  = total_leaves - borrow;

    if (cost) {
        if (extra > cost->temp_anc_needed) cost->temp_anc_needed = extra;
        if (cost->scratch_bits_needed < n_tgt + 1u) cost->scratch_bits_needed = n_tgt + 1u;
    }

    uint32_t *hmr_bits = (uint32_t *)malloc((size_t)n_tgt * sizeof(uint32_t));
    uint32_t base_b = next_cbit ? *next_cbit : 0u;
    for (uint32_t j = 0u; j < n_tgt; ++j) {
        hmr_bits[j] = base_b + j;
        qk_hmr(kb, cost, qk__q(tgt, j), hmr_bits[j]);
    }
    uint32_t pbit = base_b + n_tgt;
    uint32_t cb_cur = pbit + 1u;

    uint32_t leaves_buf[256];
    for (uint32_t i = 0u; i < borrow; ++i) leaves_buf[i] = qk__q(tgt, i);
    for (uint32_t i = 0u; i < extra; ++i)  leaves_buf[borrow + i] = qk__q(extra_anc, i);

    /* addr_msb in big-endian order from addr_lsb */
    uint32_t addr_msb[16];
    for (uint32_t i = 0u; i < n_addr; ++i) addr_msb[i] = qk__q(addr_lsb, n_addr - 1u - i);

    const uint32_t *hi_lv = leaves_buf;
    const uint32_t *lo_lv = leaves_buf + m_cnt;
    qk_onehot_compute(kb, cost, addr_msb, nm, hi_lv);
    qk_onehot_compute(kb, cost, addr_msb + nm, nl, lo_lv);

    for (uint32_t a = 0u; a < num_rows; ++a) {
        bool any = false;
        qk_bit_store0(kb, cost, pbit);
        for (uint32_t j = 0u; j < n_tgt; ++j) {
            uint32_t c = col_lo + j;
            if ((table[(size_t)a * val_words + (c >> 6u)] >> (c & 63u)) & 1u) {
                qk_bit_invert_if(kb, cost, pbit, hmr_bits[j]);
                any = true;
            }
        }
        if (any) qk_cz_if(kb, cost, hi_lv[a >> nl], lo_lv[a & (l_cnt - 1u)], pbit);
    }

    qk_onehot_uncompute_mbuc(kb, cost, addr_msb + nm, nl, lo_lv, &cb_cur);
    qk_onehot_uncompute_mbuc(kb, cost, addr_msb, nm, hi_lv, &cb_cur);

    if (next_cbit) *next_cbit = cb_cur;
    free(hmr_bits);
    return true;
}

/* Kernel 4C: Virtual-Slot-0 Differential Select-Swap QROAM (Emit-or-Count Zero-Discrepancy) */
static inline bool qk__vslot_table_bit(
    const uint64_t *table, uint32_t num_rows, uint16_t val_words, uint32_t r, uint32_t col
) {
    if (r >= num_rows) return false;
    return ((table[(size_t)r * val_words + (col >> 6u)] >> (col & 63u)) & UINT64_C(1)) != 0u;
}

static inline void qk__qrom_vslot_walk(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint64_t *pack_tab,
    uint32_t num_rows,
    uint32_t pack_stride,
    uint32_t w_act,
    uint32_t b_ss,
    uint32_t n_blocks,
    const uint32_t *addr_msb,
    uint32_t hi_bits,
    const uint32_t *all_slots,
    const uint32_t *work_anc,
    uint32_t *next_cbit,
    uint32_t ctrl,
    uint32_t lvl_msb,
    uint32_t lo_blk,
    uint32_t hi_blk
) {
    if (lo_blk >= n_blocks) return;
    if (lvl_msb >= hi_bits) {
        for (uint32_t d = 0u; d < b_ss; ++d) {
            uint32_t r = lo_blk * b_ss + d;
            if (r >= num_rows) break;
            for (uint32_t w = 0u; w < w_act; ++w) {
                if ((pack_tab[(size_t)r * pack_stride + (w >> 6u)] >> (w & 63u)) & UINT64_C(1)) {
                    uint32_t qt = all_slots[(size_t)d * w_act + w];
                    if (ctrl == UINT32_MAX) qk_x(kb, cost, qt);
                    else                    qk_cx(kb, cost, ctrl, qt);
                }
            }
        }
        return;
    }

    uint32_t mid = lo_blk + ((hi_blk - lo_blk) >> 1u);
    uint32_t qbit = addr_msb[lvl_msb];
    if (mid >= n_blocks) {
        qk__qrom_vslot_walk(kb, cost, pack_tab, num_rows, pack_stride, w_act, b_ss, n_blocks,
                            addr_msb, hi_bits, all_slots, work_anc, next_cbit, ctrl, lvl_msb + 1u, lo_blk, mid);
        return;
    }
    if (ctrl == UINT32_MAX) {
        qk_x(kb, cost, qbit);
        qk__qrom_vslot_walk(kb, cost, pack_tab, num_rows, pack_stride, w_act, b_ss, n_blocks,
                            addr_msb, hi_bits, all_slots, work_anc, next_cbit, qbit, lvl_msb + 1u, lo_blk, mid);
        qk_x(kb, cost, qbit);
        qk__qrom_vslot_walk(kb, cost, pack_tab, num_rows, pack_stride, w_act, b_ss, n_blocks,
                            addr_msb, hi_bits, all_slots, work_anc, next_cbit, qbit, lvl_msb + 1u, mid, hi_blk);
        return;
    }

    uint32_t c_child = qk__q(work_anc, lvl_msb - 1u);
    uint32_t sb = next_cbit ? (*next_cbit)++ : 0u;
    if (cost && lvl_msb > cost->temp_anc_needed) cost->temp_anc_needed = lvl_msb;
    if (cost && cost->scratch_bits_needed < 1u) cost->scratch_bits_needed = 1u;

    qk_x(kb, cost, qbit);
    qk_and_compute(kb, cost, ctrl, qbit, c_child);
    qk_x(kb, cost, qbit);
    qk__qrom_vslot_walk(kb, cost, pack_tab, num_rows, pack_stride, w_act, b_ss, n_blocks,
                        addr_msb, hi_bits, all_slots, work_anc, next_cbit, c_child, lvl_msb + 1u, lo_blk, mid);
    qk_cx(kb, cost, ctrl, c_child);
    qk__qrom_vslot_walk(kb, cost, pack_tab, num_rows, pack_stride, w_act, b_ss, n_blocks,
                        addr_msb, hi_bits, all_slots, work_anc, next_cbit, c_child, lvl_msb + 1u, mid, hi_blk);
    qk_and_uncompute(kb, cost, ctrl, qbit, c_child, sb);
}

static inline void qk_qroam_vslot0(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint64_t *table,
    uint32_t num_rows,
    uint16_t val_words,
    uint32_t col_lo,
    uint32_t col_hi,
    const uint32_t *addr_lsb,
    uint32_t addr_bits,
    const uint32_t *tgt,
    uint32_t b_sh,
    const uint32_t *anc_pool,
    uint32_t anc_pool_len,
    uint32_t *next_cbit
) {
    if (num_rows == 0u || col_lo >= col_hi || addr_bits == 0u) return;
    if (b_sh == 0u) {
        qk_lut_rom_qrom(kb, cost, table, num_rows, val_words, col_lo, col_hi, addr_lsb, addr_bits, tgt, anc_pool, next_cbit);
        return;
    }
    if (b_sh > addr_bits) b_sh = addr_bits;
    if (b_sh > 6u) b_sh = 6u;

    uint32_t n_cols = col_hi - col_lo;
    uint32_t *act_cols = (uint32_t *)malloc((size_t)n_cols * sizeof(uint32_t));
    uint32_t *act_tgts = (uint32_t *)malloc((size_t)n_cols * sizeof(uint32_t));
    if (!act_cols || !act_tgts) abort();
    uint32_t w_act = 0u;
    for (uint32_t c = col_lo; c < col_hi; ++c) {
        bool nz = false;
        for (uint32_t r = 0u; r < num_rows; ++r) {
            if (qk__vslot_table_bit(table, num_rows, val_words, r, c)) { nz = true; break; }
        }
        if (nz) {
            act_cols[w_act] = c;
            act_tgts[w_act] = qk__q(tgt, c - col_lo);
            w_act++;
        }
    }
    if (w_act == 0u) {
        free(act_cols);
        free(act_tgts);
        return;
    }

    uint32_t b_ss = UINT32_C(1) << b_sh;
    uint32_t hi_bits = addr_bits - b_sh;
    uint32_t n_blocks = (num_rows + b_ss - 1u) / b_ss;
    uint32_t addr_msb[32];
    for (uint32_t i = 0u; i < hi_bits; ++i) addr_msb[i] = qk__q(addr_lsb, addr_bits - 1u - i);

    uint32_t pack_stride = (w_act + 63u) / 64u;
    uint32_t total_rows_pack = n_blocks * b_ss;
    uint64_t *pack_tab = (uint64_t *)calloc((size_t)total_rows_pack * pack_stride, sizeof(uint64_t));
    uint32_t slot_total = (b_ss - 1u) * w_act;
    uint8_t *slot_ever_nz = (uint8_t *)calloc(slot_total, 1u);
    if (!pack_tab || !slot_ever_nz) abort();

    for (uint32_t q = 0u; q < n_blocks; ++q) {
        uint32_t r0 = q * b_ss;
        for (uint32_t w = 0u; w < w_act; ++w) {
            bool b0 = qk__vslot_table_bit(table, num_rows, val_words, r0, act_cols[w]);
            if (b0) pack_tab[(size_t)r0 * pack_stride + (w >> 6u)] |= (UINT64_C(1) << (w & 63u));
            for (uint32_t d = 1u; d < b_ss; ++d) {
                uint32_t rd = r0 + d;
                bool bd = qk__vslot_table_bit(table, num_rows, val_words, rd, act_cols[w]);
                if (bd ^ b0) {
                    pack_tab[(size_t)rd * pack_stride + (w >> 6u)] |= (UINT64_C(1) << (w & 63u));
                    slot_ever_nz[(d - 1u) * w_act + w] = 1u;
                }
            }
        }
    }

    uint32_t *slot_qs = (uint32_t *)malloc((size_t)slot_total * sizeof(uint32_t));
    uint32_t *all_slots = (uint32_t *)malloc((size_t)b_ss * w_act * sizeof(uint32_t));
    if (!slot_qs || !all_slots) abort();
    uint32_t nz_slots = 0u;
    for (uint32_t idx = 0u; idx < slot_total; ++idx) {
        slot_qs[idx] = slot_ever_nz[idx] ? qk__q(anc_pool, nz_slots++) : UINT32_MAX;
    }
    const uint32_t *work_anc = anc_pool ? (anc_pool + nz_slots) : NULL;
    for (uint32_t w = 0u; w < w_act; ++w) all_slots[w] = act_tgts[w];
    for (uint32_t idx = 0u; idx < slot_total; ++idx) all_slots[w_act + idx] = slot_qs[idx];

    /* Step A: Unary walk over hi_bits loading T[q*B+0] into tgt and D[q,d] into slot_qs */
    qk__qrom_vslot_walk(kb, cost, pack_tab, total_rows_pack, pack_stride, w_act, b_ss, n_blocks,
                        addr_msb, hi_bits, all_slots, work_anc, next_cbit,
                        UINT32_MAX, 0u, 0u, (hi_bits < 31u) ? (UINT32_C(1) << hi_bits) : n_blocks);

    /* Step B: Non-Zero One-Hot Gather over addr_lsb[0 .. b_sh-1] */
    if (b_sh == 1u) {
        uint32_t l0 = qk__q(addr_lsb, 0u);
        for (uint32_t w = 0u; w < w_act; ++w) {
            if (slot_ever_nz[w]) qk_ccx(kb, cost, l0, slot_qs[w], act_tgts[w]);
        }
    } else {
        uint32_t r_msb[16];
        for (uint32_t i = 0u; i < b_sh; ++i) r_msb[i] = qk__q(addr_lsb, b_sh - 1u - i);
        uint32_t step0 = UINT32_C(1) << (b_sh - 2u);
        qk_and_compute(kb, cost, r_msb[0], r_msb[1], qk__q(work_anc, 3u * step0 - 1u));
        qk_cx(kb, cost, r_msb[0], qk__q(work_anc, 2u * step0 - 1u));
        qk_cx(kb, cost, qk__q(work_anc, 3u * step0 - 1u), qk__q(work_anc, 2u * step0 - 1u));
        qk_cx(kb, cost, r_msb[1], qk__q(work_anc, step0 - 1u));
        qk_cx(kb, cost, qk__q(work_anc, 3u * step0 - 1u), qk__q(work_anc, step0 - 1u));

        for (uint32_t lvl = 2u; lvl < b_sh; ++lvl) {
            uint32_t step = UINT32_C(1) << (b_sh - 1u - lvl);
            uint32_t m_prev = UINT32_C(1) << lvl;
            uint32_t bit = r_msb[lvl];
            for (uint32_t k = 1u; k < m_prev; ++k) {
                uint32_t d0 = 2u * k * step, d1 = (2u * k + 1u) * step;
                qk_and_compute(kb, cost, qk__q(work_anc, d0 - 1u), bit, qk__q(work_anc, d1 - 1u));
                qk_cx(kb, cost, qk__q(work_anc, d1 - 1u), qk__q(work_anc, d0 - 1u));
            }
            qk_cx(kb, cost, bit, qk__q(work_anc, step - 1u));
            for (uint32_t k = 1u; k < m_prev; ++k) {
                qk_cx(kb, cost, qk__q(work_anc, (2u * k + 1u) * step - 1u), qk__q(work_anc, step - 1u));
            }
        }

        for (uint32_t d = 1u; d < b_ss; ++d) {
            uint32_t hq = qk__q(work_anc, d - 1u);
            for (uint32_t w = 0u; w < w_act; ++w) {
                uint32_t idx = (d - 1u) * w_act + w;
                if (slot_ever_nz[idx]) qk_ccx(kb, cost, hq, slot_qs[idx], act_tgts[w]);
            }
        }

        for (uint32_t l_rev = b_sh - 1u; l_rev >= 2u; --l_rev) {
            uint32_t step = UINT32_C(1) << (b_sh - 1u - l_rev);
            uint32_t m_prev = UINT32_C(1) << l_rev;
            uint32_t bit = r_msb[l_rev];
            for (uint32_t k = 1u; k < m_prev; ++k) {
                qk_cx(kb, cost, qk__q(work_anc, (2u * k + 1u) * step - 1u), qk__q(work_anc, step - 1u));
            }
            qk_cx(kb, cost, bit, qk__q(work_anc, step - 1u));
            for (uint32_t k = 1u; k < m_prev; ++k) {
                uint32_t d0 = 2u * k * step, d1 = (2u * k + 1u) * step;
                uint32_t sb = next_cbit ? (*next_cbit)++ : 0u;
                qk_cx(kb, cost, qk__q(work_anc, d1 - 1u), qk__q(work_anc, d0 - 1u));
                qk_and_uncompute(kb, cost, qk__q(work_anc, d0 - 1u), bit, qk__q(work_anc, d1 - 1u), sb);
            }
        }
        uint32_t sb0 = next_cbit ? (*next_cbit)++ : 0u;
        qk_cx(kb, cost, qk__q(work_anc, 3u * step0 - 1u), qk__q(work_anc, step0 - 1u));
        qk_cx(kb, cost, r_msb[1], qk__q(work_anc, step0 - 1u));
        qk_cx(kb, cost, qk__q(work_anc, 3u * step0 - 1u), qk__q(work_anc, 2u * step0 - 1u));
        qk_cx(kb, cost, r_msb[0], qk__q(work_anc, 2u * step0 - 1u));
        qk_and_uncompute(kb, cost, r_msb[0], r_msb[1], qk__q(work_anc, 3u * step0 - 1u), sb0);
    }

    /* Step C: HMR all active slot qubits + 1-Bit Split-Unary Phase Fixup */
    uint32_t *m_slots = (uint32_t *)malloc((size_t)slot_total * sizeof(uint32_t));
    if (!m_slots) abort();
    uint32_t cb_cur = next_cbit ? *next_cbit : 0u;
    for (uint32_t idx = 0u; idx < slot_total; ++idx) {
        if (slot_ever_nz[idx]) {
            m_slots[idx] = cb_cur++;
            qk_hmr(kb, cost, slot_qs[idx], m_slots[idx]);
        } else {
            m_slots[idx] = UINT32_MAX;
        }
    }
    uint32_t pbit = cb_cur++;
    if (next_cbit) *next_cbit = cb_cur;

    if (hi_bits == 0u) {
        bool any = false;
        qk_bit_store0(kb, cost, pbit);
        for (uint32_t d = 1u; d < b_ss; ++d) {
            for (uint32_t w = 0u; w < w_act; ++w) {
                if ((pack_tab[(size_t)d * pack_stride + (w >> 6u)] >> (w & 63u)) & UINT64_C(1)) {
                    qk_bit_invert_if(kb, cost, pbit, m_slots[(d - 1u) * w_act + w]);
                    any = true;
                }
            }
        }
        if (any) qk_neg_if(kb, cost, pbit);
    } else if (hi_bits == 1u) {
        for (uint32_t q = 0u; q < n_blocks; ++q) {
            bool any = false;
            qk_bit_store0(kb, cost, pbit);
            for (uint32_t d = 1u; d < b_ss; ++d) {
                uint32_t rd = q * b_ss + d;
                for (uint32_t w = 0u; w < w_act; ++w) {
                    if ((pack_tab[(size_t)rd * pack_stride + (w >> 6u)] >> (w & 63u)) & UINT64_C(1)) {
                        qk_bit_invert_if(kb, cost, pbit, m_slots[(d - 1u) * w_act + w]);
                        any = true;
                    }
                }
            }
            if (any) {
                if (q == 0u) qk_x(kb, cost, addr_msb[0]);
                qk_z_if(kb, cost, addr_msb[0], pbit);
                if (q == 0u) qk_x(kb, cost, addr_msb[0]);
            }
        }
    } else {
        uint32_t nh0 = hi_bits / 2u, nl0 = hi_bits - nh0;
        uint32_t need0 = (UINT32_C(1) << nh0) + (UINT32_C(1) << nl0);
        if (anc_pool_len >= need0) {
            const uint32_t *hi_lv = anc_pool;
            const uint32_t *lo_lv = anc_pool ? (anc_pool + (UINT32_C(1) << nh0)) : NULL;
            qk_onehot_compute(kb, cost, addr_msb, nh0, hi_lv);
            qk_onehot_compute(kb, cost, addr_msb + nh0, nl0, lo_lv);
            uint32_t lo_mask = (UINT32_C(1) << nl0) - 1u;
            for (uint32_t q = 0u; q < n_blocks; ++q) {
                bool any = false;
                qk_bit_store0(kb, cost, pbit);
                for (uint32_t d = 1u; d < b_ss; ++d) {
                    uint32_t rd = q * b_ss + d;
                    for (uint32_t w = 0u; w < w_act; ++w) {
                        if ((pack_tab[(size_t)rd * pack_stride + (w >> 6u)] >> (w & 63u)) & UINT64_C(1)) {
                            qk_bit_invert_if(kb, cost, pbit, m_slots[(d - 1u) * w_act + w]);
                            any = true;
                        }
                    }
                }
                if (any) qk_cz_if(kb, cost, qk__q(hi_lv, q >> nl0), qk__q(lo_lv, q & lo_mask), pbit);
            }
            qk_onehot_uncompute_mbuc(kb, cost, addr_msb + nh0, nl0, lo_lv, next_cbit);
            qk_onehot_uncompute_mbuc(kb, cost, addr_msb, nh0, hi_lv, next_cbit);
        } else {
            uint32_t rem1 = hi_bits - 1u;
            uint32_t nm1 = rem1 / 2u, nl1 = rem1 - nm1;
            uint32_t need1 = 1u + (UINT32_C(1) << nm1) + (UINT32_C(1) << nl1);
            if (anc_pool_len >= need1) {
                uint32_t anc_m = qk__q(anc_pool, 0u);
                const uint32_t *mid_lv = anc_pool ? (anc_pool + 1u) : NULL;
                const uint32_t *lo_lv  = mid_lv ? (mid_lv + (UINT32_C(1) << nm1)) : NULL;
                qk_onehot_compute(kb, cost, addr_msb + 1u, nm1, mid_lv);
                qk_onehot_compute(kb, cost, addr_msb + 1u + nm1, nl1, lo_lv);
                uint32_t m_cnt = UINT32_C(1) << nm1, l_cnt = UINT32_C(1) << nl1;
                for (uint32_t m = 0u; m < m_cnt; ++m) {
                    uint32_t q0_base = (m << nl1);
                    uint32_t q1_base = (UINT32_C(1) << rem1) | (m << nl1);
                    if (q0_base >= n_blocks) break;
                    qk_x(kb, cost, addr_msb[0]);
                    qk_and_compute(kb, cost, addr_msb[0], qk__q(mid_lv, m), anc_m);
                    qk_x(kb, cost, addr_msb[0]);
                    for (uint32_t l = 0u; l < l_cnt && (q0_base + l) < n_blocks; ++l) {
                        uint32_t q = q0_base + l;
                        bool any = false;
                        qk_bit_store0(kb, cost, pbit);
                        for (uint32_t d = 1u; d < b_ss; ++d) {
                            uint32_t rd = q * b_ss + d;
                            for (uint32_t w = 0u; w < w_act; ++w) {
                                if ((pack_tab[(size_t)rd * pack_stride + (w >> 6u)] >> (w & 63u)) & UINT64_C(1)) {
                                    qk_bit_invert_if(kb, cost, pbit, m_slots[(d - 1u) * w_act + w]);
                                    any = true;
                                }
                            }
                        }
                        if (any) qk_cz_if(kb, cost, anc_m, qk__q(lo_lv, l), pbit);
                    }
                    uint32_t sb_m = next_cbit ? (*next_cbit)++ : 0u;
                    if (q1_base < n_blocks) {
                        qk_cx(kb, cost, qk__q(mid_lv, m), anc_m);
                        for (uint32_t l = 0u; l < l_cnt && (q1_base + l) < n_blocks; ++l) {
                            uint32_t q = q1_base + l;
                            bool any = false;
                            qk_bit_store0(kb, cost, pbit);
                            for (uint32_t d = 1u; d < b_ss; ++d) {
                                uint32_t rd = q * b_ss + d;
                                for (uint32_t w = 0u; w < w_act; ++w) {
                                    if ((pack_tab[(size_t)rd * pack_stride + (w >> 6u)] >> (w & 63u)) & UINT64_C(1)) {
                                        qk_bit_invert_if(kb, cost, pbit, m_slots[(d - 1u) * w_act + w]);
                                        any = true;
                                    }
                                }
                            }
                            if (any) qk_cz_if(kb, cost, anc_m, qk__q(lo_lv, l), pbit);
                        }
                        qk_and_uncompute(kb, cost, addr_msb[0], qk__q(mid_lv, m), anc_m, sb_m);
                    } else {
                        qk_x(kb, cost, addr_msb[0]);
                        qk_and_uncompute(kb, cost, addr_msb[0], qk__q(mid_lv, m), anc_m, sb_m);
                        qk_x(kb, cost, addr_msb[0]);
                    }
                }
                qk_onehot_uncompute_mbuc(kb, cost, addr_msb + 1u + nm1, nl1, lo_lv, next_cbit);
                qk_onehot_uncompute_mbuc(kb, cost, addr_msb + 1u, nm1, mid_lv, next_cbit);
            } else {
                uint32_t rem2 = hi_bits - 2u;
                uint32_t nm2 = rem2 / 2u, nl2 = rem2 - nm2;
                uint32_t p01 = qk__q(anc_pool, 0u), anc_m = qk__q(anc_pool, 1u);
                const uint32_t *mid_lv = anc_pool ? (anc_pool + 2u) : NULL;
                const uint32_t *lo_lv  = mid_lv ? (mid_lv + (UINT32_C(1) << nm2)) : NULL;
                qk_onehot_compute(kb, cost, addr_msb + 2u, nm2, mid_lv);
                qk_onehot_compute(kb, cost, addr_msb + 2u + nm2, nl2, lo_lv);
                uint32_t m_cnt = UINT32_C(1) << nm2, l_cnt = UINT32_C(1) << nl2;
                for (uint32_t tv = 0u; tv < 4u; ++tv) {
                    uint32_t q_top = tv << rem2;
                    if (q_top >= n_blocks) break;
                    if ((tv & 2u) == 0u) qk_x(kb, cost, addr_msb[0]);
                    if ((tv & 1u) == 0u) qk_x(kb, cost, addr_msb[1]);
                    qk_and_compute(kb, cost, addr_msb[0], addr_msb[1], p01);
                    for (uint32_t m = 0u; m < m_cnt; ++m) {
                        uint32_t q_base = q_top | (m << nl2);
                        if (q_base >= n_blocks) break;
                        qk_and_compute(kb, cost, p01, qk__q(mid_lv, m), anc_m);
                        for (uint32_t l = 0u; l < l_cnt && (q_base + l) < n_blocks; ++l) {
                            uint32_t q = q_base + l;
                            bool any = false;
                            qk_bit_store0(kb, cost, pbit);
                            for (uint32_t d = 1u; d < b_ss; ++d) {
                                uint32_t rd = q * b_ss + d;
                                for (uint32_t w = 0u; w < w_act; ++w) {
                                    if ((pack_tab[(size_t)rd * pack_stride + (w >> 6u)] >> (w & 63u)) & UINT64_C(1)) {
                                        qk_bit_invert_if(kb, cost, pbit, m_slots[(d - 1u) * w_act + w]);
                                        any = true;
                                    }
                                }
                            }
                            if (any) qk_cz_if(kb, cost, anc_m, qk__q(lo_lv, l), pbit);
                        }
                        uint32_t sbm = next_cbit ? (*next_cbit)++ : 0u;
                        qk_and_uncompute(kb, cost, p01, qk__q(mid_lv, m), anc_m, sbm);
                    }
                    uint32_t sbp = next_cbit ? (*next_cbit)++ : 0u;
                    qk_and_uncompute(kb, cost, addr_msb[0], addr_msb[1], p01, sbp);
                    if ((tv & 1u) == 0u) qk_x(kb, cost, addr_msb[1]);
                    if ((tv & 2u) == 0u) qk_x(kb, cost, addr_msb[0]);
                }
                qk_onehot_uncompute_mbuc(kb, cost, addr_msb + 2u + nm2, nl2, lo_lv, next_cbit);
                qk_onehot_uncompute_mbuc(kb, cost, addr_msb + 2u, nm2, mid_lv, next_cbit);
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

static inline void qk_qroam_vslot0_resources(
    const uint64_t *table,
    uint32_t num_rows,
    uint16_t val_words,
    uint32_t col_lo,
    uint32_t col_hi,
    uint32_t addr_bits,
    uint32_t b_sh,
    uint32_t *out_anc,
    uint32_t *out_cbits
) {
    if (out_anc) *out_anc = 0u;
    if (out_cbits) *out_cbits = 0u;
    if (!table || num_rows == 0u || col_lo >= col_hi || addr_bits == 0u) return;
    if (b_sh == 0u) {
        if (out_anc) *out_anc = addr_bits;
        if (out_cbits) *out_cbits = (addr_bits < 16u) ? (UINT32_C(1) << addr_bits) : 4096u;
        return;
    }
    if (b_sh > addr_bits) b_sh = addr_bits;
    if (b_sh > 6u) b_sh = 6u;

    uint32_t n_cols = col_hi - col_lo;
    uint32_t *act_cols = (uint32_t *)malloc((size_t)n_cols * sizeof(uint32_t));
    if (!act_cols) abort();
    uint32_t w_act = 0u;
    for (uint32_t c = col_lo; c < col_hi; ++c) {
        for (uint32_t r = 0u; r < num_rows; ++r) {
            if (qk__vslot_table_bit(table, num_rows, val_words, r, c)) {
                act_cols[w_act++] = c;
                break;
            }
        }
    }
    if (w_act == 0u) { free(act_cols); return; }

    uint32_t b_ss = UINT32_C(1) << b_sh;
    uint32_t hi_bits = addr_bits - b_sh;
    uint32_t n_blocks = (num_rows + b_ss - 1u) / b_ss;
    uint32_t nz_slots = 0u;
    for (uint32_t d = 1u; d < b_ss; ++d) {
        for (uint32_t w = 0u; w < w_act; ++w) {
            uint32_t c = act_cols[w];
            for (uint32_t q = 0u; q < n_blocks; ++q) {
                uint32_t r0 = q * b_ss, rd = r0 + d;
                bool b0 = qk__vslot_table_bit(table, num_rows, val_words, r0, c);
                bool bd = qk__vslot_table_bit(table, num_rows, val_words, rd, c);
                if (b0 ^ bd) { nz_slots++; break; }
            }
        }
    }
    free(act_cols);

    uint32_t spine_q   = (hi_bits <= 1u) ? 0u : (hi_bits - 1u);
    uint32_t gather_q  = (b_sh >= 2u) ? (b_ss - 1u) : 0u;
    uint32_t peak_q    = nz_slots + ((spine_q > gather_q) ? spine_q : gather_q);
    if (hi_bits == 2u && peak_q < 4u) peak_q = 4u;
    else if (hi_bits == 3u && peak_q < 5u) peak_q = 5u;
    else if (hi_bits >= 4u) {
        uint32_t rem2 = hi_bits - 2u;
        uint32_t nm2 = rem2 / 2u, nl2 = rem2 - nm2;
        uint32_t need2 = 2u + (UINT32_C(1) << nm2) + (UINT32_C(1) << nl2);
        if (peak_q < need2) peak_q = need2;
    }

    qir_kernel_cost_t kc;
    memset(&kc, 0, sizeof(kc));
    uint32_t cb_cnt = 0u;
    qk_qroam_vslot0(NULL, &kc, table, num_rows, val_words, col_lo, col_hi,
                    NULL, addr_bits, NULL, b_sh, NULL, peak_q, &cb_cnt);
    if (out_anc) *out_anc = peak_q;
    if (out_cbits) *out_cbits = cb_cnt;
}

/* ============================================================================
 * Kernel 5: Exact Bilinear GF(2^12) Rank-42 & GF(64) Rank-15 & Tower Inversion
 * ============================================================================ */
static inline void qk_gf12_mul_r42_compute(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *qa,
    const uint32_t *qb,
    const uint32_t *qout,
    uint32_t anc_prod,
    uint32_t *next_cbit,
    const qir_gf12_rank42_t *plan
) {
    if (cost && cost->temp_anc_needed < 1u) cost->temp_anc_needed = 1u;
    if (cost && cost->scratch_bits_needed < 1u) cost->scratch_bits_needed = 1u;
    uint32_t sb = next_cbit ? *next_cbit : 0u;
    for (uint32_t step = 0u; step < QIR_GF12_RANK; ++step) {
        uint32_t p = plan->order[step];
        uint16_t lm = plan->l_mask[p], rm = plan->r_mask[p], wm = plan->w_mask[p];
        if (!lm || !rm || !wm) continue;
        uint32_t piv_a = (uint32_t)__builtin_ctz((unsigned)lm);
        uint32_t piv_b = (uint32_t)__builtin_ctz((unsigned)rm);
        for (uint32_t i = piv_a + 1u; i < 12u; ++i) {
            if ((lm >> i) & 1u) qk_cx(kb, cost, qk__q(qa, i), qk__q(qa, piv_a));
        }
        for (uint32_t j = piv_b + 1u; j < 12u; ++j) {
            if ((rm >> j) & 1u) qk_cx(kb, cost, qk__q(qb, j), qk__q(qb, piv_b));
        }
        qk_and_compute(kb, cost, qk__q(qa, piv_a), qk__q(qb, piv_b), anc_prod);
        for (uint32_t k = 0u; k < 12u; ++k) {
            if ((wm >> k) & 1u) qk_cx(kb, cost, anc_prod, qk__q(qout, k));
        }
        qk_and_uncompute(kb, cost, qk__q(qa, piv_a), qk__q(qb, piv_b), anc_prod, sb);
        for (uint32_t j = 11u; j > piv_b; --j) {
            if ((rm >> j) & 1u) qk_cx(kb, cost, qk__q(qb, j), qk__q(qb, piv_b));
        }
        for (uint32_t i = 11u; i > piv_a; --i) {
            if ((lm >> i) & 1u) qk_cx(kb, cost, qk__q(qa, i), qk__q(qa, piv_a));
        }
    }
    if (next_cbit) (*next_cbit)++;
}

static inline void qk_gf12_mul_r42_uncompute_hmr(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *qa,
    const uint32_t *qb,
    const uint32_t *qout,
    uint32_t *next_cbit,
    const qir_gf12_rank42_t *plan
) {
    if (cost && cost->scratch_bits_needed < 13u) cost->scratch_bits_needed = 13u;
    uint32_t base = next_cbit ? *next_cbit : 0u;
    uint32_t mbits[12];
    for (uint32_t k = 0u; k < 12u; ++k) {
        mbits[k] = base + k;
        qk_hmr(kb, cost, qk__q(qout, k), mbits[k]);
    }
    uint32_t pbit = base + 12u;
    for (uint32_t step = 0u; step < QIR_GF12_RANK; ++step) {
        uint32_t p = plan->order[step];
        uint16_t lm = plan->l_mask[p], rm = plan->r_mask[p], wm = plan->w_mask[p];
        if (!lm || !rm || !wm) continue;
        qk_bit_store0(kb, cost, pbit);
        for (uint32_t k = 0u; k < 12u; ++k) {
            if ((wm >> k) & 1u) qk_bit_invert_if(kb, cost, pbit, mbits[k]);
        }
        uint32_t piv_a = (uint32_t)__builtin_ctz((unsigned)lm);
        uint32_t piv_b = (uint32_t)__builtin_ctz((unsigned)rm);
        for (uint32_t i = piv_a + 1u; i < 12u; ++i) {
            if ((lm >> i) & 1u) qk_cx(kb, cost, qk__q(qa, i), qk__q(qa, piv_a));
        }
        for (uint32_t j = piv_b + 1u; j < 12u; ++j) {
            if ((rm >> j) & 1u) qk_cx(kb, cost, qk__q(qb, j), qk__q(qb, piv_b));
        }
        qk_cz_if(kb, cost, qk__q(qa, piv_a), qk__q(qb, piv_b), pbit);
        for (uint32_t j = 11u; j > piv_b; --j) {
            if ((rm >> j) & 1u) qk_cx(kb, cost, qk__q(qb, j), qk__q(qb, piv_b));
        }
        for (uint32_t i = 11u; i > piv_a; --i) {
            if ((lm >> i) & 1u) qk_cx(kb, cost, qk__q(qa, i), qk__q(qa, piv_a));
        }
    }
    if (next_cbit) *next_cbit = base + 13u;
}

static inline void qk_gf64_mul_r15_compute(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *qa,
    const uint32_t *qb,
    const uint32_t *qout,
    uint32_t anc_prod,
    uint32_t *next_cbit,
    const qir_gf64_rank15_t *plan
) {
    if (cost && cost->temp_anc_needed < 1u) cost->temp_anc_needed = 1u;
    if (cost && cost->scratch_bits_needed < 1u) cost->scratch_bits_needed = 1u;
    uint32_t sb = next_cbit ? *next_cbit : 0u;
    for (uint32_t p = 0u; p < QIR_GF64_RANK; ++p) {
        uint8_t lm = plan->l_mask[p], rm = plan->r_mask[p], wm = plan->w_mask[p];
        if (!lm || !rm || !wm) continue;
        uint32_t piv_a = (uint32_t)__builtin_ctz((unsigned)lm);
        uint32_t piv_b = (uint32_t)__builtin_ctz((unsigned)rm);
        for (uint32_t i = piv_a + 1u; i < 6u; ++i) {
            if ((lm >> i) & 1u) qk_cx(kb, cost, qk__q(qa, i), qk__q(qa, piv_a));
        }
        for (uint32_t j = piv_b + 1u; j < 6u; ++j) {
            if ((rm >> j) & 1u) qk_cx(kb, cost, qk__q(qb, j), qk__q(qb, piv_b));
        }
        qk_and_compute(kb, cost, qk__q(qa, piv_a), qk__q(qb, piv_b), anc_prod);
        for (uint32_t k = 0u; k < 6u; ++k) {
            if ((wm >> k) & 1u) qk_cx(kb, cost, anc_prod, qk__q(qout, k));
        }
        qk_and_uncompute(kb, cost, qk__q(qa, piv_a), qk__q(qb, piv_b), anc_prod, sb);
        for (uint32_t j = 5u; j > piv_b; --j) {
            if ((rm >> j) & 1u) qk_cx(kb, cost, qk__q(qb, j), qk__q(qb, piv_b));
        }
        for (uint32_t i = 5u; i > piv_a; --i) {
            if ((lm >> i) & 1u) qk_cx(kb, cost, qk__q(qa, i), qk__q(qa, piv_a));
        }
    }
    if (next_cbit) (*next_cbit)++;
}

static inline void qk_gf64_mul_r15_uncompute_hmr(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *qa,
    const uint32_t *qb,
    const uint32_t *qout,
    uint32_t *next_cbit,
    const qir_gf64_rank15_t *plan
) {
    if (cost && cost->scratch_bits_needed < 7u) cost->scratch_bits_needed = 7u;
    uint32_t base = next_cbit ? *next_cbit : 0u;
    uint32_t mbits[6];
    for (uint32_t k = 0u; k < 6u; ++k) {
        mbits[k] = base + k;
        qk_hmr(kb, cost, qk__q(qout, k), mbits[k]);
    }
    uint32_t pbit = base + 6u;
    for (uint32_t p = 0u; p < QIR_GF64_RANK; ++p) {
        uint8_t lm = plan->l_mask[p], rm = plan->r_mask[p], wm = plan->w_mask[p];
        if (!lm || !rm || !wm) continue;
        qk_bit_store0(kb, cost, pbit);
        for (uint32_t k = 0u; k < 6u; ++k) {
            if ((wm >> k) & 1u) qk_bit_invert_if(kb, cost, pbit, mbits[k]);
        }
        uint32_t piv_a = (uint32_t)__builtin_ctz((unsigned)lm);
        uint32_t piv_b = (uint32_t)__builtin_ctz((unsigned)rm);
        for (uint32_t i = piv_a + 1u; i < 6u; ++i) {
            if ((lm >> i) & 1u) qk_cx(kb, cost, qk__q(qa, i), qk__q(qa, piv_a));
        }
        for (uint32_t j = piv_b + 1u; j < 6u; ++j) {
            if ((rm >> j) & 1u) qk_cx(kb, cost, qk__q(qb, j), qk__q(qb, piv_b));
        }
        qk_cz_if(kb, cost, qk__q(qa, piv_a), qk__q(qb, piv_b), pbit);
        for (uint32_t j = 5u; j > piv_b; --j) {
            if ((rm >> j) & 1u) qk_cx(kb, cost, qk__q(qb, j), qk__q(qb, piv_b));
        }
        for (uint32_t i = 5u; i > piv_a; --i) {
            if ((lm >> i) & 1u) qk_cx(kb, cost, qk__q(qa, i), qk__q(qa, piv_a));
        }
    }
    if (next_cbit) *next_cbit = base + 7u;
}

/* GF(8) = GF(2)[beta]/(beta^3 + beta + 1) Rank-6 Bilinear Multiplier & 3-CCX Inverse */
static inline void qk_gf8_mul_r6_compute(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *qa,
    const uint32_t *qb,
    const uint32_t *qout,
    uint32_t anc_prod,
    uint32_t *next_cbit
) {
    if (cost && cost->temp_anc_needed < 1u) cost->temp_anc_needed = 1u;
    if (cost && cost->scratch_bits_needed < 1u) cost->scratch_bits_needed = 1u;
    static const uint8_t LM[6] = {1u, 2u, 4u, 3u, 5u, 6u};
    static const uint8_t RM[6] = {1u, 2u, 4u, 3u, 5u, 6u};
    static const uint8_t WM[6] = {7u, 5u, 1u, 2u, 4u, 3u};
    uint32_t sb = next_cbit ? *next_cbit : 0u;
    for (uint32_t p = 0u; p < 6u; ++p) {
        uint8_t lm = LM[p], rm = RM[p], wm = WM[p];
        uint32_t piv_a = (uint32_t)__builtin_ctz((unsigned)lm);
        uint32_t piv_b = (uint32_t)__builtin_ctz((unsigned)rm);
        for (uint32_t i = piv_a + 1u; i < 3u; ++i) {
            if ((lm >> i) & 1u) qk_cx(kb, cost, qk__q(qa, i), qk__q(qa, piv_a));
        }
        for (uint32_t j = piv_b + 1u; j < 3u; ++j) {
            if ((rm >> j) & 1u) qk_cx(kb, cost, qk__q(qb, j), qk__q(qb, piv_b));
        }
        qk_and_compute(kb, cost, qk__q(qa, piv_a), qk__q(qb, piv_b), anc_prod);
        for (uint32_t k = 0u; k < 3u; ++k) {
            if ((wm >> k) & 1u) qk_cx(kb, cost, anc_prod, qk__q(qout, k));
        }
        qk_and_uncompute(kb, cost, qk__q(qa, piv_a), qk__q(qb, piv_b), anc_prod, sb);
        for (uint32_t j = 2u; j > piv_b; --j) {
            if ((rm >> j) & 1u) qk_cx(kb, cost, qk__q(qb, j), qk__q(qb, piv_b));
        }
        for (uint32_t i = 2u; i > piv_a; --i) {
            if ((lm >> i) & 1u) qk_cx(kb, cost, qk__q(qa, i), qk__q(qa, piv_a));
        }
    }
    if (next_cbit) (*next_cbit)++;
}

static inline void qk_gf8_mul_r6_uncompute_hmr(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *qa,
    const uint32_t *qb,
    const uint32_t *qout,
    uint32_t *next_cbit
) {
    if (cost && cost->scratch_bits_needed < 4u) cost->scratch_bits_needed = 4u;
    static const uint8_t LM[6] = {1u, 2u, 4u, 3u, 5u, 6u};
    static const uint8_t RM[6] = {1u, 2u, 4u, 3u, 5u, 6u};
    static const uint8_t WM[6] = {7u, 5u, 1u, 2u, 4u, 3u};
    uint32_t base = next_cbit ? *next_cbit : 0u;
    uint32_t mbits[3];
    for (uint32_t k = 0u; k < 3u; ++k) {
        mbits[k] = base + k;
        qk_hmr(kb, cost, qk__q(qout, k), mbits[k]);
    }
    uint32_t pbit = base + 3u;
    for (uint32_t p = 0u; p < 6u; ++p) {
        uint8_t lm = LM[p], rm = RM[p], wm = WM[p];
        qk_bit_store0(kb, cost, pbit);
        for (uint32_t k = 0u; k < 3u; ++k) {
            if ((wm >> k) & 1u) qk_bit_invert_if(kb, cost, pbit, mbits[k]);
        }
        uint32_t piv_a = (uint32_t)__builtin_ctz((unsigned)lm);
        uint32_t piv_b = (uint32_t)__builtin_ctz((unsigned)rm);
        for (uint32_t i = piv_a + 1u; i < 3u; ++i) {
            if ((lm >> i) & 1u) qk_cx(kb, cost, qk__q(qa, i), qk__q(qa, piv_a));
        }
        for (uint32_t j = piv_b + 1u; j < 3u; ++j) {
            if ((rm >> j) & 1u) qk_cx(kb, cost, qk__q(qb, j), qk__q(qb, piv_b));
        }
        qk_cz_if(kb, cost, qk__q(qa, piv_a), qk__q(qb, piv_b), pbit);
        for (uint32_t j = 2u; j > piv_b; --j) {
            if ((rm >> j) & 1u) qk_cx(kb, cost, qk__q(qb, j), qk__q(qb, piv_b));
        }
        for (uint32_t i = 2u; i > piv_a; --i) {
            if ((lm >> i) & 1u) qk_cx(kb, cost, qk__q(qa, i), qk__q(qa, piv_a));
        }
    }
    if (next_cbit) *next_cbit = base + 4u;
}

static inline void qk_gf8_sq_xor(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *qa,
    const uint32_t *qdst
) {
    qk_cx(kb, cost, qk__q(qa, 0u), qk__q(qdst, 0u));
    qk_cx(kb, cost, qk__q(qa, 1u), qk__q(qdst, 2u));
    qk_cx(kb, cost, qk__q(qa, 2u), qk__q(qdst, 1u));
    qk_cx(kb, cost, qk__q(qa, 2u), qk__q(qdst, 2u));
}

static inline void qk_gf8_inv_alg_compute(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *qa,
    const uint32_t *qout,
    uint32_t anc_prod,
    uint32_t *next_cbit
) {
    if (cost && cost->temp_anc_needed < 1u) cost->temp_anc_needed = 1u;
    if (cost && cost->scratch_bits_needed < 1u) cost->scratch_bits_needed = 1u;
    uint32_t sb = next_cbit ? *next_cbit : 0u;
    uint32_t a0 = qk__q(qa, 0u), a1 = qk__q(qa, 1u), a2 = qk__q(qa, 2u);
    uint32_t y0 = qk__q(qout, 0u), y1 = qk__q(qout, 1u), y2 = qk__q(qout, 2u);

    qk_and_compute(kb, cost, a1, a2, anc_prod);
    qk_cx(kb, cost, anc_prod, y0);
    qk_and_uncompute(kb, cost, a1, a2, anc_prod, sb);

    qk_and_compute(kb, cost, a0, a1, anc_prod);
    qk_cx(kb, cost, anc_prod, y1);
    qk_and_uncompute(kb, cost, a0, a1, anc_prod, sb);

    qk_and_compute(kb, cost, a0, a2, anc_prod);
    qk_cx(kb, cost, anc_prod, y2);
    qk_and_uncompute(kb, cost, a0, a2, anc_prod, sb);

    qk_cx(kb, cost, a0, y0);
    qk_cx(kb, cost, a1, y0);
    qk_cx(kb, cost, a2, y0);
    qk_cx(kb, cost, a2, y1);
    qk_cx(kb, cost, a1, y2);
    qk_cx(kb, cost, a2, y2);
    if (next_cbit) (*next_cbit)++;
}

static inline void qk_gf8_inv_alg_uncompute_hmr(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *qa,
    const uint32_t *qout,
    uint32_t *next_cbit
) {
    if (cost && cost->scratch_bits_needed < 3u) cost->scratch_bits_needed = 3u;
    uint32_t base = next_cbit ? *next_cbit : 0u;
    uint32_t h0 = base, h1 = base + 1u, h2 = base + 2u;
    uint32_t a0 = qk__q(qa, 0u), a1 = qk__q(qa, 1u), a2 = qk__q(qa, 2u);
    qk_hmr(kb, cost, qk__q(qout, 0u), h0);
    qk_hmr(kb, cost, qk__q(qout, 1u), h1);
    qk_hmr(kb, cost, qk__q(qout, 2u), h2);

    qk_cz_if(kb, cost, a1, a2, h0);
    qk_cz_if(kb, cost, a0, a1, h1);
    qk_cz_if(kb, cost, a0, a2, h2);

    qk_z_if(kb, cost, a0, h0);
    qk_z_if(kb, cost, a1, h0);
    qk_z_if(kb, cost, a1, h2);
    qk_z_if(kb, cost, a2, h0);
    qk_z_if(kb, cost, a2, h1);
    qk_z_if(kb, cost, a2, h2);
    if (next_cbit) *next_cbit = base + 3u;
}

/* GF(2^12) 3-Level Algebraic Tower Inversion: GF(((2^3)^2)^2)
 * Uses 19 temp ancillas unconditional (66 CCX compute / 36 CCX HMR uncompute),
 * or 22 temp ancillas when controlled (69 CCX compute / 39 CCX HMR uncompute). */
typedef struct {
    qcc_bv256_t     to_tower[12];   /* Maps std12 -> (x0, x1) in GF(64)^2 */
    qcc_bv256_t     from_tower[12]; /* Inverse of to_tower */
    qcc_bv256_t     nu_sq_mat[6];   /* Maps x1 -> nu * (x1^2) in composite GF(64) */
    uint32_t        gf64_poly;      /* 0 => composite GF(8)[u]/(u^2+u+1) */
    qir_gf64_rank15_t gf64_plan;
} qir_gf12_tower_inv_plan_t;

static inline bool qir_build_gf12_tower_inv(uint32_t poly12, qir_gf12_tower_inv_plan_t *out) {
    memset(out, 0, sizeof(*out));
    out->gf64_poly = 0u; /* Composite GF(8)[u]/(u^2+u+1) with GF(8)=0xB */
    if (!qir_build_gf64_rank15(0u, &out->gf64_plan)) return false;

    /* 1. Find beta in GF(2^12) with beta^3 + beta + 1 = 0, and u with u^2 + u + 1 = 0 */
    uint16_t beta = 0u, u_elem = 0u;
    for (uint16_t cand = 2u; cand < 4096u; ++cand) {
        uint16_t c2 = (uint16_t)qir__gf2m_mul_u32(cand, cand, 12u, poly12);
        if (!u_elem && (c2 ^ cand ^ 1u) == 0u) u_elem = cand;
        uint16_t c3 = (uint16_t)qir__gf2m_mul_u32(c2, cand, 12u, poly12);
        if (!beta && (c3 ^ cand ^ 1u) == 0u) beta = cand;
        if (u_elem && beta) break;
    }
    if (!beta || !u_elem) return false;

    uint16_t beta2 = (uint16_t)qir__gf2m_mul_u32(beta, beta, 12u, poly12);
    uint16_t gf64_basis[6] = {
        1u, beta, beta2,
        u_elem,
        (uint16_t)qir__gf2m_mul_u32(u_elem, beta, 12u, poly12),
        (uint16_t)qir__gf2m_mul_u32(u_elem, beta2, 12u, poly12)
    };

    /* 2. Find gamma in GF(2^12) \ GF(64) such that gamma^2 + gamma = nu in GF(64) */
    uint16_t gamma = 0u;
    uint8_t nu64 = 0u;
    for (uint16_t g = 2u; g < 4096u && !gamma; ++g) {
        uint16_t g2  = (uint16_t)qir__gf2m_mul_u32(g, g, 12u, poly12);
        uint16_t sum = g2 ^ g;
        for (uint8_t c64 = 1u; c64 < 64u; ++c64) {
            uint16_t val = 0u;
            for (uint32_t b = 0u; b < 6u; ++b) if ((c64 >> b) & 1u) val ^= gf64_basis[b];
            if (val == sum) {
                bool in_sub = false;
                for (uint8_t d64 = 0u; d64 < 64u; ++d64) {
                    uint16_t v2 = 0u;
                    for (uint32_t b = 0u; b < 6u; ++b) if ((d64 >> b) & 1u) v2 ^= gf64_basis[b];
                    if (v2 == g) { in_sub = true; break; }
                }
                if (!in_sub) { gamma = g; nu64 = c64; break; }
            }
        }
    }
    if (!gamma) return false;

    uint16_t basis12[12];
    for (uint32_t i = 0u; i < 6u; ++i) {
        basis12[i]      = gf64_basis[i];
        basis12[6u + i] = (uint16_t)qir__gf2m_mul_u32(gf64_basis[i], gamma, 12u, poly12);
    }

    uint16_t mat[12], inv[12];
    for (uint32_t r = 0u; r < 12u; ++r) {
        uint16_t row = 0u;
        for (uint32_t c = 0u; c < 12u; ++c) if ((basis12[c] >> r) & 1u) row |= (uint16_t)(1u << c);
        mat[r] = row;
        inv[r] = (uint16_t)(1u << r);
    }
    for (uint32_t col = 0u; col < 12u; ++col) {
        uint32_t piv = col;
        while (piv < 12u && (((mat[piv] >> col) & 1u) == 0u)) ++piv;
        if (piv == 12u) return false;
        if (piv != col) {
            uint16_t tm = mat[col]; mat[col] = mat[piv]; mat[piv] = tm;
            uint16_t ti = inv[col]; inv[col] = inv[piv]; inv[piv] = ti;
        }
        for (uint32_t r = 0u; r < 12u; ++r) {
            if (r != col && ((mat[r] >> col) & 1u)) { mat[r] ^= mat[col]; inv[r] ^= inv[col]; }
        }
    }
    for (uint32_t r = 0u; r < 12u; ++r) {
        out->to_tower[r].w[0] = inv[r];
        uint16_t frow = 0u;
        for (uint32_t c = 0u; c < 12u; ++c) if ((basis12[c] >> r) & 1u) frow |= (uint16_t)(1u << c);
        out->from_tower[r].w[0] = frow;
    }

    /* nu * x1^2 in composite GF(64) linear rows */
    for (uint32_t c = 0u; c < 6u; ++c) {
        uint8_t b  = (uint8_t)(1u << c);
        uint8_t b2 = (uint8_t)qir__gf2m_mul_u32(b, b, 6u, 0u);
        uint8_t v  = (uint8_t)qir__gf2m_mul_u32(nu64, b2, 6u, 0u);
        for (uint32_t r = 0u; r < 6u; ++r) {
            if ((v >> r) & 1u) qcc_bv_set(out->nu_sq_mat[r].w, c);
        }
    }
    return true;
}

static inline void qk_gf12_inv_tower(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *qa,    /* [12], transformed in-place then restored cleanly */
    const uint32_t *qout,  /* [12], starts clean 0 */
    uint32_t ctrl_q,       /* UINT32_MAX if unconditional, else control qubit */
    const uint32_t *anc,   /* [19] or [22] if ctrl_q != UINT32_MAX */
    uint32_t *next_cbit,
    const qir_gf12_tower_inv_plan_t *plan
) {
    bool has_ctrl = (ctrl_q != UINT32_MAX);
    uint32_t need_anc = has_ctrl ? 22u : 19u;
    if (cost && cost->temp_anc_needed < need_anc) cost->temp_anc_needed = need_anc;
    if (cost && cost->scratch_bits_needed < 48u) cost->scratch_bits_needed = 48u;

    /* 1. Convert input qa into tower basis (x0 = qa[0..5], x1 = qa[6..11]) in-place */
    qk_linf2_inplace(kb, cost, qa, 12u, plan->to_tower, NULL);

    uint32_t qx0[6], qx1[6], qy0[6], qy1[6];
    uint32_t delta[6], d_norm[3], d_inv[3], dinv[6], d_inv_eff[3];
    for (uint32_t i = 0u; i < 6u; ++i) {
        qx0[i]   = qk__q(qa, i);
        qx1[i]   = qk__q(qa, 6u + i);
        qy0[i]   = qk__q(qout, i);
        qy1[i]   = qk__q(qout, 6u + i);
        delta[i] = qk__q(anc, i);
        dinv[i]  = qk__q(anc, 12u + i);
        qk_r(kb, cost, delta[i]);
        qk_r(kb, cost, dinv[i]);
    }
    for (uint32_t i = 0u; i < 3u; ++i) {
        d_norm[i] = qk__q(anc, 6u + i);
        d_inv[i]  = qk__q(anc, 9u + i);
        qk_r(kb, cost, d_norm[i]);
        qk_r(kb, cost, d_inv[i]);
        if (has_ctrl) {
            d_inv_eff[i] = qk__q(anc, 19u + i);
            qk_r(kb, cost, d_inv_eff[i]);
        } else {
            d_inv_eff[i] = d_inv[i];
        }
    }
    uint32_t anc_prod = qk__q(anc, 18u);

    /* 2. x0' = x0 ^ x1; delta = x0'*x1 + x0'^2 + nu*x1^2 */
    for (uint32_t i = 0u; i < 6u; ++i) qk_cx(kb, cost, qx1[i], qx0[i]);
    qk_gf64_mul_r15_compute(kb, cost, qx0, qx1, delta, anc_prod, next_cbit, &plan->gf64_plan);
    for (uint32_t c = 0u; c < 6u; ++c) {
        uint8_t sq = (uint8_t)qir__gf2m_mul_u32(1u << c, 1u << c, 6u, 0u);
        for (uint32_t r = 0u; r < 6u; ++r) if ((sq >> r) & 1u) qk_cx(kb, cost, qx0[c], delta[r]);
    }
    qk_linf2_out_of_place(kb, cost, qx1, 6u, delta, 6u, plan->nu_sq_mat, NULL);

    /* 3. Algebraic GF(64) inverse: delta = (d0, d1) over GF(8)[u]/(u^2+u+1) */
    for (uint32_t i = 0u; i < 3u; ++i) qk_cx(kb, cost, delta[3u + i], delta[i]); /* d0' = d0 ^ d1 */
    qk_gf8_mul_r6_compute(kb, cost, delta, delta + 3u, d_norm, anc_prod, next_cbit);
    qk_gf8_sq_xor(kb, cost, delta, d_norm);
    qk_gf8_sq_xor(kb, cost, delta + 3u, d_norm);
    qk_gf8_inv_alg_compute(kb, cost, d_norm, d_inv, anc_prod, next_cbit);

    const uint32_t *d_use = d_inv;
    if (has_ctrl) {
        for (uint32_t i = 0u; i < 3u; ++i) qk_and_compute(kb, cost, ctrl_q, d_inv[i], d_inv_eff[i]);
        d_use = d_inv_eff;
    }

    qk_gf8_mul_r6_compute(kb, cost, delta, d_use, dinv, anc_prod, next_cbit);
    qk_gf8_mul_r6_compute(kb, cost, delta + 3u, d_use, dinv + 3u, anc_prod, next_cbit);

    /* 4. y0 = x0' * dinv, y1 = x1 * dinv */
    qk_gf64_mul_r15_compute(kb, cost, qx0, dinv, qy0, anc_prod, next_cbit, &plan->gf64_plan);
    qk_gf64_mul_r15_compute(kb, cost, qx1, dinv, qy1, anc_prod, next_cbit, &plan->gf64_plan);

    /* 5. Uncompute dinv, d_inv_eff, d_inv, d_norm, delta via 0-CCX HMR duals */
    qk_gf8_mul_r6_uncompute_hmr(kb, cost, delta + 3u, d_use, dinv + 3u, next_cbit);
    qk_gf8_mul_r6_uncompute_hmr(kb, cost, delta, d_use, dinv, next_cbit);
    if (has_ctrl) {
        uint32_t sb_eff = next_cbit ? (*next_cbit)++ : 0u;
        for (uint32_t i = 3u; i > 0u; --i) {
            qk_and_uncompute(kb, cost, ctrl_q, d_inv[i - 1u], d_inv_eff[i - 1u], sb_eff);
        }
    }
    qk_gf8_inv_alg_uncompute_hmr(kb, cost, d_norm, d_inv, next_cbit);
    qk_gf8_sq_xor(kb, cost, delta + 3u, d_norm);
    qk_gf8_sq_xor(kb, cost, delta, d_norm);
    qk_gf8_mul_r6_uncompute_hmr(kb, cost, delta, delta + 3u, d_norm, next_cbit);
    for (uint32_t i = 0u; i < 3u; ++i) qk_cx(kb, cost, delta[3u + i], delta[i]);

    qk_linf2_out_of_place(kb, cost, qx1, 6u, delta, 6u, plan->nu_sq_mat, NULL);
    for (uint32_t c = 0u; c < 6u; ++c) {
        uint8_t sq = (uint8_t)qir__gf2m_mul_u32(1u << c, 1u << c, 6u, 0u);
        for (uint32_t r = 0u; r < 6u; ++r) if ((sq >> r) & 1u) qk_cx(kb, cost, qx0[c], delta[r]);
    }
    qk_gf64_mul_r15_uncompute_hmr(kb, cost, qx0, qx1, delta, next_cbit, &plan->gf64_plan);

    /* 6. Restore x0 = x0' ^ x1 and convert qa, qout back to std12 */
    for (uint32_t i = 0u; i < 6u; ++i) qk_cx(kb, cost, qx1[i], qx0[i]);
    qk_linf2_inplace(kb, cost, qa, 12u, plan->from_tower, NULL);
    qk_linf2_inplace(kb, cost, qout, 12u, plan->from_tower, NULL);
}

static inline void qk_gf12_inv_tower_uncompute(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *qa,    /* [12], transformed in-place then restored cleanly */
    const uint32_t *qout,  /* [12], cleared back to |0>^12 */
    uint32_t ctrl_q,       /* UINT32_MAX if unconditional, else control qubit */
    const uint32_t *anc,   /* [19] or [22] if ctrl_q != UINT32_MAX */
    uint32_t *next_cbit,
    const qir_gf12_tower_inv_plan_t *plan
) {
    bool has_ctrl = (ctrl_q != UINT32_MAX);
    uint32_t need_anc = has_ctrl ? 22u : 19u;
    if (cost && cost->temp_anc_needed < need_anc) cost->temp_anc_needed = need_anc;
    if (cost && cost->scratch_bits_needed < 48u) cost->scratch_bits_needed = 48u;

    qk_linf2_inplace(kb, cost, qa, 12u, plan->to_tower, NULL);
    qk_linf2_inplace(kb, cost, qout, 12u, plan->to_tower, NULL);

    uint32_t qx0[6], qx1[6], qy0[6], qy1[6];
    uint32_t delta[6], d_norm[3], d_inv[3], dinv[6], d_inv_eff[3];
    for (uint32_t i = 0u; i < 6u; ++i) {
        qx0[i]   = qk__q(qa, i);
        qx1[i]   = qk__q(qa, 6u + i);
        qy0[i]   = qk__q(qout, i);
        qy1[i]   = qk__q(qout, 6u + i);
        delta[i] = qk__q(anc, i);
        dinv[i]  = qk__q(anc, 12u + i);
        qk_r(kb, cost, delta[i]);
        qk_r(kb, cost, dinv[i]);
    }
    for (uint32_t i = 0u; i < 3u; ++i) {
        d_norm[i] = qk__q(anc, 6u + i);
        d_inv[i]  = qk__q(anc, 9u + i);
        qk_r(kb, cost, d_norm[i]);
        qk_r(kb, cost, d_inv[i]);
        if (has_ctrl) {
            d_inv_eff[i] = qk__q(anc, 19u + i);
            qk_r(kb, cost, d_inv_eff[i]);
        } else {
            d_inv_eff[i] = d_inv[i];
        }
    }
    uint32_t anc_prod = qk__q(anc, 18u);

    for (uint32_t i = 0u; i < 6u; ++i) qk_cx(kb, cost, qx1[i], qx0[i]);
    qk_gf64_mul_r15_compute(kb, cost, qx0, qx1, delta, anc_prod, next_cbit, &plan->gf64_plan);
    for (uint32_t c = 0u; c < 6u; ++c) {
        uint8_t sq = (uint8_t)qir__gf2m_mul_u32(1u << c, 1u << c, 6u, 0u);
        for (uint32_t r = 0u; r < 6u; ++r) if ((sq >> r) & 1u) qk_cx(kb, cost, qx0[c], delta[r]);
    }
    qk_linf2_out_of_place(kb, cost, qx1, 6u, delta, 6u, plan->nu_sq_mat, NULL);

    for (uint32_t i = 0u; i < 3u; ++i) qk_cx(kb, cost, delta[3u + i], delta[i]);
    qk_gf8_mul_r6_compute(kb, cost, delta, delta + 3u, d_norm, anc_prod, next_cbit);
    qk_gf8_sq_xor(kb, cost, delta, d_norm);
    qk_gf8_sq_xor(kb, cost, delta + 3u, d_norm);
    qk_gf8_inv_alg_compute(kb, cost, d_norm, d_inv, anc_prod, next_cbit);

    const uint32_t *d_use = d_inv;
    if (has_ctrl) {
        for (uint32_t i = 0u; i < 3u; ++i) qk_and_compute(kb, cost, ctrl_q, d_inv[i], d_inv_eff[i]);
        d_use = d_inv_eff;
    }

    qk_gf8_mul_r6_compute(kb, cost, delta, d_use, dinv, anc_prod, next_cbit);
    qk_gf8_mul_r6_compute(kb, cost, delta + 3u, d_use, dinv + 3u, anc_prod, next_cbit);

    /* Uncompute qout (qy0, qy1) with 0 CCX via HMR dual! */
    qk_gf64_mul_r15_uncompute_hmr(kb, cost, qx0, dinv, qy0, next_cbit, &plan->gf64_plan);
    qk_gf64_mul_r15_uncompute_hmr(kb, cost, qx1, dinv, qy1, next_cbit, &plan->gf64_plan);

    qk_gf8_mul_r6_uncompute_hmr(kb, cost, delta + 3u, d_use, dinv + 3u, next_cbit);
    qk_gf8_mul_r6_uncompute_hmr(kb, cost, delta, d_use, dinv, next_cbit);
    if (has_ctrl) {
        uint32_t sb_eff = next_cbit ? (*next_cbit)++ : 0u;
        for (uint32_t i = 3u; i > 0u; --i) {
            qk_and_uncompute(kb, cost, ctrl_q, d_inv[i - 1u], d_inv_eff[i - 1u], sb_eff);
        }
    }
    qk_gf8_inv_alg_uncompute_hmr(kb, cost, d_norm, d_inv, next_cbit);
    qk_gf8_sq_xor(kb, cost, delta + 3u, d_norm);
    qk_gf8_sq_xor(kb, cost, delta, d_norm);
    qk_gf8_mul_r6_uncompute_hmr(kb, cost, delta, delta + 3u, d_norm, next_cbit);
    for (uint32_t i = 0u; i < 3u; ++i) qk_cx(kb, cost, delta[3u + i], delta[i]);

    qk_linf2_out_of_place(kb, cost, qx1, 6u, delta, 6u, plan->nu_sq_mat, NULL);
    for (uint32_t c = 0u; c < 6u; ++c) {
        uint8_t sq = (uint8_t)qir__gf2m_mul_u32(1u << c, 1u << c, 6u, 0u);
        for (uint32_t r = 0u; r < 6u; ++r) if ((sq >> r) & 1u) qk_cx(kb, cost, qx0[c], delta[r]);
    }
    qk_gf64_mul_r15_uncompute_hmr(kb, cost, qx0, qx1, delta, next_cbit, &plan->gf64_plan);

    for (uint32_t i = 0u; i < 6u; ++i) qk_cx(kb, cost, qx1[i], qx0[i]);
    qk_linf2_inplace(kb, cost, qa, 12u, plan->from_tower, NULL);
}

/* ============================================================================
 * Kernel 6: Fixed-Point Multiplier, Constant-Mul & Triangle Squarer (Compute & Uncompute)
 * ============================================================================ */
static inline void qk_mul_fxp_ex(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *qa, uint32_t wa, const uint64_t zero_a[QIR_WORDS],
    const uint32_t *qb, uint32_t wb, const uint64_t zero_b[QIR_WORDS],
    const uint32_t *low_q, /* Optional persistent [low_junk] clean aux qubits */
    const uint32_t *qout, uint32_t wout,
    uint32_t shift,
    uint32_t approx_keep_bits,
    const uint32_t *cuccaro_anc, /* [acc_bits] clean temp */
    uint32_t *next_cbit,
    bool uncomp
) {
    uint32_t full_acc_bits = shift + wout;
    uint32_t min_prod_bit  = 0u;
    if (approx_keep_bits > 0u && approx_keep_bits < full_acc_bits) {
        min_prod_bit = full_acc_bits - approx_keep_bits;
        if (min_prod_bit > shift) min_prod_bit = shift;
    }
    uint32_t acc_bits = full_acc_bits - min_prod_bit;
    uint32_t low_junk = shift - min_prod_bit;
    if (acc_bits == 0u) return;
    if (cost && cost->temp_anc_needed < acc_bits) cost->temp_anc_needed = acc_bits;

    uint32_t *acc_q = (uint32_t *)malloc((size_t)acc_bits * sizeof(uint32_t));
    if (!acc_q) abort();
    for (uint32_t j = 0u; j < low_junk; ++j) acc_q[j] = low_q ? low_q[j] : j;
    for (uint32_t j = 0u; j < wout; ++j) acc_q[low_junk + j] = qk__q(qout, j);

    if (!uncomp) {
        for (uint32_t ia = 0u; ia < wa; ++ia) {
            if (zero_a && qcc_bv_get(zero_a, ia)) continue;
            if (ia >= full_acc_bits) break;
            uint32_t start_acc = (ia > min_prod_bit) ? (ia - min_prod_bit) : 0u;
            uint32_t b_start   = (ia < min_prod_bit) ? (min_prod_bit - ia) : 0u;
            if (b_start >= wb) continue;
            uint32_t avail_b = wb - b_start;
            uint32_t tgt_len = acc_bits - start_acc;
            uint32_t k_in    = (avail_b < tgt_len) ? avail_b : tgt_len;
            while (k_in > 0u && zero_b && qcc_bv_get(zero_b, b_start + k_in - 1u)) k_in--;
            if (k_in == 0u) continue;

            uint32_t *b_slice = (uint32_t *)malloc((size_t)k_in * sizeof(uint32_t));
            if (!b_slice) abort();
            for (uint32_t t = 0u; t < k_in; ++t) b_slice[t] = qk__q(qb, b_start + t);
            qk_ctrl_iadd_cuccaro(kb, cost, qk__q(qa, ia), b_slice, k_in, acc_q + start_acc, tgt_len, cuccaro_anc, next_cbit);
            free(b_slice);
        }
    } else {
        for (int32_t iia = (int32_t)wa - 1; iia >= 0; --iia) {
            uint32_t ia = (uint32_t)iia;
            if (zero_a && qcc_bv_get(zero_a, ia)) continue;
            if (ia >= full_acc_bits) continue;
            uint32_t start_acc = (ia > min_prod_bit) ? (ia - min_prod_bit) : 0u;
            uint32_t b_start   = (ia < min_prod_bit) ? (min_prod_bit - ia) : 0u;
            if (b_start >= wb) continue;
            uint32_t avail_b = wb - b_start;
            uint32_t tgt_len = acc_bits - start_acc;
            uint32_t k_in    = (avail_b < tgt_len) ? avail_b : tgt_len;
            while (k_in > 0u && zero_b && qcc_bv_get(zero_b, b_start + k_in - 1u)) k_in--;
            if (k_in == 0u) continue;

            uint32_t *b_slice = (uint32_t *)malloc((size_t)k_in * sizeof(uint32_t));
            if (!b_slice) abort();
            for (uint32_t t = 0u; t < k_in; ++t) b_slice[t] = qk__q(qb, b_start + t);
            qk_ctrl_isub_cuccaro(kb, cost, qk__q(qa, ia), b_slice, k_in, acc_q + start_acc, tgt_len, cuccaro_anc, next_cbit);
            free(b_slice);
        }
    }
    free(acc_q);
}

static inline void qk_mul_fxp(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *qa, uint32_t wa, const uint64_t zero_a[QIR_WORDS],
    const uint32_t *qb, uint32_t wb, const uint64_t zero_b[QIR_WORDS],
    const uint32_t *qout, uint32_t wout,
    uint32_t shift,
    uint32_t approx_keep_bits,
    const uint32_t *anc,
    uint32_t *next_cbit
) {
    uint32_t full_acc_bits = shift + wout;
    uint32_t min_prod_bit  = 0u;
    if (approx_keep_bits > 0u && approx_keep_bits < full_acc_bits) {
        min_prod_bit = full_acc_bits - approx_keep_bits;
        if (min_prod_bit > shift) min_prod_bit = shift;
    }
    uint32_t low_junk = shift - min_prod_bit;
    uint32_t acc_bits = full_acc_bits - min_prod_bit;
    if (cost && cost->temp_anc_needed < low_junk + acc_bits) cost->temp_anc_needed = low_junk + acc_bits;
    for (uint32_t j = 0u; j < low_junk; ++j) qk_r(kb, cost, qk__q(anc, j));
    qk_mul_fxp_ex(kb, cost, qa, wa, zero_a, qb, wb, zero_b, anc, qout, wout, shift, approx_keep_bits,
                  anc ? (anc + low_junk) : NULL, next_cbit, false);
}

static inline void qk_mul_const_fxp(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *qa, uint32_t wa, const uint64_t zero_a[QIR_WORDS],
    const uint64_t c_words[QIR_WORDS], uint32_t wc,
    const uint32_t *low_q,
    const uint32_t *qout, uint32_t wout,
    uint32_t shift,
    uint32_t approx_keep_bits,
    const uint32_t *gidney_anc,
    uint32_t *next_cbit,
    bool uncomp
) {
    uint32_t full_acc_bits = shift + wout;
    uint32_t min_prod_bit  = 0u;
    if (approx_keep_bits > 0u && approx_keep_bits < full_acc_bits) {
        min_prod_bit = full_acc_bits - approx_keep_bits;
        if (min_prod_bit > shift) min_prod_bit = shift;
    }
    uint32_t acc_bits = full_acc_bits - min_prod_bit;
    uint32_t low_junk = shift - min_prod_bit;
    if (acc_bits == 0u) return;
    if (cost && cost->temp_anc_needed < acc_bits) cost->temp_anc_needed = acc_bits;

    uint32_t *acc_q = (uint32_t *)malloc((size_t)acc_bits * sizeof(uint32_t));
    if (!acc_q) abort();
    for (uint32_t j = 0u; j < low_junk; ++j) acc_q[j] = low_q ? low_q[j] : j;
    for (uint32_t j = 0u; j < wout; ++j) acc_q[low_junk + j] = qk__q(qout, j);

    uint32_t lim = (wc < full_acc_bits) ? wc : full_acc_bits;
    for (uint32_t step = 0u; step < lim; ++step) {
        uint32_t ib = uncomp ? (lim - 1u - step) : step;
        if (!qcc_bv_get(c_words, ib)) continue;
        uint32_t start_acc = (ib > min_prod_bit) ? (ib - min_prod_bit) : 0u;
        uint32_t a_start   = (ib < min_prod_bit) ? (min_prod_bit - ib) : 0u;
        if (a_start >= wa) continue;
        uint32_t avail_a = wa - a_start;
        uint32_t tgt_len = acc_bits - start_acc;
        uint32_t k_in    = (avail_a < tgt_len) ? avail_a : tgt_len;
        while (k_in > 0u && zero_a && qcc_bv_get(zero_a, a_start + k_in - 1u)) k_in--;
        uint32_t lead = 0u;
        while (lead < k_in && zero_a && qcc_bv_get(zero_a, a_start + lead)) lead++;
        a_start   += lead;
        start_acc += lead;
        k_in      -= lead;
        tgt_len   -= lead;
        if (k_in == 0u || tgt_len == 0u) continue;

        uint32_t *a_slice = (uint32_t *)malloc((size_t)k_in * sizeof(uint32_t));
        if (!a_slice) abort();
        for (uint32_t t = 0u; t < k_in; ++t) a_slice[t] = qk__q(qa, a_start + t);
        if (!uncomp) {
            qk_gidney_iadd_asym(kb, cost, a_slice, k_in, acc_q + start_acc, tgt_len, gidney_anc, next_cbit);
        } else {
            qk_gidney_isub_asym(kb, cost, a_slice, k_in, acc_q + start_acc, tgt_len, gidney_anc, next_cbit);
        }
        free(a_slice);
    }
    free(acc_q);
}

/* R9 Algebraic Disjoint-Control Triangle Squarer: 0 a_copy qubits! */
static inline void qk_sqr_fxp_tri(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *qa, uint32_t wa, const uint64_t zero_a[QIR_WORDS],
    const uint32_t *low_q,
    const uint32_t *qout, uint32_t wout,
    uint32_t shift,
    uint32_t approx_keep_bits,
    const uint32_t *cuccaro_anc,
    uint32_t *next_cbit,
    bool uncomp
) {
    uint32_t full_acc_bits = shift + wout;
    uint32_t min_prod_bit  = 0u;
    if (approx_keep_bits > 0u && approx_keep_bits < full_acc_bits) {
        min_prod_bit = full_acc_bits - approx_keep_bits;
        if (min_prod_bit > shift) min_prod_bit = shift;
    }
    uint32_t acc_bits = full_acc_bits - min_prod_bit;
    uint32_t low_junk = shift - min_prod_bit;
    if (acc_bits == 0u) return;
    if (cost && cost->temp_anc_needed < acc_bits) cost->temp_anc_needed = acc_bits;

    uint32_t *acc_q = (uint32_t *)malloc((size_t)acc_bits * sizeof(uint32_t));
    if (!acc_q) abort();
    for (uint32_t j = 0u; j < low_junk; ++j) acc_q[j] = low_q ? low_q[j] : j;
    for (uint32_t j = 0u; j < wout; ++j) acc_q[low_junk + j] = qk__q(qout, j);

    if (!uncomp) {
        for (uint32_t i = 0u; i < wa; ++i) {
            uint32_t pos = 2u * i;
            if (pos >= min_prod_bit && pos < full_acc_bits && !(zero_a && qcc_bv_get(zero_a, i))) {
                qk_cx(kb, cost, qk__q(qa, i), acc_q[pos - min_prod_bit]);
            }
        }
        for (uint32_t i = 0u; i + 1u < wa; ++i) {
            if (zero_a && qcc_bv_get(zero_a, i)) continue;
            uint32_t pos = 2u * i + 2u;
            if (pos >= full_acc_bits) break;
            uint32_t start_acc = (pos > min_prod_bit) ? (pos - min_prod_bit) : 0u;
            uint32_t t_start   = (pos < min_prod_bit) ? (min_prod_bit - pos) : 0u;
            uint32_t off = i + 1u + t_start;
            if (off >= wa) continue;
            uint32_t avail   = wa - off;
            uint32_t tgt_len = acc_bits - start_acc;
            uint32_t k_in    = (avail < tgt_len) ? avail : tgt_len;
            while (k_in > 0u && zero_a && qcc_bv_get(zero_a, off + k_in - 1u)) k_in--;
            if (k_in == 0u) continue;
            uint32_t *slice_q = (uint32_t *)malloc((size_t)k_in * sizeof(uint32_t));
            if (!slice_q) abort();
            for (uint32_t t = 0u; t < k_in; ++t) slice_q[t] = qk__q(qa, off + t);
            qk_ctrl_iadd_cuccaro(kb, cost, qk__q(qa, i), slice_q, k_in, acc_q + start_acc, tgt_len, cuccaro_anc, next_cbit);
            free(slice_q);
        }
    } else {
        for (int32_t ii = (int32_t)wa - 2; ii >= 0; --ii) {
            uint32_t i = (uint32_t)ii;
            if (zero_a && qcc_bv_get(zero_a, i)) continue;
            uint32_t pos = 2u * i + 2u;
            if (pos >= full_acc_bits) continue;
            uint32_t start_acc = (pos > min_prod_bit) ? (pos - min_prod_bit) : 0u;
            uint32_t t_start   = (pos < min_prod_bit) ? (min_prod_bit - pos) : 0u;
            uint32_t off = i + 1u + t_start;
            if (off >= wa) continue;
            uint32_t avail   = wa - off;
            uint32_t tgt_len = acc_bits - start_acc;
            uint32_t k_in    = (avail < tgt_len) ? avail : tgt_len;
            while (k_in > 0u && zero_a && qcc_bv_get(zero_a, off + k_in - 1u)) k_in--;
            if (k_in == 0u) continue;
            uint32_t *slice_q = (uint32_t *)malloc((size_t)k_in * sizeof(uint32_t));
            if (!slice_q) abort();
            for (uint32_t t = 0u; t < k_in; ++t) slice_q[t] = qk__q(qa, off + t);
            qk_ctrl_isub_cuccaro(kb, cost, qk__q(qa, i), slice_q, k_in, acc_q + start_acc, tgt_len, cuccaro_anc, next_cbit);
            free(slice_q);
        }
        for (int32_t ii = (int32_t)wa - 1; ii >= 0; --ii) {
            uint32_t i = (uint32_t)ii;
            uint32_t pos = 2u * i;
            if (pos >= min_prod_bit && pos < full_acc_bits && !(zero_a && qcc_bv_get(zero_a, i))) {
                qk_cx(kb, cost, qk__q(qa, i), acc_q[pos - min_prod_bit]);
            }
        }
    }
    free(acc_q);
}

/* ============================================================================
 * Kernel 7: Phase-Polynomial Emitter (with Linear Phase_CSA_Split Support)
 * ============================================================================ */
static inline void qk_emit_phase_terms(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *flat_in_q,
    uint16_t total_in_bits,
    const qir_phase_term_t *terms,
    uint32_t num_terms,
    const uint32_t *anc, /* needed only if degree >= 4 AND monomials exist */
    uint32_t *next_cbit
) {
    (void)anc;
    (void)next_cbit;
    for (uint32_t t = 0u; t < num_terms; ++t) {
        const qir_phase_term_t *pt = &terms[t];
        if (pt->angle == 0u) continue;
        uint32_t qubits[QIR_MAX_VAL_BITS];
        uint32_t nq = 0u;
        for (uint16_t b = 0u; b < total_in_bits; ++b) {
            if (qcc_bv_get(pt->mask, b)) qubits[nq++] = qk__q(flat_in_q, b);
        }
        if (nq == 0u) continue;

        flow_u128 pi_turn = (flow_u128)1u << 127u;
        if (!pt->is_xor) {
            if (nq == 1u) {
                qk_z_pow_u128(kb, cost, qubits[0], pt->angle);
            } else if (nq == 2u && pt->angle == pi_turn) {
                qk_cz(kb, cost, qubits[0], qubits[1]);
            } else if (nq == 3u && pt->angle == pi_turn) {
                qk_ccz(kb, cost, qubits[0], qubits[1], qubits[2]);
            } else if (nq >= 2u) {
                if (cost && (nq - 1u) > cost->temp_anc_needed) cost->temp_anc_needed = nq - 1u;
                if (cost && cost->scratch_bits_needed < 1u) cost->scratch_bits_needed = 1u;
                uint32_t sb = next_cbit ? *next_cbit : 0u;
                uint32_t a0 = qk__q(anc, 0u);
                qk_and_compute(kb, cost, qubits[0], qubits[1], a0);
                for (uint32_t i = 2u; i < nq; ++i) {
                    qk_and_compute(kb, cost, qk__q(anc, i - 2u), qubits[i], qk__q(anc, i - 1u));
                }
                qk_z_pow_u128(kb, cost, qk__q(anc, nq - 2u), pt->angle);
                for (uint32_t i = nq - 1u; i >= 2u; --i) {
                    qk_and_uncompute(kb, cost, qk__q(anc, i - 2u), qubits[i], qk__q(anc, i - 1u), sb);
                }
                qk_and_uncompute(kb, cost, qubits[0], qubits[1], a0, sb);
                if (next_cbit) (*next_cbit)++;
            }
        } else {
            /* XOR phase gadget: CNOT ladder into qubits[0], Z_POW, uncompute CNOT ladder (0 CCX!) */
            for (uint32_t i = 1u; i < nq; ++i) qk_cx(kb, cost, qubits[i], qubits[0]);
            qk_z_pow_u128(kb, cost, qubits[0], pt->angle);
            for (uint32_t i = nq - 1u; i >= 1u; --i) qk_cx(kb, cost, qubits[i], qubits[0]);
        }
    }
}

/* ============================================================================
 * Kernel 8: X-Basis Measurement (HMR) & Classical Boolean Vector Kernels
 * ============================================================================ */
static inline void qk_hmr_meas(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *q_in,
    const uint32_t *b_out,
    uint32_t w
) {
    for (uint32_t i = 0u; i < w; ++i) {
        qk_hmr(kb, cost, qk__q(q_in, i), qk__q(b_out, i));
    }
}

static inline void qk_bit_op(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    uint8_t mode, /* 0 = XOR, 1 = AND, 2 = NOT, 3 = OR_REDUCE */
    const uint32_t *a_b,
    const uint32_t *b_b,
    const uint32_t *out_b,
    uint32_t w
) {
    if (mode == 3u) {
        uint32_t dst = qk__q(out_b, 0u);
        qk_bit_store0(kb, cost, dst);
        for (uint32_t i = 0u; i < w; ++i) {
            qk_push_cond(kb, cost, qk__q(a_b, i));
            qk_bit_store1(kb, cost, dst);
            qk_pop_cond(kb, cost);
        }
        return;
    }
    for (uint32_t i = 0u; i < w; ++i) {
        uint32_t dst = qk__q(out_b, i);
        uint32_t sa  = qk__q(a_b, i);
        if (mode == 0u) {
            uint32_t sb = qk__q(b_b, i);
            qk_bit_store0(kb, cost, dst);
            qk_bit_invert_if(kb, cost, dst, sa);
            qk_bit_invert_if(kb, cost, dst, sb);
        } else if (mode == 1u) {
            uint32_t sb = qk__q(b_b, i);
            qk_bit_store0(kb, cost, dst);
            qk_push_cond(kb, cost, sa);
            qk_bit_invert_if(kb, cost, dst, sb);
            qk_pop_cond(kb, cost);
        } else {
            qk_bit_store1(kb, cost, dst);
            qk_bit_invert_if(kb, cost, dst, sa);
        }
    }
}

/* ============================================================================
 * Kernel 9: Bijective Barrel Shifter (SHL_VAR / SHR_VAR), EQ_ZERO & PHASE_DOT
 * ============================================================================ */
static inline uint32_t qk__gcd_u32(uint32_t a, uint32_t b) {
    while (b != 0u) { uint32_t t = a % b; a = b; b = t; }
    return a ? a : 1u;
}

static inline void qk__rot_left_cswap(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    uint32_t ctrl,
    const uint32_t *buf,
    uint32_t W,
    uint32_t s,
    bool reverse
) {
    if (W <= 1u) return;
    s %= W;
    if (s == 0u) return;
    uint32_t g = qk__gcd_u32(s, W);
    uint32_t L = W / g;
    for (uint32_t c = 0u; c < g; ++c) {
        if (!reverse) {
            for (uint32_t t = L - 1u; t >= 1u; --t) {
                uint32_t i0 = (c + (t - 1u) * s) % W;
                uint32_t i1 = (c + t * s) % W;
                qk_cswap(kb, cost, ctrl, buf[i0], buf[i1]);
            }
        } else {
            for (uint32_t t = 1u; t < L; ++t) {
                uint32_t i0 = (c + (t - 1u) * s) % W;
                uint32_t i1 = (c + t * s) % W;
                qk_cswap(kb, cost, ctrl, buf[i0], buf[i1]);
            }
        }
    }
}

static inline void qk_shl_var(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *in_q, uint32_t w_in,
    const uint32_t *sh_q, uint32_t w_sh,
    const uint32_t *out_q, uint32_t out_w,
    const uint32_t *aux_q,
    uint32_t max_shift,
    bool uncomp
) {
    uint32_t W = (w_in + max_shift > out_w) ? (w_in + max_shift) : out_w;
    if (W == 0u || W > 512u) return;
    uint32_t buf[512];
    for (uint32_t i = 0u; i < out_w; ++i) buf[i] = qk__q(out_q, i);
    for (uint32_t i = out_w; i < W; ++i) buf[i] = qk__q(aux_q, i - out_w);

    if (!uncomp) {
        for (uint32_t i = 0u; i < w_in && i < W; ++i) qk_cx(kb, cost, qk__q(in_q, i), buf[i]);
        for (uint32_t k = 0u; k < w_sh && k < 30u; ++k) {
            uint32_t step = UINT32_C(1) << k;
            if (step > max_shift) break;
            if (W > 1u && (step % W) != 0u) {
                qk__rot_left_cswap(kb, cost, qk__q(sh_q, k), buf, W, step % W, false);
            }
        }
    } else {
        for (int32_t ks = (int32_t)w_sh - 1; ks >= 0; --ks) {
            uint32_t k = (uint32_t)ks;
            if (k >= 30u) continue;
            uint32_t step = UINT32_C(1) << k;
            if (step > max_shift) continue;
            if (W > 1u && (step % W) != 0u) {
                qk__rot_left_cswap(kb, cost, qk__q(sh_q, k), buf, W, step % W, true);
            }
        }
        for (uint32_t i = 0u; i < w_in && i < W; ++i) qk_cx(kb, cost, qk__q(in_q, i), buf[i]);
    }
}

static inline void qk_shr_var(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *in_q, uint32_t w_in,
    const uint32_t *sh_q, uint32_t w_sh,
    const uint32_t *out_q, uint32_t out_w,
    const uint32_t *aux_q,
    uint32_t max_shift,
    bool uncomp
) {
    uint32_t base_w = (w_in > out_w) ? w_in : out_w;
    uint32_t W = base_w + max_shift;
    if (W == 0u || W > 512u) return;
    uint32_t buf[512];
    for (uint32_t i = 0u; i < max_shift; ++i) buf[i] = qk__q(aux_q, i);
    for (uint32_t i = 0u; i < out_w; ++i) buf[max_shift + i] = qk__q(out_q, i);
    for (uint32_t i = out_w; i < base_w; ++i) buf[max_shift + i] = qk__q(aux_q, max_shift + (i - out_w));

    if (!uncomp) {
        for (uint32_t i = 0u; i < w_in && (max_shift + i) < W; ++i) {
            qk_cx(kb, cost, qk__q(in_q, i), buf[max_shift + i]);
        }
        for (uint32_t k = 0u; k < w_sh && k < 30u; ++k) {
            uint32_t step = UINT32_C(1) << k;
            if (step > max_shift) break;
            if (W > 1u && (step % W) != 0u) {
                uint32_t rot = (W - (step % W)) % W;
                qk__rot_left_cswap(kb, cost, qk__q(sh_q, k), buf, W, rot, false);
            }
        }
    } else {
        for (int32_t ks = (int32_t)w_sh - 1; ks >= 0; --ks) {
            uint32_t k = (uint32_t)ks;
            if (k >= 30u) continue;
            uint32_t step = UINT32_C(1) << k;
            if (step > max_shift) continue;
            if (W > 1u && (step % W) != 0u) {
                uint32_t rot = (W - (step % W)) % W;
                qk__rot_left_cswap(kb, cost, qk__q(sh_q, k), buf, W, rot, true);
            }
        }
        for (uint32_t i = 0u; i < w_in && (max_shift + i) < W; ++i) {
            qk_cx(kb, cost, qk__q(in_q, i), buf[max_shift + i]);
        }
    }
}

static inline void qk_eq_zero_compute(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *in_q,
    uint32_t w,
    const uint32_t *out_q,
    const uint32_t *anc,
    uint32_t *next_cbit
) {
    uint32_t qo = qk__q(out_q, 0u);
    if (w == 0u) { qk_x(kb, cost, qo); return; }
    if (w == 1u) {
        uint32_t q0 = qk__q(in_q, 0u);
        qk_x(kb, cost, q0);
        qk_cx(kb, cost, q0, qo);
        qk_x(kb, cost, q0);
        return;
    }
    if (w == 2u) {
        uint32_t q0 = qk__q(in_q, 0u), q1 = qk__q(in_q, 1u);
        qk_x(kb, cost, q0);
        qk_x(kb, cost, q1);
        qk_ccx(kb, cost, q0, q1, qo);
        qk_x(kb, cost, q1);
        qk_x(kb, cost, q0);
        return;
    }
    if (cost && cost->temp_anc_needed < w - 2u) cost->temp_anc_needed = w - 2u;
    if (cost && cost->scratch_bits_needed < 2u) cost->scratch_bits_needed = 2u;
    uint32_t sbit = next_cbit ? *next_cbit : 0u;
    if (next_cbit) (*next_cbit)++;

    for (uint32_t i = 0u; i < w; ++i) qk_x(kb, cost, qk__q(in_q, i));
    qk_and_compute(kb, cost, qk__q(in_q, 0u), qk__q(in_q, 1u), qk__q(anc, 0u));
    for (uint32_t i = 2u; i + 1u < w; ++i) {
        qk_and_compute(kb, cost, qk__q(anc, i - 2u), qk__q(in_q, i), qk__q(anc, i - 1u));
    }
    qk_ccx(kb, cost, qk__q(anc, w - 3u), qk__q(in_q, w - 1u), qo);
    for (uint32_t i = w - 2u; i >= 2u; --i) {
        qk_and_uncompute(kb, cost, qk__q(anc, i - 2u), qk__q(in_q, i), qk__q(anc, i - 1u), sbit);
    }
    qk_and_uncompute(kb, cost, qk__q(in_q, 0u), qk__q(in_q, 1u), qk__q(anc, 0u), sbit);
    for (uint32_t i = w; i > 0u; --i) qk_x(kb, cost, qk__q(in_q, i - 1u));
}

static inline void qk_eq_zero_uncompute(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *in_q,
    uint32_t w,
    const uint32_t *out_q,
    const uint32_t *anc,
    uint32_t *next_cbit
) {
    uint32_t qo = qk__q(out_q, 0u);
    if (w == 0u) { qk_x(kb, cost, qo); return; }
    if (w == 1u) {
        uint32_t q0 = qk__q(in_q, 0u);
        qk_x(kb, cost, q0);
        qk_cx(kb, cost, q0, qo);
        qk_x(kb, cost, q0);
        return;
    }
    if (w == 2u) {
        if (cost && cost->scratch_bits_needed < 1u) cost->scratch_bits_needed = 1u;
        uint32_t hbit = next_cbit ? (*next_cbit)++ : 0u;
        uint32_t q0 = qk__q(in_q, 0u), q1 = qk__q(in_q, 1u);
        qk_x(kb, cost, q0);
        qk_x(kb, cost, q1);
        qk_hmr(kb, cost, qo, hbit);
        qk_cz_if(kb, cost, q0, q1, hbit);
        qk_x(kb, cost, q1);
        qk_x(kb, cost, q0);
        return;
    }
    if (cost && cost->temp_anc_needed < w - 2u) cost->temp_anc_needed = w - 2u;
    if (cost && cost->scratch_bits_needed < 2u) cost->scratch_bits_needed = 2u;
    uint32_t base = next_cbit ? *next_cbit : 0u;
    uint32_t sbit = base, hbit = base + 1u;
    if (next_cbit) *next_cbit = base + 2u;

    for (uint32_t i = 0u; i < w; ++i) qk_x(kb, cost, qk__q(in_q, i));
    qk_and_compute(kb, cost, qk__q(in_q, 0u), qk__q(in_q, 1u), qk__q(anc, 0u));
    for (uint32_t i = 2u; i + 1u < w; ++i) {
        qk_and_compute(kb, cost, qk__q(anc, i - 2u), qk__q(in_q, i), qk__q(anc, i - 1u));
    }
    qk_hmr(kb, cost, qo, hbit);
    qk_cz_if(kb, cost, qk__q(anc, w - 3u), qk__q(in_q, w - 1u), hbit);
    for (uint32_t i = w - 2u; i >= 2u; --i) {
        qk_and_uncompute(kb, cost, qk__q(anc, i - 2u), qk__q(in_q, i), qk__q(anc, i - 1u), sbit);
    }
    qk_and_uncompute(kb, cost, qk__q(in_q, 0u), qk__q(in_q, 1u), qk__q(anc, 0u), sbit);
    for (uint32_t i = w; i > 0u; --i) qk_x(kb, cost, qk__q(in_q, i - 1u));
}

static inline void qk_phase_dot(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *q_in,
    const uint32_t *c_b,
    uint32_t w
) {
    for (uint32_t i = 0u; i < w; ++i) {
        qk_z_if(kb, cost, qk__q(q_in, i), qk__q(c_b, i));
    }
}

static inline void qk_phase_dot_gf2m(
    qcc_kmb_circuit_t *kb,
    qir_kernel_cost_t *cost,
    const uint32_t *qa,
    const uint32_t *qb,
    const uint32_t *c_b,
    uint16_t deg,
    uint32_t irred_poly,
    const qir_gf12_rank42_t *r42,
    const qir_gf64_rank15_t *r15,
    uint32_t *next_cbit
) {
    if (cost && cost->scratch_bits_needed < 1u) cost->scratch_bits_needed = 1u;
    uint32_t pbit = next_cbit ? (*next_cbit)++ : 0u;
    if (deg == 12u && r42) {
        for (uint32_t p = 0u; p < QIR_GF12_RANK; ++p) {
            uint16_t lm = r42->l_mask[p], rm = r42->r_mask[p], wm = r42->w_mask[p];
            if (!lm || !rm || !wm) continue;
            qk_bit_store0(kb, cost, pbit);
            for (uint32_t k = 0u; k < 12u; ++k) {
                if ((wm >> k) & 1u) qk_bit_invert_if(kb, cost, pbit, qk__q(c_b, k));
            }
            uint32_t piv_a = (uint32_t)__builtin_ctz((unsigned)lm);
            uint32_t piv_b = (uint32_t)__builtin_ctz((unsigned)rm);
            for (uint32_t i = piv_a + 1u; i < 12u; ++i) {
                if ((lm >> i) & 1u) qk_cx(kb, cost, qk__q(qa, i), qk__q(qa, piv_a));
            }
            for (uint32_t j = piv_b + 1u; j < 12u; ++j) {
                if ((rm >> j) & 1u) qk_cx(kb, cost, qk__q(qb, j), qk__q(qb, piv_b));
            }
            qk_cz_if(kb, cost, qk__q(qa, piv_a), qk__q(qb, piv_b), pbit);
            for (uint32_t j = 11u; j > piv_b; --j) {
                if ((rm >> j) & 1u) qk_cx(kb, cost, qk__q(qb, j), qk__q(qb, piv_b));
            }
            for (uint32_t i = 11u; i > piv_a; --i) {
                if ((lm >> i) & 1u) qk_cx(kb, cost, qk__q(qa, i), qk__q(qa, piv_a));
            }
        }
        return;
    }
    if (deg == 6u && r15) {
        for (uint32_t p = 0u; p < QIR_GF64_RANK; ++p) {
            uint8_t lm = r15->l_mask[p], rm = r15->r_mask[p], wm = r15->w_mask[p];
            if (!lm || !rm || !wm) continue;
            qk_bit_store0(kb, cost, pbit);
            for (uint32_t k = 0u; k < 6u; ++k) {
                if ((wm >> k) & 1u) qk_bit_invert_if(kb, cost, pbit, qk__q(c_b, k));
            }
            uint32_t piv_a = (uint32_t)__builtin_ctz((unsigned)lm);
            uint32_t piv_b = (uint32_t)__builtin_ctz((unsigned)rm);
            for (uint32_t i = piv_a + 1u; i < 6u; ++i) {
                if ((lm >> i) & 1u) qk_cx(kb, cost, qk__q(qa, i), qk__q(qa, piv_a));
            }
            for (uint32_t j = piv_b + 1u; j < 6u; ++j) {
                if ((rm >> j) & 1u) qk_cx(kb, cost, qk__q(qb, j), qk__q(qb, piv_b));
            }
            qk_cz_if(kb, cost, qk__q(qa, piv_a), qk__q(qb, piv_b), pbit);
            for (uint32_t j = 5u; j > piv_b; --j) {
                if ((rm >> j) & 1u) qk_cx(kb, cost, qk__q(qb, j), qk__q(qb, piv_b));
            }
            for (uint32_t i = 5u; i > piv_a; --i) {
                if ((lm >> i) & 1u) qk_cx(kb, cost, qk__q(qa, i), qk__q(qa, piv_a));
            }
        }
        return;
    }
    for (uint16_t j = 0u; j < deg; ++j) {
        for (uint16_t i = 0u; i < deg; ++i) {
            uint32_t elem = qir__gf2m_mul_u32(UINT32_C(1) << i, UINT32_C(1) << j, deg, irred_poly);
            if (elem == 0u) continue;
            qk_bit_store0(kb, cost, pbit);
            for (uint16_t k = 0u; k < deg; ++k) {
                if ((elem >> k) & 1u) qk_bit_invert_if(kb, cost, pbit, qk__q(c_b, k));
            }
            qk_cz_if(kb, cost, qk__q(qa, i), qk__q(qb, j), pbit);
        }
    }
}

#ifdef __cplusplus
}
#endif

#endif /* QIR_KERNELS_H */
