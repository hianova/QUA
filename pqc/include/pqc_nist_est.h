/*
 * pqc_nist_est.h — NIST PQC / Post-Quantum Security Resource Estimator (MAXDEPTH 2^40/2^64/2^96)
 * Part of QUA / PQC Quantum Cryptanalysis Suite.
 * Pure C17 & C++20 compatible, zero <math.h>, zero external dependencies.
 */
#ifndef PQC_NIST_EST_H
#define PQC_NIST_EST_H

#include "qcc_kmb.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t maxdepth_log2;           /* 40, 64, or 96 */
    uint32_t log2_instances_q16;      /* ceil-continuous log2(P) in Q16.16 */
    uint32_t log2_total_qubits_q16;   /* log2(P * Q_oracle) in Q16.16 */
    uint32_t log2_total_toffoli_q16;  /* log2(P * R_inst * CCX_iter) in Q16.16 */
    uint32_t log2_total_t_q16;        /* log2(P * R_inst * T_eq_iter) in Q16.16 */
    uint32_t log2_g_times_d_q16;      /* log2(total_T * actual_depth) in Q16.16 */
} pqc_maxdepth_row_t;

typedef struct {
    char               name[64];
    char               family[48];
    uint32_t           search_bits;
    qcc_kmb_stats_t    oracle_stats;
    uint64_t           t_eq_gidney;
    uint64_t           t_eq_unitary;
    uint64_t           per_iter_toffoli;
    uint64_t           per_iter_t_eq;
    uint64_t           per_iter_depth;
    uint32_t           log2_serial_iters_q16;
    uint32_t           log2_serial_toffoli_q16;
    uint32_t           log2_serial_t_q16;
    uint32_t           log2_serial_depth_q16;
    pqc_maxdepth_row_t maxdepth[3];
    uint32_t           nist_category_equiv; /* 0 (< Cat 1), 1 (AES-128), 3 (AES-192), 5 (AES-256) */
} pqc_nist_report_t;

/* Exact zero-<math.h> binary log2(x) in Q16.16 fixed-point using 128-bit squaring */
static inline uint32_t pqc_log2_q16_u128(flow_u128 x) {
    if (x <= 1u) return 0u;
    uint64_t hi = (uint64_t)(x >> 64u);
    uint64_t lo = (uint64_t)x;
    uint32_t e = hi ? (64u + 63u - (uint32_t)__builtin_clzll((unsigned long long)hi))
                    : (63u - (uint32_t)__builtin_clzll((unsigned long long)lo));

    flow_u128 m = (e <= 62u) ? (x << (62u - e)) : (x >> (e - 62u));
    const flow_u128 two_q62 = (flow_u128)1u << 63u;
    uint32_t frac = 0u;
    for (uint32_t i = 0u; i < 24u; ++i) {
        m = (m * m) >> 62u;
        frac <<= 1u;
        if (m >= two_q62) {
            frac |= 1u;
            m >>= 1u;
        }
    }
    uint32_t rounded = (frac + 128u) >> 8u; /* Round Q24 -> Q16 */
    return (e << 16u) + rounded;
}

static inline uint32_t pqc_log2_q16_u64(uint64_t x) {
    return pqc_log2_q16_u128((flow_u128)x);
}

static inline int32_t pqc_log2_pi_over_4_q16(void) {
    /* log2(pi / 4) computed purely in integer arithmetic (~ -0.348504 * 65536 = -22839) */
    uint32_t l_pi = pqc_log2_q16_u64(UINT64_C(3141592653589793238));
    uint32_t l_4  = pqc_log2_q16_u64(UINT64_C(4000000000000000000));
    return (int32_t)l_pi - (int32_t)l_4;
}

static inline void pqc_format_pow2_q16(uint32_t log2_q16, char *buf, size_t cap) {
    uint32_t int_p = log2_q16 >> 16u;
    uint32_t dec_p = ((log2_q16 & 0xFFFFu) * 100u + 32768u) >> 16u;
    if (dec_p >= 100u) {
        int_p++;
        dec_p -= 100u;
    }
    snprintf(buf, cap, "2^%u.%02u", (unsigned)int_p, (unsigned)dec_p);
}

static inline uint32_t pqc_classify_nist_cat(uint32_t log2_t_md64_q16) {
    /* NIST Category thresholds against AES-128/192/256 T-gate complexity under MAXDEPTH = 2^64 */
    if (log2_t_md64_q16 >= (214u << 16u)) return 5u;
    if (log2_t_md64_q16 >= (150u << 16u)) return 3u;
    if (log2_t_md64_q16 >= (84u  << 16u)) return 1u;
    return 0u;
}

/* Compute full Grover serial & Zalka MAXDEPTH parallel NIST resource report */
static inline void pqc_nist_estimate(
    const char *name,
    const char *family,
    uint32_t search_bits,
    const qcc_kmb_stats_t *st,
    pqc_nist_report_t *out
) {
    memset(out, 0, sizeof(*out));
    if (name)   snprintf(out->name, sizeof(out->name), "%s", name);
    if (family) snprintf(out->family, sizeof(out->family), "%s", family);
    out->search_bits  = search_bits;
    out->oracle_stats = *st;

    uint64_t extra_rot_t = st->non_clifford_z_pow * 40u;
    out->t_eq_gidney  = 4u * st->toffoli + st->t_like + extra_rot_t;
    out->t_eq_unitary = 7u * st->toffoli + st->t_like + extra_rot_t;

    /* Diffusion operator D_diff over k search bits uses (k - 2) Toffolis and O(log k) depth */
    uint64_t diff_ccx = (search_bits >= 3u) ? (uint64_t)(search_bits - 2u) : 0u;
    uint64_t diff_t_depth = 0u;
    if (search_bits >= 3u) {
        uint32_t lg = 31u - (uint32_t)__builtin_clz((unsigned)search_bits);
        diff_t_depth = (uint64_t)(2u * lg + 1u);
    }

    uint64_t base_nd = st->toffoli_depth ? st->toffoli_depth : st->circuit_depth;
    if (base_nd == 0u) base_nd = st->toffoli ? st->toffoli : 1u;

    out->per_iter_toffoli = st->toffoli + diff_ccx;
    out->per_iter_t_eq    = out->t_eq_gidney + 4u * diff_ccx;
    out->per_iter_depth   = base_nd + diff_t_depth;

    int32_t half_k_q16 = (int32_t)(search_bits << 15u);
    int32_t ser_it_s   = half_k_q16 + pqc_log2_pi_over_4_q16();
    out->log2_serial_iters_q16 = (ser_it_s > 0) ? (uint32_t)ser_it_s : 0u;

    uint32_t l_ccx = pqc_log2_q16_u64(out->per_iter_toffoli ? out->per_iter_toffoli : 1u);
    uint32_t l_t   = pqc_log2_q16_u64(out->per_iter_t_eq ? out->per_iter_t_eq : 1u);
    uint32_t l_d   = pqc_log2_q16_u64(out->per_iter_depth ? out->per_iter_depth : 1u);
    uint32_t l_q   = pqc_log2_q16_u64(st->num_qubits ? st->num_qubits : 1u);

    out->log2_serial_toffoli_q16 = out->log2_serial_iters_q16 + l_ccx;
    out->log2_serial_t_q16       = out->log2_serial_iters_q16 + l_t;
    out->log2_serial_depth_q16   = out->log2_serial_iters_q16 + l_d;

    static const uint32_t MD[3] = {40u, 64u, 96u};
    for (uint32_t mi = 0u; mi < 3u; ++mi) {
        uint32_t m_log2 = MD[mi];
        uint32_t m_q16  = m_log2 << 16u;
        uint32_t lp_q16 = 0u;
        if (out->log2_serial_depth_q16 > m_q16) {
            lp_q16 = 2u * (out->log2_serial_depth_q16 - m_q16);
        }
        uint32_t half_lp_q16 = lp_q16 >> 1u;
        pqc_maxdepth_row_t *row = &out->maxdepth[mi];
        row->maxdepth_log2          = m_log2;
        row->log2_instances_q16     = lp_q16;
        row->log2_total_qubits_q16  = lp_q16 + l_q;
        row->log2_total_toffoli_q16 = out->log2_serial_toffoli_q16 + half_lp_q16;
        row->log2_total_t_q16       = out->log2_serial_t_q16 + half_lp_q16;
        uint32_t actual_d_q16       = (out->log2_serial_depth_q16 > m_q16) ? m_q16 : out->log2_serial_depth_q16;
        row->log2_g_times_d_q16     = row->log2_total_t_q16 + actual_d_q16;
    }

    out->nist_category_equiv = pqc_classify_nist_cat(out->maxdepth[1].log2_total_t_q16);
}

#ifdef __cplusplus
}
#endif

#endif /* PQC_NIST_EST_H */
