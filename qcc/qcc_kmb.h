/*
 * qcc_kmb.h — Bidirectional Kickmix v2 (.kmb) Codec, Unified Cursor & Exact Stats
 * Part of QCC (Quantum Circuit Compiler) for FLOW v7.0.1.
 * Pure C17, zero external dependencies.
 */
#ifndef QCC_KMB_H
#define QCC_KMB_H

#include "qcc_arena.h"

#ifdef __cplusplus
extern "C" {
#endif

static const uint8_t QCC_KMB_MAGIC[16] = {
    0xd7, 0x50, 0xc7, 0xd5, 0xc3, 0x29, 0xd3, 0x26,
    0xe3, 0xcc, 0x9f, 0x68, 0x34, 0xf2, 0xb8, 0xbf
};

/* Kickmix OpType encoding bits (matches kickmix/circuit/op_type.h) */
enum {
    QCC_OP_FLAG_1Q   = 0x10u,
    QCC_OP_FLAG_2Q   = 0x20u,
    QCC_OP_FLAG_3Q   = 0x30u,
    QCC_OP_FLAG_1B   = 0x40u,
    QCC_OP_FLAG_COND = 0x80u
};

typedef enum {
    QCC_OP_NEG                  = 0x00u,
    QCC_OP_POP_CONDITION        = 0x01u,
    QCC_OP_INFO_TAG             = 0x06u,
    QCC_OP_DEBUG_PRINT_EMPTY    = 0x07u,

    QCC_OP_BIT_STORE0           = 0x00u | QCC_OP_FLAG_1B,
    QCC_OP_BIT_STORE1           = 0x01u | QCC_OP_FLAG_1B,
    QCC_OP_BIT_INVERT           = 0x02u | QCC_OP_FLAG_1B,
    QCC_OP_DEBUG_PRINT_C        = 0x07u | QCC_OP_FLAG_1B,

    QCC_OP_Z                    = 0x00u | QCC_OP_FLAG_1Q,
    QCC_OP_X                    = 0x01u | QCC_OP_FLAG_1Q,
    QCC_OP_R                    = 0x02u | QCC_OP_FLAG_1Q,
    QCC_OP_Z_POW                = 0x03u | QCC_OP_FLAG_1Q,
    QCC_OP_DEBUG_PRINT_Q        = 0x07u | QCC_OP_FLAG_1Q,

    QCC_OP_HMR                  = 0x00u | QCC_OP_FLAG_1Q | QCC_OP_FLAG_1B,

    QCC_OP_CZ                   = 0x00u | QCC_OP_FLAG_2Q,
    QCC_OP_CX                   = 0x01u | QCC_OP_FLAG_2Q,
    QCC_OP_SWAP                 = 0x02u | QCC_OP_FLAG_2Q,

    QCC_OP_CCZ                  = 0x00u | QCC_OP_FLAG_3Q,
    QCC_OP_CCX                  = 0x01u | QCC_OP_FLAG_3Q,
    QCC_OP_CSWAP                = 0x02u | QCC_OP_FLAG_3Q,

    QCC_OP_NEG_IF               = 0x00u | QCC_OP_FLAG_COND,
    QCC_OP_PUSH_CONDITION       = 0x01u | QCC_OP_FLAG_COND,
    QCC_OP_DEBUG_PRINT_EMPTY_IF = 0x07u | QCC_OP_FLAG_COND,

    QCC_OP_BIT_STORE0_IF        = 0x00u | QCC_OP_FLAG_1B | QCC_OP_FLAG_COND,
    QCC_OP_BIT_STORE1_IF        = 0x01u | QCC_OP_FLAG_1B | QCC_OP_FLAG_COND,
    QCC_OP_BIT_INVERT_IF        = 0x02u | QCC_OP_FLAG_1B | QCC_OP_FLAG_COND,
    QCC_OP_DEBUG_PRINT_C_IF     = 0x07u | QCC_OP_FLAG_1B | QCC_OP_FLAG_COND,

    QCC_OP_Z_IF                 = 0x00u | QCC_OP_FLAG_1Q | QCC_OP_FLAG_COND,
    QCC_OP_X_IF                 = 0x01u | QCC_OP_FLAG_1Q | QCC_OP_FLAG_COND,
    QCC_OP_R_IF                 = 0x02u | QCC_OP_FLAG_1Q | QCC_OP_FLAG_COND,
    QCC_OP_Z_POW_IF             = 0x03u | QCC_OP_FLAG_1Q | QCC_OP_FLAG_COND,
    QCC_OP_DEBUG_PRINT_Q_IF     = 0x07u | QCC_OP_FLAG_1Q | QCC_OP_FLAG_COND,

    QCC_OP_HMR_IF               = 0x00u | QCC_OP_FLAG_1Q | QCC_OP_FLAG_1B | QCC_OP_FLAG_COND,

    QCC_OP_CZ_IF                = 0x00u | QCC_OP_FLAG_2Q | QCC_OP_FLAG_COND,
    QCC_OP_CX_IF                = 0x01u | QCC_OP_FLAG_2Q | QCC_OP_FLAG_COND,
    QCC_OP_SWAP_IF              = 0x02u | QCC_OP_FLAG_2Q | QCC_OP_FLAG_COND,

    QCC_OP_CCZ_IF               = 0x00u | QCC_OP_FLAG_3Q | QCC_OP_FLAG_COND,
    QCC_OP_CCX_IF               = 0x01u | QCC_OP_FLAG_3Q | QCC_OP_FLAG_COND,
    QCC_OP_CSWAP_IF             = 0x02u | QCC_OP_FLAG_3Q | QCC_OP_FLAG_COND
} qcc_kmb_op_t;

typedef enum {
    QCC_PKT_REGISTER_TARGETS = 0,
    QCC_PKT_OPERATION_TYPES  = 1,
    QCC_PKT_Q0_3Q            = 2,
    QCC_PKT_Q1_3Q            = 3,
    QCC_PKT_Q2_3Q            = 4,
    QCC_PKT_Q0_2Q            = 5,
    QCC_PKT_Q1_2Q            = 6,
    QCC_PKT_Q0_1Q            = 7,
    QCC_PKT_B0_1B            = 8,
    QCC_PKT_B0_COND          = 9,
    QCC_PKT_ANGLE            = 10,
    QCC_PKT_INFO_TAG_DATA    = 11
} qcc_kmb_pkt_t;

typedef struct {
    char      name[64];
    uint32_t  name_len;
    uint32_t *targets; /* Bit 31 set => qubit (target & 0x7FFFFFFF), else bit */
    uint32_t  len;
    uint32_t  cap;
} qcc_kmb_reg_t;

typedef struct {
    uint32_t num_qubits;
    uint32_t num_bits;

    qcc_kmb_reg_t *regs;
    uint32_t       num_regs;
    uint32_t       cap_regs;

    uint8_t  *op_types;  size_t num_ops, cap_ops;
    uint32_t *q3_0, *q3_1, *q3_2; size_t num_3q, cap_3q;
    uint32_t *q2_0, *q2_1;        size_t num_2q, cap_2q;
    uint32_t *q1_0;               size_t num_1q, cap_1q;
    uint32_t *b1_0;               size_t num_1b, cap_1b;
    uint32_t *bcond;              size_t num_cond, cap_cond;
    qcc_angle128_t *angles;       size_t num_angles, cap_angles;

    uint64_t *info_tag_end_offsets; size_t num_info_tags, cap_info_tags;
    char     *info_tag_chars;       size_t num_info_tag_chars, cap_info_tag_chars;

    /* Exact online counters (kept in sync by emitters and recomputed by reader) */
    uint64_t count_toffoli; /* CCX + CCZ + CSWAP (+ _IF) */
    uint64_t count_t_gate;  /* Clifford hierarchy level-3 Z_POW (+ _IF) */
    uint64_t count_z_pow;   /* Total Z_POW (+ _IF) */
    uint64_t count_hmr;     /* HMR (+ _IF) */
} qcc_kmb_circuit_t;

typedef struct {
    uint32_t num_qubits;
    uint32_t num_bits;
    uint64_t num_ops;
    uint64_t toffoli;      /* CCX + CCZ + CSWAP (+ _IF) */
    uint64_t t_like;       /* Clifford level-3 Z_POW */
    uint64_t magic;        /* CCX+CCZ+CSWAP + level-3 Z_POW + non-Clifford Z_POW */
    uint64_t z_pow;
    uint64_t non_clifford_z_pow; /* level != 1 && level != 2 (excluding exact Clifford Z/S/I) */
    uint64_t hmr;
    uint64_t cnot;
    uint64_t cz;
    uint64_t swap;
    uint64_t r_count;
    uint64_t reaction_depth;
} qcc_kmb_stats_t;

typedef struct {
    size_t op_idx;
    size_t i3q;
    size_t i2q;
    size_t i1q;
    size_t i1b;
    size_t icond;
    size_t iang;
    size_t itag;
    size_t tag_char_prev;
} qcc_kmb_cursor_t;

typedef struct {
    uint8_t        op;
    uint8_t        n_q;
    uint8_t        has_b;
    uint8_t        has_cond;
    uint8_t        has_angle;
    uint8_t        has_info_tag;
    uint32_t       q[3];
    uint32_t       b;
    uint32_t       cond;
    qcc_angle128_t angle;
    const char    *info_tag;
    size_t         info_tag_len;
} qcc_kmb_op_view_t;

#define QCC_KMB_GROW(ptr, len, cap, type) do { \
    if ((len) >= (cap)) { \
        size_t _nc = (cap) ? (cap) * 2u : 256u; \
        void *_np = realloc((ptr), _nc * sizeof(type)); \
        if (!_np) abort(); \
        (ptr) = (type *)_np; \
        (cap) = _nc; \
    } \
} while (0)

static inline void qcc_kmb_init(qcc_kmb_circuit_t *c) {
    memset(c, 0, sizeof(*c));
}

static inline void qcc_kmb_free(qcc_kmb_circuit_t *c) {
    if (!c) return;
    for (uint32_t i = 0u; i < c->num_regs; ++i) free(c->regs[i].targets);
    free(c->regs);
    free(c->op_types);
    free(c->q3_0);
    free(c->q3_1);
    free(c->q3_2);
    free(c->q2_0);
    free(c->q2_1);
    free(c->q1_0);
    free(c->b1_0);
    free(c->bcond);
    free(c->angles);
    free(c->info_tag_end_offsets);
    free(c->info_tag_chars);
    memset(c, 0, sizeof(*c));
}

static inline void qcc_kmb_clear(qcc_kmb_circuit_t *c) {
    if (!c) return;
    for (uint32_t i = 0u; i < c->num_regs; ++i) free(c->regs[i].targets);
    c->num_qubits          = 0u;
    c->num_bits            = 0u;
    c->num_regs            = 0u;
    c->num_ops             = 0u;
    c->num_3q              = 0u;
    c->num_2q              = 0u;
    c->num_1q              = 0u;
    c->num_1b              = 0u;
    c->num_cond            = 0u;
    c->num_angles          = 0u;
    c->num_info_tags       = 0u;
    c->num_info_tag_chars  = 0u;
    c->count_toffoli       = 0u;
    c->count_t_gate        = 0u;
    c->count_z_pow         = 0u;
    c->count_hmr           = 0u;
}

static inline bool qcc_kmb_clone(qcc_kmb_circuit_t *dst, const qcc_kmb_circuit_t *src) {
    qcc_kmb_init(dst);
    dst->num_qubits          = src->num_qubits;
    dst->num_bits            = src->num_bits;
    dst->num_regs            = src->num_regs;
    dst->cap_regs            = src->num_regs;
    dst->num_ops             = src->num_ops;
    dst->cap_ops             = src->num_ops;
    dst->num_3q              = src->num_3q;
    dst->cap_3q              = src->num_3q;
    dst->num_2q              = src->num_2q;
    dst->cap_2q              = src->num_2q;
    dst->num_1q              = src->num_1q;
    dst->cap_1q              = src->num_1q;
    dst->num_1b              = src->num_1b;
    dst->cap_1b              = src->num_1b;
    dst->num_cond            = src->num_cond;
    dst->cap_cond            = src->num_cond;
    dst->num_angles          = src->num_angles;
    dst->cap_angles          = src->num_angles;
    dst->num_info_tags       = src->num_info_tags;
    dst->cap_info_tags       = src->num_info_tags;
    dst->num_info_tag_chars  = src->num_info_tag_chars;
    dst->cap_info_tag_chars  = src->num_info_tag_chars;
    dst->count_toffoli       = src->count_toffoli;
    dst->count_t_gate        = src->count_t_gate;
    dst->count_z_pow         = src->count_z_pow;
    dst->count_hmr           = src->count_hmr;

    if (src->num_regs) {
        dst->regs = (qcc_kmb_reg_t *)calloc(src->num_regs, sizeof(qcc_kmb_reg_t));
        if (!dst->regs) goto fail;
        for (uint32_t i = 0u; i < src->num_regs; ++i) {
            dst->regs[i] = src->regs[i];
            dst->regs[i].cap = src->regs[i].len;
            if (src->regs[i].len) {
                dst->regs[i].targets = (uint32_t *)malloc(src->regs[i].len * sizeof(uint32_t));
                if (!dst->regs[i].targets) goto fail;
                memcpy(dst->regs[i].targets, src->regs[i].targets, src->regs[i].len * sizeof(uint32_t));
            } else {
                dst->regs[i].targets = NULL;
            }
        }
    }
    if (src->num_ops) {
        dst->op_types = (uint8_t *)malloc(src->num_ops);
        if (!dst->op_types) goto fail;
        memcpy(dst->op_types, src->op_types, src->num_ops);
    }
    if (src->num_3q) {
        dst->q3_0 = (uint32_t *)malloc(src->num_3q * sizeof(uint32_t));
        dst->q3_1 = (uint32_t *)malloc(src->num_3q * sizeof(uint32_t));
        dst->q3_2 = (uint32_t *)malloc(src->num_3q * sizeof(uint32_t));
        if (!dst->q3_0 || !dst->q3_1 || !dst->q3_2) goto fail;
        memcpy(dst->q3_0, src->q3_0, src->num_3q * sizeof(uint32_t));
        memcpy(dst->q3_1, src->q3_1, src->num_3q * sizeof(uint32_t));
        memcpy(dst->q3_2, src->q3_2, src->num_3q * sizeof(uint32_t));
    }
    if (src->num_2q) {
        dst->q2_0 = (uint32_t *)malloc(src->num_2q * sizeof(uint32_t));
        dst->q2_1 = (uint32_t *)malloc(src->num_2q * sizeof(uint32_t));
        if (!dst->q2_0 || !dst->q2_1) goto fail;
        memcpy(dst->q2_0, src->q2_0, src->num_2q * sizeof(uint32_t));
        memcpy(dst->q2_1, src->q2_1, src->num_2q * sizeof(uint32_t));
    }
    if (src->num_1q) {
        dst->q1_0 = (uint32_t *)malloc(src->num_1q * sizeof(uint32_t));
        if (!dst->q1_0) goto fail;
        memcpy(dst->q1_0, src->q1_0, src->num_1q * sizeof(uint32_t));
    }
    if (src->num_1b) {
        dst->b1_0 = (uint32_t *)malloc(src->num_1b * sizeof(uint32_t));
        if (!dst->b1_0) goto fail;
        memcpy(dst->b1_0, src->b1_0, src->num_1b * sizeof(uint32_t));
    }
    if (src->num_cond) {
        dst->bcond = (uint32_t *)malloc(src->num_cond * sizeof(uint32_t));
        if (!dst->bcond) goto fail;
        memcpy(dst->bcond, src->bcond, src->num_cond * sizeof(uint32_t));
    }
    if (src->num_angles) {
        dst->angles = (qcc_angle128_t *)malloc(src->num_angles * sizeof(qcc_angle128_t));
        if (!dst->angles) goto fail;
        memcpy(dst->angles, src->angles, src->num_angles * sizeof(qcc_angle128_t));
    }
    if (src->num_info_tags) {
        dst->info_tag_end_offsets = (uint64_t *)malloc(src->num_info_tags * sizeof(uint64_t));
        if (!dst->info_tag_end_offsets) goto fail;
        memcpy(dst->info_tag_end_offsets, src->info_tag_end_offsets, src->num_info_tags * sizeof(uint64_t));
    }
    if (src->num_info_tag_chars) {
        dst->info_tag_chars = (char *)malloc(src->num_info_tag_chars);
        if (!dst->info_tag_chars) goto fail;
        memcpy(dst->info_tag_chars, src->info_tag_chars, src->num_info_tag_chars);
    }
    return true;
fail:
    qcc_kmb_free(dst);
    return false;
}

static inline void qcc_kmb_touch_q(qcc_kmb_circuit_t *c, uint32_t q) {
    if (q + 1u > c->num_qubits) c->num_qubits = q + 1u;
}

static inline void qcc_kmb_touch_b(qcc_kmb_circuit_t *c, uint32_t bit) {
    if (bit + 1u > c->num_bits) c->num_bits = bit + 1u;
}

static inline uint32_t qcc_kmb_add_quantum_reg(
    qcc_kmb_circuit_t *c, const char *name, const uint32_t *qubits, uint32_t n
) {
    if (c->num_regs >= c->cap_regs) {
        uint32_t nc = c->cap_regs ? c->cap_regs * 2u : 16u;
        c->regs = (qcc_kmb_reg_t *)realloc(c->regs, (size_t)nc * sizeof(qcc_kmb_reg_t));
        if (!c->regs) abort();
        c->cap_regs = nc;
    }
    uint32_t idx = c->num_regs++;
    qcc_kmb_reg_t *r = &c->regs[idx];
    memset(r, 0, sizeof(*r));
    size_t sl = name ? strlen(name) : 0u;
    if (sl > 63u) sl = 63u;
    if (sl) memcpy(r->name, name, sl);
    r->name_len = (uint32_t)sl;
    r->targets  = n ? (uint32_t *)malloc((size_t)n * sizeof(uint32_t)) : NULL;
    r->len      = n;
    r->cap      = n;
    for (uint32_t i = 0u; i < n; ++i) {
        qcc_kmb_touch_q(c, qubits[i]);
        r->targets[i] = UINT32_C(0x80000000) | qubits[i];
    }
    return idx;
}

static inline uint32_t qcc_kmb_add_classical_reg(
    qcc_kmb_circuit_t *c, const char *name, const uint32_t *bits, uint32_t n
) {
    if (c->num_regs >= c->cap_regs) {
        uint32_t nc = c->cap_regs ? c->cap_regs * 2u : 16u;
        c->regs = (qcc_kmb_reg_t *)realloc(c->regs, (size_t)nc * sizeof(qcc_kmb_reg_t));
        if (!c->regs) abort();
        c->cap_regs = nc;
    }
    uint32_t idx = c->num_regs++;
    qcc_kmb_reg_t *r = &c->regs[idx];
    memset(r, 0, sizeof(*r));
    size_t sl = name ? strlen(name) : 0u;
    if (sl > 63u) sl = 63u;
    if (sl) memcpy(r->name, name, sl);
    r->name_len = (uint32_t)sl;
    r->targets  = n ? (uint32_t *)malloc((size_t)n * sizeof(uint32_t)) : NULL;
    r->len      = n;
    r->cap      = n;
    for (uint32_t i = 0u; i < n; ++i) {
        qcc_kmb_touch_b(c, bits[i]);
        r->targets[i] = bits[i] & UINT32_C(0x7FFFFFFF);
    }
    return idx;
}

static inline void qcc_kmb_push_op_raw(qcc_kmb_circuit_t *c, uint8_t op) {
    QCC_KMB_GROW(c->op_types, c->num_ops, c->cap_ops, uint8_t);
    c->op_types[c->num_ops++] = op;
}

static inline void qcc_kmb_push_3q(qcc_kmb_circuit_t *c, uint32_t q0, uint32_t q1, uint32_t q2) {
    qcc_kmb_touch_q(c, q0);
    qcc_kmb_touch_q(c, q1);
    qcc_kmb_touch_q(c, q2);
    if (c->num_3q >= c->cap_3q) {
        size_t nc = c->cap_3q ? c->cap_3q * 2u : 256u;
        c->q3_0 = (uint32_t *)realloc(c->q3_0, nc * sizeof(uint32_t));
        c->q3_1 = (uint32_t *)realloc(c->q3_1, nc * sizeof(uint32_t));
        c->q3_2 = (uint32_t *)realloc(c->q3_2, nc * sizeof(uint32_t));
        if (!c->q3_0 || !c->q3_1 || !c->q3_2) abort();
        c->cap_3q = nc;
    }
    c->q3_0[c->num_3q] = q0;
    c->q3_1[c->num_3q] = q1;
    c->q3_2[c->num_3q] = q2;
    c->num_3q++;
}

static inline void qcc_kmb_push_2q(qcc_kmb_circuit_t *c, uint32_t q0, uint32_t q1) {
    qcc_kmb_touch_q(c, q0);
    qcc_kmb_touch_q(c, q1);
    if (c->num_2q >= c->cap_2q) {
        size_t nc = c->cap_2q ? c->cap_2q * 2u : 256u;
        c->q2_0 = (uint32_t *)realloc(c->q2_0, nc * sizeof(uint32_t));
        c->q2_1 = (uint32_t *)realloc(c->q2_1, nc * sizeof(uint32_t));
        if (!c->q2_0 || !c->q2_1) abort();
        c->cap_2q = nc;
    }
    c->q2_0[c->num_2q] = q0;
    c->q2_1[c->num_2q] = q1;
    c->num_2q++;
}

static inline void qcc_kmb_push_1q(qcc_kmb_circuit_t *c, uint32_t q0) {
    qcc_kmb_touch_q(c, q0);
    QCC_KMB_GROW(c->q1_0, c->num_1q, c->cap_1q, uint32_t);
    c->q1_0[c->num_1q++] = q0;
}

static inline void qcc_kmb_push_1b(qcc_kmb_circuit_t *c, uint32_t b0) {
    qcc_kmb_touch_b(c, b0);
    QCC_KMB_GROW(c->b1_0, c->num_1b, c->cap_1b, uint32_t);
    c->b1_0[c->num_1b++] = b0;
}

static inline void qcc_kmb_push_cond(qcc_kmb_circuit_t *c, uint32_t cond) {
    qcc_kmb_touch_b(c, cond);
    QCC_KMB_GROW(c->bcond, c->num_cond, c->cap_cond, uint32_t);
    c->bcond[c->num_cond++] = cond;
}

static inline void qcc_kmb_push_angle(qcc_kmb_circuit_t *c, uint64_t lo, uint64_t hi) {
    QCC_KMB_GROW(c->angles, c->num_angles, c->cap_angles, qcc_angle128_t);
    c->angles[c->num_angles].lo = lo;
    c->angles[c->num_angles].hi = hi;
    c->num_angles++;
}

static inline void qcc_kmb_push_info_tag(qcc_kmb_circuit_t *c, const char *tag, size_t len) {
    QCC_KMB_GROW(c->info_tag_end_offsets, c->num_info_tags, c->cap_info_tags, uint64_t);
    while (c->num_info_tag_chars + len > c->cap_info_tag_chars) {
        size_t nc = c->cap_info_tag_chars ? c->cap_info_tag_chars * 2u : 256u;
        char *np = (char *)realloc(c->info_tag_chars, nc);
        if (!np) abort();
        c->info_tag_chars = np;
        c->cap_info_tag_chars = nc;
    }
    if (len > 0u) {
        memcpy(c->info_tag_chars + c->num_info_tag_chars, tag, len);
        c->num_info_tag_chars += len;
    }
    c->info_tag_end_offsets[c->num_info_tags++] = (uint64_t)c->num_info_tag_chars;
}

/* Unified Low-Level View Emitter */
static inline void qcc_kmb_emit_view(qcc_kmb_circuit_t *c, const qcc_kmb_op_view_t *v) {
    uint8_t op = v->op;
    uint8_t qmask = op & QCC_OP_FLAG_3Q;
    qcc_kmb_push_op_raw(c, op);

    if (qmask == QCC_OP_FLAG_3Q) {
        qcc_kmb_push_3q(c, v->q[0], v->q[1], v->q[2]);
        uint8_t base = op & 0x0Fu;
        if (base == 0u || base == 1u || base == 2u) c->count_toffoli++;
    } else if (qmask == QCC_OP_FLAG_2Q) {
        qcc_kmb_push_2q(c, v->q[0], v->q[1]);
    } else if (qmask == QCC_OP_FLAG_1Q) {
        qcc_kmb_push_1q(c, v->q[0]);
    }

    if (op & QCC_OP_FLAG_1B) {
        qcc_kmb_push_1b(c, v->b);
        if ((op & ~QCC_OP_FLAG_COND) == QCC_OP_HMR) c->count_hmr++;
    }
    if (op & QCC_OP_FLAG_COND) {
        qcc_kmb_push_cond(c, v->cond);
    }
    if (op == QCC_OP_Z_POW || op == QCC_OP_Z_POW_IF) {
        qcc_kmb_push_angle(c, v->angle.lo, v->angle.hi);
        c->count_z_pow++;
        flow_u128 u = qcc_u128_from_angle(v->angle);
        if (qcc_u128_clifford_level(u) == 3u) c->count_t_gate++;
    }
    if (op == QCC_OP_INFO_TAG) {
        qcc_kmb_push_info_tag(c, v->info_tag, v->info_tag_len);
    }
}

/* Convenience Gate Emitters */
static inline void qcc_kmb_neg(qcc_kmb_circuit_t *c) { qcc_kmb_push_op_raw(c, QCC_OP_NEG); }
static inline void qcc_kmb_neg_if(qcc_kmb_circuit_t *c, uint32_t cond) {
    qcc_kmb_push_op_raw(c, QCC_OP_NEG_IF);
    qcc_kmb_push_cond(c, cond);
}
static inline void qcc_kmb_push_condition(qcc_kmb_circuit_t *c, uint32_t cond) {
    qcc_kmb_push_op_raw(c, QCC_OP_PUSH_CONDITION);
    qcc_kmb_push_cond(c, cond);
}
static inline void qcc_kmb_pop_condition(qcc_kmb_circuit_t *c) {
    qcc_kmb_push_op_raw(c, QCC_OP_POP_CONDITION);
}
static inline void qcc_kmb_info_tag(qcc_kmb_circuit_t *c, const char *s) {
    qcc_kmb_push_op_raw(c, QCC_OP_INFO_TAG);
    qcc_kmb_push_info_tag(c, s ? s : "", s ? strlen(s) : 0u);
}

static inline void qcc_kmb_bit_store0(qcc_kmb_circuit_t *c, uint32_t b) {
    qcc_kmb_push_op_raw(c, QCC_OP_BIT_STORE0);
    qcc_kmb_push_1b(c, b);
}
static inline void qcc_kmb_bit_store1(qcc_kmb_circuit_t *c, uint32_t b) {
    qcc_kmb_push_op_raw(c, QCC_OP_BIT_STORE1);
    qcc_kmb_push_1b(c, b);
}
static inline void qcc_kmb_bit_invert(qcc_kmb_circuit_t *c, uint32_t b) {
    qcc_kmb_push_op_raw(c, QCC_OP_BIT_INVERT);
    qcc_kmb_push_1b(c, b);
}

static inline void qcc_kmb_x(qcc_kmb_circuit_t *c, uint32_t q) {
    qcc_kmb_push_op_raw(c, QCC_OP_X);
    qcc_kmb_push_1q(c, q);
}
static inline void qcc_kmb_z(qcc_kmb_circuit_t *c, uint32_t q) {
    qcc_kmb_push_op_raw(c, QCC_OP_Z);
    qcc_kmb_push_1q(c, q);
}
static inline void qcc_kmb_r(qcc_kmb_circuit_t *c, uint32_t q) {
    qcc_kmb_push_op_raw(c, QCC_OP_R);
    qcc_kmb_push_1q(c, q);
}
static inline void qcc_kmb_x_if(qcc_kmb_circuit_t *c, uint32_t q, uint32_t cond) {
    qcc_kmb_push_op_raw(c, QCC_OP_X_IF);
    qcc_kmb_push_1q(c, q);
    qcc_kmb_push_cond(c, cond);
}
static inline void qcc_kmb_z_if(qcc_kmb_circuit_t *c, uint32_t q, uint32_t cond) {
    qcc_kmb_push_op_raw(c, QCC_OP_Z_IF);
    qcc_kmb_push_1q(c, q);
    qcc_kmb_push_cond(c, cond);
}
static inline void qcc_kmb_r_if(qcc_kmb_circuit_t *c, uint32_t q, uint32_t cond) {
    qcc_kmb_push_op_raw(c, QCC_OP_R_IF);
    qcc_kmb_push_1q(c, q);
    qcc_kmb_push_cond(c, cond);
}
static inline void qcc_kmb_hmr(qcc_kmb_circuit_t *c, uint32_t q, uint32_t out_bit) {
    qcc_kmb_push_op_raw(c, QCC_OP_HMR);
    qcc_kmb_push_1q(c, q);
    qcc_kmb_push_1b(c, out_bit);
    c->count_hmr++;
}
static inline void qcc_kmb_hmr_if(qcc_kmb_circuit_t *c, uint32_t q, uint32_t out_bit, uint32_t cond) {
    qcc_kmb_push_op_raw(c, QCC_OP_HMR_IF);
    qcc_kmb_push_1q(c, q);
    qcc_kmb_push_1b(c, out_bit);
    qcc_kmb_push_cond(c, cond);
    c->count_hmr++;
}
static inline void qcc_kmb_cx(qcc_kmb_circuit_t *c, uint32_t ctrl, uint32_t tgt) {
    qcc_kmb_push_op_raw(c, QCC_OP_CX);
    qcc_kmb_push_2q(c, ctrl, tgt);
}
static inline void qcc_kmb_cz(qcc_kmb_circuit_t *c, uint32_t q0, uint32_t q1) {
    qcc_kmb_push_op_raw(c, QCC_OP_CZ);
    qcc_kmb_push_2q(c, q0, q1);
}
static inline void qcc_kmb_swap(qcc_kmb_circuit_t *c, uint32_t q0, uint32_t q1) {
    qcc_kmb_push_op_raw(c, QCC_OP_SWAP);
    qcc_kmb_push_2q(c, q0, q1);
}
static inline void qcc_kmb_cx_if(qcc_kmb_circuit_t *c, uint32_t ctrl, uint32_t tgt, uint32_t cond) {
    qcc_kmb_push_op_raw(c, QCC_OP_CX_IF);
    qcc_kmb_push_2q(c, ctrl, tgt);
    qcc_kmb_push_cond(c, cond);
}
static inline void qcc_kmb_cz_if(qcc_kmb_circuit_t *c, uint32_t q0, uint32_t q1, uint32_t cond) {
    qcc_kmb_push_op_raw(c, QCC_OP_CZ_IF);
    qcc_kmb_push_2q(c, q0, q1);
    qcc_kmb_push_cond(c, cond);
}
static inline void qcc_kmb_swap_if(qcc_kmb_circuit_t *c, uint32_t q0, uint32_t q1, uint32_t cond) {
    qcc_kmb_push_op_raw(c, QCC_OP_SWAP_IF);
    qcc_kmb_push_2q(c, q0, q1);
    qcc_kmb_push_cond(c, cond);
}
static inline void qcc_kmb_ccx(qcc_kmb_circuit_t *c, uint32_t c0, uint32_t c1, uint32_t tgt) {
    qcc_kmb_push_op_raw(c, QCC_OP_CCX);
    qcc_kmb_push_3q(c, c0, c1, tgt);
    c->count_toffoli++;
}
static inline void qcc_kmb_ccz(qcc_kmb_circuit_t *c, uint32_t q0, uint32_t q1, uint32_t q2) {
    qcc_kmb_push_op_raw(c, QCC_OP_CCZ);
    qcc_kmb_push_3q(c, q0, q1, q2);
    c->count_toffoli++;
}
static inline void qcc_kmb_cswap(qcc_kmb_circuit_t *c, uint32_t ctrl, uint32_t q0, uint32_t q1) {
    qcc_kmb_push_op_raw(c, QCC_OP_CSWAP);
    qcc_kmb_push_3q(c, ctrl, q0, q1);
    c->count_toffoli++;
}
static inline void qcc_kmb_ccx_if(qcc_kmb_circuit_t *c, uint32_t c0, uint32_t c1, uint32_t tgt, uint32_t cond) {
    qcc_kmb_push_op_raw(c, QCC_OP_CCX_IF);
    qcc_kmb_push_3q(c, c0, c1, tgt);
    qcc_kmb_push_cond(c, cond);
    c->count_toffoli++;
}
static inline void qcc_kmb_ccz_if(qcc_kmb_circuit_t *c, uint32_t q0, uint32_t q1, uint32_t q2, uint32_t cond) {
    qcc_kmb_push_op_raw(c, QCC_OP_CCZ_IF);
    qcc_kmb_push_3q(c, q0, q1, q2);
    qcc_kmb_push_cond(c, cond);
    c->count_toffoli++;
}
static inline void qcc_kmb_cswap_if(qcc_kmb_circuit_t *c, uint32_t ctrl, uint32_t q0, uint32_t q1, uint32_t cond) {
    qcc_kmb_push_op_raw(c, QCC_OP_CSWAP_IF);
    qcc_kmb_push_3q(c, ctrl, q0, q1);
    qcc_kmb_push_cond(c, cond);
    c->count_toffoli++;
}
static inline void qcc_kmb_z_pow(qcc_kmb_circuit_t *c, uint32_t q, uint64_t lo, uint64_t hi) {
    if ((lo | hi) == 0u) return;
    if (lo == 0u && hi == UINT64_C(0x8000000000000000)) {
        qcc_kmb_z(c, q);
        return;
    }
    qcc_kmb_push_op_raw(c, QCC_OP_Z_POW);
    qcc_kmb_push_1q(c, q);
    qcc_kmb_push_angle(c, lo, hi);
    c->count_z_pow++;
    if (lo == 0u && ((hi & UINT64_C(0x3FFFFFFFFFFFFFFF)) == UINT64_C(0x2000000000000000))) {
        c->count_t_gate++;
    }
}
static inline void qcc_kmb_z_pow_raw(qcc_kmb_circuit_t *c, uint32_t q, uint64_t lo, uint64_t hi) {
    qcc_kmb_push_op_raw(c, QCC_OP_Z_POW);
    qcc_kmb_push_1q(c, q);
    qcc_kmb_push_angle(c, lo, hi);
    c->count_z_pow++;
    flow_u128 u = ((flow_u128)hi << 64u) | (flow_u128)lo;
    if (qcc_u128_clifford_level(u) == 3u) c->count_t_gate++;
}
static inline void qcc_kmb_z_pow_if(qcc_kmb_circuit_t *c, uint32_t q, uint64_t lo, uint64_t hi, uint32_t cond) {
    qcc_kmb_push_op_raw(c, QCC_OP_Z_POW_IF);
    qcc_kmb_push_1q(c, q);
    qcc_kmb_push_cond(c, cond);
    qcc_kmb_push_angle(c, lo, hi);
    c->count_z_pow++;
    flow_u128 u = ((flow_u128)hi << 64u) | (flow_u128)lo;
    if (qcc_u128_clifford_level(u) == 3u) c->count_t_gate++;
}
static inline void qcc_kmb_z_pow_f64(qcc_kmb_circuit_t *c, uint32_t q, double full_turns) {
    double f = full_turns - (double)(int64_t)full_turns;
    if (f < 0.0) f += 1.0;
    if (f == 0.0) return;
    double scaled = f * 18446744073709551616.0;
    uint64_t hi = (uint64_t)scaled;
    double frac = scaled - (double)hi;
    uint64_t lo = (uint64_t)(frac * 18446744073709551616.0);
    qcc_kmb_z_pow(c, q, lo, hi);
}
static inline void qcc_kmb_z_pow_u128(qcc_kmb_circuit_t *c, uint32_t q, flow_u128 full_turns) {
    qcc_kmb_z_pow(c, q, (uint64_t)full_turns, (uint64_t)(full_turns >> 64u));
}

/* Gidney Temporary-AND Compute + HMR Uncomputation */
static inline void qcc_kmb_and_compute(qcc_kmb_circuit_t *c, uint32_t c0, uint32_t c1, uint32_t anc) {
    qcc_kmb_r(c, anc);
    qcc_kmb_ccx(c, c0, c1, anc);
}
static inline void qcc_kmb_and_uncompute(qcc_kmb_circuit_t *c, uint32_t c0, uint32_t c1, uint32_t anc, uint32_t scratch_bit) {
    qcc_kmb_hmr(c, anc, scratch_bit);
    qcc_kmb_cz_if(c, c0, c1, scratch_bit);
}

static inline void qcc_kmb_emit_gidney_iadd(
    qcc_kmb_circuit_t *c,
    const uint32_t *lhs,
    const uint32_t *target,
    const uint32_t *clean_anc,
    uint32_t n,
    uint32_t scratch_bit
) {
    if (n == 0u) return;
    if (n == 1u) {
        qcc_kmb_cx(c, lhs[0], target[0]);
        return;
    }
    qcc_kmb_and_compute(c, lhs[0], target[0], clean_anc[0]);
    for (uint32_t i = 1u; i + 1u < n; ++i) {
        qcc_kmb_cx(c, clean_anc[i - 1u], lhs[i]);
        qcc_kmb_cx(c, clean_anc[i - 1u], target[i]);
        qcc_kmb_and_compute(c, lhs[i], target[i], clean_anc[i]);
        qcc_kmb_cx(c, clean_anc[i - 1u], clean_anc[i]);
    }
    qcc_kmb_cx(c, clean_anc[n - 2u], target[n - 1u]);
    qcc_kmb_cx(c, lhs[n - 1u], target[n - 1u]);
    for (uint32_t k = n - 2u; k >= 1u; --k) {
        qcc_kmb_cx(c, clean_anc[k - 1u], clean_anc[k]);
        qcc_kmb_and_uncompute(c, lhs[k], target[k], clean_anc[k], scratch_bit);
        qcc_kmb_cx(c, clean_anc[k - 1u], lhs[k]);
        qcc_kmb_cx(c, lhs[k], target[k]);
    }
    qcc_kmb_and_uncompute(c, lhs[0], target[0], clean_anc[0], scratch_bit);
    qcc_kmb_cx(c, lhs[0], target[0]);
}

/* Unified Cursor Iterator */
static inline void qcc_kmb_cursor_init(qcc_kmb_cursor_t *cur) {
    memset(cur, 0, sizeof(*cur));
}

static inline bool qcc_kmb_cursor_next(const qcc_kmb_circuit_t *c, qcc_kmb_cursor_t *cur, qcc_kmb_op_view_t *out) {
    if (cur->op_idx >= c->num_ops) return false;
    memset(out, 0, sizeof(*out));
    uint8_t op = c->op_types[cur->op_idx++];
    out->op = op;
    uint8_t qmask = op & QCC_OP_FLAG_3Q;
    if (qmask == QCC_OP_FLAG_3Q) {
        if (cur->i3q >= c->num_3q) return false;
        out->n_q  = 3u;
        out->q[0] = c->q3_0[cur->i3q];
        out->q[1] = c->q3_1[cur->i3q];
        out->q[2] = c->q3_2[cur->i3q];
        cur->i3q++;
    } else if (qmask == QCC_OP_FLAG_2Q) {
        if (cur->i2q >= c->num_2q) return false;
        out->n_q  = 2u;
        out->q[0] = c->q2_0[cur->i2q];
        out->q[1] = c->q2_1[cur->i2q];
        cur->i2q++;
    } else if (qmask == QCC_OP_FLAG_1Q) {
        if (cur->i1q >= c->num_1q) return false;
        out->n_q  = 1u;
        out->q[0] = c->q1_0[cur->i1q];
        cur->i1q++;
    }
    if (op & QCC_OP_FLAG_1B) {
        if (cur->i1b >= c->num_1b) return false;
        out->has_b = 1u;
        out->b     = c->b1_0[cur->i1b++];
    }
    if (op & QCC_OP_FLAG_COND) {
        if (cur->icond >= c->num_cond) return false;
        out->has_cond = 1u;
        out->cond     = c->bcond[cur->icond++];
    }
    if (op == QCC_OP_Z_POW || op == QCC_OP_Z_POW_IF) {
        if (cur->iang >= c->num_angles) return false;
        out->has_angle = 1u;
        out->angle     = c->angles[cur->iang++];
    }
    if (op == QCC_OP_INFO_TAG) {
        if (cur->itag >= c->num_info_tags) return false;
        uint64_t end_off = c->info_tag_end_offsets[cur->itag++];
        if (end_off < cur->tag_char_prev || end_off > c->num_info_tag_chars) return false;
        out->has_info_tag = 1u;
        out->info_tag     = c->info_tag_chars + cur->tag_char_prev;
        out->info_tag_len = (size_t)(end_off - cur->tag_char_prev);
        cur->tag_char_prev = (size_t)end_off;
    }
    return true;
}

/* Exact Kickmix-matching Validation & Statistics (including compute_reaction_depth) */
static inline bool qcc_kmb_recompute_stats(const qcc_kmb_circuit_t *c, qcc_kmb_stats_t *out) {
    memset(out, 0, sizeof(*out));
    uint32_t max_q_plus1 = 0u;
    uint32_t max_b_plus1 = 0u;
    for (uint32_t r = 0u; r < c->num_regs; ++r) {
        for (uint32_t k = 0u; k < c->regs[r].len; ++k) {
            uint32_t t = c->regs[r].targets[k];
            if (t & UINT32_C(0x80000000)) {
                uint32_t q = t & UINT32_C(0x7FFFFFFF);
                if (q + 1u > max_q_plus1) max_q_plus1 = q + 1u;
            } else {
                if (t + 1u > max_b_plus1) max_b_plus1 = t + 1u;
            }
        }
    }
    for (size_t i = 0u; i < c->num_3q; ++i) {
        if (c->q3_0[i] + 1u > max_q_plus1) max_q_plus1 = c->q3_0[i] + 1u;
        if (c->q3_1[i] + 1u > max_q_plus1) max_q_plus1 = c->q3_1[i] + 1u;
        if (c->q3_2[i] + 1u > max_q_plus1) max_q_plus1 = c->q3_2[i] + 1u;
    }
    for (size_t i = 0u; i < c->num_2q; ++i) {
        if (c->q2_0[i] + 1u > max_q_plus1) max_q_plus1 = c->q2_0[i] + 1u;
        if (c->q2_1[i] + 1u > max_q_plus1) max_q_plus1 = c->q2_1[i] + 1u;
    }
    for (size_t i = 0u; i < c->num_1q; ++i) {
        if (c->q1_0[i] + 1u > max_q_plus1) max_q_plus1 = c->q1_0[i] + 1u;
    }
    for (size_t i = 0u; i < c->num_1b; ++i) {
        if (c->b1_0[i] + 1u > max_b_plus1) max_b_plus1 = c->b1_0[i] + 1u;
    }
    for (size_t i = 0u; i < c->num_cond; ++i) {
        if (c->bcond[i] + 1u > max_b_plus1) max_b_plus1 = c->bcond[i] + 1u;
    }

    uint32_t nq = c->num_qubits > max_q_plus1 ? c->num_qubits : max_q_plus1;
    uint32_t nb = c->num_bits > max_b_plus1 ? c->num_bits : max_b_plus1;
    out->num_qubits = max_q_plus1;
    out->num_bits   = max_b_plus1;
    out->num_ops    = c->num_ops;

    uint64_t *bit_layer = nb ? (uint64_t *)calloc(nb, sizeof(uint64_t)) : NULL;
    size_t cond_cap = 64u, cond_len = 0u;
    uint64_t *cond_stack = (uint64_t *)calloc(cond_cap, sizeof(uint64_t));
    if ((nb && !bit_layer) || !cond_stack) {
        free(bit_layer);
        free(cond_stack);
        return false;
    }

    uint64_t cur_cond_depth = 0u;
    qcc_kmb_cursor_t cur;
    qcc_kmb_cursor_init(&cur);
    qcc_kmb_op_view_t v;
    while (qcc_kmb_cursor_next(c, &cur, &v)) {
        uint8_t base = v.op & (uint8_t)(~QCC_OP_FLAG_COND);
        if (v.op == QCC_OP_PUSH_CONDITION) {
            if (cond_len >= cond_cap) {
                cond_cap *= 2u;
                cond_stack = (uint64_t *)realloc(cond_stack, cond_cap * sizeof(uint64_t));
                if (!cond_stack) abort();
            }
            cond_stack[cond_len++] = bit_layer[v.cond];
            if (bit_layer[v.cond] > cur_cond_depth) cur_cond_depth = bit_layer[v.cond];
            continue;
        }
        if (v.op == QCC_OP_POP_CONDITION) {
            if (cond_len > 0u) cond_len--;
            cur_cond_depth = 0u;
            for (size_t k = 0u; k < cond_len; ++k) {
                if (cond_stack[k] > cur_cond_depth) cur_cond_depth = cond_stack[k];
            }
            continue;
        }

        uint64_t d = cur_cond_depth;
        if (v.has_cond && bit_layer[v.cond] > d) d = bit_layer[v.cond];

        if (base == QCC_OP_BIT_STORE0 || base == QCC_OP_BIT_STORE1) {
            bit_layer[v.b] = d;
        } else if (base == QCC_OP_BIT_INVERT) {
            if (d > bit_layer[v.b]) bit_layer[v.b] = d;
        } else if (base == QCC_OP_HMR) {
            out->hmr++;
            uint64_t next_d = d + 1u;
            bit_layer[v.b] = next_d;
            if (next_d > out->reaction_depth) out->reaction_depth = next_d;
        } else if (base == QCC_OP_CCX || base == QCC_OP_CCZ || base == QCC_OP_CSWAP) {
            out->toffoli++;
            out->magic++;
        } else if (base == QCC_OP_CX) {
            out->cnot++;
        } else if (base == QCC_OP_CZ) {
            out->cz++;
        } else if (base == QCC_OP_SWAP) {
            out->swap++;
        } else if (base == QCC_OP_R) {
            out->r_count++;
        } else if (base == QCC_OP_Z_POW) {
            out->z_pow++;
            flow_u128 ang = qcc_u128_from_angle(v.angle);
            uint8_t lvl = qcc_u128_clifford_level(ang);
            if (lvl == 3u) {
                out->t_like++;
                out->magic++;
            } else if (lvl > 2u) {
                out->non_clifford_z_pow++;
                out->magic++;
            }
        }
    }

    (void)nq;
    bool ok = (cur.op_idx == c->num_ops &&
               cur.i3q == c->num_3q &&
               cur.i2q == c->num_2q &&
               cur.i1q == c->num_1q &&
               cur.i1b == c->num_1b &&
               cur.icond == c->num_cond &&
               cur.iang == c->num_angles &&
               cur.itag == c->num_info_tags);
    free(bit_layer);
    free(cond_stack);
    return ok;
}

/* Binary .kmb v2 Writer & Reader */
static inline void qcc__write_pkt_hdr(FILE *f, uint64_t pkt_type, uint64_t payload_bytes) {
    fwrite(&pkt_type, sizeof(uint64_t), 1, f);
    fwrite(&payload_bytes, sizeof(uint64_t), 1, f);
}

static inline void qcc__write_u32_pkt(FILE *f, uint64_t pkt_type, const uint32_t *arr, size_t count) {
    qcc__write_pkt_hdr(f, pkt_type, count * sizeof(uint32_t));
    if (count > 0u) fwrite(arr, sizeof(uint32_t), count, f);
}

static inline int qcc_kmb_write(const qcc_kmb_circuit_t *c, const char *path) {
    qcc_kmb_stats_t st;
    if (!qcc_kmb_recompute_stats(c, &st)) return -2;

    FILE *f = fopen(path, "wb");
    if (!f) return -1;

    uint32_t version = 2u;
    uint32_t nq = st.num_qubits;
    uint32_t nb = st.num_bits;
    uint64_t num_ops  = (uint64_t)c->num_ops;
    uint64_t num_3q   = (uint64_t)c->num_3q;
    uint64_t num_2q   = (uint64_t)c->num_2q;
    uint64_t num_1q   = (uint64_t)c->num_1q;
    uint64_t num_1b   = (uint64_t)c->num_1b;
    uint64_t num_cond = (uint64_t)c->num_cond;
    uint64_t num_ang  = (uint64_t)c->num_angles;
    uint64_t num_tags = (uint64_t)c->num_info_tags;
    uint64_t num_tchr = (uint64_t)c->num_info_tag_chars;

    fwrite(QCC_KMB_MAGIC, 1, 16, f);
    fwrite(&version, sizeof(uint32_t), 1, f);
    fwrite(&nq, sizeof(uint32_t), 1, f);
    fwrite(&nb, sizeof(uint32_t), 1, f);
    fwrite(&c->num_regs, sizeof(uint32_t), 1, f);
    fwrite(&num_ops, sizeof(uint64_t), 1, f);
    fwrite(&num_3q, sizeof(uint64_t), 1, f);
    fwrite(&num_2q, sizeof(uint64_t), 1, f);
    fwrite(&num_1q, sizeof(uint64_t), 1, f);
    fwrite(&num_1b, sizeof(uint64_t), 1, f);
    fwrite(&num_cond, sizeof(uint64_t), 1, f);
    fwrite(&num_ang, sizeof(uint64_t), 1, f);
    fwrite(&num_tags, sizeof(uint64_t), 1, f);
    fwrite(&num_tchr, sizeof(uint64_t), 1, f);

    for (uint32_t i = 0u; i < c->num_regs; ++i) {
        const qcc_kmb_reg_t *r = &c->regs[i];
        uint64_t pay = (uint64_t)8u + (uint64_t)r->name_len + (uint64_t)r->len * sizeof(uint32_t);
        qcc__write_pkt_hdr(f, QCC_PKT_REGISTER_TARGETS, pay);
        fwrite(&r->name_len, sizeof(uint32_t), 1, f);
        fwrite(&r->len, sizeof(uint32_t), 1, f);
        if (r->name_len > 0u) fwrite(r->name, 1, r->name_len, f);
        if (r->len > 0u) fwrite(r->targets, sizeof(uint32_t), r->len, f);
    }

    qcc__write_pkt_hdr(f, QCC_PKT_OPERATION_TYPES, num_ops);
    if (num_ops > 0u) fwrite(c->op_types, 1, num_ops, f);

    qcc__write_u32_pkt(f, QCC_PKT_Q0_3Q, c->q3_0, c->num_3q);
    qcc__write_u32_pkt(f, QCC_PKT_Q1_3Q, c->q3_1, c->num_3q);
    qcc__write_u32_pkt(f, QCC_PKT_Q2_3Q, c->q3_2, c->num_3q);
    qcc__write_u32_pkt(f, QCC_PKT_Q0_2Q, c->q2_0, c->num_2q);
    qcc__write_u32_pkt(f, QCC_PKT_Q1_2Q, c->q2_1, c->num_2q);
    qcc__write_u32_pkt(f, QCC_PKT_Q0_1Q, c->q1_0, c->num_1q);
    qcc__write_u32_pkt(f, QCC_PKT_B0_1B, c->b1_0, c->num_1b);
    qcc__write_u32_pkt(f, QCC_PKT_B0_COND, c->bcond, c->num_cond);

    qcc__write_pkt_hdr(f, QCC_PKT_ANGLE, num_ang * sizeof(qcc_angle128_t));
    if (num_ang > 0u) fwrite(c->angles, sizeof(qcc_angle128_t), c->num_angles, f);

    qcc__write_pkt_hdr(f, QCC_PKT_INFO_TAG_DATA, (num_tags + 1u) * sizeof(uint64_t) + num_tchr);
    fwrite(&num_tags, sizeof(uint64_t), 1, f);
    if (num_tags > 0u) fwrite(c->info_tag_end_offsets, sizeof(uint64_t), c->num_info_tags, f);
    if (num_tchr > 0u) fwrite(c->info_tag_chars, 1, c->num_info_tag_chars, f);

    int rc = ferror(f) ? -3 : 0;
    fclose(f);
    return rc;
}

static inline bool qcc__read_exact(FILE *f, void *dst, size_t n) {
    return n == 0u || fread(dst, 1, n, f) == n;
}

static inline bool qcc__read_pkt_expect(FILE *f, uint64_t expected_type, uint64_t expected_bytes) {
    uint64_t pt = 0u, pb = 0u;
    if (!qcc__read_exact(f, &pt, sizeof(pt)) || !qcc__read_exact(f, &pb, sizeof(pb))) return false;
    return pt == expected_type && pb == expected_bytes;
}

static inline bool qcc__read_u32_pkt(FILE *f, uint64_t pkt_type, uint32_t **out_arr, size_t count) {
    if (!qcc__read_pkt_expect(f, pkt_type, count * sizeof(uint32_t))) return false;
    if (count == 0u) {
        *out_arr = NULL;
        return true;
    }
    *out_arr = (uint32_t *)malloc(count * sizeof(uint32_t));
    if (!*out_arr) return false;
    return qcc__read_exact(f, *out_arr, count * sizeof(uint32_t));
}

static inline int qcc_kmb_read(qcc_kmb_circuit_t *out, const char *path) {
    qcc_kmb_init(out);
    FILE *f = fopen(path, "rb");
    if (!f) return -1;

    uint8_t magic[16];
    uint32_t version = 0u, nq = 0u, nb = 0u, nregs = 0u;
    uint64_t nops = 0u, n3q = 0u, n2q = 0u, n1q = 0u, n1b = 0u, ncond = 0u, nang = 0u, ntags = 0u, ntchr = 0u;

    if (!qcc__read_exact(f, magic, 16u) || memcmp(magic, QCC_KMB_MAGIC, 16u) != 0) goto fail;
    if (!qcc__read_exact(f, &version, 4u) || version != 2u) goto fail;
    if (!qcc__read_exact(f, &nq, 4u) || !qcc__read_exact(f, &nb, 4u) || !qcc__read_exact(f, &nregs, 4u)) goto fail;
    if (!qcc__read_exact(f, &nops, 8u) || !qcc__read_exact(f, &n3q, 8u) || !qcc__read_exact(f, &n2q, 8u) ||
        !qcc__read_exact(f, &n1q, 8u) || !qcc__read_exact(f, &n1b, 8u) || !qcc__read_exact(f, &ncond, 8u) ||
        !qcc__read_exact(f, &nang, 8u) || !qcc__read_exact(f, &ntags, 8u) || !qcc__read_exact(f, &ntchr, 8u)) {
        goto fail;
    }

    out->num_qubits          = nq;
    out->num_bits            = nb;
    out->num_regs            = nregs;
    out->cap_regs            = nregs;
    out->num_ops             = (size_t)nops;
    out->cap_ops             = (size_t)nops;
    out->num_3q              = (size_t)n3q;
    out->cap_3q              = (size_t)n3q;
    out->num_2q              = (size_t)n2q;
    out->cap_2q              = (size_t)n2q;
    out->num_1q              = (size_t)n1q;
    out->cap_1q              = (size_t)n1q;
    out->num_1b              = (size_t)n1b;
    out->cap_1b              = (size_t)n1b;
    out->num_cond            = (size_t)ncond;
    out->cap_cond            = (size_t)ncond;
    out->num_angles          = (size_t)nang;
    out->cap_angles          = (size_t)nang;
    out->num_info_tags       = (size_t)ntags;
    out->cap_info_tags       = (size_t)ntags;
    out->num_info_tag_chars  = (size_t)ntchr;
    out->cap_info_tag_chars  = (size_t)ntchr;

    if (nregs > 0u) {
        out->regs = (qcc_kmb_reg_t *)calloc(nregs, sizeof(qcc_kmb_reg_t));
        if (!out->regs) goto fail;
        for (uint32_t i = 0u; i < nregs; ++i) {
            uint64_t pt = 0u, pb = 0u;
            uint32_t name_len = 0u, rlen = 0u;
            if (!qcc__read_exact(f, &pt, 8u) || !qcc__read_exact(f, &pb, 8u) || pt != QCC_PKT_REGISTER_TARGETS) goto fail;
            if (!qcc__read_exact(f, &name_len, 4u) || !qcc__read_exact(f, &rlen, 4u)) goto fail;
            if (8u + (uint64_t)name_len + (uint64_t)rlen * 4u != pb) goto fail;

            qcc_kmb_reg_t *r = &out->regs[i];
            if (name_len > 0u) {
                char *tmp = (char *)malloc(name_len);
                if (!tmp || !qcc__read_exact(f, tmp, name_len)) {
                    free(tmp);
                    goto fail;
                }
                size_t copy_n = name_len > 63u ? 63u : name_len;
                memcpy(r->name, tmp, copy_n);
                r->name_len = (uint32_t)copy_n;
                free(tmp);
            }
            r->len = rlen;
            r->cap = rlen;
            if (rlen > 0u) {
                r->targets = (uint32_t *)malloc((size_t)rlen * sizeof(uint32_t));
                if (!r->targets || !qcc__read_exact(f, r->targets, (size_t)rlen * 4u)) goto fail;
            }
        }
    }

    if (!qcc__read_pkt_expect(f, QCC_PKT_OPERATION_TYPES, nops)) goto fail;
    if (nops > 0u) {
        out->op_types = (uint8_t *)malloc((size_t)nops);
        if (!out->op_types || !qcc__read_exact(f, out->op_types, (size_t)nops)) goto fail;
    }

    if (!qcc__read_u32_pkt(f, QCC_PKT_Q0_3Q, &out->q3_0, out->num_3q)) goto fail;
    if (!qcc__read_u32_pkt(f, QCC_PKT_Q1_3Q, &out->q3_1, out->num_3q)) goto fail;
    if (!qcc__read_u32_pkt(f, QCC_PKT_Q2_3Q, &out->q3_2, out->num_3q)) goto fail;
    if (!qcc__read_u32_pkt(f, QCC_PKT_Q0_2Q, &out->q2_0, out->num_2q)) goto fail;
    if (!qcc__read_u32_pkt(f, QCC_PKT_Q1_2Q, &out->q2_1, out->num_2q)) goto fail;
    if (!qcc__read_u32_pkt(f, QCC_PKT_Q0_1Q, &out->q1_0, out->num_1q)) goto fail;
    if (!qcc__read_u32_pkt(f, QCC_PKT_B0_1B, &out->b1_0, out->num_1b)) goto fail;
    if (!qcc__read_u32_pkt(f, QCC_PKT_B0_COND, &out->bcond, out->num_cond)) goto fail;

    if (!qcc__read_pkt_expect(f, QCC_PKT_ANGLE, nang * sizeof(qcc_angle128_t))) goto fail;
    if (nang > 0u) {
        out->angles = (qcc_angle128_t *)malloc(out->num_angles * sizeof(qcc_angle128_t));
        if (!out->angles || !qcc__read_exact(f, out->angles, out->num_angles * sizeof(qcc_angle128_t))) goto fail;
    }

    if (!qcc__read_pkt_expect(f, QCC_PKT_INFO_TAG_DATA, (ntags + 1u) * sizeof(uint64_t) + ntchr)) goto fail;
    {
        uint64_t tag_cnt2 = 0u;
        if (!qcc__read_exact(f, &tag_cnt2, sizeof(uint64_t)) || tag_cnt2 != ntags) goto fail;
        if (ntags > 0u) {
            out->info_tag_end_offsets = (uint64_t *)malloc(out->num_info_tags * sizeof(uint64_t));
            if (!out->info_tag_end_offsets ||
                !qcc__read_exact(f, out->info_tag_end_offsets, out->num_info_tags * sizeof(uint64_t))) {
                goto fail;
            }
        }
        if (ntchr > 0u) {
            out->info_tag_chars = (char *)malloc(out->num_info_tag_chars);
            if (!out->info_tag_chars || !qcc__read_exact(f, out->info_tag_chars, out->num_info_tag_chars)) goto fail;
        }
    }

    if (fgetc(f) != EOF) goto fail;
    fclose(f);

    qcc_kmb_stats_t st;
    if (!qcc_kmb_recompute_stats(out, &st)) {
        qcc_kmb_free(out);
        return -3;
    }
    out->count_toffoli = st.toffoli;
    out->count_t_gate  = st.t_like;
    out->count_z_pow   = st.z_pow;
    out->count_hmr     = st.hmr;
    return 0;

fail:
    fclose(f);
    qcc_kmb_free(out);
    return -2;
}

#ifdef __cplusplus
}
#endif

#endif /* QCC_KMB_H */
