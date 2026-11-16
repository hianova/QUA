/*
 * qua_kmb.h — Thin Compatibility Adapter over Universal Compiler Core (../qcc/qcc_kmb.h)
 * Target: Google Quantum AI — Verifiable Quantum Circuit Compilation (Kickmix v2)
 * Preserves existing `qua_*` builder types & inline helpers for benchmark synthesizers
 * (`qua_ele.h`, `qua_bitop.h`, `qua_invsqrt.h`, `qua_dqi.h`, `qua_dfthcsos.h`, `qua_synth.c`).
 */
#ifndef QUA_KMB_H
#define QUA_KMB_H

#if defined(__has_include) && __has_include("qcc/qcc_kmb.h")
#include "qcc/qcc_kmb.h"
#else
#include "../../qcc/qcc_kmb.h"
#endif

#define QUA_KMB_MAGIC         QCC_KMB_MAGIC
#define QUA_OP_HAS_1Q         QCC_OP_FLAG_1Q
#define QUA_OP_HAS_2Q         QCC_OP_FLAG_2Q
#define QUA_OP_HAS_3Q         QCC_OP_FLAG_3Q
#define QUA_OP_HAS_1B         QCC_OP_FLAG_1B
#define QUA_OP_HAS_COND       QCC_OP_FLAG_COND

#define QUA_OP_NEG            QCC_OP_NEG
#define QUA_OP_POP_CONDITION  QCC_OP_POP_CONDITION
#define QUA_OP_INFO_TAG       QCC_OP_INFO_TAG
#define QUA_OP_BIT_STORE0     QCC_OP_BIT_STORE0
#define QUA_OP_BIT_STORE1     QCC_OP_BIT_STORE1
#define QUA_OP_BIT_INVERT     QCC_OP_BIT_INVERT
#define QUA_OP_Z              QCC_OP_Z
#define QUA_OP_X              QCC_OP_X
#define QUA_OP_R              QCC_OP_R
#define QUA_OP_Z_POW          QCC_OP_Z_POW
#define QUA_OP_HMR            QCC_OP_HMR
#define QUA_OP_CZ             QCC_OP_CZ
#define QUA_OP_CX             QCC_OP_CX
#define QUA_OP_SWAP           QCC_OP_SWAP
#define QUA_OP_CCZ            QCC_OP_CCZ
#define QUA_OP_CCX            QCC_OP_CCX
#define QUA_OP_CSWAP          QCC_OP_CSWAP
#define QUA_OP_NEG_IF         QCC_OP_NEG_IF
#define QUA_OP_PUSH_CONDITION QCC_OP_PUSH_CONDITION
#define QUA_OP_BIT_STORE0_IF  QCC_OP_BIT_STORE0_IF
#define QUA_OP_BIT_STORE1_IF  QCC_OP_BIT_STORE1_IF
#define QUA_OP_BIT_INVERT_IF  QCC_OP_BIT_INVERT_IF
#define QUA_OP_Z_IF           QCC_OP_Z_IF
#define QUA_OP_X_IF           QCC_OP_X_IF
#define QUA_OP_R_IF           QCC_OP_R_IF
#define QUA_OP_Z_POW_IF       QCC_OP_Z_POW_IF
#define QUA_OP_HMR_IF         QCC_OP_HMR_IF
#define QUA_OP_CZ_IF          QCC_OP_CZ_IF
#define QUA_OP_CX_IF          QCC_OP_CX_IF
#define QUA_OP_SWAP_IF        QCC_OP_SWAP_IF
#define QUA_OP_CCZ_IF         QCC_OP_CCZ_IF
#define QUA_OP_CCX_IF         QCC_OP_CCX_IF
#define QUA_OP_CSWAP_IF       QCC_OP_CSWAP_IF

#define QUA_PKT_REGISTER_TARGETS QCC_PKT_REGISTER_TARGETS
#define QUA_PKT_OPERATION_TYPES  QCC_PKT_OPERATION_TYPES
#define QUA_PKT_Q0_3Q            QCC_PKT_Q0_3Q
#define QUA_PKT_Q1_3Q            QCC_PKT_Q1_3Q
#define QUA_PKT_Q2_3Q            QCC_PKT_Q2_3Q
#define QUA_PKT_Q0_2Q            QCC_PKT_Q0_2Q
#define QUA_PKT_Q1_2Q            QCC_PKT_Q1_2Q
#define QUA_PKT_Q0_1Q            QCC_PKT_Q0_1Q
#define QUA_PKT_B0_1B            QCC_PKT_B0_1B
#define QUA_PKT_B0_COND          QCC_PKT_B0_COND
#define QUA_PKT_ANGLE            QCC_PKT_ANGLE
#define QUA_PKT_INFO_TAG         QCC_PKT_INFO_TAG_DATA

typedef qcc_kmb_op_t     qua_op_t;
typedef qcc_kmb_pkt_t    qua_pkt_type_t;
typedef qcc_angle128_t   qua_angle128_t;
typedef qcc_kmb_reg_t    qua_reg_t;
typedef qcc_kmb_circuit_t qua_builder_t;

#define QUA_GROW(ptr, len, cap, elem_sz) do { \
    if ((len) >= (cap)) { \
        size_t _nc = (cap) ? (cap) * 2u : 256u; \
        void *_np = realloc((ptr), _nc * (elem_sz)); \
        if (!_np) abort(); \
        (ptr) = _np; \
        (cap) = _nc; \
    } \
} while (0)

static inline void qua_builder_init(qua_builder_t *b) { qcc_kmb_init(b); }
static inline void qua_builder_free(qua_builder_t *b) { qcc_kmb_free(b); }
static inline void qua_touch_q(qua_builder_t *b, uint32_t q) { qcc_kmb_touch_q(b, q); }
static inline void qua_touch_b(qua_builder_t *b, uint32_t bit) { qcc_kmb_touch_b(b, bit); }
static inline void qua_push_op(qua_builder_t *b, uint8_t op) { qcc_kmb_push_op_raw(b, op); }

static inline uint32_t qua_add_quantum_reg(qua_builder_t *b, const char *name, const uint32_t *qubits, uint32_t n) {
    return qcc_kmb_add_quantum_reg(b, name, qubits, n);
}
static inline uint32_t qua_add_classical_reg(qua_builder_t *b, const char *name, const uint32_t *bits, uint32_t n) {
    return qcc_kmb_add_classical_reg(b, name, bits, n);
}

static inline void qua_z(qua_builder_t *b, uint32_t q) { qcc_kmb_z(b, q); }
static inline void qua_z_if(qua_builder_t *b, uint32_t q, uint32_t cond_bit) { qcc_kmb_z_if(b, q, cond_bit); }
static inline void qua_x(qua_builder_t *b, uint32_t q) { qcc_kmb_x(b, q); }
static inline void qua_r(qua_builder_t *b, uint32_t q) { qcc_kmb_r(b, q); }
static inline void qua_cx(qua_builder_t *b, uint32_t c, uint32_t t) { qcc_kmb_cx(b, c, t); }
static inline void qua_cz(qua_builder_t *b, uint32_t q0, uint32_t q1) { qcc_kmb_cz(b, q0, q1); }
static inline void qua_cz_if(qua_builder_t *b, uint32_t q0, uint32_t q1, uint32_t cond_bit) { qcc_kmb_cz_if(b, q0, q1, cond_bit); }
static inline void qua_ccx(qua_builder_t *b, uint32_t c0, uint32_t c1, uint32_t t) { qcc_kmb_ccx(b, c0, c1, t); }
static inline void qua_cswap(qua_builder_t *b, uint32_t c, uint32_t x, uint32_t y) { qcc_kmb_cswap(b, c, x, y); }
static inline void qua_hmr(qua_builder_t *b, uint32_t q, uint32_t out_bit) { qcc_kmb_hmr(b, q, out_bit); }
static inline void qua_bit_store0(qua_builder_t *b, uint32_t bit) { qcc_kmb_bit_store0(b, bit); }
static inline void qua_bit_store1_if(qua_builder_t *b, uint32_t bit, uint32_t cond_bit) {
    qcc_kmb_touch_b(b, bit); qcc_kmb_touch_b(b, cond_bit);
    qcc_kmb_push_op_raw(b, QCC_OP_BIT_STORE1_IF);
    qcc_kmb_push_1b(b, bit);
    qcc_kmb_push_cond(b, cond_bit);
}
static inline void qua_bit_invert(qua_builder_t *b, uint32_t bit) { qcc_kmb_bit_invert(b, bit); }
static inline void qua_bit_invert_if(qua_builder_t *b, uint32_t bit, uint32_t cond_bit) {
    qcc_kmb_touch_b(b, bit); qcc_kmb_touch_b(b, cond_bit);
    qcc_kmb_push_op_raw(b, QCC_OP_BIT_INVERT_IF);
    qcc_kmb_push_1b(b, bit);
    qcc_kmb_push_cond(b, cond_bit);
}
static inline void qua_push_condition(qua_builder_t *b, uint32_t cond_bit) { qcc_kmb_push_condition(b, cond_bit); }
static inline void qua_pop_condition(qua_builder_t *b) { qcc_kmb_pop_condition(b); }

static inline void qua_z_pow_turns_u128(qua_builder_t *b, uint32_t q, flow_u128 turns_q128) {
    qcc_kmb_z_pow(b, q, (uint64_t)turns_q128, (uint64_t)(turns_q128 >> 64));
}

static inline void qua_and_compute(qua_builder_t *b, uint32_t c0, uint32_t c1, uint32_t anc) {
    qcc_kmb_and_compute(b, c0, c1, anc);
}
static inline void qua_and_uncompute(qua_builder_t *b, uint32_t c0, uint32_t c1, uint32_t anc, uint32_t scratch_bit) {
    qcc_kmb_and_uncompute(b, c0, c1, anc, scratch_bit);
}

static inline void qua_emit_gidney_iadd(
    qua_builder_t *b,
    const uint32_t *lhs,
    const uint32_t *target,
    const uint32_t *clean_anc,
    uint32_t n,
    uint32_t scratch_bit
) {
    qcc_kmb_emit_gidney_iadd(b, lhs, target, clean_anc, n, scratch_bit);
}

static inline int qua_builder_write_kmb(const qua_builder_t *b, const char *path) {
    return qcc_kmb_write(b, path);
}

#endif /* QUA_KMB_H */
