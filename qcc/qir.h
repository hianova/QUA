/*
 * qir.h — Executable Semantic SSA Intermediate Representation (QIR)
 * Part of QCC (Quantum Circuit Compiler) for FLOW v7.0.1.
 * Pure C17, zero external dependencies, zero indirect function calls.
 */
#ifndef QIR_H
#define QIR_H

#include "qcc_kmb.h"

#ifdef __cplusplus
extern "C" {
#endif

#define QIR_MAX_VAL_BITS 256u
#define QIR_WORDS        4u /* 4 x 64-bit = 256 bits max per SSA value */

typedef enum {
    QIR_REP_STD       = 0, /* Standard binary representation (LSB = bit 0) */
    QIR_REP_CS        = 1, /* Carry-Save redundant pair (outs[0]=sum, outs[1]=carry) */
    QIR_REP_RNS       = 2, /* Residue Number System */
    QIR_REP_LIN_GAUGE = 3, /* Linear gauge M * x (invertible F2 basis) */
    QIR_REP_CLASSICAL = 4  /* Classical bit register residing in .kmb cbits */
} qir_rep_t;

enum {
    QIR_FLAG_APPROX          = 1u << 0, /* Allows bounded truncation / approximate LUT */
    QIR_FLAG_CLEAN_ANC       = 1u << 1, /* Provably clean |0> output / ancilla contract */
    QIR_FLAG_PHASE_ONLY      = 1u << 2, /* Pure diagonal phase operation */
    QIR_FLAG_ZERO_MASK_KNOWN = 1u << 3  /* zero_mask[] proven valid */
};

typedef enum {
    QIR_OP_CONST = 0,
    QIR_OP_SLICE,
    QIR_OP_CONCAT,
    QIR_OP_LINF2,
    QIR_OP_AND_VEC,
    QIR_OP_ADD,
    QIR_OP_SUB,
    QIR_OP_CS_ADD3,
    QIR_OP_MUL_FXP,
    QIR_OP_SQR_FXP,
    QIR_OP_GF2M_ADD,
    QIR_OP_GF2M_MUL,
    QIR_OP_GF2M_INV,
    QIR_OP_LUT_ROM,
    QIR_OP_CLZ_NORM,
    QIR_OP_CTRL_SEL,
    QIR_OP_PHASE_POLY,
    QIR_OP_OPAQUE_KMB,
    QIR_OP_HMR_MEAS,
    QIR_OP_BIT_OP,
    QIR_OP_SHL_VAR,
    QIR_OP_SHR_VAR,
    QIR_OP_EQ_ZERO,
    QIR_OP_PHASE_DOT,
    QIR_OP_COUNT
} qir_op_kind_t;

typedef struct {
    uint32_t  id;
    uint16_t  width;
    qir_rep_t rep;
    uint64_t  zero_mask[QIR_WORDS];  /* bit k = 1 => provably constant 0 */
    uint64_t  const_mask[QIR_WORDS]; /* bit k = 1 => provably constant (0 or 1) */
    uint64_t  const_val[QIR_WORDS];  /* value when const_mask bit k = 1 */
    uint32_t  flags;
    uint32_t  def_op;                /* UINT32_MAX if module input */
    const char *name;                /* Optional debug name (arena-owned) */
} qir_val_t;

typedef struct {
    uint64_t  mask[QIR_WORDS]; /* Bits over flattened concatenated input bits */
    uint8_t   is_xor;          /* 0 = AND monomial prod(x_i), 1 = XOR parity xor(x_i) */
    flow_u128 angle;           /* Full-turn phase angle */
} qir_phase_term_t;

typedef struct {
    qcc_bv256_t *rows;       /* [out_bits], bit j in row i => xor term from input bit j */
    qcc_bv256_t  const_vec;  /* affine constant C: y = M*x ^ C */
    uint16_t     out_bits;
    uint16_t     in_bits;
} qir_linf2_payload_t;

typedef struct {
    uint16_t  key_bits;
    uint16_t  val_bits;
    uint16_t  val_words;
    uint32_t  num_entries;
    uint64_t *table; /* [num_entries * val_words] */
} qir_lut_payload_t;

typedef struct {
    uint16_t frac_a;
    uint16_t frac_b;
    uint16_t frac_out;
    uint8_t  keep_exact;
} qir_mul_payload_t;

typedef struct {
    uint16_t m;
    uint32_t irred_poly;
} qir_gf2m_payload_t;

typedef struct {
    uint32_t          num_terms;
    uint16_t          total_in_bits;
    qir_phase_term_t *terms;
} qir_phase_poly_payload_t;

typedef struct {
    qcc_kmb_circuit_t *kmb;       /* Owned or arena-tracked subcircuit */
    uint32_t           num_q;     /* Active mapped qubits */
    uint32_t          *in_q_idx;  /* KMB qubit index corresponding to flattened input bit i */
    uint32_t          *out_q_idx; /* KMB qubit index for flattened output bit i */
    uint16_t           in_bits;
    uint16_t           out_bits;
    uint8_t            owns_kmb;
} qir_opaque_kmb_payload_t;

typedef struct {
    uint32_t      id;
    qir_op_kind_t kind;
    uint32_t      flags;
    uint32_t     *ins;
    uint16_t      num_ins;
    uint32_t      outs[4];
    uint8_t       num_outs;
    uint32_t      pred_val;       /* UINT32_MAX if unconditional, else 1-bit QIR_REP_CLASSICAL value */
    uint32_t      approx_tol_q16; /* Max acceptable error in Q16 LSBs when QIR_FLAG_APPROX set */
    union {
        struct {
            uint64_t val[QIR_WORDS];
        } cst;
        struct {
            uint16_t lo;
            uint16_t hi;
        } slice;
        qir_linf2_payload_t      linf2;
        qir_lut_payload_t        lut;
        qir_mul_payload_t        mul;
        qir_gf2m_payload_t       gf2m;
        struct {
            uint16_t mant_bits;
        } clz;
        struct {
            uint8_t mode; /* 0 = XOR, 1 = AND, 2 = NOT, 3 = OR_REDUCE */
        } bit_op;
        struct {
            uint16_t max_shift;
        } shift;
        qir_phase_poly_payload_t phase_poly;
        qir_opaque_kmb_payload_t opaque;
    } u;
} qir_op_t;

typedef struct {
    qcc_arena_t *arena;
    qir_val_t   *vals;
    uint32_t     num_vals;
    uint32_t     cap_vals;

    qir_op_t    *ops;
    uint32_t     num_ops;
    uint32_t     cap_ops;

    uint32_t    *inputs;
    uint32_t     num_inputs;
    uint32_t     cap_inputs;

    uint32_t    *outputs;
    uint32_t     num_outputs;
    uint32_t     cap_outputs;

    uint32_t    *classical_outputs;
    uint32_t     num_classical_outputs;
    uint32_t     cap_classical_outputs;

    bool         require_clean_scratch;
    const char  *name;
} qir_module_t;

typedef struct {
    uint64_t word_lanes[QIR_WORDS]; /* 64 parallel states across up to 256 bits */
} qir_state64_t;

typedef struct {
    qir_state64_t *vals;       /* [module->num_vals] */
    flow_u128      phase[64];  /* Diagonal global phase per lane */
    uint8_t        has_hmr_rand;
} qir_eval64_out_t;

static inline const char *qir_op_kind_name(qir_op_kind_t k) {
    switch (k) {
        case QIR_OP_CONST:      return "const";
        case QIR_OP_SLICE:      return "slice";
        case QIR_OP_CONCAT:     return "concat";
        case QIR_OP_LINF2:      return "linf2";
        case QIR_OP_AND_VEC:    return "and_vec";
        case QIR_OP_ADD:        return "add";
        case QIR_OP_SUB:        return "sub";
        case QIR_OP_CS_ADD3:    return "cs_add3";
        case QIR_OP_MUL_FXP:    return "mul_fxp";
        case QIR_OP_SQR_FXP:    return "sqr_fxp";
        case QIR_OP_GF2M_ADD:   return "gf2m_add";
        case QIR_OP_GF2M_MUL:   return "gf2m_mul";
        case QIR_OP_GF2M_INV:   return "gf2m_inv";
        case QIR_OP_LUT_ROM:    return "lut_rom";
        case QIR_OP_CLZ_NORM:   return "clz_norm";
        case QIR_OP_CTRL_SEL:   return "ctrl_sel";
        case QIR_OP_PHASE_POLY: return "phase_poly";
        case QIR_OP_OPAQUE_KMB: return "opaque_kmb";
        case QIR_OP_HMR_MEAS:   return "hmr_meas";
        case QIR_OP_BIT_OP:     return "bit_op";
        case QIR_OP_SHL_VAR:    return "shl_var";
        case QIR_OP_SHR_VAR:    return "shr_var";
        case QIR_OP_EQ_ZERO:    return "eq_zero";
        case QIR_OP_PHASE_DOT:  return "phase_dot";
        default:                return "unknown";
    }
}

static inline void qir_module_init(qir_module_t *m, qcc_arena_t *arena, const char *name) {
    memset(m, 0, sizeof(*m));
    m->arena = arena;
    m->name  = name ? qcc_arena_strdup(arena, name) : qcc_arena_strdup(arena, "qir_mod");
}

static inline void qir_module_free(qir_module_t *m) {
    if (!m) return;
    for (uint32_t i = 0u; i < m->num_ops; ++i) {
        if (m->ops[i].kind == QIR_OP_OPAQUE_KMB && m->ops[i].u.opaque.owns_kmb && m->ops[i].u.opaque.kmb) {
            qcc_kmb_free(m->ops[i].u.opaque.kmb);
            free(m->ops[i].u.opaque.kmb);
            m->ops[i].u.opaque.kmb = NULL;
        }
    }
    free(m->vals);
    free(m->ops);
    free(m->inputs);
    free(m->outputs);
    free(m->classical_outputs);
    memset(m, 0, sizeof(*m));
}

static inline uint32_t qir__new_val(qir_module_t *m, uint16_t width, qir_rep_t rep, uint32_t def_op, const char *name) {
    if (width > QIR_MAX_VAL_BITS) width = QIR_MAX_VAL_BITS;
    if (m->num_vals >= m->cap_vals) {
        uint32_t nc = m->cap_vals ? m->cap_vals * 2u : 64u;
        m->vals = (qir_val_t *)realloc(m->vals, (size_t)nc * sizeof(qir_val_t));
        if (!m->vals) abort();
        m->cap_vals = nc;
    }
    uint32_t id = m->num_vals++;
    qir_val_t *v = &m->vals[id];
    memset(v, 0, sizeof(*v));
    v->id     = id;
    v->width  = width;
    v->rep    = rep;
    v->def_op = def_op;
    v->name   = name ? qcc_arena_strdup(m->arena, name) : NULL;
    return id;
}

static inline uint32_t qir_add_input(qir_module_t *m, uint16_t width, const char *name) {
    uint32_t vid = qir__new_val(m, width, QIR_REP_STD, UINT32_MAX, name);
    if (m->num_inputs >= m->cap_inputs) {
        uint32_t nc = m->cap_inputs ? m->cap_inputs * 2u : 16u;
        m->inputs = (uint32_t *)realloc(m->inputs, (size_t)nc * sizeof(uint32_t));
        if (!m->inputs) abort();
        m->cap_inputs = nc;
    }
    m->inputs[m->num_inputs++] = vid;
    return vid;
}

static inline uint32_t qir_add_classical_input(qir_module_t *m, uint16_t width, const char *name) {
    uint32_t vid = qir__new_val(m, width, QIR_REP_CLASSICAL, UINT32_MAX, name);
    if (m->num_inputs >= m->cap_inputs) {
        uint32_t nc = m->cap_inputs ? m->cap_inputs * 2u : 16u;
        m->inputs = (uint32_t *)realloc(m->inputs, (size_t)nc * sizeof(uint32_t));
        if (!m->inputs) abort();
        m->cap_inputs = nc;
    }
    m->inputs[m->num_inputs++] = vid;
    return vid;
}

static inline void qir_add_output(qir_module_t *m, uint32_t val_id) {
    if (m->num_outputs >= m->cap_outputs) {
        uint32_t nc = m->cap_outputs ? m->cap_outputs * 2u : 16u;
        m->outputs = (uint32_t *)realloc(m->outputs, (size_t)nc * sizeof(uint32_t));
        if (!m->outputs) abort();
        m->cap_outputs = nc;
    }
    m->outputs[m->num_outputs++] = val_id;
}

static inline void qir_add_classical_output(qir_module_t *m, uint32_t val_id) {
    if (m->num_classical_outputs >= m->cap_classical_outputs) {
        uint32_t nc = m->cap_classical_outputs ? m->cap_classical_outputs * 2u : 16u;
        m->classical_outputs = (uint32_t *)realloc(m->classical_outputs, (size_t)nc * sizeof(uint32_t));
        if (!m->classical_outputs) abort();
        m->cap_classical_outputs = nc;
    }
    m->classical_outputs[m->num_classical_outputs++] = val_id;
}

static inline qir_op_t *qir__new_op(
    qir_module_t *m,
    qir_op_kind_t kind,
    const uint32_t *ins,
    uint16_t num_ins,
    uint32_t flags
) {
    if (m->num_ops >= m->cap_ops) {
        uint32_t nc = m->cap_ops ? m->cap_ops * 2u : 64u;
        m->ops = (qir_op_t *)realloc(m->ops, (size_t)nc * sizeof(qir_op_t));
        if (!m->ops) abort();
        m->cap_ops = nc;
    }
    uint32_t op_id = m->num_ops++;
    qir_op_t *op = &m->ops[op_id];
    memset(op, 0, sizeof(*op));
    op->id       = op_id;
    op->kind     = kind;
    op->flags    = flags;
    op->pred_val = UINT32_MAX;
    op->num_ins  = num_ins;
    if (num_ins > 0u) {
        op->ins = (uint32_t *)qcc_arena_xalloc(m->arena, (size_t)num_ins * sizeof(uint32_t), 8u);
        memcpy(op->ins, ins, (size_t)num_ins * sizeof(uint32_t));
    }
    return op;
}

/* SSA Builder Constructors */
static inline uint32_t qir_emit_const(qir_module_t *m, uint16_t width, const uint64_t val[QIR_WORDS], const char *name) {
    qir_op_t *op = qir__new_op(m, QIR_OP_CONST, NULL, 0u, 0u);
    uint32_t vid = qir__new_val(m, width, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    qcc_bv256_t mask = qcc_bv_mask_n(width);
    for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
        uint64_t v = val ? (val[w] & mask.w[w]) : 0u;
        op->u.cst.val[w]        = v;
        m->vals[vid].const_mask[w] = mask.w[w];
        m->vals[vid].const_val[w]  = v;
        m->vals[vid].zero_mask[w]  = (~v) & mask.w[w];
    }
    m->vals[vid].flags |= QIR_FLAG_ZERO_MASK_KNOWN;
    return vid;
}

static inline uint32_t qir_emit_const_u64(qir_module_t *m, uint16_t width, uint64_t v, const char *name) {
    uint64_t words[QIR_WORDS] = {v, 0u, 0u, 0u};
    return qir_emit_const(m, width, words, name);
}

static inline uint32_t qir_emit_slice(qir_module_t *m, uint32_t x, uint16_t lo, uint16_t hi, const char *name) {
    uint16_t w = (hi > lo) ? (uint16_t)(hi - lo) : 0u;
    uint32_t ins[1] = {x};
    qir_op_t *op = qir__new_op(m, QIR_OP_SLICE, ins, 1u, 0u);
    op->u.slice.lo = lo;
    op->u.slice.hi = hi;
    uint32_t vid = qir__new_val(m, w, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_concat(qir_module_t *m, const uint32_t *ins, uint16_t num_ins, const char *name) {
    uint32_t w = 0u;
    for (uint16_t i = 0u; i < num_ins; ++i) w += m->vals[ins[i]].width;
    if (w > QIR_MAX_VAL_BITS) w = QIR_MAX_VAL_BITS;
    qir_op_t *op = qir__new_op(m, QIR_OP_CONCAT, ins, num_ins, 0u);
    uint32_t vid = qir__new_val(m, (uint16_t)w, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_linf2(
    qir_module_t *m,
    const uint32_t *ins,
    uint16_t num_ins,
    uint16_t out_bits,
    const qcc_bv256_t *rows,
    const qcc_bv256_t *const_vec,
    const char *name
) {
    uint16_t in_bits = 0u;
    for (uint16_t i = 0u; i < num_ins; ++i) in_bits = (uint16_t)(in_bits + m->vals[ins[i]].width);
    qir_op_t *op = qir__new_op(m, QIR_OP_LINF2, ins, num_ins, 0u);
    op->u.linf2.out_bits = out_bits;
    op->u.linf2.in_bits  = in_bits;
    op->u.linf2.rows     = out_bits ? (qcc_bv256_t *)qcc_arena_xalloc_zero(m->arena, (size_t)out_bits * sizeof(qcc_bv256_t), 8u) : NULL;
    if (rows && out_bits > 0u) memcpy(op->u.linf2.rows, rows, (size_t)out_bits * sizeof(qcc_bv256_t));
    if (const_vec) op->u.linf2.const_vec = *const_vec;
    uint32_t vid = qir__new_val(m, out_bits, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_and_vec(qir_module_t *m, uint32_t a, uint32_t b, uint32_t flags, const char *name) {
    uint16_t w = m->vals[a].width < m->vals[b].width ? m->vals[a].width : m->vals[b].width;
    uint32_t ins[2] = {a, b};
    qir_op_t *op = qir__new_op(m, QIR_OP_AND_VEC, ins, 2u, flags);
    uint32_t vid = qir__new_val(m, w, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_add(qir_module_t *m, uint32_t a, uint32_t b, uint32_t flags, const char *name) {
    uint16_t w = m->vals[a].width > m->vals[b].width ? m->vals[a].width : m->vals[b].width;
    uint32_t ins[2] = {a, b};
    qir_op_t *op = qir__new_op(m, QIR_OP_ADD, ins, 2u, flags);
    uint32_t vid = qir__new_val(m, w, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_sub(qir_module_t *m, uint32_t a, uint32_t b, uint32_t flags, const char *name) {
    uint16_t w = m->vals[a].width > m->vals[b].width ? m->vals[a].width : m->vals[b].width;
    uint32_t ins[2] = {a, b};
    qir_op_t *op = qir__new_op(m, QIR_OP_SUB, ins, 2u, flags);
    uint32_t vid = qir__new_val(m, w, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline void qir_emit_cs_add3(
    qir_module_t *m, uint32_t a, uint32_t b, uint32_t c,
    uint32_t *out_s, uint32_t *out_c, const char *name_s, const char *name_c
) {
    uint16_t w = m->vals[a].width;
    if (m->vals[b].width > w) w = m->vals[b].width;
    if (m->vals[c].width > w) w = m->vals[c].width;
    uint32_t ins[3] = {a, b, c};
    qir_op_t *op = qir__new_op(m, QIR_OP_CS_ADD3, ins, 3u, 0u);
    uint32_t vs  = qir__new_val(m, w, QIR_REP_CS, op->id, name_s);
    uint32_t vc  = qir__new_val(m, w, QIR_REP_CS, op->id, name_c);
    op->outs[0]  = vs;
    op->outs[1]  = vc;
    op->num_outs = 2u;
    if (out_s) *out_s = vs;
    if (out_c) *out_c = vc;
}

static inline uint32_t qir_emit_mul_fxp(
    qir_module_t *m, uint32_t a, uint32_t b,
    uint16_t out_width, uint16_t frac_a, uint16_t frac_b, uint16_t frac_out,
    uint32_t flags, const char *name
) {
    uint32_t ins[2] = {a, b};
    qir_op_t *op = qir__new_op(m, QIR_OP_MUL_FXP, ins, 2u, flags);
    op->u.mul.frac_a     = frac_a;
    op->u.mul.frac_b     = frac_b;
    op->u.mul.frac_out   = frac_out;
    op->u.mul.keep_exact = (flags & QIR_FLAG_APPROX) ? 0u : 1u;
    uint32_t vid = qir__new_val(m, out_width, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_sqr_fxp(
    qir_module_t *m, uint32_t a,
    uint16_t out_width, uint16_t frac_a, uint16_t frac_out,
    uint32_t flags, const char *name
) {
    uint32_t ins[1] = {a};
    qir_op_t *op = qir__new_op(m, QIR_OP_SQR_FXP, ins, 1u, flags);
    op->u.mul.frac_a     = frac_a;
    op->u.mul.frac_b     = frac_a;
    op->u.mul.frac_out   = frac_out;
    op->u.mul.keep_exact = (flags & QIR_FLAG_APPROX) ? 0u : 1u;
    uint32_t vid = qir__new_val(m, out_width, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_gf2m_add(qir_module_t *m, uint32_t a, uint32_t b, uint16_t deg, const char *name) {
    uint32_t ins[2] = {a, b};
    qir_op_t *op = qir__new_op(m, QIR_OP_GF2M_ADD, ins, 2u, 0u);
    op->u.gf2m.m          = deg;
    op->u.gf2m.irred_poly = 0u;
    uint32_t vid = qir__new_val(m, deg, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_gf2m_mul(
    qir_module_t *m, uint32_t a, uint32_t b, uint16_t deg, uint32_t irred_poly, const char *name
) {
    uint32_t ins[2] = {a, b};
    qir_op_t *op = qir__new_op(m, QIR_OP_GF2M_MUL, ins, 2u, 0u);
    op->u.gf2m.m          = deg;
    op->u.gf2m.irred_poly = irred_poly;
    uint32_t vid = qir__new_val(m, deg, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_gf2m_inv(
    qir_module_t *m, uint32_t a, uint16_t deg, uint32_t irred_poly, const char *name
) {
    uint32_t ins[1] = {a};
    qir_op_t *op = qir__new_op(m, QIR_OP_GF2M_INV, ins, 1u, 0u);
    op->u.gf2m.m          = deg;
    op->u.gf2m.irred_poly = irred_poly;
    uint32_t vid = qir__new_val(m, deg, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_lut_rom(
    qir_module_t *m, uint32_t key_val,
    uint16_t val_bits, uint32_t num_entries, const uint64_t *table,
    uint32_t flags, const char *name
) {
    uint16_t val_words = (uint16_t)((val_bits + 63u) / 64u);
    if (val_words == 0u) val_words = 1u;
    uint32_t ins[1] = {key_val};
    qir_op_t *op = qir__new_op(m, QIR_OP_LUT_ROM, ins, 1u, flags);
    op->u.lut.key_bits    = m->vals[key_val].width;
    op->u.lut.val_bits    = val_bits;
    op->u.lut.val_words   = val_words;
    op->u.lut.num_entries = num_entries;
    size_t n_words = (size_t)num_entries * (size_t)val_words;
    op->u.lut.table = n_words ? (uint64_t *)qcc_arena_xalloc_zero(m->arena, n_words * sizeof(uint64_t), 8u) : NULL;
    if (table && n_words > 0u) memcpy(op->u.lut.table, table, n_words * sizeof(uint64_t));
    uint32_t vid = qir__new_val(m, val_bits, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline void qir_emit_clz_norm(
    qir_module_t *m, uint32_t x, uint16_t mant_bits,
    uint32_t *out_mant, uint32_t *out_clz, const char *name_mant, const char *name_clz
) {
    uint16_t in_w = m->vals[x].width;
    if (mant_bits == 0u || mant_bits > in_w) mant_bits = in_w;
    uint16_t clz_w = (uint16_t)qcc_ceil_log2_u32((uint32_t)in_w + 1u);
    if (clz_w == 0u) clz_w = 1u;
    uint32_t ins[1] = {x};
    qir_op_t *op = qir__new_op(m, QIR_OP_CLZ_NORM, ins, 1u, 0u);
    op->u.clz.mant_bits = mant_bits;
    uint32_t vm = qir__new_val(m, mant_bits, QIR_REP_STD, op->id, name_mant);
    uint32_t vc = qir__new_val(m, clz_w, QIR_REP_STD, op->id, name_clz);
    op->outs[0]  = vm;
    op->outs[1]  = vc;
    op->num_outs = 2u;
    if (out_mant) *out_mant = vm;
    if (out_clz)  *out_clz  = vc;
}

static inline uint32_t qir_emit_ctrl_sel(
    qir_module_t *m, uint32_t sel_bit, uint32_t v_false, uint32_t v_true, uint32_t flags, const char *name
) {
    uint16_t w = m->vals[v_false].width > m->vals[v_true].width ? m->vals[v_false].width : m->vals[v_true].width;
    uint32_t ins[3] = {sel_bit, v_false, v_true};
    qir_op_t *op = qir__new_op(m, QIR_OP_CTRL_SEL, ins, 3u, flags);
    uint32_t vid = qir__new_val(m, w, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_phase_poly(
    qir_module_t *m, const uint32_t *ins, uint16_t num_ins,
    const qir_phase_term_t *terms, uint32_t num_terms, const char *name
) {
    uint16_t total_bits = 0u;
    for (uint16_t i = 0u; i < num_ins; ++i) total_bits = (uint16_t)(total_bits + m->vals[ins[i]].width);
    qir_op_t *op = qir__new_op(m, QIR_OP_PHASE_POLY, ins, num_ins, QIR_FLAG_PHASE_ONLY);
    op->u.phase_poly.num_terms     = num_terms;
    op->u.phase_poly.total_in_bits = total_bits;
    op->u.phase_poly.terms = num_terms
        ? (qir_phase_term_t *)qcc_arena_xalloc(m->arena, (size_t)num_terms * sizeof(qir_phase_term_t), 16u)
        : NULL;
    if (terms && num_terms > 0u) memcpy(op->u.phase_poly.terms, terms, (size_t)num_terms * sizeof(qir_phase_term_t));
    uint32_t vid = qir__new_val(m, 0u, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_opaque_kmb(
    qir_module_t *m, const uint32_t *ins, uint16_t num_ins,
    uint16_t out_bits, const qcc_kmb_circuit_t *kmb,
    const uint32_t *in_q_idx, uint16_t in_bits,
    const uint32_t *out_q_idx, const char *name
) {
    qir_op_t *op = qir__new_op(m, QIR_OP_OPAQUE_KMB, ins, num_ins, 0u);
    qcc_kmb_circuit_t *cloned = (qcc_kmb_circuit_t *)malloc(sizeof(qcc_kmb_circuit_t));
    if (!cloned || !qcc_kmb_clone(cloned, kmb)) abort();
    op->u.opaque.kmb      = cloned;
    op->u.opaque.owns_kmb = 1u;
    op->u.opaque.num_q    = cloned->num_qubits;
    op->u.opaque.in_bits  = in_bits;
    op->u.opaque.out_bits = out_bits;
    op->u.opaque.in_q_idx = in_bits ? (uint32_t *)qcc_arena_xalloc(m->arena, (size_t)in_bits * sizeof(uint32_t), 4u) : NULL;
    if (in_bits && in_q_idx) memcpy(op->u.opaque.in_q_idx, in_q_idx, (size_t)in_bits * sizeof(uint32_t));
    op->u.opaque.out_q_idx = out_bits ? (uint32_t *)qcc_arena_xalloc(m->arena, (size_t)out_bits * sizeof(uint32_t), 4u) : NULL;
    if (out_bits && out_q_idx) memcpy(op->u.opaque.out_q_idx, out_q_idx, (size_t)out_bits * sizeof(uint32_t));

    uint32_t vid = qir__new_val(m, out_bits, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_hmr_meas(qir_module_t *m, uint32_t qval, const char *name) {
    uint16_t w = m->vals[qval].width;
    uint32_t ins[1] = {qval};
    qir_op_t *op = qir__new_op(m, QIR_OP_HMR_MEAS, ins, 1u, 0u);
    uint32_t vid = qir__new_val(m, w, QIR_REP_CLASSICAL, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_bit_op(
    qir_module_t *m, uint8_t mode, uint32_t a, uint32_t b, const char *name
) {
    uint16_t w = (mode == 3u) ? 1u : m->vals[a].width;
    uint32_t ins[2] = {a, b};
    uint16_t n_in = (mode == 2u || mode == 3u) ? 1u : 2u;
    if (mode < 2u && m->vals[b].width > w) w = m->vals[b].width;
    qir_op_t *op = qir__new_op(m, QIR_OP_BIT_OP, ins, n_in, 0u);
    op->u.bit_op.mode = mode;
    uint32_t vid = qir__new_val(m, w, QIR_REP_CLASSICAL, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_gf2m_inv_ctrl(
    qir_module_t *m, uint32_t a, uint32_t ctrl_q, uint16_t deg, uint32_t irred_poly, const char *name
) {
    uint32_t ins[2] = {a, ctrl_q};
    qir_op_t *op = qir__new_op(m, QIR_OP_GF2M_INV, ins, 2u, 0u);
    op->u.gf2m.m          = deg;
    op->u.gf2m.irred_poly = irred_poly;
    uint32_t vid = qir__new_val(m, deg, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_shl_var(
    qir_module_t *m, uint32_t val, uint32_t shift_amt, uint16_t out_w, uint16_t max_shift, const char *name
) {
    uint32_t ins[2] = {val, shift_amt};
    qir_op_t *op = qir__new_op(m, QIR_OP_SHL_VAR, ins, 2u, 0u);
    op->u.shift.max_shift = max_shift;
    uint32_t vid = qir__new_val(m, out_w, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_shr_var(
    qir_module_t *m, uint32_t val, uint32_t shift_amt, uint16_t out_w, uint16_t max_shift, const char *name
) {
    uint32_t ins[2] = {val, shift_amt};
    qir_op_t *op = qir__new_op(m, QIR_OP_SHR_VAR, ins, 2u, 0u);
    op->u.shift.max_shift = max_shift;
    uint32_t vid = qir__new_val(m, out_w, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_eq_zero(
    qir_module_t *m, uint32_t val, const char *name
) {
    uint32_t ins[1] = {val};
    qir_op_t *op = qir__new_op(m, QIR_OP_EQ_ZERO, ins, 1u, 0u);
    uint32_t vid = qir__new_val(m, 1u, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_phase_dot(
    qir_module_t *m, uint32_t q_val, uint32_t c_val, const char *name
) {
    uint32_t ins[2] = {q_val, c_val};
    qir_op_t *op = qir__new_op(m, QIR_OP_PHASE_DOT, ins, 2u, QIR_FLAG_PHASE_ONLY);
    uint32_t vid = qir__new_val(m, 0u, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline uint32_t qir_emit_phase_dot_gf2m(
    qir_module_t *m, uint32_t qa, uint32_t qb, uint32_t c_val, uint16_t deg, uint32_t irred_poly, const char *name
) {
    uint32_t ins[3] = {qa, qb, c_val};
    qir_op_t *op = qir__new_op(m, QIR_OP_PHASE_DOT, ins, 3u, QIR_FLAG_PHASE_ONLY);
    op->u.gf2m.m          = deg;
    op->u.gf2m.irred_poly = irred_poly;
    uint32_t vid = qir__new_val(m, 0u, QIR_REP_STD, op->id, name);
    op->outs[0]  = vid;
    op->num_outs = 1u;
    return vid;
}

static inline void qir_set_op_pred(qir_module_t *m, uint32_t op_id, uint32_t pred_val) {
    if (m && op_id < m->num_ops) m->ops[op_id].pred_val = pred_val;
}

static inline uint16_t qir__bit_width_u64(uint64_t v) {
    uint16_t w = 0u;
    while (v > 0u) { v >>= 1u; w++; }
    return w ? w : 1u;
}

/* Declarative piecewise 2-Jet Taylor expander for invsqrt (fn_kind=0) and recip (fn_kind=1). */
static inline uint32_t qir_emit_approx_fn(
    qir_module_t *m,
    uint32_t x_val,
    uint8_t fn_kind,
    uint16_t k_seg,
    uint16_t deg,
    uint16_t out_w,
    uint16_t frac_out,
    const char *name
) {
    uint16_t iw = m->vals[x_val].width;
    if (iw < 2u) return qir_emit_const_u64(m, out_w ? out_w : 1u, 0u, name);
    if (k_seg < 2u) k_seg = 2u;
    if (k_seg > 8u) k_seg = 8u;
    if (deg < 1u) deg = 1u;
    if (deg > 2u) deg = 2u;
    uint16_t k_mant = (fn_kind == 0u) ? (uint16_t)(k_seg - 1u) : k_seg;
    if (k_mant > iw - 1u) k_mant = (uint16_t)(iw - 1u);
    uint16_t eff_k_seg = (fn_kind == 0u) ? (uint16_t)(k_mant + 1u) : k_mant;
    uint16_t u_bits = (uint16_t)((iw - 1u) - k_mant);
    uint16_t frac_int = (uint16_t)(frac_out + 2u);
    if (frac_int > 24u) frac_int = 24u;

    uint32_t mant = 0u, clz = 0u;
    qir_emit_clz_norm(m, x_val, iw, &mant, &clz, "afn_mant", "afn_clz");
    uint16_t clz_w = m->vals[clz].width;
    uint32_t iw_m1 = qir_emit_const_u64(m, clz_w, (uint64_t)(iw - 1u), NULL);
    uint32_t exp_val = qir_emit_sub(m, iw_m1, clz, 0u, "afn_exp");

    uint32_t key_val = UINT32_MAX;
    if (fn_kind == 0u) {
        uint32_t p_bit = qir_emit_slice(m, exp_val, 0u, 1u, "afn_p");
        if (k_mant > 0u) {
            uint32_t seg_m = qir_emit_slice(m, mant, u_bits, (uint16_t)(iw - 1u), "afn_seg");
            uint32_t c_ins[2] = {seg_m, p_bit};
            key_val = qir_emit_concat(m, c_ins, 2u, "afn_key");
        } else {
            key_val = p_bit;
        }
    } else {
        key_val = (k_mant > 0u) ? qir_emit_slice(m, mant, u_bits, (uint16_t)(iw - 1u), "afn_key")
                                : qir_emit_const_u64(m, 1u, 0u, "afn_key0");
    }

    uint32_t n_ent = UINT32_C(1) << eff_k_seg;
    uint64_t c0_arr[256] = {0}, c1_arr[256] = {0}, c2_arr[256] = {0};
    uint64_t max_c0 = 1u, max_c1 = 0u, max_c2 = 0u;
    uint64_t S = UINT64_C(1) << k_mant;
    uint64_t D = 2u * S;
    uint32_t eff_deg = (u_bits == 0u) ? 0u : deg;

    for (uint32_t e = 0u; e < n_ent; ++e) {
        uint32_t p = (fn_kind == 0u) ? (e >> k_mant) : 0u;
        uint32_t s = e & (uint32_t)(S - 1u);
        uint64_t M_den = D + 2u * (uint64_t)s + 1u;
        uint64_t g0 = 0u, t1 = 0u, c2_t = 0u, c0v = 0u, c1v = 0u, c2v = 0u;
        if (fn_kind == 0u) {
            flow_u128 num = ((flow_u128)D) << (2u * frac_int);
            flow_u128 den = (flow_u128)M_den << p;
            g0 = flow_isqrt_round_u128((num + (den >> 1u)) / den);
            t1 = (g0 + M_den / 2u) / M_den;
            flow_u128 dm2 = (flow_u128)M_den * (flow_u128)M_den;
            c2_t = (uint64_t)(((flow_u128)3u * g0 + dm2) / ((flow_u128)2u * dm2));
        } else {
            flow_u128 num = ((flow_u128)D) << frac_int;
            g0 = (uint64_t)((num + (M_den >> 1u)) / M_den);
            t1 = (2u * g0 + M_den / 2u) / M_den;
            flow_u128 dm2 = (flow_u128)M_den * (flow_u128)M_den;
            c2_t = (uint64_t)(((flow_u128)4u * g0 + dm2 / 2u) / dm2);
        }
        if (eff_deg == 0u) {
            /* Direct endpoint value at u=0 for zero sub-bits */
            c0v = (fn_kind == 0u)
                ? flow_isqrt_round_u128(((((flow_u128)S) << (2u * frac_int)) + ((((flow_u128)(S + s)) << p) >> 1u)) / (((flow_u128)(S + s)) << p))
                : (uint64_t)((((flow_u128)S << frac_int) + ((S + s) >> 1u)) / (S + s));
        } else if (eff_deg == 1u) {
            c1v = t1;
            c0v = g0 + (t1 + 1u) / 2u + (c2_t + 4u) / 8u;
        } else {
            c2v = c2_t;
            c1v = t1 + c2v;
            c0v = g0 + (t1 + 1u) / 2u + (c2v + 2u) / 4u;
        }
        /* Pre-add half-LSB rounding bias before final truncation from frac_int to frac_out */
        if (frac_int > frac_out) c0v += UINT64_C(1) << (frac_int - frac_out - 1u);
        c0_arr[e] = c0v;
        c1_arr[e] = c1v;
        c2_arr[e] = c2v;
        if (c0v > max_c0) max_c0 = c0v;
        if (c1v > max_c1) max_c1 = c1v;
        if (c2v > max_c2) max_c2 = c2v;
    }

    uint16_t w_c0 = qir__bit_width_u64(max_c0);
    uint16_t w_c1 = (eff_deg >= 1u) ? qir__bit_width_u64(max_c1) : 0u;
    uint16_t w_c2 = (eff_deg >= 2u && max_c2 > 0u) ? qir__bit_width_u64(max_c2) : 0u;
    uint16_t lut_w = (uint16_t)(w_c0 + w_c1 + w_c2);
    if (lut_w > QIR_MAX_VAL_BITS) {
        w_c2 = 0u;
        lut_w = (uint16_t)(w_c0 + w_c1);
    }
    uint16_t vw = (uint16_t)((lut_w + 63u) / 64u);
    if (vw == 0u) vw = 1u;
    uint64_t *tbl = (uint64_t *)calloc((size_t)n_ent * vw, sizeof(uint64_t));
    if (!tbl) abort();
    for (uint32_t e = 0u; e < n_ent; ++e) {
        uint64_t words[QIR_WORDS] = {0};
        for (uint16_t b = 0u; b < w_c0; ++b) {
            if ((c0_arr[e] >> b) & 1u) qcc_bv_set(words, b);
        }
        for (uint16_t b = 0u; b < w_c1; ++b) {
            if ((c1_arr[e] >> b) & 1u) qcc_bv_set(words, (uint32_t)(w_c0 + b));
        }
        for (uint16_t b = 0u; b < w_c2; ++b) {
            if ((c2_arr[e] >> b) & 1u) qcc_bv_set(words, (uint32_t)(w_c0 + w_c1 + b));
        }
        for (uint16_t w = 0u; w < vw && w < QIR_WORDS; ++w) tbl[(size_t)e * vw + w] = words[w];
    }
    uint32_t lut_val = qir_emit_lut_rom(m, key_val, lut_w, n_ent, tbl, 0u, "afn_lut");
    free(tbl);

    uint32_t c0_val = qir_emit_slice(m, lut_val, 0u, w_c0, "afn_c0");
    uint32_t cur_y = c0_val;
    if (eff_deg >= 1u && u_bits > 0u) {
        uint32_t u_val = qir_emit_slice(m, mant, 0u, u_bits, "afn_u");
        uint32_t c1_val = qir_emit_slice(m, lut_val, w_c0, (uint16_t)(w_c0 + w_c1), "afn_c1");
        uint32_t p1 = qir_emit_mul_fxp(m, c1_val, u_val, w_c1, frac_int, u_bits, frac_int, QIR_FLAG_APPROX, "afn_p1");
        m->ops[m->vals[p1].def_op].approx_tol_q16 = (2u << 16u);
        cur_y = qir_emit_sub(m, c0_val, p1, 0u, "afn_lin");
        if (w_c2 > 0u && u_bits >= 2u) {
            uint32_t c2_val = qir_emit_slice(m, lut_val, (uint16_t)(w_c0 + w_c1), lut_w, "afn_c2");
            uint16_t u2_w = (u_bits > 8u) ? 8u : u_bits;
            uint32_t u2 = qir_emit_sqr_fxp(m, u_val, u2_w, u_bits, u2_w, QIR_FLAG_APPROX, "afn_u2");
            m->ops[m->vals[u2].def_op].approx_tol_q16 = (2u << 16u);
            uint32_t p2 = qir_emit_mul_fxp(m, c2_val, u2, w_c2, frac_int, u2_w, frac_int, QIR_FLAG_APPROX, "afn_p2");
            m->ops[m->vals[p2].def_op].approx_tol_q16 = (2u << 16u);
            cur_y = qir_emit_add(m, cur_y, p2, 0u, "afn_quad");
        }
    }

    /* Barrel right-shift by h = (exp >> 1) for invsqrt or exp for recip */
    uint16_t sh_lo = (fn_kind == 0u) ? 1u : 0u;
    uint16_t cur_w = m->vals[cur_y].width;
    if (clz_w > sh_lo) {
        uint32_t h_val = (sh_lo == 0u) ? exp_val : qir_emit_slice(m, exp_val, 1u, clz_w, "afn_h");
        uint16_t h_w = m->vals[h_val].width;
        for (uint16_t b = 0u; b < h_w; ++b) {
            uint32_t sh_amt = UINT32_C(1) << b;
            if (sh_amt >= cur_w) break;
            uint32_t sb = qir_emit_slice(m, h_val, b, (uint16_t)(b + 1u), NULL);
            uint32_t hi_part = qir_emit_slice(m, cur_y, (uint16_t)sh_amt, cur_w, NULL);
            uint32_t z_hi = qir_emit_const_u64(m, (uint16_t)sh_amt, 0u, NULL);
            uint32_t cins[2] = {hi_part, z_hi};
            uint32_t shifted = qir_emit_concat(m, cins, 2u, NULL);
            cur_y = qir_emit_ctrl_sel(m, sb, cur_y, shifted, 0u, "afn_rsh");
        }
    }

    uint16_t drop = (frac_int > frac_out) ? (uint16_t)(frac_int - frac_out) : 0u;
    if (cur_w >= drop + out_w) {
        return qir_emit_slice(m, cur_y, drop, (uint16_t)(drop + out_w), name);
    }
    uint32_t tail = (cur_w > drop) ? qir_emit_slice(m, cur_y, drop, cur_w, NULL)
                                   : qir_emit_const_u64(m, out_w, 0u, name);
    if (cur_w <= drop) return tail;
    uint16_t rem_z = (uint16_t)(out_w - (cur_w - drop));
    uint32_t pz = qir_emit_const_u64(m, rem_z, 0u, NULL);
    uint32_t oins[2] = {tail, pz};
    return qir_emit_concat(m, oins, 2u, name);
}

typedef void (*qir_phase_target_fn_t)(uint64_t x, uint64_t *out_target_u64, uint64_t *out_tol_u64, void *ctx);

static inline double qir__sqrt_newton_d(double v) {
    if (v <= 0.0) return 0.0;
    double x = (v >= 1.0) ? v : 1.0;
    for (int i = 0; i < 40; ++i) {
        double nxt = 0.5 * (x + v / x);
        if (nxt == x) break;
        x = nxt;
    }
    return x;
}

static inline uint64_t qir__isqrt_floor_u128(flow_u128 n) {
    if (n == 0u) return 0u;
    uint64_t hi = (uint64_t)(n >> 64u);
    uint32_t msb = hi ? (64u + 63u - (uint32_t)__builtin_clzll((unsigned long long)hi))
                      : (63u - (uint32_t)__builtin_clzll((unsigned long long)n));
    uint32_t sh = (msb >> 1u) + 1u;
    flow_u128 x = (sh >= 64u) ? (flow_u128)UINT64_MAX : ((flow_u128)1u << sh);
    for (;;) {
        flow_u128 nx = (x + n / x) >> 1u;
        if (nx >= x) break;
        x = nx;
    }
    return (uint64_t)x;
}

static inline void qir__builtin_invsqrt_ref(uint64_t x, uint64_t *out_target_u64, uint64_t *out_tol_u64, void *user_ctx) {
    (void)user_ctx;
    if (x == 1u) { *out_target_u64 = 0ULL; *out_tol_u64 = UINT64_C(1) << 40u; return; }
    if (x == 2u) { *out_target_u64 = UINT64_C(13043817825332782212); *out_tol_u64 = UINT64_C(777472127993); return; }
    if (x == 3u) { *out_target_u64 = UINT64_C(10650232656628343401); *out_tol_u64 = UINT64_C(634803334273); return; }
    if (x == 4u) { *out_target_u64 = UINT64_C(1) << 63u; *out_tol_u64 = UINT64_C(1) << 39u; return; }
    flow_u128 top  = ((flow_u128)1 << 66u);
    flow_u128 q_hi = top / x;
    flow_u128 rem  = top % x;
    flow_u128 q128 = (q_hi << 64u) | ((rem << 64u) / x);
    uint64_t r = qir__isqrt_floor_u128(q128);
    *out_target_u64 = (r + 1ULL) >> 1u;
    *out_tol_u64    = r >> 25u;
}

static inline void qir__solve_3x3_d(double A[3][3], const double b[3], double out_x[3]) {
    double M[3][4];
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) M[r][c] = A[r][c];
        M[r][3] = b[r];
    }
    for (int col = 0; col < 3; ++col) {
        int piv = col;
        double best = M[col][col] < 0.0 ? -M[col][col] : M[col][col];
        for (int r = col + 1; r < 3; ++r) {
            double v = M[r][col] < 0.0 ? -M[r][col] : M[r][col];
            if (v > best) { best = v; piv = r; }
        }
        if (piv != col) {
            for (int c = col; c < 4; ++c) { double t = M[col][c]; M[col][c] = M[piv][c]; M[piv][c] = t; }
        }
        double inv_p = (M[col][col] != 0.0) ? (1.0 / M[col][col]) : 0.0;
        for (int c = col; c < 4; ++c) M[col][c] *= inv_p;
        for (int r = 0; r < 3; ++r) {
            if (r == col) continue;
            double f = M[r][col];
            for (int c = col; c < 4; ++c) M[r][c] -= f * M[col][c];
        }
    }
    out_x[0] = M[0][3]; out_x[1] = M[1][3]; out_x[2] = M[2][3];
}

/* General Declarative Dyadic Piecewise Quadratic Phase-Oracle Builder */
static inline bool qir_emit_dyadic_phase_fn(
    qir_module_t *m,
    uint32_t x_val,
    uint64_t x_min,
    uint64_t x_max,
    uint16_t k_mant,
    uint16_t frac_bits,
    uint16_t poly_bits,
    qir_phase_target_fn_t target_fn,
    void *user_ctx
) {
    if (!m || x_val >= m->num_vals || frac_bits == 0u || frac_bits > 63u || poly_bits == 0u || poly_bits > 60u) return false;
    uint16_t iw = m->vals[x_val].width;
    if (iw < 4u || iw > 63u || k_mant < 1u || k_mant > 10u || k_mant >= iw) return false;
    if (x_min == 0u) x_min = 1u;
    uint64_t max_in = (UINT64_C(1) << iw) - 1u;
    if (x_max == 0u || x_max > max_in) x_max = max_in;
    if (x_min > x_max) return false;

    uint16_t max_sh   = (uint16_t)((iw - 1u) >> 1u);
    uint16_t c_scale  = (frac_bits > max_sh) ? (uint16_t)(frac_bits - max_sh) : 0u;
    if (c_scale >= poly_bits) return false;
    uint32_t n_charts = UINT32_C(1) << (k_mant + 1u);
    uint64_t S        = UINT64_C(1) << k_mant;
    uint16_t w_dx     = (uint16_t)((iw - 1u > k_mant) ? ((iw - 1u - k_mant > 14u) ? 14u : (iw - 1u - k_mant)) : 0u);
    uint16_t w_dx2    = (w_dx > 8u) ? 8u : w_dx;
    uint16_t sh_dx2   = (w_dx >= w_dx2) ? (uint16_t)(2u * w_dx - w_dx2) : 0u;
    uint16_t sh1      = w_dx;
    uint16_t sh2      = (uint16_t)(w_dx2 + 2u);
    uint32_t guard    = 4u;
    uint32_t minp_sq  = (sh_dx2 > guard) ? (sh_dx2 - guard) : 0u;
    uint32_t minp_1   = (sh1 > guard) ? (sh1 - guard) : 0u;
    uint64_t poly_msk = (UINT64_C(1) << poly_bits) - 1u;
    uint64_t frac_msk = (UINT64_C(1) << frac_bits) - 1u;
    uint64_t dx_msk   = w_dx ? ((UINT64_C(1) << w_dx) - 1u) : 0u;

    uint64_t *c0_arr = (uint64_t *)calloc(n_charts, sizeof(uint64_t));
    uint64_t *c1_arr = (uint64_t *)calloc(n_charts, sizeof(uint64_t));
    uint64_t *c2_arr = (uint64_t *)calloc(n_charts, sizeof(uint64_t));
    if (!c0_arr || !c1_arr || !c2_arr) abort();

    static const double CHEB_U[3] = {
        0.06698729810778067, 0.5, 0.9330127018922193
    };
    const double inv_sqrt2 = 0.70710678118654752440;
    const double two_cscale = (double)(UINT64_C(1) << c_scale);

    for (uint32_t e = 0u; e < n_charts; ++e) {
        uint32_t p = e >> k_mant;
        uint32_t s = e & (uint32_t)(S - 1u);
        double scale = two_cscale * ((p == 0u) ? inv_sqrt2 : 1.0);
        double A[3][3], b[3], sol[3];
        for (int k = 0; k < 3; ++k) {
            double u = CHEB_U[k];
            double m_val = 1.0 + ((double)s + u) / (double)S;
            double gv;
            if (!target_fn || target_fn == qir__builtin_invsqrt_ref) {
                gv = scale / qir__sqrt_newton_d(m_val);
            } else {
                uint64_t x_s = (uint64_t)(m_val * (double)(UINT64_C(1) << (iw - 1u - p)) + 0.5);
                if (x_s < x_min) x_s = x_min;
                if (x_s > x_max) x_s = x_max;
                uint64_t tu = 0u, to = 0u;
                target_fn(x_s, &tu, &to, user_ctx);
                gv = ((double)tu / 18446744073709551616.0) * (double)(UINT64_C(1) << (frac_bits - (p >> 1u)));
            }
            A[k][0] = 1.0; A[k][1] = -u; A[k][2] = u * u;
            b[k] = gv;
        }
        qir__solve_3x3_d(A, b, sol);
        int64_t ic2 = (int64_t)(sol[2] * 4.0 + 0.5);
        int64_t ic1 = (int64_t)(sol[1] + 0.5);
        int64_t ic0 = (int64_t)(sol[0] + 0.5);
        if (ic2 < 0) ic2 = 0;
        if (ic1 < 0) ic1 = 0;
        c2_arr[e] = (uint64_t)ic2;
        c1_arr[e] = (uint64_t)ic1;
        c0_arr[e] = (uint64_t)ic0 & poly_msk;
    }

    /* Coordinate minimax integer centering of c0_arr[e] over exact domain [x_min, x_max] */
    double *min_dp = (double *)malloc(n_charts * sizeof(double));
    double *max_dp = (double *)malloc(n_charts * sizeof(double));
    double *sum_dp = (double *)calloc(n_charts, sizeof(double));
    uint32_t *cnt_ch = (uint32_t *)calloc(n_charts, sizeof(uint32_t));
    double *min_wt = (double *)malloc(n_charts * sizeof(double));
    uint8_t *has_ch = (uint8_t *)calloc(n_charts, 1u);
    if (!min_dp || !max_dp || !sum_dp || !cnt_ch || !min_wt || !has_ch) abort();
    for (uint32_t e = 0u; e < n_charts; ++e) {
        min_dp[e] = 1e30; max_dp[e] = -1e30; min_wt[e] = 1e30;
    }

    uint64_t span = x_max - x_min + 1u;
    uint64_t step = (span <= 4194304u) ? 1u : (span / 2000000u + 1u);
    for (uint64_t x = x_min; x <= x_max; x += step) {
        uint32_t msb = 63u - (uint32_t)__builtin_clzll((unsigned long long)x);
        uint32_t lz  = (iw - 1u) - msb;
        uint32_t sh  = lz >> 1u;
        uint32_t p   = lz & 1u;
        uint64_t mant = x << lz;
        uint32_t e   = (p << k_mant) | (uint32_t)((mant >> (iw - 1u - k_mant)) & (S - 1u));
        uint64_t dx  = w_dx ? ((mant >> (iw - 1u - k_mant - w_dx)) & dx_msk) : 0u;

        uint64_t dx2 = 0u;
        for (uint16_t i = 0u; i < w_dx; ++i) {
            if (((dx >> i) & 1u) == 0u) continue;
            uint32_t pdiag = 2u * i;
            if (pdiag >= minp_sq) dx2 += UINT64_C(1) << (pdiag - minp_sq);
            for (uint16_t j = (uint16_t)(i + 1u); j < w_dx; ++j) {
                if ((dx >> j) & 1u) {
                    uint32_t pcross = (uint32_t)i + (uint32_t)j + 1u;
                    if (pcross >= minp_sq) dx2 += UINT64_C(1) << (pcross - minp_sq);
                }
            }
        }
        dx2 >>= (sh_dx2 - minp_sq);
        dx2 &= ((UINT64_C(1) << w_dx2) - 1u);

        uint64_t c1v = c1_arr[e];
        uint64_t t1  = 0u;
        for (uint16_t i = 0u; i < w_dx; ++i) {
            if (((dx >> i) & 1u) == 0u) continue;
            uint32_t lo_c = (minp_1 > i) ? (minp_1 - i) : 0u;
            t1 += ((c1v >> lo_c) << (lo_c + i - minp_1));
        }
        t1 >>= (sh1 - minp_1);

        uint64_t t2   = (c2_arr[e] * dx2) >> sh2;
        uint64_t poly = (c0_arr[e] - t1 + t2) & poly_msk;
        uint64_t ph_u = ((poly << sh) & frac_msk) << (64u - frac_bits);
        uint64_t targ = 0u, tol = 1u;
        if (!target_fn) qir__builtin_invsqrt_ref(x, &targ, &tol, NULL);
        else            target_fn(x, &targ, &tol, user_ctx);

        int64_t sdiff = (int64_t)(ph_u - targ);
        double unit_scale = (double)(UINT64_C(1) << (64u - frac_bits + sh));
        double dp = (double)sdiff / unit_scale;
        double wt = (double)tol / unit_scale;
        has_ch[e] = 1u;
        if (dp < min_dp[e]) min_dp[e] = dp;
        if (dp > max_dp[e]) max_dp[e] = dp;
        if (wt < min_wt[e]) min_wt[e] = wt;
        sum_dp[e] += dp;
        cnt_ch[e]++;
    }
    for (uint32_t e = 0u; e < n_charts; ++e) {
        if (!has_ch[e]) continue;
        double mid = 0.5 * (min_dp[e] + max_dp[e]);
        int64_t delta = (mid >= 0.0) ? (int64_t)(mid + 0.5) : -(int64_t)(-mid + 0.5);
        int64_t best_d = delta;
        double best_r = 1e30;
        for (int64_t cand_d = delta - 2; cand_d <= delta + 2; ++cand_d) {
            double e0 = min_dp[e] - (double)cand_d; if (e0 < 0.0) e0 = -e0;
            double e1 = max_dp[e] - (double)cand_d; if (e1 < 0.0) e1 = -e1;
            double mx = (e0 > e1) ? e0 : e1;
            double r = (min_wt[e] > 0.0) ? (mx / min_wt[e]) : mx;
            if (r < best_r) { best_r = r; best_d = cand_d; }
        }
        c0_arr[e] = (uint64_t)((int64_t)c0_arr[e] - best_d) & poly_msk;
    }
    free(min_dp); free(max_dp); free(sum_dp); free(cnt_ch); free(min_wt); free(has_ch);

    uint64_t max_c1 = 1u, max_c2 = 0u;
    for (uint32_t e = 0u; e < n_charts; ++e) {
        if (c1_arr[e] > max_c1) max_c1 = c1_arr[e];
        if (c2_arr[e] > max_c2) max_c2 = c2_arr[e];
    }
    uint16_t w_c1 = qir__bit_width_u64(max_c1);
    uint16_t w_c2 = max_c2 ? qir__bit_width_u64(max_c2) : 0u;
    uint16_t w_t1 = qir__bit_width_u64((max_c1 * (dx_msk ? dx_msk : 1u)) >> sh1);
    if (w_t1 > poly_bits) w_t1 = poly_bits;
    uint64_t max_dx2 = w_dx2 ? ((UINT64_C(1) << w_dx2) - 1u) : 0u;
    uint16_t w_t2 = max_c2 ? qir__bit_width_u64((max_c2 * (max_dx2 ? max_dx2 : 1u)) >> sh2) : 0u;
    if (w_t2 > poly_bits) w_t2 = poly_bits;

    uint32_t mant = 0u, clz = 0u;
    qir_emit_clz_norm(m, x_val, iw, &mant, &clz, "dp_mant", "dp_clz");
    uint16_t clz_w = m->vals[clz].width;
    uint32_t clz_lsb  = qir_emit_slice(m, clz, 0u, 1u, "dp_p");
    uint32_t clz_half = (clz_w > 1u) ? qir_emit_slice(m, clz, 1u, clz_w, "dp_sh")
                                     : qir_emit_const_u64(m, 1u, 0u, "dp_sh0");
    uint32_t seg = qir_emit_slice(m, mant, (uint16_t)(iw - 1u - k_mant), (uint16_t)(iw - 1u), "dp_seg");
    uint32_t kins[2] = {seg, clz_lsb};
    uint32_t key = qir_emit_concat(m, kins, 2u, "dp_key");

    uint32_t dx = (w_dx > 0u)
        ? qir_emit_slice(m, mant, (uint16_t)(iw - 1u - k_mant - w_dx), (uint16_t)(iw - 1u - k_mant), "dp_dx")
        : qir_emit_const_u64(m, 1u, 0u, "dp_dx0");

    uint32_t t2 = UINT32_MAX;
    if (w_dx2 > 0u && w_c2 > 0u && w_t2 > 0u) {
        uint32_t dx2 = qir_emit_sqr_fxp(m, dx, w_dx2, w_dx, w_dx2, QIR_FLAG_APPROX, "dp_dx2");
        m->ops[m->vals[dx2].def_op].approx_tol_q16 = (guard << 16u);
        m->ops[m->vals[dx2].def_op].u.mul.keep_exact = 2u;
        uint32_t c2_v = qir_emit_lut_rom(m, key, w_c2, n_charts, c2_arr, 0u, "dp_c2");
        t2 = qir_emit_mul_fxp(m, c2_v, dx2, w_t2, sh2, 0u, 0u, 0u, "dp_t2");
    }

    uint32_t c1_v = qir_emit_lut_rom(m, key, w_c1, n_charts, c1_arr, 0u, "dp_c1");
    uint32_t t1 = qir_emit_mul_fxp(m, c1_v, dx, w_t1, sh1, 0u, 0u, QIR_FLAG_APPROX, "dp_t1");
    m->ops[m->vals[t1].def_op].approx_tol_q16 = (guard << 16u);
    m->ops[m->vals[t1].def_op].u.mul.keep_exact = 2u;

    uint32_t c0_v = qir_emit_lut_rom(m, key, poly_bits, n_charts, c0_arr, 0u, "dp_c0");
    uint32_t diff = qir_emit_sub(m, c0_v, t1, 0u, "dp_diff");
    m->vals[diff].width = poly_bits;
    uint32_t poly = diff;
    if (t2 != UINT32_MAX) {
        poly = qir_emit_add(m, diff, t2, 0u, "dp_poly");
        m->vals[poly].width = poly_bits;
    }

    uint32_t shifted = qir_emit_shl_var(m, poly, clz_half, frac_bits, max_sh, "dp_shifted");
    qir_phase_term_t *terms = (qir_phase_term_t *)calloc(frac_bits, sizeof(qir_phase_term_t));
    if (!terms) abort();
    for (uint16_t b = 0u; b < frac_bits; ++b) {
        qcc_bv_set(terms[b].mask, b);
        terms[b].is_xor = 0u;
        uint32_t exp_b = (uint32_t)(128u - frac_bits + b);
        terms[b].angle = (exp_b < 128u) ? ((flow_u128)1u << exp_b) : 0u;
    }
    (void)qir_emit_phase_poly(m, &shifted, 1u, terms, frac_bits, "dp_phase");
    m->require_clean_scratch = true;

    free(terms);
    free(c0_arr); free(c1_arr); free(c2_arr);
    return true;
}

/* ============================================================================
 * Bit-Sliced 64-Lane Helpers & Abstract Interpretation (zero_mask propagation)
 * ============================================================================ */

static inline bool qir_bv_get_word(const uint64_t w[QIR_WORDS], uint32_t bit) {
    return qcc_bv_get(w, bit);
}

static inline void qir__extract_bit_lane(const qir_state64_t *st, uint32_t bit, uint64_t *out_lane) {
    if (bit >= QIR_MAX_VAL_BITS) {
        *out_lane = 0u;
        return;
    }
    *out_lane = st->word_lanes[0]; /* placeholder unused */
    (void)bit;
}

/* Note: in 64-lane representation for multi-bit values:
 * Storing [bit][lane] or [lane][word] -> wait! If each value has up to 256 bits,
 * storing 256 uint64_t bit-lanes per value (64 parallel instances in bit b!) OR
 * storing 64 instances of qcc_bv256_t makes variable-width arithmetic simple & fast!
 * Let's check size: 64 instances * 32 bytes = 2 KB per SSA value -> for 256 ops = 512 KB.
 * And for 64-lane GF(2)/AND/XOR/CSA operations, [lane][word] lets 64-bit integers add
 * naturally with native `+`, `-`, `*`, `__builtin_clzll`, while 256-bit shifts/GF2M are 4-word loops!
 */
typedef struct {
    uint64_t lanes[64][QIR_WORDS];
} qir_val_vec64_t;

static inline void qir__u256_add(const uint64_t a[QIR_WORDS], const uint64_t b[QIR_WORDS], uint16_t w, uint64_t out[QIR_WORDS]) {
    uint64_t carry = 0u;
    for (uint32_t i = 0u; i < QIR_WORDS; ++i) {
        __uint128_t s = (__uint128_t)a[i] + (__uint128_t)b[i] + carry;
        out[i] = (uint64_t)s;
        carry  = (uint64_t)(s >> 64u);
    }
    qcc_bv256_t m = qcc_bv_mask_n(w);
    for (uint32_t i = 0u; i < QIR_WORDS; ++i) out[i] &= m.w[i];
}

static inline void qir__u256_sub(const uint64_t a[QIR_WORDS], const uint64_t b[QIR_WORDS], uint16_t w, uint64_t out[QIR_WORDS]) {
    uint64_t borrow = 0u;
    for (uint32_t i = 0u; i < QIR_WORDS; ++i) {
        __uint128_t av = a[i];
        __uint128_t bv = (__uint128_t)b[i] + borrow;
        if (av >= bv) {
            out[i] = (uint64_t)(av - bv);
            borrow = 0u;
        } else {
            out[i] = (uint64_t)((((__uint128_t)1u) << 64u) + av - bv);
            borrow = 1u;
        }
    }
    qcc_bv256_t m = qcc_bv_mask_n(w);
    for (uint32_t i = 0u; i < QIR_WORDS; ++i) out[i] &= m.w[i];
}

static inline void qir__u256_mul_shift(
    const uint64_t a[QIR_WORDS], const uint64_t b[QIR_WORDS],
    uint32_t shift, uint16_t out_w, uint64_t out[QIR_WORDS]
) {
    uint64_t prod[8] = {0};
    for (uint32_t i = 0u; i < 4u; ++i) {
        if (a[i] == 0u) continue;
        uint64_t carry = 0u;
        for (uint32_t j = 0u; j < 4u; ++j) {
            __uint128_t cur = (__uint128_t)prod[i + j] + (__uint128_t)a[i] * (__uint128_t)b[j] + carry;
            prod[i + j] = (uint64_t)cur;
            carry       = (uint64_t)(cur >> 64u);
        }
        prod[i + 4u] += carry;
    }
    memset(out, 0, sizeof(uint64_t) * QIR_WORDS);
    uint32_t ws = shift / 64u;
    uint32_t bs = shift % 64u;
    for (uint32_t i = 0u; i < QIR_WORDS; ++i) {
        uint64_t lo = (i + ws < 8u) ? prod[i + ws] : 0u;
        uint64_t hi = (i + ws + 1u < 8u) ? prod[i + ws + 1u] : 0u;
        out[i] = (bs == 0u) ? lo : ((lo >> bs) | (hi << (64u - bs)));
    }
    qcc_bv256_t m = qcc_bv_mask_n(out_w);
    for (uint32_t i = 0u; i < QIR_WORDS; ++i) out[i] &= m.w[i];
}

static inline uint32_t qir__gf8_mul_0xB(uint32_t a, uint32_t b) {
    uint32_t res = 0u;
    a &= 7u;
    b &= 7u;
    for (uint32_t i = 0u; i < 3u; ++i) {
        if ((b >> i) & 1u) res ^= a;
        a <<= 1u;
        if (a & 8u) a ^= 0xBu;
    }
    return res & 7u;
}

static inline uint32_t qir__gf2m_mul_u32(uint32_t a, uint32_t b, uint16_t deg, uint32_t irred_poly) {
    if (deg == 6u && irred_poly == 0u) {
        uint32_t x0 = a & 7u, x1 = (a >> 3u) & 7u;
        uint32_t y0 = b & 7u, y1 = (b >> 3u) & 7u;
        uint32_t p0 = qir__gf8_mul_0xB(x0, y0);
        uint32_t p1 = qir__gf8_mul_0xB(x1, y1);
        uint32_t pm = qir__gf8_mul_0xB(x0 ^ x1, y0 ^ y1);
        return (p0 ^ p1) | ((pm ^ p0) << 3u);
    }
    uint32_t res = 0u;
    uint32_t top = UINT32_C(1) << deg;
    uint32_t mask = top - 1u;
    a &= mask;
    b &= mask;
    for (uint16_t i = 0u; i < deg; ++i) {
        if ((b >> i) & 1u) res ^= a;
        a <<= 1u;
        if (a & top) a ^= irred_poly;
    }
    return res & mask;
}

static inline uint32_t qir__gf2m_inv_u32(uint32_t a, uint16_t deg, uint32_t irred_poly) {
    uint32_t mask = (deg == 32u) ? UINT32_MAX : ((UINT32_C(1) << deg) - 1u);
    a &= mask;
    if (a == 0u) return 0u;
    uint32_t exp = mask - 1u; /* 2^deg - 2 */
    uint32_t res = 1u;
    uint32_t base = a;
    while (exp > 0u) {
        if (exp & 1u) res = qir__gf2m_mul_u32(res, base, deg, irred_poly);
        base = qir__gf2m_mul_u32(base, base, deg, irred_poly);
        exp >>= 1u;
    }
    return res & mask;
}

/* Forward abstract interpretation pass for zero_mask & const propagation */
static inline void qir_recompute_zero_masks(qir_module_t *m) {
    for (uint32_t oi = 0u; oi < m->num_ops; ++oi) {
        const qir_op_t *op = &m->ops[oi];
        if (op->num_outs == 0u) continue;
        qir_val_t *ov = &m->vals[op->outs[0]];
        qcc_bv256_t mask = qcc_bv_mask_n(ov->width);
        uint64_t zm[QIR_WORDS] = {0}, cm[QIR_WORDS] = {0}, cv[QIR_WORDS] = {0};

        switch (op->kind) {
            case QIR_OP_CONST:
                for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                    cm[w] = mask.w[w];
                    cv[w] = op->u.cst.val[w] & mask.w[w];
                    zm[w] = (~cv[w]) & mask.w[w];
                }
                break;
            case QIR_OP_SLICE: {
                const qir_val_t *iv = &m->vals[op->ins[0]];
                for (uint16_t b = 0u; b < ov->width; ++b) {
                    uint16_t sb = (uint16_t)(op->u.slice.lo + b);
                    if (sb < iv->width) {
                        if (qcc_bv_get(iv->zero_mask, sb)) qcc_bv_set(zm, b);
                        if (qcc_bv_get(iv->const_mask, sb)) qcc_bv_set(cm, b);
                        if (qcc_bv_get(iv->const_val, sb)) qcc_bv_set(cv, b);
                    } else {
                        qcc_bv_set(zm, b);
                        qcc_bv_set(cm, b);
                    }
                }
                break;
            }
            case QIR_OP_CONCAT: {
                uint32_t pos = 0u;
                for (uint16_t i = 0u; i < op->num_ins; ++i) {
                    const qir_val_t *iv = &m->vals[op->ins[i]];
                    for (uint16_t b = 0u; b < iv->width && pos < QIR_MAX_VAL_BITS; ++b, ++pos) {
                        if (qcc_bv_get(iv->zero_mask, b)) qcc_bv_set(zm, pos);
                        if (qcc_bv_get(iv->const_mask, b)) qcc_bv_set(cm, pos);
                        if (qcc_bv_get(iv->const_val, b)) qcc_bv_set(cv, pos);
                    }
                }
                break;
            }
            case QIR_OP_LINF2: {
                uint64_t flat_cm[QIR_WORDS] = {0}, flat_cv[QIR_WORDS] = {0};
                uint32_t pos = 0u;
                for (uint16_t i = 0u; i < op->num_ins; ++i) {
                    const qir_val_t *iv = &m->vals[op->ins[i]];
                    for (uint16_t b = 0u; b < iv->width && pos < QIR_MAX_VAL_BITS; ++b, ++pos) {
                        if (qcc_bv_get(iv->const_mask, b)) qcc_bv_set(flat_cm, pos);
                        if (qcc_bv_get(iv->const_val, b)) qcc_bv_set(flat_cv, pos);
                    }
                }
                for (uint16_t r = 0u; r < ov->width && r < op->u.linf2.out_bits; ++r) {
                    const qcc_bv256_t *row = &op->u.linf2.rows[r];
                    uint32_t parity = 0u;
                    bool all_const = true;
                    for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                        if ((row->w[w] & ~flat_cm[w]) != 0u) all_const = false;
                        parity ^= (uint32_t)__builtin_popcountll((unsigned long long)(row->w[w] & flat_cv[w]));
                    }
                    uint32_t cbit = qcc_bv_get(op->u.linf2.const_vec.w, r) ? 1u : 0u;
                    if (all_const) {
                        qcc_bv_set(cm, r);
                        if ((parity ^ cbit) & 1u) {
                            qcc_bv_set(cv, r);
                        } else {
                            qcc_bv_set(zm, r);
                        }
                    }
                }
                break;
            }
            case QIR_OP_AND_VEC: {
                const qir_val_t *a = &m->vals[op->ins[0]];
                const qir_val_t *b = &m->vals[op->ins[1]];
                for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                    zm[w] = (a->zero_mask[w] | b->zero_mask[w]) & mask.w[w];
                    cm[w] = zm[w] | ((a->const_mask[w] & b->const_mask[w]) & mask.w[w]);
                    cv[w] = (a->const_val[w] & b->const_val[w]) & cm[w];
                }
                break;
            }
            case QIR_OP_ADD:
            case QIR_OP_SUB: {
                /* Leading trailing-zero & leading-zero bounds: if both operands have low k bits zero, sum/diff has low k bits zero */
                const qir_val_t *a = &m->vals[op->ins[0]];
                const qir_val_t *b = &m->vals[op->ins[1]];
                for (uint16_t bit = 0u; bit < ov->width; ++bit) {
                    if (qcc_bv_get(a->zero_mask, bit) && qcc_bv_get(b->zero_mask, bit)) {
                        qcc_bv_set(zm, bit);
                        qcc_bv_set(cm, bit);
                    } else {
                        break;
                    }
                }
                break;
            }
            case QIR_OP_CS_ADD3: {
                const qir_val_t *a = &m->vals[op->ins[0]];
                const qir_val_t *b = &m->vals[op->ins[1]];
                const qir_val_t *c = &m->vals[op->ins[2]];
                qir_val_t *vs = &m->vals[op->outs[0]];
                qir_val_t *vc = &m->vals[op->outs[1]];
                qcc_bv256_t ms = qcc_bv_mask_n(vs->width);
                qcc_bv256_t mc = qcc_bv_mask_n(vc->width);
                for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                    vs->zero_mask[w]  = (a->zero_mask[w] & b->zero_mask[w] & c->zero_mask[w]) & ms.w[w];
                    uint64_t pair_zero = (a->zero_mask[w] & b->zero_mask[w]) |
                                         (b->zero_mask[w] & c->zero_mask[w]) |
                                         (a->zero_mask[w] & c->zero_mask[w]);
                    vc->zero_mask[w]  = (pair_zero << 1u) & mc.w[w];
                    qcc_bv_set(vc->zero_mask, 0u);
                    if (w > 0u && ((a->zero_mask[w - 1u] & b->zero_mask[w - 1u]) |
                                   (b->zero_mask[w - 1u] & c->zero_mask[w - 1u]) |
                                   (a->zero_mask[w - 1u] & c->zero_mask[w - 1u])) >> 63u) {
                        vc->zero_mask[w] |= UINT64_C(1);
                    }
                    vc->zero_mask[w] &= mc.w[w];
                    vs->const_mask[w] = vs->zero_mask[w];
                    vc->const_mask[w] = vc->zero_mask[w];
                }
                vs->flags |= QIR_FLAG_ZERO_MASK_KNOWN;
                vc->flags |= QIR_FLAG_ZERO_MASK_KNOWN;
                continue;
            }
            case QIR_OP_GF2M_ADD: {
                const qir_val_t *a = &m->vals[op->ins[0]];
                const qir_val_t *b = &m->vals[op->ins[1]];
                for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                    zm[w] = (a->zero_mask[w] & b->zero_mask[w]) & mask.w[w];
                    cm[w] = (a->const_mask[w] & b->const_mask[w]) & mask.w[w];
                    cv[w] = (a->const_val[w] ^ b->const_val[w]) & cm[w];
                    zm[w] = (zm[w] | (cm[w] & ~cv[w])) & mask.w[w];
                }
                break;
            }
            case QIR_OP_LUT_ROM: {
                const qir_lut_payload_t *lut = &op->u.lut;
                uint64_t or_w[QIR_WORDS] = {0}, and_w[QIR_WORDS] = {0};
                bool full_dom = (lut->key_bits < 31u && lut->num_entries >= (UINT32_C(1) << lut->key_bits));
                if (full_dom && lut->num_entries > 0u) {
                    for (uint16_t w = 0u; w < lut->val_words && w < QIR_WORDS; ++w) and_w[w] = mask.w[w];
                }
                for (uint32_t e = 0u; e < lut->num_entries; ++e) {
                    for (uint16_t w = 0u; w < lut->val_words && w < QIR_WORDS; ++w) {
                        uint64_t v = lut->table[(size_t)e * lut->val_words + w] & mask.w[w];
                        or_w[w] |= v;
                        if (full_dom) and_w[w] &= v;
                    }
                }
                for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                    zm[w] = (~or_w[w]) & mask.w[w];
                    cm[w] = zm[w] | (and_w[w] & mask.w[w]);
                    cv[w] = and_w[w] & mask.w[w];
                }
                break;
            }
            case QIR_OP_SHL_VAR: {
                const qir_val_t *iv = &m->vals[op->ins[0]];
                uint16_t top_b = (uint16_t)(iv->width + op->u.shift.max_shift);
                for (uint16_t b = top_b; b < ov->width; ++b) {
                    qcc_bv_set(zm, b);
                    qcc_bv_set(cm, b);
                }
                break;
            }
            case QIR_OP_SHR_VAR: {
                const qir_val_t *iv = &m->vals[op->ins[0]];
                uint16_t hi_nz = iv->width;
                while (hi_nz > 0u && qcc_bv_get(iv->zero_mask, hi_nz - 1u)) hi_nz--;
                for (uint16_t b = hi_nz; b < ov->width; ++b) {
                    qcc_bv_set(zm, b);
                    qcc_bv_set(cm, b);
                }
                break;
            }
            case QIR_OP_GF2M_INV: {
                const qir_val_t *a = &m->vals[op->ins[0]];
                bool a_zero = true;
                for (uint16_t b = 0u; b < a->width; ++b) {
                    if (!qcc_bv_get(a->zero_mask, b)) { a_zero = false; break; }
                }
                bool ctrl_zero = (op->num_ins >= 2u && qcc_bv_get(m->vals[op->ins[1]].zero_mask, 0u));
                if (a_zero || ctrl_zero) {
                    for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                        zm[w] = mask.w[w];
                        cm[w] = mask.w[w];
                    }
                }
                break;
            }
            default:
                break;
        }
        if (op->pred_val != UINT32_MAX) {
            /* Inactive predicate lanes stay 0; known-0 stays sound, known-1 does not */
            for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                cm[w] = zm[w];
                cv[w] = 0u;
            }
        }
        memcpy(ov->zero_mask, zm, sizeof(zm));
        memcpy(ov->const_mask, cm, sizeof(cm));
        memcpy(ov->const_val, cv, sizeof(cv));
        ov->flags |= QIR_FLAG_ZERO_MASK_KNOWN;
    }
}

/* ============================================================================
 * Pure C17 64-Lane Semantic IR Evaluator (zero indirect function calls)
 * ============================================================================ */

static inline bool qir_eval_64(
    const qir_module_t *m,
    const qir_val_vec64_t *input_lanes, /* [m->num_inputs], indexed by input ord */
    uint64_t seed,
    qcc_arena_t *scratch,
    qir_val_vec64_t **out_val_vecs,     /* Optional [m->num_vals] allocated in scratch */
    flow_u128 out_phase[64]
) {
    if (!m || !scratch) return false;
    qir_val_vec64_t *vecs = (qir_val_vec64_t *)qcc_arena_alloc_zero(
        scratch, (size_t)m->num_vals * sizeof(qir_val_vec64_t), 16u);
    if (!vecs) return false;
    flow_u128 phase[64];
    memset(phase, 0, sizeof(phase));

    for (uint32_t i = 0u; i < m->num_inputs; ++i) {
        uint32_t vid = m->inputs[i];
        qcc_bv256_t mask = qcc_bv_mask_n(m->vals[vid].width);
        if (input_lanes) {
            for (uint32_t l = 0u; l < 64u; ++l) {
                for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                    vecs[vid].lanes[l][w] = input_lanes[i].lanes[l][w] & mask.w[w];
                }
            }
        }
    }

    uint64_t rng = seed ? seed : UINT64_C(0x123456789ABCDEF0);

    for (uint32_t oi = 0u; oi < m->num_ops; ++oi) {
        const qir_op_t *op = &m->ops[oi];
        uint64_t pred_mask = UINT64_MAX;
        if (op->pred_val != UINT32_MAX && op->pred_val < m->num_vals) {
            pred_mask = 0u;
            for (uint32_t l = 0u; l < 64u; ++l) {
                if (vecs[op->pred_val].lanes[l][0] & UINT64_C(1)) pred_mask |= (UINT64_C(1) << l);
            }
        }
        switch (op->kind) {
            case QIR_OP_CONST: {
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                qcc_bv256_t mask = qcc_bv_mask_n(m->vals[op->outs[0]].width);
                for (uint32_t l = 0u; l < 64u; ++l) {
                    for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                        dst->lanes[l][w] = op->u.cst.val[w] & mask.w[w];
                    }
                }
                break;
            }
            case QIR_OP_SLICE: {
                const qir_val_vec64_t *src = &vecs[op->ins[0]];
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                uint16_t lo = op->u.slice.lo;
                uint16_t ow = m->vals[op->outs[0]].width;
                qcc_bv256_t mask = qcc_bv_mask_n(ow);
                uint32_t ws = lo / 64u;
                uint32_t bs = lo % 64u;
                for (uint32_t l = 0u; l < 64u; ++l) {
                    for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                        uint64_t v0 = (w + ws < QIR_WORDS) ? src->lanes[l][w + ws] : 0u;
                        uint64_t v1 = (bs != 0u && w + ws + 1u < QIR_WORDS) ? src->lanes[l][w + ws + 1u] : 0u;
                        dst->lanes[l][w] = ((bs == 0u) ? v0 : ((v0 >> bs) | (v1 << (64u - bs)))) & mask.w[w];
                    }
                }
                break;
            }
            case QIR_OP_CONCAT: {
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                for (uint32_t l = 0u; l < 64u; ++l) {
                    memset(dst->lanes[l], 0, sizeof(dst->lanes[l]));
                    uint32_t pos = 0u;
                    for (uint16_t i = 0u; i < op->num_ins; ++i) {
                        uint32_t iv = op->ins[i];
                        uint16_t iw = m->vals[iv].width;
                        for (uint16_t b = 0u; b < iw && pos < QIR_MAX_VAL_BITS; ++b, ++pos) {
                            if (qcc_bv_get(vecs[iv].lanes[l], b)) qcc_bv_set(dst->lanes[l], pos);
                        }
                    }
                }
                break;
            }
            case QIR_OP_LINF2: {
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                uint16_t ow = m->vals[op->outs[0]].width;
                for (uint32_t l = 0u; l < 64u; ++l) {
                    uint64_t flat[QIR_WORDS] = {0};
                    uint32_t pos = 0u;
                    for (uint16_t i = 0u; i < op->num_ins; ++i) {
                        uint32_t iv = op->ins[i];
                        uint16_t iw = m->vals[iv].width;
                        for (uint16_t b = 0u; b < iw && pos < QIR_MAX_VAL_BITS; ++b, ++pos) {
                            if (qcc_bv_get(vecs[iv].lanes[l], b)) qcc_bv_set(flat, pos);
                        }
                    }
                    memset(dst->lanes[l], 0, sizeof(dst->lanes[l]));
                    for (uint16_t r = 0u; r < ow && r < op->u.linf2.out_bits; ++r) {
                        uint32_t p = 0u;
                        for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                            p ^= (uint32_t)__builtin_popcountll((unsigned long long)(op->u.linf2.rows[r].w[w] & flat[w]));
                        }
                        if (qcc_bv_get(op->u.linf2.const_vec.w, r)) p ^= 1u;
                        if (p & 1u) qcc_bv_set(dst->lanes[l], r);
                    }
                }
                break;
            }
            case QIR_OP_AND_VEC: {
                const qir_val_vec64_t *a = &vecs[op->ins[0]];
                const qir_val_vec64_t *b = &vecs[op->ins[1]];
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                qcc_bv256_t mask = qcc_bv_mask_n(m->vals[op->outs[0]].width);
                for (uint32_t l = 0u; l < 64u; ++l) {
                    for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                        dst->lanes[l][w] = a->lanes[l][w] & b->lanes[l][w] & mask.w[w];
                    }
                }
                break;
            }
            case QIR_OP_ADD: {
                const qir_val_vec64_t *a = &vecs[op->ins[0]];
                const qir_val_vec64_t *b = &vecs[op->ins[1]];
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                uint16_t ow = m->vals[op->outs[0]].width;
                for (uint32_t l = 0u; l < 64u; ++l) {
                    qir__u256_add(a->lanes[l], b->lanes[l], ow, dst->lanes[l]);
                }
                break;
            }
            case QIR_OP_SUB: {
                const qir_val_vec64_t *a = &vecs[op->ins[0]];
                const qir_val_vec64_t *b = &vecs[op->ins[1]];
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                uint16_t ow = m->vals[op->outs[0]].width;
                for (uint32_t l = 0u; l < 64u; ++l) {
                    qir__u256_sub(a->lanes[l], b->lanes[l], ow, dst->lanes[l]);
                }
                break;
            }
            case QIR_OP_CS_ADD3: {
                const qir_val_vec64_t *a = &vecs[op->ins[0]];
                const qir_val_vec64_t *b = &vecs[op->ins[1]];
                const qir_val_vec64_t *c = &vecs[op->ins[2]];
                qir_val_vec64_t *vs = &vecs[op->outs[0]];
                qir_val_vec64_t *vc = &vecs[op->outs[1]];
                qcc_bv256_t ms = qcc_bv_mask_n(m->vals[op->outs[0]].width);
                qcc_bv256_t mc = qcc_bv_mask_n(m->vals[op->outs[1]].width);
                for (uint32_t l = 0u; l < 64u; ++l) {
                    uint64_t carry_in = 0u;
                    for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                        uint64_t av = a->lanes[l][w], bv = b->lanes[l][w], cv = c->lanes[l][w];
                        uint64_t s = av ^ bv ^ cv;
                        uint64_t maj = (av & bv) | (bv & cv) | (av & cv);
                        vs->lanes[l][w] = s & ms.w[w];
                        vc->lanes[l][w] = ((maj << 1u) | carry_in) & mc.w[w];
                        carry_in = maj >> 63u;
                    }
                }
                break;
            }
            case QIR_OP_MUL_FXP: {
                const qir_val_vec64_t *a = &vecs[op->ins[0]];
                const qir_val_vec64_t *b = &vecs[op->ins[1]];
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                uint16_t ow = m->vals[op->outs[0]].width;
                int32_t shift_s = (int32_t)op->u.mul.frac_a + (int32_t)op->u.mul.frac_b - (int32_t)op->u.mul.frac_out;
                uint32_t shift = shift_s > 0 ? (uint32_t)shift_s : 0u;
                for (uint32_t l = 0u; l < 64u; ++l) {
                    qir__u256_mul_shift(a->lanes[l], b->lanes[l], shift, ow, dst->lanes[l]);
                }
                break;
            }
            case QIR_OP_SQR_FXP: {
                const qir_val_vec64_t *a = &vecs[op->ins[0]];
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                uint16_t ow = m->vals[op->outs[0]].width;
                int32_t shift_s = 2 * (int32_t)op->u.mul.frac_a - (int32_t)op->u.mul.frac_out;
                uint32_t shift = shift_s > 0 ? (uint32_t)shift_s : 0u;
                for (uint32_t l = 0u; l < 64u; ++l) {
                    qir__u256_mul_shift(a->lanes[l], a->lanes[l], shift, ow, dst->lanes[l]);
                }
                break;
            }
            case QIR_OP_GF2M_ADD: {
                const qir_val_vec64_t *a = &vecs[op->ins[0]];
                const qir_val_vec64_t *b = &vecs[op->ins[1]];
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                qcc_bv256_t mask = qcc_bv_mask_n(op->u.gf2m.m);
                for (uint32_t l = 0u; l < 64u; ++l) {
                    for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                        dst->lanes[l][w] = (a->lanes[l][w] ^ b->lanes[l][w]) & mask.w[w];
                    }
                }
                break;
            }
            case QIR_OP_GF2M_MUL: {
                const qir_val_vec64_t *a = &vecs[op->ins[0]];
                const qir_val_vec64_t *b = &vecs[op->ins[1]];
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                uint16_t deg = op->u.gf2m.m;
                uint32_t poly = op->u.gf2m.irred_poly;
                for (uint32_t l = 0u; l < 64u; ++l) {
                    dst->lanes[l][0] = qir__gf2m_mul_u32((uint32_t)a->lanes[l][0], (uint32_t)b->lanes[l][0], deg, poly);
                }
                break;
            }
            case QIR_OP_GF2M_INV: {
                const qir_val_vec64_t *a = &vecs[op->ins[0]];
                const qir_val_vec64_t *ctrl = (op->num_ins >= 2u) ? &vecs[op->ins[1]] : NULL;
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                uint16_t deg = op->u.gf2m.m;
                uint32_t poly = op->u.gf2m.irred_poly;
                for (uint32_t l = 0u; l < 64u; ++l) {
                    if (ctrl && ((ctrl->lanes[l][0] & 1u) == 0u)) {
                        dst->lanes[l][0] = 0u;
                    } else {
                        dst->lanes[l][0] = qir__gf2m_inv_u32((uint32_t)a->lanes[l][0], deg, poly);
                    }
                }
                break;
            }
            case QIR_OP_LUT_ROM: {
                const qir_val_vec64_t *k = &vecs[op->ins[0]];
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                uint16_t vw = op->u.lut.val_words;
                qcc_bv256_t mask = qcc_bv_mask_n(op->u.lut.val_bits);
                uint64_t kmask = (op->u.lut.key_bits >= 64u) ? UINT64_MAX : ((UINT64_C(1) << op->u.lut.key_bits) - 1u);
                for (uint32_t l = 0u; l < 64u; ++l) {
                    uint64_t addr = k->lanes[l][0] & kmask;
                    if (addr < op->u.lut.num_entries) {
                        for (uint16_t w = 0u; w < vw && w < QIR_WORDS; ++w) {
                            dst->lanes[l][w] = op->u.lut.table[(size_t)addr * vw + w] & mask.w[w];
                        }
                    }
                }
                break;
            }
            case QIR_OP_CLZ_NORM: {
                const qir_val_vec64_t *x = &vecs[op->ins[0]];
                qir_val_vec64_t *vm = &vecs[op->outs[0]];
                qir_val_vec64_t *vc = &vecs[op->outs[1]];
                uint16_t iw = m->vals[op->ins[0]].width;
                uint16_t mw = m->vals[op->outs[0]].width;
                for (uint32_t l = 0u; l < 64u; ++l) {
                    uint32_t lz = iw;
                    for (int32_t b = (int32_t)iw - 1; b >= 0; --b) {
                        if (qcc_bv_get(x->lanes[l], (uint32_t)b)) {
                            lz = (uint32_t)(iw - 1 - b);
                            break;
                        }
                    }
                    vc->lanes[l][0] = lz;
                    /* Shift left by lz, then extract top mw bits of [0..iw-1] */
                    uint64_t shifted[QIR_WORDS] = {0};
                    for (uint32_t b = 0u; b + lz < iw; ++b) {
                        if (qcc_bv_get(x->lanes[l], b)) qcc_bv_set(shifted, b + lz);
                    }
                    uint16_t drop = (iw > mw) ? (uint16_t)(iw - mw) : 0u;
                    for (uint16_t b = 0u; b < mw; ++b) {
                        if (qcc_bv_get(shifted, drop + b)) qcc_bv_set(vm->lanes[l], b);
                    }
                }
                break;
            }
            case QIR_OP_CTRL_SEL: {
                const qir_val_vec64_t *sel = &vecs[op->ins[0]];
                const qir_val_vec64_t *vf  = &vecs[op->ins[1]];
                const qir_val_vec64_t *vt  = &vecs[op->ins[2]];
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                qcc_bv256_t mask = qcc_bv_mask_n(m->vals[op->outs[0]].width);
                for (uint32_t l = 0u; l < 64u; ++l) {
                    uint64_t s = sel->lanes[l][0] & 1u;
                    for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                        dst->lanes[l][w] = (s ? vt->lanes[l][w] : vf->lanes[l][w]) & mask.w[w];
                    }
                }
                break;
            }
            case QIR_OP_PHASE_POLY: {
                for (uint32_t l = 0u; l < 64u; ++l) {
                    if (((pred_mask >> l) & 1u) == 0u) continue;
                    uint64_t flat[QIR_WORDS] = {0};
                    uint32_t pos = 0u;
                    for (uint16_t i = 0u; i < op->num_ins; ++i) {
                        uint32_t iv = op->ins[i];
                        uint16_t iw = m->vals[iv].width;
                        for (uint16_t b = 0u; b < iw && pos < QIR_MAX_VAL_BITS; ++b, ++pos) {
                            if (qcc_bv_get(vecs[iv].lanes[l], b)) qcc_bv_set(flat, pos);
                        }
                    }
                    for (uint32_t t = 0u; t < op->u.phase_poly.num_terms; ++t) {
                        const qir_phase_term_t *pt = &op->u.phase_poly.terms[t];
                        uint32_t fire = 0u;
                        if (pt->is_xor) {
                            uint32_t pop = 0u;
                            for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                                pop ^= (uint32_t)__builtin_popcountll((unsigned long long)(flat[w] & pt->mask[w]));
                            }
                            fire = pop & 1u;
                        } else {
                            fire = 1u;
                            for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                                if ((flat[w] & pt->mask[w]) != pt->mask[w]) {
                                    fire = 0u;
                                    break;
                                }
                            }
                        }
                        if (fire) phase[l] += pt->angle;
                    }
                }
                break;
            }
            case QIR_OP_OPAQUE_KMB: {
                const qcc_kmb_circuit_t *k = op->u.opaque.kmb;
                if (!k) break;
                uint32_t nq = k->num_qubits;
                uint32_t nb = k->num_bits;
                uint64_t *qstate = nq ? (uint64_t *)calloc(nq, sizeof(uint64_t)) : NULL;
                uint64_t *bstate = nb ? (uint64_t *)calloc(nb, sizeof(uint64_t)) : NULL;
                uint64_t *cond_stack = (uint64_t *)calloc(64u, sizeof(uint64_t));
                if ((nq && !qstate) || (nb && !bstate) || !cond_stack) {
                    free(qstate); free(bstate); free(cond_stack);
                    return false;
                }
                size_t cond_len = 0u, cond_cap = 64u;
                /* Flatten IR inputs into mapped qubits */
                uint32_t pos = 0u;
                for (uint16_t i = 0u; i < op->num_ins; ++i) {
                    uint32_t iv = op->ins[i];
                    uint16_t iw = m->vals[iv].width;
                    for (uint16_t b = 0u; b < iw && pos < op->u.opaque.in_bits; ++b, ++pos) {
                        uint32_t qid = op->u.opaque.in_q_idx[pos];
                        if (qid < nq) {
                            uint64_t lane_mask = 0u;
                            for (uint32_t l = 0u; l < 64u; ++l) {
                                if (qcc_bv_get(vecs[iv].lanes[l], b)) lane_mask |= (UINT64_C(1) << l);
                            }
                            qstate[qid] = lane_mask;
                        }
                    }
                }
                uint64_t cond_all = pred_mask;
                qcc_kmb_cursor_t cur;
                qcc_kmb_cursor_init(&cur);
                qcc_kmb_op_view_t v;
                while (qcc_kmb_cursor_next(k, &cur, &v)) {
                    if (v.op == QCC_OP_PUSH_CONDITION) {
                        if (cond_len >= cond_cap) {
                            cond_cap *= 2u;
                            uint64_t *nc = (uint64_t *)realloc(cond_stack, cond_cap * sizeof(uint64_t));
                            if (!nc) {
                                free(qstate); free(bstate); free(cond_stack);
                                return false;
                            }
                            cond_stack = nc;
                        }
                        cond_stack[cond_len++] = bstate[v.cond];
                        cond_all &= bstate[v.cond];
                        continue;
                    }
                    if (v.op == QCC_OP_POP_CONDITION) {
                        if (cond_len > 0u) cond_len--;
                        cond_all = pred_mask;
                        for (size_t s = 0u; s < cond_len; ++s) cond_all &= cond_stack[s];
                        continue;
                    }
                    uint64_t active = cond_all & (v.has_cond ? bstate[v.cond] : UINT64_MAX);
                    uint8_t base = v.op & (uint8_t)(~QCC_OP_FLAG_COND);
                    if (base == QCC_OP_NEG) {
                        for (uint32_t l = 0u; l < 64u; ++l) {
                            if ((active >> l) & 1u) phase[l] += ((flow_u128)1u << 127u);
                        }
                    } else if (base == QCC_OP_BIT_STORE0) {
                        bstate[v.b] &= ~active;
                    } else if (base == QCC_OP_BIT_STORE1) {
                        bstate[v.b] |= active;
                    } else if (base == QCC_OP_BIT_INVERT) {
                        bstate[v.b] ^= active;
                    } else if (base == QCC_OP_X) {
                        qstate[v.q[0]] ^= active;
                    } else if (base == QCC_OP_Z) {
                        uint64_t hit = active & qstate[v.q[0]];
                        for (uint32_t l = 0u; l < 64u; ++l) {
                            if ((hit >> l) & 1u) phase[l] += ((flow_u128)1u << 127u);
                        }
                    } else if (base == QCC_OP_R) {
                        uint64_t rbits = qcc_splitmix64(&rng);
                        uint64_t hit = active & qstate[v.q[0]] & rbits;
                        for (uint32_t l = 0u; l < 64u; ++l) {
                            if ((hit >> l) & 1u) phase[l] += ((flow_u128)1u << 127u);
                        }
                        qstate[v.q[0]] &= ~active;
                    } else if (base == QCC_OP_Z_POW) {
                        flow_u128 ang = qcc_u128_from_angle(v.angle);
                        uint64_t hit = active & qstate[v.q[0]];
                        for (uint32_t l = 0u; l < 64u; ++l) {
                            if ((hit >> l) & 1u) phase[l] += ang;
                        }
                    } else if (base == QCC_OP_HMR) {
                        uint64_t rbits = qcc_splitmix64(&rng);
                        uint64_t flip = active & rbits & qstate[v.q[0]];
                        for (uint32_t l = 0u; l < 64u; ++l) {
                            if ((flip >> l) & 1u) phase[l] += ((flow_u128)1u << 127u);
                        }
                        bstate[v.b] = (bstate[v.b] & ~active) | (rbits & active);
                        qstate[v.q[0]] &= ~active;
                    } else if (base == QCC_OP_CX) {
                        qstate[v.q[1]] ^= (active & qstate[v.q[0]]);
                    } else if (base == QCC_OP_CZ) {
                        uint64_t hit = active & qstate[v.q[0]] & qstate[v.q[1]];
                        for (uint32_t l = 0u; l < 64u; ++l) {
                            if ((hit >> l) & 1u) phase[l] += ((flow_u128)1u << 127u);
                        }
                    } else if (base == QCC_OP_SWAP) {
                        uint64_t diff = (qstate[v.q[0]] ^ qstate[v.q[1]]) & active;
                        qstate[v.q[0]] ^= diff;
                        qstate[v.q[1]] ^= diff;
                    } else if (base == QCC_OP_CCX) {
                        qstate[v.q[2]] ^= (active & qstate[v.q[0]] & qstate[v.q[1]]);
                    } else if (base == QCC_OP_CCZ) {
                        uint64_t hit = active & qstate[v.q[0]] & qstate[v.q[1]] & qstate[v.q[2]];
                        for (uint32_t l = 0u; l < 64u; ++l) {
                            if ((hit >> l) & 1u) phase[l] += ((flow_u128)1u << 127u);
                        }
                    } else if (base == QCC_OP_CSWAP) {
                        uint64_t act2 = active & qstate[v.q[0]];
                        uint64_t diff = (qstate[v.q[1]] ^ qstate[v.q[2]]) & act2;
                        qstate[v.q[1]] ^= diff;
                        qstate[v.q[2]] ^= diff;
                    }
                }
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                for (uint16_t b = 0u; b < op->u.opaque.out_bits && b < QIR_MAX_VAL_BITS; ++b) {
                    uint32_t qid = op->u.opaque.out_q_idx[b];
                    uint64_t msk = (qid < nq) ? qstate[qid] : 0u;
                    for (uint32_t l = 0u; l < 64u; ++l) {
                        if ((msk >> l) & 1u) qcc_bv_set(dst->lanes[l], b);
                    }
                }
                free(qstate);
                free(bstate);
                free(cond_stack);
                break;
            }
            case QIR_OP_HMR_MEAS: {
                const qir_val_vec64_t *src = &vecs[op->ins[0]];
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                uint16_t w = m->vals[op->outs[0]].width;
                for (uint32_t l = 0u; l < 64u; ++l) memset(dst->lanes[l], 0, sizeof(dst->lanes[l]));
                for (uint16_t b = 0u; b < w; ++b) {
                    uint64_t rbits = qcc_splitmix64(&rng);
                    for (uint32_t l = 0u; l < 64u; ++l) {
                        if (((pred_mask >> l) & 1u) == 0u) continue;
                        uint64_t rb = (rbits >> l) & 1u;
                        uint64_t qb = qcc_bv_get(src->lanes[l], b) ? 1u : 0u;
                        if (rb) qcc_bv_set(dst->lanes[l], b);
                        if (rb & qb) phase[l] += ((flow_u128)1u << 127u);
                    }
                }
                break;
            }
            case QIR_OP_BIT_OP: {
                const qir_val_vec64_t *a = &vecs[op->ins[0]];
                const qir_val_vec64_t *b = (op->num_ins > 1u) ? &vecs[op->ins[1]] : NULL;
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                qcc_bv256_t mask = qcc_bv_mask_n(m->vals[op->outs[0]].width);
                uint8_t mode = op->u.bit_op.mode;
                for (uint32_t l = 0u; l < 64u; ++l) {
                    if (mode == 3u) {
                        uint64_t any = 0u;
                        for (uint32_t w = 0u; w < QIR_WORDS; ++w) any |= a->lanes[l][w];
                        dst->lanes[l][0] = any ? 1u : 0u;
                    } else {
                        for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                            uint64_t av = a->lanes[l][w], bv = b ? b->lanes[l][w] : 0u;
                            uint64_t rv = (mode == 0u) ? (av ^ bv) : ((mode == 1u) ? (av & bv) : (~av));
                            dst->lanes[l][w] = rv & mask.w[w];
                        }
                    }
                }
                break;
            }
            case QIR_OP_SHL_VAR: {
                const qir_val_vec64_t *val = &vecs[op->ins[0]];
                const qir_val_vec64_t *sh  = &vecs[op->ins[1]];
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                uint16_t iw = m->vals[op->ins[0]].width;
                uint16_t ow = m->vals[op->outs[0]].width;
                for (uint32_t l = 0u; l < 64u; ++l) {
                    memset(dst->lanes[l], 0, sizeof(dst->lanes[l]));
                    uint64_t amt = sh->lanes[l][0];
                    for (uint16_t b = 0u; b < iw; ++b) {
                        if (qcc_bv_get(val->lanes[l], b) && amt + b < ow) {
                            qcc_bv_set(dst->lanes[l], (uint32_t)(amt + b));
                        }
                    }
                }
                break;
            }
            case QIR_OP_SHR_VAR: {
                const qir_val_vec64_t *val = &vecs[op->ins[0]];
                const qir_val_vec64_t *sh  = &vecs[op->ins[1]];
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                uint16_t iw = m->vals[op->ins[0]].width;
                uint16_t ow = m->vals[op->outs[0]].width;
                for (uint32_t l = 0u; l < 64u; ++l) {
                    memset(dst->lanes[l], 0, sizeof(dst->lanes[l]));
                    uint64_t amt = sh->lanes[l][0];
                    for (uint16_t b = 0u; b < ow; ++b) {
                        if (amt + b < iw && qcc_bv_get(val->lanes[l], (uint32_t)(amt + b))) {
                            qcc_bv_set(dst->lanes[l], b);
                        }
                    }
                }
                break;
            }
            case QIR_OP_EQ_ZERO: {
                const qir_val_vec64_t *a = &vecs[op->ins[0]];
                qir_val_vec64_t *dst = &vecs[op->outs[0]];
                for (uint32_t l = 0u; l < 64u; ++l) {
                    uint64_t any = 0u;
                    for (uint32_t w = 0u; w < QIR_WORDS; ++w) any |= a->lanes[l][w];
                    dst->lanes[l][0] = any ? 0u : 1u;
                }
                break;
            }
            case QIR_OP_PHASE_DOT: {
                for (uint32_t l = 0u; l < 64u; ++l) {
                    if (((pred_mask >> l) & 1u) == 0u) continue;
                    uint32_t par = 0u;
                    if (op->num_ins == 2u) {
                        for (uint32_t w = 0u; w < QIR_WORDS; ++w) {
                            par ^= (uint32_t)__builtin_popcountll((unsigned long long)(vecs[op->ins[0]].lanes[l][w] & vecs[op->ins[1]].lanes[l][w]));
                        }
                    } else if (op->num_ins >= 3u) {
                        uint32_t p = qir__gf2m_mul_u32((uint32_t)vecs[op->ins[0]].lanes[l][0],
                                                       (uint32_t)vecs[op->ins[1]].lanes[l][0],
                                                       op->u.gf2m.m, op->u.gf2m.irred_poly);
                        par = (uint32_t)__builtin_popcount((unsigned)(p & (uint32_t)vecs[op->ins[2]].lanes[l][0]));
                    }
                    if (par & 1u) phase[l] += ((flow_u128)1u << 127u);
                }
                break;
            }
            default:
                return false;
        }
        if (pred_mask != UINT64_MAX) {
            for (uint8_t ko = 0u; ko < op->num_outs; ++ko) {
                uint32_t ovid = op->outs[ko];
                for (uint32_t l = 0u; l < 64u; ++l) {
                    if (((pred_mask >> l) & 1u) == 0u) {
                        memset(vecs[ovid].lanes[l], 0, sizeof(vecs[ovid].lanes[l]));
                    }
                }
            }
        }
    }

    if (out_val_vecs) *out_val_vecs = vecs;
    if (out_phase) memcpy(out_phase, phase, sizeof(phase));
    return true;
}

#ifdef __cplusplus
}
#endif

#endif /* QIR_H */
