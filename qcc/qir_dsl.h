/*
 * qir_dsl.h — S-Expression (.qir) Text DSL Printer & Parser
 * Part of QCC (Quantum Circuit Compiler) for FLOW v7.0.1.
 * Pure C17, zero external dependencies.
 */
#ifndef QIR_DSL_H
#define QIR_DSL_H

#include "qir.h"
#include <ctype.h>

#ifdef __cplusplus
extern "C" {
#endif

static inline void qir__fprint_bv_hex(FILE *f, const uint64_t w[QIR_WORDS], uint16_t nbits) {
    uint32_t top = nbits ? (uint32_t)((nbits - 1u) / 64u) : 0u;
    if (top >= QIR_WORDS) top = QIR_WORDS - 1u;
    while (top > 0u && w[top] == 0u) top--;
    fprintf(f, "0x%llx", (unsigned long long)w[top]);
    for (int32_t i = (int32_t)top - 1; i >= 0; --i) {
        fprintf(f, "%016llx", (unsigned long long)w[i]);
    }
}

static inline bool qir_print_module(const qir_module_t *m, FILE *f) {
    if (!m || !f) return false;
    fprintf(f, "(qir \"%s\"\n", m->name ? m->name : "qir_mod");
    if (m->require_clean_scratch) fprintf(f, "  (clean_anc)\n");

    for (uint32_t i = 0u; i < m->num_inputs; ++i) {
        uint32_t vid = m->inputs[i];
        const qir_val_t *v = &m->vals[vid];
        char pfx = (v->rep == QIR_REP_CLASSICAL) ? 'c' : 'u';
        fprintf(f, "  (input %%%s %c%u)\n", v->name ? v->name : "in", pfx, (unsigned)v->width);
    }

    for (uint32_t oi = 0u; oi < m->num_ops; ++oi) {
        const qir_op_t *op = &m->ops[oi];
        fprintf(f, "  (let ");
        if (op->num_outs == 1u) {
            const qir_val_t *ov = &m->vals[op->outs[0]];
            char pfx = (ov->rep == QIR_REP_CLASSICAL) ? 'c' : 'u';
            if (ov->name) fprintf(f, "%%%s %c%u ", ov->name, pfx, (unsigned)ov->width);
            else fprintf(f, "%%v%u %c%u ", ov->id, pfx, (unsigned)ov->width);
        } else {
            fprintf(f, "(");
            for (uint8_t k = 0u; k < op->num_outs; ++k) {
                const qir_val_t *ov = &m->vals[op->outs[k]];
                char pfx = (ov->rep == QIR_REP_CLASSICAL) ? 'c' : 'u';
                if (k > 0u) fprintf(f, " ");
                if (ov->name) fprintf(f, "%%%s %c%u", ov->name, pfx, (unsigned)ov->width);
                else fprintf(f, "%%v%u %c%u", ov->id, pfx, (unsigned)ov->width);
            }
            fprintf(f, ") ");
        }

        const char *kname = qir_op_kind_name(op->kind);
        if (op->kind == QIR_OP_BIT_OP) {
            kname = (op->u.bit_op.mode == 0u) ? "bit_xor" : ((op->u.bit_op.mode == 1u) ? "bit_and" : "bit_not");
        }
        fprintf(f, "(%s", kname);

        if (op->kind == QIR_OP_CONST) {
            fprintf(f, " ");
            qir__fprint_bv_hex(f, op->u.cst.val, m->vals[op->outs[0]].width);
        } else if (op->kind == QIR_OP_LINF2) {
            fprintf(f, " (");
            for (uint16_t i = 0u; i < op->num_ins; ++i) {
                const qir_val_t *iv = &m->vals[op->ins[i]];
                if (i > 0u) fprintf(f, " ");
                if (iv->name) fprintf(f, "%%%s", iv->name);
                else fprintf(f, "%%v%u", iv->id);
            }
            fprintf(f, ") :c ");
            qir__fprint_bv_hex(f, op->u.linf2.const_vec.w, op->u.linf2.out_bits);
            fprintf(f, " :rows (");
            for (uint16_t r = 0u; r < op->u.linf2.out_bits; ++r) {
                if (r > 0u) fprintf(f, " ");
                qir__fprint_bv_hex(f, op->u.linf2.rows[r].w, op->u.linf2.in_bits);
            }
            fprintf(f, ")");
        } else if (op->kind == QIR_OP_PHASE_POLY) {
            fprintf(f, " (");
            for (uint16_t i = 0u; i < op->num_ins; ++i) {
                const qir_val_t *iv = &m->vals[op->ins[i]];
                if (i > 0u) fprintf(f, " ");
                if (iv->name) fprintf(f, "%%%s", iv->name);
                else fprintf(f, "%%v%u", iv->id);
            }
            fprintf(f, ")");
            for (uint32_t t = 0u; t < op->u.phase_poly.num_terms; ++t) {
                const qir_phase_term_t *pt = &op->u.phase_poly.terms[t];
                qcc_angle128_t a = qcc_angle_from_u128(pt->angle);
                fprintf(f, " (%s ", pt->is_xor ? "xor" : "and");
                qir__fprint_bv_hex(f, pt->mask, op->u.phase_poly.total_in_bits);
                fprintf(f, " 0x%016llx 0x%016llx)", (unsigned long long)a.hi, (unsigned long long)a.lo);
            }
        } else {
            for (uint16_t i = 0u; i < op->num_ins; ++i) {
                const qir_val_t *iv = &m->vals[op->ins[i]];
                if (iv->name) fprintf(f, " %%%s", iv->name);
                else fprintf(f, " %%v%u", iv->id);
            }
            if (op->kind == QIR_OP_SLICE) {
                fprintf(f, " %u %u", (unsigned)op->u.slice.lo, (unsigned)op->u.slice.hi);
            } else if (op->kind == QIR_OP_MUL_FXP) {
                fprintf(f, " :frac %u %u %u", (unsigned)op->u.mul.frac_a,
                        (unsigned)op->u.mul.frac_b, (unsigned)op->u.mul.frac_out);
            } else if (op->kind == QIR_OP_SQR_FXP) {
                fprintf(f, " :frac %u %u", (unsigned)op->u.mul.frac_a, (unsigned)op->u.mul.frac_out);
            } else if (op->kind == QIR_OP_GF2M_ADD) {
                fprintf(f, " %u", (unsigned)op->u.gf2m.m);
            } else if (op->kind == QIR_OP_GF2M_MUL || op->kind == QIR_OP_GF2M_INV) {
                fprintf(f, " %u 0x%x", (unsigned)op->u.gf2m.m, (unsigned)op->u.gf2m.irred_poly);
            } else if (op->kind == QIR_OP_CLZ_NORM) {
                fprintf(f, " %u", (unsigned)op->u.clz.mant_bits);
            } else if (op->kind == QIR_OP_LUT_ROM) {
                fprintf(f, " (");
                for (uint32_t e = 0u; e < op->u.lut.num_entries; ++e) {
                    if (e > 0u) fprintf(f, " ");
                    const uint64_t *ew = &op->u.lut.table[(size_t)e * op->u.lut.val_words];
                    uint64_t tmp[QIR_WORDS] = {0};
                    for (uint16_t w = 0u; w < op->u.lut.val_words && w < QIR_WORDS; ++w) tmp[w] = ew[w];
                    qir__fprint_bv_hex(f, tmp, op->u.lut.val_bits);
                }
                fprintf(f, ")");
            } else if (op->kind == QIR_OP_OPAQUE_KMB && op->u.opaque.kmb) {
                fprintf(f, " :ops %zu :q %u", op->u.opaque.kmb->num_ops, (unsigned)op->u.opaque.num_q);
            }
        }

        if (op->flags & QIR_FLAG_APPROX) fprintf(f, " :approx %u", (unsigned)op->approx_tol_q16);
        if (op->flags & QIR_FLAG_CLEAN_ANC) fprintf(f, " :clean_anc");
        if (op->pred_val != UINT32_MAX && op->pred_val < m->num_vals) {
            const qir_val_t *pv = &m->vals[op->pred_val];
            if (pv->name) fprintf(f, " :if %%%s", pv->name);
            else fprintf(f, " :if %%v%u", pv->id);
        }
        fprintf(f, "))\n");
    }

    for (uint32_t i = 0u; i < m->num_outputs; ++i) {
        const qir_val_t *ov = &m->vals[m->outputs[i]];
        if (ov->name) fprintf(f, "  (output %%%s)\n", ov->name);
        else fprintf(f, "  (output %%v%u)\n", ov->id);
    }
    for (uint32_t i = 0u; i < m->num_classical_outputs; ++i) {
        const qir_val_t *ov = &m->vals[m->classical_outputs[i]];
        if (ov->name) fprintf(f, "  (output_bits %%%s)\n", ov->name);
        else fprintf(f, "  (output_bits %%v%u)\n", ov->id);
    }
    fprintf(f, ")\n");
    return true;
}

static inline bool qir_write_file(const qir_module_t *m, const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) return false;
    bool ok = qir_print_module(m, f);
    fclose(f);
    return ok;
}

/* ============================================================================
 * S-Expression Parser
 * ============================================================================ */

typedef struct {
    const char *src;
    size_t      pos;
    size_t      len;
} qir__lex_t;

typedef struct {
    const char **names;
    uint32_t    *ids;
    uint32_t     count;
    uint32_t     cap;
} qir__symtab_t;

static inline void qir__lex_skip(qir__lex_t *L) {
    while (L->pos < L->len) {
        char c = L->src[L->pos];
        if (isspace((unsigned char)c)) {
            L->pos++;
        } else if (c == ';') {
            while (L->pos < L->len && L->src[L->pos] != '\n') L->pos++;
        } else {
            break;
        }
    }
}

static inline int qir__lex_peek(qir__lex_t *L) {
    qir__lex_skip(L);
    return (L->pos < L->len) ? (unsigned char)L->src[L->pos] : -1;
}

static inline bool qir__lex_expect(qir__lex_t *L, char ch) {
    qir__lex_skip(L);
    if (L->pos >= L->len || L->src[L->pos] != ch) return false;
    L->pos++;
    return true;
}

static inline bool qir__lex_token(qir__lex_t *L, char *buf, size_t cap) {
    qir__lex_skip(L);
    if (L->pos >= L->len) return false;
    char c = L->src[L->pos];
    if (c == '(' || c == ')') return false;
    size_t n = 0u;
    if (c == '"') {
        L->pos++;
        while (L->pos < L->len && L->src[L->pos] != '"') {
            if (n + 1u < cap) buf[n++] = L->src[L->pos];
            L->pos++;
        }
        if (L->pos < L->len && L->src[L->pos] == '"') L->pos++;
        buf[n] = '\0';
        return true;
    }
    while (L->pos < L->len) {
        c = L->src[L->pos];
        if (isspace((unsigned char)c) || c == '(' || c == ')' || c == ';') break;
        if (n + 1u < cap) buf[n++] = c;
        L->pos++;
    }
    buf[n] = '\0';
    return n > 0u;
}

static inline bool qir__parse_bv256(const char *s, uint64_t out[QIR_WORDS]) {
    memset(out, 0, sizeof(uint64_t) * QIR_WORDS);
    if (!s || !*s) return false;
    if (s[0] == '0' && (s[1] == 'x' || s[1] == 'X')) {
        const char *h = s + 2;
        size_t len = strlen(h);
        if (len == 0u) return false;
        size_t pos = 0u;
        for (int32_t i = (int32_t)len - 1; i >= 0; --i) {
            char c = h[i];
            uint64_t nib = 0u;
            if (c >= '0' && c <= '9') nib = (uint64_t)(c - '0');
            else if (c >= 'a' && c <= 'f') nib = (uint64_t)(c - 'a' + 10);
            else if (c >= 'A' && c <= 'F') nib = (uint64_t)(c - 'A' + 10);
            else return false;
            size_t bit = pos * 4u;
            if (bit < 256u) out[bit >> 6u] |= (nib << (bit & 63u));
            pos++;
        }
        return true;
    }
    char *end = NULL;
    unsigned long long v = strtoull(s, &end, 10);
    if (!end || *end != '\0') return false;
    out[0] = (uint64_t)v;
    return true;
}

static inline bool qir__parse_width(const char *tok, uint16_t *out_w) {
    if (!tok || (tok[0] != 'u' && tok[0] != 'c')) return false;
    char *end = NULL;
    unsigned long w = strtoul(tok + 1, &end, 10);
    if (!end || *end != '\0' || w > QIR_MAX_VAL_BITS) return false;
    *out_w = (uint16_t)w;
    return true;
}

static inline void qir__sym_put(qir_module_t *m, qir__symtab_t *st, const char *tok, uint32_t id) {
    const char *name = (tok[0] == '%') ? (tok + 1) : tok;
    if (st->count >= st->cap) {
        uint32_t nc = st->cap ? st->cap * 2u : 64u;
        st->names = (const char **)realloc(st->names, nc * sizeof(const char *));
        st->ids   = (uint32_t *)realloc(st->ids, nc * sizeof(uint32_t));
        st->cap   = nc;
    }
    st->names[st->count] = qcc_arena_strdup(m->arena, name);
    st->ids[st->count]   = id;
    st->count++;
}

static inline bool qir__sym_get(const qir__symtab_t *st, const char *tok, uint32_t *out_id) {
    const char *name = (tok[0] == '%') ? (tok + 1) : tok;
    for (int32_t i = (int32_t)st->count - 1; i >= 0; --i) {
        if (strcmp(st->names[i], name) == 0) {
            *out_id = st->ids[i];
            return true;
        }
    }
    return false;
}

static inline bool qir_parse_str(qir_module_t *m, qcc_arena_t *arena, const char *src, size_t len) {
    qir_module_init(m, arena, "parsed");
    qir__lex_t L = {src, 0u, len};
    qir__symtab_t sym = {0};
    char tok[256];

    if (!qir__lex_expect(&L, '(')) goto fail;
    if (!qir__lex_token(&L, tok, sizeof(tok)) || strcmp(tok, "qir") != 0) goto fail;
    if (qir__lex_peek(&L) == '"') {
        if (!qir__lex_token(&L, tok, sizeof(tok))) goto fail;
        m->name = qcc_arena_strdup(arena, tok);
    }

    while (qir__lex_peek(&L) == '(') {
        if (!qir__lex_expect(&L, '(')) goto fail;
        if (!qir__lex_token(&L, tok, sizeof(tok))) goto fail;

        if (strcmp(tok, "clean_anc") == 0) {
            m->require_clean_scratch = true;
            if (!qir__lex_expect(&L, ')')) goto fail;
        } else if (strcmp(tok, "input") == 0) {
            char vtok[128], wtok[32];
            if (!qir__lex_token(&L, vtok, sizeof(vtok)) || !qir__lex_token(&L, wtok, sizeof(wtok))) goto fail;
            uint16_t w = 0u;
            if (!qir__parse_width(wtok, &w)) goto fail;
            const char *clean = (vtok[0] == '%') ? (vtok + 1) : vtok;
            uint32_t vid = qir_add_input(m, w, clean);
            qir__sym_put(m, &sym, vtok, vid);
            if (!qir__lex_expect(&L, ')')) goto fail;
        } else if (strcmp(tok, "output") == 0) {
            char vtok[128];
            if (!qir__lex_token(&L, vtok, sizeof(vtok))) goto fail;
            uint32_t vid = 0u;
            if (!qir__sym_get(&sym, vtok, &vid)) goto fail;
            qir_add_output(m, vid);
            if (!qir__lex_expect(&L, ')')) goto fail;
        } else if (strcmp(tok, "output_bits") == 0) {
            char vtok[128];
            if (!qir__lex_token(&L, vtok, sizeof(vtok))) goto fail;
            uint32_t vid = 0u;
            if (!qir__sym_get(&sym, vtok, &vid)) goto fail;
            qir_add_classical_output(m, vid);
            if (!qir__lex_expect(&L, ')')) goto fail;
        } else if (strcmp(tok, "let") == 0) {
            char out_names[4][128];
            uint16_t out_widths[4] = {0};
            uint8_t n_outs = 0u;
            if (qir__lex_peek(&L) == '(') {
                if (!qir__lex_expect(&L, '(')) goto fail;
                while (qir__lex_peek(&L) != ')') {
                    char wtok[32];
                    if (n_outs >= 4u) goto fail;
                    if (!qir__lex_token(&L, out_names[n_outs], 128u) || !qir__lex_token(&L, wtok, sizeof(wtok))) goto fail;
                    if (!qir__parse_width(wtok, &out_widths[n_outs])) goto fail;
                    n_outs++;
                }
                if (!qir__lex_expect(&L, ')')) goto fail;
            } else {
                char wtok[32];
                if (!qir__lex_token(&L, out_names[0], 128u) || !qir__lex_token(&L, wtok, sizeof(wtok))) goto fail;
                if (!qir__parse_width(wtok, &out_widths[0])) goto fail;
                n_outs = 1u;
            }

            if (!qir__lex_expect(&L, '(')) goto fail;
            char opname[64];
            if (!qir__lex_token(&L, opname, sizeof(opname))) goto fail;

            uint32_t out_ids[4] = {0};
            uint32_t flags = 0u;
            uint32_t approx_tol = 0u;

            if (strcmp(opname, "const") == 0) {
                uint64_t val[QIR_WORDS];
                if (!qir__lex_token(&L, tok, sizeof(tok)) || !qir__parse_bv256(tok, val)) goto fail;
                out_ids[0] = qir_emit_const(m, out_widths[0], val, out_names[0] + (out_names[0][0] == '%'));
            } else if (strcmp(opname, "slice") == 0) {
                char vtok[128], lo_s[32], hi_s[32];
                uint32_t iv = 0u;
                if (!qir__lex_token(&L, vtok, sizeof(vtok)) || !qir__sym_get(&sym, vtok, &iv)) goto fail;
                if (!qir__lex_token(&L, lo_s, sizeof(lo_s)) || !qir__lex_token(&L, hi_s, sizeof(hi_s))) goto fail;
                out_ids[0] = qir_emit_slice(m, iv, (uint16_t)strtoul(lo_s, NULL, 10), (uint16_t)strtoul(hi_s, NULL, 10),
                                            out_names[0] + (out_names[0][0] == '%'));
            } else if (strcmp(opname, "concat") == 0) {
                uint32_t ins[64];
                uint16_t ni = 0u;
                while (qir__lex_peek(&L) != ')' && qir__lex_peek(&L) != ':') {
                    if (!qir__lex_token(&L, tok, sizeof(tok))) goto fail;
                    if (!qir__sym_get(&sym, tok, &ins[ni])) goto fail;
                    ni++;
                }
                out_ids[0] = qir_emit_concat(m, ins, ni, out_names[0] + (out_names[0][0] == '%'));
            } else if (strcmp(opname, "and_vec") == 0 || strcmp(opname, "add") == 0 || strcmp(opname, "sub") == 0) {
                char a_t[128], b_t[128];
                uint32_t ia = 0u, ib = 0u;
                if (!qir__lex_token(&L, a_t, sizeof(a_t)) || !qir__sym_get(&sym, a_t, &ia)) goto fail;
                if (!qir__lex_token(&L, b_t, sizeof(b_t)) || !qir__sym_get(&sym, b_t, &ib)) goto fail;
                const char *nm = out_names[0] + (out_names[0][0] == '%');
                if (strcmp(opname, "and_vec") == 0) out_ids[0] = qir_emit_and_vec(m, ia, ib, 0u, nm);
                else if (strcmp(opname, "add") == 0) out_ids[0] = qir_emit_add(m, ia, ib, 0u, nm);
                else out_ids[0] = qir_emit_sub(m, ia, ib, 0u, nm);
                m->vals[out_ids[0]].width = out_widths[0];
            } else if (strcmp(opname, "cs_add3") == 0) {
                char a_t[128], b_t[128], c_t[128];
                uint32_t ia = 0u, ib = 0u, ic = 0u;
                if (!qir__lex_token(&L, a_t, sizeof(a_t)) || !qir__sym_get(&sym, a_t, &ia)) goto fail;
                if (!qir__lex_token(&L, b_t, sizeof(b_t)) || !qir__sym_get(&sym, b_t, &ib)) goto fail;
                if (!qir__lex_token(&L, c_t, sizeof(c_t)) || !qir__sym_get(&sym, c_t, &ic)) goto fail;
                if (n_outs != 2u) goto fail;
                qir_emit_cs_add3(m, ia, ib, ic, &out_ids[0], &out_ids[1],
                                 out_names[0] + (out_names[0][0] == '%'),
                                 out_names[1] + (out_names[1][0] == '%'));
            } else if (strcmp(opname, "ctrl_sel") == 0) {
                char s_t[128], f_t[128], t_t[128];
                uint32_t is = 0u, ifv = 0u, itv = 0u;
                if (!qir__lex_token(&L, s_t, sizeof(s_t)) || !qir__sym_get(&sym, s_t, &is)) goto fail;
                if (!qir__lex_token(&L, f_t, sizeof(f_t)) || !qir__sym_get(&sym, f_t, &ifv)) goto fail;
                if (!qir__lex_token(&L, t_t, sizeof(t_t)) || !qir__sym_get(&sym, t_t, &itv)) goto fail;
                out_ids[0] = qir_emit_ctrl_sel(m, is, ifv, itv, 0u, out_names[0] + (out_names[0][0] == '%'));
                m->vals[out_ids[0]].width = out_widths[0];
            } else if (strcmp(opname, "mul_fxp") == 0 || strcmp(opname, "sqr_fxp") == 0) {
                bool is_sqr = (strcmp(opname, "sqr_fxp") == 0);
                char a_t[128], b_t[128];
                uint32_t ia = 0u, ib = 0u;
                if (!qir__lex_token(&L, a_t, sizeof(a_t)) || !qir__sym_get(&sym, a_t, &ia)) goto fail;
                if (!is_sqr) {
                    if (!qir__lex_token(&L, b_t, sizeof(b_t)) || !qir__sym_get(&sym, b_t, &ib)) goto fail;
                }
                uint16_t fa = 0u, fb = 0u, fo = 0u;
                if (qir__lex_peek(&L) == ':') {
                    /* consume :frac if present */
                }
                const char *nm = out_names[0] + (out_names[0][0] == '%');
                out_ids[0] = is_sqr
                    ? qir_emit_sqr_fxp(m, ia, out_widths[0], fa, fo, 0u, nm)
                    : qir_emit_mul_fxp(m, ia, ib, out_widths[0], fa, fb, fo, 0u, nm);
            } else if (strcmp(opname, "gf2m_add") == 0) {
                char a_t[128], b_t[128], deg_s[32];
                uint32_t ia = 0u, ib = 0u;
                if (!qir__lex_token(&L, a_t, sizeof(a_t)) || !qir__sym_get(&sym, a_t, &ia)) goto fail;
                if (!qir__lex_token(&L, b_t, sizeof(b_t)) || !qir__sym_get(&sym, b_t, &ib)) goto fail;
                if (!qir__lex_token(&L, deg_s, sizeof(deg_s))) goto fail;
                out_ids[0] = qir_emit_gf2m_add(m, ia, ib, (uint16_t)strtoul(deg_s, NULL, 10),
                                               out_names[0] + (out_names[0][0] == '%'));
            } else if (strcmp(opname, "gf2m_mul") == 0) {
                char a_t[128], b_t[128], deg_s[32], p_s[64];
                uint32_t ia = 0u, ib = 0u;
                if (!qir__lex_token(&L, a_t, sizeof(a_t)) || !qir__sym_get(&sym, a_t, &ia)) goto fail;
                if (!qir__lex_token(&L, b_t, sizeof(b_t)) || !qir__sym_get(&sym, b_t, &ib)) goto fail;
                if (!qir__lex_token(&L, deg_s, sizeof(deg_s)) || !qir__lex_token(&L, p_s, sizeof(p_s))) goto fail;
                out_ids[0] = qir_emit_gf2m_mul(m, ia, ib, (uint16_t)strtoul(deg_s, NULL, 10),
                                               (uint32_t)strtoul(p_s, NULL, 0),
                                               out_names[0] + (out_names[0][0] == '%'));
            } else if (strcmp(opname, "gf2m_inv") == 0) {
                char a_t[128], deg_s[32], p_s[64];
                uint32_t ia = 0u;
                if (!qir__lex_token(&L, a_t, sizeof(a_t)) || !qir__sym_get(&sym, a_t, &ia)) goto fail;
                if (!qir__lex_token(&L, deg_s, sizeof(deg_s)) || !qir__lex_token(&L, p_s, sizeof(p_s))) goto fail;
                out_ids[0] = qir_emit_gf2m_inv(m, ia, (uint16_t)strtoul(deg_s, NULL, 10),
                                               (uint32_t)strtoul(p_s, NULL, 0),
                                               out_names[0] + (out_names[0][0] == '%'));
            } else if (strcmp(opname, "clz_norm") == 0) {
                char a_t[128], mb_s[32];
                uint32_t ia = 0u;
                if (!qir__lex_token(&L, a_t, sizeof(a_t)) || !qir__sym_get(&sym, a_t, &ia)) goto fail;
                if (!qir__lex_token(&L, mb_s, sizeof(mb_s))) goto fail;
                if (n_outs != 2u) goto fail;
                qir_emit_clz_norm(m, ia, (uint16_t)strtoul(mb_s, NULL, 10),
                                  &out_ids[0], &out_ids[1],
                                  out_names[0] + (out_names[0][0] == '%'),
                                  out_names[1] + (out_names[1][0] == '%'));
            } else if (strcmp(opname, "lut_rom") == 0) {
                char k_t[128];
                uint32_t ik = 0u;
                if (!qir__lex_token(&L, k_t, sizeof(k_t)) || !qir__sym_get(&sym, k_t, &ik)) goto fail;
                if (!qir__lex_expect(&L, '(')) goto fail;
                uint16_t vw = (uint16_t)((out_widths[0] + 63u) / 64u);
                if (vw == 0u) vw = 1u;
                size_t cap = 64u, cnt = 0u;
                uint64_t *tbl = (uint64_t *)calloc(cap * vw, sizeof(uint64_t));
                while (qir__lex_peek(&L) != ')') {
                    uint64_t v[QIR_WORDS];
                    if (!qir__lex_token(&L, tok, sizeof(tok)) || !qir__parse_bv256(tok, v)) {
                        free(tbl);
                        goto fail;
                    }
                    if (cnt >= cap) {
                        cap *= 2u;
                        tbl = (uint64_t *)realloc(tbl, cap * vw * sizeof(uint64_t));
                    }
                    for (uint16_t w = 0u; w < vw && w < QIR_WORDS; ++w) tbl[cnt * vw + w] = v[w];
                    cnt++;
                }
                if (!qir__lex_expect(&L, ')')) {
                    free(tbl);
                    goto fail;
                }
                out_ids[0] = qir_emit_lut_rom(m, ik, out_widths[0], (uint32_t)cnt, tbl, 0u,
                                              out_names[0] + (out_names[0][0] == '%'));
                free(tbl);
            } else if (strcmp(opname, "linf2") == 0) {
                uint32_t ins[64];
                uint16_t ni = 0u;
                if (!qir__lex_expect(&L, '(')) goto fail;
                while (qir__lex_peek(&L) != ')') {
                    if (!qir__lex_token(&L, tok, sizeof(tok)) || !qir__sym_get(&sym, tok, &ins[ni])) goto fail;
                    ni++;
                }
                if (!qir__lex_expect(&L, ')')) goto fail;
                qcc_bv256_t cvec = qcc_bv_zero();
                qcc_bv256_t *rows = (qcc_bv256_t *)qcc_arena_alloc_zero(arena, (size_t)out_widths[0] * sizeof(qcc_bv256_t), 8u);
                uint32_t linf2_pred = UINT32_MAX;
                while (qir__lex_peek(&L) == ':') {
                    if (!qir__lex_token(&L, tok, sizeof(tok))) goto fail;
                    if (strcmp(tok, ":c") == 0) {
                        if (!qir__lex_token(&L, tok, sizeof(tok)) || !qir__parse_bv256(tok, cvec.w)) goto fail;
                    } else if (strcmp(tok, ":rows") == 0) {
                        if (!qir__lex_expect(&L, '(')) goto fail;
                        uint16_t r = 0u;
                        while (qir__lex_peek(&L) != ')') {
                            if (!qir__lex_token(&L, tok, sizeof(tok))) goto fail;
                            if (r < out_widths[0]) {
                                if (!qir__parse_bv256(tok, rows[r].w)) goto fail;
                            }
                            r++;
                        }
                        if (!qir__lex_expect(&L, ')')) goto fail;
                    } else if (strcmp(tok, ":approx") == 0) {
                        if (!qir__lex_token(&L, tok, sizeof(tok))) goto fail;
                        flags |= QIR_FLAG_APPROX;
                        approx_tol = (uint32_t)strtoul(tok, NULL, 0);
                    } else if (strcmp(tok, ":if") == 0) {
                        char pt[128];
                        uint32_t pid = 0u;
                        if (!qir__lex_token(&L, pt, sizeof(pt)) || !qir__sym_get(&sym, pt, &pid)) goto fail;
                        linf2_pred = pid;
                    } else {
                        goto fail;
                    }
                }
                out_ids[0] = qir_emit_linf2(m, ins, ni, out_widths[0], rows, &cvec,
                                            out_names[0] + (out_names[0][0] == '%'));
                if (linf2_pred != UINT32_MAX) {
                    qir_set_op_pred(m, m->vals[out_ids[0]].def_op, linf2_pred);
                }
            } else if (strcmp(opname, "hmr_meas") == 0) {
                char q_t[128];
                uint32_t iq = 0u;
                if (!qir__lex_token(&L, q_t, sizeof(q_t)) || !qir__sym_get(&sym, q_t, &iq)) goto fail;
                out_ids[0] = qir_emit_hmr_meas(m, iq, out_names[0] + (out_names[0][0] == '%'));
            } else if (strcmp(opname, "bit_xor") == 0 || strcmp(opname, "bit_and") == 0 || strcmp(opname, "bit_not") == 0) {
                uint8_t mode = (strcmp(opname, "bit_xor") == 0) ? 0u : ((strcmp(opname, "bit_and") == 0) ? 1u : 2u);
                char a_t[128], b_t[128];
                uint32_t ia = 0u, ib = 0u;
                if (!qir__lex_token(&L, a_t, sizeof(a_t)) || !qir__sym_get(&sym, a_t, &ia)) goto fail;
                if (mode != 2u) {
                    if (!qir__lex_token(&L, b_t, sizeof(b_t)) || !qir__sym_get(&sym, b_t, &ib)) goto fail;
                }
                out_ids[0] = qir_emit_bit_op(m, mode, ia, ib, out_names[0] + (out_names[0][0] == '%'));
            } else if (strcmp(opname, "approx_fn") == 0) {
                char x_t[128];
                uint32_t ix = 0u;
                if (!qir__lex_token(&L, x_t, sizeof(x_t)) || !qir__sym_get(&sym, x_t, &ix)) goto fail;
                uint8_t fn_kind = 0u;
                uint16_t k_seg = 4u, deg = 2u;
                uint16_t frac_out = (out_widths[0] > 1u) ? (uint16_t)(out_widths[0] - 1u) : 1u;
                while (qir__lex_peek(&L) == ':') {
                    if (!qir__lex_token(&L, tok, sizeof(tok))) goto fail;
                    if (strcmp(tok, ":fn") == 0) {
                        char fns[32];
                        if (!qir__lex_token(&L, fns, sizeof(fns))) goto fail;
                        fn_kind = (strcmp(fns, "recip") == 0) ? 1u : 0u;
                    } else if (strcmp(tok, ":seg") == 0) {
                        char vs[32];
                        if (!qir__lex_token(&L, vs, sizeof(vs))) goto fail;
                        k_seg = (uint16_t)strtoul(vs, NULL, 10);
                    } else if (strcmp(tok, ":deg") == 0) {
                        char vs[32];
                        if (!qir__lex_token(&L, vs, sizeof(vs))) goto fail;
                        deg = (uint16_t)strtoul(vs, NULL, 10);
                    } else if (strcmp(tok, ":frac") == 0) {
                        char vs[32];
                        if (!qir__lex_token(&L, vs, sizeof(vs))) goto fail;
                        frac_out = (uint16_t)strtoul(vs, NULL, 10);
                    } else {
                        goto fail;
                    }
                }
                out_ids[0] = qir_emit_approx_fn(m, ix, fn_kind, k_seg, deg, out_widths[0], frac_out,
                                                out_names[0] + (out_names[0][0] == '%'));
            } else if (strcmp(opname, "phase_poly") == 0) {
                uint32_t ins[64];
                uint16_t ni = 0u;
                if (!qir__lex_expect(&L, '(')) goto fail;
                while (qir__lex_peek(&L) != ')') {
                    if (!qir__lex_token(&L, tok, sizeof(tok)) || !qir__sym_get(&sym, tok, &ins[ni])) goto fail;
                    ni++;
                }
                if (!qir__lex_expect(&L, ')')) goto fail;
                size_t cap = 32u, nterms = 0u;
                qir_phase_term_t *terms = (qir_phase_term_t *)calloc(cap, sizeof(qir_phase_term_t));
                while (qir__lex_peek(&L) == '(') {
                    if (!qir__lex_expect(&L, '(')) { free(terms); goto fail; }
                    if (nterms >= cap) {
                        cap *= 2u;
                        terms = (qir_phase_term_t *)realloc(terms, cap * sizeof(qir_phase_term_t));
                    }
                    qir_phase_term_t *pt = &terms[nterms];
                    memset(pt, 0, sizeof(*pt));
                    char kind_s[16], mask_s[128];
                    if (!qir__lex_token(&L, kind_s, sizeof(kind_s)) || !qir__lex_token(&L, mask_s, sizeof(mask_s))) {
                        free(terms); goto fail;
                    }
                    pt->is_xor = (strcmp(kind_s, "xor") == 0) ? 1u : 0u;
                    if (!qir__parse_bv256(mask_s, pt->mask)) { free(terms); goto fail; }
                    if (qir__lex_peek(&L) == ':') {
                        char kw[32], arg[64];
                        if (!qir__lex_token(&L, kw, sizeof(kw)) || !qir__lex_token(&L, arg, sizeof(arg))) {
                            free(terms); goto fail;
                        }
                        if (strcmp(kw, ":turns") == 0) {
                            double f = strtod(arg, NULL);
                            f = f - (double)(int64_t)f;
                            if (f < 0.0) f += 1.0;
                            double sc = f * 18446744073709551616.0;
                            uint64_t hi = (uint64_t)sc;
                            uint64_t lo = (uint64_t)((sc - (double)hi) * 18446744073709551616.0);
                            pt->angle = ((flow_u128)hi << 64u) | (flow_u128)lo;
                        }
                    } else {
                        char hi_s[64], lo_s[64];
                        if (!qir__lex_token(&L, hi_s, sizeof(hi_s)) || !qir__lex_token(&L, lo_s, sizeof(lo_s))) {
                            free(terms); goto fail;
                        }
                        uint64_t hi = (uint64_t)strtoull(hi_s, NULL, 0);
                        uint64_t lo = (uint64_t)strtoull(lo_s, NULL, 0);
                        pt->angle = ((flow_u128)hi << 64u) | (flow_u128)lo;
                    }
                    if (!qir__lex_expect(&L, ')')) { free(terms); goto fail; }
                    nterms++;
                }
                out_ids[0] = qir_emit_phase_poly(m, ins, ni, terms, (uint32_t)nterms,
                                                 out_names[0] + (out_names[0][0] == '%'));
                free(terms);
            } else {
                goto fail;
            }

            /* Parse optional trailing attributes :frac, :approx, :clean_anc, :if */
            qir_op_t *last_op = &m->ops[m->num_ops - 1u];
            while (qir__lex_peek(&L) == ':') {
                if (!qir__lex_token(&L, tok, sizeof(tok))) goto fail;
                if (strcmp(tok, ":approx") == 0) {
                    char tol_s[64];
                    if (!qir__lex_token(&L, tol_s, sizeof(tol_s))) goto fail;
                    flags |= QIR_FLAG_APPROX;
                    approx_tol = (uint32_t)strtoul(tol_s, NULL, 0);
                } else if (strcmp(tok, ":clean_anc") == 0) {
                    flags |= QIR_FLAG_CLEAN_ANC;
                } else if (strcmp(tok, ":if") == 0) {
                    char pt[128];
                    uint32_t pid = 0u;
                    if (!qir__lex_token(&L, pt, sizeof(pt)) || !qir__sym_get(&sym, pt, &pid)) goto fail;
                    qir_set_op_pred(m, last_op->id, pid);
                } else if (strcmp(tok, ":frac") == 0) {
                    char f0[32], f1[32], f2[32];
                    if (!qir__lex_token(&L, f0, sizeof(f0)) || !qir__lex_token(&L, f1, sizeof(f1))) goto fail;
                    if (last_op->kind == QIR_OP_MUL_FXP) {
                        if (!qir__lex_token(&L, f2, sizeof(f2))) goto fail;
                        last_op->u.mul.frac_a   = (uint16_t)strtoul(f0, NULL, 10);
                        last_op->u.mul.frac_b   = (uint16_t)strtoul(f1, NULL, 10);
                        last_op->u.mul.frac_out = (uint16_t)strtoul(f2, NULL, 10);
                    } else if (last_op->kind == QIR_OP_SQR_FXP) {
                        last_op->u.mul.frac_a   = (uint16_t)strtoul(f0, NULL, 10);
                        last_op->u.mul.frac_b   = last_op->u.mul.frac_a;
                        last_op->u.mul.frac_out = (uint16_t)strtoul(f1, NULL, 10);
                    }
                } else {
                    goto fail;
                }
            }
            last_op->flags |= flags;
            if (flags & QIR_FLAG_APPROX) {
                last_op->approx_tol_q16 = approx_tol;
                if (last_op->kind == QIR_OP_MUL_FXP || last_op->kind == QIR_OP_SQR_FXP) {
                    last_op->u.mul.keep_exact = 0u;
                }
            }

            if (!qir__lex_expect(&L, ')')) goto fail; /* end op */
            if (!qir__lex_expect(&L, ')')) goto fail; /* end let */
            for (uint8_t k = 0u; k < n_outs; ++k) {
                qir__sym_put(m, &sym, out_names[k], out_ids[k]);
            }
        } else {
            goto fail;
        }
    }

    if (!qir__lex_expect(&L, ')')) goto fail;
    free((void *)sym.names);
    free(sym.ids);
    qir_recompute_zero_masks(m);
    return true;

fail:
    free((void *)sym.names);
    free(sym.ids);
    qir_module_free(m);
    return false;
}

static inline bool qir_parse_file(qir_module_t *m, qcc_arena_t *arena, const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) return false;
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (sz < 0) {
        fclose(f);
        return false;
    }
    char *buf = (char *)malloc((size_t)sz + 1u);
    if (!buf) {
        fclose(f);
        return false;
    }
    size_t nr = fread(buf, 1, (size_t)sz, f);
    fclose(f);
    buf[nr] = '\0';
    bool ok = qir_parse_str(m, arena, buf, nr);
    free(buf);
    return ok;
}

#ifdef __cplusplus
}
#endif

#endif /* QIR_DSL_H */
