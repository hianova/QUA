/*
 * qir_fir.h — Unified FLOW_IR (.fir) -> QIR Bridge & QTT Linear Gate
 * Enforces Flow.Morphism Multiplicity == FLOW_IR_MULT_LINEAR before lowering
 * into QCC's reversible QIR -> .kmb fault-tolerant compiler.
 */
#ifndef QIR_FIR_H
#define QIR_FIR_H

#include "flow_ir.h"
#include "qir.h"

#ifdef __cplusplus
extern "C" {
#endif

static inline bool qir_load_fir_bytes(qir_module_t *m, qcc_arena_t *arena,
                                      const uint8_t *buf, size_t len) {
    if (!m || !arena || !buf) {
        return false;
    }
    flow_ir_prog_t *prog = (flow_ir_prog_t *)calloc(1u, sizeof(flow_ir_prog_t));
    if (!prog) {
        return false;
    }
    if (!flow_ir_load_bytes(prog, buf, len)) {
        free(prog);
        return false;
    }
    if (!flow_ir_admits_unitary(prog)) {
        fprintf(stderr,
                "qir_fir: QTT error — morphism '%s' has multiplicity=%u (quantum requires linear)\n",
                prog->name, (unsigned)prog->mult);
        free(prog);
        return false;
    }

    qir_module_init(m, arena, prog->name);
    m->require_clean_scratch = (prog->flags & FLOW_IR_FLAG_CLEAN) != 0u;

    uint32_t vmap[FLOW_IR_MAX_VALS];
    for (uint32_t i = 0u; i < FLOW_IR_MAX_VALS; ++i) {
        vmap[i] = UINT32_MAX;
    }

    for (uint32_t i = 0u; i < prog->num_inputs; ++i) {
        uint32_t vid = prog->inputs[i];
        const flow_ir_val_t *fv = &prog->vals[vid];
        vmap[vid] = fv->is_classical
            ? qir_add_classical_input(m, fv->width, fv->name)
            : qir_add_input(m, fv->width, fv->name);
    }

    bool ok = true;
    for (uint32_t i = 0u; i < prog->num_ops && ok; ++i) {
        const flow_ir_op_t *op = &prog->ops[i];
        uint32_t oid           = op->out_val;
        if (oid >= prog->num_vals) {
            ok = false;
            break;
        }
        uint16_t ow       = prog->vals[oid].width;
        const char *vname = prog->vals[oid].name;
        uint16_t k        = op->kind;
        uint32_t qv       = UINT32_MAX;

        for (uint16_t j = 0u; j < op->num_ins; ++j) {
            if (op->ins[j] >= oid || vmap[op->ins[j]] == UINT32_MAX) {
                ok = false;
                break;
            }
        }
        if (!ok) {
            break;
        }

        if (k == FLOW_IR_OP_CONST_U64) {
            qv = qir_emit_const_u64(m, ow, op->p0, vname);
        } else if (k == FLOW_IR_OP_SLICE && op->num_ins == 1u) {
            qv = qir_emit_slice(m, vmap[op->ins[0]], (uint16_t)op->p0, (uint16_t)(op->p0 + ow), vname);
        } else if (k == FLOW_IR_OP_CONCAT && op->num_ins >= 1u) {
            uint32_t cins[FLOW_IR_MAX_INS];
            for (uint16_t j = 0u; j < op->num_ins; ++j) {
                cins[j] = vmap[op->ins[j]];
            }
            qv = qir_emit_concat(m, cins, op->num_ins, vname);
        } else if (k == FLOW_IR_OP_XOR_VEC && op->num_ins == 2u) {
            qv = qir_emit_gf2m_add(m, vmap[op->ins[0]], vmap[op->ins[1]], ow, vname);
        } else if (k == FLOW_IR_OP_AND_VEC && op->num_ins == 2u) {
            qv = qir_emit_and_vec(m, vmap[op->ins[0]], vmap[op->ins[1]], 0u, vname);
        } else if (k == FLOW_IR_OP_NOT_VEC && op->num_ins == 1u) {
            uint32_t ones = qir_emit_const_u64(m, ow, flow_ir_mask_u64(ow), NULL);
            qv = qir_emit_gf2m_add(m, vmap[op->ins[0]], ones, ow, vname);
        } else if (k == FLOW_IR_OP_OR_REDUCE && op->num_ins == 1u) {
            uint32_t ez  = qir_emit_eq_zero(m, vmap[op->ins[0]], NULL);
            uint32_t one = qir_emit_const_u64(m, 1u, UINT64_C(1), NULL);
            qv = qir_emit_gf2m_add(m, ez, one, 1u, vname);
        } else if (k == FLOW_IR_OP_ADD && op->num_ins == 2u) {
            qv = qir_emit_add(m, vmap[op->ins[0]], vmap[op->ins[1]], 0u, vname);
        } else if (k == FLOW_IR_OP_SUB && op->num_ins == 2u) {
            qv = qir_emit_sub(m, vmap[op->ins[0]], vmap[op->ins[1]], 0u, vname);
        } else if (k == FLOW_IR_OP_MUL_FXP && op->num_ins == 2u) {
            qv = qir_emit_mul_fxp(m, vmap[op->ins[0]], vmap[op->ins[1]], ow,
                                  (uint16_t)op->p0, (uint16_t)op->p1, (uint16_t)op->p2, 0u, vname);
        } else if (k == FLOW_IR_OP_SQR_FXP && op->num_ins == 1u) {
            qv = qir_emit_sqr_fxp(m, vmap[op->ins[0]], ow,
                                  (uint16_t)op->p0, (uint16_t)op->p2, 0u, vname);
        } else if (k == FLOW_IR_OP_GF2M_MUL && op->num_ins == 2u) {
            qv = qir_emit_gf2m_mul(m, vmap[op->ins[0]], vmap[op->ins[1]],
                                   (uint16_t)op->p0, (uint32_t)op->p1, vname);
        } else if (k == FLOW_IR_OP_GF2M_INV && op->num_ins == 1u) {
            qv = qir_emit_gf2m_inv(m, vmap[op->ins[0]], (uint16_t)op->p0, (uint32_t)op->p1, vname);
        } else if (k == FLOW_IR_OP_GF2M_INV && op->num_ins >= 2u) {
            qv = qir_emit_gf2m_inv_ctrl(m, vmap[op->ins[0]], vmap[op->ins[1]],
                                        (uint16_t)op->p0, (uint32_t)op->p1, vname);
        } else if (k == FLOW_IR_OP_LUT_ROM && op->num_ins == 1u) {
            qv = qir_emit_lut_rom(m, vmap[op->ins[0]], ow, op->lut_len,
                                  &prog->lut_pool[op->lut_off], 0u, vname);
        } else if (k == FLOW_IR_OP_CTRL_SEL && op->num_ins == 3u) {
            qv = qir_emit_ctrl_sel(m, vmap[op->ins[0]], vmap[op->ins[1]], vmap[op->ins[2]],
                                   0u, vname);
        } else if (k == FLOW_IR_OP_SHL_VAR && op->num_ins == 2u) {
            uint16_t sw  = prog->vals[op->ins[1]].width;
            uint16_t msh = op->p0 ? (uint16_t)op->p0
                                  : (uint16_t)((sw < 16u) ? ((1u << sw) - 1u) : 65535u);
            qv = qir_emit_shl_var(m, vmap[op->ins[0]], vmap[op->ins[1]], ow, msh, vname);
        } else if (k == FLOW_IR_OP_SHR_VAR && op->num_ins == 2u) {
            uint16_t sw  = prog->vals[op->ins[1]].width;
            uint16_t msh = op->p0 ? (uint16_t)op->p0
                                  : (uint16_t)((sw < 16u) ? ((1u << sw) - 1u) : 65535u);
            qv = qir_emit_shr_var(m, vmap[op->ins[0]], vmap[op->ins[1]], ow, msh, vname);
        } else if (k == FLOW_IR_OP_EQ_ZERO && op->num_ins == 1u) {
            qv = qir_emit_eq_zero(m, vmap[op->ins[0]], vname);
        } else {
            ok = false;
            break;
        }
        vmap[oid] = qv;
    }

    if (!ok) {
        qir_module_free(m);
        free(prog);
        return false;
    }

    for (uint32_t i = 0u; i < prog->num_outputs; ++i) {
        uint32_t oid = prog->outputs[i];
        if (oid >= prog->num_vals || vmap[oid] == UINT32_MAX) {
            qir_module_free(m);
            free(prog);
            return false;
        }
        if (prog->vals[oid].is_classical) {
            qir_add_classical_output(m, vmap[oid]);
        } else {
            qir_add_output(m, vmap[oid]);
        }
    }

    qir_recompute_zero_masks(m);
    free(prog);
    return true;
}

static inline bool qir_load_fir_file(qir_module_t *m, qcc_arena_t *arena, const char *path) {
    if (!m || !arena || !path) {
        return false;
    }
    FILE *fp = fopen(path, "rb");
    if (!fp) {
        return false;
    }
    uint8_t *buf = (uint8_t *)calloc(1u, FLOW_IR_MAX_FILE_BYTES);
    if (!buf) {
        fclose(fp);
        return false;
    }
    size_t n  = fread(buf, 1u, FLOW_IR_MAX_FILE_BYTES, fp);
    int extra = fgetc(fp);
    fclose(fp);
    if (extra != EOF || n == 0u) {
        free(buf);
        return false;
    }
    bool ok = qir_load_fir_bytes(m, arena, buf, n);
    free(buf);
    return ok;
}

#ifdef __cplusplus
}
#endif

#endif /* QIR_FIR_H */
