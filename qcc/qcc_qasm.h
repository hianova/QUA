/*
 * qcc_qasm.h — Pure C17 / C++20 OpenQASM 3.0 Exporter for .kmb v2 Circuits
 * Part of QCC v3 (QUA). Zero external dependencies.
 */
#ifndef QCC_QASM_H
#define QCC_QASM_H

#include "qcc_kmb.h"

#ifdef __cplusplus
extern "C" {
#endif

static inline void qcc__qasm_emit_cond_prefix(FILE *out, const uint32_t *stack, size_t stack_len,
                                              bool has_cond, uint32_t cond_b, bool *out_wrapped) {
    if (stack_len == 0u && !has_cond) {
        *out_wrapped = false;
        return;
    }
    *out_wrapped = true;
    fprintf(out, "if (");
    bool first = true;
    for (size_t i = 0u; i < stack_len; ++i) {
        fprintf(out, "%sc[%u]", first ? "" : " && ", (unsigned)stack[i]);
        first = false;
    }
    if (has_cond) {
        fprintf(out, "%sc[%u]", first ? "" : " && ", (unsigned)cond_b);
    }
    fprintf(out, ") { ");
}

static inline void qcc__qasm_emit_z_pow(FILE *out, uint32_t q0, qcc_angle128_t angle) {
    flow_u128 ang = qcc_u128_from_angle(angle);
    if (ang == 0u) {
        return;
    }
    const flow_u128 t1 = (flow_u128)1u << 125u;
    if (ang == ((flow_u128)1u << 127u)) {
        fprintf(out, "z q[%u];", (unsigned)q0);
    } else if (ang == ((flow_u128)1u << 126u)) {
        fprintf(out, "s q[%u];", (unsigned)q0);
    } else if (ang == ((flow_u128)3u << 126u)) {
        fprintf(out, "sdg q[%u];", (unsigned)q0);
    } else if (ang == t1) {
        fprintf(out, "t q[%u];", (unsigned)q0);
    } else if (ang == (t1 * 7u)) {
        fprintf(out, "tdg q[%u];", (unsigned)q0);
    } else {
        /* 2 * signed_fraction(turns) * pi */
        uint64_t hi_u   = (uint64_t)(ang >> 64u);
        int64_t hi_s    = (int64_t)hi_u;
        double turns    = (double)hi_s * (1.0 / 18446744073709551616.0);
        uint64_t lo_u   = (uint64_t)ang;
        if (lo_u != 0u) {
            turns += ((double)lo_u * (1.0 / 18446744073709551616.0)) * (1.0 / 18446744073709551616.0);
        }
        fprintf(out, "rz(%.15g * pi) q[%u];", 2.0 * turns, (unsigned)q0);
    }
}

static inline bool qcc_kmb_export_openqasm3(const qcc_kmb_circuit_t *c, FILE *out,
                                            uint64_t *out_lines) {
    if (!c || !out) {
        return false;
    }
    qcc_kmb_stats_t ks;
    if (!qcc_kmb_recompute_stats(c, &ks)) {
        return false;
    }
    uint32_t nq = ks.num_qubits > c->num_qubits ? ks.num_qubits : c->num_qubits;
    uint32_t nb = ks.num_bits > c->num_bits ? ks.num_bits : c->num_bits;
    if (nq == 0u) {
        nq = 1u;
    }

    uint64_t lines = 0u;
    fprintf(out, "OPENQASM 3.0;\ninclude \"stdgates.inc\";\nqubit[%u] q;\n", (unsigned)nq);
    lines += 3u;
    if (nb > 0u) {
        fprintf(out, "bit[%u] c;\n", (unsigned)nb);
        lines += 1u;
    }

    size_t cond_cap = 64u, cond_len = 0u;
    uint32_t *cond_stack = (uint32_t *)malloc(cond_cap * sizeof(uint32_t));
    if (!cond_stack) {
        return false;
    }

    qcc_kmb_cursor_t cur;
    qcc_kmb_cursor_init(&cur);
    qcc_kmb_op_view_t v;
    while (qcc_kmb_cursor_next(c, &cur, &v)) {
        if (v.op == QCC_OP_PUSH_CONDITION) {
            if (cond_len >= cond_cap) {
                cond_cap *= 2u;
                uint32_t *ns = (uint32_t *)realloc(cond_stack, cond_cap * sizeof(uint32_t));
                if (!ns) {
                    free(cond_stack);
                    return false;
                }
                cond_stack = ns;
            }
            cond_stack[cond_len++] = v.cond;
            continue;
        }
        if (v.op == QCC_OP_POP_CONDITION) {
            if (cond_len > 0u) {
                cond_len--;
            }
            continue;
        }

        uint8_t base = v.op & (uint8_t)(~QCC_OP_FLAG_COND);
        if (base == QCC_OP_INFO_TAG || base == QCC_OP_DEBUG_PRINT_EMPTY ||
            base == QCC_OP_DEBUG_PRINT_C || base == QCC_OP_DEBUG_PRINT_Q) {
            continue;
        }

        bool wrapped = false;
        qcc__qasm_emit_cond_prefix(out, cond_stack, cond_len, v.has_cond != 0u, v.cond, &wrapped);
        switch (base) {
            case QCC_OP_NEG:
                fprintf(out, "gphase(pi);");
                break;
            case QCC_OP_BIT_STORE0:
                fprintf(out, "c[%u] = 0;", (unsigned)v.b);
                break;
            case QCC_OP_BIT_STORE1:
                fprintf(out, "c[%u] = 1;", (unsigned)v.b);
                break;
            case QCC_OP_BIT_INVERT:
                fprintf(out, "c[%u] = ~c[%u];", (unsigned)v.b, (unsigned)v.b);
                break;
            case QCC_OP_X:
                fprintf(out, "x q[%u];", (unsigned)v.q[0]);
                break;
            case QCC_OP_Z:
                fprintf(out, "z q[%u];", (unsigned)v.q[0]);
                break;
            case QCC_OP_R:
                fprintf(out, "reset q[%u];", (unsigned)v.q[0]);
                break;
            case QCC_OP_Z_POW:
                qcc__qasm_emit_z_pow(out, v.q[0], v.angle);
                break;
            case QCC_OP_HMR:
                fprintf(out, "h q[%u]; c[%u] = measure q[%u]; reset q[%u];",
                        (unsigned)v.q[0], (unsigned)v.b, (unsigned)v.q[0], (unsigned)v.q[0]);
                break;
            case QCC_OP_CX:
                fprintf(out, "cx q[%u], q[%u];", (unsigned)v.q[0], (unsigned)v.q[1]);
                break;
            case QCC_OP_CZ:
                fprintf(out, "cz q[%u], q[%u];", (unsigned)v.q[0], (unsigned)v.q[1]);
                break;
            case QCC_OP_SWAP:
                fprintf(out, "swap q[%u], q[%u];", (unsigned)v.q[0], (unsigned)v.q[1]);
                break;
            case QCC_OP_CCX:
                fprintf(out, "ctrl(2) @ x q[%u], q[%u], q[%u];",
                        (unsigned)v.q[0], (unsigned)v.q[1], (unsigned)v.q[2]);
                break;
            case QCC_OP_CCZ:
                fprintf(out, "ctrl(2) @ z q[%u], q[%u], q[%u];",
                        (unsigned)v.q[0], (unsigned)v.q[1], (unsigned)v.q[2]);
                break;
            case QCC_OP_CSWAP:
                fprintf(out, "ctrl @ swap q[%u], q[%u], q[%u];",
                        (unsigned)v.q[0], (unsigned)v.q[1], (unsigned)v.q[2]);
                break;
            default:
                break;
        }
        fprintf(out, "%s\n", wrapped ? " }" : "");
        lines += 1u;
    }

    free(cond_stack);
    if (out_lines) {
        *out_lines = lines;
    }
    return true;
}

static inline bool qcc_kmb_write_openqasm3(const qcc_kmb_circuit_t *c, const char *path,
                                           uint64_t *out_lines) {
    if (!c || !path) {
        return false;
    }
    FILE *fp = fopen(path, "w");
    if (!fp) {
        return false;
    }
    bool ok = qcc_kmb_export_openqasm3(c, fp, out_lines);
    if (fclose(fp) != 0) {
        ok = false;
    }
    return ok;
}

#ifdef __cplusplus
}
#endif

#endif /* QCC_QASM_H */
