#ifndef QCTRL_DECODER_RING0_H
#define QCTRL_DECODER_RING0_H

/* ============================================================================
 * FLOW-QCtrl Ring-0 Pure-GPR Rotated Surface Code Decoder (d = 3, 5, 7)
 * ============================================================================
 * Pure integer freestanding kernel compatible with C17 & C++20:
 *   - Zero FPU / SIMD instructions (-mgeneral-regs-only clean)
 *   - Zero udiv / sdiv / __divti3 (power-of-two shifts & bitmasks only)
 *   - Zero heap allocation & zero external libc symbols (no memset/memcpy)
 *   - Zero indirect branches / jump tables (-fno-jump-tables clean)
 *   - Exact Bitmask-DP MWPM for syndrome weight k <= 8 (provably exact for
 *     all weight <= (d-1)/2 errors on d = 3, 5, 7) + greedy Ring-0 fallback.
 *   - Syndrome-Driven EWMA+CUSUM Drift Monitor (zero extra probe overhead).
 * ============================================================================ */

#include "flow.h"

#ifdef __cplusplus
extern "C" {
#endif

#define QCTRL_DEC_MAX_D     7U
#define QCTRL_DEC_MAX_DATA  49U
#define QCTRL_DEC_MAX_STAB  24U
#define QCTRL_DEC_DP_MAX_K  8U
#define QCTRL_DEC_INF_DIST  0x3FFFU

typedef struct {
    uint32_t d;
    uint32_t n_data;
    uint32_t n_stab; /* per Pauli type: (d^2 - 1)/2 <= 24 */
    uint64_t z_stab_mask[QCTRL_DEC_MAX_STAB]; /* detects X errors */
    uint64_t x_stab_mask[QCTRL_DEC_MAX_STAB]; /* detects Z errors */
    uint64_t log_x_mask;                      /* logical X = column c=0 */
    uint64_t log_z_mask;                      /* logical Z = row r=0 */

    /* Matching graph shortest paths (g=0: X-err / Z-stab, g=1: Z-err / X-stab) */
    uint8_t  bdist[2][QCTRL_DEC_MAX_STAB];
    uint64_t bpath_mask[2][QCTRL_DEC_MAX_STAB];
    uint8_t  pair_dist[2][QCTRL_DEC_MAX_STAB][QCTRL_DEC_MAX_STAB];
    uint64_t pair_mask[2][QCTRL_DEC_MAX_STAB][QCTRL_DEC_MAX_STAB];

    /* O(1) lookup table for d=3 (16 syndromes -> 9-bit data correction) */
    uint16_t lut_d3[2][16];
} qctrl_sc_code_t;

typedef struct {
    int64_t base_rate_q16;  /* baseline syndrome defect fraction in Q16.16 */
    int64_t ewma_rate_q16;  /* EWMA defect rate in Q16.16 (alpha = 1/8) */
    int64_t cusum_q16;      /* one-sided positive CUSUM in Q16.16 */
    int64_t thresh_q16;     /* recalibration trigger threshold */
    uint32_t rounds;
    uint32_t triggers;
} qctrl_synd_monitor_t;

/* Pure-GPR SWAR popcount (zero SIMD cnt / libc call) */
FLOW_INL uint32_t qctrl__popcnt64_gpr(uint64_t x) {
    x = x - ((x >> 1U) & UINT64_C(0x5555555555555555));
    x = (x & UINT64_C(0x3333333333333333)) + ((x >> 2U) & UINT64_C(0x3333333333333333));
    x = (x + (x >> 4U)) & UINT64_C(0x0F0F0F0F0F0F0F0F);
    return (uint32_t)((x * UINT64_C(0x0101010101010101)) >> 56U);
}

FLOW_INL uint32_t qctrl__parity64_gpr(uint64_t x) {
    x ^= x >> 32U;
    x ^= x >> 16U;
    x ^= x >> 8U;
    x ^= x >> 4U;
    return (uint32_t)((UINT64_C(0x6996) >> (x & 0xFU)) & 1U);
}

FLOW_INL uint32_t qctrl__ctz64_gpr(uint64_t x) {
    uint32_t c = 0U;
    if (x == 0ULL) {
        return 64U;
    }
    if ((x & UINT64_C(0xFFFFFFFF)) == 0ULL) { x >>= 32U; c += 32U; }
    if ((x & UINT64_C(0xFFFF)) == 0ULL)     { x >>= 16U; c += 16U; }
    if ((x & UINT64_C(0xFF)) == 0ULL)       { x >>= 8U;  c += 8U;  }
    if ((x & UINT64_C(0xF)) == 0ULL)        { x >>= 4U;  c += 4U;  }
    if ((x & UINT64_C(0x3)) == 0ULL)        { x >>= 2U;  c += 2U;  }
    if ((x & UINT64_C(0x1)) == 0ULL)        { c += 1U; }
    return c;
}

FLOW_INL void qctrl__build_graph_apsp(qctrl_sc_code_t *code, uint32_t g,
                                      const uint64_t *stab_masks) {
    uint32_t n_s = code->n_stab;
    uint32_t n_q = code->n_data;
    uint8_t adj_cnt[QCTRL_DEC_MAX_STAB][QCTRL_DEC_MAX_STAB];
    uint64_t edge_qmask[QCTRL_DEC_MAX_STAB][QCTRL_DEC_MAX_STAB];
    uint32_t u, v, q, s;

    for (u = 0U; u < QCTRL_DEC_MAX_STAB; ++u) {
        code->bdist[g][u]      = 0x7FU;
        code->bpath_mask[g][u] = 0ULL;
        for (v = 0U; v < QCTRL_DEC_MAX_STAB; ++v) {
            adj_cnt[u][v]          = 0U;
            edge_qmask[u][v]       = 0ULL;
            code->pair_dist[g][u][v] = (u == v) ? 0U : 0x7FU;
            code->pair_mask[g][u][v] = 0ULL;
        }
    }

    /* Map each data qubit to an interior edge (u <-> v) or a boundary edge (u <-> bnd) */
    for (q = 0U; q < n_q; ++q) {
        uint64_t qbit = UINT64_C(1) << q;
        uint32_t hits[2] = {0U, 0U};
        uint32_t nh = 0U;
        for (s = 0U; s < n_s; ++s) {
            if ((stab_masks[s] & qbit) != 0ULL) {
                if (nh < 2U) {
                    hits[nh] = s;
                }
                nh += 1U;
            }
        }
        if (nh == 1U) {
            uint32_t su = hits[0];
            if (code->bdist[g][su] > 1U) {
                code->bdist[g][su]      = 1U;
                code->bpath_mask[g][su] = qbit;
            }
        } else if (nh == 2U) {
            uint32_t su = hits[0];
            uint32_t sv = hits[1];
            adj_cnt[su][sv]    = 1U;
            adj_cnt[sv][su]    = 1U;
            edge_qmask[su][sv] = qbit;
            edge_qmask[sv][su] = qbit;
        }
    }

    /* Unweighted BFS from each source stabilizer u <= 23 -> exact shortest path & XOR mask */
    for (u = 0U; u < n_s; ++u) {
        uint32_t queue[QCTRL_DEC_MAX_STAB];
        uint32_t qh = 0U, qt = 0U;
        uint8_t dist[QCTRL_DEC_MAX_STAB];
        uint64_t pmask[QCTRL_DEC_MAX_STAB];
        for (v = 0U; v < n_s; ++v) {
            dist[v]  = 0x7FU;
            pmask[v] = 0ULL;
        }
        dist[u]      = 0U;
        pmask[u]     = 0ULL;
        queue[qt++]  = u;
        while (qh < qt) {
            uint32_t cur = queue[qh++];
            uint8_t nd   = (uint8_t)(dist[cur] + 1U);
            for (v = 0U; v < n_s; ++v) {
                if (adj_cnt[cur][v] != 0U && dist[v] == 0x7FU) {
                    dist[v]     = nd;
                    pmask[v]    = pmask[cur] ^ edge_qmask[cur][v];
                    queue[qt++] = v;
                }
            }
        }
        for (v = 0U; v < n_s; ++v) {
            code->pair_dist[g][u][v] = dist[v];
            code->pair_mask[g][u][v] = pmask[v];
        }
    }

    /* Relax boundary distance over all shortest paths u -> v -> boundary */
    for (u = 0U; u < n_s; ++u) {
        uint8_t best_d  = code->bdist[g][u];
        uint64_t best_m = code->bpath_mask[g][u];
        for (v = 0U; v < n_s; ++v) {
            if (code->bdist[g][v] == 1U && code->pair_dist[g][u][v] < 0x7FU) {
                uint8_t cand = (uint8_t)(code->pair_dist[g][u][v] + 1U);
                if (cand < best_d) {
                    best_d = cand;
                    best_m = code->pair_mask[g][u][v] ^ code->bpath_mask[g][v];
                }
            }
        }
        code->bdist[g][u]      = best_d;
        code->bpath_mask[g][u] = best_m;
    }
}

FLOW_INL uint64_t qctrl__solve_mwpm_ring0(const qctrl_sc_code_t *code, uint32_t synd, uint32_t g) {
    if (synd == 0U) {
        return 0ULL;
    }
    uint32_t verts[QCTRL_DEC_MAX_STAB];
    uint32_t k = 0U;
    uint32_t tmp = synd;
    while (tmp != 0U) {
        uint32_t bit = qctrl__ctz64_gpr((uint64_t)tmp);
        tmp &= (tmp - 1U);
        if (bit < code->n_stab && k < QCTRL_DEC_MAX_STAB) {
            verts[k++] = bit;
        }
    }
    if (k == 0U) {
        return 0ULL;
    }

    /* Exact Bitmask-DP MWPM for up to k=8 triggered stabilizers (2^8 = 256 states) */
    if (k <= QCTRL_DEC_DP_MAX_K) {
        uint32_t n_states = UINT32_C(1) << k;
        uint16_t dp[256];
        uint8_t choice[256];
        uint32_t m;
        dp[0]     = 0U;
        choice[0] = 0xFFU;
        for (m = 1U; m < n_states; ++m) {
            uint32_t i      = qctrl__ctz64_gpr((uint64_t)m);
            uint32_t rest_i = m ^ (UINT32_C(1) << i);
            uint16_t best   = (uint16_t)((uint16_t)code->bdist[g][verts[i]] + dp[rest_i]);
            uint8_t best_c  = (uint8_t)i; /* self == match to virtual boundary */
            uint32_t j;
            for (j = i + 1U; j < k; ++j) {
                if ((rest_i & (UINT32_C(1) << j)) != 0U) {
                    uint32_t rest_ij = rest_i ^ (UINT32_C(1) << j);
                    uint16_t cand =
                        (uint16_t)((uint16_t)code->pair_dist[g][verts[i]][verts[j]] + dp[rest_ij]);
                    if (cand < best) {
                        best   = cand;
                        best_c = (uint8_t)((i << 4U) | j | 0x80U);
                    }
                }
            }
            dp[m]     = best;
            choice[m] = best_c;
        }
        uint64_t corr = 0ULL;
        uint32_t cur  = n_states - 1U;
        while (cur != 0U) {
            uint8_t ch = choice[cur];
            if ((ch & 0x80U) == 0U) {
                uint32_t i = (uint32_t)ch;
                corr ^= code->bpath_mask[g][verts[i]];
                cur  ^= (UINT32_C(1) << i);
            } else {
                uint32_t i = (uint32_t)((ch >> 4U) & 0x7U);
                uint32_t j = (uint32_t)(ch & 0x7U);
                corr ^= code->pair_mask[g][verts[i]][verts[j]];
                cur  ^= (UINT32_C(1) << i) | (UINT32_C(1) << j);
            }
        }
        return corr;
    }

    /* Deterministic Pure-GPR Greedy Pairing + Boundary Fallback for heavy k > 8 bursts */
    uint8_t used[QCTRL_DEC_MAX_STAB];
    uint32_t rem = k;
    uint32_t i, j;
    for (i = 0U; i < k; ++i) {
        used[i] = 0U;
    }
    uint64_t corr = 0ULL;
    while (rem >= 2U) {
        uint16_t best_w  = QCTRL_DEC_INF_DIST;
        uint32_t best_i  = 0U;
        uint32_t best_j  = 0U;
        bool best_bnd    = false;
        for (i = 0U; i < k; ++i) {
            if (used[i] != 0U) { continue; }
            for (j = i + 1U; j < k; ++j) {
                if (used[j] != 0U) { continue; }
                uint16_t pw = (uint16_t)code->pair_dist[g][verts[i]][verts[j]];
                uint16_t bw = (uint16_t)((uint16_t)code->bdist[g][verts[i]] +
                                         (uint16_t)code->bdist[g][verts[j]]);
                if (bw < pw && bw < best_w) {
                    best_w   = bw;
                    best_i   = i;
                    best_j   = j;
                    best_bnd = true;
                } else if (pw < best_w) {
                    best_w   = pw;
                    best_i   = i;
                    best_j   = j;
                    best_bnd = false;
                }
            }
        }
        used[best_i] = 1U;
        used[best_j] = 1U;
        rem -= 2U;
        if (best_bnd) {
            corr ^= code->bpath_mask[g][verts[best_i]] ^ code->bpath_mask[g][verts[best_j]];
        } else {
            corr ^= code->pair_mask[g][verts[best_i]][verts[best_j]];
        }
    }
    for (i = 0U; i < k; ++i) {
        if (used[i] == 0U) {
            corr ^= code->bpath_mask[g][verts[i]];
        }
    }
    return corr;
}

FLOW_INL bool qctrl_sc_init_ring0(qctrl_sc_code_t *code, uint32_t d) {
    uint32_t i, j, r, c;
    if (!code || (d != 3U && d != 5U && d != 7U)) {
        return false;
    }
    code->d          = d;
    code->n_data     = d * d;
    code->n_stab     = (code->n_data - 1U) >> 1U;
    code->log_x_mask = 0ULL;
    code->log_z_mask = 0ULL;
    for (i = 0U; i < QCTRL_DEC_MAX_STAB; ++i) {
        code->z_stab_mask[i] = 0ULL;
        code->x_stab_mask[i] = 0ULL;
    }

    for (r = 0U; r < d; ++r) {
        code->log_x_mask |= (UINT64_C(1) << (r * d));
    }
    for (c = 0U; c < d; ++c) {
        code->log_z_mask |= (UINT64_C(1) << c);
    }

    /* Rotated surface code Z-stabilizers (detect X) */
    uint32_t nz = 0U;
    for (r = 0U; r + 1U < d; ++r) {
        for (c = 0U; c + 1U < d; ++c) {
            if (((r + c) & 1U) == 0U && nz < code->n_stab) {
                uint64_t m = (UINT64_C(1) << (r * d + c)) |
                             (UINT64_C(1) << (r * d + (c + 1U))) |
                             (UINT64_C(1) << ((r + 1U) * d + c)) |
                             (UINT64_C(1) << ((r + 1U) * d + (c + 1U)));
                code->z_stab_mask[nz++] = m;
            }
        }
    }
    for (r = 1U; r + 1U < d; r += 2U) {
        if (nz < code->n_stab) {
            code->z_stab_mask[nz++] = (UINT64_C(1) << (r * d)) | (UINT64_C(1) << ((r + 1U) * d));
        }
    }
    for (r = 0U; r + 1U < d; r += 2U) {
        if (nz < code->n_stab) {
            code->z_stab_mask[nz++] =
                (UINT64_C(1) << (r * d + (d - 1U))) | (UINT64_C(1) << ((r + 1U) * d + (d - 1U)));
        }
    }

    /* Rotated surface code X-stabilizers (detect Z) */
    uint32_t nx = 0U;
    for (r = 0U; r + 1U < d; ++r) {
        for (c = 0U; c + 1U < d; ++c) {
            if (((r + c) & 1U) == 1U && nx < code->n_stab) {
                uint64_t m = (UINT64_C(1) << (r * d + c)) |
                             (UINT64_C(1) << (r * d + (c + 1U))) |
                             (UINT64_C(1) << ((r + 1U) * d + c)) |
                             (UINT64_C(1) << ((r + 1U) * d + (c + 1U)));
                code->x_stab_mask[nx++] = m;
            }
        }
    }
    for (c = 0U; c + 1U < d; c += 2U) {
        if (nx < code->n_stab) {
            code->x_stab_mask[nx++] = (UINT64_C(1) << c) | (UINT64_C(1) << (c + 1U));
        }
    }
    for (c = 1U; c + 1U < d; c += 2U) {
        if (nx < code->n_stab) {
            code->x_stab_mask[nx++] =
                (UINT64_C(1) << ((d - 1U) * d + c)) | (UINT64_C(1) << ((d - 1U) * d + c + 1U));
        }
    }
    if (nz != code->n_stab || nx != code->n_stab) {
        return false;
    }

    qctrl__build_graph_apsp(code, 0U, code->z_stab_mask);
    qctrl__build_graph_apsp(code, 1U, code->x_stab_mask);

    for (i = 0U; i < 2U; ++i) {
        for (j = 0U; j < 16U; ++j) {
            code->lut_d3[i][j] = 0U;
        }
    }
    if (d == 3U) {
        uint32_t s;
        for (i = 0U; i < 2U; ++i) {
            for (s = 1U; s < 16U; ++s) {
                code->lut_d3[i][s] = (uint16_t)qctrl__solve_mwpm_ring0(code, s, i);
            }
        }
    }
    return true;
}

FLOW_INL uint32_t qctrl_sc_syndrome_ring0(const qctrl_sc_code_t *code, uint64_t err_mask,
                                          bool is_z_err) {
    const uint64_t *masks = is_z_err ? code->x_stab_mask : code->z_stab_mask;
    uint32_t n_s          = code->n_stab;
    uint32_t synd         = 0U;
    uint32_t s;
    for (s = 0U; s < n_s; ++s) {
        if (qctrl__parity64_gpr(err_mask & masks[s]) != 0U) {
            synd |= (UINT32_C(1) << s);
        }
    }
    return synd;
}

FLOW_INL uint64_t qctrl_sc_decode_ring0(const qctrl_sc_code_t *code, uint32_t synd, bool is_z_err) {
    uint32_t g = is_z_err ? 1U : 0U;
    if (synd == 0U) {
        return 0ULL;
    }
    if (code->d == 3U) {
        return (uint64_t)code->lut_d3[g][synd & 0xFU];
    }
    return qctrl__solve_mwpm_ring0(code, synd, g);
}

/* Ring-0 Syndrome-Driven EWMA + CUSUM Drift Monitor (zero division, zero extra probes) */
FLOW_INL void qctrl_synd_monitor_reset(qctrl_synd_monitor_t *mon, int64_t base_rate_q16,
                                       int64_t thresh_q16) {
    mon->base_rate_q16 = base_rate_q16;
    mon->ewma_rate_q16 = base_rate_q16;
    mon->cusum_q16     = 0;
    mon->thresh_q16    = thresh_q16;
    mon->rounds        = 0U;
    mon->triggers      = 0U;
}

FLOW_INL bool qctrl_synd_monitor_step_ring0(qctrl_synd_monitor_t *mon, uint32_t synd_x,
                                            uint32_t synd_z, uint32_t n_stab_per_type) {
    uint32_t defects = qctrl__popcnt64_gpr((uint64_t)synd_x) + qctrl__popcnt64_gpr((uint64_t)synd_z);
    /* Scale defects by 1/(2*n_stab_per_type) via fixed reciprocal shift 2^16 / denom */
    uint32_t total_stab = n_stab_per_type << 1U;
    int64_t rate_q16    = 0;
    if (total_stab == 8U) {
        rate_q16 = (int64_t)defects << 13U;
    } else if (total_stab == 24U) {
        rate_q16 = (int64_t)defects * INT64_C(2731); /* ~65536 / 24 */
    } else if (total_stab == 48U) {
        rate_q16 = (int64_t)defects * INT64_C(1365); /* ~65536 / 48 */
    } else {
        rate_q16 = (int64_t)defects * INT64_C(2731);
    }

    /* EWMA alpha = 1/8 = 3-bit right shift */
    mon->ewma_rate_q16 = mon->ewma_rate_q16 + flow_ashr_i64(rate_q16 - mon->ewma_rate_q16, 3U);
    /* CUSUM drift allowance: 1.5x baseline + 0.005 slack (~328 in Q16.16) */
    int64_t allow_q16  = mon->base_rate_q16 + flow_ashr_i64(mon->base_rate_q16, 1U) + INT64_C(328);
    int64_t inc_q16    = mon->ewma_rate_q16 - allow_q16;
    int64_t nxt_c      = mon->cusum_q16 + inc_q16;
    mon->cusum_q16     = flow_sel_i64(nxt_c > 0, nxt_c, 0);
    mon->rounds += 1U;

    if (mon->cusum_q16 >= mon->thresh_q16) {
        mon->cusum_q16 = 0;
        mon->triggers += 1U;
        return true;
    }
    return false;
}

#ifdef __cplusplus
}
#endif

#endif /* QCTRL_DECODER_RING0_H */
