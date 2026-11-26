#ifndef QCTRL_DECODER3D_RING0_H
#define QCTRL_DECODER3D_RING0_H

/* ============================================================================
 * FLOW-QCtrl Ring-0 Pure-GPR 3D Space-Time Surface Code Decoder (d = 3, 5, 7)
 * ============================================================================
 * Extends `qctrl_decoder_ring0.h` to phenomenological 3D space-time decoding
 * across `d` noisy measurement rounds + 1 final ideal readout slice (`<= 8`
 * slices total), remaining 100% Pure-GPR freestanding:
 *   - Zero FPU / SIMD instructions (-mgeneral-regs-only clean)
 *   - Zero udiv / sdiv / libc calls / heap / indirect jump tables
 *   - Exact Safe-Component Partition (`d3d(i,j) < bdist(i) + bdist(j)`)
 *   - Exact Bitmask-DP MWPM per component up to k <= 10 (<= 1024 states)
 *   - Pure-GPR Savings Greedy + 2-Opt Local Swap fallback for dense k > 10
 * ============================================================================ */

#include "qctrl_decoder_ring0.h"

#ifdef __cplusplus
extern "C" {
#endif

#define QCTRL_DEC3D_MAX_ROUNDS 8U
#define QCTRL_DEC3D_MAX_DETS   192U
#define QCTRL_DEC3D_DP_MAX_K   10U
#define QCTRL_DEC3D_MATE_BND   0xFFU

FLOW_INL void qctrl_sc_detectors3d_ring0(const uint32_t *meas_synd, uint32_t n_slices,
                                         uint32_t *out_det) {
    uint32_t t;
    if (!meas_synd || !out_det || n_slices == 0U) {
        return;
    }
    out_det[0] = meas_synd[0];
    for (t = 1U; t < n_slices; ++t) {
        out_det[t] = meas_synd[t] ^ meas_synd[t - 1U];
    }
}

FLOW_INL uint16_t qctrl__dist3d_ring0(const qctrl_sc_code_t *code, uint32_t g,
                                      uint8_t si, uint8_t ti, uint8_t sj, uint8_t tj) {
    uint16_t dt = (ti >= tj) ? (uint16_t)(ti - tj) : (uint16_t)(tj - ti);
    return (uint16_t)((uint16_t)code->pair_dist[g][si][sj] + dt);
}

FLOW_INL uint64_t qctrl__solve_comp_dp3d_ring0(const qctrl_sc_code_t *code, uint32_t g,
                                               const uint8_t *cs, const uint8_t *ct, uint32_t k) {
    uint32_t n_states = UINT32_C(1) << k;
    uint16_t dp[1024];
    uint16_t choice[1024];
    uint32_t m;
    dp[0]     = 0U;
    choice[0] = 0xFFFFU;
    for (m = 1U; m < n_states; ++m) {
        uint32_t i      = qctrl__ctz64_gpr((uint64_t)m);
        uint32_t rest_i = m ^ (UINT32_C(1) << i);
        uint8_t si      = cs[i];
        uint8_t ti      = ct[i];
        uint16_t best   = (uint16_t)((uint16_t)code->bdist[g][si] + dp[rest_i]);
        uint16_t best_c = (uint16_t)i;
        uint32_t j;
        for (j = i + 1U; j < k; ++j) {
            if ((rest_i & (UINT32_C(1) << j)) != 0U) {
                uint32_t rest_ij = rest_i ^ (UINT32_C(1) << j);
                uint16_t d3d     = qctrl__dist3d_ring0(code, g, si, ti, cs[j], ct[j]);
                uint16_t cand    = (uint16_t)(d3d + dp[rest_ij]);
                if (cand < best) {
                    best   = cand;
                    best_c = (uint16_t)(UINT32_C(0x8000) | (i << 4U) | j);
                }
            }
        }
        dp[m]     = best;
        choice[m] = best_c;
    }
    uint64_t corr = 0ULL;
    uint32_t cur  = n_states - 1U;
    while (cur != 0U) {
        uint16_t ch = choice[cur];
        if ((ch & 0x8000U) == 0U) {
            uint32_t i = (uint32_t)(ch & 0xFU);
            corr ^= code->bpath_mask[g][cs[i]];
            cur  ^= (UINT32_C(1) << i);
        } else {
            uint32_t i = (uint32_t)((ch >> 4U) & 0xFU);
            uint32_t j = (uint32_t)(ch & 0xFU);
            corr ^= code->pair_mask[g][cs[i]][cs[j]];
            cur  ^= (UINT32_C(1) << i) | (UINT32_C(1) << j);
        }
    }
    return corr;
}

FLOW_INL uint16_t qctrl__mate_cost3d(const qctrl_sc_code_t *code, uint32_t g,
                                     const uint8_t *cs, const uint8_t *ct,
                                     const uint8_t *mate, uint32_t u) {
    uint8_t mv = mate[u];
    if (mv == QCTRL_DEC3D_MATE_BND) {
        return (uint16_t)code->bdist[g][cs[u]];
    }
    return qctrl__dist3d_ring0(code, g, cs[u], ct[u], cs[mv], ct[mv]);
}

FLOW_INL uint64_t qctrl__solve_comp_2opt3d_ring0(const qctrl_sc_code_t *code, uint32_t g,
                                                 const uint8_t *cs, const uint8_t *ct, uint32_t k) {
    uint8_t mate[QCTRL_DEC3D_MAX_DETS];
    uint32_t i, j, pass;
    for (i = 0U; i < k; ++i) {
        mate[i] = QCTRL_DEC3D_MATE_BND;
    }

    /* Phase 1: Max-savings greedy initialization over pairs with d3d(i,j) < bdist(i)+bdist(j) */
    for (;;) {
        int32_t best_gain = 0;
        uint32_t best_i = 0U, best_j = 0U;
        for (i = 0U; i < k; ++i) {
            if (mate[i] != QCTRL_DEC3D_MATE_BND) { continue; }
            uint16_t bi = (uint16_t)code->bdist[g][cs[i]];
            for (j = i + 1U; j < k; ++j) {
                if (mate[j] != QCTRL_DEC3D_MATE_BND) { continue; }
                uint16_t pw   = qctrl__dist3d_ring0(code, g, cs[i], ct[i], cs[j], ct[j]);
                uint16_t bw   = (uint16_t)(bi + (uint16_t)code->bdist[g][cs[j]]);
                int32_t gain  = (int32_t)bw - (int32_t)pw;
                if (gain > best_gain) {
                    best_gain = gain;
                    best_i    = i;
                    best_j    = j;
                }
            }
        }
        if (best_gain <= 0) {
            break;
        }
        mate[best_i] = (uint8_t)best_j;
        mate[best_j] = (uint8_t)best_i;
    }

    /* Phase 2: Pure-GPR 2-Opt local refinement (swap pairs or boundary reassignment) */
    for (pass = 0U; pass < 4U; ++pass) {
        bool improved = false;
        for (i = 0U; i < k; ++i) {
            for (j = i + 1U; j < k; ++j) {
                if (mate[i] == (uint8_t)j) {
                    uint16_t w_cur = qctrl__dist3d_ring0(code, g, cs[i], ct[i], cs[j], ct[j]);
                    uint16_t w_bnd = (uint16_t)((uint16_t)code->bdist[g][cs[i]] +
                                                (uint16_t)code->bdist[g][cs[j]]);
                    if (w_bnd < w_cur) {
                        mate[i] = QCTRL_DEC3D_MATE_BND;
                        mate[j] = QCTRL_DEC3D_MATE_BND;
                        improved = true;
                    }
                    continue;
                }
                uint8_t mi = mate[i];
                uint8_t mj = mate[j];
                uint16_t cur_w = 0U;
                cur_w = (uint16_t)(cur_w + ((mi == QCTRL_DEC3D_MATE_BND)
                            ? (uint16_t)code->bdist[g][cs[i]]
                            : qctrl__dist3d_ring0(code, g, cs[i], ct[i], cs[mi], ct[mi])));
                cur_w = (uint16_t)(cur_w + ((mj == QCTRL_DEC3D_MATE_BND)
                            ? (uint16_t)code->bdist[g][cs[j]]
                            : qctrl__dist3d_ring0(code, g, cs[j], ct[j], cs[mj], ct[mj])));
                uint16_t w_ij = qctrl__dist3d_ring0(code, g, cs[i], ct[i], cs[j], ct[j]);
                uint16_t w_mimj = 0U;
                if (mi != QCTRL_DEC3D_MATE_BND && mj != QCTRL_DEC3D_MATE_BND) {
                    w_mimj = qctrl__dist3d_ring0(code, g, cs[mi], ct[mi], cs[mj], ct[mj]);
                    uint16_t b_mimj = (uint16_t)((uint16_t)code->bdist[g][cs[mi]] +
                                                 (uint16_t)code->bdist[g][cs[mj]]);
                    if (b_mimj < w_mimj) {
                        /* Reconnect i<->j and split mi, mj to boundary if cheaper */
                        if ((uint16_t)(w_ij + b_mimj) < cur_w) {
                            mate[i]  = (uint8_t)j;
                            mate[j]  = (uint8_t)i;
                            mate[mi] = QCTRL_DEC3D_MATE_BND;
                            mate[mj] = QCTRL_DEC3D_MATE_BND;
                            improved = true;
                            continue;
                        }
                    }
                    if ((uint16_t)(w_ij + w_mimj) < cur_w) {
                        mate[i]  = (uint8_t)j;
                        mate[j]  = (uint8_t)i;
                        mate[mi] = mj;
                        mate[mj] = mi;
                        improved = true;
                    }
                } else if (mi != QCTRL_DEC3D_MATE_BND && mj == QCTRL_DEC3D_MATE_BND) {
                    uint16_t new_w = (uint16_t)(w_ij + (uint16_t)code->bdist[g][cs[mi]]);
                    if (new_w < cur_w) {
                        mate[i]  = (uint8_t)j;
                        mate[j]  = (uint8_t)i;
                        mate[mi] = QCTRL_DEC3D_MATE_BND;
                        improved = true;
                    }
                } else if (mi == QCTRL_DEC3D_MATE_BND && mj != QCTRL_DEC3D_MATE_BND) {
                    uint16_t new_w = (uint16_t)(w_ij + (uint16_t)code->bdist[g][cs[mj]]);
                    if (new_w < cur_w) {
                        mate[i]  = (uint8_t)j;
                        mate[j]  = (uint8_t)i;
                        mate[mj] = QCTRL_DEC3D_MATE_BND;
                        improved = true;
                    }
                }
            }
        }
        if (!improved) {
            break;
        }
    }

    uint64_t corr = 0ULL;
    for (i = 0U; i < k; ++i) {
        if (mate[i] == QCTRL_DEC3D_MATE_BND) {
            corr ^= code->bpath_mask[g][cs[i]];
        } else if (i < (uint32_t)mate[i]) {
            corr ^= code->pair_mask[g][cs[i]][cs[mate[i]]];
        }
    }
    return corr;
}

FLOW_INL uint64_t qctrl_sc_decode3d_ring0(const qctrl_sc_code_t *code, const uint32_t *det,
                                          uint32_t n_slices, bool is_z_err) {
    uint8_t det_s[QCTRL_DEC3D_MAX_DETS];
    uint8_t det_t[QCTRL_DEC3D_MAX_DETS];
    uint32_t k = 0U;
    uint32_t t;
    uint32_t g = is_z_err ? 1U : 0U;

    if (!code || !det || n_slices == 0U) {
        return 0ULL;
    }
    if (n_slices > QCTRL_DEC3D_MAX_ROUNDS) {
        n_slices = QCTRL_DEC3D_MAX_ROUNDS;
    }

    for (t = 0U; t < n_slices; ++t) {
        uint32_t word = det[t];
        while (word != 0U) {
            uint32_t s = qctrl__ctz64_gpr((uint64_t)word);
            word &= (word - 1U);
            if (s < code->n_stab && k < QCTRL_DEC3D_MAX_DETS) {
                det_s[k] = (uint8_t)s;
                det_t[k] = (uint8_t)t;
                k += 1U;
            }
        }
    }
    if (k == 0U) {
        return 0ULL;
    }
    if (k <= QCTRL_DEC3D_DP_MAX_K) {
        return qctrl__solve_comp_dp3d_ring0(code, g, det_s, det_t, k);
    }

    /* Exact Safe-Component BFS partitioning: edge (i,j) iff d3d(i,j) < bdist(i) + bdist(j) */
    uint8_t vis[QCTRL_DEC3D_MAX_DETS];
    uint32_t i;
    for (i = 0U; i < k; ++i) {
        vis[i] = 0U;
    }
    uint64_t corr = 0ULL;
    for (i = 0U; i < k; ++i) {
        if (vis[i] != 0U) { continue; }
        uint8_t queue[QCTRL_DEC3D_MAX_DETS];
        uint8_t cs[QCTRL_DEC3D_MAX_DETS];
        uint8_t ct[QCTRL_DEC3D_MAX_DETS];
        uint32_t qh = 0U, qt = 0U;
        vis[i]      = 1U;
        queue[qt++] = (uint8_t)i;
        while (qh < qt) {
            uint32_t u = (uint32_t)queue[qh++];
            uint8_t su = det_s[u];
            uint8_t tu = det_t[u];
            uint16_t bu = (uint16_t)code->bdist[g][su];
            cs[qh - 1U] = su;
            ct[qh - 1U] = tu;
            uint32_t v;
            for (v = 0U; v < k; ++v) {
                if (vis[v] != 0U) { continue; }
                uint8_t sv  = det_s[v];
                uint8_t tv  = det_t[v];
                uint16_t d3 = qctrl__dist3d_ring0(code, g, su, tu, sv, tv);
                if (d3 < (uint16_t)(bu + (uint16_t)code->bdist[g][sv])) {
                    vis[v]      = 1U;
                    queue[qt++] = (uint8_t)v;
                }
            }
        }
        if (qt == 1U) {
            corr ^= code->bpath_mask[g][cs[0]];
        } else if (qt <= QCTRL_DEC3D_DP_MAX_K) {
            corr ^= qctrl__solve_comp_dp3d_ring0(code, g, cs, ct, qt);
        } else {
            corr ^= qctrl__solve_comp_2opt3d_ring0(code, g, cs, ct, qt);
        }
    }
    return corr;
}

#ifdef __cplusplus
}
#endif

#endif /* QCTRL_DECODER3D_RING0_H */
