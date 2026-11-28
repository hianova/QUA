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
 *   - Fast exact Bitmask-DP MWPM per component up to k <= 8 (<= 256 states)
 *   - Exact Iterative Integer Edmonds' Blossom MWM (k <= 48) via Savings
 *     Reduction (w2 = 2 * max(0, b_i + b_j - d3d(i,j)), exact bit-shift duals)
 *   - Pure-GPR Savings Greedy + 2-Opt fallback only for extreme bursts k > 48
 * ============================================================================ */

#include "qctrl_decoder_ring0.h"

#ifdef __cplusplus
extern "C" {
#endif

#define QCTRL_DEC3D_MAX_ROUNDS 8U
#define QCTRL_DEC3D_MAX_DETS   192U
#define QCTRL_DEC3D_DP_MAX_K   8U
#define QCTRL_DEC3D_MATE_BND   0xFFU

#define QCTRL_BLOS_MAX_N       48U
#define QCTRL_BLOS3D_MAX_K     QCTRL_BLOS_MAX_N
#define QCTRL_BLOS_MAX_V       74U

#if defined(__GNUC__) || defined(__clang__)
#define QCTRL_DEC3D_NOINL static __attribute__((noinline, unused))
#else
#define QCTRL_DEC3D_NOINL static
#endif

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

FLOW_INL uint64_t qctrl__solve_comp_dp3d_ring0_ex(const qctrl_sc_code_t *code, uint32_t g,
                                                  const uint8_t *cs, const uint8_t *ct,
                                                  uint32_t k, uint32_t *out_cost) {
    uint32_t n_states = UINT32_C(1) << k;
    uint16_t dp[256];
    uint16_t choice[256];
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
    if (out_cost != NULL) {
        *out_cost = (uint32_t)dp[n_states - 1U];
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

FLOW_INL uint64_t qctrl__solve_comp_dp3d_ring0(const qctrl_sc_code_t *code, uint32_t g,
                                               const uint8_t *cs, const uint8_t *ct, uint32_t k) {
    return qctrl__solve_comp_dp3d_ring0_ex(code, g, cs, ct, k, NULL);
}

/* ============================================================================
 * Iterative Integer Edmonds' Blossom MWM (Pure-GPR, Zero Recursion, Zero Heap)
 * ============================================================================ */

typedef struct {
    uint8_t  g_u[QCTRL_BLOS_MAX_V][QCTRL_BLOS_MAX_V];
    uint8_t  g_v[QCTRL_BLOS_MAX_V][QCTRL_BLOS_MAX_V];
    uint8_t  g_w[QCTRL_BLOS_MAX_V][QCTRL_BLOS_MAX_V];
    uint8_t  flower_from[QCTRL_BLOS_MAX_V][QCTRL_BLOS_MAX_N + 1U];
    uint8_t  flower[QCTRL_BLOS_MAX_V][QCTRL_BLOS_MAX_V];
    uint8_t  flen[QCTRL_BLOS_MAX_V];
    int32_t  lab[QCTRL_BLOS_MAX_V];
    uint8_t  match[QCTRL_BLOS_MAX_V];
    uint8_t  slack[QCTRL_BLOS_MAX_V];
    uint8_t  st[QCTRL_BLOS_MAX_V];
    uint8_t  pa[QCTRL_BLOS_MAX_V];
    int8_t   S[QCTRL_BLOS_MAX_V];
    uint32_t vis[QCTRL_BLOS_MAX_V];
    uint8_t  q[256];
    uint32_t n;
    uint32_t n_x;
    uint32_t vis_timer;
    uint32_t qh;
    uint32_t qt;
} qctrl_blos3d_ws_t;

FLOW_INL int32_t qctrl__blos_edelta(const qctrl_blos3d_ws_t *ws, uint32_t x, uint32_t y) {
    return ws->lab[ws->g_u[x][y]] + ws->lab[ws->g_v[x][y]] - (int32_t)ws->g_w[x][y];
}

FLOW_INL void qctrl__blos_update_slack(qctrl_blos3d_ws_t *ws, uint32_t u, uint32_t x) {
    if (ws->slack[x] == 0U ||
        qctrl__blos_edelta(ws, u, x) < qctrl__blos_edelta(ws, (uint32_t)ws->slack[x], x)) {
        ws->slack[x] = (uint8_t)u;
    }
}

FLOW_INL void qctrl__blos_set_slack(qctrl_blos3d_ws_t *ws, uint32_t x) {
    uint32_t u;
    ws->slack[x] = 0U;
    for (u = 1U; u <= ws->n; ++u) {
        if (ws->g_w[u][x] > 0U && ws->st[u] != (uint8_t)x && ws->S[ws->st[u]] == 0) {
            qctrl__blos_update_slack(ws, u, x);
        }
    }
}

FLOW_INL void qctrl__blos_qpush(qctrl_blos3d_ws_t *ws, uint32_t root) {
    uint8_t stk[QCTRL_BLOS_MAX_V * 2U];
    uint32_t sp = 0U;
    stk[sp++]   = (uint8_t)root;
    while (sp > 0U) {
        uint32_t x = (uint32_t)stk[--sp];
        if (x <= ws->n) {
            if (ws->qt < 256U) {
                ws->q[ws->qt++] = (uint8_t)x;
            }
        } else {
            uint32_t len = (uint32_t)ws->flen[x];
            uint32_t i;
            /* Push in reverse so children leave stack in forward cycle order */
            for (i = len; i > 0U; --i) {
                if (sp < QCTRL_BLOS_MAX_V * 2U) {
                    stk[sp++] = ws->flower[x][i - 1U];
                }
            }
        }
    }
}

FLOW_INL void qctrl__blos_set_st(qctrl_blos3d_ws_t *ws, uint32_t root, uint8_t b) {
    uint8_t stk[QCTRL_BLOS_MAX_V * 2U];
    uint32_t sp = 0U;
    stk[sp++]   = (uint8_t)root;
    while (sp > 0U) {
        uint32_t x = (uint32_t)stk[--sp];
        ws->st[x]  = b;
        if (x > ws->n) {
            uint32_t len = (uint32_t)ws->flen[x];
            uint32_t i;
            for (i = 0U; i < len; ++i) {
                if (sp < QCTRL_BLOS_MAX_V * 2U) {
                    stk[sp++] = ws->flower[x][i];
                }
            }
        }
    }
}

FLOW_INL uint32_t qctrl__blos_get_pr(qctrl_blos3d_ws_t *ws, uint32_t b, uint8_t xr) {
    uint32_t len = (uint32_t)ws->flen[b];
    uint32_t pr  = 0U;
    while (pr < len && ws->flower[b][pr] != xr) {
        pr += 1U;
    }
    if ((pr & 1U) != 0U && len > 1U) {
        uint32_t lo = 1U;
        uint32_t hi = len - 1U;
        while (lo < hi) {
            uint8_t tmp       = ws->flower[b][lo];
            ws->flower[b][lo] = ws->flower[b][hi];
            ws->flower[b][hi] = tmp;
            lo += 1U;
            hi -= 1U;
        }
        return len - pr;
    }
    return pr;
}

FLOW_INL void qctrl__blos_set_match(qctrl_blos3d_ws_t *ws, uint32_t u_in, uint32_t v_in) {
    uint8_t su_stk[QCTRL_BLOS_MAX_V * 2U];
    uint8_t sv_stk[QCTRL_BLOS_MAX_V * 2U];
    uint32_t sp = 0U;
    su_stk[sp]  = (uint8_t)u_in;
    sv_stk[sp]  = (uint8_t)v_in;
    sp += 1U;
    while (sp > 0U) {
        sp -= 1U;
        uint32_t u    = (uint32_t)su_stk[sp];
        uint32_t v    = (uint32_t)sv_stk[sp];
        uint8_t eu    = ws->g_u[u][v];
        uint8_t ev    = ws->g_v[u][v];
        ws->match[u]  = ev;
        if (u > ws->n) {
            uint8_t xr   = ws->flower_from[u][eu];
            uint32_t pr  = qctrl__blos_get_pr(ws, u, xr);
            uint32_t len = (uint32_t)ws->flen[u];
            uint32_t i;
            if (sp < QCTRL_BLOS_MAX_V * 2U) {
                su_stk[sp] = xr;
                sv_stk[sp] = (uint8_t)v;
                sp += 1U;
            }
            for (i = 0U; i < pr; ++i) {
                if (sp < QCTRL_BLOS_MAX_V * 2U) {
                    su_stk[sp] = ws->flower[u][i];
                    sv_stk[sp] = ws->flower[u][i ^ 1U];
                    sp += 1U;
                }
            }
            if (pr > 0U && pr < len) {
                uint8_t rot[QCTRL_BLOS_MAX_V];
                for (i = 0U; i < len; ++i) {
                    uint32_t src = i + pr;
                    if (src >= len) {
                        src -= len;
                    }
                    rot[i] = ws->flower[u][src];
                }
                for (i = 0U; i < len; ++i) {
                    ws->flower[u][i] = rot[i];
                }
            }
        }
    }
}

FLOW_INL void qctrl__blos_augment(qctrl_blos3d_ws_t *ws, uint32_t u, uint32_t v) {
    for (;;) {
        uint32_t xnv = (uint32_t)ws->st[ws->match[u]];
        qctrl__blos_set_match(ws, u, v);
        if (xnv == 0U) {
            return;
        }
        uint32_t p_xnv = (uint32_t)ws->st[ws->pa[xnv]];
        qctrl__blos_set_match(ws, xnv, p_xnv);
        u = p_xnv;
        v = xnv;
    }
}

FLOW_INL uint32_t qctrl__blos_lca(qctrl_blos3d_ws_t *ws, uint32_t u, uint32_t v) {
    ws->vis_timer += 1U;
    uint32_t tm = ws->vis_timer;
    while (u != 0U || v != 0U) {
        if (u != 0U) {
            if (ws->vis[u] == tm) {
                return u;
            }
            ws->vis[u] = tm;
            u          = (uint32_t)ws->st[ws->match[u]];
            if (u != 0U) {
                u = (uint32_t)ws->st[ws->pa[u]];
            }
        }
        {
            uint32_t tmp = u;
            u            = v;
            v            = tmp;
        }
    }
    return 0U;
}

FLOW_INL void qctrl__blos_add_blossom(qctrl_blos3d_ws_t *ws, uint32_t u, uint32_t lca, uint32_t v) {
    uint32_t b = ws->n + 1U;
    while (b <= ws->n_x && ws->st[b] != 0U) {
        b += 1U;
    }
    if (b > ws->n_x) {
        if (ws->n_x + 1U >= QCTRL_BLOS_MAX_V) {
            return;
        }
        ws->n_x += 1U;
        b = ws->n_x;
    }

    ws->lab[b]   = 0;
    ws->S[b]     = 0;
    ws->match[b] = ws->match[lca];
    ws->pa[b]    = 0U;
    ws->slack[b] = 0U;
    ws->vis[b]   = 0U;
    ws->flen[b]  = 0U;
    ws->flower[b][ws->flen[b]++] = (uint8_t)lca;

    uint32_t x = u;
    while (x != lca && (uint32_t)ws->flen[b] + 2U < QCTRL_BLOS_MAX_V) {
        uint32_t y                   = (uint32_t)ws->st[ws->match[x]];
        ws->flower[b][ws->flen[b]++] = (uint8_t)x;
        ws->flower[b][ws->flen[b]++] = (uint8_t)y;
        qctrl__blos_qpush(ws, y);
        x = (uint32_t)ws->st[ws->pa[y]];
    }
    /* Reverse flower[b][1 .. flen[b]-1] */
    if (ws->flen[b] > 2U) {
        uint32_t lo = 1U;
        uint32_t hi = (uint32_t)ws->flen[b] - 1U;
        while (lo < hi) {
            uint8_t tmp       = ws->flower[b][lo];
            ws->flower[b][lo] = ws->flower[b][hi];
            ws->flower[b][hi] = tmp;
            lo += 1U;
            hi -= 1U;
        }
    }
    x = v;
    while (x != lca && (uint32_t)ws->flen[b] + 2U < QCTRL_BLOS_MAX_V) {
        uint32_t y                   = (uint32_t)ws->st[ws->match[x]];
        ws->flower[b][ws->flen[b]++] = (uint8_t)x;
        ws->flower[b][ws->flen[b]++] = (uint8_t)y;
        qctrl__blos_qpush(ws, y);
        x = (uint32_t)ws->st[ws->pa[y]];
    }

    qctrl__blos_set_st(ws, b, (uint8_t)b);

    uint32_t i;
    for (i = 1U; i <= ws->n_x; ++i) {
        ws->g_u[b][i] = (uint8_t)b;
        ws->g_v[b][i] = (uint8_t)i;
        ws->g_w[b][i] = 0U;
        ws->g_u[i][b] = (uint8_t)i;
        ws->g_v[i][b] = (uint8_t)b;
        ws->g_w[i][b] = 0U;
    }
    for (i = 1U; i <= ws->n; ++i) {
        ws->flower_from[b][i] = 0U;
    }

    uint32_t len = (uint32_t)ws->flen[b];
    for (i = 0U; i < len; ++i) {
        uint32_t xs = (uint32_t)ws->flower[b][i];
        uint32_t j;
        for (j = 1U; j <= ws->n_x; ++j) {
            if (ws->g_w[xs][j] > 0U &&
                (ws->g_w[b][j] == 0U ||
                 qctrl__blos_edelta(ws, xs, j) < qctrl__blos_edelta(ws, b, j))) {
                ws->g_u[b][j] = ws->g_u[xs][j];
                ws->g_v[b][j] = ws->g_v[xs][j];
                ws->g_w[b][j] = ws->g_w[xs][j];
                ws->g_u[j][b] = ws->g_u[j][xs];
                ws->g_v[j][b] = ws->g_v[j][xs];
                ws->g_w[j][b] = ws->g_w[j][xs];
            }
        }
        for (j = 1U; j <= ws->n; ++j) {
            if (ws->flower_from[xs][j] != 0U) {
                ws->flower_from[b][j] = (uint8_t)xs;
            }
        }
    }
    qctrl__blos_set_slack(ws, b);
}

FLOW_INL void qctrl__blos_expand_blossom(qctrl_blos3d_ws_t *ws, uint32_t b) {
    uint32_t len = (uint32_t)ws->flen[b];
    uint32_t i;
    for (i = 0U; i < len; ++i) {
        uint32_t c = (uint32_t)ws->flower[b][i];
        qctrl__blos_set_st(ws, c, (uint8_t)c);
    }
    uint8_t xr  = ws->flower_from[b][ws->g_u[b][ws->pa[b]]];
    uint32_t pr = qctrl__blos_get_pr(ws, b, xr);
    for (i = 0U; i < pr; i += 2U) {
        uint32_t xs    = (uint32_t)ws->flower[b][i];
        uint32_t xns   = (uint32_t)ws->flower[b][i + 1U];
        ws->pa[xs]     = ws->g_u[xns][xs];
        ws->S[xs]      = 1;
        ws->S[xns]     = 0;
        ws->slack[xs]  = 0U;
        qctrl__blos_set_slack(ws, xns);
        qctrl__blos_qpush(ws, xns);
    }
    ws->S[xr]  = 1;
    ws->pa[xr] = ws->pa[b];
    for (i = pr + 1U; i < len; ++i) {
        uint32_t xs = (uint32_t)ws->flower[b][i];
        ws->S[xs]   = -1;
        qctrl__blos_set_slack(ws, xs);
    }
    ws->st[b] = 0U;
}

FLOW_INL bool qctrl__blos_on_found_edge(qctrl_blos3d_ws_t *ws, uint8_t eu, uint8_t ev) {
    uint32_t u = (uint32_t)ws->st[eu];
    uint32_t v = (uint32_t)ws->st[ev];
    if (ws->S[v] == -1) {
        ws->pa[v] = eu;
        ws->S[v]  = 1;
        uint32_t nu = (uint32_t)ws->st[ws->match[v]];
        ws->slack[v]  = 0U;
        ws->slack[nu] = 0U;
        ws->S[nu]     = 0;
        qctrl__blos_qpush(ws, nu);
    } else if (ws->S[v] == 0) {
        uint32_t lca = qctrl__blos_lca(ws, u, v);
        if (lca == 0U) {
            qctrl__blos_augment(ws, u, v);
            qctrl__blos_augment(ws, v, u);
            return true;
        } else {
            qctrl__blos_add_blossom(ws, u, lca, v);
        }
    }
    return false;
}

FLOW_INL bool qctrl__blos_matching_step(qctrl_blos3d_ws_t *ws) {
    uint32_t i;
    for (i = 1U; i <= ws->n_x; ++i) {
        ws->S[i]     = -1;
        ws->slack[i] = 0U;
    }
    ws->qh = 0U;
    ws->qt = 0U;
    for (i = 1U; i <= ws->n_x; ++i) {
        if (ws->st[i] == (uint8_t)i && ws->match[i] == 0U) {
            ws->pa[i] = 0U;
            ws->S[i]  = 0;
            qctrl__blos_qpush(ws, i);
        }
    }
    if (ws->qt == 0U) {
        return false;
    }

    uint32_t iter;
    for (iter = 0U; iter < 512U; ++iter) {
        while (ws->qh < ws->qt) {
            uint32_t u = (uint32_t)ws->q[ws->qh++];
            if (ws->S[ws->st[u]] == 1) {
                continue;
            }
            uint32_t v;
            for (v = 1U; v <= ws->n; ++v) {
                if (ws->g_w[u][v] > 0U && ws->st[u] != ws->st[v]) {
                    if (qctrl__blos_edelta(ws, u, v) == 0) {
                        if (qctrl__blos_on_found_edge(ws, ws->g_u[u][v], ws->g_v[u][v])) {
                            return true;
                        }
                    } else {
                        qctrl__blos_update_slack(ws, u, (uint32_t)ws->st[v]);
                    }
                }
            }
        }
        ws->qh = 0U;
        ws->qt = 0U;

        int32_t d = INT32_C(0x3FFFFFFF);
        uint32_t u, b, x;
        for (u = 1U; u <= ws->n; ++u) {
            if (ws->S[ws->st[u]] == 0 && ws->lab[u] < d) {
                d = ws->lab[u];
            }
        }
        for (b = ws->n + 1U; b <= ws->n_x; ++b) {
            if (ws->st[b] == (uint8_t)b && ws->S[b] == 1) {
                int32_t half_b = (int32_t)((uint32_t)ws->lab[b] >> 1U);
                if (half_b < d) {
                    d = half_b;
                }
            }
        }
        for (x = 1U; x <= ws->n_x; ++x) {
            if (ws->st[x] == (uint8_t)x && ws->slack[x] != 0U) {
                int32_t ed = qctrl__blos_edelta(ws, (uint32_t)ws->slack[x], x);
                if (ws->S[x] == -1) {
                    if (ed < d) {
                        d = ed;
                    }
                } else if (ws->S[x] == 0) {
                    int32_t half_ed = (int32_t)((uint32_t)ed >> 1U);
                    if (half_ed < d) {
                        d = half_ed;
                    }
                }
            }
        }

        for (u = 1U; u <= ws->n; ++u) {
            if (ws->S[ws->st[u]] == 0) {
                if (ws->lab[u] <= d) {
                    return false;
                }
                ws->lab[u] -= d;
            } else if (ws->S[ws->st[u]] == 1) {
                ws->lab[u] += d;
            }
        }
        for (b = ws->n + 1U; b <= ws->n_x; ++b) {
            if (ws->st[b] == (uint8_t)b) {
                if (ws->S[b] == 0) {
                    ws->lab[b] += (d << 1);
                } else if (ws->S[b] == 1) {
                    ws->lab[b] -= (d << 1);
                }
            }
        }

        for (x = 1U; x <= ws->n_x; ++x) {
            if (ws->st[x] == (uint8_t)x && ws->slack[x] != 0U &&
                ws->st[ws->slack[x]] != (uint8_t)x &&
                qctrl__blos_edelta(ws, (uint32_t)ws->slack[x], x) == 0) {
                if (qctrl__blos_on_found_edge(ws, ws->g_u[ws->slack[x]][x],
                                              ws->g_v[ws->slack[x]][x])) {
                    return true;
                }
            }
        }
        for (b = ws->n + 1U; b <= ws->n_x; ++b) {
            if (ws->st[b] == (uint8_t)b && ws->S[b] == 1 && ws->lab[b] == 0) {
                qctrl__blos_expand_blossom(ws, b);
            }
        }
    }
    return false;
}

QCTRL_DEC3D_NOINL uint64_t qctrl__solve_comp_blossom3d_ring0(const qctrl_sc_code_t *code, uint32_t g,
                                                             const uint8_t *cs, const uint8_t *ct,
                                                             uint32_t k, uint32_t *out_cost) {
    if (!code || !cs || !ct || k == 0U || k > QCTRL_BLOS3D_MAX_K) {
        if (out_cost != NULL) {
            *out_cost = 0U;
        }
        return 0ULL;
    }

    qctrl_blos3d_ws_t ws;
    uint32_t i, j;
    uint8_t w_max  = 0U;

    ws.n         = k;
    ws.n_x       = k;
    ws.vis_timer = 0U;
    ws.qh        = 0U;
    ws.qt        = 0U;

    for (i = 0U; i < (QCTRL_BLOS_MAX_N + 1U); ++i) {
        ws.lab[i]   = 0;
        ws.match[i] = 0U;
        ws.slack[i] = 0U;
        ws.st[i]    = (uint8_t)i;
        ws.pa[i]    = 0U;
        ws.S[i]     = -1;
        ws.vis[i]   = 0U;
        ws.flen[i]  = 0U;
        for (j = 0U; j < (QCTRL_BLOS_MAX_N + 1U); ++j) {
            ws.g_u[i][j] = (uint8_t)i;
            ws.g_v[i][j] = (uint8_t)j;
            ws.g_w[i][j] = 0U;
        }
        for (j = 0U; j < (QCTRL_BLOS_MAX_N + 1U); ++j) {
            ws.flower_from[i][j] = (i > 0U && i == j) ? (uint8_t)i : 0U;
        }
    }

    for (i = 0U; i < k; ++i) {
        uint16_t bi = (uint16_t)code->bdist[g][cs[i]];
        for (j = i + 1U; j < k; ++j) {
            uint16_t d3 = qctrl__dist3d_ring0(code, g, cs[i], ct[i], cs[j], ct[j]);
            uint16_t bw = (uint16_t)(bi + (uint16_t)code->bdist[g][cs[j]]);
            if (d3 < bw) {
                uint32_t gain = (uint32_t)(bw - d3);
                uint32_t w2   = gain << 1U;
                uint8_t w2_u8 = (w2 > 254U) ? 254U : (uint8_t)w2;
                ws.g_w[i + 1U][j + 1U] = w2_u8;
                ws.g_w[j + 1U][i + 1U] = w2_u8;
                if (w2_u8 > w_max) {
                    w_max = w2_u8;
                }
            }
        }
    }

    if (w_max == 0U) {
        uint64_t corr0 = 0ULL;
        uint32_t cost0 = 0U;
        for (i = 0U; i < k; ++i) {
            corr0 ^= code->bpath_mask[g][cs[i]];
            cost0 += (uint32_t)code->bdist[g][cs[i]];
        }
        if (out_cost != NULL) {
            *out_cost = cost0;
        }
        return corr0;
    }

    for (i = 1U; i <= QCTRL_BLOS_MAX_N; ++i) {
        if (i <= k) {
            ws.lab[i] = (int32_t)w_max;
        }
    }

    while (qctrl__blos_matching_step(&ws)) {
    }

    uint64_t corr = 0ULL;
    uint32_t cost = 0U;
    for (i = 0U; i < k; ++i) {
        uint32_t m1 = (uint32_t)ws.match[i + 1U];
        if (m1 == 0U) {
            corr ^= code->bpath_mask[g][cs[i]];
            cost += (uint32_t)code->bdist[g][cs[i]];
        } else {
            uint32_t mj = m1 - 1U;
            if (i < mj) {
                corr ^= code->pair_mask[g][cs[i]][cs[mj]];
                cost += (uint32_t)qctrl__dist3d_ring0(code, g, cs[i], ct[i], cs[mj], ct[mj]);
            }
        }
    }
    if (out_cost != NULL) {
        *out_cost = cost;
    }
    return corr;
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
                uint16_t pw  = qctrl__dist3d_ring0(code, g, cs[i], ct[i], cs[j], ct[j]);
                uint16_t bw  = (uint16_t)(bi + (uint16_t)code->bdist[g][cs[j]]);
                int32_t gain = (int32_t)bw - (int32_t)pw;
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
                        mate[i]  = QCTRL_DEC3D_MATE_BND;
                        mate[j]  = QCTRL_DEC3D_MATE_BND;
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
                uint16_t w_ij   = qctrl__dist3d_ring0(code, g, cs[i], ct[i], cs[j], ct[j]);
                uint16_t w_mimj = 0U;
                if (mi != QCTRL_DEC3D_MATE_BND && mj != QCTRL_DEC3D_MATE_BND) {
                    w_mimj = qctrl__dist3d_ring0(code, g, cs[mi], ct[mi], cs[mj], ct[mj]);
                    uint16_t b_mimj = (uint16_t)((uint16_t)code->bdist[g][cs[mi]] +
                                                 (uint16_t)code->bdist[g][cs[mj]]);
                    if (b_mimj < w_mimj) {
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
    if (k == 1U) {
        return code->bpath_mask[g][det_s[0]];
    }
    if (k == 2U) {
        uint16_t d3 = qctrl__dist3d_ring0(code, g, det_s[0], det_t[0], det_s[1], det_t[1]);
        uint16_t bw = (uint16_t)((uint16_t)code->bdist[g][det_s[0]] +
                                 (uint16_t)code->bdist[g][det_s[1]]);
        return (d3 < bw)
                   ? code->pair_mask[g][det_s[0]][det_s[1]]
                   : (code->bpath_mask[g][det_s[0]] ^ code->bpath_mask[g][det_s[1]]);
    }
    if (k <= 4U) {
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
            uint32_t u  = (uint32_t)queue[qh++];
            uint8_t su  = det_s[u];
            uint8_t tu  = det_t[u];
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
        } else if (qt == 2U) {
            uint16_t d3 = qctrl__dist3d_ring0(code, g, cs[0], ct[0], cs[1], ct[1]);
            uint16_t bw = (uint16_t)((uint16_t)code->bdist[g][cs[0]] +
                                     (uint16_t)code->bdist[g][cs[1]]);
            corr ^= (d3 < bw)
                        ? code->pair_mask[g][cs[0]][cs[1]]
                        : (code->bpath_mask[g][cs[0]] ^ code->bpath_mask[g][cs[1]]);
        } else if (qt <= QCTRL_DEC3D_DP_MAX_K) {
            corr ^= qctrl__solve_comp_dp3d_ring0(code, g, cs, ct, qt);
        } else if (qt <= QCTRL_BLOS_MAX_N) {
            corr ^= qctrl__solve_comp_blossom3d_ring0(code, g, cs, ct, qt, NULL);
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
