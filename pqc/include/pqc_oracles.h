/*
 * pqc_oracles.h — Declarative QIR PQC / NIST Cryptanalytic Oracles & Classical KAT References
 * Part of QUA / PQC Quantum Cryptanalysis Suite.
 * Pure C17 & C++20 compatible, zero <math.h>, zero raw qcc_kmb_* gate emission.
 */
#ifndef PQC_ORACLES_H
#define PQC_ORACLES_H

#include "qir.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*pqc_f2_affine_fn_t)(const uint8_t *in_bytes, uint8_t *out_bytes, const void *ctx);

static inline uint8_t pqc__get_bit(const uint8_t *buf, uint32_t idx) {
    return (uint8_t)((buf[idx >> 3u] >> (idx & 7u)) & 1u);
}

static inline void pqc__set_bit(uint8_t *buf, uint32_t idx, uint8_t bit) {
    if (bit & 1u) buf[idx >> 3u] |= (uint8_t)(1u << (idx & 7u));
    else          buf[idx >> 3u] &= (uint8_t)~(1u << (idx & 7u));
}

static inline uint32_t pqc__rotr32(uint32_t x, uint32_t n) {
    return (n == 0u) ? x : ((x >> n) | (x << (32u - n)));
}

static inline uint64_t pqc__rotr64(uint64_t x, uint32_t n) {
    return (n == 0u) ? x : ((x >> n) | (x << (64u - n)));
}

static inline uint16_t pqc__rotl16(uint16_t x, uint32_t n) {
    n &= 15u;
    return (n == 0u) ? x : (uint16_t)(((x << n) | (x >> (16u - n))) & 0xFFFFu);
}

static inline uint32_t pqc__load_be32(const uint8_t *p) {
    return ((uint32_t)p[0] << 24u) | ((uint32_t)p[1] << 16u) |
           ((uint32_t)p[2] << 8u)  |  (uint32_t)p[3];
}

static inline void pqc__store_be32(uint8_t *p, uint32_t v) {
    p[0] = (uint8_t)(v >> 24u);
    p[1] = (uint8_t)(v >> 16u);
    p[2] = (uint8_t)(v >> 8u);
    p[3] = (uint8_t)v;
}

static inline uint32_t pqc__load_le32(const uint8_t *p) {
    return ((uint32_t)p[0])       | ((uint32_t)p[1] << 8u) |
           ((uint32_t)p[2] << 16u) | ((uint32_t)p[3] << 24u);
}

static inline void pqc__store_le32(uint8_t *p, uint32_t v) {
    p[0] = (uint8_t)v;
    p[1] = (uint8_t)(v >> 8u);
    p[2] = (uint8_t)(v >> 16u);
    p[3] = (uint8_t)(v >> 24u);
}

static inline void pqc__unpack_u32_le(uint32_t x, uint8_t out4[4]) {
    pqc__store_le32(out4, x);
}

static inline uint32_t pqc__pack_u32_le(const uint8_t in4[4]) {
    return pqc__load_le32(in4);
}

static inline void pqc__unpack_u64_le(uint64_t x, uint8_t out8[8]) {
    for (uint32_t i = 0u; i < 8u; ++i) out8[i] = (uint8_t)(x >> (8u * i));
}

static inline uint64_t pqc__pack_u64_le(const uint8_t in8[8]) {
    uint64_t r = 0u;
    for (uint32_t i = 0u; i < 8u; ++i) r |= ((uint64_t)in8[i] << (8u * i));
    return r;
}

/* Build an exact QIR_OP_LINF2 node from any classical F2-affine callback (in_bits, out_bits <= 256) */
static inline uint32_t pqc_emit_linf2_from_fn(
    qir_module_t *m,
    const uint32_t *ins,
    uint16_t num_ins,
    uint16_t in_bits,
    uint16_t out_bits,
    pqc_f2_affine_fn_t fn,
    const void *ctx,
    const char *name
) {
    uint8_t in_b[32], base_out[32], cur_out[32];
    memset(in_b, 0, sizeof(in_b));
    memset(base_out, 0, sizeof(base_out));
    memset(cur_out, 0, sizeof(cur_out));
    fn(in_b, base_out, ctx);

    qcc_bv256_t *rows = (qcc_bv256_t *)calloc(out_bits ? out_bits : 1u, sizeof(qcc_bv256_t));
    qcc_bv256_t cvec = qcc_bv_zero();
    for (uint16_t r = 0u; r < out_bits; ++r) {
        if (pqc__get_bit(base_out, r)) qcc_bv_set(cvec.w, r);
    }
    for (uint16_t c = 0u; c < in_bits; ++c) {
        memset(in_b, 0, sizeof(in_b));
        memset(cur_out, 0, sizeof(cur_out));
        pqc__set_bit(in_b, c, 1u);
        fn(in_b, cur_out, ctx);
        for (uint16_t r = 0u; r < out_bits; ++r) {
            uint8_t bit = (uint8_t)(pqc__get_bit(cur_out, r) ^ pqc__get_bit(base_out, r));
            if (bit) qcc_bv_set(rows[r].w, c);
        }
    }
    uint32_t vid = qir_emit_linf2(m, ins, num_ins, out_bits, rows, &cvec, name);
    free(rows);
    return vid;
}

typedef struct {
    const uint8_t *target_bytes;
    uint16_t n_bits;
} pqc__xor_const_ctx_t;

static inline void pqc__fn_xor_const(const uint8_t *in_b, uint8_t *out_b, const void *ctx_p) {
    const pqc__xor_const_ctx_t *c = (const pqc__xor_const_ctx_t *)ctx_p;
    uint32_t nb = (uint32_t)((c->n_bits + 7u) / 8u);
    for (uint32_t i = 0u; i < nb; ++i) {
        out_b[i] = (uint8_t)(in_b[i] ^ (c->target_bytes ? c->target_bytes[i] : 0u));
    }
}

/* Flip phase by pi iff val == target_bytes (using 1 eq_zero ancilla + Z gate) */
static inline void pqc_emit_grover_eq_target(
    qir_module_t *m,
    uint32_t val,
    const uint8_t *target_bytes,
    uint16_t n_bits
) {
    uint32_t diff = val;
    bool all_zero = true;
    uint32_t nb = (uint32_t)((n_bits + 7u) / 8u);
    if (target_bytes) {
        for (uint32_t i = 0u; i < nb; ++i) {
            if (target_bytes[i] != 0u) { all_zero = false; break; }
        }
    }
    if (!all_zero) {
        pqc__xor_const_ctx_t ctx;
        ctx.target_bytes = target_bytes;
        ctx.n_bits       = n_bits;
        diff = pqc_emit_linf2_from_fn(m, &val, 1u, n_bits, n_bits, pqc__fn_xor_const, &ctx, "diff");
    }
    uint32_t is_zero = qir_emit_eq_zero(m, diff, "is_sol");
    qir_phase_term_t term;
    memset(&term, 0, sizeof(term));
    qcc_bv_set(term.mask, 0u);
    term.is_xor = 0u;
    term.angle  = (flow_u128)1u << 127u; /* exact pi phase flip */
    (void)qir_emit_phase_poly(m, &is_zero, 1u, &term, 1u, "grover_pi");
}

/* ============================================================================
 * Family 1: FIPS-197 AES-128 / AES-192 / AES-256 (NIST Categories 1, 3, 5)
 * ============================================================================ */

static inline uint8_t pqc__aes_xtime(uint8_t a) {
    return (uint8_t)(((a << 1u) ^ ((a & 0x80u) ? 0x1Bu : 0u)) & 0xFFu);
}

static inline uint8_t pqc__aes_sbox_affine(uint8_t b) {
    uint8_t res = 0x63u;
    for (uint32_t i = 0u; i < 8u; ++i) {
        uint8_t bit = (uint8_t)(((b >> i) ^ (b >> ((i + 4u) & 7u)) ^
                                 (b >> ((i + 5u) & 7u)) ^ (b >> ((i + 6u) & 7u)) ^
                                 (b >> ((i + 7u) & 7u))) & 1u);
        res ^= (uint8_t)(bit << i);
    }
    return res;
}

static inline uint8_t pqc__aes_sbox(uint8_t x) {
    uint8_t inv = (uint8_t)qir__gf2m_inv_u32(x, 8u, 0x11Bu);
    return pqc__aes_sbox_affine(inv);
}

static inline void pqc_ref_aes_encrypt(
    const uint8_t pt[16],
    const uint8_t *key,
    uint32_t key_bits,
    uint32_t rounds,
    uint8_t out_ct[16]
) {
    static const uint8_t RCON[16] = {
        0x00u, 0x01u, 0x02u, 0x04u, 0x08u, 0x10u, 0x20u, 0x40u,
        0x80u, 0x1Bu, 0x36u, 0x6Cu, 0xD8u, 0xABu, 0x4Du, 0x9Au
    };
    uint32_t nk = key_bits / 32u;
    uint8_t w[240];
    memset(w, 0, sizeof(w));
    memcpy(w, key, nk * 4u);
    uint32_t total_words = (rounds + 1u) * 4u;
    for (uint32_t i = nk; i < total_words; ++i) {
        uint8_t t[4];
        memcpy(t, &w[4u * (i - 1u)], 4u);
        if (i % nk == 0u) {
            uint8_t z0 = pqc__aes_sbox(t[1]) ^ RCON[i / nk];
            t[1] = pqc__aes_sbox(t[2]);
            t[2] = pqc__aes_sbox(t[3]);
            t[3] = pqc__aes_sbox(t[0]);
            t[0] = z0;
        } else if (nk == 8u && (i % nk) == 4u) {
            for (uint32_t b = 0u; b < 4u; ++b) t[b] = pqc__aes_sbox(t[b]);
        }
        for (uint32_t b = 0u; b < 4u; ++b) {
            w[4u * i + b] = (uint8_t)(w[4u * (i - nk) + b] ^ t[b]);
        }
    }

    uint8_t st[16];
    for (uint32_t i = 0u; i < 16u; ++i) st[i] = (uint8_t)(pt[i] ^ w[i]);
    for (uint32_t r = 1u; r <= rounds; ++r) {
        uint8_t sb[16];
        for (uint32_t i = 0u; i < 16u; ++i) sb[i] = pqc__aes_sbox(st[i]);
        uint8_t sr[16];
        for (uint32_t c = 0u; c < 4u; ++c) {
            for (uint32_t row = 0u; row < 4u; ++row) {
                sr[4u * c + row] = sb[4u * ((c + row) & 3u) + row];
            }
        }
        if (r < rounds) {
            for (uint32_t c = 0u; c < 4u; ++c) {
                uint8_t a0 = sr[4u * c], a1 = sr[4u * c + 1u], a2 = sr[4u * c + 2u], a3 = sr[4u * c + 3u];
                uint8_t t = (uint8_t)(a0 ^ a1 ^ a2 ^ a3);
                st[4u * c + 0u] = (uint8_t)(a0 ^ t ^ pqc__aes_xtime((uint8_t)(a0 ^ a1)) ^ w[16u * r + 4u * c + 0u]);
                st[4u * c + 1u] = (uint8_t)(a1 ^ t ^ pqc__aes_xtime((uint8_t)(a1 ^ a2)) ^ w[16u * r + 4u * c + 1u]);
                st[4u * c + 2u] = (uint8_t)(a2 ^ t ^ pqc__aes_xtime((uint8_t)(a2 ^ a3)) ^ w[16u * r + 4u * c + 2u]);
                st[4u * c + 3u] = (uint8_t)(a3 ^ t ^ pqc__aes_xtime((uint8_t)(a3 ^ a0)) ^ w[16u * r + 4u * c + 3u]);
            }
        } else {
            for (uint32_t i = 0u; i < 16u; ++i) st[i] = (uint8_t)(sr[i] ^ w[16u * r + i]);
        }
    }
    memcpy(out_ct, st, 16u);
}

typedef struct {
    uint8_t pt[16];
    uint8_t do_mix;
    uint8_t rcon;
    uint8_t rk0_words; /* For KeySchedule linear layer */
} pqc__aes_fn_ctx_t;

static inline void pqc__aes_fn_init_ark(const uint8_t *in_b, uint8_t *out_b, const void *ctx_p) {
    const pqc__aes_fn_ctx_t *c = (const pqc__aes_fn_ctx_t *)ctx_p;
    for (uint32_t i = 0u; i < 16u; ++i) out_b[i] = (uint8_t)(in_b[i] ^ c->pt[i]);
}

/* Input: in_b[0..15] = inv(state), in_b[16..31] = rk[r]; Output: out_b[0..15] */
static inline void pqc__aes_fn_round_ark(const uint8_t *in_b, uint8_t *out_b, const void *ctx_p) {
    const pqc__aes_fn_ctx_t *c = (const pqc__aes_fn_ctx_t *)ctx_p;
    uint8_t sb[16], sr[16];
    for (uint32_t i = 0u; i < 16u; ++i) sb[i] = pqc__aes_sbox_affine(in_b[i]);
    for (uint32_t col = 0u; col < 4u; ++col) {
        for (uint32_t row = 0u; row < 4u; ++row) {
            sr[4u * col + row] = sb[4u * ((col + row) & 3u) + row];
        }
    }
    if (c->do_mix) {
        for (uint32_t col = 0u; col < 4u; ++col) {
            uint8_t a0 = sr[4u * col], a1 = sr[4u * col + 1u], a2 = sr[4u * col + 2u], a3 = sr[4u * col + 3u];
            uint8_t t = (uint8_t)(a0 ^ a1 ^ a2 ^ a3);
            out_b[4u * col + 0u] = (uint8_t)(a0 ^ t ^ pqc__aes_xtime((uint8_t)(a0 ^ a1)) ^ in_b[16u + 4u * col + 0u]);
            out_b[4u * col + 1u] = (uint8_t)(a1 ^ t ^ pqc__aes_xtime((uint8_t)(a1 ^ a2)) ^ in_b[16u + 4u * col + 1u]);
            out_b[4u * col + 2u] = (uint8_t)(a2 ^ t ^ pqc__aes_xtime((uint8_t)(a2 ^ a3)) ^ in_b[16u + 4u * col + 2u]);
            out_b[4u * col + 3u] = (uint8_t)(a3 ^ t ^ pqc__aes_xtime((uint8_t)(a3 ^ a0)) ^ in_b[16u + 4u * col + 3u]);
        }
    } else {
        for (uint32_t i = 0u; i < 16u; ++i) out_b[i] = (uint8_t)(sr[i] ^ in_b[16u + i]);
    }
}

/* KeySchedule step: in_b[0..4*Nk-1] = prev_block, in_b[4*Nk..4*Nk+3] = inv(t[0..3]) */
static inline void pqc__aes_fn_keysched(const uint8_t *in_b, uint8_t *out_b, const void *ctx_p) {
    const pqc__aes_fn_ctx_t *c = (const pqc__aes_fn_ctx_t *)ctx_p;
    uint32_t nk_bytes = (uint32_t)c->rk0_words * 4u;
    uint8_t sub[4];
    for (uint32_t b = 0u; b < 4u; ++b) sub[b] = pqc__aes_sbox_affine(in_b[nk_bytes + b]);
    sub[0] ^= c->rcon;
    for (uint32_t b = 0u; b < 4u; ++b) out_b[b] = (uint8_t)(in_b[b] ^ sub[b]);
    for (uint32_t w = 1u; w < c->rk0_words; ++w) {
        for (uint32_t b = 0u; b < 4u; ++b) {
            out_b[4u * w + b] = (uint8_t)(in_b[4u * w + b] ^ out_b[4u * (w - 1u) + b]);
        }
    }
}

static inline bool pqc_build_aes_qir(
    qir_module_t *m,
    uint32_t key_bits,
    uint32_t rounds,
    const uint8_t pt[16],
    const uint8_t target_ct[16],
    bool phase_oracle
) {
    if (key_bits != 128u && key_bits != 192u && key_bits != 256u) return false;
    if (rounds == 0u || rounds > 14u) return false;
    static const uint8_t RCON[16] = {
        0x00u, 0x01u, 0x02u, 0x04u, 0x08u, 0x10u, 0x20u, 0x40u,
        0x80u, 0x1Bu, 0x36u, 0x6Cu, 0xD8u, 0xABu, 0x4Du, 0x9Au
    };

    uint32_t key_in = qir_add_input(m, (uint16_t)key_bits, "key");
    if (phase_oracle) qir_add_output(m, key_in);

    uint32_t nk = key_bits / 32u;
    uint32_t total_words = (rounds + 1u) * 4u;
    uint32_t word_vals[64];
    memset(word_vals, 0, sizeof(word_vals));

    for (uint32_t i = 0u; i < nk; ++i) {
        word_vals[i] = qir_emit_slice(m, key_in, (uint16_t)(i * 32u), (uint16_t)((i + 1u) * 32u), "kw");
    }

    uint32_t w_pos = nk;
    uint32_t rc_idx = 1u;
    while (w_pos < total_words) {
        if (nk == 4u || nk == 6u) {
            uint32_t prev_blk = qir_emit_concat(m, &word_vals[w_pos - nk], (uint16_t)nk, "kblk");
            uint32_t last_w   = word_vals[w_pos - 1u];
            uint32_t sub_in[4];
            sub_in[0] = qir_emit_slice(m, last_w, 8u,  16u, "ks0");
            sub_in[1] = qir_emit_slice(m, last_w, 16u, 24u, "ks1");
            sub_in[2] = qir_emit_slice(m, last_w, 24u, 32u, "ks2");
            sub_in[3] = qir_emit_slice(m, last_w, 0u,  8u,  "ks3");
            uint32_t sub_inv[4];
            for (uint32_t b = 0u; b < 4u; ++b) {
                sub_inv[b] = qir_emit_gf2m_inv(m, sub_in[b], 8u, 0x11Bu, "kinv");
            }
            uint32_t sub_cat = qir_emit_concat(m, sub_inv, 4u, "ksub");
            uint32_t lins[2] = {prev_blk, sub_cat};
            pqc__aes_fn_ctx_t kctx;
            memset(&kctx, 0, sizeof(kctx));
            kctx.rcon      = RCON[rc_idx++];
            kctx.rk0_words = (uint8_t)nk;
            uint32_t next_blk = pqc_emit_linf2_from_fn(m, lins, 2u, (uint16_t)(nk * 32u + 32u),
                                                       (uint16_t)(nk * 32u), pqc__aes_fn_keysched, &kctx, "knext");
            for (uint32_t j = 0u; j < nk && w_pos < total_words; ++j, ++w_pos) {
                word_vals[w_pos] = qir_emit_slice(m, next_blk, (uint16_t)(j * 32u), (uint16_t)((j + 1u) * 32u), "kw");
            }
        } else {
            /* AES-256: generate 4-word half-blocks (128b + 32b = 160b <= 256b) */
            bool is_rot = (((w_pos - nk) & 7u) == 0u);
            uint32_t base_blk = qir_emit_concat(m, &word_vals[w_pos - 8u], 4u, "khalf");
            uint32_t last_w   = word_vals[w_pos - 1u];
            uint32_t sub_in[4];
            if (is_rot) {
                sub_in[0] = qir_emit_slice(m, last_w, 8u,  16u, "ks0");
                sub_in[1] = qir_emit_slice(m, last_w, 16u, 24u, "ks1");
                sub_in[2] = qir_emit_slice(m, last_w, 24u, 32u, "ks2");
                sub_in[3] = qir_emit_slice(m, last_w, 0u,  8u,  "ks3");
            } else {
                sub_in[0] = qir_emit_slice(m, last_w, 0u,  8u,  "ks0");
                sub_in[1] = qir_emit_slice(m, last_w, 8u,  16u, "ks1");
                sub_in[2] = qir_emit_slice(m, last_w, 16u, 24u, "ks2");
                sub_in[3] = qir_emit_slice(m, last_w, 24u, 32u, "ks3");
            }
            uint32_t sub_inv[4];
            for (uint32_t b = 0u; b < 4u; ++b) {
                sub_inv[b] = qir_emit_gf2m_inv(m, sub_in[b], 8u, 0x11Bu, "kinv");
            }
            uint32_t sub_cat = qir_emit_concat(m, sub_inv, 4u, "ksub");
            uint32_t lins[2] = {base_blk, sub_cat};
            pqc__aes_fn_ctx_t kctx;
            memset(&kctx, 0, sizeof(kctx));
            kctx.rcon      = is_rot ? RCON[rc_idx++] : 0u;
            kctx.rk0_words = 4u;
            uint32_t next_blk = pqc_emit_linf2_from_fn(m, lins, 2u, 160u, 128u, pqc__aes_fn_keysched, &kctx, "knext");
            for (uint32_t j = 0u; j < 4u && w_pos < total_words; ++j, ++w_pos) {
                word_vals[w_pos] = qir_emit_slice(m, next_blk, (uint16_t)(j * 32u), (uint16_t)((j + 1u) * 32u), "kw");
            }
        }
    }

    pqc__aes_fn_ctx_t rctx;
    memset(&rctx, 0, sizeof(rctx));
    memcpy(rctx.pt, pt, 16u);
    uint32_t rk0 = qir_emit_concat(m, &word_vals[0], 4u, "rk0");
    uint32_t state = pqc_emit_linf2_from_fn(m, &rk0, 1u, 128u, 128u, pqc__aes_fn_init_ark, &rctx, "st0");

    for (uint32_t r = 1u; r <= rounds; ++r) {
        uint32_t inv_b[16];
        for (uint32_t b = 0u; b < 16u; ++b) {
            uint32_t sl = qir_emit_slice(m, state, (uint16_t)(b * 8u), (uint16_t)((b + 1u) * 8u), "sb");
            inv_b[b] = qir_emit_gf2m_inv(m, sl, 8u, 0x11Bu, "inv");
        }
        uint32_t inv_cat = qir_emit_concat(m, inv_b, 16u, "inv128");
        uint32_t rkr = qir_emit_concat(m, &word_vals[4u * r], 4u, "rkr");
        uint32_t rins[2] = {inv_cat, rkr};
        rctx.do_mix = (uint8_t)((r < rounds) ? 1u : 0u);
        state = pqc_emit_linf2_from_fn(m, rins, 2u, 256u, 128u, pqc__aes_fn_round_ark, &rctx,
                                       (!phase_oracle && r == rounds) ? "ct" : "st");
    }

    if (phase_oracle) pqc_emit_grover_eq_target(m, state, target_ct, 128u);
    else              qir_add_output(m, state);
    return true;
}

/* ============================================================================
 * Family 2: Hash & Sponge (SHA-256, Ascon-p[12], Keccak-f[400])
 * ============================================================================ */

static const uint32_t PQC_SHA256_IV[8] = {
    0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
    0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u
};

static const uint32_t PQC_SHA256_K[64] = {
    0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
    0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u, 0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
    0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
    0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
    0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u, 0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
    0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
    0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
    0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u, 0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u
};

static inline void pqc_ref_sha256_bytes(const uint8_t *msg, size_t len, uint8_t digest[32]) {
    uint32_t h[8];
    memcpy(h, PQC_SHA256_IV, sizeof(h));
    size_t padded_len = ((len + 9u + 63u) / 64u) * 64u;
    uint8_t *buf = (uint8_t *)calloc(padded_len, 1u);
    if (len > 0u) memcpy(buf, msg, len);
    buf[len] = 0x80u;
    uint64_t bit_len = (uint64_t)len * 8u;
    for (uint32_t b = 0u; b < 8u; ++b) buf[padded_len - 1u - b] = (uint8_t)(bit_len >> (8u * b));

    for (size_t off = 0u; off < padded_len; off += 64u) {
        uint32_t w[64];
        for (uint32_t i = 0u; i < 16u; ++i) w[i] = pqc__load_be32(buf + off + 4u * i);
        for (uint32_t i = 16u; i < 64u; ++i) {
            uint32_t s0 = pqc__rotr32(w[i - 15u], 7u) ^ pqc__rotr32(w[i - 15u], 18u) ^ (w[i - 15u] >> 3u);
            uint32_t s1 = pqc__rotr32(w[i - 2u], 17u) ^ pqc__rotr32(w[i - 2u], 19u) ^ (w[i - 2u] >> 10u);
            w[i] = w[i - 16u] + s0 + w[i - 7u] + s1;
        }
        uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
        for (uint32_t i = 0u; i < 64u; ++i) {
            uint32_t s1 = pqc__rotr32(e, 6u) ^ pqc__rotr32(e, 11u) ^ pqc__rotr32(e, 25u);
            uint32_t ch = (e & f) ^ ((~e) & g);
            uint32_t t1 = hh + s1 + ch + PQC_SHA256_K[i] + w[i];
            uint32_t s0 = pqc__rotr32(a, 2u) ^ pqc__rotr32(a, 13u) ^ pqc__rotr32(a, 22u);
            uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
            uint32_t t2 = s0 + maj;
            hh = g; g = f; f = e; e = d + t1;
            d = c; c = b; b = a; a = t1 + t2;
        }
        h[0] += a; h[1] += b; h[2] += c; h[3] += d;
        h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
    }
    free(buf);
    for (uint32_t i = 0u; i < 8u; ++i) pqc__store_be32(digest + 4u * i, h[i]);
}

/* Fixed-size 256-bit preimage block matching quantum circuit bit layout (8 little-endian 32-bit words) */
static inline void pqc_ref_sha256_msg256(const uint8_t msg32[32], uint32_t rounds, uint8_t digest32[32]) {
    uint32_t w[64];
    memset(w, 0, sizeof(w));
    for (uint32_t i = 0u; i < 8u; ++i) w[i] = pqc__load_le32(msg32 + 4u * i);
    w[8] = 0x80000000u;
    w[15] = 256u;
    for (uint32_t i = 16u; i < rounds; ++i) {
        uint32_t s0 = pqc__rotr32(w[i - 15u], 7u) ^ pqc__rotr32(w[i - 15u], 18u) ^ (w[i - 15u] >> 3u);
        uint32_t s1 = pqc__rotr32(w[i - 2u], 17u) ^ pqc__rotr32(w[i - 2u], 19u) ^ (w[i - 2u] >> 10u);
        w[i] = w[i - 16u] + s0 + w[i - 7u] + s1;
    }
    uint32_t h[8];
    memcpy(h, PQC_SHA256_IV, sizeof(h));
    uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
    for (uint32_t i = 0u; i < rounds; ++i) {
        uint32_t s1 = pqc__rotr32(e, 6u) ^ pqc__rotr32(e, 11u) ^ pqc__rotr32(e, 25u);
        uint32_t ch = (e & f) ^ ((~e) & g);
        uint32_t t1 = hh + s1 + ch + PQC_SHA256_K[i] + w[i];
        uint32_t s0 = pqc__rotr32(a, 2u) ^ pqc__rotr32(a, 13u) ^ pqc__rotr32(a, 22u);
        uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
        uint32_t t2 = s0 + maj;
        hh = g; g = f; f = e; e = d + t1;
        d = c; c = b; b = a; a = t1 + t2;
    }
    h[0] += a; h[1] += b; h[2] += c; h[3] += d;
    h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
    for (uint32_t i = 0u; i < 8u; ++i) pqc__store_le32(digest32 + 4u * i, h[i]);
}

static inline void pqc__sha_fn_s0(const uint8_t *in_b, uint8_t *out_b, const void *ctx_p) {
    (void)ctx_p;
    uint32_t x = pqc__pack_u32_le(in_b);
    uint32_t y = pqc__rotr32(x, 7u) ^ pqc__rotr32(x, 18u) ^ (x >> 3u);
    pqc__unpack_u32_le(y, out_b);
}
static inline void pqc__sha_fn_s1(const uint8_t *in_b, uint8_t *out_b, const void *ctx_p) {
    (void)ctx_p;
    uint32_t x = pqc__pack_u32_le(in_b);
    uint32_t y = pqc__rotr32(x, 17u) ^ pqc__rotr32(x, 19u) ^ (x >> 10u);
    pqc__unpack_u32_le(y, out_b);
}
static inline void pqc__sha_fn_S0(const uint8_t *in_b, uint8_t *out_b, const void *ctx_p) {
    (void)ctx_p;
    uint32_t x = pqc__pack_u32_le(in_b);
    uint32_t y = pqc__rotr32(x, 2u) ^ pqc__rotr32(x, 13u) ^ pqc__rotr32(x, 22u);
    pqc__unpack_u32_le(y, out_b);
}
static inline void pqc__sha_fn_S1(const uint8_t *in_b, uint8_t *out_b, const void *ctx_p) {
    (void)ctx_p;
    uint32_t x = pqc__pack_u32_le(in_b);
    uint32_t y = pqc__rotr32(x, 6u) ^ pqc__rotr32(x, 11u) ^ pqc__rotr32(x, 25u);
    pqc__unpack_u32_le(y, out_b);
}

static inline bool pqc_build_sha256_qir(
    qir_module_t *m,
    uint32_t rounds,
    const uint8_t target_digest[32],
    bool phase_oracle
) {
    if (rounds == 0u || rounds > 64u) return false;
    uint32_t preimage = qir_add_input(m, 256u, "preimage");
    if (phase_oracle) qir_add_output(m, preimage);

    uint32_t w[64];
    bool     w_is_const[64];
    uint32_t w_const_val[64];
    memset(w_is_const, 0, sizeof(w_is_const));
    memset(w_const_val, 0, sizeof(w_const_val));

    for (uint32_t i = 0u; i < 8u; ++i) {
        w[i] = qir_emit_slice(m, preimage, (uint16_t)(i * 32u), (uint16_t)((i + 1u) * 32u), "w");
    }
    w_is_const[8] = true; w_const_val[8] = 0x80000000u;
    for (uint32_t i = 9u; i < 15u; ++i) { w_is_const[i] = true; w_const_val[i] = 0u; }
    w_is_const[15] = true; w_const_val[15] = 256u;

    for (uint32_t i = 16u; i < rounds; ++i) {
        uint32_t s0 = w_is_const[i - 15u] ? UINT32_MAX : pqc_emit_linf2_from_fn(m, &w[i - 15u], 1u, 32u, 32u, pqc__sha_fn_s0, NULL, "s0");
        uint32_t s1 = w_is_const[i - 2u]  ? UINT32_MAX : pqc_emit_linf2_from_fn(m, &w[i - 2u],  1u, 32u, 32u, pqc__sha_fn_s1, NULL, "s1");
        uint32_t c_acc = 0u;
        uint32_t qterms[4], nq = 0u;
        if (w_is_const[i - 16u]) {
            c_acc += w_const_val[i - 16u];
        } else {
            qterms[nq++] = w[i - 16u];
        }
        if (w_is_const[i - 15u]) {
            uint32_t xv = w_const_val[i - 15u];
            c_acc += (pqc__rotr32(xv, 7u) ^ pqc__rotr32(xv, 18u) ^ (xv >> 3u));
        } else {
            qterms[nq++] = s0;
        }
        if (w_is_const[i - 7u]) {
            c_acc += w_const_val[i - 7u];
        } else {
            qterms[nq++] = w[i - 7u];
        }
        if (w_is_const[i - 2u]) {
            uint32_t xv = w_const_val[i - 2u];
            c_acc += (pqc__rotr32(xv, 17u) ^ pqc__rotr32(xv, 19u) ^ (xv >> 10u));
        } else {
            qterms[nq++] = s1;
        }
        if (nq == 0u) {
            w_is_const[i]  = true;
            w_const_val[i] = c_acc;
        } else {
            uint32_t acc = qterms[0];
            for (uint32_t t = 1u; t < nq; ++t) acc = qir_emit_add(m, acc, qterms[t], 0u, "wsum");
            if (c_acc != 0u) {
                uint32_t kc = qir_emit_const_u64(m, 32u, c_acc, "wc");
                acc = qir_emit_add(m, acc, kc, 0u, "wsumc");
            }
            w[i] = acc;
        }
    }

    uint32_t st[8];
    uint32_t st_c[8];
    bool     st_is_c[8];
    for (uint32_t i = 0u; i < 8u; ++i) {
        st_c[i]    = PQC_SHA256_IV[i];
        st_is_c[i] = true;
        st[i]      = UINT32_MAX;
    }

    for (uint32_t i = 0u; i < rounds; ++i) {
        /* Compute t1 = hh + S1(e) + Ch(e,f,g) + K[i] + w[i] */
        uint32_t c_sum = PQC_SHA256_K[i];
        uint32_t terms1[4], n1 = 0u;
        if (st_is_c[7]) c_sum += st_c[7]; else terms1[n1++] = st[7];
        if (w_is_const[i]) c_sum += w_const_val[i]; else terms1[n1++] = w[i];

        if (st_is_c[4]) {
            uint32_t ev = st_c[4];
            c_sum += (pqc__rotr32(ev, 6u) ^ pqc__rotr32(ev, 11u) ^ pqc__rotr32(ev, 25u));
        } else {
            terms1[n1++] = pqc_emit_linf2_from_fn(m, &st[4], 1u, 32u, 32u, pqc__sha_fn_S1, NULL, "S1");
        }

        /* Ch(e,f,g) = g ^ (e & (f ^ g)) */
        if (st_is_c[4] && st_is_c[5] && st_is_c[6]) {
            c_sum += ((st_c[4] & st_c[5]) ^ ((~st_c[4]) & st_c[6]));
        } else {
            uint32_t ev = st_is_c[4] ? qir_emit_const_u64(m, 32u, st_c[4], "ec") : st[4];
            uint32_t fv = st_is_c[5] ? qir_emit_const_u64(m, 32u, st_c[5], "fc") : st[5];
            uint32_t gv = st_is_c[6] ? qir_emit_const_u64(m, 32u, st_c[6], "gc") : st[6];
            uint32_t fg = qir_emit_gf2m_add(m, fv, gv, 32u, "fg");
            uint32_t a1 = qir_emit_and_vec(m, ev, fg, 0u, "ch_a");
            uint32_t ch = qir_emit_gf2m_add(m, gv, a1, 32u, "ch");
            terms1[n1++] = ch;
        }

        uint32_t t1 = UINT32_MAX;
        bool t1_is_c = (n1 == 0u);
        uint32_t t1_c = c_sum;
        if (!t1_is_c) {
            t1 = terms1[0];
            for (uint32_t t = 1u; t < n1; ++t) t1 = qir_emit_add(m, t1, terms1[t], 0u, "t1");
            if (c_sum != 0u) {
                uint32_t kc = qir_emit_const_u64(m, 32u, c_sum, "c1");
                t1 = qir_emit_add(m, t1, kc, 0u, "t1c");
            }
        }

        /* Compute t2 = S0(a) + Maj(a,b,c), Maj(a,b,c) = a ^ ((a^b)&(a^c)) */
        uint32_t c_sum2 = 0u;
        uint32_t terms2[2], n2 = 0u;
        if (st_is_c[0]) {
            uint32_t av = st_c[0];
            c_sum2 += (pqc__rotr32(av, 2u) ^ pqc__rotr32(av, 13u) ^ pqc__rotr32(av, 22u));
        } else {
            terms2[n2++] = pqc_emit_linf2_from_fn(m, &st[0], 1u, 32u, 32u, pqc__sha_fn_S0, NULL, "S0");
        }
        if (st_is_c[0] && st_is_c[1] && st_is_c[2]) {
            c_sum2 += ((st_c[0] & st_c[1]) ^ (st_c[0] & st_c[2]) ^ (st_c[1] & st_c[2]));
        } else {
            uint32_t av = st_is_c[0] ? qir_emit_const_u64(m, 32u, st_c[0], "ac") : st[0];
            uint32_t bv = st_is_c[1] ? qir_emit_const_u64(m, 32u, st_c[1], "bc") : st[1];
            uint32_t cv = st_is_c[2] ? qir_emit_const_u64(m, 32u, st_c[2], "cc") : st[2];
            uint32_t ab = qir_emit_gf2m_add(m, av, bv, 32u, "ab");
            uint32_t ac = qir_emit_gf2m_add(m, av, cv, 32u, "ac");
            uint32_t m1 = qir_emit_and_vec(m, ab, ac, 0u, "mj_a");
            terms2[n2++] = qir_emit_gf2m_add(m, av, m1, 32u, "maj");
        }

        uint32_t new_e = UINT32_MAX, new_a = UINT32_MAX;
        bool new_e_c = false, new_a_c = false;
        uint32_t new_e_cv = 0u, new_a_cv = 0u;

        if (t1_is_c && st_is_c[3]) {
            new_e_c = true; new_e_cv = st_c[3] + t1_c;
        } else if (t1_is_c) {
            uint32_t tc = qir_emit_const_u64(m, 32u, t1_c, "t1k");
            new_e = qir_emit_add(m, st[3], tc, 0u, "ne");
        } else if (st_is_c[3]) {
            uint32_t dc = qir_emit_const_u64(m, 32u, st_c[3], "dc");
            new_e = qir_emit_add(m, t1, dc, 0u, "ne");
        } else {
            new_e = qir_emit_add(m, st[3], t1, 0u, "ne");
        }

        if (t1_is_c && n2 == 0u) {
            new_a_c = true; new_a_cv = t1_c + c_sum2;
        } else {
            uint32_t at[3], na = 0u;
            if (!t1_is_c) at[na++] = t1;
            else c_sum2 += t1_c;
            for (uint32_t t = 0u; t < n2; ++t) at[na++] = terms2[t];
            uint32_t aa = at[0];
            for (uint32_t t = 1u; t < na; ++t) aa = qir_emit_add(m, aa, at[t], 0u, "na");
            if (c_sum2 != 0u) {
                uint32_t kc = qir_emit_const_u64(m, 32u, c_sum2, "c2");
                aa = qir_emit_add(m, aa, kc, 0u, "nac");
            }
            new_a = aa;
        }

        st[7] = st[6]; st_c[7] = st_c[6]; st_is_c[7] = st_is_c[6];
        st[6] = st[5]; st_c[6] = st_c[5]; st_is_c[6] = st_is_c[5];
        st[5] = st[4]; st_c[5] = st_c[4]; st_is_c[5] = st_is_c[4];
        st[4] = new_e; st_c[4] = new_e_cv; st_is_c[4] = new_e_c;
        st[3] = st[2]; st_c[3] = st_c[2]; st_is_c[3] = st_is_c[2];
        st[2] = st[1]; st_c[2] = st_c[1]; st_is_c[2] = st_is_c[1];
        st[1] = st[0]; st_c[1] = st_c[0]; st_is_c[1] = st_is_c[0];
        st[0] = new_a; st_c[0] = new_a_cv; st_is_c[0] = new_a_c;
    }

    uint32_t out_w[8];
    for (uint32_t i = 0u; i < 8u; ++i) {
        if (st_is_c[i]) {
            out_w[i] = qir_emit_const_u64(m, 32u, PQC_SHA256_IV[i] + st_c[i], "hc");
        } else {
            uint32_t iv = qir_emit_const_u64(m, 32u, PQC_SHA256_IV[i], "iv");
            out_w[i] = qir_emit_add(m, iv, st[i], 0u, "hout");
        }
    }
    uint32_t digest = qir_emit_concat(m, out_w, 8u, "digest");
    if (phase_oracle) pqc_emit_grover_eq_target(m, digest, target_digest, 256u);
    else              qir_add_output(m, digest);
    return true;
}

/* --- NIST SP 800-232 Ascon-p[12] --- */
static inline void pqc_ref_ascon_p(uint64_t s[5], uint32_t rounds) {
    static const uint8_t RC[12] = {
        0xF0u, 0xE1u, 0xD2u, 0xC3u, 0xB4u, 0xA5u,
        0x96u, 0x87u, 0x78u, 0x69u, 0x5Au, 0x4Bu
    };
    uint32_t start = 12u - rounds;
    for (uint32_t r = start; r < 12u; ++r) {
        s[2] ^= (uint64_t)RC[r];
        s[0] ^= s[4]; s[4] ^= s[3]; s[2] ^= s[1];
        uint64_t t0 = (~s[0]) & s[1];
        uint64_t t1 = (~s[1]) & s[2];
        uint64_t t2 = (~s[2]) & s[3];
        uint64_t t3 = (~s[3]) & s[4];
        uint64_t t4 = (~s[4]) & s[0];
        s[0] ^= t1; s[1] ^= t2; s[2] ^= t3; s[3] ^= t4; s[4] ^= t0;
        s[1] ^= s[0]; s[0] ^= s[4]; s[3] ^= s[2]; s[2] = ~s[2];
        s[0] ^= pqc__rotr64(s[0], 19u) ^ pqc__rotr64(s[0], 28u);
        s[1] ^= pqc__rotr64(s[1], 61u) ^ pqc__rotr64(s[1], 39u);
        s[2] ^= pqc__rotr64(s[2], 1u)  ^ pqc__rotr64(s[2], 6u);
        s[3] ^= pqc__rotr64(s[3], 10u) ^ pqc__rotr64(s[3], 17u);
        s[4] ^= pqc__rotr64(s[4], 7u)  ^ pqc__rotr64(s[4], 41u);
    }
}

typedef struct {
    uint32_t lane_idx;
    uint64_t not_mask;
} pqc__ascon_diff_ctx_t;

static inline void pqc__ascon_fn_sigma(const uint8_t *in_b, uint8_t *out_b, const void *ctx_p) {
    static const uint8_t R0[5] = {19u, 61u, 1u, 10u, 7u};
    static const uint8_t R1[5] = {28u, 39u, 6u, 17u, 41u};
    const pqc__ascon_diff_ctx_t *c = (const pqc__ascon_diff_ctx_t *)ctx_p;
    uint64_t x = pqc__pack_u64_le(in_b) ^ c->not_mask;
    uint64_t y = x ^ pqc__rotr64(x, R0[c->lane_idx]) ^ pqc__rotr64(x, R1[c->lane_idx]);
    pqc__unpack_u64_le(y, out_b);
}

static inline bool pqc_build_ascon_qir(
    qir_module_t *m,
    uint32_t rounds,
    const uint64_t nonce_iv[3],
    const uint64_t target_tag[2],
    bool phase_oracle
) {
    if (rounds == 0u || rounds > 12u) return false;
    static const uint8_t RC[12] = {
        0xF0u, 0xE1u, 0xD2u, 0xC3u, 0xB4u, 0xA5u,
        0x96u, 0x87u, 0x78u, 0x69u, 0x5Au, 0x4Bu
    };
    uint32_t key_in = qir_add_input(m, 128u, "key");
    if (phase_oracle) qir_add_output(m, key_in);

    uint32_t x[5];
    x[0] = qir_emit_const_u64(m, 64u, nonce_iv[0], "iv");
    x[1] = qir_emit_slice(m, key_in, 0u, 64u, "k0");
    x[2] = qir_emit_slice(m, key_in, 64u, 128u, "k1");
    x[3] = qir_emit_const_u64(m, 64u, nonce_iv[1], "n0");
    x[4] = qir_emit_const_u64(m, 64u, nonce_iv[2], "n1");

    uint32_t start = 12u - rounds;
    for (uint32_t r = start; r < 12u; ++r) {
        uint32_t rc = qir_emit_const_u64(m, 64u, RC[r], "rc");
        x[2] = qir_emit_gf2m_add(m, x[2], rc, 64u, "x2c");
        x[0] = qir_emit_gf2m_add(m, x[0], x[4], 64u, "p0");
        x[4] = qir_emit_gf2m_add(m, x[4], x[3], 64u, "p4");
        x[2] = qir_emit_gf2m_add(m, x[2], x[1], 64u, "p2");

        /* (~x_i) & x_{i+1} = x_{i+1} ^ (x_i & x_{i+1}) */
        uint32_t ands[5], ts[5];
        for (uint32_t i = 0u; i < 5u; ++i) {
            uint32_t nxt = (i + 1u < 5u) ? (i + 1u) : 0u;
            ands[i] = qir_emit_and_vec(m, x[i], x[nxt], 0u, "aa");
            ts[i]   = qir_emit_gf2m_add(m, x[nxt], ands[i], 64u, "at");
        }
        uint32_t y[5];
        y[0] = qir_emit_gf2m_add(m, x[0], ts[1], 64u, "y0");
        y[1] = qir_emit_gf2m_add(m, x[1], ts[2], 64u, "y1");
        y[2] = qir_emit_gf2m_add(m, x[2], ts[3], 64u, "y2");
        y[3] = qir_emit_gf2m_add(m, x[3], ts[4], 64u, "y3");
        y[4] = qir_emit_gf2m_add(m, x[4], ts[0], 64u, "y4");

        y[1] = qir_emit_gf2m_add(m, y[1], y[0], 64u, "z1");
        y[0] = qir_emit_gf2m_add(m, y[0], y[4], 64u, "z0");
        y[3] = qir_emit_gf2m_add(m, y[3], y[2], 64u, "z3");

        for (uint32_t i = 0u; i < 5u; ++i) {
            pqc__ascon_diff_ctx_t dctx;
            dctx.lane_idx = i;
            dctx.not_mask = (i == 2u) ? UINT64_MAX : 0u;
            x[i] = pqc_emit_linf2_from_fn(m, &y[i], 1u, 64u, 64u, pqc__ascon_fn_sigma, &dctx, "ax");
        }
    }

    uint32_t out_lanes[2] = {x[0], x[1]};
    uint32_t out128 = qir_emit_concat(m, out_lanes, 2u, "ascon_out");
    if (phase_oracle) {
        uint8_t tgt[16];
        pqc__unpack_u64_le(target_tag[0], tgt);
        pqc__unpack_u64_le(target_tag[1], tgt + 8);
        pqc_emit_grover_eq_target(m, out128, tgt, 128u);
    } else {
        qir_add_output(m, out128);
    }
    return true;
}

/* --- SHA-3 Keccak-f[400] (25 x 16-bit lanes, 20 rounds) --- */
static const uint16_t PQC_KECCAK_F400_RC[20] = {
    0x0001u, 0x8082u, 0x808Au, 0x8000u, 0x808Bu, 0x0001u, 0x8081u, 0x8009u,
    0x008Au, 0x0088u, 0x8009u, 0x000Au, 0x808Bu, 0x008Bu, 0x8089u, 0x8003u,
    0x8002u, 0x0080u, 0x800Au, 0x000Au
};

static const uint8_t PQC_KECCAK_ROT[25] = {
     0u,  1u, 62u, 28u, 27u,
    36u, 44u,  6u, 55u, 20u,
     3u, 10u, 43u, 25u, 39u,
    41u, 45u, 15u, 21u,  8u,
    18u,  2u, 61u, 56u, 14u
};

static inline void pqc_ref_keccak_f400(uint16_t st[25], uint32_t rounds) {
    for (uint32_t r = 0u; r < rounds; ++r) {
        uint16_t c[5], d[5], b[25];
        for (uint32_t x = 0u; x < 5u; ++x) {
            c[x] = (uint16_t)(st[x] ^ st[x + 5u] ^ st[x + 10u] ^ st[x + 15u] ^ st[x + 20u]);
        }
        for (uint32_t x = 0u; x < 5u; ++x) {
            d[x] = (uint16_t)(c[(x + 4u) % 5u] ^ pqc__rotl16(c[(x + 1u) % 5u], 1u));
        }
        for (uint32_t y = 0u; y < 5u; ++y) {
            for (uint32_t x = 0u; x < 5u; ++x) {
                uint32_t idx = y * 5u + x;
                uint32_t new_x = y;
                uint32_t new_y = (2u * x + 3u * y) % 5u;
                b[new_y * 5u + new_x] = pqc__rotl16((uint16_t)(st[idx] ^ d[x]), PQC_KECCAK_ROT[idx] & 15u);
            }
        }
        for (uint32_t y = 0u; y < 5u; ++y) {
            for (uint32_t x = 0u; x < 5u; ++x) {
                uint16_t v0 = b[y * 5u + x];
                uint16_t v1 = b[y * 5u + (x + 1u) % 5u];
                uint16_t v2 = b[y * 5u + (x + 2u) % 5u];
                st[y * 5u + x] = (uint16_t)(v0 ^ ((uint16_t)(~v1) & v2));
            }
        }
        st[0] ^= PQC_KECCAK_F400_RC[r % 20u];
    }
}

static inline void pqc__keccak_fn_parity16(const uint8_t *in_b, uint8_t *out_b, const void *ctx_p) {
    (void)ctx_p;
    uint16_t acc = 0u;
    for (uint32_t k = 0u; k < 5u; ++k) {
        uint16_t v = (uint16_t)in_b[2u * k] | ((uint16_t)in_b[2u * k + 1u] << 8u);
        acc ^= v;
    }
    out_b[0] = (uint8_t)(acc & 0xFFu);
    out_b[1] = (uint8_t)(acc >> 8u);
}

static inline void pqc__keccak_fn_drow(const uint8_t *in_b, uint8_t *out_b, const void *ctx_p) {
    (void)ctx_p;
    uint16_t c[5];
    for (uint32_t x = 0u; x < 5u; ++x) c[x] = (uint16_t)in_b[2u * x] | ((uint16_t)in_b[2u * x + 1u] << 8u);
    for (uint32_t x = 0u; x < 5u; ++x) {
        uint16_t dx = (uint16_t)(c[(x + 4u) % 5u] ^ pqc__rotl16(c[(x + 1u) % 5u], 1u));
        out_b[2u * x]      = (uint8_t)(dx & 0xFFu);
        out_b[2u * x + 1u] = (uint8_t)(dx >> 8u);
    }
}

static inline void pqc__keccak_fn_theta_rot(const uint8_t *in_b, uint8_t *out_b, const void *ctx_p) {
    uint32_t sh = *(const uint32_t *)ctx_p;
    uint16_t a = (uint16_t)in_b[0] | ((uint16_t)in_b[1] << 8u);
    uint16_t d = (uint16_t)in_b[2] | ((uint16_t)in_b[3] << 8u);
    uint16_t y = pqc__rotl16((uint16_t)(a ^ d), sh);
    out_b[0] = (uint8_t)(y & 0xFFu);
    out_b[1] = (uint8_t)(y >> 8u);
}

static inline bool pqc_build_keccak_f400_qir(
    qir_module_t *m,
    uint32_t rounds,
    const uint16_t target_rate[8],
    bool phase_oracle
) {
    if (rounds == 0u || rounds > 20u) return false;
    uint32_t seed_in = qir_add_input(m, 128u, "seed");
    if (phase_oracle) qir_add_output(m, seed_in);

    uint32_t st[25];
    for (uint32_t i = 0u; i < 8u; ++i) {
        st[i] = qir_emit_slice(m, seed_in, (uint16_t)(i * 16u), (uint16_t)((i + 1u) * 16u), "ks");
    }
    for (uint32_t i = 8u; i < 25u; ++i) {
        st[i] = qir_emit_const_u64(m, 16u, 0u, "kz");
    }

    for (uint32_t r = 0u; r < rounds; ++r) {
        uint32_t c_lanes[5];
        for (uint32_t x = 0u; x < 5u; ++x) {
            uint32_t col[5] = {st[x], st[x + 5u], st[x + 10u], st[x + 15u], st[x + 20u]};
            c_lanes[x] = pqc_emit_linf2_from_fn(m, col, 5u, 80u, 16u, pqc__keccak_fn_parity16, NULL, "kc");
        }
        uint32_t d_all = pqc_emit_linf2_from_fn(m, c_lanes, 5u, 80u, 80u, pqc__keccak_fn_drow, NULL, "kd");
        uint32_t d[5];
        for (uint32_t x = 0u; x < 5u; ++x) {
            d[x] = qir_emit_slice(m, d_all, (uint16_t)(x * 16u), (uint16_t)((x + 1u) * 16u), "kdx");
        }
        uint32_t b[25];
        for (uint32_t y = 0u; y < 5u; ++y) {
            for (uint32_t x = 0u; x < 5u; ++x) {
                uint32_t idx   = y * 5u + x;
                uint32_t new_x = y;
                uint32_t new_y = (2u * x + 3u * y) % 5u;
                uint32_t sh    = PQC_KECCAK_ROT[idx] & 15u;
                uint32_t ins[2] = {st[idx], d[x]};
                b[new_y * 5u + new_x] = pqc_emit_linf2_from_fn(m, ins, 2u, 32u, 16u, pqc__keccak_fn_theta_rot, &sh, "kb");
            }
        }
        for (uint32_t y = 0u; y < 5u; ++y) {
            for (uint32_t x = 0u; x < 5u; ++x) {
                uint32_t b0 = b[y * 5u + x];
                uint32_t b1 = b[y * 5u + (x + 1u) % 5u];
                uint32_t b2 = b[y * 5u + (x + 2u) % 5u];
                uint32_t a12 = qir_emit_and_vec(m, b1, b2, 0u, "ka");
                uint32_t t1  = qir_emit_gf2m_add(m, b2, a12, 16u, "kt");
                st[y * 5u + x] = qir_emit_gf2m_add(m, b0, t1, 16u, "kx");
            }
        }
        uint32_t rc = qir_emit_const_u64(m, 16u, PQC_KECCAK_F400_RC[r % 20u], "krc");
        st[0] = qir_emit_gf2m_add(m, st[0], rc, 16u, "kiota");
    }

    uint32_t out128 = qir_emit_concat(m, st, 8u, "keccak_rate");
    if (phase_oracle) {
        uint8_t tgt[16];
        for (uint32_t i = 0u; i < 8u; ++i) {
            tgt[2u * i]      = (uint8_t)(target_rate[i] & 0xFFu);
            tgt[2u * i + 1u] = (uint8_t)(target_rate[i] >> 8u);
        }
        pqc_emit_grover_eq_target(m, out128, tgt, 128u);
    } else {
        qir_add_output(m, out128);
    }
    return true;
}

/* ============================================================================
 * Family 3: Code-Based PQC (Classic McEliece GF(2^12) Goppa + ISD Syndrome/Weight)
 * ============================================================================ */

static inline uint16_t pqc_ref_mceliece_goppa_inv(uint16_t alpha, const uint16_t coeffs[9], uint32_t poly12) {
    /* Monic degree-8 Horner: g(alpha) = (((...((coeffs[8]*a + coeffs[7])*a + ...))) */
    uint16_t acc = (uint16_t)(coeffs[8] & 0xFFFu);
    for (int32_t i = 7; i >= 0; --i) {
        acc = (uint16_t)(qir__gf2m_mul_u32(acc, alpha, 12u, poly12) ^ (coeffs[i] & 0xFFFu));
    }
    return (uint16_t)qir__gf2m_inv_u32(acc, 12u, poly12);
}

static inline bool pqc_build_mceliece_goppa_qir(
    qir_module_t *m,
    const uint16_t coeffs[9],
    uint32_t poly12,
    uint16_t target_inv,
    bool phase_oracle
) {
    uint32_t alpha = qir_add_input(m, 12u, "alpha");
    if (phase_oracle) qir_add_output(m, alpha);

    uint32_t acc = qir_emit_const_u64(m, 12u, coeffs[8] & 0xFFFu, "c8");
    for (int32_t i = 7; i >= 0; --i) {
        uint32_t prod = qir_emit_gf2m_mul(m, acc, alpha, 12u, poly12, "gp");
        if ((coeffs[i] & 0xFFFu) != 0u) {
            uint32_t ci = qir_emit_const_u64(m, 12u, coeffs[i] & 0xFFFu, "ci");
            acc = qir_emit_gf2m_add(m, prod, ci, 12u, "gs");
        } else {
            acc = prod;
        }
    }
    uint32_t inv = qir_emit_gf2m_inv(m, acc, 12u, poly12, "g_inv");
    if (phase_oracle) {
        uint8_t tgt[2] = { (uint8_t)(target_inv & 0xFFu), (uint8_t)((target_inv >> 8u) & 0x0Fu) };
        pqc_emit_grover_eq_target(m, inv, tgt, 12u);
    } else {
        qir_add_output(m, inv);
    }
    return true;
}

static inline void pqc_ref_mceliece_isd(
    uint64_t err_vec,
    const uint64_t h_rows[24],
    uint32_t *out_syn,
    uint8_t *out_wt
) {
    uint32_t syn = 0u;
    for (uint32_t r = 0u; r < 24u; ++r) {
        if (__builtin_parityll((unsigned long long)(h_rows[r] & err_vec))) {
            syn |= (UINT32_C(1) << r);
        }
    }
    *out_syn = syn;
    *out_wt  = (uint8_t)__builtin_popcountll((unsigned long long)err_vec);
}

static inline bool pqc_build_mceliece_isd_qir(
    qir_module_t *m,
    const uint64_t h_rows[24],
    uint32_t target_syn,
    uint8_t target_wt,
    bool phase_oracle
) {
    uint32_t e_in = qir_add_input(m, 64u, "err_vec");
    if (phase_oracle) qir_add_output(m, e_in);

    qcc_bv256_t rows[24];
    for (uint32_t r = 0u; r < 24u; ++r) {
        rows[r] = qcc_bv_zero();
        rows[r].w[0] = h_rows[r];
    }
    qcc_bv256_t zero_c = qcc_bv_zero();
    uint32_t syn = qir_emit_linf2(m, &e_in, 1u, 24u, rows, &zero_c, "syn");

    /* Reversible popcount adder tree over 64 bits -> 7-bit Hamming weight */
    uint32_t lvl[64];
    for (uint32_t i = 0u; i < 64u; ++i) lvl[i] = qir_emit_slice(m, e_in, (uint16_t)i, (uint16_t)(i + 1u), "eb");
    uint32_t count = 64u;
    while (count > 1u) {
        uint32_t nxt = count >> 1u;
        uint32_t zero1 = qir_emit_const_u64(m, 1u, 0u, "z1");
        for (uint32_t j = 0u; j < nxt; ++j) {
            uint32_t ext_a[2] = {lvl[2u * j], zero1};
            uint32_t ext_b[2] = {lvl[2u * j + 1u], zero1};
            uint32_t pa = qir_emit_concat(m, ext_a, 2u, "pa");
            uint32_t pb = qir_emit_concat(m, ext_b, 2u, "pb");
            lvl[j] = qir_emit_add(m, pa, pb, 0u, "pcnt");
        }
        count = nxt;
    }
    uint32_t wt = lvl[0]; /* 7 bits */

    if (phase_oracle) {
        uint32_t pair[2] = {syn, wt};
        uint32_t chk = qir_emit_concat(m, pair, 2u, "isd_chk");
        uint8_t tgt[4] = {
            (uint8_t)(target_syn & 0xFFu),
            (uint8_t)((target_syn >> 8u) & 0xFFu),
            (uint8_t)((target_syn >> 16u) & 0xFFu),
            (uint8_t)(target_wt & 0x7Fu)
        };
        pqc_emit_grover_eq_target(m, chk, tgt, 31u);
    } else {
        qir_add_output(m, syn);
        qir_add_output(m, wt);
    }
    return true;
}

/* ============================================================================
 * Family 4: Lattice-Based PQC (FIPS-203 ML-KEM / Kyber q = 3329 Module-LWE Search)
 * ============================================================================ */

#define PQC_KYBER_Q 3329u

static inline uint32_t pqc_ref_mlkem_lwe(
    uint32_t s_packed24,        /* 8 x 3-bit signed secrets in {-2,-1,0,1,2} encoded as val + 2 (0..4) */
    const uint16_t a_mat[4][8], /* Entries in [0, 3328] */
    const uint16_t b_vec[4],    /* Entries in [0, 3328] */
    uint16_t out_res[4]
) {
    uint16_t worst = 0u;
    for (uint32_t i = 0u; i < 4u; ++i) {
        int32_t acc = 0;
        for (uint32_t j = 0u; j < 8u; ++j) {
            uint32_t code = (s_packed24 >> (3u * j)) & 7u;
            int32_t sj = (code <= 4u) ? ((int32_t)code - 2) : 0;
            acc += (int32_t)a_mat[i][j] * sj;
        }
        int32_t diff = ((int32_t)b_vec[i] - acc) % (int32_t)PQC_KYBER_Q;
        if (diff < 0) diff += (int32_t)PQC_KYBER_Q;
        out_res[i] = (uint16_t)diff;
        int32_t centered = (diff > (int32_t)(PQC_KYBER_Q / 2u)) ? ((int32_t)PQC_KYBER_Q - diff) : diff;
        if ((uint16_t)centered > worst) worst = (uint16_t)centered;
    }
    return worst;
}

/* Reversible (a + b) mod 3329 -> 12-bit clean residue in [0, 3328] */
static inline uint32_t pqc__emit_mod3329_add(qir_module_t *m, uint32_t a12, uint32_t b12) {
    uint32_t z1 = qir_emit_const_u64(m, 1u, 0u, "z1");
    uint32_t ap[2] = {a12, z1};
    uint32_t bp[2] = {b12, z1};
    uint32_t a13 = qir_emit_concat(m, ap, 2u, "a13");
    uint32_t b13 = qir_emit_concat(m, bp, 2u, "b13");
    uint32_t sum = qir_emit_add(m, a13, b13, 0u, "sum13");
    uint32_t qc  = qir_emit_const_u64(m, 13u, PQC_KYBER_Q, "q13");
    uint32_t sub = qir_emit_sub(m, sum, qc, 0u, "sub13");
    uint32_t neg = qir_emit_slice(m, sub, 12u, 13u, "neg");
    uint32_t restored = qir_emit_add(m, sub, qc, 0u, "addq");
    uint32_t sel = qir_emit_ctrl_sel(m, neg, sub, restored, 0u, "modsum");
    return qir_emit_slice(m, sel, 0u, 12u, "r12");
}

static inline void pqc__fn_flip1(const uint8_t *in_b, uint8_t *out_b, const void *ctx_p) {
    (void)ctx_p;
    out_b[0] = (uint8_t)((in_b[0] & 1u) ^ 1u);
}

static inline bool pqc_build_mlkem_lwe_qir(
    qir_module_t *m,
    const uint16_t a_mat[4][8],
    const uint16_t b_vec[4],
    uint16_t eta,
    bool phase_oracle
) {
    if (eta >= 100u) return false;
    uint32_t s_in = qir_add_input(m, 24u, "s_cbd");
    if (phase_oracle) qir_add_output(m, s_in);

    uint32_t res_vals[4];
    uint32_t bad_bits[4];
    for (uint32_t i = 0u; i < 4u; ++i) {
        uint32_t acc = UINT32_MAX;
        for (uint32_t j = 0u; j < 8u; ++j) {
            uint32_t sj = qir_emit_slice(m, s_in, (uint16_t)(3u * j), (uint16_t)(3u * j + 3u), "sj");
            uint64_t lut8[8];
            for (uint32_t code = 0u; code < 8u; ++code) {
                int32_t sv = (code <= 4u) ? ((int32_t)code - 2) : 0;
                int32_t term = (-(int32_t)a_mat[i][j] * sv) % (int32_t)PQC_KYBER_Q;
                if (term < 0) term += (int32_t)PQC_KYBER_Q;
                lut8[code] = (uint64_t)term;
            }
            uint32_t prod = qir_emit_lut_rom(m, sj, 12u, 8u, lut8, 0u, "as_j");
            acc = (j == 0u) ? prod : pqc__emit_mod3329_add(m, acc, prod);
        }
        uint32_t b_c  = qir_emit_const_u64(m, 12u, b_vec[i] % PQC_KYBER_Q, "bi");
        uint32_t err  = pqc__emit_mod3329_add(m, acc, b_c);
        res_vals[i] = err;

        if (phase_oracle) {
            uint32_t eta_c = qir_emit_const_u64(m, 12u, eta, "eta");
            uint32_t sh    = pqc__emit_mod3329_add(m, err, eta_c);
            uint32_t z1    = qir_emit_const_u64(m, 1u, 0u, "z1");
            uint32_t sp[2] = {sh, z1};
            uint32_t sh13  = qir_emit_concat(m, sp, 2u, "sh13");
            uint32_t win_c = qir_emit_const_u64(m, 13u, (uint64_t)(2u * eta + 1u), "win");
            uint32_t cmp   = qir_emit_sub(m, sh13, win_c, 0u, "cmp");
            uint32_t okbit = qir_emit_slice(m, cmp, 12u, 13u, "okb");
            bad_bits[i] = pqc_emit_linf2_from_fn(m, &okbit, 1u, 1u, 1u, pqc__fn_flip1, NULL, "bad");
        }
    }

    if (phase_oracle) {
        uint32_t chk = qir_emit_concat(m, bad_bits, 4u, "lwe_bad");
        uint8_t zero_tgt[1] = {0u};
        pqc_emit_grover_eq_target(m, chk, zero_tgt, 4u);
    } else {
        uint32_t out48 = qir_emit_concat(m, res_vals, 4u, "lwe_res");
        qir_add_output(m, out48);
    }
    return true;
}

/* ============================================================================
 * Family 5: Multivariate Quadratic PQC (MQ n=24, m=24 Core — UOV / MAYO Style)
 * ============================================================================ */

typedef struct {
    uint32_t quad_coeff[24][12]; /* For eq k, bit i of quad_coeff[k][d-1] is Q_{k, i, (i+d)%24} */
    uint32_t lin_coeff[24];      /* 24-bit linear terms L_k */
    uint32_t const_bits;         /* Bit k is affine constant c_k */
} pqc_mq24_sys_t;

static inline uint32_t pqc_ref_mq24(uint32_t x24, const pqc_mq24_sys_t *sys) {
    uint32_t out = sys->const_bits;
    for (uint32_t k = 0u; k < 24u; ++k) {
        uint32_t acc = __builtin_parity((unsigned)(sys->lin_coeff[k] & x24));
        for (uint32_t d = 1u; d <= 12u; ++d) {
            uint32_t x_rot = ((x24 >> d) | (x24 << (24u - d))) & UINT32_C(0xFFFFFF);
            uint32_t p_d   = x24 & x_rot;
            uint32_t mask  = (d == 12u) ? UINT32_C(0xFFF) : UINT32_C(0xFFFFFF);
            acc ^= (uint32_t)__builtin_parity((unsigned)(sys->quad_coeff[k][d - 1u] & p_d & mask));
        }
        if (acc & 1u) out ^= (UINT32_C(1) << k);
    }
    return out;
}

static inline void pqc__mq_fn_rot24(const uint8_t *in_b, uint8_t *out_b, const void *ctx_p) {
    uint32_t d = *(const uint32_t *)ctx_p;
    uint32_t x = (uint32_t)in_b[0] | ((uint32_t)in_b[1] << 8u) | ((uint32_t)in_b[2] << 16u);
    uint32_t y = ((x >> d) | (x << (24u - d))) & UINT32_C(0xFFFFFF);
    out_b[0] = (uint8_t)(y & 0xFFu);
    out_b[1] = (uint8_t)((y >> 8u) & 0xFFu);
    out_b[2] = (uint8_t)((y >> 16u) & 0xFFu);
}

static inline bool pqc_build_mq24_qir(
    qir_module_t *m,
    const pqc_mq24_sys_t *sys,
    uint32_t target_y24,
    bool phase_oracle
) {
    uint32_t x_in = qir_add_input(m, 24u, "x");
    if (phase_oracle) qir_add_output(m, x_in);

    uint32_t prods[12];
    for (uint32_t d = 1u; d <= 12u; ++d) {
        uint32_t xr = pqc_emit_linf2_from_fn(m, &x_in, 1u, 24u, 24u, pqc__mq_fn_rot24, &d, "xr");
        prods[d - 1u] = qir_emit_and_vec(m, x_in, xr, 0u, "pd");
    }

    /* Group 1: x_in (24b) + prods[0..8] (9*24=216b) -> 240b <= 256b */
    uint32_t g1_ins[10];
    g1_ins[0] = x_in;
    for (uint32_t i = 0u; i < 9u; ++i) g1_ins[1u + i] = prods[i];
    qcc_bv256_t rows1[24], rows2[24];
    qcc_bv256_t c1 = qcc_bv_zero(), c2 = qcc_bv_zero();
    for (uint32_t k = 0u; k < 24u; ++k) {
        rows1[k] = qcc_bv_zero();
        rows2[k] = qcc_bv_zero();
        rows1[k].w[0] = sys->lin_coeff[k] & UINT64_C(0xFFFFFF);
        for (uint32_t d = 1u; d <= 9u; ++d) {
            uint32_t off = 24u * d;
            uint64_t qc  = sys->quad_coeff[k][d - 1u] & UINT64_C(0xFFFFFF);
            for (uint32_t b = 0u; b < 24u; ++b) {
                if ((qc >> b) & 1u) qcc_bv_set(rows1[k].w, off + b);
            }
        }
        rows2[k].w[0] = (uint64_t)(sys->quad_coeff[k][9] & UINT32_C(0xFFFFFF)) |
                        ((uint64_t)(sys->quad_coeff[k][10] & UINT32_C(0xFFFFFF)) << 24u) |
                        ((uint64_t)(sys->quad_coeff[k][11] & UINT32_C(0xFFF)) << 48u);
        if ((sys->const_bits >> k) & 1u) qcc_bv_set(c1.w, k);
    }
    uint32_t p12_lo = qir_emit_slice(m, prods[11], 0u, 12u, "p12h");
    uint32_t g2_ins[3] = {prods[9], prods[10], p12_lo};
    uint32_t y1 = qir_emit_linf2(m, g1_ins, 10u, 24u, rows1, &c1, "mq1");
    uint32_t y2 = qir_emit_linf2(m, g2_ins, 3u,  24u, rows2, &c2, "mq2");
    uint32_t y  = qir_emit_gf2m_add(m, y1, y2, 24u, "mq_y");

    if (phase_oracle) {
        uint8_t tgt[3] = {
            (uint8_t)(target_y24 & 0xFFu),
            (uint8_t)((target_y24 >> 8u) & 0xFFu),
            (uint8_t)((target_y24 >> 16u) & 0xFFu)
        };
        pqc_emit_grover_eq_target(m, y, tgt, 24u);
    } else {
        qir_add_output(m, y);
    }
    return true;
}

#ifdef __cplusplus
}
#endif

#endif /* PQC_ORACLES_H */
