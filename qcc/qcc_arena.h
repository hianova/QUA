/*
 * qcc_arena.h — Zero-Fragmentation Linear Arena, 256-Bit Bitvectors & 128-Bit Phase Primitives
 * Part of QCC (Quantum Circuit Compiler) for FLOW v7.0.1.
 * Pure C17, zero external dependencies.
 */
#ifndef QCC_ARENA_H
#define QCC_ARENA_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "flow.h"

#ifdef __cplusplus
extern "C" {
#endif

#define QCC_ARENA_DEFAULT_SLAB_BYTES (256u * 1024u)

typedef struct qcc_arena_slab_s {
    struct qcc_arena_slab_s *next;
    size_t                   cap;
    size_t                   offset;
    uint8_t                 *data;
} qcc_arena_slab_t;

typedef struct {
    qcc_arena_slab_t *head;
    qcc_arena_slab_t *cur;
    size_t            default_slab_cap;
    size_t            total_allocated;
    /* Compatibility view of active slab matching plan.md fields */
    uint8_t          *base;
    size_t            cap;
    size_t            offset;
} qcc_arena_t;

typedef struct {
    qcc_arena_slab_t *slab;
    size_t            offset;
} qcc_arena_mark_t;

static inline qcc_arena_slab_t *qcc__arena_new_slab(size_t cap) {
    qcc_arena_slab_t *s = (qcc_arena_slab_t *)malloc(sizeof(qcc_arena_slab_t));
    if (!s) return NULL;
    s->data = (uint8_t *)malloc(cap);
    if (!s->data) {
        free(s);
        return NULL;
    }
    s->next   = NULL;
    s->cap    = cap;
    s->offset = 0u;
    return s;
}

static inline bool qcc_arena_init(qcc_arena_t *a, size_t initial_cap) {
    memset(a, 0, sizeof(*a));
    if (initial_cap == 0u) initial_cap = QCC_ARENA_DEFAULT_SLAB_BYTES;
    a->default_slab_cap = initial_cap;
    qcc_arena_slab_t *s = qcc__arena_new_slab(initial_cap);
    if (!s) return false;
    a->head   = s;
    a->cur    = s;
    a->base   = s->data;
    a->cap    = s->cap;
    a->offset = 0u;
    return true;
}

static inline void qcc_arena_free(qcc_arena_t *a) {
    if (!a) return;
    qcc_arena_slab_t *s = a->head;
    while (s) {
        qcc_arena_slab_t *nxt = s->next;
        free(s->data);
        free(s);
        s = nxt;
    }
    memset(a, 0, sizeof(*a));
}

static inline void qcc_arena_reset(qcc_arena_t *a) {
    if (!a || !a->head) return;
    qcc_arena_slab_t *s = a->head;
    while (s) {
        s->offset = 0u;
        s = s->next;
    }
    a->cur             = a->head;
    a->base            = a->head->data;
    a->cap             = a->head->cap;
    a->offset          = 0u;
    a->total_allocated = 0u;
}

static inline void *qcc_arena_alloc(qcc_arena_t *a, size_t bytes, size_t align) {
    if (!a) return NULL;
    if (align == 0u) align = 8u;
    if (bytes == 0u) bytes = 1u;
    if (!a->cur) {
        if (!qcc_arena_init(a, QCC_ARENA_DEFAULT_SLAB_BYTES)) return NULL;
    }

    size_t mask = align - 1u;
    for (;;) {
        qcc_arena_slab_t *s = a->cur;
        size_t aligned = (s->offset + mask) & ~mask;
        if (aligned <= s->cap && bytes <= s->cap - aligned) {
            void *ptr = s->data + aligned;
            s->offset = aligned + bytes;
            a->base   = s->data;
            a->cap    = s->cap;
            a->offset = s->offset;
            a->total_allocated += bytes;
            return ptr;
        }
        if (s->next) {
            a->cur        = s->next;
            a->cur->offset = 0u;
            continue;
        }
        size_t need = bytes + align;
        size_t ncap = a->default_slab_cap;
        if (need > ncap) ncap = need;
        qcc_arena_slab_t *ns = qcc__arena_new_slab(ncap);
        if (!ns) return NULL;
        s->next = ns;
        a->cur  = ns;
    }
}

static inline void *qcc_arena_alloc_zero(qcc_arena_t *a, size_t bytes, size_t align) {
    void *p = qcc_arena_alloc(a, bytes, align);
    if (p && bytes > 0u) memset(p, 0, bytes);
    return p;
}

static inline void *qcc_xmalloc(size_t n) {
    void *p = malloc(n ? n : 1u);
    if (!p) abort();
    return p;
}

static inline void *qcc_xcalloc(size_t cnt, size_t sz) {
    void *p = calloc(cnt ? cnt : 1u, sz ? sz : 1u);
    if (!p) abort();
    return p;
}

static inline void *qcc_xrealloc(void *ptr, size_t n) {
    void *p = realloc(ptr, n ? n : 1u);
    if (!p) abort();
    return p;
}

static inline void *qcc_arena_xalloc(qcc_arena_t *a, size_t bytes, size_t align) {
    void *p = qcc_arena_alloc(a, bytes, align);
    if (!p) abort();
    return p;
}

static inline void *qcc_arena_xalloc_zero(qcc_arena_t *a, size_t bytes, size_t align) {
    void *p = qcc_arena_alloc_zero(a, bytes, align);
    if (!p) abort();
    return p;
}

static inline char *qcc_arena_strdup(qcc_arena_t *a, const char *s) {
    if (!s) return NULL;
    size_t n = strlen(s) + 1u;
    char *p = (char *)qcc_arena_alloc(a, n, 1u);
    if (p) memcpy(p, s, n);
    return p;
}

static inline qcc_arena_mark_t qcc_arena_Save(const qcc_arena_t *a) {
    qcc_arena_mark_t m = {a ? a->cur : NULL, a && a->cur ? a->cur->offset : 0u};
    return m;
}

static inline void qcc_arena_Restore(qcc_arena_t *a, qcc_arena_mark_t mark) {
    if (!a || !mark.slab) return;
    qcc_arena_slab_t *s = mark.slab->next;
    while (s) {
        s->offset = 0u;
        s = s->next;
    }
    mark.slab->offset = mark.offset;
    a->cur            = mark.slab;
    a->base           = mark.slab->data;
    a->cap            = mark.slab->cap;
    a->offset         = mark.offset;
}

/* ============================================================================
 * 256-Bit Bitvector Helpers (up to 256 qubits/wires per register)
 * ============================================================================ */

typedef struct {
    uint64_t w[4];
} qcc_bv256_t;

static inline qcc_bv256_t qcc_bv_zero(void) {
    qcc_bv256_t r = {{0u, 0u, 0u, 0u}};
    return r;
}

static inline qcc_bv256_t qcc_bv_mask_n(uint32_t n) {
    qcc_bv256_t r = qcc_bv_zero();
    for (uint32_t i = 0u; i < 4u; ++i) {
        if (n >= 64u * (i + 1u)) {
            r.w[i] = UINT64_MAX;
        } else if (n > 64u * i) {
            uint32_t rem = n - 64u * i;
            r.w[i] = (UINT64_C(1) << rem) - 1u;
        }
    }
    return r;
}

static inline bool qcc_bv_get(const uint64_t w[4], uint32_t bit) {
    if (bit >= 256u) return false;
    return ((w[bit >> 6u] >> (bit & 63u)) & UINT64_C(1)) != 0u;
}

static inline void qcc_bv_set(uint64_t w[4], uint32_t bit) {
    if (bit < 256u) w[bit >> 6u] |= (UINT64_C(1) << (bit & 63u));
}

static inline void qcc_bv_clr(uint64_t w[4], uint32_t bit) {
    if (bit < 256u) w[bit >> 6u] &= ~(UINT64_C(1) << (bit & 63u));
}

static inline bool qcc_bv_is_zero(const uint64_t w[4]) {
    return (w[0] | w[1] | w[2] | w[3]) == 0u;
}

static inline uint32_t qcc_bv_popcount(const uint64_t w[4]) {
    return (uint32_t)(__builtin_popcountll((unsigned long long)w[0]) +
                      __builtin_popcountll((unsigned long long)w[1]) +
                      __builtin_popcountll((unsigned long long)w[2]) +
                      __builtin_popcountll((unsigned long long)w[3]));
}

/* ============================================================================
 * 128-Bit Full-Turn Phase & Deterministic PRNG Helpers
 * ============================================================================ */

typedef struct {
    uint64_t lo; /* words[0] / 2^128 full-turns */
    uint64_t hi; /* words[1] / 2^64  full-turns */
} qcc_angle128_t;

static inline flow_u128 qcc_u128_from_angle(qcc_angle128_t a) {
    return ((flow_u128)a.hi << 64u) | (flow_u128)a.lo;
}

static inline qcc_angle128_t qcc_angle_from_u128(flow_u128 turns) {
    qcc_angle128_t a;
    a.lo = (uint64_t)turns;
    a.hi = (uint64_t)(turns >> 64u);
    return a;
}

static inline uint8_t qcc_u128_clifford_level(flow_u128 turns) {
    uint64_t lo = (uint64_t)turns;
    uint64_t hi = (uint64_t)(turns >> 64u);
    if (lo != 0u) return (uint8_t)(128u - (uint32_t)__builtin_ctzll((unsigned long long)lo));
    if (hi != 0u) return (uint8_t)(64u - (uint32_t)__builtin_ctzll((unsigned long long)hi));
    return 0u;
}

static inline bool qcc_u128_is_mult_180(flow_u128 turns) {
    uint64_t lo = (uint64_t)turns;
    uint64_t hi = (uint64_t)(turns >> 64u);
    return (lo == 0u) && ((hi & UINT64_C(0x7FFFFFFFFFFFFFFF)) == 0u);
}

static inline bool qcc_u128_is_mult_90(flow_u128 turns) {
    uint64_t lo = (uint64_t)turns;
    uint64_t hi = (uint64_t)(turns >> 64u);
    return (lo == 0u) && ((hi & UINT64_C(0x3FFFFFFFFFFFFFFF)) == 0u);
}

static inline bool qcc_u128_is_mult_45(flow_u128 turns) {
    uint64_t lo = (uint64_t)turns;
    uint64_t hi = (uint64_t)(turns >> 64u);
    return (lo == 0u) && ((hi & UINT64_C(0x1FFFFFFFFFFFFFFF)) == 0u);
}

static inline uint64_t qcc_splitmix64(uint64_t *state) {
    uint64_t z = (*state += UINT64_C(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30u)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27u)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31u);
}

static inline uint32_t qcc_ceil_log2_u32(uint32_t v) {
    if (v <= 1u) return 0u;
    uint32_t b = 0u;
    uint32_t x = v - 1u;
    while (x > 0u) {
        x >>= 1u;
        ++b;
    }
    return b;
}

#ifdef __cplusplus
}
#endif

#endif /* QCC_ARENA_H */
