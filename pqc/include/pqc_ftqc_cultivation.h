/*
 * pqc_ftqc_cultivation.h — Full-Stack FTQC Co-Design Estimator:
 *   qctrl Ring-0 Calibration x Magic State Cultivation (Gidney-Shutty-Jones)
 *   x Surface Code Logical Routing & 15-to-1 Distillation Comparison.
 * Pure C17 & C++20 compatible, zero <math.h>, pure integer Q16.16 / ppb math.
 */
#ifndef PQC_FTQC_CULTIVATION_H
#define PQC_FTQC_CULTIVATION_H

#include "pqc_nist_est.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char     name[64];
    uint32_t f_2q_ppm;         /* e.g. 996000 for 99.60% */
    uint32_t t1_us;            /* e.g. 530 us */
    uint32_t t2_us;            /* e.g. 300 us */
    uint32_t patch_rms_mrad;   /* from qctrl bench40 */
    uint32_t seam_rms_mrad;    /* from qctrl bench40 */
    uint32_t acts_per_6700;    /* calibration interventions */
    uint32_t ring0_ns_sweep;   /* 0 for host/offline, ~12000 for Ring-0 */
} pqc_hw_cal_profile_t;

typedef struct {
    uint64_t p_incoh_ppb;         /* 2Q + T1/T2 base physical error in ppb (1e-9) */
    uint64_t p_coh_zz_ppb;        /* coherent ZZ error floor (phi/2)^2 * 2.67 in ppb */
    uint64_t p_eff_ppb;           /* combined effective physical cycle error in ppb */
    bool     below_sc_threshold;  /* p_eff < 1.0e-2 (10,000,000 ppb) */
    uint32_t d_data;              /* minimal odd rotated SC distance achieving <=1% block err */
    uint64_t data_phys_qubits;    /* 2 * d^2 * Q_logical */
    bool     cult_viable;         /* true if Cultivation achieves <= 1%/(2*N_T) error budget */
    uint32_t cult_d;              /* cultivation distance 3 or 5 */
    uint64_t cult_eps_ppb;        /* cultivated |T> residual infidelity in ppb */
    uint32_t cult_p_surv_pct;     /* post-selection survival probability (%) */
    uint64_t cult_factory_qubits; /* qubits for parallelized Cultivation factories */
    uint64_t cult_spacetime_qcyc; /* total physical qubit-cycles with Cultivation */
    uint64_t dist_factory_qubits; /* qubits for 15-to-1 distillation factory */
    uint64_t dist_spacetime_qcyc; /* total physical qubit-cycles with 15-to-1 Distillation */
    uint32_t vol_reduction_q8;    /* (dist_spacetime / cult_spacetime) * 256 in Q8 */
} pqc_ftqc_design_t;

typedef struct {
    uint32_t ring0_cal_sweep_ns;
    uint32_t decoder_d3_ns_x100;
    uint32_t decoder_d5_ns_x100;
    uint32_t decoder_d7_ns_x100;
    uint32_t lambda_d3_to_d5_x100;
    uint32_t lambda_d5_to_d7_x100;
} pqc_decoder_summary_t;

/* Pure integer computation of p_eff in ppb (1 = 10^-9) */
static inline void pqc_ftqc_eval_noise_ppb(
    const pqc_hw_cal_profile_t *prof,
    uint64_t *out_incoh_ppb,
    uint64_t *out_coh_ppb,
    uint64_t *out_eff_ppb
) {
    /* 1. 2Q gate infidelity + T1/T2 idle per 1us cycle (exact match to qctrl_breakeven.h) */
    uint64_t p_2q_ppb   = (prof->f_2q_ppm < 1000000u) ? ((uint64_t)(1000000u - prof->f_2q_ppm) * 1000ULL) : 0ULL;
    uint64_t p_t1t2_ppb = (prof->t2_us > 0u ? (500000000ULL / (uint64_t)prof->t2_us) : 0ULL) +
                          (prof->t1_us > 0u ? (250000000ULL / (uint64_t)prof->t1_us) : 0ULL);
    uint64_t p_incoh    = p_2q_ppb + p_t1t2_ppb;

    /* 2. Coherent ZZ error floor: ~ 2.67 couplers/data * (patch_rms / 2)^2 */
    uint64_t mrad       = (uint64_t)prof->patch_rms_mrad;
    uint64_t p_coh      = (mrad * mrad * 2670ULL) / 4ULL;

    if (out_incoh_ppb) *out_incoh_ppb = p_incoh;
    if (out_coh_ppb)   *out_coh_ppb   = p_coh;
    if (out_eff_ppb)   *out_eff_ppb   = p_incoh + p_coh;
}

/* Fixed-point log2 of logical error rate eps_L(d) ≈ 0.1 * (p_eff / p_th)^((d+1)/2)
 * where p_th = 0.01 (10,000,000 ppb). Returns negative log2(1 / eps_L) in Q16.16 (pos = safer).
 */
static inline uint32_t pqc_ftqc_neg_log2_eps_l_q16(uint64_t p_eff_ppb, uint32_t d) {
    const uint64_t p_th_ppb = 10000000ULL; /* 1.0% threshold */
    if (p_eff_ppb == 0ULL) return 200u << 16u;
    if (p_eff_ppb >= p_th_ppb) return 0u; /* Above threshold: distance gives 0 gain */

    /* log2(p_th / p_eff) in Q16.16 via ratio scaled by 2^32 */
    uint64_t ratio_q32 = (p_th_ppb << 32u) / p_eff_ppb;
    uint32_t l_ratio_q16 = pqc_log2_q16_u64(ratio_q32) - (32u << 16u);
    uint32_t exp_pow = (d + 1u) >> 1u;
    /* -log2(0.1) ≈ 3.32193 in Q16.16 = 217706 */
    return 217706u + exp_pow * l_ratio_q16;
}

/* Pure integer 100-bit-safe estimation of total qubits and qubit-cycles for Cultivation vs 15-to-1 */
static inline void pqc_ftqc_codesign(
    const pqc_hw_cal_profile_t *prof,
    uint32_t logical_qubits,
    uint64_t t_count_circuit,
    uint64_t t_depth_circuit,
    pqc_ftqc_design_t *out
) {
    memset(out, 0, sizeof(*out));
    pqc_ftqc_eval_noise_ppb(prof, &out->p_incoh_ppb, &out->p_coh_zz_ppb, &out->p_eff_ppb);

    const uint64_t p_th_ppb = 10000000ULL;
    out->below_sc_threshold = (out->p_eff_ppb < p_th_ppb);

    uint64_t q_l   = logical_qubits ? (uint64_t)logical_qubits : 1ULL;
    uint64_t n_t   = t_count_circuit ? t_count_circuit : 1ULL;
    uint64_t d_t   = t_depth_circuit ? t_depth_circuit : 1ULL;

    /* Required -log2(eps_L_round) >= log2(200 * q_l * d_t) so total memory error <= 0.5% */
    uint64_t vol_q64 = (q_l <= (UINT64_MAX / 200ULL / (d_t ? d_t : 1ULL)))
                         ? (200ULL * q_l * d_t)
                         : UINT64_MAX;
    uint32_t req_neg_log2_q16 = pqc_log2_q16_u64(vol_q64);

    if (!out->below_sc_threshold) {
        out->d_data = 0u; /* Cannot stabilize at any distance without calibration */
        return;
    }

    uint32_t best_d = 35u;
    for (uint32_t d = 3u; d <= 35u; d += 2u) {
        if (pqc_ftqc_neg_log2_eps_l_q16(out->p_eff_ppb, d) >= req_neg_log2_q16) {
            best_d = d;
            break;
        }
    }
    out->d_data = best_d;

    /* Logical block footprint including lattice surgery routing overhead (2 * d^2 per logical qubit) */
    uint64_t phys_per_lq = 2ULL * (uint64_t)best_d * (uint64_t)best_d;
    out->data_phys_qubits = phys_per_lq * q_l;

    /* Cultivation model (Gidney-Shutty-Jones 2024 empirical scaling anchored at p=1e-3):
     *   d_c = 3: eps_cult ≈ 4e-7 * (p_eff / 1e-3)^4, footprint ~ 4*d_c^2 = 36q
     *   d_c = 5: eps_cult ≈ 2e-9 * (p_eff / 1e-3)^6, footprint ~ 4*d_c^2 = 100q
     *   Threshold cutoff for post-selection stability: p_eff <= 6.0e-3 (d_c=5) / 5.5e-3 (d_c=3)
     */
    const uint64_t p_ref_ppb = 1000000ULL; /* 1e-3 */
    flow_u128 r_q12 = ((flow_u128)out->p_eff_ppb << 12u) / (flow_u128)p_ref_ppb;
    flow_u128 r2    = (r_q12 * r_q12) >> 12u;
    flow_u128 r4    = (r2 * r2) >> 12u;
    flow_u128 r6    = (r4 * r2) >> 12u;
    uint64_t eps_c3_ppb = (uint64_t)((400ULL * r4) >> 12u);
    uint64_t eps_c5_ppb = (uint64_t)((2ULL * r6) >> 12u);
    if (eps_c3_ppb == 0ULL) eps_c3_ppb = 1ULL;
    if (eps_c5_ppb == 0ULL) eps_c5_ppb = 1ULL;

    /* 15-to-1 Distillation baseline (Litinski T-factory at d_dist = max(3, (d_data/2) | 1)) */
    uint32_t d_dist = (best_d / 2u);
    if (d_dist < 3u) d_dist = 3u;
    if ((d_dist & 1u) == 0u) d_dist++;
    uint64_t dist_unit_q   = 15ULL * (uint64_t)d_dist * (uint64_t)d_dist;
    uint64_t dist_unit_cyc = 6ULL * (uint64_t)d_dist;

    uint64_t total_exec_cyc = d_t * (uint64_t)best_d;
    uint64_t par_dist = (n_t * dist_unit_cyc + total_exec_cyc - 1ULL) / total_exec_cyc;
    if (par_dist < 1ULL) par_dist = 1ULL;
    if (par_dist > 16ULL) par_dist = 16ULL;

    out->dist_factory_qubits = par_dist * dist_unit_q;
    uint64_t dist_cyc_tot    = total_exec_cyc;
    if (n_t * dist_unit_cyc > par_dist * total_exec_cyc) {
        dist_cyc_tot = (n_t * dist_unit_cyc) / par_dist;
    }
    out->dist_spacetime_qcyc = (out->data_phys_qubits + out->dist_factory_qubits) * dist_cyc_tot;

    /* Required per-T error <= 0.5% / n_t (i.e. eps_T * 200 * n_t <= 1e9 ppb) */
    uint64_t max_eps_t_ppb = (n_t <= (1000000000ULL / 200ULL)) ? (1000000000ULL / (200ULL * n_t)) : 1ULL;
    if (max_eps_t_ppb == 0ULL) max_eps_t_ppb = 1ULL;

    uint32_t chosen_dc = 0u;
    if (eps_c3_ppb <= max_eps_t_ppb && out->p_eff_ppb <= 5500000ULL) {
        chosen_dc = 3u;
        out->cult_eps_ppb = eps_c3_ppb;
    } else if (eps_c5_ppb <= max_eps_t_ppb && out->p_eff_ppb <= 6500000ULL) {
        chosen_dc = 5u;
        out->cult_eps_ppb = eps_c5_ppb;
    }

    if (chosen_dc == 0u) {
        out->cult_viable         = false;
        out->cult_d              = 0u;
        out->cult_factory_qubits = out->dist_factory_qubits;
        out->cult_spacetime_qcyc = out->dist_spacetime_qcyc;
        out->vol_reduction_q8    = 256u;
        return;
    }

    out->cult_viable = true;
    out->cult_d      = chosen_dc;
    uint64_t penalty = ((chosen_dc == 3u ? 95ULL : 140ULL) * out->p_eff_ppb) / 10000000ULL;
    uint32_t p_surv  = (penalty < 92ULL) ? (uint32_t)(100ULL - penalty) : 8u;
    out->cult_p_surv_pct = p_surv;

    /* Cultivated state unit footprint ~ 4 * d_c^2, expected cycles = (4 * d_c * 100) / p_surv */
    uint64_t cult_unit_q   = 4ULL * (uint64_t)chosen_dc * (uint64_t)chosen_dc;
    uint64_t cult_unit_cyc = ((4ULL * (uint64_t)chosen_dc * 100ULL) + (uint64_t)p_surv - 1ULL) / (uint64_t)p_surv;
    if (cult_unit_cyc < 6ULL) cult_unit_cyc = 6ULL;

    uint64_t par_cult = (n_t * cult_unit_cyc + total_exec_cyc - 1ULL) / total_exec_cyc;
    if (par_cult < 1ULL) par_cult = 1ULL;
    if (par_cult > 32ULL) par_cult = 32ULL;

    out->cult_factory_qubits = par_cult * cult_unit_q;
    uint64_t cult_cyc_tot    = total_exec_cyc;
    if (n_t * cult_unit_cyc > par_cult * total_exec_cyc) {
        cult_cyc_tot = (n_t * cult_unit_cyc) / par_cult;
    }
    out->cult_spacetime_qcyc = (out->data_phys_qubits + out->cult_factory_qubits) * cult_cyc_tot;
    if (out->cult_spacetime_qcyc > 0ULL) {
        uint64_t red = (out->dist_spacetime_qcyc * 256ULL) / out->cult_spacetime_qcyc;
        out->vol_reduction_q8 = (uint32_t)(red > 65535ULL ? 65535ULL : red);
    } else {
        out->vol_reduction_q8 = 256u;
    }
}

#ifdef __cplusplus
}
#endif

#endif /* PQC_FTQC_CULTIVATION_H */
