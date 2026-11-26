#define QCTRL_EXPORT_RING0 1
#include "qctrl_ring0.h"
#include "qctrl_decoder_ring0.h"

bool qctrl_edge_step_ring0_export(qctrl_edge_state_t *st, const qctrl_cfg_t *cfg, int32_t km,
                                  int32_t k0, int32_t kp, int32_t n_shots) {
    return qctrl_edge_step_ring0(st, cfg, km, k0, kp, n_shots);
}

uint32_t qctrl_grid40_color_sweep_ring0(qctrl_edge_state_t states[QCTRL_NUM_EDGES],
                                        const qctrl_edge_desc_t edges[QCTRL_NUM_EDGES],
                                        const qctrl_cfg_t *cfg, uint8_t target_color,
                                        const int32_t km[QCTRL_NUM_EDGES],
                                        const int32_t k0[QCTRL_NUM_EDGES],
                                        const int32_t kp[QCTRL_NUM_EDGES], int32_t n_shots) {
    uint32_t acts = 0U;
    uint32_t e;
    for (e = 0U; e < QCTRL_NUM_EDGES; ++e) {
        if (edges[e].color == target_color) {
            if (qctrl_edge_step_ring0(&states[e], cfg, km[e], k0[e], kp[e], n_shots)) {
                acts += 1U;
            }
        }
    }
    return acts;
}

bool qctrl_sc_init_ring0_export(qctrl_sc_code_t *code, uint32_t d) {
    return qctrl_sc_init_ring0(code, d);
}

uint64_t qctrl_sc_decode_ring0_export(const qctrl_sc_code_t *code, uint32_t synd, bool is_z_err) {
    return qctrl_sc_decode_ring0(code, synd, is_z_err);
}

bool qctrl_synd_monitor_step_ring0_export(qctrl_synd_monitor_t *mon, uint32_t synd_x,
                                          uint32_t synd_z, uint32_t n_stab_per_type) {
    return qctrl_synd_monitor_step_ring0(mon, synd_x, synd_z, n_stab_per_type);
}
