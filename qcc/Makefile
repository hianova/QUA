export PATH := $(HOME)/.elan/bin:/opt/homebrew/opt/llvm/bin:/opt/homebrew/bin:/usr/local/bin:$(PATH)

FLOW_DIR ?= ../../FLOW
CC       := cc
CXX      ?= c++
LEAN     ?= lean
CFLAGS   := -std=c17 -O3 -Wall -Wextra -Werror -I$(FLOW_DIR)/core -I$(FLOW_DIR)/sdk
BUILD    := build
TARGET   := $(BUILD)/qcc
HDRS     := qcc_arena.h qcc_kmb.h qir.h qir_lift.h qir_rewrite.h \
            qir_kernels.h qir_codegen.h qir_opt.h qir_pass_kmb.h qir_verify.h qir_fir.h

.PHONY: all clean check-headers test lean-fir examples

all: $(TARGET)

$(BUILD):
	mkdir -p $(BUILD)

$(TARGET): qcc_main.c $(HDRS) $(FLOW_DIR)/core/flow.h $(FLOW_DIR)/core/searu.h $(FLOW_DIR)/core/flow_opt_params.h $(FLOW_DIR)/sdk/flow_ir.h | $(BUILD)
	$(CC) $(CFLAGS) -o $@ qcc_main.c

check-headers:
	$(CC) -std=c17 -Wall -Wextra -Werror -Wno-unused-function -I$(FLOW_DIR)/core -I$(FLOW_DIR)/sdk -I. -fsyntax-only qir_verify.h qir_pass_kmb.h qir_fir.h
	$(CXX) -std=c++20 -Wall -Wextra -Werror -Wno-unused-function -iquote $(FLOW_DIR)/core -iquote $(FLOW_DIR)/sdk -iquote . -x c++-header -fsyntax-only qir_verify.h qir_pass_kmb.h qir_fir.h

test: $(TARGET) check-headers
	./$(TARGET) selftest

lean-fir: $(TARGET)
	@if grep -Eq '^[[:space:]]*(theorem|lemma|axiom)[[:space:]]' ../Qua.lean; then \
		echo "ERROR: Qua.lean violates SSOT pure-config rule (found theorem)!" >&2; exit 1; \
	fi
	@qua_lines="$$(wc -l < ../Qua.lean | tr -d '[:space:]')"; \
	if [ "$$qua_lines" -gt 120 ]; then echo "ERROR: Qua.lean exceeds 120 lines ($$qua_lines)" >&2; exit 1; fi
	$(MAKE) -C "$(FLOW_DIR)" build/lean/Flow.olean build/flow
	LEAN_PATH="$(FLOW_DIR)/build/lean" $(LEAN) -DwarningAsError=true --run ../Qua.lean $(BUILD)
	"$(FLOW_DIR)/build/flow" fir-synth $(BUILD)/gf256_tower_spec.fir $(BUILD)/gf256_tower_a32.fir --name gf256_tower_a32
	"$(FLOW_DIR)/build/flow" fir-synth $(BUILD)/gf256_tower_spec.fir $(BUILD)/gf256_tower_d4.fir --name gf256_tower_d4 --max-depth 2
	"$(FLOW_DIR)/build/flow" fir-synth $(BUILD)/lut5_bijection_spec.fir $(BUILD)/lut5_bijection_opt.fir --name lut5_bijection_opt
	./$(TARGET) opt $(BUILD)/gf256_sbox_mac.fir -o $(BUILD)/gf256_sbox_mac.kmb --clean-anc
	./$(TARGET) verify $(BUILD)/gf256_sbox_mac.kmb --ref $(BUILD)/gf256_sbox_mac.fir --states 4096
	./$(TARGET) opt $(BUILD)/rev_alu_fxp.fir -o $(BUILD)/rev_alu_fxp.kmb --clean-anc
	./$(TARGET) verify $(BUILD)/rev_alu_fxp.kmb --ref $(BUILD)/rev_alu_fxp.fir --states 4096
	./$(TARGET) opt $(BUILD)/gf256_tower_spec.fir -o $(BUILD)/gf256_tower_spec.kmb --clean-anc
	./$(TARGET) verify $(BUILD)/gf256_tower_spec.kmb --ref $(BUILD)/gf256_tower_spec.fir --states 4096
	./$(TARGET) opt $(BUILD)/gf256_tower_a32.fir -o $(BUILD)/gf256_tower_a32.kmb --clean-anc
	./$(TARGET) verify $(BUILD)/gf256_tower_a32.kmb --ref $(BUILD)/gf256_tower_spec.fir --states 4096
	./$(TARGET) opt $(BUILD)/gf256_tower_d4.fir -o $(BUILD)/gf256_tower_d4.kmb --clean-anc
	./$(TARGET) verify $(BUILD)/gf256_tower_d4.kmb --ref $(BUILD)/gf256_tower_spec.fir --states 4096
	./$(TARGET) opt $(BUILD)/lut5_bijection_spec.fir -o $(BUILD)/lut5_bijection_spec.kmb --clean-anc
	./$(TARGET) verify $(BUILD)/lut5_bijection_spec.kmb --ref $(BUILD)/lut5_bijection_spec.fir --states 4096
	./$(TARGET) opt $(BUILD)/lut5_bijection_opt.fir -o $(BUILD)/lut5_bijection_opt.kmb --clean-anc
	./$(TARGET) verify $(BUILD)/lut5_bijection_opt.kmb --ref $(BUILD)/lut5_bijection_spec.fir --states 4096
	$(MAKE) -C ../pqc FLOW_DIR=$(abspath $(FLOW_DIR)) build/pqc
	../pqc/build/pqc analyze $(BUILD)/gf256_tower_spec.kmb --bits 8 --name GF256_TOWER_SPEC
	../pqc/build/pqc analyze $(BUILD)/gf256_tower_a32.kmb --bits 8 --name GF256_TOWER_SEARU_A32
	../pqc/build/pqc analyze $(BUILD)/gf256_tower_d4.kmb --bits 8 --name GF256_TOWER_SEARU_D4
	../pqc/build/pqc analyze $(BUILD)/lut5_bijection_spec.kmb --bits 5 --name LUT5_BIJ_SPEC
	../pqc/build/pqc analyze $(BUILD)/lut5_bijection_opt.kmb --bits 5 --name LUT5_BIJ_SEARU_OPT
	@if ./$(TARGET) opt $(BUILD)/cartesian_reject.fir -o $(BUILD)/bad_cart.kmb 2>/dev/null; then \
		echo "ERROR: cartesian_reject.fir was NOT rejected by quantum compiler!" >&2; exit 1; \
	fi
	@for f in gf256_sbox_mac.fir gf256_tower_spec.fir gf256_tower_a32.fir gf256_tower_d4.fir lut5_bijection_spec.fir lut5_bijection_opt.fir; do \
		if [ -f "../../ELE/build/$$f" ]; then \
			cmp -s "$(BUILD)/$$f" "../../ELE/build/$$f" || \
			(echo "ERROR: Cross-repo SSOT mismatch between ELE and QUA $$f!" >&2 && exit 1); \
			echo "  [PASS] Cross-repo byte-identical SSOT verified: ELE $$f == QUA $$f"; \
		fi; \
	done

examples: $(TARGET) lean-fir
	./$(TARGET) opt $(BUILD)/adder_csa_phase.fir -o $(BUILD)/adder_csa_phase.kmb --clean-anc
	./$(TARGET) opt $(BUILD)/gf12_arith.fir -o $(BUILD)/gf12_arith.kmb --clean-anc
	./$(TARGET) opt $(BUILD)/wide_lut_cap.fir -o $(BUILD)/wide_lut_cap.kmb --clean-anc
	./$(TARGET) verify $(BUILD)/adder_csa_phase.kmb --ref $(BUILD)/adder_csa_phase.fir
	./$(TARGET) verify $(BUILD)/gf12_arith.kmb --ref $(BUILD)/gf12_arith.fir
	./$(TARGET) verify $(BUILD)/wide_lut_cap.kmb --ref $(BUILD)/wide_lut_cap.fir
	./$(TARGET) verify $(BUILD)/wide_lut_cap.kmb

clean:
	rm -rf $(BUILD)
