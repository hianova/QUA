export PATH := $(HOME)/.elan/bin:/opt/homebrew/opt/llvm/bin:/opt/homebrew/bin:/usr/local/bin:$(PATH)

FLOW_DIR ?= ../../FLOW
CC       := cc
CXX      ?= c++
LEAN_BIN ?= lean
LEAN     ?= false
Z3       ?= $(shell command -v z3 2>/dev/null || true)
FLOW_CLI := $(FLOW_DIR)/build/flow
CFLAGS   := -std=c17 -O3 -Wall -Wextra -Werror -I$(FLOW_DIR)/core -I$(FLOW_DIR)/sdk
BUILD    := build
TARGET   := $(BUILD)/qcc
HDRS     := qcc_arena.h qcc_kmb.h qcc_qasm.h qir.h qir_lift.h qir_rewrite.h \
            qir_kernels.h qir_codegen.h qir_opt.h qir_pass_kmb.h qir_verify.h qir_fir.h qir_qaoa_mb.h

IR_MODULES := gf256_sbox_mac rev_alu_fxp cartesian_reject adder_csa_phase \
              gf12_arith wide_lut_cap gf256_tower_spec lut5_bijection_spec

.PHONY: all clean check-headers test fir-artifacts lean-fir examples

all: $(TARGET)

$(BUILD):
	mkdir -p $(BUILD)

$(TARGET): qcc_main.c $(HDRS) $(FLOW_DIR)/core/flow.h $(FLOW_DIR)/core/searu.h $(FLOW_DIR)/core/flow_opt_params.h $(FLOW_DIR)/sdk/flow_ir_wire.h $(FLOW_DIR)/sdk/flow_ir.h | $(BUILD)
	$(CC) $(CFLAGS) -o $@ qcc_main.c

check-headers:
	$(CC) -std=c17 -Wall -Wextra -Werror -Wno-unused-function -I$(FLOW_DIR)/core -I$(FLOW_DIR)/sdk -I. -fsyntax-only qir_verify.h qir_pass_kmb.h qir_fir.h qir_qaoa_mb.h qcc_qasm.h
	$(CXX) -std=c++20 -Wall -Wextra -Werror -Wno-unused-function -iquote $(FLOW_DIR)/core -iquote $(FLOW_DIR)/sdk -iquote . -x c++-header -fsyntax-only qir_verify.h qir_pass_kmb.h qir_fir.h qir_qaoa_mb.h qcc_qasm.h

test: $(TARGET) check-headers
	./$(TARGET) selftest

fir-artifacts: $(TARGET) | $(BUILD)
ifeq ($(LEAN),true)
	@$(MAKE) lean-fir FLOW_DIR="$(FLOW_DIR)" LEAN_BIN="$(LEAN_BIN)"
else
	$(MAKE) -C "$(FLOW_DIR)" build/flow
	@for m in $(IR_MODULES); do \
		"$(FLOW_CLI)" fir asm "../ir/$$m.fasm" "$(BUILD)/$$m.fir"; \
	done
	"$(FLOW_CLI)" fir-synth $(BUILD)/gf256_tower_spec.fir $(BUILD)/gf256_tower_a32.fir --name gf256_tower_a32
	"$(FLOW_CLI)" fir-synth $(BUILD)/gf256_tower_spec.fir $(BUILD)/gf256_tower_d4.fir --name gf256_tower_d4 --max-depth 2
	"$(FLOW_CLI)" fir-synth $(BUILD)/lut5_bijection_spec.fir $(BUILD)/lut5_bijection_opt.fir --name lut5_bijection_opt
	@if [ -n "$(Z3)" ]; then \
		for pair in "gf256_tower_spec gf256_tower_a32" "gf256_tower_spec gf256_tower_d4" "lut5_bijection_spec lut5_bijection_opt"; do \
			set -- $$pair; \
			"$(FLOW_CLI)" fir equiv "$(BUILD)/$$1.fir" "$(BUILD)/$$2.fir" "$(BUILD)/equiv_$$2.smt2"; \
			res=$$($(Z3) "$(BUILD)/equiv_$$2.smt2"); \
			if [ "$$res" != "unsat" ]; then echo "ERROR: Z3 QF_BV equiv failed on $$2: $$res" >&2; exit 1; fi; \
			echo "  [Z3 PASS] $$2 == $$1 (unsat)"; \
		done; \
	fi
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
endif

lean-fir: $(TARGET)
	@if grep -Eq '^[[:space:]]*(theorem|lemma|axiom)[[:space:]]' ../Qua.lean; then \
		echo "ERROR: Qua.lean violates SSOT pure-config rule (found theorem)!" >&2; exit 1; \
	fi
	@qua_lines="$$(wc -l < ../Qua.lean | tr -d '[:space:]')"; \
	if [ "$$qua_lines" -gt 120 ]; then echo "ERROR: Qua.lean exceeds 120 lines ($$qua_lines)" >&2; exit 1; fi
	$(MAKE) -C "$(FLOW_DIR)" build/lean/Flow.olean build/flow
	LEAN_PATH="$(FLOW_DIR)/build/lean" $(LEAN_BIN) -DwarningAsError=true --run ../Qua.lean $(BUILD)
	@mkdir -p $(BUILD)/fasm_parity
	@for m in gf256_sbox_mac rev_alu_fxp adder_csa_phase gf12_arith wide_lut_cap gf256_tower_spec lut5_bijection_spec; do \
		"$(FLOW_CLI)" fir asm "../ir/$$m.fasm" "$(BUILD)/fasm_parity/$$m.fir"; \
		cmp -s "$(BUILD)/$$m.fir" "$(BUILD)/fasm_parity/$$m.fir" || \
			{ echo "ERROR: .fasm <-> Qua.lean byte mismatch on $$m.fir!" >&2; exit 1; }; \
	done
	@echo "  [PASS] .fasm <-> Qua.lean 7/7 linear .fir byte-identical parity"
	"$(FLOW_CLI)" fir-synth $(BUILD)/gf256_tower_spec.fir $(BUILD)/gf256_tower_a32.fir --name gf256_tower_a32
	"$(FLOW_CLI)" fir-synth $(BUILD)/gf256_tower_spec.fir $(BUILD)/gf256_tower_d4.fir --name gf256_tower_d4 --max-depth 2
	"$(FLOW_CLI)" fir-synth $(BUILD)/lut5_bijection_spec.fir $(BUILD)/lut5_bijection_opt.fir --name lut5_bijection_opt
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
			{ echo "ERROR: Cross-repo SSOT mismatch between ELE and QUA $$f!" >&2; exit 1; }; \
			echo "  [PASS] Cross-repo byte-identical SSOT verified: ELE $$f == QUA $$f"; \
		fi; \
	done

examples: $(TARGET) fir-artifacts
	./$(TARGET) opt $(BUILD)/adder_csa_phase.fir -o $(BUILD)/adder_csa_phase.kmb --clean-anc
	./$(TARGET) opt $(BUILD)/gf12_arith.fir -o $(BUILD)/gf12_arith.kmb --clean-anc
	./$(TARGET) opt $(BUILD)/wide_lut_cap.fir -o $(BUILD)/wide_lut_cap.kmb --clean-anc
	./$(TARGET) verify $(BUILD)/adder_csa_phase.kmb --ref $(BUILD)/adder_csa_phase.fir
	./$(TARGET) verify $(BUILD)/gf12_arith.kmb --ref $(BUILD)/gf12_arith.fir
	./$(TARGET) verify $(BUILD)/wide_lut_cap.kmb --ref $(BUILD)/wide_lut_cap.fir
	./$(TARGET) verify $(BUILD)/wide_lut_cap.kmb
	./$(TARGET) emit-bench $(BUILD)

clean:
	rm -rf $(BUILD)
