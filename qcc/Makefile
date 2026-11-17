export PATH := $(HOME)/.elan/bin:/opt/homebrew/opt/llvm/bin:/opt/homebrew/bin:/usr/local/bin:$(PATH)

FLOW_DIR ?= ../../FLOW
CC       := cc
CXX      ?= c++
LEAN     ?= lean
CFLAGS   := -std=c17 -O3 -Wall -Wextra -Werror -I$(FLOW_DIR)
BUILD    := build
TARGET   := $(BUILD)/qcc
HDRS     := qcc_arena.h qcc_kmb.h qir.h qir_lift.h qir_rewrite.h \
            qir_kernels.h qir_codegen.h qir_opt.h qir_pass_kmb.h qir_verify.h qir_fir.h

.PHONY: all clean check-headers test lean-fir examples

all: $(TARGET)

$(BUILD):
	mkdir -p $(BUILD)

$(TARGET): qcc_main.c $(HDRS) $(FLOW_DIR)/flow.h $(FLOW_DIR)/searu.h $(FLOW_DIR)/flow_opt_params.h $(FLOW_DIR)/flow_ir.h | $(BUILD)
	$(CC) $(CFLAGS) -o $@ qcc_main.c

check-headers:
	$(CC) -std=c17 -Wall -Wextra -Werror -Wno-unused-function -I$(FLOW_DIR) -I. -fsyntax-only qir_verify.h qir_pass_kmb.h qir_fir.h
	$(CXX) -std=c++20 -Wall -Wextra -Werror -Wno-unused-function -iquote $(FLOW_DIR) -iquote . -x c++-header -fsyntax-only qir_verify.h qir_pass_kmb.h qir_fir.h

test: $(TARGET) check-headers
	./$(TARGET) selftest

lean-fir: $(TARGET)
	@if grep -Eq '^[[:space:]]*(theorem|lemma|axiom)[[:space:]]' ../Qua.lean; then \
		echo "ERROR: Qua.lean violates SSOT pure-config rule (found theorem)!" >&2; exit 1; \
	fi
	@qua_lines="$$(wc -l < ../Qua.lean | tr -d '[:space:]')"; \
	if [ "$$qua_lines" -gt 120 ]; then echo "ERROR: Qua.lean exceeds 120 lines ($$qua_lines)" >&2; exit 1; fi
	$(MAKE) -C "$(FLOW_DIR)" build/lean/Flow.olean
	LEAN_PATH="$(FLOW_DIR)/build/lean" $(LEAN) -DwarningAsError=true --run ../Qua.lean $(BUILD)
	./$(TARGET) opt $(BUILD)/gf256_sbox_mac.fir -o $(BUILD)/gf256_sbox_mac.kmb --clean-anc
	./$(TARGET) verify $(BUILD)/gf256_sbox_mac.kmb --ref $(BUILD)/gf256_sbox_mac.fir --states 4096
	./$(TARGET) opt $(BUILD)/rev_alu_fxp.fir -o $(BUILD)/rev_alu_fxp.kmb --clean-anc
	./$(TARGET) verify $(BUILD)/rev_alu_fxp.kmb --ref $(BUILD)/rev_alu_fxp.fir --states 4096
	@if ./$(TARGET) opt $(BUILD)/cartesian_reject.fir -o $(BUILD)/bad_cart.kmb 2>/dev/null; then \
		echo "ERROR: cartesian_reject.fir was NOT rejected by quantum compiler!" >&2; exit 1; \
	fi
	@if [ -f ../../ELE/build/gf256_sbox_mac.fir ]; then \
		cmp -s $(BUILD)/gf256_sbox_mac.fir ../../ELE/build/gf256_sbox_mac.fir || \
		(echo "ERROR: Cross-repo SSOT mismatch between ELE and QUA gf256_sbox_mac.fir!" >&2 && exit 1); \
		echo "  [PASS] Cross-repo byte-identical SSOT verified: ELE gf256_sbox_mac.fir == QUA gf256_sbox_mac.fir"; \
	fi

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
