FLOW_DIR ?= ../../FLOW
CC       := cc
CFLAGS   := -std=c17 -O3 -Wall -Wextra -Werror -I$(FLOW_DIR)
BUILD    := build
TARGET   := $(BUILD)/qcc
HDRS     := qcc_arena.h qcc_kmb.h qir.h qir_dsl.h qir_lift.h qir_rewrite.h \
            qir_kernels.h qir_codegen.h qir_opt.h qir_pass_kmb.h qir_verify.h

.PHONY: all clean test examples

all: $(TARGET)

$(BUILD):
	mkdir -p $(BUILD)

$(TARGET): qcc_main.c $(HDRS) $(FLOW_DIR)/flow.h $(FLOW_DIR)/searu.h $(FLOW_DIR)/flow_opt_params.h | $(BUILD)
	$(CC) $(CFLAGS) -o $@ qcc_main.c

test: $(TARGET)
	./$(TARGET) selftest

examples: $(TARGET)
	./$(TARGET) opt examples/adder_csa_phase.qir -o $(BUILD)/adder_csa_phase.kmb
	./$(TARGET) opt examples/invsqrt_poly.qir -o $(BUILD)/invsqrt_poly.kmb --max-q 64
	./$(TARGET) opt examples/gf12_arith.qir -o $(BUILD)/gf12_arith.kmb
	./$(TARGET) opt examples/wide_lut_cap.qir -o $(BUILD)/wide_lut_cap.kmb --max-q 40
	./$(TARGET) opt examples/invsqrt_phase_oracle.qir -o $(BUILD)/invsqrt_phase_oracle.kmb
	./$(TARGET) verify $(BUILD)/adder_csa_phase.kmb --ref examples/adder_csa_phase.qir
	./$(TARGET) verify $(BUILD)/gf12_arith.kmb --ref examples/gf12_arith.qir
	./$(TARGET) verify $(BUILD)/invsqrt_phase_oracle.kmb --ref examples/invsqrt_phase_oracle.qir

clean:
	rm -rf $(BUILD)
