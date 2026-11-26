export PATH := $(HOME)/.elan/bin:/opt/homebrew/opt/llvm/bin:/opt/homebrew/bin:/usr/local/bin:$(PATH)

FLOW_DIR ?= ../../FLOW
CC       ?= clang
CXX      ?= clang++
PYTHON   ?= python3

CFLAGS   := -std=c17 -O3 -Wall -Wextra -Werror -I$(FLOW_DIR)/core
RING0_CF := -std=c17 -ffreestanding -nostdlib -fno-builtin -fno-stack-protector \
            -mgeneral-regs-only -fno-jump-tables -DQCTRL_EXPORT_RING0=1 \
            -Wall -Wextra -Werror -O2 -I$(FLOW_DIR)/core

BUILD    := build
TARGET   := $(BUILD)/qctrl
RING0_O  := $(BUILD)/qctrl_ring0.o

.PHONY: all ring0 check-headers test demo clean

all: $(TARGET) ring0 check-headers

$(BUILD):
	mkdir -p $(BUILD)

$(RING0_O): qctrl_ring0.c qctrl_ring0.h qctrl_decoder_ring0.h qctrl_decoder3d_ring0.h $(FLOW_DIR)/core/flow.h $(FLOW_DIR)/core/flow_opt_params.h | $(BUILD)
	$(CC) $(RING0_CF) -c qctrl_ring0.c -o $@

ring0: $(RING0_O)
	@echo "=== [qctrl Ring-0 Gate] Pure-GPR Freestanding Audit (0 undef, 0 udiv, 0 FPU/SIMD, 0 br/blr) ==="
	@test -z "$$(nm -u $(RING0_O))" || (echo "ERROR: Undefined external symbols in $(RING0_O):" && nm -u $(RING0_O) && exit 1)
	@! objdump -d $(RING0_O) | grep -Eq '[[:space:]](udiv|sdiv|divq?|idivq?|blr|br)[[:space:]]' || \
		(echo "ERROR: Division or indirect branch/call leaked into $(RING0_O)!" && objdump -d $(RING0_O) && exit 1)
	@! objdump -d $(RING0_O) | grep -Eq '([[:space:]]|%)(fadd|fsub|fmul|fdiv|scvtf|ucvtf|fcvtzs|fcvtzu|fcmp|fmov|addsd|subsd|mulsd|divsd|cvtsi2sd|cvttsd2si|xmm[0-9]|ymm[0-9]|zmm[0-9])' || \
		(echo "ERROR: FPU/SIMD instruction leaked into $(RING0_O)!" && objdump -d $(RING0_O) && exit 1)
	@echo "  [PASS] $(RING0_O): zero undefined symbols, zero udiv/sdiv/blr, zero FPU/SIMD"

check-headers:
	$(CC) -std=c17 -Wall -Wextra -Werror -Wno-unused-function -I$(FLOW_DIR)/core -I. -fsyntax-only qctrl_ring0.h qctrl_breakeven.h qctrl_decoder_ring0.h qctrl_decoder3d_ring0.h
	$(CXX) -std=c++20 -Wall -Wextra -Werror -Wno-unused-function -iquote $(FLOW_DIR)/core -iquote . -x c++-header -fsyntax-only qctrl_ring0.h qctrl_breakeven.h qctrl_decoder_ring0.h qctrl_decoder3d_ring0.h
	@echo "  [PASS] qctrl_ring0.h + qctrl_breakeven.h + qctrl_decoder_ring0.h + qctrl_decoder3d_ring0.h clean under C17 + C++20"

$(TARGET): qctrl_sim40.c qctrl_ring0.h qctrl_breakeven.h qctrl_decoder_ring0.h qctrl_decoder3d_ring0.h $(FLOW_DIR)/core/flow.h $(FLOW_DIR)/core/flow_opt_params.h | $(BUILD)
	$(CC) $(CFLAGS) qctrl_sim40.c -lm -o $@

test: $(TARGET) ring0 check-headers
	./$(TARGET) selftest

demo: test
	./$(TARGET) bench40 --out $(BUILD)/bench40.json
	./$(TARGET) breakeven --in $(BUILD)/bench40.json --report DEMO_REPORT.md --ftqc-out $(BUILD)/qctrl_ftqc_profile.json
	./$(TARGET) decode-bench --profile-in $(BUILD)/qctrl_ftqc_profile.json --profile-out $(BUILD)/qctrl_ftqc_profile.json --report DEMO_REPORT.md
	./$(TARGET) decode3d-bench --profile-in $(BUILD)/qctrl_ftqc_profile.json --profile-out $(BUILD)/qctrl_ftqc_profile.json --report DEMO_REPORT.md

clean:
	rm -rf $(BUILD)
