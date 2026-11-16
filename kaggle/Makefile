# Google Quantum AI — Verifiable Quantum Circuit Compilation (./kaggle)
# Native Apple Silicon Toolchain: C17 (FLOW v7.0.1 + SEARU) + C++20 Kickmix + Python 3.14

FLOW_DIR ?= ../../FLOW
QCC_DIR  ?= ../qcc

CC       := cc
CXX      := c++
CFLAGS   := -std=c17 -O3 -Wall -Wextra -Werror -mcpu=native -I.. -I$(FLOW_DIR) -Isrc
CXXFLAGS := -std=c++20 -O3 -fno-strict-aliasing -mcpu=native
PYTHON   := ./.venv/bin/python
PYTEST   := ./.venv/bin/pytest
KAGGLE   := ./.venv/bin/kaggle
COMP     := google-verifiable-quantum-circuit-compilation

BUILD_DIR     := build
DATA_DIR      := data
COMP_ROOT     := $(DATA_DIR)/$(COMP)
BASELINE_DIR  := $(COMP_ROOT)/sample_submission
VERIFIERS_DIR := $(COMP_ROOT)/verifiers
SUB_DIR       := submissions
STAGE_DIR     := $(SUB_DIR)/stage

.PHONY: all build-synth build-bootstrap build-kickmix test verify-qcc data inspect score-baseline stage-init stage-sync score score-full pack clean

all: build-synth build-bootstrap test

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

build-synth: $(BUILD_DIR) src/qua_synth.c src/qua_kmb.h src/qua_ele.h src/qua_bitop.h src/qua_invsqrt.h src/qua_dqi.h src/qua_dfthcsos.h $(QCC_DIR)/qcc_kmb.h $(QCC_DIR)/qcc_arena.h $(FLOW_DIR)/flow.h $(FLOW_DIR)/searu.h $(FLOW_DIR)/flow_opt_params.h
	$(CC) $(CFLAGS) src/qua_synth.c -o $(BUILD_DIR)/qua_synth

build-bootstrap: $(BUILD_DIR) dfthcsos_opt/dft_bootstrap.c $(FLOW_DIR)/flow.h $(FLOW_DIR)/searu.h $(FLOW_DIR)/flow_opt_params.h
	$(CC) $(CFLAGS) -DACCELERATE_NEW_LAPACK dfthcsos_opt/dft_bootstrap.c -framework Accelerate -o $(BUILD_DIR)/dft_bootstrap

build-kickmix: $(BUILD_DIR)
	cd vendor/Qualtran/kickmix && \
	SRCS=$$(find src -name "*.cc" ! -name "*.test.cc" ! -name "*.perf.cc" ! -name "*.pybind.cc") && \
	$(CXX) $(CXXFLAGS) -Isrc $$SRCS -o ../../../$(BUILD_DIR)/kickmix
	./.venv/bin/pip install --no-build-isolation --quiet ./vendor/Qualtran/kickmix

test: build-synth
	mkdir -p $(BUILD_DIR)/dfthcsos
	test -f $(BUILD_DIR)/dfthcsos/words.bin || cp dfthcsos_opt/words.bin $(BUILD_DIR)/dfthcsos/words.bin
	test -f $(BUILD_DIR)/dfthcsos/submission.npz || cp dfthcsos_opt/submission.npz $(BUILD_DIR)/dfthcsos/submission.npz
	./$(BUILD_DIR)/qua_synth $(BUILD_DIR)
	./$(BUILD_DIR)/kickmix count_operations --groups --in $(BUILD_DIR)/selftest_gf2_12_mul.kmb
	$(PYTHON) scripts/verify_harness.py --selftest $(BUILD_DIR)
	$(PYTEST) vendor/Qualtran/kickmix/src -k "not test_builder_mark_errors" -q

verify-qcc: build-synth
	$(MAKE) -C $(QCC_DIR) FLOW_DIR=$(abspath $(FLOW_DIR)) all test examples
	$(PYTHON) scripts/verify_harness.py --kmb $(QCC_DIR)/build/adder_csa_phase.kmb $(QCC_DIR)/build/invsqrt_poly.kmb $(QCC_DIR)/build/gf12_arith.kmb $(QCC_DIR)/build/wide_lut_cap.kmb $(QCC_DIR)/build/invsqrt_phase_oracle.kmb
	$(QCC_DIR)/build/qcc verify $(BUILD_DIR)/selftest_adder16.kmb
	$(QCC_DIR)/build/qcc verify $(BUILD_DIR)/selftest_gf2_12_mul.kmb
	$(QCC_DIR)/build/qcc verify $(BUILD_DIR)/selftest_linf2.kmb
	$(QCC_DIR)/build/qcc verify $(BUILD_DIR)/selftest_zpow128.kmb

data:
	mkdir -p $(DATA_DIR)
	$(KAGGLE) competitions download -c $(COMP) -p $(DATA_DIR)
	unzip -o $(DATA_DIR)/$(COMP).zip -d $(COMP_ROOT)
	$(PYTHON) scripts/verify_harness.py --inspect-data $(DATA_DIR)

inspect:
	$(PYTHON) scripts/verify_harness.py --inspect-data $(DATA_DIR)

score-baseline:
	$(PYTHON) scripts/verify_harness.py --score $(BASELINE_DIR) $(VERIFIERS_DIR)

stage-init:
	mkdir -p $(STAGE_DIR)/dfthcsos
	ln -sf $(abspath $(BASELINE_DIR)/invsqrt.kmb) $(STAGE_DIR)/invsqrt.kmb
	ln -sf $(abspath $(BASELINE_DIR)/dqi.kmb) $(STAGE_DIR)/dqi.kmb
	for f in $(abspath $(BASELINE_DIR)/dfthcsos)/*; do ln -sf "$$f" $(STAGE_DIR)/dfthcsos/$$(basename "$$f"); done

stage-sync:
	mkdir -p $(STAGE_DIR)/dfthcsos
	rm -f $(STAGE_DIR)/dfthcsos/angle_k_*.kmb $(STAGE_DIR)/dfthcsos/angle_k_*.kmx
	ln -sf $(abspath $(BUILD_DIR)/invsqrt.kmb) $(STAGE_DIR)/invsqrt.kmb
	ln -sf $(abspath $(BUILD_DIR)/dqi.kmb) $(STAGE_DIR)/dqi.kmb
	ln -sf $(abspath $(BUILD_DIR)/dfthcsos/submission.npz) $(STAGE_DIR)/dfthcsos/submission.npz
	for f in $(abspath $(BUILD_DIR)/dfthcsos)/angle_k_*.kmb; do ln -sf "$$f" $(STAGE_DIR)/dfthcsos/$$(basename "$$f"); done

score:
	$(PYTHON) scripts/verify_harness.py --score $(STAGE_DIR) $(VERIFIERS_DIR)

score-full:
	$(PYTHON) scripts/verify_harness.py --score $(STAGE_DIR) $(VERIFIERS_DIR) --full

pack:
	$(PYTHON) scripts/verify_harness.py --pack $(STAGE_DIR) $(SUB_DIR)/submission.zip

clean:
	rm -rf $(BUILD_DIR)/*.o $(BUILD_DIR)/qua_synth $(BUILD_DIR)/selftest_*.kmb
