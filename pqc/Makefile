FLOW_DIR ?= ../../FLOW
QCC_DIR  ?= ../qcc

CC       := cc
CXX      ?= c++
CFLAGS   := -std=c17 -O3 -Wall -Wextra -Werror -I$(FLOW_DIR) -I$(QCC_DIR) -Iinclude
BUILD    := build
TARGET   := $(BUILD)/pqc
HDRS     := include/pqc_nist_est.h include/pqc_oracles.h

.PHONY: all check-headers test bench clean

all: $(TARGET)

$(BUILD):
	mkdir -p $(BUILD)

$(TARGET): src/pqc_main.c $(HDRS) $(wildcard $(QCC_DIR)/*.h) $(FLOW_DIR)/flow.h $(FLOW_DIR)/searu.h $(FLOW_DIR)/flow_opt_params.h | $(BUILD)
	$(CC) $(CFLAGS) src/pqc_main.c -o $@

check-headers:
	$(CC) -std=c17 -Wall -Wextra -Werror -Wno-unused-function -I$(FLOW_DIR) -I$(QCC_DIR) -Iinclude -fsyntax-only $(HDRS)
	$(CXX) -std=c++20 -Wall -Wextra -Werror -Wno-unused-function -I$(FLOW_DIR) -I$(QCC_DIR) -Iinclude -x c++-header -fsyntax-only $(HDRS)

test: $(TARGET) check-headers
	./$(TARGET) selftest --emit-sample $(BUILD)/selftest_aes128_fwd.kmb
	test -x $(QCC_DIR)/build/qcc || $(MAKE) -C $(QCC_DIR) FLOW_DIR=$(abspath $(FLOW_DIR)) all
	$(QCC_DIR)/build/qcc verify $(BUILD)/selftest_aes128_fwd.kmb --states 256

bench: $(TARGET)
	./$(TARGET) bench $(BUILD) --json $(BUILD)/pqc_report.json --md $(BUILD)/pqc_report.md

clean:
	rm -rf $(BUILD)
