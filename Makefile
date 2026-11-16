FLOW_DIR ?= ../FLOW

.PHONY: all qcc kaggle test clean

all: qcc
	$(MAKE) -C kaggle FLOW_DIR=$(abspath $(FLOW_DIR)) build-synth build-bootstrap

qcc:
	$(MAKE) -C qcc FLOW_DIR=$(abspath $(FLOW_DIR)) all test examples

kaggle:
	$(MAKE) -C kaggle FLOW_DIR=$(abspath $(FLOW_DIR)) all

test:
	$(MAKE) -C kaggle FLOW_DIR=$(abspath $(FLOW_DIR)) verify-qcc test

clean:
	$(MAKE) -C qcc clean
	$(MAKE) -C kaggle clean
