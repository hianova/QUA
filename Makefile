FLOW_DIR ?= ../FLOW

.PHONY: all qcc pqc qctrl kaggle test clean

all: qcc pqc
	$(MAKE) -C kaggle FLOW_DIR=$(abspath $(FLOW_DIR)) build-synth build-bootstrap

qcc:
	$(MAKE) -C qcc FLOW_DIR=$(abspath $(FLOW_DIR)) all test examples

pqc:
	$(MAKE) -C pqc FLOW_DIR=$(abspath $(FLOW_DIR)) all test

qctrl:
	$(MAKE) -C qctrl FLOW_DIR=$(abspath $(FLOW_DIR)) test demo

kaggle:
	$(MAKE) -C kaggle FLOW_DIR=$(abspath $(FLOW_DIR)) all

test:
	$(MAKE) -C pqc FLOW_DIR=$(abspath $(FLOW_DIR)) test
	$(MAKE) -C kaggle FLOW_DIR=$(abspath $(FLOW_DIR)) verify-qcc test

clean:
	$(MAKE) -C qcc clean
	$(MAKE) -C pqc clean
	$(MAKE) -C kaggle clean
