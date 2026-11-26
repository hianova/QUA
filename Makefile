FLOW_DIR ?= ../FLOW
LEAN     ?= false

.PHONY: all qcc pqc qctrl kaggle ftqc test clean

all: qcc pqc qctrl
	$(MAKE) -C kaggle FLOW_DIR=$(abspath $(FLOW_DIR)) build-synth build-bootstrap

qcc:
	$(MAKE) -C qcc FLOW_DIR=$(abspath $(FLOW_DIR)) LEAN=$(LEAN) all test examples

pqc:
	$(MAKE) -C pqc FLOW_DIR=$(abspath $(FLOW_DIR)) all test

qctrl:
	$(MAKE) -C qctrl FLOW_DIR=$(abspath $(FLOW_DIR)) test demo

kaggle:
	$(MAKE) -C kaggle FLOW_DIR=$(abspath $(FLOW_DIR)) all

ftqc: qcc pqc qctrl
	$(MAKE) -C pqc FLOW_DIR=$(abspath $(FLOW_DIR)) ftqc

test:
	$(MAKE) -C qcc FLOW_DIR=$(abspath $(FLOW_DIR)) LEAN=$(LEAN) test
	$(MAKE) -C qctrl FLOW_DIR=$(abspath $(FLOW_DIR)) test demo
	$(MAKE) -C pqc FLOW_DIR=$(abspath $(FLOW_DIR)) test ftqc
	$(MAKE) -C kaggle FLOW_DIR=$(abspath $(FLOW_DIR)) verify-qcc test

clean:
	$(MAKE) -C qcc clean
	$(MAKE) -C pqc clean
	$(MAKE) -C qctrl clean
	$(MAKE) -C kaggle clean
