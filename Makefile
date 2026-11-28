export PATH := $(HOME)/.elan/bin:$(PATH)

FLOW_DIR    ?= ../FLOW
LEAN        ?= false
LEAN_BIN    ?= lean
LEANCHECKER ?= leanchecker

.PHONY: all qcc pqc qctrl kaggle ftqc test lean-qecc clean

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

lean-qecc:
	@python3 -c 'import re,sys; s=open("Qua/QeccCliffordT.lean",encoding="utf-8").read(); s=re.sub(r"/-.*?-/", "", s, flags=re.S); s=re.sub(r"--.*", "", s); m=re.findall(r"\b(?:sorry|axiom|native_decide|set_option)\b", s); sys.exit("ERROR: forbidden token(s) in Qua/QeccCliffordT.lean: " + str(m)) if m else None'
	@mkdir -p build/lean/Qua
	$(LEAN_BIN) -DwarningAsError=true -R . -o build/lean/Qua/QeccCliffordT.olean Qua/QeccCliffordT.lean
	LEAN_PATH=build/lean $(LEANCHECKER) Qua.QeccCliffordT
	@python3 -c 'import re,subprocess,os; src=open("Qua/QeccCliffordT.lean",encoding="utf-8").read(); thms=re.findall(r"^(?:theorem|lemma)\s+([A-Za-z0-9_]+)", src, re.M); assert thms, "no public theorems found"; chk="build/lean/_check_qecc_axioms.lean"; open(chk,"w").write("import Qua.QeccCliffordT\n" + "\n".join(f"#print axioms Qua.QeccCliffordT.{t}" for t in thms) + "\n"); env=dict(os.environ, LEAN_PATH="build/lean"); r=subprocess.run(["$(LEAN_BIN)", "-DwarningAsError=true", chk], env=env, capture_output=True, text=True); os.remove(chk); print(r.stdout.strip()); assert r.returncode==0, r.stderr; assert "sorryAx" not in r.stdout and "sorryAx" not in r.stderr, "sorryAx detected"; print(f"  [PASS] Qua.QeccCliffordT: {len(thms)} public theorems verified (leanchecker + #print axioms 0 sorryAx)")'

clean:
	$(MAKE) -C qcc clean
	$(MAKE) -C pqc clean
	$(MAKE) -C qctrl clean
	$(MAKE) -C kaggle clean
	rm -rf build
