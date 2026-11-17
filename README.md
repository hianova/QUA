# `QUA` — QCC v3 通用量子電路編譯器、PQC / NIST 量子密碼分析資源編譯器與 Kaggle 驗證工作區

本儲存庫整合純 C17（C++20 標頭相容）容錯量子電路編譯器 **QCC v3**、**PQC / NIST 安全基準量子密碼分析資源編譯器**、統一 `FLOW v8.0.0` SSOT 介面（[`Qua.lean`](Qua.lean) + [`qcc/qir_fir.h`](qcc/qir_fir.h)）與 Google Quantum AI 競賽驗證管線，並透過相鄰目錄 `../FLOW`（依 `FLOW_DIR` 指定）引用 `Flow.lean`、`flow.h`、`searu.h`、`flow_opt_params.h` 與 `flow_ir.h`。

---

## 目錄結構

- [`Qua.lean`](Qua.lean)：純參數配置檔（`import Flow`，`0` 本地定理、`0` 自訂序列化），以 `Flow.Morphism.Builder` 宣告並透過 `emitFirChecked` 產出 `.linear` 量子算子圖（`gf256_sbox_mac.fir`、`rev_alu_fxp.fir`、`adder_csa_phase.fir`、`gf12_arith.fir`、`wide_lut_cap.fir`）與 `.cartesian` 拒斥見證（`cartesian_reject.fir`）
- [`qcc/`](qcc/README.md)：**QCC v3 通用容錯量子電路編譯器**
  - `.kmb` v2 雙向編解碼、關鍵路徑 `D_toff` / `D_circ` 統計（[`qcc_kmb.h`](qcc/qcc_kmb.h)）、Arena 與 256-bit / 128-bit 基礎元件（[`qcc_arena.h`](qcc/qcc_arena.h)）
  - SSA QIR（[`qir.h`](qcc/qir.h)）、統一 `.fir` QTT 線性閘口橋接器（[`qir_fir.h`](qcc/qir_fir.h)，拒斥 `.cartesian` 並降階 `.linear` 至 QIR）、無尺度 `>256` qubit 反編譯提升（[`qir_lift.h`](qcc/qir_lift.h)）
  - E-Graph 飽和重寫 `R1`–`R11`（含 $\text{GF}(2^8)$ Rank-27、[`qir_rewrite.h`](qcc/qir_rewrite.h)）、容錯核心庫（[`qir_kernels.h`](qcc/qir_kernels.h)）、實體配置降階（[`qir_codegen.h`](qcc/qir_codegen.h)）、`SEARU` × `FLOW` 2-Jet 排程最佳化與 65,536-qubit 配置（[`qir_opt.h`](qcc/qir_opt.h)）、`.kmb` 通道與統一編譯 API（[`qir_pass_kmb.h`](qcc/qir_pass_kmb.h)）、雙種子等價驗證（[`qir_verify.h`](qcc/qir_verify.h)）及 CLI / `T1`–`T22` 測試套件 + `lean-fir` 跨界同圖雙投影驗證（`ELE` `.sv`/`.gds` vs `QUA` `.kmb`，[`qcc_main.c`](qcc/qcc_main.c)）
- [`pqc/`](pqc/README.md)：**PQC / NIST 安全基準量子密碼分析資源編譯器**
  - [`include/pqc_oracles.h`](pqc/include/pqc_oracles.h)：純宣告式 SSA QIR Oracle 建構器與獨立經典參考實作，涵蓋 **AES-128/192/256**、**SHA-256 / Ascon-p[12] / Keccak-f[400]**、**Classic McEliece ($\text{GF}(2^{12})$ Goppa & ISD)**、**ML-KEM ($q=3329$ LWE)** 與 **MQ ($n=24, m=24$)**（零手動 `.kmb` 閘發射）
  - [`include/pqc_nist_est.h`](pqc/include/pqc_nist_est.h)：純整數 Q16.16 `log2` 與 Zalka / JNRV `MAXDEPTH` ($2^{40} / 2^{64} / 2^{96}$) 容錯資源估算引擎（Gidney `4 * CCX` 與酉 `7 * CCX` 等效 T 閘數、DAG `T-depth` 與 `Circ-depth`）
  - [`src/pqc_main.c`](pqc/src/pqc_main.c)：`selftest`（FIPS-197 / FIPS-180-4 / SP 800-232 KAT + 64-lane 雙種子 + 植入解相位驗證）、`bench`（JSON/Markdown 報告輸出）與 `analyze`
- [`kaggle/`](kaggle/README.md)：**Google Quantum AI — Verifiable Quantum Circuit Compilation**
  - [`src/qua_synth.c`](kaggle/src/qua_synth.c)：純薄層 QCC 呼叫端（無特製閘級啟發式、無直接 `.kmb` 閘發射），僅宣告 QIR 問題規格並由 `qcc_compile_module_to_file` 編譯輸出 8 項基礎原語自測電路、`invsqrt` 與 `dfthcsos` 電路
  - [`dfthcsos_opt/dft_bootstrap.c`](kaggle/dfthcsos_opt/dft_bootstrap.c)：Apple Silicon `Accelerate` 張量分解與 Givens 量化器
  - [`scripts/verify_harness.py`](kaggle/scripts/verify_harness.py)：原語 Fuzz 驗證、Kickmix 雙種子模擬、官方驗證器執行與提交打包腳本

---

## 依賴設定

- 預設透過 `FLOW_DIR ?= ../FLOW`（於子目錄執行時為 `../../FLOW`）引用隔壁目錄的 `flow.h`、`searu.h` 與 `flow_opt_params.h`，若 `FLOW` 位於其他位置可於執行 `make` 時覆寫 `FLOW_DIR=/path/to/FLOW`。

---

## 常用指令

```bash
cd /Users/kuangtalin/Documents/QUA
make qcc              # 檢查 C17/C++20 標頭相容性，執行 qcc T1..T21 自測與 5 個範例驗證
make pqc              # 編譯 pqc 並執行 KAT + 64-lane 雙種子自測與 qcc verify
make -C pqc bench     # 產出全系列 PQC/NIST 基準 .kmb、build/pqc_report.json 與 build/pqc_report.md
make all              # 編譯 qcc、pqc 與 kaggle 執行檔
make test             # 執行 pqc 自測、qcc 聯動驗證 (verify-qcc)、8 項 selftest Fuzz 與 Kickmix pytest
make clean            # 清除編譯產物
```

進一步的 Kaggle 評分與打包指令：

```bash
make -C kaggle stage-init stage-sync score pack
```
