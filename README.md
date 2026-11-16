# `QUA` — QCC v3 通用量子電路編譯器與 Kaggle Verifiable Quantum Circuit Compilation 工作區

本儲存庫整合純 C17（C++20 標頭相容）容錯量子電路編譯器 **QCC v3** 與 Google Quantum AI 競賽驗證管線，並透過相鄰目錄 `../FLOW`（不複製標頭檔，依 `FLOW_DIR` 指定）引用 `flow.h`、`searu.h` 與 `flow_opt_params.h`。

---

## 目錄結構

- [`qcc/`](qcc/README.md)：**QCC v3 通用容錯量子電路編譯器**
  - `.kmb` v2 雙向編解碼（[`qcc_kmb.h`](qcc/qcc_kmb.h)）、Arena 與 256-bit / 128-bit 基礎元件（[`qcc_arena.h`](qcc/qcc_arena.h)）
  - SSA QIR（[`qir.h`](qcc/qir.h) / [`qir_dsl.h`](qcc/qir_dsl.h)）、無尺度 `>256` qubit 反編譯提升（[`qir_lift.h`](qcc/qir_lift.h)）
  - E-Graph 飽和重寫 `R1`–`R11`（[`qir_rewrite.h`](qcc/qir_rewrite.h)）、容錯核心庫（[`qir_kernels.h`](qcc/qir_kernels.h)）、實體配置降階（[`qir_codegen.h`](qcc/qir_codegen.h)）、`SEARU` × `FLOW` 2-Jet 排程最佳化與零垃圾級聯反計算（[`qir_opt.h`](qcc/qir_opt.h)）、`.kmb` 通道與統一編譯 API（[`qir_pass_kmb.h`](qcc/qir_pass_kmb.h)）、雙種子等價驗證（[`qir_verify.h`](qcc/qir_verify.h)）及 CLI / `T1`–`T20` 測試套件（[`qcc_main.c`](qcc/qcc_main.c)）
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
make qcc     # 檢查 C17/C++20 標頭相容性，執行 qcc T1..T20 自測與 5 個範例驗證
make all     # 編譯 qcc、kaggle/build/qua_synth 與 kaggle/build/dft_bootstrap
make test    # 執行 qcc 聯動驗證 (verify-qcc)、8 項 selftest Fuzz 與 Kickmix pytest
make clean   # 清除編譯產物
```

進一步的 Kaggle 評分與打包指令：

```bash
make -C kaggle stage-init stage-sync score pack
```
