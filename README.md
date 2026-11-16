# `qua` — QCC 通用量子電路編譯器與 Kaggle Verifiable Quantum Circuit Compilation 工作區

本儲存庫由原 `FLOW` 專案抽離獨立而成，整合了純 C17 容錯量子電路編譯器與 Google Quantum AI 競賽驗證管線，並透過相鄰的 [`../FLOW`](file:///Users/kuangtalin/Documents/FLOW)（單一來源，不複製標頭檔）提供 `flow.h`、`searu.h` 與 `flow_opt_params.h`。

## 目錄結構

- [`qcc/`](file:///Users/kuangtalin/Documents/qua/qcc/Makefile)：**QCC v2 通用量子電路編譯器**
  - `.kmb` v2 雙向編解碼 (`qcc_kmb.h`)、Arena 記憶體管理 (`qcc_arena.h`)
  - SSA QIR (`qir.h` / `qir_dsl.h`)、無尺度 >256 Qubit `.kmb` 提升 (`qir_lift.h`)
  - E-Graph 飽和重寫 (`qir_rewrite.h`)、QROAM / LinF2 / GF(2^12) / MulConstAdd / TriSqr 核心 (`qir_kernels.h`)、Symbiotic SEARU × FLOW 2-Jet 最佳化與零垃圾串聯反計算 (`qir_opt.h`)、`qir_pass_kmb.h` 與雙種子驗證器 (`qir_verify.h`)
- [`kaggle/`](file:///Users/kuangtalin/Documents/qua/kaggle/README.md)：**Google Quantum AI — Verifiable Quantum Circuit Compilation**
  - `src/`：三大 Benchmark（`invsqrt`、`dqi`、`dfthcsos`）純 C17 合成器（`qua_synth.c`）
  - `dfthcsos_opt/`：`Accelerate` 張量分解與 Givens 量化最佳化器（`dft_bootstrap.c`）
  - `scripts/verify_harness.py` & `vendor/Qualtran/kickmix`：Fuzz 驗證、Kickmix 雙種子模擬與官方評分打包工具

## 依賴設定

- 預設透過 `FLOW_DIR ?= ../FLOW` 引用隔壁目錄的 [`flow.h`](file:///Users/kuangtalin/Documents/FLOW/flow.h)、[`searu.h`](file:///Users/kuangtalin/Documents/FLOW/searu.h) 與 [`flow_opt_params.h`](file:///Users/kuangtalin/Documents/FLOW/flow_opt_params.h)（若 `FLOW` 位於其他路徑，可於指令列覆寫 `make FLOW_DIR=/path/to/FLOW ...`）。

## 常用指令

```bash
cd /Users/kuangtalin/Documents/qua
make qcc              # 建置並執行 qcc 全部 15 項自測與 5 個範例驗證
make all              # 建置 qcc、qua_synth 與 dft_bootstrap
make test             # 執行 kaggle 與 qcc 聯動驗證 (verify-qcc) 及 3,402 項 Kickmix pytest
make -C kaggle score  # 執行三大官方 Verifier 快速模式並計算幾何平均壓縮倍率
make -C kaggle pack   # 驗證所有競賽合約並產出 kaggle/submissions/submission.zip
```
