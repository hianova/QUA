# `QUA` — 全端容錯量子架構（QCC v3 × Ring-0 QEC 控制 × PQC × Magic State Cultivation）

本儲存庫整合四個互鎖的純 C17 / C++20 容錯量子子系統，並支援 **`FLOW v10.0.0` Zero-Lean 雙軌介面**（`../FLOW` via `FLOW_DIR`）：預設透過原生 [`ir/*.fasm`](ir) + Z3 `QF_BV` miter 完成無 Lean 依賴建置與驗證（`LEAN=false`），同時保留 [`Qua.lean`](Qua.lean) 之位元組完全一致檢查。

---

## 目錄結構

- [`qctrl/`](qctrl/DEMO_REPORT.md)：**Ring-0 40 位元超導晶片自主校準、Break-Even 與實時 2D / 3D 表面碼解碼引擎**
  - [`qctrl_ring0.h`](qctrl/qctrl_ring0.h) / [`qctrl_decoder_ring0.h`](qctrl/qctrl_decoder_ring0.h) / [`qctrl_decoder3d_ring0.h`](qctrl/qctrl_decoder3d_ring0.h) / [`qctrl_ring0.c`](qctrl/qctrl_ring0.c)：Freestanding Pure-GPR 4 色並行校準核＋**實時 2D/3D 旋轉表面碼解碼器（`d = 3, 5, 7` Bitmask-DP MWPM ＋ 3D 時空 Safe-Component BFS ＋ 精確迭代整數 Edmonds' Blossom `k <= 48` ＋ EWMA/CUSUM 漂移監控）**，`0` 浮點/SIMD、`0` 除法、`0` heap、`0` 間接跳躍；校準全掃 ~12 µs（`~178 ns/edge`），2D 解碼 `< 8 ns/rnd`、3D 解碼 `3.9–34.0 ns/rnd`（`Λ_3D(3→5) = 4.02×`、`Λ_3D(5→7) = 5.71×`）並通過 2D/3D 全權重 `≤ (d-1)/2` 窮舉 `0` 邏輯誤判及 `500/500` DP 等價驗證。
  - [`qctrl_breakeven.h`](qctrl/qctrl_breakeven.h) / [`qctrl_sim40.c`](qctrl/qctrl_sim40.c)：純 C17 精確 512 維 FWHT `[[9, 1, 3]]` 表面碼驗證、`decode-bench` 與 `decode3d-bench`，輸出 [`qctrl/DEMO_REPORT.md`](qctrl/DEMO_REPORT.md) 與 `build/qctrl_ftqc_profile.json`。
- [`qcc/`](qcc/README.md)：**QCC v3 通用容錯量子電路編譯器**
  - `.kmb` v2（[`qcc_kmb.h`](qcc/qcc_kmb.h)）、OpenQASM 3.0 匯出器（[`qcc_qasm.h`](qcc/qcc_qasm.h)，支援 `qcc export-qasm`）、SSA QIR（[`qir.h`](qcc/qir.h)）、`.fir` QTT 線性閘口（[`qir_fir.h`](qcc/qir_fir.h)）、E-Graph 重寫 `R1`–`R11`（[`qir_rewrite.h`](qcc/qir_rewrite.h)）、`SEARU` 2-Jet 排程與乾淨反計算驗證（[`qir_opt.h`](qcc/qir_opt.h)、[`qir_verify.h`](qcc/qir_verify.h)）。
  - [`qir_qaoa_mb.h`](qcc/qir_qaoa_mb.h)：宣告式 **QAOA MaxCut-12（Base vs Popcount-LUT vs Wallace 3:2）**、**2D 4×4 多體 Ising Trotter**、**Plaquette LGT 四體微擾** 與 **3-SAT Grover Oracle** 基準生成器，配合 `T1`–`T26` 雙種子自測與 `qcc emit-bench`（同時輸出 6 組 `.kmb` + `.qasm`）。
- [`pqc/`](pqc/README.md)：**PQC / NIST 安全基準與全端 FTQC 共設計引擎**
  - [`include/pqc_oracles.h`](pqc/include/pqc_oracles.h) / [`include/pqc_nist_est.h`](pqc/include/pqc_nist_est.h)：AES-128/192/256、SHA-256、Ascon、Keccak-f[400]、Classic McEliece、ML-KEM、MQ Oracle 與純整數 `Q16.16` `MAXDEPTH` 資源估算。
  - [`include/pqc_ftqc_cultivation.h`](pqc/include/pqc_ftqc_cultivation.h)：整合 `qctrl` 實測物理誤差＋實時解碼延遲＋表面碼距與 **Magic State Cultivation vs 15-to-1 Distillation** 時空體積評估（純整數運算，零 `<math.h>`），透過 `pqc ftqc-suite` 對 8 大基準電路產出動態結論完全一致的 [`FTQC_CODESIGN_REPORT.md`](FTQC_CODESIGN_REPORT.md) 與 JSON（Gen-2 解鎖 **5/8** 基準）。
- [`ir/`](ir) & [`Qua.lean`](Qua.lean)：`FLOW v10.0.0` 雙軌規格（8 份 `.fasm` + SEARU synthesis + Z3 `QF_BV` equivalence 與 `Qua.lean` bit-exact 核對）。
- [`ARCHITECTURE_DISCOVERY_REPORT.md`](ARCHITECTURE_DISCOVERY_REPORT.md) & [`FTQC_CODESIGN_REPORT.md`](FTQC_CODESIGN_REPORT.md)：跨層相變分析、經濟斷崖線（$T_{\text{eq}} \lesssim 250\text{–}300$）、與全端資源實測戰報。
- [`kaggle/`](kaggle/README.md)：Google Quantum AI — Verifiable Quantum Circuit Compilation 本地工作區。

---

## 常用指令

```bash
make qcc LEAN=false   # 免 Lean 預設建置：.fasm -> SEARU synth -> Z3 equiv -> qcc T1..T23 + emit-bench
make qctrl            # Ring-0 Pure-GPR 稽核 + 純 C17 512 維 Break-Even 戰報
make pqc ftqc         # PQC KAT/雙種子測試 + 全端 Cultivation 共設計報告 (FTQC_CODESIGN_REPORT.md)
make -C qcc lean-fir  # 驗證 Qua.lean 與 .fasm 位元組完全一致 + 跨 ELE 同源檢驗
make clean            # 清除 qcc / pqc / qctrl / kaggle 編譯產物
```

---

## CI 驗證範圍

- GitHub Actions CI（`ubuntu-latest` + `macos-latest`）同時驗證 **Zero-Lean 預設管線（`LEAN=false`）**：`qcc`（`T1`–`T23` + `.fasm` + Z3 miter）、`qctrl`（Ring-0 稽核 + `selftest` + `breakeven`）、`pqc` + `ftqc` 閉環報告，以及選配 `lean-fir` 與 `Qua.lean` 等價校驗。
