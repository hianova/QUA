# `./kaggle` — Google Quantum AI: Verifiable Quantum Circuit Compilation

本目錄整合 [Google Verifiable Quantum Circuit Compilation](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/data) 競賽的驗證與張量前處理環境。電路編譯已全面收斂為**無特製閘級啟發式的薄層呼叫器** [`src/qua_synth.c`](src/qua_synth.c)，直接呼叫相鄰通用量子編譯器 [`../qcc/`](../qcc/README.md)（`qcc_compile_module_to_file`）進行合成與驗證，搭配 Apple Silicon `Accelerate` 張量量化器 [`dfthcsos_opt/dft_bootstrap.c`](dfthcsos_opt/dft_bootstrap.c) 與 [`scripts/verify_harness.py`](scripts/verify_harness.py)。

---

## 1. 架構與職責邊界

- [`src/qua_synth.c`](src/qua_synth.c)：純 C17 薄層驅動程式（編譯選項 `-std=c17 -O3 -Wall -Wextra -Werror`）。不含有任何 `qcc_kmb_*` 閘發射呼叫（僅使用 `qcc_kmb_stats_t` 列印統計），僅以宣告式 `qir_*` 描述問題並委託 QCC v3 編譯為 `.kmb` v2 電路，包含：
  - `selftest`：8 項基礎原語自測電路
  - `invsqrt`：呼叫 `qir_emit_dyadic_phase_fn` 建立 22-qubit $1/\sqrt{x}$ 相位預言機
  - `dfthcsos`：讀取張量量化表 [`dfthcsos_opt/words.bin`](dfthcsos_opt/words.bin)，組建批次（`J=5`, `K=15`）`LUT_ROM + CONCAT + GF2M_ADD + SLICE` QIR 並於 `max_qubits=128` 限制下編譯出 `angle_k_0..14.kmb`
- `dqi` 目前不於 [`src/qua_synth.c`](src/qua_synth.c) 合成；使用 `make stage-init && make stage-sync` 打包完整三題提交檔時，若 `build/dqi.kmb` 不存在會保留官方基準 `sample_submission/dqi.kmb` 以通過整體打包合約。
- [`dfthcsos_opt/dft_bootstrap.c`](dfthcsos_opt/dft_bootstrap.c)：以 `ACCELERATE_NEW_LAPACK` 進行 DF-THC 張量分解與 Givens 座標量化，產出 [`dfthcsos_opt/words.bin`](dfthcsos_opt/words.bin) 與 [`dfthcsos_opt/submission.npz`](dfthcsos_opt/submission.npz)。
- [`scripts/verify_harness.py`](scripts/verify_harness.py)：執行 8 項原語 Kickmix 雙種子模擬驗證、通用 `.kmb` 交叉檢驗、官方評分與 `submissions/submission.zip` 打包。

---

## 2. QCC v3 實測結果

### 官方驗證器（`invsqrt` 與 `dfthcsos`）

| Benchmark | 官方基準成本 | QCC v3 實測成本 | 壓縮倍率 | QCC v3 實測規格 |
| :--- | ---: | ---: | ---: | :--- |
| **`invsqrt`** | `2,076,774.23` | **`1,239,630.15`** | **`1.6753×`** | `PASS`（`3,133,452` shots 全域覆蓋），`Q=255`，`CCX=4,630`，`T=1`，`rot=34`，`C_eq=4,861.29`，`D_react=1,370`，`worst_error_over_tolerance=0.4667` |
| **`dfthcsos`** | `3,355,599,588.32` | **`2,039,968,295.19`** | **`1.6449×`** | `PASS`（`12,784` shots，`gap_in_domain=true`），15 張批次電路（`J=5`, `brot=10`, `rows=3,196`），`angles_qubits=121 <= 128`，`angles_toffoli_per_pass=26,093`，總 `T=60,916`，`Q=680` |

### 8 項基礎原語自測電路（`qua_synth selftest` + `verify_harness.py --selftest`）

| 電路 | Qubits | 實測關鍵指標 | Kickmix Fuzz 驗證 |
| :--- | ---: | :--- | :--- |
| `selftest_adder16.kmb` | `47` | `CCX=15, HMR=15` | `4,096` shots `PASS` |
| `selftest_gf2_12_mul.kmb` | `37` | `CCX=42, HMR=42`（Chudnovsky Rank-42） | `4,096` shots `PASS` |
| `selftest_gf64_mul.kmb` | `19` | `CCX=15, HMR=15`（Chudnovsky Rank-15） | `4,096` pairs `PASS` |
| `selftest_qroam_vslot0.kmb` | `23` | `CCX=35, HMR=35` | `2,368` shots `PASS` |
| `selftest_zpow128.kmb` | `4` | `Z_POW=4, T_LIKE=2` | 128-bit 精度 `PASS` |
| `selftest_hmr_lut.kmb` | `36` | `CCX=47, HMR=59` | `4,096` shots `PASS` |
| `selftest_linf2.kmb` | `12` | `CCX=0, @clean=0` | `4,096` states `PASS` |
| `selftest_cuccaro.kmb` | `34` | `CCX=46, @clean=1` | `4,096` shots `PASS` |

---

## 3. 常用指令

```bash
cd /Users/kuangtalin/Documents/QUA/kaggle
make all                  # 編譯 qua_synth、dft_bootstrap 並執行 test
make test                 # 執行 qua_synth、8 項 Kickmix Fuzz 驗證與 Kickmix pytest
make verify-qcc           # 聯動驗證 ../qcc T1..T21、5 個範例與 selftest 雙種子檢查
make stage-init           # 從 sample_submission 初始化暫存目錄
make stage-sync           # 同步 build/invsqrt.kmb 與 dfthcsos 到 stage（缺 dqi 時沿用基準）
make score                # 執行官方 Verifier 快速模式並投影分數
make score-full           # 執行官方 Verifier 完整規模
make pack                 # 檢查競賽合約並產出 submissions/submission.zip
make clean                # 清除編譯與合成暫存檔
```
