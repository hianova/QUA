# `./kaggle` — Google Quantum AI: Verifiable Quantum Circuit Compilation

本工作區針對 Kaggle 競賽 **[Google Verifiable Quantum Circuit Compilation](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/data)** 建立，已全面重構為純 C17 容錯量子電路合成管線（[`src/qua_synth.c`](file:///Users/kuangtalin/Documents/qua/kaggle/src/qua_synth.c) 與標頭模組），搭配 Apple Silicon `Accelerate` 張量壓縮器 [`dft_bootstrap.c`](file:///Users/kuangtalin/Documents/qua/kaggle/dfthcsos_opt/dft_bootstrap.c)、通用量子編譯器 [`qcc/`](file:///Users/kuangtalin/Documents/qua/qcc/Makefile) 及 [`verify_harness.py`](file:///Users/kuangtalin/Documents/qua/kaggle/scripts/verify_harness.py) 驗證打包管線，並產出已通過合約驗證之 [`submissions/submission.zip`](file:///Users/kuangtalin/Documents/qua/kaggle/submissions/submission.zip)（`58,031,820 B`、`10` 個檔案）。

---

## 1. 官方驗證成績總覽（Geometric Mean Score = **`5.8230×`**）

競賽總分為三大 Benchmark 相對於官方 `sample_submission` 成本壓縮比的**幾何平均數（越高越好）**：

$$
\text{Score} = \left( \prod_{b=1}^{3} \frac{\text{cost}_{b}^{\text{ref}}}{\text{cost}_{b}} \right)^{1/3} = (9.1098 \times 12.5800 \times 1.7228)^{1/3} = \mathbf{5.8230\times}
$$

| Benchmark | 官方基準成本 $\text{cost}_b^{\text{ref}}$ | 實測驗證成本 $\text{cost}_b$ | 提升倍率 | 合成架構與實測規格 |
| :--- | ---: | ---: | ---: | :--- |
| **1. `invsqrt`**<br>([`qua_invsqrt.h`](file:///Users/kuangtalin/Documents/qua/kaggle/src/qua_invsqrt.h)) | `2,076,774.23`<br>(`Q=260, CCX=7,735, rot=37`) | **`227,972.06`**<br>(`Q=99, CCX=2,044, T=1, rot=33, C_eq=2,302.75`) | **`9.1098×`**<br>(`6,266,904` shots `PASS`) | • `K_m=8` 分段 2-Jet 查表多項式（`w=(28,15,6)`, `w_res=28`, `b_sh_sq=1, b_sh_main=1, zprod=1`）<br>• 130-bit 精確整數參考 (`isqrt`) 與雙階段區間 minimax 打磨：`D_react=251`，`worst_error_over_tolerance = 0.996960`（`worst_x = 2,897,863`） |
| **2. `dqi`**<br>([`qua_dqi.h`](file:///Users/kuangtalin/Documents/qua/kaggle/src/qua_dqi.h)) | `24,439,470,855,968.22`<br>(`Q=1,893, T=5,725,110`) | **`1,942,717,385,781.65`**<br>(競賽 `102,400`-shot / $\log_2 \delta=-128$ 投影；`4,096`-shot 有限樣本上界 `4,867,586,158,613.31`) | **`12.5800×`**<br>(`4,096` shots `0` 錯誤 `PASS`) | • Half-EEA + Vandermonde + Toeplitz (`toep_base=3`) + Chien-Forney 解碼管線（`R1=23, S1=24, K_MIN=23, HALT_START=17`）<br>• Rank-42 $\text{GF}(2^{12})$ / Rank-15 $\text{GF}(64)$ 原生雙線性乘法器與區域 Clifford/CCZ 窺孔消除（`-618,260`）：`Q=1,677, T_avg=579,876.9, T_max=579,956, HMR=703,559, D_react=33,971` |
| **3. `dfthcsos`**<br>([`qua_dfthcsos.h`](file:///Users/kuangtalin/Documents/qua/kaggle/src/qua_dfthcsos.h)<br>[`dft_bootstrap.c`](file:///Users/kuangtalin/Documents/qua/kaggle/dfthcsos_opt/dft_bootstrap.c)) | `3,355,599,588.32`<br>(`Q=730, T=93,237, λ=49.301`) | **`1,947,735,832.79`**<br>(`Q=740, T=53,446.0, λ_eff=49.2474`) | **`1.7228×`**<br>(`204,544` shots `PASS`) | • `Accelerate` 張量壓縮量化 `(R=52, B=60, C=37, brot=10)` 與二進位表 [`words.bin`](file:///Users/kuangtalin/Documents/qua/kaggle/dfthcsos_opt/words.bin) / [`submission.npz`](file:///Users/kuangtalin/Documents/qua/kaggle/dfthcsos_opt/submission.npz)<br>• Virtual-Slot-0 QROAM (`J=11, K=7`, `borrow_clean=1`) 產出 7 張電路：`angles_Q=126 <= 128, angles_T=44,716, D_react=77` |

---

## 2. 專案結構與核心檔案

- [Makefile](file:///Users/kuangtalin/Documents/qua/kaggle/Makefile)：編譯 `build/qua_synth` 與 `build/dft_bootstrap`，並驅動自測、驗證與打包。
- **C17 核心合成器 (`src/`)**：
  - [`qua_synth.c`](file:///Users/kuangtalin/Documents/qua/kaggle/src/qua_synth.c)：CLI 總驅動，輸出 8 個基礎原語自測電路與三大 Benchmark 提交電路。
  - [`qua_kmb.h`](file:///Users/kuangtalin/Documents/qua/kaggle/src/qua_kmb.h)：`.kmb` v2 編碼器與相位輔助函式（包裝 [`qcc_kmb.h`](file:///Users/kuangtalin/Documents/qua/qcc/qcc_kmb.h)）。
  - [`qua_bitop.h`](file:///Users/kuangtalin/Documents/qua/kaggle/src/qua_bitop.h)：Bilevel SEARU × 2-Jet 搜尋、Virtual-Slot-0 QROAM、`HMR_Dual` parity LUT、`LinF2_Fuse` 與 Cuccaro 單輔助受控加法器。
  - [`qua_ele.h`](file:///Users/kuangtalin/Documents/qua/kaggle/src/qua_ele.h)：SSA pebbling、異質位寬裝箱、Rank-42 $\text{GF}(2^{12})$ / Rank-15 $\text{GF}(64)$ 乘法、區域窺孔優化與 `qubit` 緊緻化。
  - [`qua_invsqrt.h`](file:///Users/kuangtalin/Documents/qua/kaggle/src/qua_invsqrt.h) / [`qua_dqi.h`](file:///Users/kuangtalin/Documents/qua/kaggle/src/qua_dqi.h) / [`qua_dfthcsos.h`](file:///Users/kuangtalin/Documents/qua/kaggle/src/qua_dfthcsos.h)：三大 Benchmark C17 搜尋與電路合成器。
- **DF-THC 張量量化與持久化資料 (`dfthcsos_opt/`)**：
  - [`dft_bootstrap.c`](file:///Users/kuangtalin/Documents/qua/kaggle/dfthcsos_opt/dft_bootstrap.c)：以 Apple `ACCELERATE_NEW_LAPACK` 進行張量分解與 Givens 座標下降量化。
  - [`words.bin`](file:///Users/kuangtalin/Documents/qua/kaggle/dfthcsos_opt/words.bin) / [`submission.npz`](file:///Users/kuangtalin/Documents/qua/kaggle/dfthcsos_opt/submission.npz)：已驗證之 `(52, 60, 37, brot=10)` 基準解。
- **驗證與打包**：
  - [`scripts/verify_harness.py`](file:///Users/kuangtalin/Documents/qua/kaggle/scripts/verify_harness.py)：C17 原語隨機 Fuzz 驗證、Kickmix 雙種子交叉驗證、官方 verifier 評分投影、`words.bin` 匯出與 [`submissions/submission.zip`](file:///Users/kuangtalin/Documents/qua/kaggle/submissions/submission.zip) 嚴格合約打包。

### 常用指令

```bash
cd /Users/kuangtalin/Documents/qua/kaggle
make all             # 編譯 qua_synth 與 dft_bootstrap，產生 selftest 與三大 Benchmark .kmb
make test            # 編譯並執行 8 項 C17 原語 Fuzz 驗證與 3,402 項 Kickmix pytest
make verify-qcc      # 聯動驗證相鄰 ../qcc 編譯產物與 Kickmix 雙種子交叉檢查
make score           # 執行三大官方 Verifier 快速模式並投影 Leaderboard 幾何平均倍數
make score-full      # 執行競賽完整規模官方 Verifier
make pack            # 檢查所有競賽合約（<=128 Qubits、dqi/dfthcsos 禁用 Z_POW 等）並產出 submissions/submission.zip
```

---

## 3. 通用量子編譯器基礎設施（[`qcc/`](file:///Users/kuangtalin/Documents/qua/qcc/Makefile)）

通用中間表示與編譯最佳化模組位於相鄰 [`../qcc/`](file:///Users/kuangtalin/Documents/qua/qcc/Makefile)：

- [`qcc_arena.h`](file:///Users/kuangtalin/Documents/qua/qcc/qcc_arena.h) / [`qcc_kmb.h`](file:///Users/kuangtalin/Documents/qua/qcc/qcc_kmb.h)：`.kmb` v2 雙向編解碼、遊標掃描與原語構造器。
- [`qir.h`](file:///Users/kuangtalin/Documents/qua/qcc/qir.h) / [`qir_dsl.h`](file:///Users/kuangtalin/Documents/qua/qcc/qir_dsl.h) / [`qir_lift.h`](file:///Users/kuangtalin/Documents/qua/qcc/qir_lift.h)：SSA 量子中間表示、`.qir` 解析與反編譯提升。
- [`qir_rewrite.h`](file:///Users/kuangtalin/Documents/qua/qcc/qir_rewrite.h) / [`qir_opt.h`](file:///Users/kuangtalin/Documents/qua/qcc/qir_opt.h) / [`qir_kernels.h`](file:///Users/kuangtalin/Documents/qua/qcc/qir_kernels.h) / [`qir_codegen.h`](file:///Users/kuangtalin/Documents/qua/qcc/qir_codegen.h) / [`qir_pass_kmb.h`](file:///Users/kuangtalin/Documents/qua/qcc/qir_pass_kmb.h) / [`qir_verify.h`](file:///Users/kuangtalin/Documents/qua/qcc/qir_verify.h)：飽和重寫、2-Jet 鬆弛與雙種子模擬驗證。
