# `qcc` — QCC v3 通用容錯量子電路編譯器

**QCC v3 (Quantum Circuit Compiler)** 是一套純 C17（且所有標頭皆相容 C++20 `-Wall -Wextra -Werror`）、零外部執行期依賴的容錯量子電路編譯器，建立於 `FLOW v8.0.0`（`core/` + `sdk/`）與 `SEARU` 之上（透過 `FLOW_DIR` 引用，無須複製標頭檔）。

---

## 1. 標頭模組架構

| 檔案 | 功能定位 |
| :--- | :--- |
| [qcc_arena.h](qcc_arena.h) | 線性 Arena 配置器、256-bit 位元向量與 128-bit 固定點相位角度輔助運算（Clifford 階層判斷） |
| [qcc_kmb.h](qcc_kmb.h) | Kickmix `.kmb` v2 二進位格式的雙向編解碼、遊標掃描、精確資源與關鍵路徑計數（`CCX`、`Z_POW`、`T_LIKE`、`HMR`、`D_react`、`D_toff`、`D_circ`）與暫存器管理 |
| [qir.h](qir.h) | SSA 量子中間表示（QIR）、64-lane 位元切片模擬器 `qir_eval_64`、zero-mask 抽象解釋與組建 API |
| [qir_fir.h](qir_fir.h) | `FLOW.Morphism` `.fir` QTT 線性閘口橋接器（嚴格驗證 `.linear` 准入並降階至 QIR） |
| [qir_lift.h](qir_lift.h) | 無尺度 `.kmb` v2 結構反編譯提升器（支援 `>256` qubit 分塊、Gidney temporary-AND、`LinF2` 與 `phase_poly` 島辨識） |
| [qir_rewrite.h](qir_rewrite.h) | E-Graph 等價飽和重寫引擎（規則 `R1`–`R11`）、256-bit 線性規範（Gauge）壓縮與 Rank-42 / Rank-27 / Rank-15 $\text{GF}(2^m)$ 代數表建構 |
| [qir_kernels.h](qir_kernels.h) | 容錯量子核心庫：Gidney temporary-AND、Cuccaro 單輔助受控加減、QROAM/PARITY-LUT、`CSWAP` 筒形位移、`0-CCX` `CZ_IF` 相位內積、3 層代數塔 $\text{GF}(2^{12})$ 反元素（`66`/`69` CCX）與 $\text{GF}(2^8)$ `GF(((2^2)^2)^2)` 代數塔反元素（`36`/`38` CCX，HMR 反計算 `18`/`20` CCX） |
| [qir_codegen.h](qir_codegen.h) | 物理暫存器配置、摺疊子圖降階、`@clean` 輔助暫存器輸出與 `.kmb` 閘序列發射 |
| [qir_opt.h](qir_opt.h) | `SEARU` × `FLOW` 2-Jet 混合最佳化、$O(N \log N)$ 最大堆積零垃圾級聯反計算、in-place/線路借用排程（支援至 65,536 qubits 與大規模電路快速重寫播種） |
| [qir_pass_kmb.h](qir_pass_kmb.h) | 閘級 `.kmb` 相位融合/自我反元素對消、`@clean` 裁剪壓縮，以及統一高階編譯 API（`qcc_compile_module_to_kmb` / `qcc_compile_module_to_file`） |
| [qir_verify.h](qir_verify.h) | 雙種子 `.kmb` / QIR 等價與無相位回踢驗證器（支援經典輸入、分離式 `init_mask`/`cmp_mask` 遮罩與容許近似誤差驗證） |
| [qcc_main.c](qcc_main.c) | CLI 工具（`opt` / `lift` / `verify` / `inspect` / `selftest`）與 `T1`–`T22` 自我測試套件（27 項測試） |

---

## 2. SSA QIR 與重寫規則（`R1`–`R11`）

### 支援的核心操作

- **算術與定點數**：`qir_emit_add`、`qir_emit_sub`、`qir_emit_clz_norm`、`qir_emit_sqr_fxp`、`qir_emit_mul_fxp`、`qir_emit_ctrl_sel`
- **有限體與線性/位元邏輯**：`qir_emit_gf2m_add`、`qir_emit_gf2m_mul`（含 Chudnovsky Rank-42 $\text{GF}(2^{12})$、Rank-27 $\text{GF}(2^8)$ 與 Rank-15 複合 $\text{GF}(64)$）、`qir_emit_gf2m_inv` / `qir_emit_gf2m_inv_ctrl`、可變位移 `qir_emit_shl_var` / `qir_emit_shr_var`、`qir_emit_eq_zero`、`qir_emit_linf2`、位元邏輯 `bit_or_reduce` 與經典暫存器 `qir_add_classical_input`
- **查表與相位**：`qir_emit_lut_rom`、`qir_emit_phase_poly`、`qir_emit_phase_dot` / `qir_emit_phase_dot_gf2m`
- **高階宣告式建構**：`approx_fn` 與 `qir_emit_dyadic_phase_fn`（`phase_fn`，含 Chebyshev 擬合、座標 minimax 定心與 128-bit 整數開方）

### E-Graph 重寫規則

- **R1 (`LinF2_Fuse`)**：將單次使用鏈結的多個 $\mathbb{F}_2$ 線性變換融合為單一 Gauss 消去線性網格。
- **R2 (`HMR_ParityLUT`)**：透過 `HMR` + 經典條件 `CZ_IF` 進行低費用查表（最多支援 13 個位址位元）與 Gidney temporary-AND 反計算。
- **R3 (`CSA_Tree`)**：將多運算元連續加法壓縮為 Carry-Save Adder 樹。
- **R4 (`Phase_AddTree_Split`)**：將饋入模組化線性 `phase_poly` 的 `ADD`/`SUB` 樹展開為 0-CCX 相位旋轉。
- **R5 (`LUT_ColSplit`) / R7 (`BasisIso_Gauge`)**：對寬輸出 QROM 進行分欄切割，並利用 256-bit $\mathbb{F}_2$ 列簡化（RREF）壓縮規範基底。
- **R6 (`GF2m_Tower_Rank`)**：將 $\text{GF}(2^{12})$ 乘法降階為 Rank-42 Chudnovsky（42 CCX）、$\text{GF}(2^8)$ 乘法降階為 Rank-27（27 CCX）、複合 $\text{GF}(64)$ 乘法降階為 Rank-15（15 CCX），並將 $\text{GF}(2^{12})$ 反元素降階為 `GF(((2^3)^2)^2)`（無條件 `66` / 受控 `69` CCX）、任意不可約多項式 $\text{GF}(2^8)$ 反元素降階為 `GF(((2^2)^2)^2)`（無條件 `36` / 受控 `38` CCX，HMR 反計算 `18` / `20` CCX）。
- **R8 / R8b (`QROAM_SelectSwap` / `ColSplit`)**：Differential Virtual-Slot-0 Select-Swap QROAM（`b_sh = 1..6`）與分欄 Select-Swap 融合。
- **R9 (`MUL_CONST_ADD` / `SQR_TRIANGLE`)**：常數乘加與對稱三角平方縮減。
- **R10a / R10b / R10c**：就地累加 `INPLACE_ARITH`、就地 Cuccaro 受控加法 `CTRL_ADD_CUCCARO`，以及 `CONCAT -> XOR_ACCUM -> SLICE` 零拷貝匯流排融合。
- **R11 (`GF2M_CONST_LINF2`)**：將常數有限體乘法化簡為 0-CCX 的 `LinF2` 線性網格。

---

## 3. 建置、CLI 用法與測試

### 建置與執行測試

```bash
cd qcc
make FLOW_DIR=../../FLOW clean all test examples
```

- `make check-headers` 會同時驗證 `cc -std=c17 -Wall -Wextra -Werror` 與 `c++ -std=c++20 -Wall -Wextra -Werror` 零警告編譯。
- `make test` 執行 `./build/qcc selftest`，驗證全部 `T1`–`T22`（27 項）測試。
- `make lean-fir examples` 從 [`Qua.lean`](../Qua.lean) 編譯並以 `FLOW fir-synth` 合成與雙種子驗證全部 `.fir` 目標：
  `gf256_sbox_mac.fir`、`rev_alu_fxp.fir`、`adder_csa_phase.fir`、`gf12_arith.fir`、`wide_lut_cap.fir`、`gf256_tower_spec.fir` $\to$ `gf256_tower_{a32,d4}.fir`（`CCX: 92 -> 32/33`, `T-depth: 88 -> 6`）、`lut5_bijection_spec.fir` $\to$ `lut5_bijection_opt.fir`（`CCX: 21 -> 5`, `T-depth: 13 -> 1`），並驗證拒斥 `.cartesian` 及與 `../ELE` 的 6 份跨專案 `.fir` 位元組一致性。

### CLI 指令

```bash
./build/qcc --help
./build/qcc opt <in.kmb|in.fir> -o <out.kmb> [--max-q N] [--metric default|qc|q2t|qt] [--iters N] [--verify-batches N] [--clean-anc] [--emit-clean-reg]
./build/qcc verify <circuit.kmb> [--ref <ref.kmb|ref.fir>] [--states N] [--in-regs N]
./build/qcc inspect <circuit.kmb>
./build/qcc selftest
```

### C17 / C++20 程式庫呼叫介面

透過 `#include "qir_pass_kmb.h"` 即可將任意組建好的 `qir_module_t` 編譯為 `.kmb` v2 電路或寫入檔案：

```c
qcc_opt_cfg_t cfg = qcc_opt_cfg_default();
cfg.require_clean_scratch = true;
qcc_kmb_stats_t st;
bool ok = qcc_compile_module_to_file(&mod, &arena, &cfg, "out.kmb", &st);
```
