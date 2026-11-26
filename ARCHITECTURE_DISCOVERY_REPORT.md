# `QUA` 全端架構演進與經濟性邊界戰報

> 本文整合 [`FTQC_CODESIGN_REPORT.md`](FTQC_CODESIGN_REPORT.md) 與 [`qctrl/DEMO_REPORT.md`](qctrl/DEMO_REPORT.md) 之純 C17 / C++20 實測數據，說明為何單層物理優化會撞牆，而將 **Ring-0 即時自主校準 × 表面碼實時解碼 × QCC E-Graph 編譯 × Magic State Cultivation** 串成閉環後，會在跨層交會處產生明確的成本相變。

---

## 一、為何 `PHY` 單層撞牆，跨層 `QUA` 卻長出可落地突破？

- **單層物理極限不可繞過**：基礎退相干（$T_1 = 530\,\mu\text{s}$、$T_2 = 300\,\mu\text{s}$）與雙量子閘保真度 $F_{2Q}$ 屬於元件物理下限；若只從單一量子電路角度對 `AES-128` 這類大算術電路微調 10%–20% 邏輯閘數，在百萬至千萬實體位元尺度下仍不具備商業意義。
- **跨層組合相變（Cross-Layer Phase Transition）**：透過跨層模型推導與實驗驗證（見 [`Time`](../Time) 對 2026–2028 前沿方向的定位：`Magic State Cultivation`、`ApproxOpt × QEC`、`Many-Body × 2Q-Gate`），當「物理校準壓低連貫串擾」、「表面碼實時解碼壓進周期時限」與「編譯器壓低 `T-eq`」同時滿足門檻時，系統會從 **高門檻不可運行／15-to-1 巨型蒸餾** 直接跳躍進 **Cultivation 輕量化視窗**，出現 **2×–4.5×** 的斷崖式實體時空體積縮減。

---

## 二、實測確立的四項關鍵工程與經濟發現

### 1. Ring-0 自主校準是開啟表面碼與 Cultivation 的二元開關
依 [`FTQC_CODESIGN_REPORT.md`](FTQC_CODESIGN_REPORT.md) 在 Gen-2 目標硬體（$F_{2Q} = 99.85\%$, $T_1 = 530\,\mu\text{s}$, $T_2 = 300\,\mu\text{s}$）之實測閉環：

| 控制策略 | `ZZ` RMS | 有效雜訊 $p_{\text{eff}}$ | 表面碼距 $d$ | 魔術態方案 | 實測結論 |
| :--- | ---: | ---: | :---: | :--- | :--- |
| `No-Cal (Pure QEC)` | `156 mrad` | `1.99%` | `ABOVE_TH` | 被阻擋（$p_{\text{eff}} \ge 1.0\%$） | 即使硬件達 `99.85%`，不壓連貫 `ZZ` 仍無法容錯 |
| `Periodic Offline (q20 ep)` | `70 mrad` | `0.69%` | `d=35` | 15-to-1 蒸餾 | 仍高於 Cultivation 雜訊窗口，需龐大蒸餾工廠 |
| `Single-Step 3sig Schmitt` | `57 mrad` | `0.58%` | `d=35` | 15-to-1 蒸餾 | 靜默誤觸 `180/6700`，無法解鎖 Cultivation |
| **`FLOW-QCtrl Ring-0`** | **`32 mrad`** | **`0.43%`** | **`d=19–29`** | **解鎖 `dc=5` Cultivation** | **5/8 基準進入 Cultivation，時空體積縮減 `1.97×–4.55×`** |

### 2. 誠實的經濟斷崖線：$T_{\text{eq}} \lesssim 250\text{–}300$ 與結構化前沿
閉環測試明確畫出了 Gen-2 硬體的真實可行域，而非盲目樂觀：
- **進入 `dc=5` Cultivation 的 5 大基準**（Gen-2 + Ring-0）：
  - `Ising4x4_Trotter`（`88 T-eq`）：`108,560 -> 14,752` 實體位元（**`4.55×`** 時空縮減）
  - `Plaquette4x4_LGT`（`169 T-eq`）：`108,560 -> 14,752` 實體位元（**`4.54×`** 時空縮減）
  - `GF256_SEARU_A32`（`128 T-eq`）：`402,560 -> 147,088` 實體位元（**`3.10×`** 時空縮減）
  - `KSAT12_Grover`（`116 T-eq`）：`213,910 -> 64,122` 實體位元（**`1.97×`** 時空縮減）
  - **Stage C 新增突破 — `QAOA_MaxCut12_Wallace`（`73 T-eq`、`T-depth = 4`）**：以 3:2 Wallace Full-Adder 欄位壓縮改寫（[`qcc/qir_qaoa_mb.h`](qcc/qir_qaoa_mb.h) `T23d_QAOA_MaxCut12_Wallace_Opt`），將邏輯位元設為 `101`（換取 `CCX` 由 `136 -> 16`、`T-eq` 由 `553 -> 73`），在 Gen-2 Ring-0 下距離降至 `d=21`、總實體位元 **`92,282`**（`89,082 + 3,200`），成功跨越門檻解鎖 `Cultivation (dc=5, P=40%)`（**`3.38×`** 時空縮減），將 Gen-2 解鎖基準數推升至 **5/8**！
- **尚未跨越 Cultivation 門檻的對照組**：
  - `QAOA_MaxCut12_Base -> PopLUT Opt` 雖藉 4-bit LUT 分塊把邏輯位元由 `151 -> 79`、`CCX` 由 `136 -> 106`、Gen-2 實體位元由 `307,982 -> 186,878`（`-39.3%`），但因 `T-eq = 433` 仍高於該雜訊預算，回退至 `15-to-1` 蒸餾；與 `Wallace 3:2` 的對照直接在系統內部證明「純編譯器 `T-eq` 壓縮即可將演算法推入 Cultivation 經濟區間」。
  - `AES128_10R_Fwd`（`43,200 T-eq`）仍需 **`10,911,175`** 顆實體位元，證明大規模密碼分析短期不具經濟性。

### 3. Pure-GPR Ring-0 即時 2D / 3D 表面碼解碼（無需 GPU / AI 推論延遲）
- **2D Code-Capacity 精確解碼（[`qctrl/qctrl_decoder_ring0.h`](qctrl/qctrl_decoder_ring0.h)）**：純整數 freestanding Pure-GPR（`0` 浮點/SIMD、`0` `udiv`/`sdiv`、`0` heap、`0` 間接跳躍）實作 `d = 3, 5, 7` 精確 Bitmask-DP MWPM：
  - 最新基準平均解碼延遲：**`2.14 ns/rnd`（`d=3`） / `3.90 ns/rnd`（`d=5`） / `7.26 ns/rnd`（`d=7`）**，比 `1,000 ns` 表面碼週期硬性時限快 **137×–467×**。
  - 窮舉 `39,966` 組 `wt <= (d-1)/2` Pauli 錯誤達成 **`0` 邏輯誤判**；`p = 0.010` 下測得指數抑錯 **`Λ_2D(3→5) = 6.03×`** 與 **`Λ_2D(5→7) = 5.80×`**；並內建 EWMA+CUSUM 症候漂移監控器（200 回靜默 `0` 誤報、TLS 突發 19 回觸發）。
- **Stage B 新增突破 — 3D 時空容錯解碼（[`qctrl/qctrl_decoder3d_ring0.h`](qctrl/qctrl_decoder3d_ring0.h)）**：同以 Pure-GPR 實作 Safe-Component BFS（`d3d < bdist(i)+bdist(j)`）＋ Bitmask-DP（`k <= 10`）＋ Savings Greedy ＋ 2-Opt，於 `d` 輪含資料與量測雜訊的時空圖即時解碼：
  - 驗證 `40,510` 組 `data+meas` 混合故障（`d=3, 5` 全窮舉 `wt <= t` ＋ `5,000` 組 `d=7` `wt=3` 三元組）達成 **`0` 邏輯誤判**。
  - 3D 單輪平均延遲：**`5.00 ns/rnd`（`15.0 ns/blk`，`d=3`） / `11.86 ns/rnd`（`59.3 ns/blk`，`d=5`） / `92.96 ns/rnd`（`650.7 ns/blk`，`d=7`）**，高階 `d=7` 仍比 `1,000 ns` 硬性時限快 **`10.8×`**。
  - **誠實縮放極限觀測**：`p_data = p_meas = 0.002` 時 `d=3/5/7` 失效率為 `1.300e-03 / 1.500e-04 / 0.000e+00`；`p = 0.006` 下測得 **`Λ_3D(3→5) = 4.00×`**（`8.000e-03 -> 2.000e-03`），但 `d=7` 為 `3.400e-03`（**`Λ_3D(5→7) = 0.59×`**）。原因在於 `7` 輪含測錯共 `511` 個故障點在高雜訊下會產生 `k > 10` 大連通叢集，觸發 Savings+2-Opt 近似路徑而非精確 Blossom，明確指出下一代 3D `d >= 7` 必須升級為純整數交替樹 Sparse Blossom。

### 4. 商業定位修正：從「賣通用量子算力」轉為「賣跨越門檻的控制＋解碼＋編譯鏟子」
短期（2026–2028）真正的商業與工程 Alpha 不在於自建百萬位元機組破解密碼，而在於任何超導晶片團隊要從 40q Break-Even 邁向百位元級 Cultivation 實驗時，都必須依賴能同時壓低連貫 `ZZ`、滿足 `< 1 µs` 即時解碼與自動壓縮 `T-eq` 的全端基礎設施。

---

## 三、深化進度與下一階段路線（Roadmap）

1. ~~**Frontier 1｜3D 時空容錯解碼（`d` 輪含測錯之純 GPR 3D MWPM）**~~：**已完成**（[`qctrl/qctrl_decoder3d_ring0.h`](qctrl/qctrl_decoder3d_ring0.h)，`40,510` 筆 `0` 誤判、`5.00–92.96 ns/rnd`、`Λ_3D(3→5)=4.00×`）。
2. ~~**Frontier 2｜QAOA MaxCut-12 Wallace 3:2 壓縮跨越 Cultivation 懸崖**~~：**已完成**（[`qcc/qir_qaoa_mb.h`](qcc/qir_qaoa_mb.h) `T23d`，`T-eq: 553 -> 73`、Gen-2 Ring-0 `92,282` Q、解鎖 Cultivation **`3.38×`**）。
3. **Frontier 3｜Pure-GPR 交替樹 Sparse Blossom（突破 3D `d >= 7` `k > 10` 高雜訊簇）**：以零 heap 固定陣列交替樹 Blossom 取代 `k > 10` Savings+2-Opt 回退路徑，在維持 Pure-GPR `< 1 µs/rnd` 硬實時性的同時恢復 `Λ_3D(5→7) > 1` 指數抑錯。
