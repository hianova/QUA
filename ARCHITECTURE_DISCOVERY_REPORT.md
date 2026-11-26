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
| **`FLOW-QCtrl Ring-0`** | **`32 mrad`** | **`0.43%`** | **`d=19–29`** | **解鎖 `dc=5` Cultivation** | **4/7 基準進入 Cultivation，時空體積縮減 `1.97×–4.55×`** |

### 2. 誠實的經濟斷崖線：$T_{\text{eq}} \lesssim 250\text{–}300$ 與結構化前沿
閉環測試明確畫出了 Gen-2 硬體的真實可行域，而非盲目樂觀：
- **進入 `dc=5` Cultivation 的 4 大基準**（Gen-2 + Ring-0）：
  - `Ising4x4_Trotter`（`88 T-eq`）：`108,560 -> 14,752` 實體位元（**`4.55×`** 時空縮減）
  - `Plaquette4x4_LGT`（`169 T-eq`）：`108,560 -> 14,752` 實體位元（**`4.54×`** 時空縮減）
  - `GF256_SEARU_A32`（`128 T-eq`）：`402,560 -> 147,088` 實體位元（**`3.10×`** 時空縮減）
  - `KSAT12_Grover`（`116 T-eq`）：`213,910 -> 64,122` 實體位元（**`1.97×`** 時空縮減）
- **尚未跨越 Cultivation 門檻的對照組**：
  - `QAOA_MaxCut12_Base -> PopLUT Opt` 雖藉 4-bit LUT 分塊把邏輯位元由 `151 -> 79`、`CCX` 由 `136 -> 106`、Gen-2 實體位元由 `307,982 -> 186,878`（`-39.3%`），但因 `T-eq = 433` 仍略高於此雜訊下的單 T 容錯預算（`~1.15e-5`），誠實回退至 `15-to-1` 蒸餾（提示需進一步用 Wallace 3:2 壓縮跨過 `250 T-eq` 線）。
  - `AES128_10R_Fwd`（`43,200 T-eq`）仍需 **`10,911,175`** 顆實體位元，證明大規模密碼分析短期不具經濟性。

### 3. Pure-GPR Ring-0 即時表面碼解碼（無需 GPU / AI 推論延遲）
[`qctrl/qctrl_decoder_ring0.h`](qctrl/qctrl_decoder_ring0.h) 以純整數 freestanding Pure-GPR（`0` 浮點/SIMD、`0` `udiv`/`sdiv`、`0` heap、`0` 間接跳躍）實作 `d = 3, 5, 7` 精確 Bitmask-DP MWPM：
- `X+Z` 雙組平均解碼延遲：**`2.12 ns`（`d=3`） / `4.22 ns`（`d=5`） / `7.38 ns`（`d=7`）**，比 `1,000 ns` 表面碼週期硬性時限快 **135×–471×**。
- 窮舉 `39,966` 組 `wt <= (d-1)/2` Pauli 錯誤達成 **`0` 邏輯誤判**；`p = 0.010` 下測得指數抑錯 **`Λ(3→5) = 6.03×`** 與 **`Λ(5→7) = 5.80×`**；並內建 EWMA+CUSUM 症候漂移監控器（200 回靜默 `0` 誤報、TLS 突發 19 回觸發）。

### 4. 商業定位修正：從「賣通用量子算力」轉為「賣跨越門檻的控制＋解碼＋編譯鏟子」
短期（2026–2028）真正的商業與工程 Alpha 不在於自建百萬位元機組破解密碼，而在於任何超導晶片團隊要從 40q Break-Even 邁向百位元級 Cultivation 實驗時，都必須依賴能同時壓低連貫 `ZZ`、滿足 `< 1 µs` 即時解碼與自動壓縮 `T-eq` 的全端基礎設施。

---

## 三、下一階段深化路線（Roadmap）

1. **Frontier 1｜3D 時空容錯解碼（`d` 輪含測錯之純 GPR 3D MWPM）**：將解碼從 2D 單輪 Code-Capacity 推進至含測量雜訊的 `d+1` 切片時空圖，驗證所有權重 `<= t` 資料/量測混合故障 `0` 誤判與微秒內延遲。
2. **Frontier 2｜QAOA MaxCut-12 Wallace 3:2 壓縮跨越 Cultivation 懸崖**：以 1-CCX 3:2 Full-Adder 與 2:2 Half-Adder 欄位壓縮替代位元寬進位加法樹，將 `QAOA_MaxCut12` 壓至 `< 250 T-eq`，在系統內部驗證「僅靠編譯器優化就將演算法推入 Cultivation 經濟區間」。
