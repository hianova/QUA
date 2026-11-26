# `pqc` — PQC / NIST 安全基準量子密碼分析資源編譯器

**`pqc`** 基於 [`../qcc`](../qcc/README.md) 與 [`../../FLOW`](../../FLOW) 建置，專攻**後量子密碼學（PQC）與 NIST 安全基準之可驗證、零垃圾（Zero-Garbage）量子密碼分析 Oracle 合成與容錯資源估算**。

所有電路均透過宣告式 SSA QIR（[`include/pqc_oracles.h`](include/pqc_oracles.h)）交由 `qcc_optimize_module` / `qir_codegen_emit` 自動最佳化與降階（零手動閘級發射），並以 `QCC` 的 `64-lane` 雙種子位元切片模擬器結合獨立經典參考實作與 KAT 已知答案向量進行端對端相位與數值驗證。

---

## 支援的 5 大密碼分析家族

| 家族 | 原語與規格 | 核心編譯技術 |
| :--- | :--- | :--- |
| **NIST 對稱加密基準** | `FIPS-197` AES-128 / AES-192 / AES-256（10 / 12 / 14 輪與完整 KeySchedule） | $\text{GF}(2^8)$ 三層代數塔反元素（36 / 18 CCX）+ `R1` 融合 ShiftRows/MixColumns/Affine（0 CCX） |
| **Hash / Sponge（SLH-DSA / SHA-3）** | `FIPS 180-4` SHA-256（64 輪）、`NIST SP 800-232` Ascon-p[12]、Keccak-f[400]（20 輪） | 1-AND `Ch`/`Maj`、`R3` Carry-Save Adder Tree、位元切片 `AND_VEC` + `LINF2` |
| **Code-Based PQC** | Classic McEliece `mceliece348864`（$\text{GF}(2^{12})$, `0x1009`）Goppa 求值與反演 + ISD | Chudnovsky Rank-42 $\text{GF}(2^{12})$ 乘（42 CCX）+ 代數塔反演（66 CCX）+ 可逆 popcount 樹 |
| **Lattice-Based PQC** | `FIPS-203` ML-KEM / Kyber（$q = 3329$）Module-LWE 有限距離搜尋 | 可逆 $\bmod\,3329$ Cuccaro 受控減加、3-bit 小秘密 QROM 乘表與無窮範數界限檢查 |
| **Multivariate Quadratic** | UOV / MAYO 類型密集二次系統（$n = 24, m = 24$ over $\mathbb{F}_2$） | 雙線性移位群 `AND_VEC`（276 CCX = 每二次單項式 1 CCX）+ `LINF2` + `EQ_ZERO` 相位標記 |

---

## 資源估算模型（[`include/pqc_nist_est.h`](include/pqc_nist_est.h)）

- 純整數 / Q16.16 定點數 `log2` 計算（零 `<math.h>` 依賴）。
- 同時回報 **Gidney 暫態 AND 等效 T 閘數**（`4 * CCX + T`）與 **酉 Toffoli 等效 T 閘數**（`7 * CCX + T`），以及由 `.kmb` 有向無環圖（DAG）精確計算的 `toffoli_depth`、`circuit_depth` 與 `reaction_depth`。
- 依 Zalka / Jaques–Naehrig–Roetteler–Virdia 平行 Grover 折衷模型輸出 `MAXDEPTH = 2^40 / 2^64 / 2^96` 下的實例數 $P$、總邏輯 Qubit、總 `T` 閘與 $G \cdot D$ 指標。

---

## 實測編譯基準與 NIST `MAXDEPTH` 資源估算（`make bench`）

以下數據由 `./build/pqc bench build --json build/pqc_report.json --md build/pqc_report.md` 實測產出（所有電路皆通過 64-lane 雙種子 + 經典參考 / KAT 驗證，`D_react = 1`）：

| Oracle / Benchmark | Family | k | Q | CCX | T-eq (Gidney) | T-depth (`D_toff`) | Circ-depth (`D_circ`) | Serial T | MD `2^40` T | MD `2^64` T | MD `2^64` Q | MD `2^96` T | NIST Cat |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `aes128_fwd` | AES-Cat1 | 128 | 4427 | 10800 | 43200 | 8703 | 64635 | `2^79.07` | `2^115.81` | `2^91.81` | `2^37.59` | `2^79.07` | Cat 1 |
| `aes128_grover` | AES-Cat1 | 128 | 4543 | 11053 | 44212 | 8956 | 65360 | `2^79.10` | `2^115.88` | `2^91.88` | `2^37.71` | `2^79.10` | Cat 1 |
| `aes192_grover` | AES-Cat3 | 192 | 5311 | 12349 | 49396 | 10654 | 77965 | `2^111.27` | `2^180.30` | `2^156.30` | `2^102.44` | `2^124.30` | Cat 3 |
| `aes256_grover` | AES-Cat5 | 256 | 6175 | 15157 | 60628 | 12351 | 90541 | `2^143.56` | `2^244.81` | `2^220.81` | `2^167.08` | `2^188.81` | Cat 5 |
| `sha256_grover` | Hash | 256 | 32319 | 39618 | 158472 | 16102 | 127832 | `2^144.93` | `2^246.56` | `2^222.56` | `2^170.24` | `2^190.56` | Cat 5 |
| `ascon128_grover` | Sponge | 128 | 20288 | 3969 | 15876 | 313 | 2709 | `2^77.65` | `2^109.66` | `2^85.66` | `2^30.33` | `2^77.65` | Cat 1 |
| `keccak400_grover` | Sponge | 128 | 35440 | 8253 | 33012 | 368 | 2020 | `2^78.68` | `2^110.92` | `2^86.92` | `2^31.58` | `2^78.68` | Cat 1 |
| `mceliece_goppa12` | Code-GF12 | 12 | 287 | 585 | 2340 | 584 | 10296 | `2^16.87` | `2^16.87` | `2^16.87` | `2^8.16` | `2^16.87` | Cat 0 |
| `mceliece_isd` | Code-ISD | 64 | 338 | 299 | 1196 | 297 | 1905 | `2^42.15` | `2^42.15` | `2^42.15` | `2^8.40` | `2^42.15` | Cat 0 |
| `mlkem_lwe3329` | Lattice | 24 | 2535 | 4735 | 18940 | 4636 | 22220 | `2^25.87` | `2^25.87` | `2^25.87` | `2^11.31` | `2^25.87` | Cat 0 |
| `mq24_grover` | MQ | 24 | 719 | 333 | 1332 | 57 | 545 | `2^22.12` | `2^22.12` | `2^22.12` | `2^9.49` | `2^22.12` | Cat 0 |

- 重點指標觀察：`aes128_grover` 透過 QCC 三層代數塔 `GF(((2^2)^2)^2)` 反元素（36 / 18 CCX）與 `R1` 線性層融合，達到 `Q = 4,543`、`CCX = 11,053`（`T-eq = 44,212`）與 `T-depth = 8,956`；`MAXDEPTH = 2^64` 下為 `2^91.88 T`、`2^37.71` 邏輯 Qubit。
- 非線性代數度低且並行的海綿結構與多變量系統（`ascon128_grover` `T-depth = 313`、`keccak400_grover` `T-depth = 368`、`mq24_grover` `T-depth = 57`）在 QCC 自動排程下展現極低關鍵路徑深度。

---

## 全端 FTQC 協同設計（[`include/pqc_ftqc_cultivation.h`](include/pqc_ftqc_cultivation.h)）

- 以純整數定點（零 `<math.h>`）將 `../qctrl/build/qctrl_ftqc_profile.json` 的實測物理連貫 `ZZ` 雜訊、`d=3/5/7` Ring-0 解碼延遲與各 `.kmb` 基準電路深度/`T-eq` 結合，對比 **Magic State Cultivation (`dc=3/5`)** 與 **15-to-1 蒸餾** 的實體時空體積。
- 執行 `make ftqc` 自動產生 [`../FTQC_CODESIGN_REPORT.md`](../FTQC_CODESIGN_REPORT.md) 與 `build/ftqc_codesign.json`，`selftest` 亦包含 `KAT_AND_FTQC_Invariants` 驗證閾值阻擋、單調性與 `N_T = 10^12` 溢位保護。

---

## 建置與執行指令

```bash
cd pqc
make all      # 編譯 build/pqc
make test     # 執行 C17/C++20 標頭檢查、KAT + FTQC 不變量 + 64-lane 雙種子自測與 qcc verify
make bench    # 編譯全系列基準並輸出 build/*.kmb、build/pqc_report.json 與 build/pqc_report.md
make ftqc     # 產出全端 Cultivation 資源共設計戰報 ../FTQC_CODESIGN_REPORT.md 與 JSON
./build/pqc analyze build/aes128_grover.kmb --bits 128 --name AES128
```
