# QUA Full-Stack FTQC Co-Design Report (`qctrl Ring-0` × `qcc` × `Magic State Cultivation`)

- **Profile Source**: `../qctrl/build/qctrl_ftqc_profile.json` (live measured qctrl bench40, Ring-0 40q sweep = `17431 ns`)
- **Fixed-Point Engine**: Pure C17 integer `pqc_ftqc_cultivation.h` (zero `<math.h>`)

## 1. Verified Logical Benchmarks (`qcc` + `pqc` Compiled `.kmb` Artifacts)

| Benchmark | Frontier Domain | Logical Q | CCX | Z_POW / T | T-eq (Cult/Gidney) | T/Circ Depth |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: |
| `QAOA_MaxCut12_Base` | ApproxOpt x QEC | 151 | 136 | 5 | 553 | 87 |
| `QAOA_MaxCut12_Opt` | PopLUT Opt QAOA | 79 | 106 | 5 | 433 | 85 |
| `QAOA_MaxCut12_Wallace` | Wallace 3:2 QAOA | 101 | 16 | 5 | 73 | 4 |
| `Ising4x4_Trotter` | Many-Body 2Q-Gate | 16 | 0 | 64 | 88 | 7 |
| `Plaquette4x4_LGT` | Lattice Gauge 4B | 16 | 0 | 58 | 169 | 12 |
| `KSAT12_Grover` | NP Combinatorial | 59 | 29 | 0 | 116 | 15 |
| `GF256_SEARU_A32` | Post-Quantum SBox | 128 | 32 | 0 | 128 | 6 |
| `AES128_10R_Fwd` | PQC AES-128 | 4427 | 10800 | 0 | 43200 | 8703 |

## 2. Ring-0 Real-Time QEC Decoder & Cycle Deadline Budget (`1,000 ns`)

| Decoder Mode & Distance | Measured Ring-0 Decode Latency | `1,000 ns` Budget Headroom | Below-Threshold Suppression `Λ` |
| :---: | :---: | :---: | :---: |
| 2D `d = 3` (Code Cap) | **`1.26 ns / round`** | 793.7x faster | baseline |
| 2D `d = 5` (Code Cap) | **`4.84 ns / round`** | 206.6x faster | **`Λ_2D(3→5) = 6.03x`** |
| 2D `d = 7` (Code Cap) | **`9.02 ns / round`** | 110.9x faster | **`Λ_2D(5→7) = 5.80x`** |
| 3D `d = 3` (`d` noisy rnds) | **`5.37 ns / round`** | 186.2x faster | baseline |
| 3D `d = 5` (`d` noisy rnds) | **`11.35 ns / round`** | 88.1x faster | **`Λ_3D(3→5) = 4.02x`** |
| 3D `d = 7` (`d` noisy rnds) | **`48.64 ns / round`** | 20.6x faster | **`Λ_3D(5→7) = 5.71x`** |

## 3. Gen-1 Baseline Hardware (`F_2Q = 99.60%`, `T1 = 530 us`, `T2 = 300 us`)

| Benchmark | Controller | `ZZ` RMS | `p_eff` | SC Dist `d` | Magic State Scheme | Total Phys Q (Data+Magic) | Spacetime Reduction |
| :--- | :--- | ---: | ---: | ---: | :--- | ---: | ---: |
| `QAOA_MaxCut12_Base` | `No-Cal (Pure QEC)` | `156 mrad` | `2.24%` | `ABOVE_TH` | Blocked (`p_eff >= 1.0%`) | - | - |
| `QAOA_MaxCut12_Base` | `Periodic Offline (q20 ep)` | `70 mrad` | `0.94%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 439310 (369950+69360) | 1.00x (baseline) |
| `QAOA_MaxCut12_Base` | `Single-Step 3sig Schmitt` | `57 mrad` | `0.83%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 439310 (369950+69360) | 1.00x (baseline) |
| `QAOA_MaxCut12_Base` | `FLOW-QCtrl Ring-0 (Huber+PS)` | `32 mrad` | `0.68%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 439310 (369950+69360) | 1.00x (baseline) |
| `QAOA_MaxCut12_Opt` | `No-Cal (Pure QEC)` | `156 mrad` | `2.24%` | `ABOVE_TH` | Blocked (`p_eff >= 1.0%`) | - | - |
| `QAOA_MaxCut12_Opt` | `Periodic Offline (q20 ep)` | `70 mrad` | `0.94%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 258575 (193550+65025) | 1.00x (baseline) |
| `QAOA_MaxCut12_Opt` | `Single-Step 3sig Schmitt` | `57 mrad` | `0.83%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 258575 (193550+65025) | 1.00x (baseline) |
| `QAOA_MaxCut12_Opt` | `FLOW-QCtrl Ring-0 (Huber+PS)` | `32 mrad` | `0.68%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 258575 (193550+65025) | 1.00x (baseline) |
| `QAOA_MaxCut12_Wallace` | `No-Cal (Pure QEC)` | `156 mrad` | `2.24%` | `ABOVE_TH` | Blocked (`p_eff >= 1.0%`) | - | - |
| `QAOA_MaxCut12_Wallace` | `Periodic Offline (q20 ep)` | `70 mrad` | `0.94%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 316810 (247450+69360) | 1.00x (baseline) |
| `QAOA_MaxCut12_Wallace` | `Single-Step 3sig Schmitt` | `57 mrad` | `0.83%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 316810 (247450+69360) | 1.00x (baseline) |
| `QAOA_MaxCut12_Wallace` | `FLOW-QCtrl Ring-0 (Huber+PS)` | `32 mrad` | `0.68%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 316810 (247450+69360) | 1.00x (baseline) |
| `Ising4x4_Trotter` | `No-Cal (Pure QEC)` | `156 mrad` | `2.24%` | `ABOVE_TH` | Blocked (`p_eff >= 1.0%`) | - | - |
| `Ising4x4_Trotter` | `Periodic Offline (q20 ep)` | `70 mrad` | `0.94%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 108560 (39200+69360) | 1.00x (baseline) |
| `Ising4x4_Trotter` | `Single-Step 3sig Schmitt` | `57 mrad` | `0.83%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 108560 (39200+69360) | 1.00x (baseline) |
| `Ising4x4_Trotter` | `FLOW-QCtrl Ring-0 (Huber+PS)` | `32 mrad` | `0.68%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 108560 (39200+69360) | 1.00x (baseline) |
| `Plaquette4x4_LGT` | `No-Cal (Pure QEC)` | `156 mrad` | `2.24%` | `ABOVE_TH` | Blocked (`p_eff >= 1.0%`) | - | - |
| `Plaquette4x4_LGT` | `Periodic Offline (q20 ep)` | `70 mrad` | `0.94%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 108560 (39200+69360) | 1.00x (baseline) |
| `Plaquette4x4_LGT` | `Single-Step 3sig Schmitt` | `57 mrad` | `0.83%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 108560 (39200+69360) | 1.00x (baseline) |
| `Plaquette4x4_LGT` | `FLOW-QCtrl Ring-0 (Huber+PS)` | `32 mrad` | `0.68%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 108560 (39200+69360) | 1.00x (baseline) |
| `KSAT12_Grover` | `No-Cal (Pure QEC)` | `156 mrad` | `2.24%` | `ABOVE_TH` | Blocked (`p_eff >= 1.0%`) | - | - |
| `KSAT12_Grover` | `Periodic Offline (q20 ep)` | `70 mrad` | `0.94%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 213910 (144550+69360) | 1.00x (baseline) |
| `KSAT12_Grover` | `Single-Step 3sig Schmitt` | `57 mrad` | `0.83%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 213910 (144550+69360) | 1.00x (baseline) |
| `KSAT12_Grover` | `FLOW-QCtrl Ring-0 (Huber+PS)` | `32 mrad` | `0.68%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 213910 (144550+69360) | 1.00x (baseline) |
| `GF256_SEARU_A32` | `No-Cal (Pure QEC)` | `156 mrad` | `2.24%` | `ABOVE_TH` | Blocked (`p_eff >= 1.0%`) | - | - |
| `GF256_SEARU_A32` | `Periodic Offline (q20 ep)` | `70 mrad` | `0.94%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 382960 (313600+69360) | 1.00x (baseline) |
| `GF256_SEARU_A32` | `Single-Step 3sig Schmitt` | `57 mrad` | `0.83%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 382960 (313600+69360) | 1.00x (baseline) |
| `GF256_SEARU_A32` | `FLOW-QCtrl Ring-0 (Huber+PS)` | `32 mrad` | `0.68%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 382960 (313600+69360) | 1.00x (baseline) |
| `AES128_10R_Fwd` | `No-Cal (Pure QEC)` | `156 mrad` | `2.24%` | `ABOVE_TH` | Blocked (`p_eff >= 1.0%`) | - | - |
| `AES128_10R_Fwd` | `Periodic Offline (q20 ep)` | `70 mrad` | `0.94%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 10911175 (10846150+65025) | 1.00x (baseline) |
| `AES128_10R_Fwd` | `Single-Step 3sig Schmitt` | `57 mrad` | `0.83%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 10911175 (10846150+65025) | 1.00x (baseline) |
| `AES128_10R_Fwd` | `FLOW-QCtrl Ring-0 (Huber+PS)` | `32 mrad` | `0.68%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 10911175 (10846150+65025) | 1.00x (baseline) |

## 4. Gen-2 Target Hardware (`F_2Q = 99.85%`, `T1 = 530 us`, `T2 = 300 us`)

| Benchmark | Controller | `ZZ` RMS | `p_eff` | SC Dist `d` | Magic State Scheme | Total Phys Q (Data+Magic) | Spacetime Reduction |
| :--- | :--- | ---: | ---: | ---: | :--- | ---: | ---: |
| `QAOA_MaxCut12_Base` | `No-Cal (Pure QEC)` | `156 mrad` | `1.99%` | `ABOVE_TH` | Blocked (`p_eff >= 1.0%`) | - | - |
| `QAOA_MaxCut12_Base` | `Periodic Offline (q20 ep)` | `70 mrad` | `0.69%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 439310 (369950+69360) | 1.00x (baseline) |
| `QAOA_MaxCut12_Base` | `Single-Step 3sig Schmitt` | `57 mrad` | `0.58%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 439310 (369950+69360) | 1.00x (baseline) |
| `QAOA_MaxCut12_Base` | `FLOW-QCtrl Ring-0 (Huber+PS)` | `32 mrad` | `0.43%` | `d=29` | 15-to-1 Distill (cult eps too hi) | 307982 (253982+54000) | 1.00x (baseline) |
| `QAOA_MaxCut12_Opt` | `No-Cal (Pure QEC)` | `156 mrad` | `1.99%` | `ABOVE_TH` | Blocked (`p_eff >= 1.0%`) | - | - |
| `QAOA_MaxCut12_Opt` | `Periodic Offline (q20 ep)` | `70 mrad` | `0.69%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 258575 (193550+65025) | 1.00x (baseline) |
| `QAOA_MaxCut12_Opt` | `Single-Step 3sig Schmitt` | `57 mrad` | `0.58%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 258575 (193550+65025) | 1.00x (baseline) |
| `QAOA_MaxCut12_Opt` | `FLOW-QCtrl Ring-0 (Huber+PS)` | `32 mrad` | `0.43%` | `d=29` | 15-to-1 Distill (cult eps too hi) | 186878 (132878+54000) | 1.00x (baseline) |
| `QAOA_MaxCut12_Wallace` | `No-Cal (Pure QEC)` | `156 mrad` | `1.99%` | `ABOVE_TH` | Blocked (`p_eff >= 1.0%`) | - | - |
| `QAOA_MaxCut12_Wallace` | `Periodic Offline (q20 ep)` | `70 mrad` | `0.69%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 316810 (247450+69360) | 1.00x (baseline) |
| `QAOA_MaxCut12_Wallace` | `Single-Step 3sig Schmitt` | `57 mrad` | `0.58%` | `d=33` | 15-to-1 Distill (cult eps too hi) | 289338 (219978+69360) | 1.00x (baseline) |
| `QAOA_MaxCut12_Wallace` | `FLOW-QCtrl Ring-0 (Huber+PS)` | `32 mrad` | `0.43%` | `d=21` | Cultivation (dc=5, P=40%) | 92282 (89082+3200) | **3.38x** |
| `Ising4x4_Trotter` | `No-Cal (Pure QEC)` | `156 mrad` | `1.99%` | `ABOVE_TH` | Blocked (`p_eff >= 1.0%`) | - | - |
| `Ising4x4_Trotter` | `Periodic Offline (q20 ep)` | `70 mrad` | `0.69%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 108560 (39200+69360) | 1.00x (baseline) |
| `Ising4x4_Trotter` | `Single-Step 3sig Schmitt` | `57 mrad` | `0.58%` | `d=29` | 15-to-1 Distill (cult eps too hi) | 80912 (26912+54000) | 1.00x (baseline) |
| `Ising4x4_Trotter` | `FLOW-QCtrl Ring-0 (Huber+PS)` | `32 mrad` | `0.43%` | `d=19` | Cultivation (dc=5, P=40%) | 14752 (11552+3200) | **4.55x** |
| `Plaquette4x4_LGT` | `No-Cal (Pure QEC)` | `156 mrad` | `1.99%` | `ABOVE_TH` | Blocked (`p_eff >= 1.0%`) | - | - |
| `Plaquette4x4_LGT` | `Periodic Offline (q20 ep)` | `70 mrad` | `0.69%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 108560 (39200+69360) | 1.00x (baseline) |
| `Plaquette4x4_LGT` | `Single-Step 3sig Schmitt` | `57 mrad` | `0.58%` | `d=31` | 15-to-1 Distill (cult eps too hi) | 84752 (30752+54000) | 1.00x (baseline) |
| `Plaquette4x4_LGT` | `FLOW-QCtrl Ring-0 (Huber+PS)` | `32 mrad` | `0.43%` | `d=19` | Cultivation (dc=5, P=40%) | 14752 (11552+3200) | **4.54x** |
| `KSAT12_Grover` | `No-Cal (Pure QEC)` | `156 mrad` | `1.99%` | `ABOVE_TH` | Blocked (`p_eff >= 1.0%`) | - | - |
| `KSAT12_Grover` | `Periodic Offline (q20 ep)` | `70 mrad` | `0.69%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 213910 (144550+69360) | 1.00x (baseline) |
| `KSAT12_Grover` | `Single-Step 3sig Schmitt` | `57 mrad` | `0.58%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 213910 (144550+69360) | 1.00x (baseline) |
| `KSAT12_Grover` | `FLOW-QCtrl Ring-0 (Huber+PS)` | `32 mrad` | `0.43%` | `d=23` | Cultivation (dc=5, P=40%) | 64122 (62422+1700) | **1.97x** |
| `GF256_SEARU_A32` | `No-Cal (Pure QEC)` | `156 mrad` | `1.99%` | `ABOVE_TH` | Blocked (`p_eff >= 1.0%`) | - | - |
| `GF256_SEARU_A32` | `Periodic Offline (q20 ep)` | `70 mrad` | `0.69%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 382960 (313600+69360) | 1.00x (baseline) |
| `GF256_SEARU_A32` | `Single-Step 3sig Schmitt` | `57 mrad` | `0.58%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 382960 (313600+69360) | 1.00x (baseline) |
| `GF256_SEARU_A32` | `FLOW-QCtrl Ring-0 (Huber+PS)` | `32 mrad` | `0.43%` | `d=23` | Cultivation (dc=5, P=40%) | 138624 (135424+3200) | **3.13x** |
| `AES128_10R_Fwd` | `No-Cal (Pure QEC)` | `156 mrad` | `1.99%` | `ABOVE_TH` | Blocked (`p_eff >= 1.0%`) | - | - |
| `AES128_10R_Fwd` | `Periodic Offline (q20 ep)` | `70 mrad` | `0.69%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 10911175 (10846150+65025) | 1.00x (baseline) |
| `AES128_10R_Fwd` | `Single-Step 3sig Schmitt` | `57 mrad` | `0.58%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 10911175 (10846150+65025) | 1.00x (baseline) |
| `AES128_10R_Fwd` | `FLOW-QCtrl Ring-0 (Huber+PS)` | `32 mrad` | `0.43%` | `d=35` | 15-to-1 Distill (cult eps too hi) | 10911175 (10846150+65025) | 1.00x (baseline) |

### Full-Stack Co-Design Takeaways (Auto-Derived from Measured Table)
- **Calibration gates Cultivation**: `No-Cal` (`156 mrad` `ZZ`) leaves `p_eff = 2.24%` (Gen-1) / `1.99%` (Gen-2) above the `1.00%` surface-code threshold (`ABOVE_TH`), blocking fault-tolerant scaling.
- **Ring-0 unlocks Cultivation where `T`-budget permits**: `FLOW-QCtrl Ring-0` (`32 mrad` `ZZ`) lowers Gen-2 `p_eff` to `0.43%`, unlocking Cultivation for **5/8** benchmarks (`QAOA_MaxCut12_Wallace` (`dc=5`, **3.38x**), `Ising4x4_Trotter` (`dc=5`, **4.55x**), `Plaquette4x4_LGT` (`dc=5`, **4.54x**), `KSAT12_Grover` (`dc=5`, **1.97x**), `GF256_SEARU_A32` (`dc=5`, **3.13x**)), while higher-budget workloads (`QAOA_MaxCut12_Base` (`553 T-eq`, `d=29`), `QAOA_MaxCut12_Opt` (`433 T-eq`, `d=29`), `AES128_10R_Fwd` (`43200 T-eq`, `d=35`)) stay on `15-to-1` distillation.
- **3-Stage Compiler QAOA Trajectory (PopLUT -> Wallace 3:2 crosses Cultivation Cliff)**: `Base -> PopLUT -> Wallace 3:2` compresses `QAOA_MaxCut12` logical footprint `Q: 151 -> 79 -> 101` and `T-eq: 553 -> 433 -> 73`, cutting Gen-2 Ring-0 physical footprint **`307982 -> 186878 -> 92282`** qubits (Cultivation UNLOCKED!).
- **Sub-microsecond Ring-0 2D & 3D decoding**: Pure-GPR MWPM latency spans `1.26–9.02 ns/rnd` (2D) and `5.37–48.64 ns/rnd` (3D noisy rounds, `d=3..7`), inside the `1,000 ns` cycle budget with `Λ_2D(5→7) = 5.80x` & `Λ_3D(5→7) = 5.71x`.
