namespace Qua.QeccCliffordT

/-!
  QUA Pure-Init Lean 4.36.0 Certificate — QECC Stabilizer/CWS & Exact Z[ζ_8, 1/2] Clifford+T
  (`zero Mathlib`, `zero sorry`, `zero axiom`, `zero native_decide`, `zero set_option`)
  =========================================================================================
  1. Exact QECC Stabilizer Code Distance Witnesses (`[[5,1,3]]`, `[[7,1,3]]`, `[[8,3,3]]`, `[[9,1,3]]`,
     `[[10,2,3]]`, `[[10,2,4]]`, `[[11,1,5]]`, `[[12,0,6]]`) with pairwise commutation, full `F_2` rank,
     exhaustive centralizer / weight-bounded checks, and explicit logical operator witnesses.
  2. Non-Additive CWS Codes Beating Stabilizer Bounds (`((5,6,2))`, `((9,12,3))`, `((10,24,3))`) with universal
     weight-1 and weight-2 soundness theorems.
  3. Exact Cyclotomic Ring `Z[ζ_8, 1/2]` (`ζ_8^4 + 1 = 0`) Clifford+T Identities (`CS` 3-T, `CCZ` 7-T,
     `Toffoli` 7-T, `Fredkin` 7-T, `RelativePhaseToffoli` 4-T), exhaustive `minT(CS) = 3` proof across
     241 Pauli rotations, and universal algebraic `minT(CCZ_phase_poly) = 7` theorem (`ccz_phase_poly_min7`).
-/

/-! ## Common Bounded Quantifiers -/

def allPow2 : Nat → (Nat → Bool) → Bool
  | 0, p => p 0
  | k + 1, p => allPow2 k p && allPow2 k (fun i => p (i + 2 ^ k))

theorem allPow2_sound : ∀ k (p : Nat → Bool), allPow2 k p = true → ∀ i, i < 2 ^ k → p i = true
  | 0, p, h, i, hi => by
    have : i = 0 := by omega
    subst this; exact h
  | k + 1, p, h, i, hi => by
    unfold allPow2 at h
    simp only [Bool.and_eq_true] at h
    by_cases hlt : i < 2 ^ k
    · exact allPow2_sound k p h.1 i hlt
    · have hi2 : i - 2 ^ k < 2 ^ k := by
        have hpow : 2 ^ (k + 1) = 2 ^ k + 2 ^ k := by rw [Nat.pow_succ]; omega
        omega
      have hr := allPow2_sound k (fun j => p (j + 2 ^ k)) h.2 (i - 2 ^ k) hi2
      have heq : (i - 2 ^ k) + 2 ^ k = i := by omega
      rwa [heq] at hr

def allBelow : Nat → (Nat → Bool) → Bool
  | 0, _ => true
  | n + 1, p => allBelow n p && p n

theorem allBelow_sound : ∀ n (p : Nat → Bool), allBelow n p = true → ∀ i, i < n → p i = true
  | 0, _, _, i, hi => by omega
  | n + 1, p, h, i, hi => by
    unfold allBelow at h
    simp only [Bool.and_eq_true] at h
    by_cases heq : i = n
    · subst heq; exact h.2
    · have hlt : i < n := by omega
      exact allBelow_sound n p h.1 i hlt

/-! ## Section 1: Exact Stabilizer Codes over `F_2^{2n}` -/

private def pop4 (x : Nat) : Nat :=
  match x with
  | 0 => 0 | 1 => 1 | 2 => 1 | 3 => 2
  | 4 => 1 | 5 => 2 | 6 => 2 | 7 => 3
  | 8 => 1 | 9 => 2 | 10 => 2 | 11 => 3
  | 12 => 2 | 13 => 3 | 14 => 3 | _ => 4

/-- Fast O(1) popcount for up to 16-bit integers via 4 nibble lookups. -/
def popcnt16 (x : Nat) : Nat :=
  pop4 (x % 16) + pop4 ((x / 16) % 16) + pop4 ((x / 256) % 16) + pop4 ((x / 4096) % 16)

def pauliWt (x z : Nat) : Nat :=
  popcnt16 (x ||| z)

def sympInner (x1 z1 x2 z2 : Nat) : Nat :=
  (popcnt16 ((x1 &&& z2) ^^^ (z1 &&& x2))) % 2

def commutesAll (gens : List (Nat × Nat)) (u v : Nat) : Bool :=
  gens.all fun g => sympInner g.1 g.2 u v == 0

def pairwiseCommute (gens : List (Nat × Nat)) : Bool :=
  gens.all fun g1 => gens.all fun g2 => sympInner g1.1 g1.2 g2.1 g2.2 == 0

def inStabRref (n : Nat) (rref : List (Nat × Nat)) (pivs : List Nat) (u v : Nat) : Bool :=
  let rec go (rs : List (Nat × Nat)) (ps : List Nat) (vec : Nat) : Nat :=
    match rs, ps with
    | r :: rs', p :: ps' =>
      let vec' := if ((vec >>> p) % 2) == 1 then vec ^^^ (r.1 ||| (r.2 <<< n)) else vec
      go rs' ps' vec'
    | _, _ => vec
  go rref pivs (u ||| (v <<< n)) == 0

def linCombPairs (basis : List (Nat × Nat)) (mask : Nat) : Nat × Nat :=
  let rec go (bs : List (Nat × Nat)) (m accX accZ : Nat) : Nat × Nat :=
    match bs with
    | [] => (accX, accZ)
    | b :: bs' =>
      if m % 2 == 1 then go bs' (m / 2) (accX ^^^ b.1) (accZ ^^^ b.2)
      else go bs' (m / 2) accX accZ
  go basis mask 0 0

/-- Centralizer verifier: checks every basis vector commutes with `gens`, and every linear
    combination `idx < 2^dimPerp` not in `S` has weight ≥ `d`. -/
def verifyCentralizerPow2 (n dimPerp d : Nat) (gens perpBasis sRref : List (Nat × Nat)) (sPivs : List Nat) : Bool :=
  (perpBasis.all fun b => commutesAll gens b.1 b.2) &&
  allPow2 dimPerp fun idx =>
    let v := linCombPairs perpBasis idx
    inStabRref n sRref sPivs v.1 v.2 || decide (pauliWt v.1 v.2 ≥ d)

def f2Rank (cols : Nat) (vecs : List Nat) : Nat :=
  let rec elimCol (c : Nat) (rows : List Nat) (rk : Nat) : List Nat × Nat :=
    match rows.findIdx? (fun r => ((r >>> c) % 2) == 1) with
    | none => (rows, rk)
    | some idx =>
      let piv := rows[idx]!
      let rest := (rows.eraseIdx idx).map (fun r => if ((r >>> c) % 2) == 1 then r ^^^ piv else r)
      (rest, rk + 1)
  let rec loop (c : Nat) (rows : List Nat) (rk : Nat) : Nat :=
    match c with
    | 0 => rk
    | c' + 1 =>
      let (rows', rk') := elimCol c' rows rk
      loop c' rows' rk'
  loop cols vecs 0

def wt1u (q op : Nat) : Nat := if op == 0 || op == 2 then 1 <<< q else 0
def wt1v (q op : Nat) : Nat := if op == 1 || op == 2 then 1 <<< q else 0
def wt2u (q1 q2 op1 op2 : Nat) : Nat := wt1u q1 op1 ||| wt1u q2 op2
def wt2v (q1 q2 op1 op2 : Nat) : Nat := wt1v q1 op1 ||| wt1v q2 op2

def stabCheckWtLt3 (n : Nat) (gens sRref : List (Nat × Nat)) (sPivs : List Nat) : Bool :=
  (allBelow n fun q => allBelow 3 fun op =>
    !(commutesAll gens (wt1u q op) (wt1v q op)) ||
    inStabRref n sRref sPivs (wt1u q op) (wt1v q op)) &&
  (allBelow n fun q1 => allBelow q1 fun q2 => allBelow 3 fun op1 => allBelow 3 fun op2 =>
    !(commutesAll gens (wt2u q1 q2 op1 op2) (wt2v q1 q2 op1 op2)) ||
    inStabRref n sRref sPivs (wt2u q1 q2 op1 op2) (wt2v q1 q2 op1 op2))

theorem stabCheckWtLt3_sound (n : Nat) (gens sRref : List (Nat × Nat)) (sPivs : List Nat)
    (h : stabCheckWtLt3 n gens sRref sPivs = true) :
    (∀ q op, q < n → op < 3 →
      commutesAll gens (wt1u q op) (wt1v q op) = true →
      inStabRref n sRref sPivs (wt1u q op) (wt1v q op) = true) ∧
    (∀ q1 q2 op1 op2, q1 < n → q2 < q1 → op1 < 3 → op2 < 3 →
      commutesAll gens (wt2u q1 q2 op1 op2) (wt2v q1 q2 op1 op2) = true →
      inStabRref n sRref sPivs (wt2u q1 q2 op1 op2) (wt2v q1 q2 op1 op2) = true) := by
  unfold stabCheckWtLt3 at h
  simp only [Bool.and_eq_true] at h
  refine ⟨?_, ?_⟩
  · intro q op hq hop hc
    have hq1 := allBelow_sound n _ h.1 q hq
    have hop1 := allBelow_sound 3 _ hq1 op hop
    simp only [Bool.or_eq_true, Bool.not_eq_true'] at hop1
    rcases hop1 with hnc | hins
    · rw [hc] at hnc; contradiction
    · exact hins
  · intro q1 q2 op1 op2 hq1 hq2 hop1 hop2 hc
    have a1 := allBelow_sound n _ h.2 q1 hq1
    have a2 := allBelow_sound q1 _ a1 q2 hq2
    have a3 := allBelow_sound 3 _ a2 op1 hop1
    have a4 := allBelow_sound 3 _ a3 op2 hop2
    simp only [Bool.or_eq_true, Bool.not_eq_true'] at a4
    rcases a4 with hnc | hins
    · rw [hc] at hnc; contradiction
    · exact hins

-- 5-qubit [[5,1,3]] [[5,1,3]]
def stab513_gens : List (Nat × Nat) := [(9, 6), (18, 12), (5, 24), (10, 17)]
def stab513_perp : List (Nat × Nat) := [(31, 0), (12, 1), (7, 2), (14, 4), (3, 8), (6, 16)]
def stab513_srref : List (Nat × Nat) := [(17, 27), (18, 12), (20, 3), (24, 29)]
def stab513_spivs : List Nat := [0, 1, 2, 3]
theorem stab513_valid :
    pairwiseCommute stab513_gens = true ∧
    f2Rank 10 (stab513_gens.map (fun g => g.1 ||| (g.2 <<< 5))) = 4 ∧
    f2Rank 10 (stab513_perp.map (fun g => g.1 ||| (g.2 <<< 5))) = 6 := by
  refine ⟨by decide, by decide, by decide⟩

theorem stab513_centralizer_dist :
    verifyCentralizerPow2 5 6 3 stab513_gens stab513_perp stab513_srref stab513_spivs = true ∧
    (commutesAll stab513_gens 19 1 = true ∧
     inStabRref 5 stab513_srref stab513_spivs 19 1 = false ∧
     pauliWt 19 1 = 3) := by
  refine ⟨by decide, by decide, by decide, by decide⟩

theorem stab513_wt_lt_3_detected :
    (∀ q op, q < 5 → op < 3 → commutesAll stab513_gens (wt1u q op) (wt1v q op) = true → inStabRref 5 stab513_srref stab513_spivs (wt1u q op) (wt1v q op) = true) ∧
    (∀ q1 q2 op1 op2, q1 < 5 → q2 < q1 → op1 < 3 → op2 < 3 → commutesAll stab513_gens (wt2u q1 q2 op1 op2) (wt2v q1 q2 op1 op2) = true → inStabRref 5 stab513_srref stab513_spivs (wt2u q1 q2 op1 op2) (wt2v q1 q2 op1 op2) = true) :=
  stabCheckWtLt3_sound 5 stab513_gens stab513_srref stab513_spivs (by decide)
-- Steane [[7,1,3]] [[7,1,3]]
def stab713_gens : List (Nat × Nat) := [(120, 0), (102, 0), (85, 0), (0, 120), (0, 102), (0, 85)]
def stab713_perp : List (Nat × Nat) := [(7, 0), (25, 0), (42, 0), (75, 0), (0, 7), (0, 25), (0, 42), (0, 75)]
def stab713_srref : List (Nat × Nat) := [(85, 0), (102, 0), (120, 0), (0, 85), (0, 102), (0, 120)]
def stab713_spivs : List Nat := [0, 1, 3, 7, 8, 10]
theorem stab713_valid :
    pairwiseCommute stab713_gens = true ∧
    f2Rank 14 (stab713_gens.map (fun g => g.1 ||| (g.2 <<< 7))) = 6 ∧
    f2Rank 14 (stab713_perp.map (fun g => g.1 ||| (g.2 <<< 7))) = 8 := by
  refine ⟨by decide, by decide, by decide⟩

theorem stab713_centralizer_dist :
    verifyCentralizerPow2 7 8 3 stab713_gens stab713_perp stab713_srref stab713_spivs = true ∧
    (commutesAll stab713_gens 7 0 = true ∧
     inStabRref 7 stab713_srref stab713_spivs 7 0 = false ∧
     pauliWt 7 0 = 3) := by
  refine ⟨by decide, by decide, by decide, by decide⟩

theorem stab713_wt_lt_3_detected :
    (∀ q op, q < 7 → op < 3 → commutesAll stab713_gens (wt1u q op) (wt1v q op) = true → inStabRref 7 stab713_srref stab713_spivs (wt1u q op) (wt1v q op) = true) ∧
    (∀ q1 q2 op1 op2, q1 < 7 → q2 < q1 → op1 < 3 → op2 < 3 → commutesAll stab713_gens (wt2u q1 q2 op1 op2) (wt2v q1 q2 op1 op2) = true → inStabRref 7 stab713_srref stab713_spivs (wt2u q1 q2 op1 op2) (wt2v q1 q2 op1 op2) = true) :=
  stabCheckWtLt3_sound 7 stab713_gens stab713_srref stab713_spivs (by decide)
-- 8-3-3 [[8,3,3]] [[8,3,3]]
def stab833_gens : List (Nat × Nat) := [(15, 153), (17, 170), (34, 85), (71, 45), (132, 75)]
def stab833_perp : List (Nat × Nat) := [(30, 0), (45, 0), (75, 0), (135, 0), (9, 3), (2, 5), (13, 9), (10, 17), (3, 33), (7, 65), (8, 129)]
def stab833_srref : List (Nat × Nat) := [(225, 51), (34, 85), (132, 75), (72, 180), (240, 153)]
def stab833_spivs : List Nat := [0, 1, 2, 3, 4]
theorem stab833_valid :
    pairwiseCommute stab833_gens = true ∧
    f2Rank 16 (stab833_gens.map (fun g => g.1 ||| (g.2 <<< 8))) = 5 ∧
    f2Rank 16 (stab833_perp.map (fun g => g.1 ||| (g.2 <<< 8))) = 11 := by
  refine ⟨by decide, by decide, by decide⟩

theorem stab833_centralizer_dist :
    verifyCentralizerPow2 8 11 3 stab833_gens stab833_perp stab833_srref stab833_spivs = true ∧
    (commutesAll stab833_gens 66 3 = true ∧
     inStabRref 8 stab833_srref stab833_spivs 66 3 = false ∧
     pauliWt 66 3 = 3) := by
  refine ⟨by decide, by decide, by decide, by decide⟩

theorem stab833_wt_lt_3_detected :
    (∀ q op, q < 8 → op < 3 → commutesAll stab833_gens (wt1u q op) (wt1v q op) = true → inStabRref 8 stab833_srref stab833_spivs (wt1u q op) (wt1v q op) = true) ∧
    (∀ q1 q2 op1 op2, q1 < 8 → q2 < q1 → op1 < 3 → op2 < 3 → commutesAll stab833_gens (wt2u q1 q2 op1 op2) (wt2v q1 q2 op1 op2) = true → inStabRref 8 stab833_srref stab833_spivs (wt2u q1 q2 op1 op2) (wt2v q1 q2 op1 op2) = true) :=
  stabCheckWtLt3_sound 8 stab833_gens stab833_srref stab833_spivs (by decide)
-- Shor [[9,1,3]] [[9,1,3]]
def stab913_gens : List (Nat × Nat) := [(0, 3), (0, 6), (0, 24), (0, 48), (0, 192), (0, 384), (63, 0), (504, 0)]
def stab913_perp : List (Nat × Nat) := [(7, 0), (56, 0), (448, 0), (0, 3), (0, 5), (0, 24), (0, 40), (0, 73), (0, 137), (0, 265)]
def stab913_srref : List (Nat × Nat) := [(455, 0), (504, 0), (0, 5), (0, 6), (0, 40), (0, 48), (0, 320), (0, 384)]
def stab913_spivs : List Nat := [0, 3, 9, 10, 12, 13, 15, 16]
theorem stab913_valid :
    pairwiseCommute stab913_gens = true ∧
    f2Rank 18 (stab913_gens.map (fun g => g.1 ||| (g.2 <<< 9))) = 8 ∧
    f2Rank 18 (stab913_perp.map (fun g => g.1 ||| (g.2 <<< 9))) = 10 := by
  refine ⟨by decide, by decide, by decide⟩

theorem stab913_centralizer_dist :
    verifyCentralizerPow2 9 10 3 stab913_gens stab913_perp stab913_srref stab913_spivs = true ∧
    (commutesAll stab913_gens 7 0 = true ∧
     inStabRref 9 stab913_srref stab913_spivs 7 0 = false ∧
     pauliWt 7 0 = 3) := by
  refine ⟨by decide, by decide, by decide, by decide⟩

theorem stab913_wt_lt_3_detected :
    (∀ q op, q < 9 → op < 3 → commutesAll stab913_gens (wt1u q op) (wt1v q op) = true → inStabRref 9 stab913_srref stab913_spivs (wt1u q op) (wt1v q op) = true) ∧
    (∀ q1 q2 op1 op2, q1 < 9 → q2 < q1 → op1 < 3 → op2 < 3 → commutesAll stab913_gens (wt2u q1 q2 op1 op2) (wt2v q1 q2 op1 op2) = true → inStabRref 9 stab913_srref stab913_spivs (wt2u q1 q2 op1 op2) (wt2v q1 q2 op1 op2) = true) :=
  stabCheckWtLt3_sound 9 stab913_gens stab913_srref stab913_spivs (by decide)
-- 10-2-3 [[10,2,3]] [[10,2,3]]
def stab1023_gens : List (Nat × Nat) := [(2, 554), (554, 552), (4, 85), (85, 81), (8, 170), (170, 162), (16, 340), (340, 324)]
def stab1023_perp : List (Nat × Nat) := [(321, 0), (672, 0), (68, 1), (136, 2), (81, 4), (162, 8), (5, 16), (10, 32), (20, 64), (40, 128), (80, 256), (34, 512)]
def stab1023_srref : List (Nat × Nat) := [(257, 277), (2, 554), (4, 85), (8, 170), (16, 340), (544, 168), (320, 69), (640, 650)]
def stab1023_spivs : List Nat := [0, 1, 2, 3, 4, 5, 6, 7]
theorem stab1023_valid :
    pairwiseCommute stab1023_gens = true ∧
    f2Rank 20 (stab1023_gens.map (fun g => g.1 ||| (g.2 <<< 10))) = 8 ∧
    f2Rank 20 (stab1023_perp.map (fun g => g.1 ||| (g.2 <<< 10))) = 12 := by
  refine ⟨by decide, by decide, by decide⟩

theorem stab1023_centralizer_dist :
    verifyCentralizerPow2 10 12 3 stab1023_gens stab1023_perp stab1023_srref stab1023_spivs = true ∧
    (commutesAll stab1023_gens 321 0 = true ∧
     inStabRref 10 stab1023_srref stab1023_spivs 321 0 = false ∧
     pauliWt 321 0 = 3) := by
  refine ⟨by decide, by decide, by decide, by decide⟩

theorem stab1023_wt_lt_3_detected :
    (∀ q op, q < 10 → op < 3 → commutesAll stab1023_gens (wt1u q op) (wt1v q op) = true → inStabRref 10 stab1023_srref stab1023_spivs (wt1u q op) (wt1v q op) = true) ∧
    (∀ q1 q2 op1 op2, q1 < 10 → q2 < q1 → op1 < 3 → op2 < 3 → commutesAll stab1023_gens (wt2u q1 q2 op1 op2) (wt2v q1 q2 op1 op2) = true → inStabRref 10 stab1023_srref stab1023_spivs (wt2u q1 q2 op1 op2) (wt2v q1 q2 op1 op2) = true) :=
  stabCheckWtLt3_sound 10 stab1023_gens stab1023_srref stab1023_spivs (by decide)
-- 10-2-4 [[10,2,4]] [[10,2,4]]
def stab1024_gens : List (Nat × Nat) := [(1023, 951), (445, 640), (681, 1017), (785, 540), (926, 609), (640, 667), (835, 216), (278, 679)]
def stab1024_perp : List (Nat × Nat) := [(341, 0), (682, 0), (189, 1), (46, 2), (225, 4), (151, 8), (198, 16), (217, 32), (91, 64), (182, 128), (133, 256), (94, 512)]
def stab1024_srref : List (Nat × Nat) := [(257, 495), (514, 990), (260, 650), (8, 605), (528, 1011), (288, 720), (64, 745), (640, 667)]
def stab1024_spivs : List Nat := [0, 1, 2, 3, 4, 5, 6, 7]
theorem stab1024_valid :
    pairwiseCommute stab1024_gens = true ∧
    f2Rank 20 (stab1024_gens.map (fun g => g.1 ||| (g.2 <<< 10))) = 8 ∧
    f2Rank 20 (stab1024_perp.map (fun g => g.1 ||| (g.2 <<< 10))) = 12 := by
  refine ⟨by decide, by decide, by decide⟩

theorem stab1024_centralizer_dist :
    verifyCentralizerPow2 10 12 4 stab1024_gens stab1024_perp stab1024_srref stab1024_spivs = true ∧
    (commutesAll stab1024_gens 147 3 = true ∧
     inStabRref 10 stab1024_srref stab1024_spivs 147 3 = false ∧
     pauliWt 147 3 = 4) := by
  refine ⟨by decide, by decide, by decide, by decide⟩
-- 11-1-5 [[11,1,5]] [[11,1,5]]
def stab1115_gens : List (Nat × Nat) := [(1603, 1518), (556, 1358), (1468, 1625), (913, 1998), (1258, 1559), (1918, 700), (434, 1839), (1219, 1909), (367, 1430), (278, 679)]
def stab1115_perp : List (Nat × Nat) := [(682, 0), (1365, 0), (379, 2), (92, 5), (450, 8), (302, 17), (396, 32), (435, 65), (182, 128), (365, 257), (267, 512), (189, 1025)]
def stab1115_srref : List (Nat × Nat) := [(1025, 1241), (2, 1175), (1028, 1980), (8, 605), (16, 1210), (1056, 2022), (64, 745), (128, 1490), (1280, 1334), (512, 1865)]
def stab1115_spivs : List Nat := [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
theorem stab1115_valid :
    pairwiseCommute stab1115_gens = true ∧
    f2Rank 22 (stab1115_gens.map (fun g => g.1 ||| (g.2 <<< 11))) = 10 ∧
    f2Rank 22 (stab1115_perp.map (fun g => g.1 ||| (g.2 <<< 11))) = 12 := by
  refine ⟨by decide, by decide, by decide⟩

theorem stab1115_centralizer_dist :
    verifyCentralizerPow2 11 12 5 stab1115_gens stab1115_perp stab1115_srref stab1115_spivs = true ∧
    (commutesAll stab1115_gens 682 0 = true ∧
     inStabRref 11 stab1115_srref stab1115_spivs 682 0 = false ∧
     pauliWt 682 0 = 5) := by
  refine ⟨by decide, by decide, by decide, by decide⟩
-- Dodeca [[12,0,6]] self-dual stabilizer state
def dodeca12_gens : List (Nat × Nat) := [(2849, 759), (1603, 1518), (3206, 3036), (2317, 1977), (539, 3954), (1078, 3813), (2156, 3531), (217, 2967), (434, 1839), (868, 3678), (1736, 3261), (3472, 2427)]
theorem dodeca12_valid_and_dist6 :
    pairwiseCommute dodeca12_gens = true ∧
    f2Rank 24 (dodeca12_gens.map (fun g => g.1 ||| (g.2 <<< 12))) = 12 ∧
    allPow2 12 (fun idx => idx == 0 || decide (pauliWt (linCombPairs dodeca12_gens idx).1 (linCombPairs dodeca12_gens idx).2 ≥ 6)) = true ∧
    (commutesAll dodeca12_gens 1025 1241 = true ∧ pauliWt 1025 1241 = 6) := by
  refine ⟨by decide, by decide, by decide, by decide, by decide⟩

/-! ## Section 2: Non-Additive CWS Codes Beating Stabilizer Bounds -/

def cwsGamma (adj : List Nat) (u : Nat) : Nat :=
  let rec go (i : Nat) (xs : List Nat) (acc : Nat) : Nat :=
    match xs with
    | [] => acc
    | row :: xs' =>
      let acc' := if ((u >>> i) % 2) == 1 then acc ^^^ row else acc
      go (i + 1) xs' acc'
  go 0 adj 0

def cwsDetects (adj cw : List Nat) (u v : Nat) : Bool :=
  let cl := v ^^^ cwsGamma adj u
  if cl == 0 then
    match cw with
    | [] => true
    | c0 :: cs =>
      let p0 := popcnt16 (u &&& c0) % 2
      cs.all fun c => popcnt16 (u &&& c) % 2 == p0
  else
    cw.all fun c1 => cw.all fun c2 => (c1 == c2) || ((c1 ^^^ c2) != cl)

def cwsDetectsFast (adj cw : List Nat) (inCw : Nat → Bool) (u v : Nat) : Bool :=
  let cl := v ^^^ cwsGamma adj u
  if cl == 0 then
    match cw with
    | [] => true
    | c0 :: cs =>
      let p0 := popcnt16 (u &&& c0) % 2
      cs.all fun c => popcnt16 (u &&& c) % 2 == p0
  else
    cw.all fun c1 => !(inCw (c1 ^^^ cl))

theorem cwsDetectsFast_sound (adj cw : List Nat) (inCw : Nat → Bool)
    (hIn : ∀ c, c ∈ cw → inCw c = true) (u v : Nat)
    (h : cwsDetectsFast adj cw inCw u v = true) :
    cwsDetects adj cw u v = true := by
  unfold cwsDetectsFast at h
  unfold cwsDetects
  by_cases hcl : (v ^^^ cwsGamma adj u) == 0
  · simp only [hcl, ↓reduceIte] at h ⊢; exact h
  · simp only [hcl, Bool.false_eq_true, ↓reduceIte] at h ⊢
    rw [List.all_eq_true] at h ⊢
    intro c1 hc1
    rw [List.all_eq_true]
    intro c2 hc2
    by_cases h12 : c1 == c2
    · simp [h12]
    · simp only [h12, Bool.false_or, bne_iff_ne, ne_eq]
      intro hxor
      have hfast := h c1 hc1
      have hc2' : c1 ^^^ (v ^^^ cwsGamma adj u) = c2 := by
        rw [← hxor, ← Nat.xor_assoc, Nat.xor_self, Nat.zero_xor]
      rw [hc2', hIn c2 hc2] at hfast
      contradiction

def cwsCheckWt1Fast (n : Nat) (adj cw : List Nat) (inCw : Nat → Bool) : Bool :=
  allBelow n fun q => allBelow 3 fun op => cwsDetectsFast adj cw inCw (wt1u q op) (wt1v q op)

def cwsCheckWt2Fast (n : Nat) (adj cw : List Nat) (inCw : Nat → Bool) : Bool :=
  allBelow n fun q1 => allBelow q1 fun q2 => allBelow 3 fun op1 => allBelow 3 fun op2 =>
    cwsDetectsFast adj cw inCw (wt2u q1 q2 op1 op2) (wt2v q1 q2 op1 op2)

theorem cwsCheckWt1Fast_sound (n : Nat) (adj cw : List Nat) (inCw : Nat → Bool)
    (hIn : ∀ c, c ∈ cw → inCw c = true) (h : cwsCheckWt1Fast n adj cw inCw = true) :
    ∀ q op, q < n → op < 3 → cwsDetects adj cw (wt1u q op) (wt1v q op) = true := by
  intro q op hq hop
  exact cwsDetectsFast_sound adj cw inCw hIn _ _ (allBelow_sound 3 _ (allBelow_sound n _ h q hq) op hop)

theorem cwsCheckWt2Fast_sound (n : Nat) (adj cw : List Nat) (inCw : Nat → Bool)
    (hIn : ∀ c, c ∈ cw → inCw c = true) (h : cwsCheckWt2Fast n adj cw inCw = true) :
    ∀ q1 q2 op1 op2, q1 < n → q2 < q1 → op1 < 3 → op2 < 3 →
      cwsDetects adj cw (wt2u q1 q2 op1 op2) (wt2v q1 q2 op1 op2) = true := by
  intro q1 q2 op1 op2 hq1 hq2 hop1 hop2
  exact cwsDetectsFast_sound adj cw inCw hIn _ _
    (allBelow_sound 3 _ (allBelow_sound 3 _ (allBelow_sound q1 _ (allBelow_sound n _ h q1 hq1) q2 hq2) op1 hop1) op2 hop2)

def cws562_adj : List Nat := [18, 5, 10, 20, 9]
def cws562_cw : List Nat := [0, 3, 12, 22, 27, 29]
def cws562_inCw (x : Nat) : Bool := match x with | 0 => true | 3 => true | 12 => true | 22 => true | 27 => true | 29 => true | _ => false
private theorem cws562_inCw_sound : ∀ c, c ∈ cws562_cw → cws562_inCw c = true := by decide
theorem cws562_beats_stabilizer_and_exact_d2 :
    cws562_cw.length = 6 ∧ cws562_cw.length > 4 ∧
    (∀ q op, q < 5 → op < 3 → cwsDetects cws562_adj cws562_cw (wt1u q op) (wt1v q op) = true) ∧
    (cwsDetects cws562_adj cws562_cw 3 0 = false ∧ pauliWt 3 0 = 2) := by
  refine ⟨rfl, by decide, cwsCheckWt1Fast_sound 5 cws562_adj cws562_cw cws562_inCw cws562_inCw_sound (by decide), by decide, by decide⟩

def cws9123_adj : List Nat := [258, 5, 10, 20, 40, 80, 160, 320, 129]
def cws9123_cw : List Nat := [0, 73, 35, 273, 439, 230, 175, 412, 469, 106, 344, 510]
def cws9123_inCw (x : Nat) : Bool := match x with | 0 => true | 73 => true | 35 => true | 273 => true | 439 => true | 230 => true | 175 => true | 412 => true | 469 => true | 106 => true | 344 => true | 510 => true | _ => false
private theorem cws9123_inCw_sound : ∀ c, c ∈ cws9123_cw → cws9123_inCw c = true := by decide
theorem cws9123_beats_stabilizer_and_exact_d3 :
    cws9123_cw.length = 12 ∧ cws9123_cw.length > 8 ∧
    (∀ q op, q < 9 → op < 3 → cwsDetects cws9123_adj cws9123_cw (wt1u q op) (wt1v q op) = true) ∧
    (∀ q1 q2 op1 op2, q1 < 9 → q2 < q1 → op1 < 3 → op2 < 3 → cwsDetects cws9123_adj cws9123_cw (wt2u q1 q2 op1 op2) (wt2v q1 q2 op1 op2) = true) ∧
    (cwsDetects cws9123_adj cws9123_cw 2 5 = false ∧ pauliWt 2 5 = 3) := by
  refine ⟨rfl, by decide, cwsCheckWt1Fast_sound 9 cws9123_adj cws9123_cw cws9123_inCw cws9123_inCw_sound (by decide), cwsCheckWt2Fast_sound 9 cws9123_adj cws9123_cw cws9123_inCw cws9123_inCw_sound (by decide), by decide, by decide⟩

def cws10243_adj : List Nat := [434, 837, 682, 372, 713, 973, 954, 885, 747, 502]
def cws10243_cw : List Nat := [0, 298, 829, 463, 843, 365, 992, 596, 169, 338, 730, 437, 31, 645, 196, 392, 545, 926, 98, 784, 694, 635, 1023, 247]
def cws10243_inCw (x : Nat) : Bool := match x with | 0 => true | 298 => true | 829 => true | 463 => true | 843 => true | 365 => true | 992 => true | 596 => true | 169 => true | 338 => true | 730 => true | 437 => true | 31 => true | 645 => true | 196 => true | 392 => true | 545 => true | 926 => true | 98 => true | 784 => true | 694 => true | 635 => true | 1023 => true | 247 => true | _ => false
private theorem cws10243_inCw_sound : ∀ c, c ∈ cws10243_cw → cws10243_inCw c = true := by decide
theorem cws10243_beats_stabilizer_and_exact_d3 :
    cws10243_cw.length = 24 ∧ cws10243_cw.length > 16 ∧
    (∀ q op, q < 10 → op < 3 → cwsDetects cws10243_adj cws10243_cw (wt1u q op) (wt1v q op) = true) ∧
    (∀ q1 q2 op1 op2, q1 < 10 → q2 < q1 → op1 < 3 → op2 < 3 → cwsDetects cws10243_adj cws10243_cw (wt2u q1 q2 op1 op2) (wt2v q1 q2 op1 op2) = true) ∧
    (cwsDetects cws10243_adj cws10243_cw 6 1 = false ∧ pauliWt 6 1 = 3) := by
  refine ⟨rfl, by decide, cwsCheckWt1Fast_sound 10 cws10243_adj cws10243_cw cws10243_inCw cws10243_inCw_sound (by decide), cwsCheckWt2Fast_sound 10 cws10243_adj cws10243_cw cws10243_inCw cws10243_inCw_sound (by decide), by decide, by decide⟩

/-! ## Section 3: Exact Cyclotomic Ring `Z[ζ_8, 1/2]` Clifford+T Identities & Min-T Bounds -/

structure Z8Int where
  a : Int
  b : Int
  c : Int
  d : Int
  deriving DecidableEq, Repr

def z8Zero : Z8Int := ⟨0, 0, 0, 0⟩
def z8One  : Z8Int := ⟨1, 0, 0, 0⟩
def z8MOne : Z8Int := ⟨-1, 0, 0, 0⟩
def z8W    : Z8Int := ⟨0, 1, 0, 0⟩
def z8I    : Z8Int := ⟨0, 0, 1, 0⟩
def z8W3   : Z8Int := ⟨0, 0, 0, 1⟩
def z8Sqrt2 : Z8Int := ⟨0, 1, 0, -1⟩

def z8Add (x y : Z8Int) : Z8Int := ⟨x.a + y.a, x.b + y.b, x.c + y.c, x.d + y.d⟩
def z8Sub (x y : Z8Int) : Z8Int := ⟨x.a - y.a, x.b - y.b, x.c - y.c, x.d - y.d⟩
def z8Neg (x : Z8Int) : Z8Int := ⟨-x.a, -x.b, -x.c, -x.d⟩
def z8Conj (x : Z8Int) : Z8Int := ⟨x.a, -x.d, -x.c, -x.b⟩
def z8Mul (x y : Z8Int) : Z8Int :=
  ⟨x.a * y.a - x.b * y.d - x.c * y.c - x.d * y.b,
   x.a * y.b + x.b * y.a - x.c * y.d - x.d * y.c,
   x.a * y.c + x.b * y.b + x.c * y.a - x.d * y.d,
   x.a * y.d + x.b * y.c + x.c * y.b + x.d * y.a⟩
def z8ScaleInt (s : Int) (x : Z8Int) : Z8Int := ⟨s * x.a, s * x.b, s * x.c, s * x.d⟩

theorem z8_cyclotomic_eq :
    z8Add (z8Mul (z8Mul z8W z8W) (z8Mul z8W z8W)) z8One = z8Zero := by
  decide

theorem z8_sqrt2_sq : z8Mul z8Sqrt2 z8Sqrt2 = ⟨2, 0, 0, 0⟩ := by
  decide

theorem z8_conj_mul_self_real_subring (x : Z8Int) :
    (z8Mul (z8Conj x) x).c = 0 ∧ (z8Mul (z8Conj x) x).b + (z8Mul (z8Conj x) x).d = 0 := by
  rcases x with ⟨a, b, c, d⟩
  simp only [z8Mul, z8Conj]
  constructor
  · have h1 : a * c = c * a := Int.mul_comm a c
    have h2 : (-d) * b = -(b * d) := by rw [Int.neg_mul, Int.mul_comm d b]
    have h3 : (-c) * a = -(c * a) := Int.neg_mul c a
    have h4 : (-b) * d = -(b * d) := Int.neg_mul b d
    rw [h1, h2, h3, h4]; omega
  · have h1 : a * b = b * a := Int.mul_comm a b
    have h2 : (-d) * a = -(a * d) := by rw [Int.neg_mul, Int.mul_comm d a]
    have h3 : (-c) * d = -(c * d) := Int.neg_mul c d
    have h4 : (-b) * c = -(b * c) := Int.neg_mul b c
    have h5 : a * d = d * a := Int.mul_comm a d
    have h6 : (-d) * c = -(c * d) := by rw [Int.neg_mul, Int.mul_comm d c]
    have h7 : (-c) * b = -(b * c) := by rw [Int.neg_mul, Int.mul_comm c b]
    have h8 : (-b) * a = -(b * a) := Int.neg_mul b a
    rw [h1, h2, h3, h4, h5, h6, h7, h8]; omega

def allFinN (n : Nat) (p : Fin n → Bool) : Bool :=
  let rec go (k : Nat) (hk : k ≤ n) : Bool :=
    match k with
    | 0 => true
    | k' + 1 => go k' (by omega) && p ⟨k', by omega⟩
  go n (Nat.le_refl n)

theorem allFinN_sound (n : Nat) (p : Fin n → Bool) (h : allFinN n p = true) : ∀ i : Fin n, p i = true := by
  have hgo : ∀ k (hk : k ≤ n), allFinN.go n p k hk = true → ∀ i : Fin n, i.val < k → p i = true := by
    intro k hk hg
    induction k with
    | zero => intro i hi; omega
    | succ k' ih =>
      unfold allFinN.go at hg
      simp only [Bool.and_eq_true] at hg
      intro i hi
      by_cases heq : i.val = k'
      · have hfin : i = ⟨k', by omega⟩ := Fin.ext heq
        rw [hfin]; exact hg.2
      · exact ih (by omega) hg.1 i (by omega)
  intro i
  exact hgo n (Nat.le_refl n) h i i.isLt

def z8MatMul {n : Nat} (A B : Fin n → Fin n → Z8Int) (i j : Fin n) : Z8Int :=
  Fin.foldl n (fun acc k => z8Add acc (z8Mul (A i k) (B k j))) z8Zero

def z8MatDagMul {n : Nat} (A B : Fin n → Fin n → Z8Int) (i j : Fin n) : Z8Int :=
  Fin.foldl n (fun acc k => z8Add acc (z8Mul (z8Conj (A k i)) (B k j))) z8Zero

def z8MatEq {n : Nat} (A B : Fin n → Fin n → Z8Int) : Bool :=
  allFinN n fun i => allFinN n fun j => A i j == B i j

theorem z8MatEq_sound {n : Nat} {A B : Fin n → Fin n → Z8Int} (h : z8MatEq A B = true) :
    ∀ i j : Fin n, A i j = B i j := by
  intro i j
  exact beq_iff_eq.mp (allFinN_sound n _ (allFinN_sound n _ h i) j)

def z8IsUnitaryScaled {n : Nat} (scale : Int) (M : Fin n → Fin n → Z8Int) : Bool :=
  allFinN n fun i => allFinN n fun j =>
    z8MatDagMul M M i j == if i == j then ⟨scale, 0, 0, 0⟩ else z8Zero

-- 2-qubit gates & exhaustive minT(CS) = 3
def gI2 (r c : Fin 2) : Z8Int := if r == c then z8One else z8Zero
def gX2 (r c : Fin 2) : Z8Int := if r != c then z8One else z8Zero
def gZ2 (r c : Fin 2) : Z8Int := if r == c then (if r.val == 0 then z8One else z8MOne) else z8Zero
def gY2 (r c : Fin 2) : Z8Int :=
  match r.val, c.val with
  | 0, 1 => z8Neg z8I
  | 1, 0 => z8I
  | _, _ => z8Zero
def gT2 (r c : Fin 2) : Z8Int := if r == c then (if r.val == 0 then z8One else z8W) else z8Zero
def gTdag2 (r c : Fin 2) : Z8Int := if r == c then (if r.val == 0 then z8One else z8Conj z8W) else z8Zero

def gKron2 (A B : Fin 2 → Fin 2 → Z8Int) (r c : Fin 4) : Z8Int :=
  z8Mul (A ⟨r.val / 2, by omega⟩ ⟨c.val / 2, by omega⟩)
        (B ⟨r.val % 2, Nat.mod_lt _ (by decide)⟩ ⟨c.val % 2, Nat.mod_lt _ (by decide)⟩)

def gCNOT01 (r c : Fin 4) : Z8Int :=
  let target : Nat := match c.val with | 2 => 3 | 3 => 2 | x => x
  if r.val == target then z8One else z8Zero

def gCS_target (r c : Fin 4) : Z8Int :=
  if r == c then (if r.val == 3 then z8I else z8One) else z8Zero

def gCS_circ3T : Fin 4 → Fin 4 → Z8Int :=
  z8MatMul (gKron2 gT2 gT2) (z8MatMul gCNOT01 (z8MatMul (gKron2 gI2 gTdag2) gCNOT01))

theorem cs_3T_exact_identity_and_unitary :
    (∀ i j : Fin 4, gCS_circ3T i j = gCS_target i j) ∧
    z8IsUnitaryScaled 1 gCS_circ3T = true := by
  exact ⟨z8MatEq_sound (by decide), by decide⟩

def pauli1q (idx : Fin 4) : Fin 2 → Fin 2 → Z8Int :=
  match idx.val with
  | 0 => gI2 | 1 => gX2 | 2 => gY2 | _ => gZ2

def pauli2q15 (p : Fin 15) : Fin 4 → Fin 4 → Z8Int :=
  let idx := p.val + 1
  gKron2 (pauli1q ⟨idx / 4, by omega⟩) (pauli1q ⟨idx % 4, Nat.mod_lt _ (by decide)⟩)

def tRotDagScaled2 (p : Fin 15) (r c : Fin 4) : Z8Int :=
  let P := pauli2q15 p
  let d := if r == c then z8One else z8Zero
  z8Add (z8Mul (z8Add z8One (z8Conj z8W)) d) (z8Mul (z8Sub z8One (z8Conj z8W)) (P r c))

def matEqScaled (s : Int) (A B : Fin 4 → Fin 4 → Z8Int) : Bool :=
  allFinN 4 fun i => allFinN 4 fun j => A i j == z8ScaleInt s (B i j)

def isSignedPauliScaled (s : Int) (C : Fin 4 → Fin 4 → Z8Int) : Bool :=
  Fin.foldl 15 (fun acc p =>
    let P := pauli2q15 p
    acc || matEqScaled s C P || matEqScaled (-s) C P) false

-- Fast concrete 4x4 matrix type `Mat4` for 2-qubit Clifford+T exhaustive search
structure Vec4 (α : Type) where
  x0 : α
  x1 : α
  x2 : α
  x3 : α
  deriving DecidableEq, Repr

def vec4Get {α : Type} (v : Vec4 α) (i : Fin 4) : α :=
  match i.val with | 0 => v.x0 | 1 => v.x1 | 2 => v.x2 | _ => v.x3

def vec4Of {α : Type} (f : Fin 4 → α) : Vec4 α := ⟨f 0, f 1, f 2, f 3⟩
abbrev Mat4 := Vec4 (Vec4 Z8Int)
def mat4Get (M : Mat4) (r c : Fin 4) : Z8Int := vec4Get (vec4Get M r) c
def mat4Of (f : Fin 4 → Fin 4 → Z8Int) : Mat4 := vec4Of fun r => vec4Of fun c => f r c
def mat4Mul (A B : Mat4) : Mat4 :=
  let dot4 (r c : Vec4 Z8Int) : Z8Int :=
    z8Add (z8Add (z8Mul r.x0 c.x0) (z8Mul r.x1 c.x1)) (z8Add (z8Mul r.x2 c.x2) (z8Mul r.x3 c.x3))
  let col (j : Fin 4) : Vec4 Z8Int := ⟨mat4Get B 0 j, mat4Get B 1 j, mat4Get B 2 j, mat4Get B 3 j⟩
  let c0 := col 0; let c1 := col 1; let c2 := col 2; let c3 := col 3
  ⟨⟨dot4 A.x0 c0, dot4 A.x0 c1, dot4 A.x0 c2, dot4 A.x0 c3⟩,
   ⟨dot4 A.x1 c0, dot4 A.x1 c1, dot4 A.x1 c2, dot4 A.x1 c3⟩,
   ⟨dot4 A.x2 c0, dot4 A.x2 c1, dot4 A.x2 c2, dot4 A.x2 c3⟩,
   ⟨dot4 A.x3 c0, dot4 A.x3 c1, dot4 A.x3 c2, dot4 A.x3 c3⟩⟩

def mat4Dag (M : Mat4) : Mat4 := mat4Of fun i j => z8Conj (mat4Get M j i)

def gCS_mat4 : Mat4 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩⟩⟩
private theorem gCS_mat4_eq : gCS_mat4 = mat4Of gCS_target := by decide
def gens2q_mat4 : List Mat4 := [
  ⟨⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩,
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩⟩,
  ⟨⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩,
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩⟩
]
def p15_mat4 : List Mat4 := [
  ⟨⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩,
  ⟨⟨⟨0, 0, 0, 0⟩, ⟨0, 0, -1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, -1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩,
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩⟩,
  ⟨⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩,
  ⟨⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩,
  ⟨⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, -1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, -1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩,
  ⟨⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩,
  ⟨⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, -1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, -1, 0⟩⟩,
   ⟨⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩,
  ⟨⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, -1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, -1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩,
  ⟨⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩,
  ⟨⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, -1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩⟩,
   ⟨⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, -1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩,
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩⟩,
  ⟨⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩,
  ⟨⟨⟨0, 0, 0, 0⟩, ⟨0, 0, -1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, -1, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩,
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩⟩
]
def tRotDag_get (p : Fin 15) : Mat4 :=
  match p.val with
  | 0 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 1⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨1, 0, 0, 1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 1⟩, ⟨1, 0, 0, -1⟩⟩⟩
  | 1 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 1, -1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, -1, 1, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 1, -1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, -1, 1, 0⟩, ⟨1, 0, 0, -1⟩⟩⟩
  | 2 => ⟨⟨⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -2⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -2⟩⟩⟩
  | 3 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 1⟩⟩,
   ⟨⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩⟩⟩
  | 4 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 1⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩⟩⟩
  | 5 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, -1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, -1, 1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 1, -1, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, -1, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩⟩⟩
  | 6 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, -1⟩⟩,
   ⟨⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩⟩⟩
  | 7 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, -1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, -1, 0⟩⟩,
   ⟨⟨0, -1, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, -1, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩⟩⟩
  | 8 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, -1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 1, -1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, -1, 1, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, -1, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩⟩⟩
  | 9 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, -1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 1⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨-1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩⟩⟩
  | 10 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, -1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, -1, 1, 0⟩⟩,
   ⟨⟨0, -1, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 1, -1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩⟩⟩
  | 11 => ⟨⟨⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -2⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -2⟩⟩⟩
  | 12 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 1⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨-1, 0, 0, -1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, -1⟩, ⟨1, 0, 0, -1⟩⟩⟩
  | 13 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 1, -1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, -1, 1, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, -1, 1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, -1, 0⟩, ⟨1, 0, 0, -1⟩⟩⟩
  | _ => ⟨⟨⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -2⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -2⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩⟩⟩
private def csM1_get (p : Fin 15) : Mat4 :=
  match p.val with
  | 0 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 1⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨1, 0, 0, 1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, -1, 1, 0⟩, ⟨0, 1, 1, 0⟩⟩⟩
  | 1 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 1, -1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, -1, 1, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 1, -1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, -1⟩, ⟨0, 1, 1, 0⟩⟩⟩
  | 2 => ⟨⟨⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -2⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 2, 0, 0⟩⟩⟩
  | 3 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 1⟩⟩,
   ⟨⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, -1, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 1, 0⟩⟩⟩
  | 4 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 1⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, -1, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 1, 0⟩⟩⟩
  | 5 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, -1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, -1, 1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 1, -1, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨-1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 1, 0⟩⟩⟩
  | 6 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, -1⟩⟩,
   ⟨⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 1, -1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 1, 0⟩⟩⟩
  | 7 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, -1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, -1, 0⟩⟩,
   ⟨⟨0, -1, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 1, 0⟩⟩⟩
  | 8 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, -1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 1, -1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, -1, 1, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨-1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 1, 0⟩⟩⟩
  | 9 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, -1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 1⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 1, -1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 1, 0⟩⟩⟩
  | 10 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, -1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, -1, 1, 0⟩⟩,
   ⟨⟨0, -1, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 1, 0⟩⟩⟩
  | 11 => ⟨⟨⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -2⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 2, 0, 0⟩⟩⟩
  | 12 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨1, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 1⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨-1, 0, 0, -1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, -1, 0⟩, ⟨0, 1, 1, 0⟩⟩⟩
  | 13 => ⟨⟨⟨1, 0, 0, -1⟩, ⟨0, 1, -1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, -1, 1, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, -1⟩, ⟨0, -1, 1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 1⟩, ⟨0, 1, 1, 0⟩⟩⟩
  | _ => ⟨⟨⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -2⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -2⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 2, 0⟩⟩⟩
private theorem p15_mat4_eq : allFinN 15 (fun p => p15_mat4.get ⟨p.val, by simp [p15_mat4]⟩ == mat4Of (pauli2q15 p)) = true := by decide
private theorem tRotDag_get_eq_all : allFinN 15 (fun p => tRotDag_get p == mat4Of (tRotDagScaled2 p)) = true := by decide
private theorem csM1_get_eq_all : allFinN 15 (fun p => mat4Mul gCS_mat4 (tRotDag_get p) == csM1_get p) = true := by decide

def mat4EqScaled (s : Int) (A B : Mat4) : Bool :=
  allFinN 4 fun i => allFinN 4 fun j => mat4Get A i j == z8ScaleInt s (mat4Get B i j)

def z8IsPauliEntry (s2 : Int) (e : Z8Int) : Bool :=
  (e.b == 0 && e.d == 0) &&
  ((e.a == 0 && e.c == 0) ||
   (e.c == 0 && (e.a == s2 || e.a == -s2)) ||
   (e.a == 0 && (e.c == s2 || e.c == -s2)))

def mat4AllPauliEntries (s2 : Int) (C : Mat4) : Bool :=
  allFinN 4 fun i => allFinN 4 fun j => z8IsPauliEntry s2 (mat4Get C i j)

def isClifford2qFast (s : Int) (M : Mat4) : Bool :=
  let Md := mat4Dag M
  let s2 := s * s
  gens2q_mat4.all fun G =>
    let C := mat4Mul M (mat4Mul G Md)
    mat4AllPauliEntries s2 C && p15_mat4.any fun P => mat4EqScaled s2 C P || mat4EqScaled (-s2) C P

def isClifford2qScaled (s : Int) (M : Fin 4 → Fin 4 → Z8Int) : Bool :=
  isClifford2qFast s (mat4Of M)

private theorem isClifford2qFast_eq_scaled (s : Int) (M : Fin 4 → Fin 4 → Z8Int) :
    isClifford2qScaled s M = isClifford2qFast s (mat4Of M) := rfl

private def csCheck1 (p1 : Fin 15) : Bool :=
  !(isClifford2qFast 2 (csM1_get p1))

private def csCheck2 (p1 p2 : Fin 15) : Bool :=
  !(isClifford2qFast 4 (mat4Mul (csM1_get p1) (tRotDag_get p2)))

private def csCheck2Row (p1 : Fin 15) : Bool := allFinN 15 (csCheck2 p1)

private theorem cs2_row0 : csCheck2Row 0 = true := by decide
private theorem cs2_row1 : csCheck2Row 1 = true := by decide
private theorem cs2_row2 : csCheck2Row 2 = true := by decide
private theorem cs2_row3 : csCheck2Row 3 = true := by decide
private theorem cs2_row4 : csCheck2Row 4 = true := by decide
private theorem cs2_row5 : csCheck2Row 5 = true := by decide
private theorem cs2_row6 : csCheck2Row 6 = true := by decide
private theorem cs2_row7 : csCheck2Row 7 = true := by decide
private theorem cs2_row8 : csCheck2Row 8 = true := by decide
private theorem cs2_row9 : csCheck2Row 9 = true := by decide
private theorem cs2_row10 : csCheck2Row 10 = true := by decide
private theorem cs2_row11 : csCheck2Row 11 = true := by decide
private theorem cs2_row12 : csCheck2Row 12 = true := by decide
private theorem cs2_row13 : csCheck2Row 13 = true := by decide
private theorem cs2_row14 : csCheck2Row 14 = true := by decide

private theorem csCheck2Row_all (p1 : Fin 15) : allFinN 15 (csCheck2 p1) = true := by
  rcases p1 with ⟨v, hv⟩
  match v, hv with
  | 0, _ => exact cs2_row0 | 1, _ => exact cs2_row1 | 2, _ => exact cs2_row2
  | 3, _ => exact cs2_row3 | 4, _ => exact cs2_row4 | 5, _ => exact cs2_row5
  | 6, _ => exact cs2_row6 | 7, _ => exact cs2_row7 | 8, _ => exact cs2_row8
  | 9, _ => exact cs2_row9 | 10, _ => exact cs2_row10 | 11, _ => exact cs2_row11
  | 12, _ => exact cs2_row12 | 13, _ => exact cs2_row13 | 14, _ => exact cs2_row14

private theorem tRotDag_mat4_eq (p : Fin 15) (r c : Fin 4) :
    mat4Get (tRotDag_get p) r c = tRotDagScaled2 p r c := by
  have hq := beq_iff_eq.mp (allFinN_sound 15 _ tRotDag_get_eq_all p)
  rw [hq]
  rcases r with ⟨rv, hr⟩; rcases c with ⟨cv, hc⟩
  match rv, hr, cv, hc with
  | 0, _, 0, _ => rfl | 0, _, 1, _ => rfl | 0, _, 2, _ => rfl | 0, _, 3, _ => rfl
  | 1, _, 0, _ => rfl | 1, _, 1, _ => rfl | 1, _, 2, _ => rfl | 1, _, 3, _ => rfl
  | 2, _, 0, _ => rfl | 2, _, 1, _ => rfl | 2, _, 2, _ => rfl | 2, _, 3, _ => rfl
  | 3, _, 0, _ => rfl | 3, _, 1, _ => rfl | 3, _, 2, _ => rfl | 3, _, 3, _ => rfl

theorem cs_minT_gt_2_exhaustive :
    isClifford2qFast 1 gCS_mat4 = false ∧
    (∀ p1 : Fin 15, isClifford2qFast 2 (mat4Mul gCS_mat4 (tRotDag_get p1)) = false) ∧
    (∀ p1 p2 : Fin 15,
      isClifford2qFast 4 (mat4Mul (mat4Mul gCS_mat4 (tRotDag_get p1)) (tRotDag_get p2)) = false) := by
  refine ⟨by decide, fun p1 => ?_, fun p1 p2 => ?_⟩
  · have hm1 := beq_iff_eq.mp (allFinN_sound 15 _ csM1_get_eq_all p1)
    rw [hm1]
    have h := allFinN_sound 15 csCheck1 (by decide) p1
    simpa [csCheck1] using h
  · have hm1 := beq_iff_eq.mp (allFinN_sound 15 _ csM1_get_eq_all p1)
    rw [hm1]
    have h := allFinN_sound 15 (csCheck2 p1) (csCheck2Row_all p1) p2
    simpa [csCheck2] using h

-- Concrete unboxed 8x8 matrix type `Mat8` for linear-time exact Clifford+T composition (no thunk blowup!)
structure Vec8 (α : Type) where
  x0 : α
  x1 : α
  x2 : α
  x3 : α
  x4 : α
  x5 : α
  x6 : α
  x7 : α
  deriving DecidableEq, Repr

def vec8Get {α : Type} (v : Vec8 α) (i : Fin 8) : α :=
  match i.val with
  | 0 => v.x0 | 1 => v.x1 | 2 => v.x2 | 3 => v.x3
  | 4 => v.x4 | 5 => v.x5 | 6 => v.x6 | _ => v.x7

def vec8Of {α : Type} (f : Fin 8 → α) : Vec8 α :=
  ⟨f 0, f 1, f 2, f 3, f 4, f 5, f 6, f 7⟩

abbrev Mat8 := Vec8 (Vec8 Z8Int)

def mat8Get (M : Mat8) (r c : Fin 8) : Z8Int :=
  vec8Get (vec8Get M r) c

def mat8Of (f : Fin 8 → Fin 8 → Z8Int) : Mat8 :=
  vec8Of fun r => vec8Of fun c => f r c

def mat8Id : Mat8 :=
  mat8Of fun r c => if r == c then z8One else z8Zero

def dotVec8 (a b : Vec8 Z8Int) : Z8Int :=
  z8Add (z8Add (z8Add (z8Mul a.x0 b.x0) (z8Mul a.x1 b.x1))
               (z8Add (z8Mul a.x2 b.x2) (z8Mul a.x3 b.x3)))
        (z8Add (z8Add (z8Mul a.x4 b.x4) (z8Mul a.x5 b.x5))
               (z8Add (z8Mul a.x6 b.x6) (z8Mul a.x7 b.x7)))

def mat8Col (B : Mat8) (j : Fin 8) : Vec8 Z8Int :=
  ⟨vec8Get B.x0 j, vec8Get B.x1 j, vec8Get B.x2 j, vec8Get B.x3 j,
   vec8Get B.x4 j, vec8Get B.x5 j, vec8Get B.x6 j, vec8Get B.x7 j⟩

def mat8MulRow (r : Vec8 Z8Int) (c0 c1 c2 c3 c4 c5 c6 c7 : Vec8 Z8Int) : Vec8 Z8Int :=
  ⟨dotVec8 r c0, dotVec8 r c1, dotVec8 r c2, dotVec8 r c3,
   dotVec8 r c4, dotVec8 r c5, dotVec8 r c6, dotVec8 r c7⟩

def mat8Mul (A B : Mat8) : Mat8 :=
  let c0 := mat8Col B 0; let c1 := mat8Col B 1; let c2 := mat8Col B 2; let c3 := mat8Col B 3
  let c4 := mat8Col B 4; let c5 := mat8Col B 5; let c6 := mat8Col B 6; let c7 := mat8Col B 7
  ⟨mat8MulRow A.x0 c0 c1 c2 c3 c4 c5 c6 c7,
   mat8MulRow A.x1 c0 c1 c2 c3 c4 c5 c6 c7,
   mat8MulRow A.x2 c0 c1 c2 c3 c4 c5 c6 c7,
   mat8MulRow A.x3 c0 c1 c2 c3 c4 c5 c6 c7,
   mat8MulRow A.x4 c0 c1 c2 c3 c4 c5 c6 c7,
   mat8MulRow A.x5 c0 c1 c2 c3 c4 c5 c6 c7,
   mat8MulRow A.x6 c0 c1 c2 c3 c4 c5 c6 c7,
   mat8MulRow A.x7 c0 c1 c2 c3 c4 c5 c6 c7⟩

def composeFromMat8 (acc : Mat8) : List Mat8 → Mat8
  | [] => acc
  | g :: gs => composeFromMat8 (mat8Mul g acc) gs

def applySeqMat8 (gs : List Mat8) : Mat8 := composeFromMat8 mat8Id gs

private theorem composeFromMat8_append (acc : Mat8) (xs ys : List Mat8) :
    composeFromMat8 acc (xs ++ ys) = composeFromMat8 (composeFromMat8 acc xs) ys := by
  induction xs generalizing acc with
  | nil => rfl
  | cons x xs ih => simp [composeFromMat8, ih]

def gCNOT3 (ctrl tgt : Nat) : Mat8 :=
  mat8Of fun r c =>
    let bc := (c.val >>> (2 - ctrl)) % 2
    let y := if bc == 1 then c.val ^^^ (1 <<< (2 - tgt)) else c.val
    if r.val == y then z8One else z8Zero

def gT3 (q : Nat) (inv : Bool) : Mat8 :=
  let ph := if inv then z8Conj z8W else z8W
  mat8Of fun r c =>
    if r == c then (if ((c.val >>> (2 - q)) % 2) == 1 then ph else z8One) else z8Zero

def gH3Int (q : Nat) : Mat8 :=
  let mask := 1 <<< (2 - q)
  mat8Of fun r c =>
    if (r.val &&& (7 ^^^ mask)) == (c.val &&& (7 ^^^ mask)) then
      if ((r.val &&& mask) != 0) && ((c.val &&& mask) != 0) then z8MOne else z8One
    else z8Zero

def cczSeq : List Mat8 :=
  [gT3 0 false, gT3 1 false, gT3 2 false,
   gCNOT3 0 1, gT3 1 true,
   gCNOT3 0 2, gT3 2 true,
   gCNOT3 1 2, gT3 2 true,
   gCNOT3 0 2, gT3 2 false,
   gCNOT3 0 2, gCNOT3 1 2, gCNOT3 0 2, gCNOT3 0 1]

def gCCZ_target (r c : Fin 8) : Z8Int :=
  if r == c then (if r.val == 7 then z8MOne else z8One) else z8Zero

def gTOF_target (r c : Fin 8) : Z8Int :=
  let y := if (c.val &&& 6) == 6 then c.val ^^^ 1 else c.val
  if r.val == y then z8One else z8Zero

def gFRED_target (r c : Fin 8) : Z8Int :=
  let b0 := (c.val >>> 2) % 2
  let b1 := (c.val >>> 1) % 2
  let b2 := c.val % 2
  let y := if b0 == 1 && b1 != b2 then c.val ^^^ 3 else c.val
  if r.val == y then z8One else z8Zero

def gCCZ_7T : Mat8 := applySeqMat8 cczSeq
def gTOF_7T_scaled2 : Mat8 := applySeqMat8 ([gH3Int 2] ++ cczSeq ++ [gH3Int 2])
def gFRED_7T_scaled2 : Mat8 := applySeqMat8 ([gCNOT3 2 1, gH3Int 2] ++ cczSeq ++ [gH3Int 2, gCNOT3 2 1])
def gRTOF_4T_scaled2 : Mat8 :=
  applySeqMat8 [gH3Int 2, gT3 2 false, gCNOT3 1 2, gT3 2 true, gCNOT3 0 2, gT3 2 false, gCNOT3 1 2, gT3 2 true, gH3Int 2]

private def mat8Id_lit : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩⟩
private theorem mat8Id_eq_lit : mat8Id = mat8Id_lit := by decide
private def gT3_0_f_lit : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩⟩⟩
private theorem gT3_0_f_eq : (gT3 0 false) = gT3_0_f_lit := by decide
private def gT3_1_f_lit : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩⟩⟩
private theorem gT3_1_f_eq : (gT3 1 false) = gT3_1_f_lit := by decide
private def gT3_2_f_lit : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩⟩⟩
private theorem gT3_2_f_eq : (gT3 2 false) = gT3_2_f_lit := by decide
private def gT3_1_t_lit : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩⟩⟩
private theorem gT3_1_t_eq : (gT3 1 true) = gT3_1_t_lit := by decide
private def gT3_2_t_lit : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩⟩⟩
private theorem gT3_2_t_eq : (gT3 2 true) = gT3_2_t_lit := by decide
private def gCNOT3_0_1_lit : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem gCNOT3_0_1_eq : (gCNOT3 0 1) = gCNOT3_0_1_lit := by decide
private def gCNOT3_0_2_lit : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem gCNOT3_0_2_eq : (gCNOT3 0 2) = gCNOT3_0_2_lit := by decide
private def gCNOT3_1_2_lit : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem gCNOT3_1_2_eq : (gCNOT3 1 2) = gCNOT3_1_2_lit := by decide
private def gCNOT3_2_1_lit : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem gCNOT3_2_1_eq : (gCNOT3 2 1) = gCNOT3_2_1_lit := by decide
private def gH3_2_lit : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩⟩
private theorem gH3_2_eq : (gH3Int 2) = gH3_2_lit := by decide

private def ccz_mid0 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩⟩⟩
private theorem ccz_step0 : composeFromMat8 mat8Id_lit [gT3_0_f_lit] = ccz_mid0 := by decide
private def ccz_mid1 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩⟩⟩
private theorem ccz_step1 : composeFromMat8 ccz_mid0 [gT3_1_f_lit] = ccz_mid1 := by decide
private def ccz_mid2 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 1⟩⟩⟩
private theorem ccz_step2 : composeFromMat8 ccz_mid1 [gT3_2_f_lit] = ccz_mid2 := by decide
private def ccz_mid3 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem ccz_step3 : composeFromMat8 ccz_mid2 [gCNOT3_0_1_lit] = ccz_mid3 := by decide
private def ccz_mid4 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem ccz_step4 : composeFromMat8 ccz_mid3 [gT3_1_t_lit] = ccz_mid4 := by decide
private def ccz_mid5 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem ccz_step5 : composeFromMat8 ccz_mid4 [gCNOT3_0_2_lit] = ccz_mid5 := by decide
private def ccz_mid6 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem ccz_step6 : composeFromMat8 ccz_mid5 [gT3_2_t_lit] = ccz_mid6 := by decide
private def ccz_mid7 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem ccz_step7 : composeFromMat8 ccz_mid6 [gCNOT3_1_2_lit] = ccz_mid7 := by decide
private def ccz_mid8 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem ccz_step8 : composeFromMat8 ccz_mid7 [gT3_2_t_lit] = ccz_mid8 := by decide
private def ccz_mid9 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem ccz_step9 : composeFromMat8 ccz_mid8 [gCNOT3_0_2_lit] = ccz_mid9 := by decide
private def ccz_mid10 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem ccz_step10 : composeFromMat8 ccz_mid9 [gT3_2_f_lit] = ccz_mid10 := by decide
private def ccz_mid11 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem ccz_step11 : composeFromMat8 ccz_mid10 [gCNOT3_0_2_lit] = ccz_mid11 := by decide
private def ccz_mid12 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem ccz_step12 : composeFromMat8 ccz_mid11 [gCNOT3_1_2_lit] = ccz_mid12 := by decide
private def ccz_mid13 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem ccz_step13 : composeFromMat8 ccz_mid12 [gCNOT3_0_2_lit] = ccz_mid13 := by decide
def ccz_lit : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩⟩
private theorem ccz_step14 : composeFromMat8 ccz_mid13 [gCNOT3_0_1_lit] = ccz_lit := by decide
private theorem ccz_split : (cczSeq) = [gT3_0_f_lit] ++ ([gT3_1_f_lit] ++ ([gT3_2_f_lit] ++ ([gCNOT3_0_1_lit] ++ ([gT3_1_t_lit] ++ ([gCNOT3_0_2_lit] ++ ([gT3_2_t_lit] ++ ([gCNOT3_1_2_lit] ++ ([gT3_2_t_lit] ++ ([gCNOT3_0_2_lit] ++ ([gT3_2_f_lit] ++ ([gCNOT3_0_2_lit] ++ ([gCNOT3_1_2_lit] ++ ([gCNOT3_0_2_lit] ++ ([gCNOT3_0_1_lit])))))))))))))) := by unfold cczSeq; rw [gCNOT3_0_1_eq, gCNOT3_0_2_eq, gCNOT3_1_2_eq, gT3_0_f_eq, gT3_1_f_eq, gT3_1_t_eq, gT3_2_f_eq, gT3_2_t_eq]; rfl
private theorem gCCZ_7T_eq_lit : gCCZ_7T = ccz_lit := by
  dsimp [gCCZ_7T, applySeqMat8]
  rw [mat8Id_eq_lit, ccz_split, composeFromMat8_append, ccz_step0, composeFromMat8_append, ccz_step1, composeFromMat8_append, ccz_step2, composeFromMat8_append, ccz_step3, composeFromMat8_append, ccz_step4, composeFromMat8_append, ccz_step5, composeFromMat8_append, ccz_step6, composeFromMat8_append, ccz_step7, composeFromMat8_append, ccz_step8, composeFromMat8_append, ccz_step9, composeFromMat8_append, ccz_step10, composeFromMat8_append, ccz_step11, composeFromMat8_append, ccz_step12, composeFromMat8_append, ccz_step13, ccz_step14]

private def tof_mid0 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩⟩
private theorem tof_step0 : composeFromMat8 mat8Id_lit [gH3_2_lit] = tof_mid0 := by decide
private def tof_mid1 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 1, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩⟩⟩
private theorem tof_step1 : composeFromMat8 tof_mid0 [gT3_0_f_lit] = tof_mid1 := by decide
private def tof_mid2 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, -1, 0⟩⟩⟩
private theorem tof_step2 : composeFromMat8 tof_mid1 [gT3_1_f_lit] = tof_mid2 := by decide
private def tof_mid3 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, -1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, -1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, -1⟩⟩⟩
private theorem tof_step3 : composeFromMat8 tof_mid2 [gT3_2_f_lit] = tof_mid3 := by decide
private def tof_mid4 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, -1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, -1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, -1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem tof_step4 : composeFromMat8 tof_mid3 [gCNOT3_0_1_lit] = tof_mid4 := by decide
private def tof_mid5 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, -1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem tof_step5 : composeFromMat8 tof_mid4 [gT3_1_t_lit] = tof_mid5 := by decide
private def tof_mid6 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, -1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem tof_step6 : composeFromMat8 tof_mid5 [gCNOT3_0_2_lit] = tof_mid6 := by decide
private def tof_mid7 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, -1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 1, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem tof_step7 : composeFromMat8 tof_mid6 [gT3_2_t_lit] = tof_mid7 := by decide
private def tof_mid8 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, -1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 1, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem tof_step8 : composeFromMat8 tof_mid7 [gCNOT3_1_2_lit] = tof_mid8 := by decide
private def tof_mid9 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, -1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem tof_step9 : composeFromMat8 tof_mid8 [gT3_2_t_lit] = tof_mid9 := by decide
private def tof_mid10 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, -1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem tof_step10 : composeFromMat8 tof_mid9 [gCNOT3_0_2_lit] = tof_mid10 := by decide
private def tof_mid11 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem tof_step11 : composeFromMat8 tof_mid10 [gT3_2_f_lit] = tof_mid11 := by decide
private def tof_mid12 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem tof_step12 : composeFromMat8 tof_mid11 [gCNOT3_0_2_lit] = tof_mid12 := by decide
private def tof_mid13 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem tof_step13 : composeFromMat8 tof_mid12 [gCNOT3_1_2_lit] = tof_mid13 := by decide
private def tof_mid14 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem tof_step14 : composeFromMat8 tof_mid13 [gCNOT3_0_2_lit] = tof_mid14 := by decide
private def tof_mid15 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩⟩
private theorem tof_step15 : composeFromMat8 tof_mid14 [gCNOT3_0_1_lit] = tof_mid15 := by decide
def tof_lit : Mat8 :=
  ⟨⟨⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem tof_step16 : composeFromMat8 tof_mid15 [gH3_2_lit] = tof_lit := by decide
private theorem tof_split : (gH3Int 2 :: (cczSeq ++ [gH3Int 2])) = [gH3_2_lit] ++ ([gT3_0_f_lit] ++ ([gT3_1_f_lit] ++ ([gT3_2_f_lit] ++ ([gCNOT3_0_1_lit] ++ ([gT3_1_t_lit] ++ ([gCNOT3_0_2_lit] ++ ([gT3_2_t_lit] ++ ([gCNOT3_1_2_lit] ++ ([gT3_2_t_lit] ++ ([gCNOT3_0_2_lit] ++ ([gT3_2_f_lit] ++ ([gCNOT3_0_2_lit] ++ ([gCNOT3_1_2_lit] ++ ([gCNOT3_0_2_lit] ++ ([gCNOT3_0_1_lit] ++ ([gH3_2_lit])))))))))))))))) := by unfold cczSeq; rw [gCNOT3_0_1_eq, gCNOT3_0_2_eq, gCNOT3_1_2_eq, gH3_2_eq, gT3_0_f_eq, gT3_1_f_eq, gT3_1_t_eq, gT3_2_f_eq, gT3_2_t_eq]; rfl
private theorem gTOF_7T_scaled2_eq_lit : gTOF_7T_scaled2 = tof_lit := by
  dsimp [gTOF_7T_scaled2, applySeqMat8]
  rw [mat8Id_eq_lit, tof_split, composeFromMat8_append, tof_step0, composeFromMat8_append, tof_step1, composeFromMat8_append, tof_step2, composeFromMat8_append, tof_step3, composeFromMat8_append, tof_step4, composeFromMat8_append, tof_step5, composeFromMat8_append, tof_step6, composeFromMat8_append, tof_step7, composeFromMat8_append, tof_step8, composeFromMat8_append, tof_step9, composeFromMat8_append, tof_step10, composeFromMat8_append, tof_step11, composeFromMat8_append, tof_step12, composeFromMat8_append, tof_step13, composeFromMat8_append, tof_step14, composeFromMat8_append, tof_step15, tof_step16]

private def fred_mid0 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem fred_step0 : composeFromMat8 mat8Id_lit [gCNOT3_2_1_lit] = fred_mid0 := by decide
private def fred_mid1 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem fred_step1 : composeFromMat8 fred_mid0 [gH3_2_lit] = fred_mid1 := by decide
private def fred_mid2 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, -1, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem fred_step2 : composeFromMat8 fred_mid1 [gT3_0_f_lit] = fred_mid2 := by decide
private def fred_mid3 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, -1, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, -1, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem fred_step3 : composeFromMat8 fred_mid2 [gT3_1_f_lit] = fred_mid3 := by decide
private def fred_mid4 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, -1, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, -1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem fred_step4 : composeFromMat8 fred_mid3 [gT3_2_f_lit] = fred_mid4 := by decide
private def fred_mid5 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, -1, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, -1, 0⟩⟩⟩
private theorem fred_step5 : composeFromMat8 fred_mid4 [gCNOT3_0_1_lit] = fred_mid5 := by decide
private def fred_mid6 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, -1, 0, 0⟩⟩⟩
private theorem fred_step6 : composeFromMat8 fred_mid5 [gT3_1_t_lit] = fred_mid6 := by decide
private def fred_mid7 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, -1, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩⟩
private theorem fred_step7 : composeFromMat8 fred_mid6 [gCNOT3_0_2_lit] = fred_mid7 := by decide
private def fred_mid8 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, -1, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩⟩⟩
private theorem fred_step8 : composeFromMat8 fred_mid7 [gT3_2_t_lit] = fred_mid8 := by decide
private def fred_mid9 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, -1, 0, 0⟩⟩⟩
private theorem fred_step9 : composeFromMat8 fred_mid8 [gCNOT3_1_2_lit] = fred_mid9 := by decide
private def fred_mid10 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩⟩
private theorem fred_step10 : composeFromMat8 fred_mid9 [gT3_2_t_lit] = fred_mid10 := by decide
private def fred_mid11 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 1⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩⟩⟩
private theorem fred_step11 : composeFromMat8 fred_mid10 [gCNOT3_0_2_lit] = fred_mid11 := by decide
private def fred_mid12 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩⟩
private theorem fred_step12 : composeFromMat8 fred_mid11 [gT3_2_f_lit] = fred_mid12 := by decide
private def fred_mid13 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩⟩
private theorem fred_step13 : composeFromMat8 fred_mid12 [gCNOT3_0_2_lit] = fred_mid13 := by decide
private def fred_mid14 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩⟩
private theorem fred_step14 : composeFromMat8 fred_mid13 [gCNOT3_1_2_lit] = fred_mid14 := by decide
private def fred_mid15 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩⟩
private theorem fred_step15 : composeFromMat8 fred_mid14 [gCNOT3_0_2_lit] = fred_mid15 := by decide
private def fred_mid16 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem fred_step16 : composeFromMat8 fred_mid15 [gCNOT3_0_1_lit] = fred_mid16 := by decide
private def fred_mid17 : Mat8 :=
  ⟨⟨⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem fred_step17 : composeFromMat8 fred_mid16 [gH3_2_lit] = fred_mid17 := by decide
def fred_lit : Mat8 :=
  ⟨⟨⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩⟩⟩
private theorem fred_step18 : composeFromMat8 fred_mid17 [gCNOT3_2_1_lit] = fred_lit := by decide
private theorem fred_split : (gCNOT3 2 1 :: gH3Int 2 :: (cczSeq ++ [gH3Int 2, gCNOT3 2 1])) = [gCNOT3_2_1_lit] ++ ([gH3_2_lit] ++ ([gT3_0_f_lit] ++ ([gT3_1_f_lit] ++ ([gT3_2_f_lit] ++ ([gCNOT3_0_1_lit] ++ ([gT3_1_t_lit] ++ ([gCNOT3_0_2_lit] ++ ([gT3_2_t_lit] ++ ([gCNOT3_1_2_lit] ++ ([gT3_2_t_lit] ++ ([gCNOT3_0_2_lit] ++ ([gT3_2_f_lit] ++ ([gCNOT3_0_2_lit] ++ ([gCNOT3_1_2_lit] ++ ([gCNOT3_0_2_lit] ++ ([gCNOT3_0_1_lit] ++ ([gH3_2_lit] ++ ([gCNOT3_2_1_lit])))))))))))))))))) := by unfold cczSeq; rw [gCNOT3_0_1_eq, gCNOT3_0_2_eq, gCNOT3_1_2_eq, gCNOT3_2_1_eq, gH3_2_eq, gT3_0_f_eq, gT3_1_f_eq, gT3_1_t_eq, gT3_2_f_eq, gT3_2_t_eq]; rfl
private theorem gFRED_7T_scaled2_eq_lit : gFRED_7T_scaled2 = fred_lit := by
  dsimp [gFRED_7T_scaled2, applySeqMat8]
  rw [mat8Id_eq_lit, fred_split, composeFromMat8_append, fred_step0, composeFromMat8_append, fred_step1, composeFromMat8_append, fred_step2, composeFromMat8_append, fred_step3, composeFromMat8_append, fred_step4, composeFromMat8_append, fred_step5, composeFromMat8_append, fred_step6, composeFromMat8_append, fred_step7, composeFromMat8_append, fred_step8, composeFromMat8_append, fred_step9, composeFromMat8_append, fred_step10, composeFromMat8_append, fred_step11, composeFromMat8_append, fred_step12, composeFromMat8_append, fred_step13, composeFromMat8_append, fred_step14, composeFromMat8_append, fred_step15, composeFromMat8_append, fred_step16, composeFromMat8_append, fred_step17, fred_step18]

private def rtof_mid0 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩⟩⟩
private theorem rtof_step0 : composeFromMat8 mat8Id_lit [gH3_2_lit] = rtof_mid0 := by decide
private def rtof_mid1 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩⟩⟩
private theorem rtof_step1 : composeFromMat8 rtof_mid0 [gT3_2_f_lit] = rtof_mid1 := by decide
private def rtof_mid2 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩⟩⟩
private theorem rtof_step2 : composeFromMat8 rtof_mid1 [gCNOT3_1_2_lit] = rtof_mid2 := by decide
private def rtof_mid3 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, -1⟩⟩⟩
private theorem rtof_step3 : composeFromMat8 rtof_mid2 [gT3_2_t_lit] = rtof_mid3 := by decide
private def rtof_mid4 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, -1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩⟩⟩
private theorem rtof_step4 : composeFromMat8 rtof_mid3 [gCNOT3_0_2_lit] = rtof_mid4 := by decide
private def rtof_mid5 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, -1⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, -1, 0⟩⟩⟩
private theorem rtof_step5 : composeFromMat8 rtof_mid4 [gT3_2_f_lit] = rtof_mid5 := by decide
private def rtof_mid6 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, -1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 1, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, -1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, -1⟩, ⟨0, 0, 0, -1⟩⟩⟩
private theorem rtof_step6 : composeFromMat8 rtof_mid5 [gCNOT3_1_2_lit] = rtof_mid6 := by decide
private def rtof_mid7 : Mat8 :=
  ⟨⟨⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨-1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨1, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 1, 0⟩, ⟨0, 0, -1, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, -1, 0⟩, ⟨0, 0, -1, 0⟩⟩⟩
private theorem rtof_step7 : composeFromMat8 rtof_mid6 [gT3_2_t_lit] = rtof_mid7 := by decide
def rtof_lit : Mat8 :=
  ⟨⟨⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨-2, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, -2, 0⟩⟩,
   ⟨⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 0, 0⟩, ⟨0, 0, 2, 0⟩, ⟨0, 0, 0, 0⟩⟩⟩
private theorem rtof_step8 : composeFromMat8 rtof_mid7 [gH3_2_lit] = rtof_lit := by decide
private theorem rtof_split : ([gH3Int 2, gT3 2 false, gCNOT3 1 2, gT3 2 true, gCNOT3 0 2, gT3 2 false, gCNOT3 1 2, gT3 2 true, gH3Int 2]) = [gH3_2_lit] ++ ([gT3_2_f_lit] ++ ([gCNOT3_1_2_lit] ++ ([gT3_2_t_lit] ++ ([gCNOT3_0_2_lit] ++ ([gT3_2_f_lit] ++ ([gCNOT3_1_2_lit] ++ ([gT3_2_t_lit] ++ ([gH3_2_lit])))))))) := by rw [gCNOT3_0_2_eq, gCNOT3_1_2_eq, gH3_2_eq, gT3_2_f_eq, gT3_2_t_eq]; rfl
private theorem gRTOF_4T_scaled2_eq_lit : gRTOF_4T_scaled2 = rtof_lit := by
  dsimp [gRTOF_4T_scaled2, applySeqMat8]
  rw [mat8Id_eq_lit, rtof_split, composeFromMat8_append, rtof_step0, composeFromMat8_append, rtof_step1, composeFromMat8_append, rtof_step2, composeFromMat8_append, rtof_step3, composeFromMat8_append, rtof_step4, composeFromMat8_append, rtof_step5, composeFromMat8_append, rtof_step6, composeFromMat8_append, rtof_step7, rtof_step8]


theorem ccz_7T_exact_identity_and_unitary :
    (∀ i j : Fin 8, mat8Get gCCZ_7T i j = gCCZ_target i j) ∧
    z8IsUnitaryScaled 1 (mat8Get gCCZ_7T) = true := by
  rw [gCCZ_7T_eq_lit]
  exact ⟨z8MatEq_sound (by decide), by decide⟩

theorem toffoli_7T_exact_identity_and_unitary :
    (∀ i j : Fin 8, mat8Get gTOF_7T_scaled2 i j = z8ScaleInt 2 (gTOF_target i j)) ∧
    z8IsUnitaryScaled 4 (mat8Get gTOF_7T_scaled2) = true := by
  rw [gTOF_7T_scaled2_eq_lit]
  exact ⟨z8MatEq_sound (by decide), by decide⟩

theorem fredkin_7T_exact_identity_and_unitary :
    (∀ i j : Fin 8, mat8Get gFRED_7T_scaled2 i j = z8ScaleInt 2 (gFRED_target i j)) ∧
    z8IsUnitaryScaled 4 (mat8Get gFRED_7T_scaled2) = true := by
  rw [gFRED_7T_scaled2_eq_lit]
  exact ⟨z8MatEq_sound (by decide), by decide⟩

private def rtofOffDiagLitOk (r c : Fin 8) : Bool :=
  (r == c) || (z8MatMul (mat8Get rtof_lit) (fun u v => z8Conj (gTOF_target v u)) r c == z8Zero)

theorem rtof_4T_unitary_and_relative_phase :
    z8IsUnitaryScaled 4 (mat8Get gRTOF_4T_scaled2) = true ∧
    (∀ i j : Fin 8, i ≠ j → z8MatMul (mat8Get gRTOF_4T_scaled2) (fun r c => z8Conj (gTOF_target c r)) i j = z8Zero) := by
  rw [gRTOF_4T_scaled2_eq_lit]
  refine ⟨by decide, fun i j hij => ?_⟩
  have h := allFinN_sound 8 (rtofOffDiagLitOk i) (allFinN_sound 8 (fun r => allFinN 8 (rtofOffDiagLitOk r)) (by decide) i) j
  unfold rtofOffDiagLitOk at h
  rw [beq_false_of_ne hij] at h
  exact beq_iff_eq.mp h

-- Universal algebraic lower bound: ANY 3-qubit CNOT+T phase polynomial realizing CCZ requires 7 odd terms (min T = 7)
def bitDot3 (a x : Nat) : Nat :=
  ((a % 2) * (x % 2) + ((a / 2) % 2) * ((x / 2) % 2) + ((a / 4) % 2) * ((x / 4) % 2)) % 2

def phasePoly3 (c : Fin 7 → Fin 8) (x : Fin 7) : Nat :=
  let xv := x.val + 1
  ((c 0).val * bitDot3 1 xv + (c 1).val * bitDot3 2 xv + (c 2).val * bitDot3 3 xv +
   (c 3).val * bitDot3 4 xv + (c 4).val * bitDot3 5 xv + (c 5).val * bitDot3 6 xv +
   (c 6).val * bitDot3 7 xv) % 8

private theorem fin7_cases (b : Fin 7) :
    b = 0 ∨ b = 1 ∨ b = 2 ∨ b = 3 ∨ b = 4 ∨ b = 5 ∨ b = 6 := by
  rcases b with ⟨v, hv⟩
  match v, hv with
  | 0, _ => exact Or.inl rfl | 1, _ => exact Or.inr (Or.inl rfl)
  | 2, _ => exact Or.inr (Or.inr (Or.inl rfl)) | 3, _ => exact Or.inr (Or.inr (Or.inr (Or.inl rfl)))
  | 4, _ => exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inl rfl))))
  | 5, _ => exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl rfl)))))
  | 6, _ => exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr rfl)))))

def phasePoly3Eval (c0 c1 c2 c3 c4 c5 c6 : Nat) (x : Fin 7) : Nat :=
  match x.val with
  | 0 => (c0 + c2 + c4 + c6) % 8
  | 1 => (c1 + c2 + c5 + c6) % 8
  | 2 => (c0 + c1 + c4 + c5) % 8
  | 3 => (c3 + c4 + c5 + c6) % 8
  | 4 => (c0 + c2 + c3 + c5) % 8
  | 5 => (c1 + c2 + c3 + c4) % 8
  | _ => (c0 + c1 + c3 + c6) % 8

private theorem phasePoly3_eq_eval (c : Fin 7 → Fin 8) (x : Fin 7) :
    phasePoly3 c x = phasePoly3Eval (c 0).val (c 1).val (c 2).val (c 3).val (c 4).val (c 5).val (c 6).val x := by
  rcases fin7_cases x with rfl | rfl | rfl | rfl | rfl | rfl | rfl
  · change ((c 0).val * 1 + (c 1).val * 0 + (c 2).val * 1 + (c 3).val * 0 + (c 4).val * 1 + (c 5).val * 0 + (c 6).val * 1) % 8 = _; simp [phasePoly3Eval]; try rfl
  · change ((c 0).val * 0 + (c 1).val * 1 + (c 2).val * 1 + (c 3).val * 0 + (c 4).val * 0 + (c 5).val * 1 + (c 6).val * 1) % 8 = _; simp [phasePoly3Eval]; try rfl
  · change ((c 0).val * 1 + (c 1).val * 1 + (c 2).val * 0 + (c 3).val * 0 + (c 4).val * 1 + (c 5).val * 1 + (c 6).val * 0) % 8 = _; simp [phasePoly3Eval]; try rfl
  · change ((c 0).val * 0 + (c 1).val * 0 + (c 2).val * 0 + (c 3).val * 1 + (c 4).val * 1 + (c 5).val * 1 + (c 6).val * 1) % 8 = _; simp [phasePoly3Eval]; try rfl
  · change ((c 0).val * 1 + (c 1).val * 0 + (c 2).val * 1 + (c 3).val * 1 + (c 4).val * 0 + (c 5).val * 1 + (c 6).val * 0) % 8 = _; simp [phasePoly3Eval]; try rfl
  · change ((c 0).val * 0 + (c 1).val * 1 + (c 2).val * 1 + (c 3).val * 1 + (c 4).val * 1 + (c 5).val * 0 + (c 6).val * 0) % 8 = _; simp [phasePoly3Eval]; try rfl
  · change ((c 0).val * 1 + (c 1).val * 1 + (c 2).val * 0 + (c 3).val * 1 + (c 4).val * 0 + (c 5).val * 0 + (c 6).val * 1) % 8 = _; simp [phasePoly3Eval]; try rfl

private theorem sum4_mod8_0004 (a b c d : Nat)
    (ha : a % 8 = 0) (hb : b % 8 = 0) (hc : c % 8 = 0) (hd : d % 8 = 4) :
    (a + b + c + d) % 8 = 4 := by simp [Nat.add_mod, ha, hb, hc, hd]

private theorem sum4_mod8_0000 (a b c d : Nat)
    (ha : a % 8 = 0) (hb : b % 8 = 0) (hc : c % 8 = 0) (hd : d % 8 = 0) :
    (a + b + c + d) % 8 = 0 := by simp [Nat.add_mod, ha, hb, hc, hd]

private theorem sum3_mod8_000 (a b c : Nat)
    (ha : a % 8 = 0) (hb : b % 8 = 0) (hc : c % 8 = 0) :
    (a + b + c) % 8 = 0 := by simp [Nat.add_mod, ha, hb, hc]

private theorem sum3_mod8_004 (a b c : Nat)
    (ha : a % 8 = 0) (hb : b % 8 = 0) (hc : c % 8 = 4) :
    (a + b + c) % 8 = 4 := by simp [Nat.add_mod, ha, hb, hc]

private theorem odd_of_sum_diff4_a (s4 s3 c : Nat) (heq : s4 = s3 + 4 * c)
    (h4 : s4 % 8 = 4) (h3 : s3 % 8 = 0) : c % 2 = 1 := by omega

private theorem odd_of_sum_diff4_b (s4 s3 c : Nat) (heq : s4 = s3 + 4 * c)
    (h4 : s4 % 8 = 0) (h3 : s3 % 8 = 4) : c % 2 = 1 := by omega

private theorem walsh_diff (c0 c1 c2 c3 c4 c5 c6 : Nat) :
    let s0 := c0 + c2 + c4 + c6
    let s1 := c1 + c2 + c5 + c6
    let s2 := c0 + c1 + c4 + c5
    let s3 := c3 + c4 + c5 + c6
    let s4 := c0 + c2 + c3 + c5
    let s5 := c1 + c2 + c3 + c4
    let s6 := c0 + c1 + c3 + c6
    (s0 + s2 + s4 + s6 = s1 + s3 + s5 + 4 * c0) ∧
    (s1 + s2 + s5 + s6 = s0 + s3 + s4 + 4 * c1) ∧
    (s0 + s1 + s4 + s5 = s2 + s3 + s6 + 4 * c2) ∧
    (s3 + s4 + s5 + s6 = s0 + s1 + s2 + 4 * c3) ∧
    (s0 + s2 + s3 + s5 = s1 + s4 + s6 + 4 * c4) ∧
    (s1 + s2 + s3 + s4 = s0 + s5 + s6 + 4 * c5) ∧
    (s0 + s1 + s3 + s6 = s2 + s4 + s5 + 4 * c6) := by
  refine ⟨by omega, by omega, by omega, by omega, by omega, by omega, by omega⟩

theorem ccz_phase_poly_min7 (c : Fin 7 → Fin 8)
    (hspec : ∀ x : Fin 7, phasePoly3 c x = if x.val == 6 then 4 else 0) :
    ∀ b : Fin 7, (c b).val % 2 = 1 := by
  have hev : ∀ x : Fin 7, phasePoly3Eval (c 0).val (c 1).val (c 2).val (c 3).val (c 4).val (c 5).val (c 6).val x =
      if x.val == 6 then 4 else 0 := fun x => by rw [← phasePoly3_eq_eval]; exact hspec x
  let s0 := (c 0).val + (c 2).val + (c 4).val + (c 6).val
  let s1 := (c 1).val + (c 2).val + (c 5).val + (c 6).val
  let s2 := (c 0).val + (c 1).val + (c 4).val + (c 5).val
  let s3 := (c 3).val + (c 4).val + (c 5).val + (c 6).val
  let s4 := (c 0).val + (c 2).val + (c 3).val + (c 5).val
  let s5 := (c 1).val + (c 2).val + (c 3).val + (c 4).val
  let s6 := (c 0).val + (c 1).val + (c 3).val + (c 6).val
  have h0 : s0 % 8 = 0 := hev 0
  have h1 : s1 % 8 = 0 := hev 1
  have h2 : s2 % 8 = 0 := hev 2
  have h3 : s3 % 8 = 0 := hev 3
  have h4 : s4 % 8 = 0 := hev 4
  have h5 : s5 % 8 = 0 := hev 5
  have h6 : s6 % 8 = 4 := hev 6
  rcases walsh_diff (c 0).val (c 1).val (c 2).val (c 3).val (c 4).val (c 5).val (c 6).val with
    ⟨e0, e1, e2, e3, e4, e5, e6⟩
  intro b
  rcases fin7_cases b with rfl | rfl | rfl | rfl | rfl | rfl | rfl
  · exact odd_of_sum_diff4_a (s0 + s2 + s4 + s6) (s1 + s3 + s5) (c 0).val e0 (sum4_mod8_0004 s0 s2 s4 s6 h0 h2 h4 h6) (sum3_mod8_000 s1 s3 s5 h1 h3 h5)
  · exact odd_of_sum_diff4_a (s1 + s2 + s5 + s6) (s0 + s3 + s4) (c 1).val e1 (sum4_mod8_0004 s1 s2 s5 s6 h1 h2 h5 h6) (sum3_mod8_000 s0 s3 s4 h0 h3 h4)
  · exact odd_of_sum_diff4_b (s0 + s1 + s4 + s5) (s2 + s3 + s6) (c 2).val e2 (sum4_mod8_0000 s0 s1 s4 s5 h0 h1 h4 h5) (sum3_mod8_004 s2 s3 s6 h2 h3 h6)
  · exact odd_of_sum_diff4_a (s3 + s4 + s5 + s6) (s0 + s1 + s2) (c 3).val e3 (sum4_mod8_0004 s3 s4 s5 s6 h3 h4 h5 h6) (sum3_mod8_000 s0 s1 s2 h0 h1 h2)
  · exact odd_of_sum_diff4_b (s0 + s2 + s3 + s5) (s1 + s4 + s6) (c 4).val e4 (sum4_mod8_0000 s0 s2 s3 s5 h0 h2 h3 h5) (sum3_mod8_004 s1 s4 s6 h1 h4 h6)
  · exact odd_of_sum_diff4_b (s1 + s2 + s3 + s4) (s0 + s5 + s6) (c 5).val e5 (sum4_mod8_0000 s1 s2 s3 s4 h1 h2 h3 h4) (sum3_mod8_004 s0 s5 s6 h0 h5 h6)
  · exact odd_of_sum_diff4_a (s0 + s1 + s3 + s6) (s2 + s4 + s5) (c 6).val e6 (sum4_mod8_0004 s0 s1 s3 s6 h0 h1 h3 h6) (sum3_mod8_000 s2 s4 s5 h2 h4 h5)

end Qua.QeccCliffordT
