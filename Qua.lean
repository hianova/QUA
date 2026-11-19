/-
  QUA — Pure Quantum Circuit Parameter Config (`import Flow`, 0 Custom Theorems)
  ==============================================================================
  SSOT compliance: all QTT linearity proofs and `.fir` Morphism IR encoders live
  in domain-free `../FLOW/Flow/*.lean`.
  This file configures QUA circuits via `Flow.Morphism.Builder`, emits two
  `.linear` quantum targets + one `.cartesian` guard witness, and relies on
  `Flow.Morphism.emitFirChecked` (`wf = true := by decide`) with 0 theorems.
-/
import Flow

set_option maxRecDepth 8192

namespace Qua

def gf256SboxMacGraph : Flow.Morphism.Graph :=
  let b := Flow.Morphism.Builder.init "gf256_sbox_mac" .linear true
  let (b, a)     := b.input "a" 8
  let (b, bv)    := b.input "b" 8
  let (b, c)     := b.input "c" 8
  let (b, sel)   := b.input "sel" 1
  let (b, cst63) := b.constU64 "cst63" 8 0x63
  let (b, invA) := b.gf2mInv "inv_a" 8 a 0x11B; let (b, prod) := b.gf2mMul "prod" 8 invA bv 0x11B
  let (b, aff)  := b.xorVec "aff" 8 prod cst63; let (b, sum) := b.add "sum" 8 c aff
  let (b, y)    := b.ctrlSel "y" 8 sel prod sum
  b.build [y]

def revAluFxpGraph : Flow.Morphism.Graph :=
  let b := Flow.Morphism.Builder.init "rev_alu_fxp" .linear true
  let (b, x) := b.input "x" 8; let (b, y) := b.input "y" 8
  let (b, sq) := b.sqrFxp "sq" 8 x 4 4 4; let (b, sum) := b.add "sum" 8 sq y
  let (b, ez) := b.eqZero "ez" x; let (b, res) := b.ctrlSel "res" 8 ez sum y
  b.build [res]

def cartesianRejectGraph : Flow.Morphism.Graph :=
  let (b, s0) := (Flow.Morphism.Builder.init "cartesian_guard" .cartesian false).input "state_in" 512
  let (b, p0) := b.stagePipe "stage0" 64 [s0] 0 17 64 0
  b.build [p0]

def adderCsaPhaseGraph : Flow.Morphism.Graph :=
  let b := Flow.Morphism.Builder.init "adder_csa_phase" .linear true
  let (b, a) := b.input "a" 8; let (b, bv) := b.input "b" 8; let (b, c) := b.input "c" 8; let (b, d) := b.input "d" 8
  let (b, s1) := b.add "s1" 8 a bv; let (b, s2) := b.add "s2" 8 s1 c; let (b, sum) := b.add "sum" 8 s2 d
  b.build [sum]

def gf12ArithGraph : Flow.Morphism.Graph :=
  let b := Flow.Morphism.Builder.init "gf12_arith" .linear true
  let (b, a) := b.input "a" 12; let (b, bv) := b.input "b" 12
  let (b, ainv) := b.gf2mInv "ainv" 12 a 0x1053; let (b, prod) := b.gf2mMul "prod" 12 ainv bv 0x1053
  let (b, res) := b.xorVec "res" 12 prod a
  b.build [res]

def wideLutCapGraph : Flow.Morphism.Graph :=
  let b := Flow.Morphism.Builder.init "wide_lut_cap" .linear true
  let (b, idx) := b.input "idx" 4
  let tblWords : List UInt64 := [
    0x0000, 0x0101, 0x0202, 0x0303, 0x0404, 0x0505, 0x0606, 0x0707,
    0x0808, 0x0909, 0x0a0a, 0x0b0b, 0x0c0c, 0x0d0d, 0x0e0e, 0x0f0f ]
  let (b, tbl) := b.op .lut_rom "tbl" 16 [idx] (lutWords := tblWords)
  b.build [tbl]

def gf256TowerSpecGraph : Flow.Morphism.Graph :=
  let (b, a)   := (Flow.Morphism.Builder.init "gf256_tower_spec" .linear true).input "a" 8
  let (b, inv) := b.gf2mQuadInv "q" 4 a 0x13 0x8
  let (b, y)   := b.slice "y" 8 inv 0
  b.build [y]

def lut5BijectionSpecGraph : Flow.Morphism.Graph :=
  let (b, x) := (Flow.Morphism.Builder.init "lut5_bijection_spec" .linear true).input "x" 5
  let tbl : List UInt64 := [
    0x04,0x0b,0x1f,0x14,0x1a,0x15,0x09,0x02,0x1b,0x05,0x08,0x12,0x1d,0x03,0x06,0x1c,
    0x1e,0x13,0x07,0x0e,0x00,0x0d,0x11,0x18,0x10,0x0c,0x01,0x19,0x16,0x0a,0x0f,0x17 ]
  let (b, y) := b.op .lut_rom "y" 5 [x] (lutWords := tbl)
  b.build [y]

def main (args : List String) : IO UInt32 := do
  let outDir : System.FilePath := match args with
    | d :: _ => System.FilePath.mk d
    | []     => System.FilePath.mk "build"
  IO.FS.createDirAll outDir
  if !(← Flow.Morphism.emitFirChecked (outDir / "gf256_sbox_mac.fir") gf256SboxMacGraph) then return 1
  if !(← Flow.Morphism.emitFirChecked (outDir / "rev_alu_fxp.fir") revAluFxpGraph) then return 1
  if !(← Flow.Morphism.emitFirChecked (outDir / "cartesian_reject.fir") cartesianRejectGraph) then return 1
  if !(← Flow.Morphism.emitFirChecked (outDir / "adder_csa_phase.fir") adderCsaPhaseGraph) then return 1
  if !(← Flow.Morphism.emitFirChecked (outDir / "gf12_arith.fir") gf12ArithGraph) then return 1
  if !(← Flow.Morphism.emitFirChecked (outDir / "wide_lut_cap.fir") wideLutCapGraph) then return 1
  if !(← Flow.Morphism.emitFirChecked (outDir / "gf256_tower_spec.fir") gf256TowerSpecGraph) then return 1
  if !(← Flow.Morphism.emitFirChecked (outDir / "lut5_bijection_spec.fir") lut5BijectionSpecGraph) then return 1
  IO.println s!"  [QUA-LEAN] Pure config emitted 7x linear .fir + 1x cartesian QTT guard into {outDir} (0 local theorems)"
  return 0

end Qua

def main (args : List String) : IO UInt32 :=
  Qua.main args
