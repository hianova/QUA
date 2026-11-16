#!/usr/bin/env python3
"""
verify_harness.py — Local Kickmix Fuzz-Verifier, Official Score Runner & Submission Packager
Google Quantum AI — Verifiable Quantum Circuit Compilation (./kaggle/qua)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import pathlib
import struct
import subprocess
import sys
import zipfile

import kickmix as km
import numpy as np

# Exact verified official sample_submission reference costs
REF_COST_INVSQRT = 2076774.233715282
REF_COST_DFTHCSOS = 3355599588.3155375


def exact_dqi_ref_cost(num_shots: int = 102400, log2_delta: float = -128.0) -> float:
    eps_plus = -math.expm1(log2_delta * math.log(2.0) / num_shots)
    rep = math.log(2.0) / -math.log(0.5 + 2.0 * math.sqrt(eps_plus))
    # Static 5,725,461 minus 351 Chien-inversion Toffolis skipped when classical E_j == 0 (p = 1/4096 across 4095 positions)
    mean_toff = 5725461.0 - 351.0 * (4095.0 / 4096.0)
    return rep * (1893.0 ** 2) * mean_toff


REF_COST_DQI = exact_dqi_ref_cost(102400, -128.0)


def unpack_shots(vals: np.ndarray, width: int) -> np.ndarray:
    shifts = np.arange(width, dtype=np.uint64)[:, None]
    return ((vals[None, :] >> shifts) & 1).astype(np.bool_)


def pack_shots(bits_2d: np.ndarray) -> np.ndarray:
    b = np.asarray(bits_2d, dtype=np.uint64)
    shifts = np.arange(b.shape[0], dtype=np.uint64)[:, None]
    return np.sum(b << shifts, axis=0, dtype=np.uint64)


def fuzz_verify_adder16(kmb_path: pathlib.Path, shots: int = 4096, seed: int = 20261115) -> None:
    circ = km.Circuit.from_file(str(kmb_path))
    stats = circ.operation_statistics()
    assert circ.num_qubits == 47, f"Expected 47 qubits, got {circ.num_qubits}"
    assert stats.num_toffoli_like == 15, f"Expected 15 Toffolis, got {stats.num_toffoli_like}"

    rng = np.random.default_rng(seed)
    sim = km.Simulator(batch_size=shots)
    sim.use_same_registers_as(circ)

    lhs = rng.integers(0, 1 << 16, size=shots, dtype=np.uint64)
    tgt = rng.integers(0, 1 << 16, size=shots, dtype=np.uint64)

    sim.write_across_shots("lhs", unpack_shots(lhs, 16))
    sim.write_across_shots("target", unpack_shots(tgt, 16))
    sim.write_across_shots("@clean", np.zeros((15, shots), dtype=np.bool_))
    sim.do(circ)

    got_lhs = pack_shots(sim.read_across_shots("lhs"))
    got_tgt = pack_shots(sim.read_across_shots("target"))
    got_cln = pack_shots(sim.read_across_shots("@clean"))
    phase_flip = np.asarray(sim.read_phase_flipped_across_shots(), dtype=bool)

    expected = (lhs + tgt) & 0xFFFF
    if not np.array_equal(got_lhs, lhs):
        raise AssertionError("lhs register mutated!")
    if not np.array_equal(got_tgt, expected):
        raise AssertionError("target sum mismatch!")
    if np.any(got_cln != 0):
        raise AssertionError("@clean ancillas not returned to |0>!")
    if np.any(phase_flip):
        raise AssertionError(f"HMR phase kickback error on {int(np.sum(phase_flip))}/{shots} shots!")

    print(f"[PASS] {kmb_path.name}: {shots} shots bit-exact & 0 phase-kickback errors (Q={circ.num_qubits}, CCX={stats.num_toffoli_like})")


def fuzz_verify_gf2_12_mul(kmb_path: pathlib.Path, shots: int = 4096, seed: int = 20261116) -> None:
    circ = km.Circuit.from_file(str(kmb_path))
    stats = circ.operation_statistics()
    assert circ.num_qubits == 37, f"Expected 37 qubits (36 IO + 1 ancilla), got {circ.num_qubits}"
    assert stats.num_toffoli_like == 42, f"Expected 42 Chudnovsky P^1(GF(2^3)) Toffolis, got {stats.num_toffoli_like}"

    field = km.GF2Field(12)  # modulus 0x1053 = x^12 + x^6 + x^4 + x + 1
    rng = np.random.default_rng(seed)
    sim = km.Simulator(batch_size=shots)
    sim.use_same_registers_as(circ)

    a = rng.integers(0, 1 << 12, size=shots, dtype=np.uint64)
    b = rng.integers(0, 1 << 12, size=shots, dtype=np.uint64)
    out_init = rng.integers(0, 1 << 12, size=shots, dtype=np.uint64)

    sim.write_across_shots("a", unpack_shots(a, 12))
    sim.write_across_shots("b", unpack_shots(b, 12))
    sim.write_across_shots("out", unpack_shots(out_init, 12))
    sim.write_across_shots("@clean", np.zeros((1, shots), dtype=np.bool_))
    sim.do(circ)

    got_a = pack_shots(sim.read_across_shots("a"))
    got_b = pack_shots(sim.read_across_shots("b"))
    got_out = pack_shots(sim.read_across_shots("out"))
    got_cln = pack_shots(sim.read_across_shots("@clean"))
    phase_flip = np.asarray(sim.read_phase_flipped_across_shots(), dtype=bool)

    expected = np.array([int(out_init[i]) ^ int(field.mul(int(a[i]), int(b[i]))) for i in range(shots)], dtype=np.uint64)

    if not np.array_equal(got_a, a) or not np.array_equal(got_b, b):
        raise AssertionError("GF(2^12) input registers a/b were not cleanly restored!")
    if not np.array_equal(got_out, expected):
        bad = int(np.argmax(got_out != expected))
        raise AssertionError(f"GF(2^12) product mismatch at shot {bad}: got {got_out[bad]:#x}, expected {expected[bad]:#x}")
    if np.any(got_cln != 0):
        raise AssertionError("@clean ancilla not returned to |0>!")
    if np.any(phase_flip):
        raise AssertionError(f"GF(2^12) HMR phase kickback error on {int(np.sum(phase_flip))}/{shots} shots!")

    print(f"[PASS] {kmb_path.name}: {shots} shots GF(2^12) exact & 0 phase errors (Q={circ.num_qubits}, CCX={stats.num_toffoli_like} vs 144 naive)")


def verify_zpow128(kmb_path: pathlib.Path) -> None:
    circ = km.Circuit.from_file(str(kmb_path))
    stats = circ.operation_statistics()
    assert stats["Z_POW"] == 4, f"Expected 4 Z_POW ops, got {stats['Z_POW']}"
    assert stats.num_t_like == 2, f"Expected 2 T-like ops, got {stats.num_t_like}"
    print(f"[PASS] {kmb_path.name}: 128-bit FixedPrecisionAngle128 parsed cleanly (levels={stats.z_pow_clifford_hierarchy_counts})")


def fuzz_verify_gf64_mul(kmb_path: pathlib.Path) -> None:
    circ = km.Circuit.from_file(str(kmb_path))
    stats = circ.operation_statistics()
    assert circ.num_qubits == 19, f"Expected 19 qubits (18 IO + 1 ancilla), got {circ.num_qubits}"
    assert stats.num_toffoli_like == 15, f"Expected 15 Chudnovsky P^1(GF(4)) Toffolis, got {stats.num_toffoli_like}"

    def gf8_mul(a: int, b: int) -> int:
        res = 0
        for _ in range(3):
            if b & 1:
                res ^= a
            hi = a & 4
            a = (a << 1) & 7
            if hi:
                a ^= 3
            b >>= 1
        return res

    def gf64_mul(x: int, y: int) -> int:
        x0, x1 = x & 7, (x >> 3) & 7
        y0, y1 = y & 7, (y >> 3) & 7
        p0 = gf8_mul(x0, y0)
        p1 = gf8_mul(x1, y1)
        pm = gf8_mul(x0 ^ x1, y0 ^ y1)
        return (p0 ^ p1) | ((pm ^ p0) << 3)

    shots = 4096
    sim = km.Simulator(batch_size=shots, seed=20261117)
    sim.use_same_registers_as(circ)
    idx = np.arange(shots, dtype=np.uint64)
    a = idx & 63
    b = (idx >> 6) & 63
    out_init = (idx * 17 + 5) & 63
    sim.write_across_shots("a", unpack_shots(a, 6))
    sim.write_across_shots("b", unpack_shots(b, 6))
    sim.write_across_shots("out", unpack_shots(out_init, 6))
    sim.write_across_shots("@clean", np.zeros((1, shots), dtype=np.bool_))
    sim.do(circ)

    got_a = pack_shots(sim.read_across_shots("a"))
    got_b = pack_shots(sim.read_across_shots("b"))
    got_out = pack_shots(sim.read_across_shots("out"))
    got_cln = pack_shots(sim.read_across_shots("@clean"))
    phase_flip = np.asarray(sim.read_phase_flipped_across_shots(), dtype=bool)
    expected = np.array([int(out_init[i]) ^ gf64_mul(int(a[i]), int(b[i])) for i in range(shots)], dtype=np.uint64)

    assert np.array_equal(got_a, a) and np.array_equal(got_b, b), "GF(64) inputs mutated!"
    assert np.array_equal(got_out, expected), "GF(64) product mismatch!"
    assert not np.any(got_cln != 0), "GF(64) ancilla not returned to |0>!"
    assert not np.any(phase_flip), "GF(64) HMR phase kickback error!"
    print(f"[PASS] {kmb_path.name}: ALL 4096 GF(64) pairs exact & 0 phase errors (Q={circ.num_qubits}, CCX={stats.num_toffoli_like} vs 18 Karatsuba)")


def fuzz_verify_qroam_vslot0(kmb_path: pathlib.Path) -> None:
    circ = km.Circuit.from_file(str(kmb_path))
    stats = circ.operation_statistics()
    repeats = 64
    num_rows = 37
    shots = num_rows * repeats
    sim = km.Simulator(batch_size=shots, seed=20261118)
    sim.use_same_registers_as(circ)
    sel = np.tile(np.arange(num_rows, dtype=np.uint64), repeats)
    tgt_init = (np.arange(shots, dtype=np.uint64) * 0x27D + 0x19) & 0xFFF
    anc_len = len(circ.register_data["@clean"])
    sim.write_across_shots("selection", unpack_shots(sel, 6))
    sim.write_across_shots("target", unpack_shots(tgt_init, 12))
    sim.write_across_shots("@clean", np.zeros((anc_len, shots), dtype=np.bool_))
    sim.do(circ)

    got_sel = pack_shots(sim.read_across_shots("selection"))
    got_tgt = pack_shots(sim.read_across_shots("target"))
    got_cln = pack_shots(sim.read_across_shots("@clean"))
    phase_flip = np.asarray(sim.read_phase_flipped_across_shots(), dtype=bool)
    expected = tgt_init ^ ((sel * 0x5BD + 0x137) & 0x7FF)

    assert np.array_equal(got_sel, sel), "QROAM selection mutated!"
    assert np.array_equal(got_tgt, expected), "QROAM target mismatch!"
    assert not np.any(got_cln != 0), "QROAM ancillas not returned to |0>!"
    assert not np.any(phase_flip), f"QROAM HMR phase kickback error on {int(np.sum(phase_flip))}/{shots} shots!"
    print(f"[PASS] {kmb_path.name}: {shots} shots Virtual-Slot-0 QROAM exact & 0 phase errors (Q={circ.num_qubits}, CCX={stats.num_toffoli_like})")


def fuzz_verify_hmr_lut(kmb_path: pathlib.Path, shots: int = 4096, seed: int = 20261119) -> None:
    circ = km.Circuit.from_file(str(kmb_path))
    stats = circ.operation_statistics()
    assert len(circ.register_data["selection"]) == 6, "Expected 6-qubit selection register"
    assert len(circ.register_data["out"]) == 12, "Expected 12-qubit out register"
    assert stats["HMR"] > 0, "Expected HMR operations in HMR_Dual Parity-LUT circuit"

    rng = np.random.default_rng(seed)
    num_rows = 37
    sel = np.concatenate([
        np.arange(num_rows, dtype=np.uint64),
        rng.integers(0, num_rows, size=shots - num_rows, dtype=np.uint64),
    ])
    out_init = rng.integers(0, 1 << 12, size=shots, dtype=np.uint64)
    anc_len = len(circ.register_data["@clean"])

    sim = km.Simulator(batch_size=shots, seed=seed)
    sim.use_same_registers_as(circ)
    sim.write_across_shots("selection", unpack_shots(sel, 6))
    sim.write_across_shots("out", unpack_shots(out_init, 12))
    sim.write_across_shots("@clean", np.zeros((anc_len, shots), dtype=np.bool_))
    sim.do(circ)

    got_sel = pack_shots(sim.read_across_shots("selection"))
    got_out = pack_shots(sim.read_across_shots("out"))
    got_cln = pack_shots(sim.read_across_shots("@clean"))
    phase_flip = np.asarray(sim.read_phase_flipped_across_shots(), dtype=bool)
    expected = out_init ^ ((sel * 0x5BD + 0x137) & 0xFFF)

    assert np.array_equal(got_sel, sel), "HMR_Dual LUT selection mutated!"
    assert np.array_equal(got_out, expected), "HMR_Dual LUT output mismatch!"
    assert not np.any(got_cln != 0), "HMR_Dual LUT @clean target/ancillas not returned to |0>!"
    assert not np.any(phase_flip), f"HMR_Dual LUT phase kickback error on {int(np.sum(phase_flip))}/{shots} shots!"
    print(f"[PASS] {kmb_path.name}: {shots} shots HMR_Dual Parity-LUT exact & 0 phase errors (Q={circ.num_qubits}, CCX={stats.num_toffoli_like}, HMR={stats['HMR']})")


def fuzz_verify_linf2(kmb_path: pathlib.Path, seed: int = 20261120) -> None:
    circ = km.Circuit.from_file(str(kmb_path))
    stats = circ.operation_statistics()
    assert circ.num_qubits == 12, f"Expected 12 qubits, got {circ.num_qubits}"
    assert stats.num_toffoli_like == 0, f"Expected 0 Toffolis, got {stats.num_toffoli_like}"
    assert "@clean" not in circ.register_data, "Expected 0 @clean ancillas for in-place GL_12(F_2)"

    def rows_to_mat(rows: list[int], w: int = 12) -> np.ndarray:
        return np.array([[(row >> c) & 1 for c in range(w)] for row in rows], dtype=np.uint8)

    l1_rows, u1_rows, l2_rows, u2_rows = [], [], [], []
    for r in range(12):
        lo_mask = (1 << r) - 1
        hi_mask = 0xFFF ^ ((1 << (r + 1)) - 1)
        l1_rows.append((1 << r) | (((r + 1) * 0x5AD) & lo_mask))
        u1_rows.append((1 << r) | (((r + 3) * 0x36B) & hi_mask))
        l2_rows.append((1 << r) | (((r + 5) * 0x6D3) & lo_mask))
        u2_rows.append((1 << r) | (((r + 2) * 0x59D) & hi_mask))

    m1 = (rows_to_mat(l1_rows) @ rows_to_mat(u1_rows)) & 1
    m2 = (rows_to_mat(u2_rows) @ rows_to_mat(l2_rows)) & 1
    m21 = (m2 @ m1) & 1
    assert not np.array_equal(m21, (m1 @ m2) & 1), "Expected non-commuting GL_12(F_2) test matrices"

    shots = 4096
    x = np.arange(shots, dtype=np.uint64)
    sim = km.Simulator(batch_size=shots, seed=seed)
    sim.use_same_registers_as(circ)
    sim.write_across_shots("x", unpack_shots(x, 12))
    sim.do(circ)

    got_x = pack_shots(sim.read_across_shots("x"))
    phase_flip = np.asarray(sim.read_phase_flipped_across_shots(), dtype=bool)
    x_bits = unpack_shots(x, 12).astype(np.uint8)
    expected = pack_shots(((m21 @ x_bits) & 1).astype(np.bool_))

    assert np.array_equal(got_x, expected), "GL_12(F_2) M2*M1 linear output mismatch!"
    assert not np.any(phase_flip), f"GL_12(F_2) phase kickback error on {int(np.sum(phase_flip))}/{shots} shots!"
    print(f"[PASS] {kmb_path.name}: ALL {shots} GL_12(F_2) states exact & 0 phase errors (Q={circ.num_qubits}, CCX={stats.num_toffoli_like})")


def fuzz_verify_cuccaro(kmb_path: pathlib.Path, shots: int = 4096, seed: int = 20261121) -> None:
    circ = km.Circuit.from_file(str(kmb_path))
    stats = circ.operation_statistics()
    assert circ.num_qubits == 34, f"Expected 34 qubits (1 ctrl + 16 lhs + 16 target + 1 @clean), got {circ.num_qubits}"
    assert len(circ.register_data["@clean"]) == 1, f"Expected 1 @clean ancilla, got {len(circ.register_data['@clean'])}"
    assert stats.num_toffoli_like == 46, f"Expected 46 Toffolis (3*16 - 2), got {stats.num_toffoli_like}"

    rng = np.random.default_rng(seed)
    ctrl = rng.integers(0, 2, size=shots, dtype=np.uint64)
    lhs = rng.integers(0, 1 << 16, size=shots, dtype=np.uint64)
    tgt = rng.integers(0, 1 << 16, size=shots, dtype=np.uint64)

    sim = km.Simulator(batch_size=shots, seed=seed)
    sim.use_same_registers_as(circ)
    sim.write_across_shots("ctrl", unpack_shots(ctrl, 1))
    sim.write_across_shots("lhs", unpack_shots(lhs, 16))
    sim.write_across_shots("target", unpack_shots(tgt, 16))
    sim.write_across_shots("@clean", np.zeros((1, shots), dtype=np.bool_))
    sim.do(circ)

    got_ctrl = pack_shots(sim.read_across_shots("ctrl"))
    got_lhs = pack_shots(sim.read_across_shots("lhs"))
    got_tgt = pack_shots(sim.read_across_shots("target"))
    got_cln = pack_shots(sim.read_across_shots("@clean"))
    phase_flip = np.asarray(sim.read_phase_flipped_across_shots(), dtype=bool)
    expected = np.where(ctrl == 1, (tgt + lhs) & 0xFFFF, tgt)

    assert np.array_equal(got_ctrl, ctrl), "Cuccaro ctrl mutated!"
    assert np.array_equal(got_lhs, lhs), "Cuccaro lhs mutated!"
    assert np.array_equal(got_tgt, expected), "Cuccaro controlled sum mismatch!"
    assert not np.any(got_cln != 0), "Cuccaro @clean ancilla not returned to |0>!"
    assert not np.any(phase_flip), f"Cuccaro phase kickback error on {int(np.sum(phase_flip))}/{shots} shots!"
    print(f"[PASS] {kmb_path.name}: {shots} shots Cuccaro 1-ancilla ctrl-iadd exact & 0 phase errors (Q={circ.num_qubits}, CCX={stats.num_toffoli_like})")


def inspect_data_dir(data_dir: pathlib.Path) -> int:
    if not data_dir.exists() or not any(data_dir.iterdir()):
        print(f"[INFO] {data_dir} is currently empty. Download Kaggle competition files via `make data` or drop the folder into {data_dir}/.")
        return 0

    print(f"=== Inspecting Competition Data in {data_dir} ===")
    for p in sorted(data_dir.rglob("*")):
        if p.is_dir():
            continue
        sz = p.stat().st_size
        rel = p.relative_to(data_dir)
        if p.suffix in {".kmb", ".kmx"}:
            c = km.Circuit.from_file(str(p))
            st = c.operation_statistics()
            c1 = c.num_qubits * (st.num_toffoli_like + st.num_t_like / 4.0)
            print(f"  [CIRCUIT] {rel} ({sz:,} B) -> Q={c.num_qubits}, Toffoli={st.num_toffoli_like:,}, T={st.num_t_like:,}, Q*C={c1:,.1f}")
        elif p.suffix == ".npz":
            with np.load(p, allow_pickle=False) as z:
                shapes = ", ".join(f"{k}:{z[k].shape}:{z[k].dtype}" for k in z.files)
            print(f"  [NPZ]     {rel} ({sz:,} B) -> {shapes}")
        else:
            print(f"  [FILE]    {rel} ({sz:,} B)")
    return 0


def run_official_score(sub_dir: pathlib.Path, verifiers_dir: pathlib.Path, full: bool = False) -> int:
    sub_dir = sub_dir.resolve()
    verifiers_dir = verifiers_dir.resolve()
    inv_f = sub_dir / "invsqrt.kmb"
    dqi_f = sub_dir / "dqi.kmb"
    thc_d = sub_dir / "dfthcsos"
    for p in (inv_f, dqi_f, thc_d):
        if not p.exists():
            print(f"[ERROR] Missing {p}", file=sys.stderr)
            return 1

    py = sys.executable

    # 1. invsqrt (3,133,452 = 196 * 15,987)
    inv_cmd = [
        py,
        str(verifiers_dir / "invsqrt" / "verify_submission_invsqrt.py"),
        str(inv_f.resolve()),
    ] + (["--shots-per-x", "64", "--batch-size", "256"] if full else ["--shots-per-x", "1", "--batch-size", "196"])
    r1 = subprocess.run(inv_cmd, cwd=str(verifiers_dir), capture_output=True, text=True)
    if r1.returncode != 0:
        print(f"[FAIL] invsqrt: {r1.stdout}\n{r1.stderr}", file=sys.stderr)
        return 1
    j1 = json.loads(r1.stdout)
    if j1.get("status") != "PASS":
        print(f"[FAIL] invsqrt: {r1.stdout}", file=sys.stderr)
        return 1
    cost_1 = float(j1["score"])

    # 2. dqi
    if full:
        dqi_cmd = [
            py,
            str(verifiers_dir / "dqi" / "verify_submission_dqi.py"),
            str(dqi_f.resolve()),
            "--max-workers", "2",
        ]
    else:
        dqi_cmd = [
            py,
            str(verifiers_dir / "dqi" / "verify_submission_dqi.py"),
            str(dqi_f.resolve()),
            "--num-shots", "256",
            "--batch-size", "256",
            "--max-workers", "1",
            "--log2-delta", "-2",
        ]
    r2 = subprocess.run(dqi_cmd, cwd=str(verifiers_dir), capture_output=True, text=True)
    if r2.returncode != 0:
        print(f"[FAIL] dqi: {r2.stdout}\n{r2.stderr}", file=sys.stderr)
        return 1
    j2 = json.loads(r2.stdout)
    if j2.get("status") != "PASS":
        print(f"[FAIL] dqi: {r2.stdout}", file=sys.stderr)
        return 1
    if full:
        cost_2 = float(j2["score"])
    else:
        # Project empirical error onto the competition 102,400-shot log2_delta=-128 Clopper-Pearson curve
        if str(verifiers_dir) not in sys.path:
            sys.path.insert(0, str(verifiers_dir))
        from dqi.verify_submission_dqi import clopper_pearson_upper_bound, repetition_overhead

        emp_err = float(j2["empirical_error_rate"])
        k_proj = int(round(emp_err * 102400.0))
        eps_plus = clopper_pearson_upper_bound(k_proj, 102400, -128.0)
        if eps_plus >= 1.0 / 16.0:
            print(f"[FAIL] dqi projected eps_plus={eps_plus:.4f} >= 1/16", file=sys.stderr)
            return 1
        rep_102k = repetition_overhead(eps_plus)
        cost_2 = rep_102k * (float(j2["num_qubits"]) ** 2) * float(j2["toffoli_count"])

    # 3. dfthcsos
    thc_cmd = [
        py,
        str(verifiers_dir / "dfthcsos" / "verify.py"),
        str(thc_d.resolve()),
        "--shots-per-address",
        "64" if full else "4",
    ]
    r3 = subprocess.run(thc_cmd, cwd=str(verifiers_dir), capture_output=True, text=True)
    if r3.returncode != 0:
        print(f"[FAIL] dfthcsos: {r3.stdout}\n{r3.stderr}", file=sys.stderr)
        return 1
    j3 = json.loads(r3.stdout)
    if j3.get("status") != "PASS":
        print(f"[FAIL] dfthcsos: {r3.stdout}", file=sys.stderr)
        return 1
    cost_3 = float(j3["score"])

    s1 = REF_COST_INVSQRT / cost_1
    s2 = REF_COST_DQI / cost_2
    s3 = REF_COST_DFTHCSOS / cost_3
    geomean = (s1 * s2 * s3) ** (1.0 / 3.0)

    mode_str = "FULL (Competition Shots)" if full else "FAST (Exact invsqrt + 256-shot DQI + 4-shot DFTHCSOS)"
    print(f"=== Official Verifier Leaderboard Summary [{mode_str}] ===")
    print(f"  [1] invsqrt  : PASS | cost = {cost_1:>18,.2f} | ref = {REF_COST_INVSQRT:>18,.2f} | gain = {s1:7.4f}x (Q={j1['num_qubits']}, C={j1['toffoli_equivalent']:.1f})")
    print(f"  [2] dqi      : PASS | cost = {cost_2:>18,.2f} | ref = {REF_COST_DQI:>18,.2f} | gain = {s2:7.4f}x (Q={j2['num_qubits']}, T={j2['toffoli_count']:,.1f})")
    print(f"  [3] dfthcsos : PASS | cost = {cost_3:>18,.2f} | ref = {REF_COST_DFTHCSOS:>18,.2f} | gain = {s3:7.4f}x (Q={j3['num_qubits']}, T={j3['toffoli_count']:,.1f}, lam={j3['lambda_eff']:.3f})")
    print(f"  -------------------------------------------------------------------------------------")
    print(f"  >>> GEOMETRIC MEAN SCORE = {geomean:.6f}x <<<")
    return 0


def verify_kmb_generic(kmb_path: pathlib.Path, shots: int = 2048) -> None:
    import re
    circ = km.Circuit.from_file(str(kmb_path))
    stats = circ.operation_statistics()
    regs = [(name, tgts) for name, tgts in circ.register_data.items() if not name.startswith("@")]
    q_regs = [(name, tgts) for name, tgts in regs if all(isinstance(e, km.q) for e in tgts)]

    # Auto-detect exact (input %name ...) list if sibling/examples .qir exists; else schema / fallback
    init_regs = None
    qir_backed = False
    for cand_qir in [
        kmb_path.with_suffix(".qir"),
        kmb_path.parent.parent / "examples" / f"{kmb_path.stem}.qir",
    ]:
        if cand_qir.exists():
            qir_ins = set(re.findall(r"\(\s*input\s+%([^\s\)]+)", cand_qir.read_text()))
            if qir_ins:
                init_regs = [(n, t) for n, t in q_regs if n in qir_ins]
                qir_backed = True
                break
    if not init_regs:
        q_names = [n for n, _ in q_regs]
        if "selection" in q_names and all(n == "selection" or n.startswith("word_") for n in q_names):
            init_regs = [(n, t) for n, t in q_regs if n == "selection"]
        else:
            init_regs = q_regs[:-1] if len(q_regs) > 1 else q_regs

    s1 = km.Simulator(batch_size=shots, seed=10101)
    s2 = km.Simulator(batch_size=shots, seed=20202)
    s1.use_same_registers_as(circ)
    s2.use_same_registers_as(circ)

    rng = np.random.default_rng(20261116)
    for rname, targets in init_regs:
        if len(targets) == 0:
            continue
        bits = rng.integers(0, 2, size=(len(targets), shots), dtype=np.bool_)
        s1.write_across_shots(rname, bits)
        s2.write_across_shots(rname, bits)

    s1.do(circ)
    s2.do(circ)

    try:
        pf1_arr = np.asarray(s1.read_phase_flipped_across_shots(), dtype=bool)
        pf2_arr = np.asarray(s2.read_phase_flipped_across_shots(), dtype=bool)
    except ValueError:
        ph1 = [s1.read_shot_phase(i) for i in range(shots)]
        ph2 = [s2.read_shot_phase(i) for i in range(shots)]
        drift_ph = sum(1 for a, b in zip(ph1, ph2) if a != b)
        if drift_ph != 0:
            raise AssertionError(f"{kmb_path.name}: HMR dual-seed phase drift on {drift_ph}/{shots} shots")
    else:
        pf1 = int(np.sum(pf1_arr))
        pf2 = int(np.sum(pf2_arr))
        if pf1 != 0 or pf2 != 0:
            raise AssertionError(f"{kmb_path.name}: Kickmix phase flip error pf1={pf1}, pf2={pf2}")

    all_reg_qs = {q.id for _, tgts in regs for q in tgts if isinstance(q, km.q)}
    scratch_qs = [km.q(i) for i in range(circ.num_qubits) if i not in all_reg_qs]
    if scratch_qs:
        sc_arr = km.array(scratch_qs)
        sc1 = s1.read_across_shots(sc_arr)
        sc2 = s2.read_across_shots(sc_arr)
        sc_drift = int(np.sum(sc1 != sc2))
        if sc_drift != 0:
            raise AssertionError(f"{kmb_path.name}: HMR dual-seed drift on {sc_drift} scratch qubit bits")
        if not qir_backed:
            for label, sc_bits in (("seed1", sc1), ("seed2", sc2)):
                dirty = sc_bits.any(axis=1)
                if dirty.any():
                    bad_q = scratch_qs[int(np.argmax(dirty))]
                    raise AssertionError(f"{kmb_path.name}: scratch/@clean qubit {bad_q} did not return to |0> ({label})")

    seed_drift = 0
    for rname, _ in regs:
        b1 = s1.read_across_shots(rname)
        b2 = s2.read_across_shots(rname)
        seed_drift += int(np.sum(b1 != b2))
    if seed_drift != 0:
        raise AssertionError(f"{kmb_path.name}: HMR dual-seed drift on {seed_drift} register bits")

    print(
        f"[KICKMIX PASS] {kmb_path.name}: Q={circ.num_qubits} B={circ.num_bits} "
        f"CCX={stats.num_toffoli_like} Z_POW={stats['Z_POW']} T_like={stats.num_t_like} "
        f"HMR={stats['HMR']} D_react={circ.reaction_depth()} ({shots} dual-seed shots, 0 phase/drift)"
    )


def export_dfthc_words(sub_dir: pathlib.Path, out_words_bin: pathlib.Path) -> int:
    verifiers_root = pathlib.Path(__file__).resolve().parent.parent / "data" / "google-verifiable-quantum-circuit-compilation" / "verifiers"
    if str(verifiers_root) not in sys.path:
        sys.path.insert(0, str(verifiers_root))
    from dfthcsos.scoreinfra.circuits import load_circuits, recover_words_from_circuits

    npz_f = sub_dir / "submission.npz"
    if not npz_f.is_file():
        print(f"[ERROR] Missing {npz_f}", file=sys.stderr)
        return 1
    with np.load(npz_f, allow_pickle=False) as z:
        rank_r, bases_b, ensemble_c = (int(x) for x in z["W"].shape)

    circuits = load_circuits(sub_dir)
    num_vectors = 76 + rank_r * bases_b
    words, layout, _ = recover_words_from_circuits(
        circuits, num_vectors=num_vectors, num_angles=75, shots_per_address=1, seed=0
    )
    if words.shape != (num_vectors, 75):
        print(f"[ERROR] Unexpected recovered shape {words.shape}, expected {(num_vectors, 75)}", file=sys.stderr)
        return 1
    out_words_bin.parent.mkdir(parents=True, exist_ok=True)
    payload = struct.pack("<4i", rank_r, bases_b, ensemble_c, layout.brot) + words.astype("<u2").tobytes()
    out_words_bin.write_bytes(payload)
    print(f"[OK] Exported {out_words_bin} (R={rank_r}, B={bases_b}, C={ensemble_c}, brot={layout.brot}, rows={num_vectors})")
    return 0


def pack_submission(stage_dir: pathlib.Path, out_zip: pathlib.Path) -> int:
    print(f"=== Packaging Submission from {stage_dir} -> {out_zip} ===")
    inv_f = stage_dir / "invsqrt.kmb"
    dqi_f = stage_dir / "dqi.kmb"
    thc_d = stage_dir / "dfthcsos"
    npz_f = thc_d / "submission.npz"
    if not (inv_f.exists() and dqi_f.exists() and thc_d.is_dir() and npz_f.exists()):
        print(
            f"[ERROR] Expected official layout in {stage_dir}: invsqrt.kmb, dqi.kmb, dfthcsos/submission.npz, dfthcsos/angle_k_*.kmb",
            file=sys.stderr,
        )
        return 1

    # Validate invsqrt register contract
    inv_c = km.Circuit.from_file(str(inv_f))
    inv_regs = inv_c.register_data
    if len(inv_regs) != inv_c.num_registers or list(inv_regs) != ["x"]:
        print(f"[REJECTED] invsqrt.kmb invalid register schema: {list(inv_regs)}", file=sys.stderr)
        return 1
    x_reg = list(inv_regs["x"])
    if len(x_reg) != 22 or not all(isinstance(e, km.q) for e in x_reg) or len(set(x_reg)) != 22:
        print("[REJECTED] invsqrt.kmb 'x' must hold 22 distinct qubits", file=sys.stderr)
        return 1

    # Validate dqi register contract & no Z_POW
    dqi_c = km.Circuit.from_file(str(dqi_f))
    if dqi_c.operation_statistics()["Z_POW"] != 0:
        print("[REJECTED] dqi.kmb contains Z_POW ops (forbidden in dqi)!", file=sys.stderr)
        return 1
    dqi_regs = dqi_c.register_data
    s_names = [f"S{i}" for i in range(70)]
    e_names = [f"E{j}" for j in range(4095)]
    if len(dqi_regs) != dqi_c.num_registers or set(dqi_regs) != set(s_names) | set(e_names):
        print("[REJECTED] dqi.kmb must declare exact registers S0..S69 and E0..E4094!", file=sys.stderr)
        return 1
    s_all = [e for name in s_names for e in dqi_regs[name]]
    e_all = [e for name in e_names for e in dqi_regs[name]]
    if len(s_all) != 70 * 12 or not all(isinstance(e, km.q) for e in s_all) or len(set(s_all)) != 70 * 12:
        print("[REJECTED] dqi.kmb S0..S69 must each hold 12 distinct qubits!", file=sys.stderr)
        return 1
    if len(e_all) != 4095 * 12 or not all(isinstance(e, km.b) for e in e_all) or len(set(e_all)) != 4095 * 12:
        print("[REJECTED] dqi.kmb E0..E4094 must each hold 12 distinct classical bits!", file=sys.stderr)
        return 1

    # Validate dfthcsos submission.npz keys & angle_k_*.kmb contract
    with np.load(npz_f, allow_pickle=False) as z:
        req_keys = {"U", "W", "H2_I", "H1_I", "Hsym_lower", "tree_edges"}
        missing = req_keys - set(z.files)
        if missing:
            print(f"[REJECTED] {npz_f} missing required keys: {sorted(missing)}", file=sys.stderr)
            return 1

    angle_files = sorted(thc_d.glob("angle_k_*.kmb"), key=lambda p: int(p.stem.split("_")[-1]))
    if not angle_files:
        print(f"[ERROR] No angle_k_*.kmb circuits found in {thc_d}", file=sys.stderr)
        return 1

    verifiers_root = pathlib.Path(__file__).resolve().parent.parent / "data" / "google-verifiable-quantum-circuit-compilation" / "verifiers"
    if str(verifiers_root) not in sys.path:
        sys.path.insert(0, str(verifiers_root))
    from dfthcsos.scoreinfra.circuits import InvalidCircuits, register_layout

    thc_circs = []
    for idx, af in enumerate(angle_files):
        if af.name != f"angle_k_{idx}.kmb":
            print(f"[REJECTED] Contiguous angle_k_<idx>.kmb required; expected angle_k_{idx}.kmb, got {af.name}", file=sys.stderr)
            return 1
        c = km.Circuit.from_file(str(af))
        st = c.operation_statistics()
        if c.num_qubits > 128:
            print(f"[REJECTED] {af.name} uses {c.num_qubits} qubits > 128 hard limit!", file=sys.stderr)
            return 1
        if st["Z_POW"] != 0:
            print(f"[REJECTED] {af.name} contains {st['Z_POW']} Z_POW ops (forbidden in dfthcsos)!", file=sys.stderr)
            return 1
        thc_circs.append(c)
    try:
        register_layout(thc_circs)
    except InvalidCircuits as exc:
        print(f"[REJECTED] dfthcsos register contract failed: {exc}", file=sys.stderr)
        return 1

    out_zip.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(out_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for f in [inv_f, dqi_f, npz_f] + angle_files:
            arcname = str(f.relative_to(stage_dir))
            zf.write(f, arcname=arcname)
    digest = hashlib.sha256(out_zip.read_bytes()).hexdigest()
    print(f"[OK] Wrote {out_zip} ({out_zip.stat().st_size:,} B, sha256={digest[:16]}...)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="Kickmix local verifier, official score runner & Kaggle packager")
    ap.add_argument("--selftest", metavar="BUILD_DIR", help="Fuzz-verify C17-synthesized selftest .kmb circuits")
    ap.add_argument("--kmb", nargs="+", metavar="FILE.kmb", help="Cross-verify arbitrary .kmb v2 circuits via Kickmix C++ simulator")
    ap.add_argument("--inspect-data", metavar="DATA_DIR", help="Inspect downloaded Kaggle baseline circuits & .npz")
    ap.add_argument("--score", nargs=2, metavar=("SUBMISSION_DIR", "VERIFIERS_DIR"), help="Run official verifiers and compute Leaderboard Score")
    ap.add_argument("--full", action="store_true", help="Use full competition shot counts (200M invsqrt, 102.4k dqi, 64/addr dfthcsos)")
    ap.add_argument("--pack", nargs=2, metavar=("STAGE_DIR", "OUT_ZIP"), help="Validate and zip submission bundle")
    ap.add_argument("--export-dfthc-words", nargs=2, metavar=("SUB_DIR", "OUT_WORDS_BIN"), help="Recover DF-THC words.bin from submission directory")
    args = ap.parse_args()

    if args.kmb:
        for p in args.kmb:
            verify_kmb_generic(pathlib.Path(p))
    if args.selftest:
        bdir = pathlib.Path(args.selftest)
        fuzz_verify_adder16(bdir / "selftest_adder16.kmb")
        fuzz_verify_gf2_12_mul(bdir / "selftest_gf2_12_mul.kmb")
        fuzz_verify_gf64_mul(bdir / "selftest_gf64_mul.kmb")
        fuzz_verify_qroam_vslot0(bdir / "selftest_qroam_vslot0.kmb")
        verify_zpow128(bdir / "selftest_zpow128.kmb")
        fuzz_verify_hmr_lut(bdir / "selftest_hmr_lut.kmb")
        fuzz_verify_linf2(bdir / "selftest_linf2.kmb")
        fuzz_verify_cuccaro(bdir / "selftest_cuccaro.kmb")
    if args.inspect_data:
        return inspect_data_dir(pathlib.Path(args.inspect_data))
    if args.score:
        return run_official_score(pathlib.Path(args.score[0]), pathlib.Path(args.score[1]), full=args.full)
    if args.pack:
        return pack_submission(pathlib.Path(args.pack[0]), pathlib.Path(args.pack[1]))
    if args.export_dfthc_words:
        return export_dfthc_words(pathlib.Path(args.export_dfthc_words[0]), pathlib.Path(args.export_dfthc_words[1]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
