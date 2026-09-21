# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""relabel_radial_ma_shells.py (supervisor decision 56): the byte-level radial-shell
relabel of a MSH 2.2 binary mesh on a synthetic strip (ring-aligned quadrangles on the
metal top face and sidewall of one edge, corner triangles, a far triangle, non-MA
surfaces and a volume element), the metal edge lines of the two-edge test data, and the
qualify plumbing that consumes the census: case_inputs.expand_radial_shells (one MA
Dielectric entry per shell ordinal, the parent labels replaced everywhere),
physics_run_parameters with the labeled LinearTol deviation, compare_matrices /
p_sequence with a reference labeled by its own interface map (per-type sums over the
shells reproduce the whole-MA participation)."""
import copy
import csv
import json
from pathlib import Path
import struct
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "qualify"))
import case_inputs  # noqa: E402
import compare_matrices  # noqa: E402
import p_sequence  # noqa: E402
import radial_ma_profile  # noqa: E402
import relabel_radial_ma_shells as relabel  # noqa: E402

RADII = [0.00025 * (2 ** k - 1) for k in range(1, 8)]
TESTDATA = HERE / "testdata" / "two-edge-8dd4bc70f183"


def msh22(names, nodes, elements):
    """A MSH 2.2 binary file: names [(dim, tag, name)], nodes {tag: xyz}, elements
    [(type, physical, elementary, nodes)] in one-element blocks (the publisher's layout)."""
    out = bytearray(b"$MeshFormat\n2.2 1 8\n" + struct.pack("<i", 1) + b"\n$EndMeshFormat\n")
    out += ("$PhysicalNames\n" + str(len(names)) + "\n" + "".join(f'{d} {t} "{n}"\n' for d, t, n in names)
            + "$EndPhysicalNames\n").encode()
    out += f"$Nodes\n{len(nodes)}\n".encode()
    for tag, xyz in nodes.items():
        out += struct.pack("<i3d", tag, *xyz)
    out += f"\n$EndNodes\n$Elements\n{len(elements)}\n".encode()
    for tag, (element_type, physical, elementary, connectivity) in enumerate(elements, start=1):
        out += struct.pack("<3i", element_type, 1, 2)
        out += struct.pack(f"<{3 + len(connectivity)}i", tag, physical, elementary, *connectivity)
    out += b"\n$EndElements\n"
    return bytes(out)


def strip_mesh():
    """MA label 6001 around the metal edge x = 0 (line (0, 0)-(0, 1), top z = 0.1, bottom
    z = 0): top-face quads spanning each ring interval [r_(k-1), r_k) at x < 0, sidewall
    quads below the top edge and above the bottom edge, two corner triangles inside
    ring 1 and ring 3, one far top-face triangle; a 3100 triangle and a tetrahedron."""
    nodes, elements, tag = {}, [], 0

    def node(xyz):
        nonlocal tag
        tag += 1
        nodes[tag] = xyz
        return tag

    bounds = [0.0] + RADII
    expected = {}
    for k in range(len(RADII)):
        a, b = bounds[k], bounds[k + 1]
        quad = [node((-a, 0.0, 0.1)), node((-b, 0.0, 0.1)), node((-b, 1.0, 0.1)), node((-a, 1.0, 0.1))]
        elements.append((3, 6001, 7, quad))
        expected[len(elements)] = ("top", k + 1)
        quad = [node((0.0, 0.0, 0.1 - a)), node((0.0, 0.0, 0.1 - b)), node((0.0, 1.0, 0.1 - b)), node((0.0, 1.0, 0.1 - a))]
        elements.append((3, 6001, 8, quad))
        expected[len(elements)] = ("top", k + 1)
        quad = [node((0.0, 0.0, a)), node((0.0, 0.0, b)), node((0.0, 1.0, b)), node((0.0, 1.0, a))]
        elements.append((3, 6001, 8, quad))
        expected[len(elements)] = ("bottom", k + 1)
    elements.append((2, 6001, 7, [node((-0.0001, 1.0, 0.1)), node((-0.0002, 1.0, 0.1)), node((-0.00015, 1.0001, 0.1))]))
    expected[len(elements)] = ("top", 1)
    elements.append((2, 6001, 7, [node((-0.001, 1.0, 0.1)), node((-0.0015, 1.0, 0.1)), node((-0.0012, 1.0004, 0.1))]))
    expected[len(elements)] = ("top", 3)
    elements.append((2, 6001, 7, [node((-0.5, 0.0, 0.1)), node((-0.9, 0.0, 0.1)), node((-0.7, 1.0, 0.1))]))
    expected[len(elements)] = ("far", 0)
    elements.append((2, 3100, 9, [node((0.5, 0.0, -0.05)), node((0.9, 0.0, -0.05)), node((0.7, 1.0, -0.05))]))
    elements.append((4, 2, 11, [node((0.5, 0.5, 0.2)), node((0.6, 0.5, 0.2)), node((0.5, 0.6, 0.2)), node((0.5, 0.5, 0.3))]))
    names = [(2, 1, "matching_surface"), (2, 3100, "surface_3100"), (2, 6001, "surface_6001"), (3, 2, "vacuum")]
    return msh22(names, nodes, elements), expected


LINES = [{"Kind": "top", "First": [0.0, 0.0], "Last": [0.0, 1.0], "Z": 0.1},
         {"Kind": "bottom", "First": [0.0, 0.0], "Last": [0.0, 1.0], "Z": 0.0}]


class RadialShellRelabelTest(unittest.TestCase):
    def test_strip_relabel_is_label_only_and_ring_aligned(self):
        data, expected = strip_mesh()
        out, mesh, shells, relabeled, ma_labels = relabel.relabel(data, lines=LINES, radii=RADII)
        self.assertEqual(ma_labels, [6001])
        label_only = relabel.assert_label_only(data, out, mesh, relabeled)
        self.assertEqual(label_only["RelabeledElements"], len(expected))
        self.assertLessEqual(label_only["DifferingBytesAfterPhysicalNames"], 8 * len(expected))
        after = relabel.read_msh22_binary(out)
        self.assertEqual(data[mesh["NodesSpan"][0]:mesh["NodesSpan"][1]], out[after["NodesSpan"][0]:after["NodesSpan"][1]])
        by_ordinal = {}
        for element_tag, (kind, ring) in expected.items():
            physical = after["Elements"][element_tag - 1][2][0]
            ordinal = physical // relabel.SHELL_LABEL_STRIDE
            self.assertEqual(physical % relabel.SHELL_LABEL_STRIDE, 6001)
            self.assertEqual(relabel.ordinal_description(ordinal, RADII)[:2], (kind, ring))
            by_ordinal.setdefault(ordinal, []).append(element_tag)
        # Every top ring holds its two quads (and the corner triangles of rings 1 / 3),
        # every bottom ring one quad, the far shell the one far triangle.
        rings = len(RADII)
        for k in range(1, rings + 1):
            self.assertEqual(len(by_ordinal[1 + k]), 2 + (k in (1, 3)))
            self.assertEqual(len(by_ordinal[1 + rings + k]), 1)
        self.assertEqual(len(by_ordinal[relabel.FAR_ORDINAL]), 1)
        # Non-MA elements and the volume element keep their tags; the parent name is gone,
        # the shell names present, the other names kept.
        self.assertEqual(after["Elements"][-1][2], (2, 11))
        self.assertEqual(after["Elements"][-2][2], (3100, 9))
        names = {(d, t): n for d, t, n in after["PhysicalNames"]}
        self.assertNotIn((2, 6001), names)
        self.assertEqual(names[(2, 26001)], "ma_shell_top_1_6001")
        self.assertEqual(names[(2, 16001)], "ma_shell_far_0_6001")
        self.assertEqual(names[(2, 96001)], "ma_shell_bottom_1_6001")
        self.assertEqual(names[(3, 2)], "vacuum")
        # Areas: ring k top = 2 x h_k x 1 um (+ corner triangles), bottom = h_k; the shells sum
        # to the parent area; ring-aligned quads never straddle, the 6001 closure is exact.
        sizes = [RADII[0]] + [b - a for a, b in zip(RADII, RADII[1:])]
        for k, h in enumerate(sizes, start=1):
            shell = shells[(1 + rings + k, 6001)]
            self.assertAlmostEqual(shell["Area"], h, places=15)
            self.assertEqual(shell["StraddlingMeasure"], 0.0)
            top = shells[(1 + k, 6001)]
            self.assertGreaterEqual(top["Area"], 2 * h - 1e-15)
        total = sum(shell["Area"] for shell in shells.values())
        # 2 R (top-face + upper sidewall quads) + R (lower sidewall quads) + the far triangle
        # (0.4 x 1 / 2) + the two corner triangles (1e-4 x 1e-4 / 2, 5e-4 x 4e-4 / 2).
        self.assertAlmostEqual(total, 3 * RADII[-1] + 0.2 + 0.5e-8 + 1e-7, places=14)
        closure = relabel.closure_check(shells)
        self.assertEqual(closure["RelativeClosure"], 0.0)

    def test_metal_edge_lines_of_the_two_edge_case(self):
        boundary = relabel.read_csv_rows(TESTDATA / "plan-view-boundary.csv")
        signature = relabel.read_csv_rows(TESTDATA / "mesh-signature.csv")
        lines = relabel.metal_edge_lines(boundary, signature, 0.1)
        # Six Physical segments (the two cut sides on the box are Continuation) x top / bottom.
        self.assertEqual(len(lines), 12)
        self.assertEqual({line["Z"] for line in lines if line["Kind"] == "top"}, {0.1})
        self.assertEqual({line["Z"] for line in lines if line["Kind"] == "bottom"}, {0.0})
        segments = {(tuple(line["First"]), tuple(line["Last"])) for line in lines if line["Kind"] == "top"}
        self.assertIn(((-1.0, 0.0), (-1.0, 1.0)), segments)
        self.assertIn(((4.0, 1.0), (0.0, 1.0)), segments)
        self.assertNotIn(((-5.0, 1.0), (-5.0, 0.0)), segments)
        # A point 0.3 nm below the top edge of conductor 2's y = 0 line is in top ring 2.
        import numpy as np
        distance, kind = relabel.edge_distances(np.asarray([[2.0, 0.0, 0.1 - 0.0003], [2.0, 0.0, 0.0001]]), lines)
        self.assertEqual([relabel.shell_ordinal(d, k, RADII) for d, k in zip(distance, kind)], [3, 1 + 7 + 1])

    def test_expand_radial_shells_and_physics_run_deviation(self):
        shells = {"Shells": [{"Label": 10000 * o + p, "Parent": p, "Ordinal": o, "Kind": "far" if o == 1 else "top", "Ring": max(o - 1, 0),
                              "InnerRadius": 0.0, "OuterRadius": None, "Area": 1.0}
                             for o in (1, 2, 3) for p in (6001, 6002)]}
        config = {"Boundaries": {"Ground": {"Attributes": [5001, 6001, 5002, 6002]},
                                 "PrescribedPotential": [{"Index": 1, "Attributes": [1], "DataFile": "a.csv"},
                                                         {"Index": 2, "Attributes": [1], "TerminalAttributes": [5002, 6002], "DataFile": "b.csv"}],
                                 "Postprocessing": {"Dielectric": [
                                     {"Index": 1, "Attributes": [6001, 6002], "Type": "MA", "EdgeAttributes": [3100], "EdgeExcludeAttributes": [1]},
                                     {"Index": 2, "Attributes": [5001, 5002], "Type": "MS", "EdgeAttributes": [3100], "EdgeExcludeAttributes": [1]}]}}}
        expanded, shell_map = case_inputs.expand_radial_shells(config, shells)
        self.assertEqual(expanded["Boundaries"]["Ground"]["Attributes"], [5001, 16001, 26001, 36001, 5002, 16002, 26002, 36002])
        self.assertEqual(expanded["Boundaries"]["PrescribedPotential"][1]["TerminalAttributes"], [5002, 16002, 26002, 36002])
        entries = expanded["Boundaries"]["Postprocessing"]["Dielectric"]
        self.assertEqual([(e["Index"], e["Type"], e["Attributes"]) for e in entries],
                         [(3, "MA", [16001, 16002]), (4, "MA", [26001, 26002]), (5, "MA", [36001, 36002]), (2, "MS", [5001, 5002])])
        self.assertEqual(sorted(shell_map), [3, 4, 5])
        self.assertEqual(shell_map[4]["Ordinal"], 2)
        self.assertEqual(shell_map[3]["BaseIndex"], 1)
        self.assertEqual(case_inputs.interface_types(expanded), {2: "MS", 3: "MA", 4: "MA", 5: "MA"})
        # A mixed entry (a parent next to a foreign attribute) fails closed.
        mixed = copy.deepcopy(config)
        mixed["Boundaries"]["Postprocessing"]["Dielectric"][0]["Attributes"] = [6001, 4001]
        with self.assertRaisesRegex(case_inputs.CaseInputError, "pure MA entry"):
            case_inputs.expand_radial_shells(mixed, shells)
        # The recipe tolerance, and the labeled per-case deviation (fail closed when unlabeled).
        manifest = {"ProductionRecipe": {"PhysicsRun": {"Order": 4, "LinearTol": 1e-10, "Provenance": "p"}},
                    "Calibration": {"GateDeviations": {"LinearTol": {"Production": 1e-10, "Calibration": 1e-8, "Cases": ["c"],
                                                                     "ProductionUse": "FORBIDDEN"}}}}
        case = {"Id": "c", "Calibration": {"PhysicsRun": {"LinearTol": 1e-8}}}
        self.assertEqual(case_inputs.physics_run_parameters(manifest)["LinearTol"], 1e-10)
        parameters = case_inputs.physics_run_parameters(manifest, case)
        self.assertEqual((parameters["Order"], parameters["LinearTol"]), (4, 1e-8))
        self.assertEqual(parameters["Deviation"]["LinearTol"]["Production"], 1e-10)
        self.assertEqual(case_inputs.physics_run_parameters(manifest, {"Id": "other", "Calibration": {}})["LinearTol"], 1e-10)
        for broken in ({"Id": "other", "Calibration": {"PhysicsRun": {"LinearTol": 1e-8}}},
                       {"Id": "c", "Calibration": {"PhysicsRun": {"LinearTol": 1e-9}}},
                       {"Id": "c", "Calibration": {"PhysicsRun": {"LinearTol": 1e-8, "Order": 5}}}):
            with self.assertRaisesRegex(case_inputs.CaseInputError, "labeled calibration-only deviation"):
                case_inputs.physics_run_parameters(manifest, broken)
        stale = copy.deepcopy(manifest)
        stale["Calibration"]["GateDeviations"]["LinearTol"]["Production"] = 1e-9
        with self.assertRaisesRegex(case_inputs.CaseInputError, "labeled calibration-only deviation"):
            case_inputs.physics_run_parameters(stale, case)

    def test_compare_and_p_sequence_sum_shells_per_type(self):
        """A run whose MA is split over shell interfaces 3..5 against a reference with one MA
        interface 1: p_MA of the run sums the shells, the per-entry rows cover the common
        MS interface only, and p_sequence reads the reference with its own map."""
        header = ["interface", "edge", "R (m)", "basis_i", "basis_j", "Q_ij (J)", "Q_ij normal (J)", "Q_ij tangential (J)",
                  "Q_total_ij (J)", "Q_total_ij normal (J)", "Q_total_ij tangential (J)"]
        sources = [1, 2]
        domain = {(1, 1): 2.0, (2, 2): 4.0, (1, 2): 0.5}
        ma = {(1, 1): 0.2, (2, 2): 0.8, (1, 2): 0.1}
        ms = {(1, 1): 0.3, (2, 2): 0.6, (1, 2): 0.05}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name, split in (("ref", None), ("run", (0.5, 0.3, 0.2))):
                directory = root / name
                directory.mkdir()
                with (directory / "domain-response-matrix.csv").open("w", newline="") as stream:
                    writer = csv.writer(stream)
                    writer.writerow(["basis_i", "basis_j", "Q_ij (J)"])
                    for (i, j), value in domain.items():
                        writer.writerow([i, j, value * (1.0 if split is None else 1.001)])
                with (directory / "surface-response-matrix.csv").open("w", newline="") as stream:
                    writer = csv.writer(stream)
                    writer.writerow(header)
                    for (i, j), value in ma.items():
                        if split is None:
                            writer.writerow([1, 1, 2e-6, i, j] + [value] * 6)
                        else:
                            for index, fraction in zip((3, 4, 5), split):
                                writer.writerow([index, 1, 2e-6, i, j] + [value * fraction * 1.01] * 6)
                    for (i, j), value in ms.items():
                        writer.writerow([2, 1, 2e-6, i, j] + [value] * 6)
            run_types = {2: "MS", 3: "MA", 4: "MA", 5: "MA"}
            reference_types = {1: "MA", 2: "MS"}
            with self.assertRaisesRegex(ValueError, "reference surface matrix interfaces"):
                compare_matrices.compare(root / "ref", root / "run", interface_types=run_types)
            comparison = compare_matrices.compare(root / "ref", root / "run", interface_types=run_types,
                                                  reference_interface_types=reference_types)
            per_source = comparison["PerSource"]
            self.assertAlmostEqual(per_source[1]["p_MA_rel"], 1.01 / 1.001 - 1.0, places=12)
            self.assertAlmostEqual(per_source[2]["p_MS_rel"], 1.0 / 1.001 - 1.0, places=12)
            self.assertEqual(comparison["InterfaceNames"], {2: "MS"})
            self.assertEqual({row["interface"] for row in comparison["Rows"] if row["kind"] == "surface"}, {2})
            summary = compare_matrices.write_comparison(root / "ref", root / "run", root / "cmp", interface_types=run_types,
                                                        reference_interface_types=reference_types)
            self.assertEqual(sorted(summary["PerSource"]), ["1", "2"])
            sequence = p_sequence.p_sequence({"main": root / "run", "ref": root / "ref", "low": None, "high": None}, sources,
                                             run_types, reference_types)
            self.assertAlmostEqual(sequence[2]["p_MA"]["vs_ref"]["main"], 1.01 / 1.001 - 1.0, places=12)
            self.assertAlmostEqual(sequence[1]["Q_MA"]["values"]["main"], 0.2 * 1.01, places=12)
            self.assertAlmostEqual(sequence[1]["Q_MA"]["values"]["ref"], 0.2, places=12)

    def test_radial_profile_fit_recovers_a_power_law(self):
        """Ring energies of f(r) = c r^alpha over the production ring set: the free fit
        returns alpha and a zero remainder; a deficient innermost ring gives the remainder
        exactly; the theoretical -2/3 fit on a -2/3 law is exact."""
        bounds = [0.0] + RADII
        for alpha in (-2.0 / 3.0, -0.4):
            shells = {k: (bounds[k - 1], bounds[k], 3.0 * radial_ma_profile.ring_energy_factor(alpha, bounds[k - 1], bounds[k]))
                      for k in range(1, 8)}
            profile = radial_ma_profile.profile_of(shells)
            for name in ("Fit2-K", "Fit2-4"):
                estimate = profile["Estimates"][name]
                self.assertAlmostEqual(estimate["Alpha"], alpha, places=4)
                self.assertLess(estimate["ResidualRMS"], 1e-6)
                self.assertAlmostEqual(estimate["Remainder"] / shells[1][2], 0.0, places=4)
            self.assertTrue(profile["CleanPowerLaw"])
            if alpha == -2.0 / 3.0:
                for name in ("Theory@2", "Theory2-K"):
                    self.assertAlmostEqual(profile["Estimates"][name]["Remainder"] / shells[1][2], 0.0, places=8)
            deficient = dict(shells)
            deficient[1] = (shells[1][0], shells[1][1], 0.8 * shells[1][2])
            profile = radial_ma_profile.profile_of(deficient)
            self.assertAlmostEqual(profile["Estimates"]["Fit2-K"]["Remainder"] / shells[1][2], 0.2, places=4)
            self.assertAlmostEqual(profile["Estimates"]["Fit2-4"]["Ring1ResolvedOverModel"], 0.8, places=4)
        # Theory@2 anchors the -2/3 law on ring 2 alone: Q_1 model = Q_2 r_1^(1/3) / (r_2^(1/3) - r_1^(1/3)).
        estimate = profile["Estimates"]["Theory@2"]
        self.assertAlmostEqual(estimate["Ring1Model"], shells[2][2] * 0.25 ** (1 / 3) / (0.75 ** (1 / 3) - 0.25 ** (1 / 3)), places=12)
        # Local slopes of a pure power law equal alpha up to the midpoint approximation
        # (the innermost ring, whose inner radius is 0, is off by more: not a fitted ring).
        slopes = radial_ma_profile.local_slopes([shells[k] for k in range(1, 8)])
        self.assertTrue(all(abs(s + 0.4) < 0.02 for s in slopes[1:]))
        self.assertLess(slopes[0], -0.5)


if __name__ == "__main__":
    unittest.main()
