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
import math
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
        self.assertLess(closure["RelativeClosure"], 1e-12)
        self.assertAlmostEqual(closure["OwnedMeasure"], closure["WholeMeasure"], places=14)
        # The closure is independent of the census: the quadrature weights are closed against
        # the cross-product areas (element_area), so a rule that lost weight fails it.
        self.assertAlmostEqual(relabel.element_area(2, [(0, 0, 0), (1, 0, 0), (0, 2, 0)]), 1.0, places=15)
        self.assertAlmostEqual(relabel.element_area(3, [(0, 0, 0), (2, 0, 0), (2, 1, 0), (0, 1, 0)]), 2.0, places=15)
        lossy = copy.deepcopy(shells)
        lossy[(1 + rings + 1, 6001)]["QuadratureMeasure"] *= 0.5
        with self.assertRaisesRegex(relabel.RelabelError, "does not close"):
            relabel.closure_check(lossy)

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

    def test_radial_profile_analyze_and_markdown(self):
        """analyze() on a synthetic shell run (one main order, 15 MA shell interfaces next to
        MS 2, a reference with one MA interface): the Consistent estimator sums the top edge's
        Theory@2 remainder and the bottom edge's Fit2-4 remainder; a source whose bottom rings
        carry no energy (no bottom fit) and a second order are handled (None guards) and
        markdown() renders both."""
        header = ["interface", "edge", "R (m)", "basis_i", "basis_j", "Q_ij (J)", "Q_ij normal (J)", "Q_ij tangential (J)",
                  "Q_total_ij (J)", "Q_total_ij normal (J)", "Q_total_ij tangential (J)"]
        bounds = [0.0] + RADII
        rings = len(RADII)
        top = {k: 3.0 * radial_ma_profile.ring_energy_factor(-2.0 / 3.0, bounds[k - 1], bounds[k]) for k in range(1, rings + 1)}
        top[1] *= 0.7  # the innermost ring resolves 0.7 of the -2/3 law
        bottom = {k: 5.0 * radial_ma_profile.ring_energy_factor(-0.4, bounds[k - 1], bounds[k]) for k in range(1, rings + 1)}
        bottom[1] *= 0.9
        far = 0.05
        shells = {"3": {"Kind": "far", "Ring": 0, "InnerRadius": RADII[-1], "OuterRadius": None}}
        for k in range(1, rings + 1):
            shells[str(3 + k)] = {"Kind": "top", "Ring": k, "InnerRadius": bounds[k - 1], "OuterRadius": bounds[k]}
            shells[str(3 + rings + k)] = {"Kind": "bottom", "Ring": k, "InnerRadius": bounds[k - 1], "OuterRadius": bounds[k]}
        interfaces = {"2": "MS", **{index: "MA" for index in shells}}
        # Source 1: both edges; source 2: no bottom energy (no bottom fit); source 3: zero trace.
        energies = {1: 10.0, 2: 20.0, 3: 5.0}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root / "reference"
            reference.mkdir()
            (root / "reference-config.json").write_text(json.dumps(
                {"Boundaries": {"Postprocessing": {"Dielectric": [{"Index": 1, "Type": "MA"}, {"Index": 2, "Type": "MS"}]}}}))
            for order, scale in (("p4", 1.0), ("p5", 1.01)):
                reducer = root / "case" / "results" / "main" / f"syn-{order}" / "reducer"
                reducer.mkdir(parents=True)
                with (reducer / "domain-response-matrix.csv").open("w", newline="") as stream:
                    writer = csv.writer(stream)
                    writer.writerow(["basis_i", "basis_j", "Q_ij (J)"])
                    for i, e in energies.items():
                        writer.writerow([i, i, e])
                with (reducer / "surface-response-matrix.csv").open("w", newline="") as stream:
                    writer = csv.writer(stream)
                    writer.writerow(header)
                    for i in energies:
                        for index, shell in shells.items():
                            if shell["Kind"] == "far":
                                value = far
                            elif shell["Kind"] == "top":
                                value = top[shell["Ring"]] * (scale if shell["Ring"] == 1 else 1.0)
                            else:
                                value = 0.0 if i == 2 else bottom[shell["Ring"]]
                            writer.writerow([int(index), 1, 2e-6, i, i] + [value] * 6)
                        writer.writerow([2, 1, 2e-6, i, i] + [0.3] * 6)
            for name in ("domain-response-matrix.csv", "surface-response-matrix.csv"):
                (reference / name).write_text((root / "case" / "results" / "main" / "syn-p4" / "reducer" / name).read_text()
                                              .replace("\n3,1,2e-06", "\n1,1,2e-06"))
            record = {"Cases": [{"Case": "syn", "Root": str(root / "case"),
                                 "Inputs": {"RadialShells": {"Interfaces": shells}, "Interfaces": interfaces},
                                 "Sources": {"ZeroTrace": [3]},
                                 "Stages": [{"Role": "main", "Prefix": "syn-p4", "Order": 4}, {"Role": "main", "Prefix": "syn-p5", "Order": 5},
                                            {"Role": "control", "Prefix": "syn-p3-control", "Order": 3}],
                                 "Reference": {"Results": str(reference), "Config": str(root / "reference-config.json")}}]}
            (root / "library-qualification.json").write_text(json.dumps(record))
            out = radial_ma_profile.analyze(root / "library-qualification.json", "syn", strongest=1)
            self.assertEqual(out["HeadlineEstimator"], "Consistent")
            self.assertEqual(out["CombinedEstimators"]["Consistent"], {"top": "Theory@2", "bottom": "Fit2-4"})
            p4 = out["Orders"]["p4"]
            self.assertEqual(p4["Sources"], [1, 2])
            self.assertEqual(p4["Strongest"], [1])
            one = p4["PerSource"]["1"]
            q_ma = sum(top.values()) + sum(bottom.values()) + far
            self.assertAlmostEqual(one["Q_MA"], q_ma, places=12)
            top_theory = one["Kinds"]["top"]["Estimates"]["Theory@2"]["Remainder"]
            bottom_fit = one["Kinds"]["bottom"]["Estimates"]["Fit2-4"]["Remainder"]
            self.assertAlmostEqual(top_theory / (top[1] / 0.7), 0.3, places=6)
            self.assertAlmostEqual(bottom_fit / (bottom[1] / 0.9), 0.1, places=3)
            self.assertAlmostEqual(one["Remainder"]["Consistent"], top_theory + bottom_fit, places=15)
            self.assertAlmostEqual(one["Deficit"]["Consistent"], (top_theory + bottom_fit) / q_ma, places=15)
            self.assertAlmostEqual(one["DeficitTop"]["Consistent"], top_theory / q_ma, places=15)
            # Theory@2 on both edges overstates the bottom remainder against the Consistent estimator.
            self.assertGreater(one["Deficit"]["Theory@2"], one["Deficit"]["Consistent"])
            self.assertLess(one["Deficit"]["Fit2-4"], one["Deficit"]["Consistent"])
            # Source 2 has no bottom energy: no bottom fit, the Consistent remainder is the top part only.
            two = p4["PerSource"]["2"]
            self.assertIsNone(two["Kinds"]["bottom"]["Estimates"]["Fit2-4"])
            self.assertAlmostEqual(two["Remainder"]["Consistent"], two["Kinds"]["top"]["Estimates"]["Theory@2"]["Remainder"], places=15)
            summary = p4["Summary"]
            self.assertEqual(summary["Alpha"]["bottom:Fit2-4"]["Fitted"], 1)
            self.assertAlmostEqual(summary["Alpha"]["bottom:Fit2-4"]["Median"], -0.4, places=3)
            above = summary["Alpha"]["bottom:Fit2-4"]["ShareAboveFloor"]
            self.assertEqual((above["Floor"], above["Sources"]), (radial_ma_profile.BOTTOM_SHARE_FLOOR, 1))
            self.assertNotIn("ShareAboveFloor", summary["Alpha"]["bottom:Theory@2"])
            self.assertAlmostEqual(summary["Deficit"]["Consistent"]["StrongestMedian"], one["Deficit"]["Consistent"], places=15)
            self.assertTrue(all(v is not None for v in summary["LocalSlopeMedians"]["bottom"]))
            self.assertTrue(all(v is None for v in p4["PerSource"]["2"]["Kinds"]["bottom"]["LocalSlopes"]))
            self.assertEqual(summary["Deficit"]["Consistent"]["At"], {})
            # The p-step: ring 1 of the top edge moved by the scale, the other rings did not.
            step = out["PStep"]["MedianRelativeStep"]
            self.assertEqual(out["PStep"]["Orders"], ["p4", "p5"])
            self.assertAlmostEqual(step["top"][0], 0.01, places=12)
            self.assertAlmostEqual(step["top"][1], 0.0, places=12)
            self.assertEqual(step["bottom"][0], 0.0)
            self.assertAlmostEqual(step["Ring1ShareOfStep"]["top"], 1.0, places=12)
            text = radial_ma_profile.markdown(out, detail=[1, 2, 9])
            self.assertIn("| **Consistent** (top Theory@2 + bottom Fit2-4) |", text)
            self.assertIn("bottom-edge alpha Fit2-4 over the 1 sources with bottom share > 20%", text)
            self.assertIn("## p5: 2 free sources", text)
            self.assertIn("## p-step p4 -> p5", text)
            self.assertIn("| 2 |", text)
            self.assertNotIn("\n| 9 |", text)
            # A bottom-less source set (every bottom fit None) still summarizes and renders.
            for order in ("p4", "p5"):
                out["Orders"][order]["PerSource"].pop("1")
                out["Orders"][order]["Sources"] = [2]
                out["Orders"][order]["Strongest"] = [2]
            self.assertIn("n/a", radial_ma_profile.markdown(out, detail=[2]))


if __name__ == "__main__":
    unittest.main()


class ProductionRadialShellsTest(unittest.TestCase):
    """Decision 61a: the shells at the placement stage (publish_rigid_coupon_mesh.
    apply_radial_shells) on a rotated placement of the synthetic strip, the parent view
    the audits read and the ownership report's per-parent sums of the shell rows."""

    def canonical_build(self, root, data, radii=RADII):
        canonical = root / "canonical.msh"
        canonical.write_bytes(data)
        census = root / "build-census.json"
        census.write_text(json.dumps({"PrismTubes": {"Section": {"RingRadii": radii}}}))
        partition = root / "canonical.msh.interface-partition.csv"
        with partition.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["attribute", "elements", "area"])
            writer.writerow([3100, 1, 0.2])
            writer.writerow([6001, 24, 3 * RADII[-1] + 0.2 + 0.5e-8 + 1e-7])
        record = {"CanonicalArtifacts": {
            "canonical-candidate-mesh": {"Path": str(canonical), "SHA256": relabel.sha256(canonical)},
            "build-census": {"Path": str(census), "SHA256": relabel.sha256(census)},
            "canonical-ownership-partition": {"Path": str(partition), "SHA256": relabel.sha256(partition)}}}
        return canonical, record

    def sources(self, root):
        boundary = root / "plan-view-boundary.csv"
        boundary.write_text("Loop,Vertex,Conductor,Plane,Hole,Class,X,Y\n"
                            "0,0,1,0.0,0,Physical,0.0,0.0\n0,1,1,0.0,0,Continuation,0.0,1.0\n"
                            "0,2,1,0.0,0,Continuation,-1.0,1.0\n0,3,1,0.0,0,Continuation,-1.0,0.0\n")
        signature = root / "mesh-signature.csv"
        signature.write_text("Index,Slot,Conductor,Px,Py,Pz,Gx,Gy,Gz,Tx,Ty,Tz,S0,S1,Nz\n"
                             "0,0,1,0.0,0.5,0.0,1.0,0.0,0.0,0.0,1.0,0.0,-0.5,0.5,1.0\n")
        process = root / "process.toml"
        process.write_text('Units = "um"\nRadius = 2.0\nMetalThickness = 0.1\nOveretch = 0.05\n')
        return boundary, signature, process

    def test_rotated_placement_shells_equal_the_identity_shells(self):
        import numpy as np
        import publish_rigid_coupon_mesh as publisher
        data, expected = strip_mesh()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            canonical, build = self.canonical_build(root, data)
            boundary, signature, process = self.sources(root)
            angle = 0.63
            matrix = np.array([[math.cos(angle), -math.sin(angle), 0, 1.2], [math.sin(angle), math.cos(angle), 0, -0.7],
                               [0, 0, 1, 0.9], [0, 0, 0, 1]])
            results = {}
            for name, transform in (("identity", np.eye(4)), ("rotated", matrix)):
                output = root / f"{name}.msh"
                publisher.transform_gmsh22(canonical, output, transform)
                parent_digest = relabel.sha256(output)
                record = publisher.apply_radial_shells(
                    output, Path(str(output) + publisher.RADIAL_SHELLS_SUFFIX), canonical_build=build, matrix=transform,
                    signature=signature, boundary=boundary, process=process, parent_digest=parent_digest,
                    canonical_digest=relabel.sha256(canonical))
                self.assertTrue(record["Applied"])
                self.assertEqual(record["RingRadii"], RADII)
                self.assertEqual(record["MAParents"], [6001])
                self.assertEqual(record["LabelOnly"]["RelabeledElements"], len(expected))
                self.assertLess(record["ParentAreaMaximumRelativeDifference"], 1e-12)
                census = json.loads(Path(record["Path"]).read_text())
                self.assertEqual(census["Mesh"]["SHA256"], relabel.sha256(output))
                self.assertEqual(census["ParentMesh"]["SHA256"], parent_digest)
                self.assertEqual(census["Placement"]["Transform"], [float(v) for row in transform for v in row])
                after = relabel.read_msh22_binary(output.read_bytes())
                results[name] = ([element[2][0] for element in after["Elements"]], census)
            # The shells are classified in source-local coordinates: the rotated placement
            # carries exactly the identity's labels, element by element, and the same areas.
            self.assertEqual(results["identity"][0], results["rotated"][0])
            for a, b in zip(results["identity"][1]["Shells"], results["rotated"][1]["Shells"]):
                self.assertEqual((a["Label"], a["Elements"]), (b["Label"], b["Elements"]))
                self.assertAlmostEqual(a["Area"], b["Area"], places=14)
            # The parent view the audits read collapses the shells to 6001 and nothing else.
            import meshio
            from mixed_mesh import parent_label_view, shell_labels, shell_parent
            shelled = meshio.read(root / "rotated.msh")
            labels = shell_labels(shelled)
            self.assertEqual(labels[0], 16001)
            self.assertEqual(sorted({shell_parent(v) for v in labels}), [6001])
            view = parent_label_view(shelled)
            surface = np.concatenate([np.asarray(v) for cell, v in zip(view.cells, view.cell_data["gmsh:physical"])
                                      if cell.type in ("triangle", "quad")])
            self.assertEqual(sorted(set(surface.tolist())), [1, 3100, 6001] if 1 in surface else [3100, 6001])
            self.assertIs(parent_label_view(view), view)
            # The census is refused where it already exists (fresh outputs only).
            with self.assertRaisesRegex(ValueError, "fresh"):
                publisher.apply_radial_shells(
                    root / "identity.msh", Path(str(root / "identity.msh") + publisher.RADIAL_SHELLS_SUFFIX),
                    canonical_build=build, matrix=np.eye(4), signature=signature, boundary=boundary, process=process,
                    parent_digest="x", canonical_digest="y")
            # A canonical build without a gmsh-build census (legacy pipeline) publishes no shell.
            legacy = {"CanonicalArtifacts": {k: v for k, v in build["CanonicalArtifacts"].items() if k != "build-census"}}
            output = root / "legacy.msh"
            publisher.transform_gmsh22(canonical, output, np.eye(4))
            record = publisher.apply_radial_shells(output, root / "legacy.json", canonical_build=legacy, matrix=np.eye(4),
                                                   signature=signature, boundary=boundary, process=process,
                                                   parent_digest="x", canonical_digest="y")
            self.assertFalse(record["Applied"])
            self.assertEqual(shell_labels(meshio.read(output)), [])
            # A gmsh-build census without a ring set fails closed.
            no_rings = root / "no-rings.json"
            no_rings.write_text(json.dumps({"PrismTubes": {}}))
            build_no_rings = json.loads(json.dumps(build))
            build_no_rings["CanonicalArtifacts"]["build-census"] = {"Path": str(no_rings), "SHA256": relabel.sha256(no_rings)}
            with self.assertRaisesRegex(ValueError, "ring set"):
                publisher.apply_radial_shells(output, root / "no-rings-census.json", canonical_build=build_no_rings,
                                              matrix=np.eye(4), signature=signature, boundary=boundary, process=process,
                                              parent_digest="x", canonical_digest="y")

    def test_ownership_report_sums_the_shell_rows_per_parent(self):
        from general_mesh_audit_producer import _ownership_report
        contract = {"CutSurfaceRoles": ["cut"],
                    "BoundaryLabels": [{"Attribute": 1, "Role": "cut"}, {"Attribute": 3100, "Role": "sa"},
                                       {"Attribute": 6001, "Role": "ma-1"}]}
        summary = ["Gauss4", "4", "120", "1.0", "1.0", "0.0", "1e-12", "0", "0", "1"]
        header = ("attribute,elements,area,ambiguous_area,ambiguous_fraction,unresolved_elements,unresolved_area,"
                  "unresolved_fraction,quadrature_rule,quadrature_order,quadrature_points,quadrature_whole_measure,"
                  "quadrature_owned_measure,quadrature_relative_closure,quadrature_closure_tolerance,quadrature_unmatched,"
                  "quadrature_overlaps,quadrature_positive_weights")
        rows = [[3100, 1, 0.2, 0.0, 0.0, 0, 0.0, 0.0], [16001, 1, 0.3, 0.0, 0.0, 0, 0.0, 0.0],
                [26001, 2, 0.1, 0.05, 0.5, 1, 0.05, 0.5], [96001, 3, 0.4, 0.0, 0.0, 0, 0.0, 0.0]]
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            partition = root / "ownership.csv"
            partition.write_text(header + "\n" + "".join(",".join(str(v) for v in row + summary) + "\n" for row in rows))
            quadrature = root / "ownership.quadrature.csv"
            quadrature.write_text("attribute,measure\n3100,0.2\n16001,0.3\n26001,0.1\n96001,0.4\n")
            report = _ownership_report(partition, quadrature, contract)
        self.assertEqual(report["ResponseOwnership"]["OwnerAttributes"], [3100, 6001])
        self.assertAlmostEqual(report["ResponseOwnership"]["OwnerMeasures"][1], 0.8, places=15)
        self.assertTrue(report["ResponseOwnership"]["NoDuplicateOrMissingOwners"])
        self.assertTrue(report["ResponseOwnership"]["Exhaustive"])
        self.assertEqual(report["PhysicalSurfaceCoverage"]["PartitionRows"], 2)
        self.assertEqual(report["PhysicalSurfaceCoverage"]["InterfaceElements"], 7)
        self.assertEqual(report["RadialShells"]["Labels"], [16001, 26001, 96001])
        self.assertEqual(report["RadialShells"]["Parents"], [6001])
        self.assertEqual(report["RadialShells"]["PartitionRows"], 4)
        self.assertEqual(report["WholeElementAmbiguityDiagnostics"]["AmbiguousRows"], 1)
        self.assertEqual(report["WholeElementAmbiguityDiagnostics"]["UnresolvedElements"], 1)


class LabelOnlyRepublicationTest(unittest.TestCase):
    """verify_label_only_republication: the re-published (shelled) identity vs the previous
    production identity - labels only, the census count, the receipt's parent digest."""

    def test_shelled_strip_vs_its_parent_is_labels_only(self):
        import hashlib
        import verify_label_only_republication as verifier
        data, expected = strip_mesh()
        out, mesh, shells, relabeled, ma_labels = relabel.relabel(data, lines=LINES, radii=RADII)
        label_only = relabel.assert_label_only(data, out, mesh, relabeled)
        census = {"Mesh": {"SHA256": hashlib.sha256(out).hexdigest()}, "LabelOnly": label_only,
                  "Shells": [{"Label": shell["Label"], "Parent": shell["Parent"]} for shell in shells.values()]}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "previous.msh").write_bytes(data)
            (root / "identity.msh").write_bytes(out)
            (root / "census.json").write_text(json.dumps(census))
            (root / "receipt.json").write_text(json.dumps({"ParentLabeledMeshSHA256": hashlib.sha256(data).hexdigest()}))
            record = verifier.verify(root / "previous.msh", root / "identity.msh", root / "census.json", root / "receipt.json")
            self.assertEqual(record["Verdict"], "labels-only")
            self.assertEqual(record["RelabeledElements"], len(expected))
            self.assertEqual(record["RelabeledPerParent"], {"6001": len(expected)})
            self.assertTrue(record["Receipt"]["ParentLabeledEqualsPrevious"])
            self.assertEqual(record["PhysicalNamesRemoved"], [6001])
            self.assertIn(26001, record["PhysicalNamesAdded"])
            # A wrong parent digest, a foreign census and a changed node fail closed.
            (root / "bad-receipt.json").write_text(json.dumps({"ParentLabeledMeshSHA256": "0" * 64}))
            with self.assertRaisesRegex(verifier.RepublicationError, "ParentLabeledMeshSHA256"):
                verifier.verify(root / "previous.msh", root / "identity.msh", root / "census.json", root / "bad-receipt.json")
            foreign = dict(census, Mesh={"SHA256": "1" * 64})
            (root / "foreign.json").write_text(json.dumps(foreign))
            with self.assertRaisesRegex(verifier.RepublicationError, "another mesh"):
                verifier.verify(root / "previous.msh", root / "identity.msh", root / "foreign.json")
            moved = bytearray(data)
            moved[mesh["NodesSpan"][0] + 40] ^= 1
            (root / "moved.msh").write_bytes(bytes(moved))
            with self.assertRaisesRegex(verifier.RepublicationError, "Nodes"):
                verifier.verify(root / "moved.msh", root / "identity.msh", root / "census.json")
            # A relabel that is not a shell of the old label (a foreign MA label) is refused.
            wrong = bytearray(out)
            after = relabel.read_msh22_binary(out)
            index = next(i for i, element in enumerate(after["Elements"]) if element[2][0] >= relabel.SHELL_LABEL_STRIDE)
            offset = after["Elements"][index][1]
            struct.pack_into("<i", wrong, offset, 26002)
            (root / "wrong.msh").write_bytes(bytes(wrong))
            wrong_census = dict(census, Mesh={"SHA256": hashlib.sha256(bytes(wrong)).hexdigest()})
            (root / "wrong-census.json").write_text(json.dumps(wrong_census))
            with self.assertRaisesRegex(verifier.RepublicationError, "not a shell of the old label"):
                verifier.verify(root / "previous.msh", root / "wrong.msh", root / "wrong-census.json")
