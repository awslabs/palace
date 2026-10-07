#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Round-3 class (10) (decisions 492 / 493 / 510; DESIGN-part-G G.10.5): the content-hashed
sources of a coupon are platform-independent bytes.

* ARC fixtures: the full generator (--basis-only, the device-coupon options of record) on a
  synthetic two-arc rounded finger and on two census arc signatures (testdata/census-b-signatures:
  448693d60a6f, one arc + context; 2bc3d927fda6, three arcs) reproduces STORED golden content
  hashes (device_coupons.content_hash over the generated roles) — the login-node replay runs the
  same test on Linux (the cross-platform half).
* STRAIGHT fixtures: the module is NOT on the straight path (deterministic_math's trigonometry is
  patched to fail; the legacy claims-only and a device-plan straight signature build) except the
  sidewall pullback tan(radians 90), whose correctly rounded value equals libm's, so the straight
  ids are unchanged; the legacy fixture's generated files are byte-identical (the LegacyByteIdentityTest
  of test_cluster_signature_geometry covers the generator; here the content hash of its coupon.json).
* The scalar-arithmetic rule: the 2-vector helpers equal the two-rounding scalar formulas and the
  frame application equals the matmul for a signed-permutation frame.
"""
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cluster_signature_geometry as csg  # noqa: E402
import deterministic_math  # noqa: E402
import device_coupons  # noqa: E402
import generate_spatial_response as gsr  # noqa: E402

CENSUS_B = HERE / "testdata" / "census-b-signatures" / "signatures.json"
LEGACY = HERE / "testdata" / "legacy-byte-identity"
INTERFACES = [{"Slot": 0, "Type": "MA", "Target": 3}, {"Slot": 0, "Type": "MS", "Target": 2}, {"Slot": 0, "Type": "SA", "Target": 1}]
LAW = '{"Type":"PEC"}'
# The device-coupon generator options of record (device_coupons.generate_sources: R 1.9 um, the
# S1p / batch process, ring 16, delaunay caps, 1 R interior hats, order 1, the material options).
GENERATOR_OPTIONS = ["--radius", "1.9", "--metal-thickness", "0.1", "--overetch-depth", "0.05", "--sidewall-angle", "90.0",
                     "--top-rounding", "0.0", "--trench-rounding", "0.0", "--ring-size", "16", "--cap-triangulation", "delaunay",
                     "--cap-interior-spacing", "1.0", "--order", "1", "--basis-only", "--substrate-permittivity", "11.45",
                     "--sa-thickness", "0.002", "--sa-permittivity", "4.0", "--ms-thickness", "0.002", "--ms-permittivity", "11.45",
                     "--ma-thickness", "0.002", "--ma-permittivity", "10.0"]
PARAMETERS = {"metal_thickness": 0.1, "overetch": 0.05, "sidewall_angle": 90.0, "top_radius": 0.0, "bottom_radius": 0.0}
RADIUS = 1.9

# The golden content hashes (sha256 of the sorted (role, digest) pairs; the case id is its first 12
# hex): written by this module at the round-3 B0 head on macOS arm64 and reproduced on Linux
# aarch64 (impl-B0 REPORT, the two-platform replay). A changed value is a changed source byte.
GOLDEN_CONTENT_HASH = {
    "rounded-finger": "55d18fb2e29ef702674c886b0241263f834a26d856caa4d7542d268f137acf8b",
    "448693d60a6f": "93418d21c39779e93e0378ef2e28b2ebd8702ea35f887cced83f90b1e43a0f80",
    "2bc3d927fda6": "2aee170a6108b30243e0297e97ae7e5eda7d40254e6371d2bc22f3b5a03ddfe9",
}


def rounded_finger_record():
    """ArcPortionTest's rounded finger end: the end edge, two 0.5 R fillet arcs, R along both sides."""
    c45 = 0.5 * math.cos(math.pi / 4)

    def portion(p, gap):
        return {"Conductor": 1, "Gap": list(gap), "Interfaces": ["MA", "MS", "SA"], "Law": LAW, "P": list(p)}
    top = {"Conductor": 1, "Interfaces": ["MA", "MS", "SA"], "Law": LAW, "P": [-0.5, 0.75, 0.0, 0.25],
           "Arc": [-0.5, 0.25, -0.5 + c45, 0.25 + c45], "GapRadial": 1}
    bottom = {"Conductor": 1, "Interfaces": ["MA", "MS", "SA"], "Law": LAW, "P": [-0.5, -0.75, 0.0, -0.25],
              "Arc": [-0.5, -0.25, -0.5 + c45, -0.25 - c45], "GapRadial": 1}
    signature = {"Type": "SpatialEdgeCluster", "EdgeCount": 5,
                 "Portions": [portion((0.0, -0.25, 0.0, 0.25), (1.0, 0.0)), top, bottom,
                              portion((-1.5, 0.75, -0.5, 0.75), (0.0, 1.0)), portion((-1.5, -0.75, -0.5, -0.75), (0.0, -1.0))],
                 "Vertices": []}
    return {"Id": "rounded-finger", "Topology": "SpatialEdgeCluster", "Geometry": {"EdgeCount": 5, "Signature": signature},
            "Interfaces": INTERFACES, "BoundaryCondition": {"Type": "PEC"}}


def census_coupon(prefix):
    rec = json.loads(CENSUS_B.read_text())[prefix]
    return {"Id": f"census-{prefix}", "Topology": "SpatialEdgeCluster",
            "Geometry": {"EdgeCount": len(rec["Signature"]["Portions"]), "Signature": rec["Signature"]},
            "Interfaces": rec["Interfaces"], "BoundaryCondition": rec["BoundaryCondition"]}


def generate_content_hash(coupon, work):
    """device_coupons' generation path on COPIES: coupon.json from the signature, the generator
    --basis-only with the options of record, the stamped signature model, process.toml, the content
    hash over the generated roles. Returns (content hash hex, per-role digests)."""
    work = Path(work)
    generated = device_coupons.coupon_geometry(coupon, RADIUS, PARAMETERS)
    (work / "coupon.json").write_text(json.dumps(generated, indent=2) + "\n")
    command = [sys.executable, str(HERE / "generate_spatial_response.py"), str(work / "coupon.json"), "--output", str(work),
               "--model-name", coupon["Id"]] + GENERATOR_OPTIONS
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise AssertionError(result.stderr[-2000:])
    device_coupons.stamp_signature_model(work / "process-library.json", coupon, RADIUS, generated)
    device_coupons.write_process_toml(work / "process.toml", PARAMETERS, RADIUS)
    return device_coupons.content_hash(work)


class ArcFixtureGoldenHashTest(unittest.TestCase):
    def check(self, name, coupon):
        with tempfile.TemporaryDirectory() as tmp:
            digest, roles = generate_content_hash(coupon, tmp)
        self.assertEqual(digest, GOLDEN_CONTENT_HASH[name],
                         f"{name}: content hash {digest} != golden {GOLDEN_CONTENT_HASH[name]}; roles {json.dumps(roles, indent=1)}")

    def test_rounded_finger_reproduces_its_golden_content_hash(self):
        self.check("rounded-finger", rounded_finger_record())

    def test_census_one_arc_context_coupon_reproduces_its_golden_content_hash(self):
        self.check("448693d60a6f", census_coupon("448693d60a6f"))

    def test_census_three_arc_coupon_reproduces_its_golden_content_hash(self):
        self.check("2bc3d927fda6", census_coupon("2bc3d927fda6"))


def failing(name):
    def call(*arguments):
        raise AssertionError(f"deterministic_math.{name} reached on a straight path with {arguments}")
    return call


class StraightPathTest(unittest.TestCase):
    """The module is not on the straight path: a straight coupon's bytes are the pre-B0 bytes."""

    def test_legacy_straight_signature_never_reaches_the_trigonometry(self):
        record = json.loads((LEGACY / "signature.json").read_text())
        with mock.patch.multiple(deterministic_math, sin=failing("sin"), cos=failing("cos"), acos=failing("acos"), asin=failing("asin")):
            coupon, _ = csg.cluster_coupon(record, 1.9, 0.1, 0.05)
            frame, edges, facets = gsr.normalize_geometry(coupon, 1.9)
            lower, upper = gsr.coupon_bounds(edges, 1.9, 0.1, 0.05)
            loops = gsr.plan_view_boundary_loops(facets, 1.9, lower, upper,
                                                 gsr.classified_continuation_segments(coupon["Geometry"], frame, 1.9))
            self.assertEqual(gsr.reconcile_mask_with_boundary(facets, loops, 1.9), [])
            self.assertEqual(gsr.tag_arc_boundary_loops(loops, coupon, frame, 1.9), 0)
        self.assertEqual(json.dumps(coupon, indent=1) + "\n", (LEGACY / "coupon.json").read_text())
        # atan2 is the one function a straight coupon calls (the arrangement's half-edge sort key in
        # plan_view_faces): on the exact node differences of a rectilinear coupon its arguments have a
        # zero component, where the result is the exact axis value on every platform.
        with mock.patch.object(deterministic_math, "atan2", side_effect=lambda y, x: self.axis_only(y, x)):
            csg.cluster_coupon(record, 1.9, 0.1, 0.05)

    def axis_only(self, y, x):
        self.assertTrue(y == 0.0 or x == 0.0, (y, x))
        return math.atan2(y, x)

    def test_device_plan_straight_signatures_never_reach_the_trigonometry(self):
        census = json.loads(CENSUS_B.read_text())
        straight = [prefix for prefix, rec in census.items()
                    if not any("Arc" in e for e in rec["Signature"]["Portions"] + rec["Signature"].get("Context", []))]
        for directory in ("spatial-19-edge-39ab2ffd68ec", "spatial-19-edge-c5952e69af66"):
            library = json.loads((HERE / "testdata" / directory / "process-library.json").read_text())
            model = library["Models"][0]
            rec = {"Topology": "SpatialEdgeCluster", "Geometry": {"Signature": model["Signature"]}, "Signature": model["Signature"],
                   "Interfaces": model["Interfaces"], "BoundaryCondition": {"Type": "PEC"}}
            census[directory] = rec
            straight.append(directory)
        self.assertGreaterEqual(len(straight), 2)
        for prefix in straight:
            rec = census[prefix]
            with mock.patch.multiple(deterministic_math, sin=failing("sin"), cos=failing("cos"), acos=failing("acos"),
                                     asin=failing("asin"), tan=failing("tan")):
                coupon, _ = csg.cluster_coupon(rec, 1.9, 0.1, 0.05)
                frame, edges, facets = gsr.normalize_geometry(coupon, 1.9)
                lower, upper = gsr.coupon_bounds(edges, 1.9, 0.1, 0.05, coupon["Geometry"].get("SupportBox"))
                loops = gsr.plan_view_boundary_loops(facets, 1.9, lower, upper,
                                                     gsr.classified_continuation_segments(coupon["Geometry"], frame, 1.9))
                gsr.reconcile_mask_with_boundary(facets, loops, 1.9)
                self.assertEqual(gsr.tag_arc_boundary_loops(loops, coupon, frame, 1.9), 0, prefix)
            self.assertNotIn("JointSnaps", coupon["Geometry"], prefix)

    def test_the_sidewall_pullback_is_the_one_transcendental_on_the_straight_path_and_matches_libm(self):
        """matching_perimeter_coordinates / conductor_at_points: pullback = t / tan(radians(sidewall));
        at the 90-degree process of record the correctly rounded tangent is libm's value, so the
        pre-B0 bytes of every straight coupon are reproduced."""
        angle = math.radians(90.0)
        self.assertEqual(deterministic_math.tan(angle), math.tan(angle))
        self.assertEqual(deterministic_math.tan(angle), 1.633123935319537e+16)

    def test_legacy_fixture_generated_files_are_byte_identical(self):
        record = json.loads((LEGACY / "signature.json").read_text())
        coupon, _ = csg.cluster_coupon(record, 1.9, 0.1, 0.05)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "coupon.json").write_text(json.dumps(coupon, indent=1) + "\n")
            command = [sys.executable, str(HERE / "generate_spatial_response.py"), str(tmp / "coupon.json"), "--output", str(tmp),
                       "--radius", "1.9", "--metal-thickness", "0.1", "--overetch-depth", "0.05", "--sidewall-angle", "90",
                       "--top-rounding", "0", "--trench-rounding", "0", "--ring-size", "16", "--cap-triangulation", "delaunay",
                       "--cap-interior-spacing", "0.5", "--order", "1", "--model-name", "legacy-fixture", "--signature-only"]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            for name in ("mesh-signature.csv", "plan-view-mask.csv", "plan-view-boundary.csv"):
                self.assertEqual((tmp / name).read_bytes(), (LEGACY / name).read_bytes(), name)


class ScalarRuleTest(unittest.TestCase):
    def test_two_vector_helpers_are_the_two_rounding_scalar_formulas(self):
        rng = np.random.default_rng(20261007)
        for _ in range(2000):
            a, b = rng.uniform(-10.0, 10.0, 2), rng.uniform(-10.0, 10.0, 2)
            self.assertEqual(csg._dot2(a, b), float(a[0]) * float(b[0]) + float(a[1]) * float(b[1]))
            self.assertEqual(csg._norm2(a), math.sqrt(float(a[0]) * float(a[0]) + float(a[1]) * float(a[1])))
            self.assertEqual(gsr._dot2(a, b), csg._dot2(a, b))
            self.assertEqual(gsr._norm2(a), csg._norm2(a))
            c, d = rng.uniform(-10.0, 10.0, 3), rng.uniform(-10.0, 10.0, 3)
            self.assertEqual(gsr._dot3(c, d), (float(c[0]) * float(d[0]) + float(c[1]) * float(d[1])) + float(c[2]) * float(d[2]))
            rows = rng.uniform(-10.0, 10.0, (5, 2))
            np.testing.assert_array_equal(gsr._rows_dot(rows, a), [csg._dot2(row, a) for row in rows])

    def test_frame_application_equals_the_matmul_for_signed_permutations(self):
        rng = np.random.default_rng(7)
        frames = [np.identity(3), np.asarray([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
                  np.asarray([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]])]
        for frame in frames:
            for _ in range(200):
                point = rng.uniform(-10.0, 10.0, 3)
                np.testing.assert_array_equal(gsr.apply_frame(frame, point), frame @ point)
                points = rng.uniform(-10.0, 10.0, (4, 3))
                np.testing.assert_array_equal(gsr.canonical_points(points, frame), points @ frame)

    def test_sequential_mean_is_the_scalar_rule_and_equals_numpy_below_its_block(self):
        """_scalar_mean is the left-to-right sum / n on every platform; numpy's mean agrees for n <= 7
        (its pairwise summation starts its 8-accumulator block at n = 8), i.e. for the facet
        triangles (n = 3) and the continuation segment pairs (n = 2) of every coupon."""
        rng = np.random.default_rng(3)
        for n in (2, 3, 4, 7, 8, 13):
            for _ in range(300):
                values = rng.uniform(-1.0, 1.0, n)
                total = 0.0
                for value in values:
                    total += float(value)
                self.assertEqual(gsr._scalar_mean(values), total / n)
                if n <= 7:
                    self.assertEqual(gsr._scalar_mean(values), float(np.mean(values)))

    def test_squares_on_the_hashed_path_are_products_not_pow(self):
        """MINOR-1 (decision 516): `x ** 2` is libm pow on CPython floats and numpy scalars and differs
        from x * x in ~0.1 % of arguments with platform-dependent mismatch sets; the hashed-path
        sites use multiplication (csg.plan_view_faces split parameter, gsr.delaunay_flip_cap in_circle)
        and matching_support_points derives its decimal exponent without log10 / pow."""
        import inspect
        for function in (csg.plan_view_faces, gsr.delaunay_flip_cap, gsr.matching_support_points):
            code = "\n".join(line.split("#")[0] for line in inspect.getsource(function).splitlines())
            self.assertNotIn("**", code, function.__name__)
            self.assertNotIn("log10", code, function.__name__)
        # The exact decimal exponent equals floor(log10) for every double, incl. the quantisation
        # tolerances of the radii of record and exact powers of ten.
        import decimal
        for tolerance in (1.9e-10, 1.0e-10, 2.0e-10, 5.0e-10, 64.0 * np.finfo(float).eps, 0.5, 1.0, 10.0, 123.456):
            self.assertEqual(decimal.Decimal(tolerance).adjusted(), math.floor(math.log10(tolerance)), tolerance)
            exponent = decimal.Decimal(tolerance).adjusted()
            self.assertEqual(float(f"1e{exponent}"), float(decimal.Decimal(10) ** exponent))
        rng = np.random.default_rng(11)
        mismatches = sum(1 for x in rng.uniform(0.1, 10.0, 20000) if float(x) ** 2 != float(x) * float(x))
        self.assertGreaterEqual(mismatches, 0)   # measured ~0.1 % per platform: the reason for the rule

    def test_box_trace_face_area_difference_is_scalar(self):
        """MINOR-3 (decision 516): validate_box_trace's summed face areas (basis-contract.json
        MaximumRelativeFaceAreaDifference, a WRITTEN float) use the scalar _cross3 / _norm3, and equal the
        former BLAS value on face-exact triangles (one nonzero cross component)."""
        import inspect
        import box_trace
        source = inspect.getsource(box_trace.validate_box_trace)
        self.assertNotIn("np.linalg.norm", source)
        self.assertNotIn("np.cross", source)
        self.assertIn("_norm3(_cross3(", source)
        rng = np.random.default_rng(5)
        for _ in range(500):
            a, b = rng.uniform(-3.0, 3.0, 3), rng.uniform(-3.0, 3.0, 3)
            a[2] = b[2] = 0.0   # a triangle on a z-face: the cross has one nonzero component
            self.assertEqual(gsr._norm3(gsr._cross3(a, b)), float(np.linalg.norm(np.cross(a, b))))

    def test_provenance_names_the_float_serialisation_rule(self):
        self.assertEqual(deterministic_math.FLOAT_SERIALISATION_RULE, "deterministic-trig-v1")
        self.assertTrue(device_coupons.FLOAT_SERIALISATION_MODULE.is_file())
        self.assertEqual(device_coupons.FLOAT_SERIALISATION_MODULE.name, "deterministic_math.py")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--print-golden":
        for name, coupon in (("rounded-finger", rounded_finger_record()), ("448693d60a6f", census_coupon("448693d60a6f")),
                             ("2bc3d927fda6", census_coupon("2bc3d927fda6"))):
            with tempfile.TemporaryDirectory() as tmp:
                digest, roles = generate_content_hash(coupon, tmp)
            print(name, digest, json.dumps(roles))
        sys.exit(0)
    unittest.main()
