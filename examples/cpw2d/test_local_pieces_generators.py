#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""U-PIECES (identification fix round DESIGN A.3 (b) / A.10 G-3 / G-6; decisions 1064 (4), 1096 (8), 1097): the Python
local-piece rule == the C++ LocalPiecesOfParallelEdges, (1) on the closed-form cases of lane B's unit test
(test/unit/test-conductorconsistency-r0a.cpp) and (2) on EVERY ParallelEdgeCluster model of the v3-j library of record
through the C++ library reader itself: a copy of v3-j whose PEC models carry Python's Edges[].LocalPiece and one
conductor reference per piece (on metal) is loaded by the Palace binary (--surface-response-preflight on the S6 J1
window of record); the reader REFUSES BY NAME a LocalPiece that is not the piece its GapDirections bound, so a load that
succeeds is the cross-check, and a deliberately wrong piece must be refused. The binary = $PALACE_EXECUTABLE (the
chain's form; skipped when unset or when the records are absent). Plus the generators' default-OFF identity: the
cluster-signature CSV and the pair / cluster library entries are byte-identical without --local-pieces."""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path

CPW2D = Path(__file__).resolve().parent
sys.path.insert(0, str(CPW2D))
sys.path.insert(0, str(CPW2D.parent / "surface_response_identification"))
import local_pieces as lp  # noqa: E402
import generate_edge_cluster_response as CLUSTER  # noqa: E402
import prepare_surface_response_coupons as PLANNER  # noqa: E402

V3J = ("/Users/simlap/bedrock-tests/coupon-accuracy-assessment-20260913/stage2-20261004/v3h/stretchfix/local/library/"
       "s2-r1p9-v3-j-local/process-library-local.json")
S6_CONFIG = ("/Users/simlap/bedrock-tests/coupon-accuracy-assessment-20260913/readiness-20261029/idfix/lane-B/local/"
             "dryrun-f8007abf2a/sct002-S6-thin-preflight-wi/lane/config.json")
PALACE = os.environ.get("PALACE_EXECUTABLE")
LOCK = Path("/tmp/idfix-heavy.lock")


class ClosedFormCases(unittest.TestCase):
    """The cases of lane B's C++ test, by value."""

    def test_cross_sections(self):
        self.assertEqual(lp.local_pieces_of_parallel_edges([0.0, 2.0, 4.0, 6.0], [1, -1, 1, -1]), [1, 2, 2, 3])
        self.assertEqual(lp.local_pieces_of_parallel_edges([0.0, 3.0, 6.0, 9.0], [1, -1, 1, -1]), [1, 2, 2, 3])
        self.assertEqual(lp.local_pieces_of_parallel_edges([0.0, 2.0, 4.0], [1, -1, 1]), [1, 2, 2])
        self.assertEqual(lp.local_pieces_of_parallel_edges([-1.0, 1.0], [1, -1]), [1, 2])
        self.assertEqual(lp.local_pieces_of_parallel_edges([-1.0, 1.0], [-1, 1]), [1, 1])
        self.assertEqual(lp.local_pieces_of_parallel_edges([0.0], [1]), [1])
        self.assertEqual(lp.local_pieces_of_parallel_edges([0.0], [-1]), [1])
        self.assertEqual(lp.local_pieces_of_parallel_edges(list(range(9)), [1, -1, 1, -1, 1, -1, 1, -1, 1]), [1, 2, 2, 3, 3, 4, 4, 5, 5])

    def test_piece_of_offset(self):
        offsets, gaps = [0.0, 2.0, 4.0, 6.0], [1, -1, 1, -1]
        pieces = lp.local_pieces_of_parallel_edges(offsets, gaps)
        for u, piece in ((3.0, 2), (-1.9, 1), (0.0, 1), (2.0, 2), (7.5, 3), (1.0, 0), (5.0, 0)):
            self.assertEqual(lp.local_piece_of_offset(offsets, gaps, pieces, u), piece, u)
        self.assertEqual(lp.local_piece_of_offset([-1.0, 1.0], [1, -1], [1, 2], -1.95), 1)
        self.assertEqual(lp.local_piece_of_offset([-1.0, 1.0], [1, -1], [1, 2], 1.95), 2)

    def test_fails_closed_on_a_non_cross_section(self):
        with self.assertRaisesRegex(ValueError, "no metal between them"):
            lp.local_pieces_of_parallel_edges([0.0, 1.0], [1, 1])
        with self.assertRaisesRegex(ValueError, "no gap between them"):
            lp.local_pieces_of_parallel_edges([0.0, 1.0], [-1, -1])
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            lp.local_pieces_of_parallel_edges([0.0, 0.0], [1, -1])
        with self.assertRaisesRegex(ValueError, "-1 or 1"):
            lp.local_pieces_of_parallel_edges([0.0], [0])


def local_piece_library(models, radius):
    """Every PEC model with Python's Edges[].LocalPiece and one reference per piece at the midpoint of the piece's metal
    interval (on metal); the other models untouched. Returns (models, the PEC count, the piece counts)."""
    out = []
    counts = {}
    for model in models:
        model = json.loads(json.dumps(model))
        if model.get("Topology") == "ParallelEdgeCluster" and model.get("Edges"):
            offsets = [float(e["Offset"]) for e in model["Edges"]]
            gaps = [int(e["GapDirection"]) for e in model["Edges"]]
            pieces = lp.local_pieces_of_parallel_edges(offsets, gaps)
            for edge, piece in zip(model["Edges"], pieces):
                edge["LocalPiece"] = piece
            intervals = lp.metal_intervals_by_piece(offsets, gaps, pieces, radius)
            references = [[0.5 * (lo + hi), 0.0, 0.0] for _, (lo, hi) in sorted(intervals.items())]
            model.pop("Reference", None)
            model.pop("ConductorReferences", None)
            model.pop("OpenContourPaths", None)
            if len(references) == 1:
                model["Reference"] = references[0]
            else:
                # the reader requires open paths connecting every reference: one synthetic path per consecutive piece
                # pair (the preflight reads no matrices; a load-time contract only)
                model["ConductorReferences"] = references
                model["OpenContourPaths"] = [{"Indices": [k + 1], "StartConductor": k + 1, "EndConductor": k + 2}
                                             for k in range(len(references) - 1)]
            counts[model["Name"]] = max(pieces)
        out.append(model)
    return out, counts


def with_heavy_lock(command, timeout):
    """One palace preflight under the lanes' heavy-process lock (decision 1062; mkdir atomic; released after)."""
    import time
    waited = 0
    while True:
        try:
            LOCK.mkdir()
            break
        except FileExistsError:
            time.sleep(5)
            waited += 5
            if waited > 600:
                raise RuntimeError("the heavy lock is held for > 10 min")
    try:
        (LOCK / "owner").write_text(f"C {os.getpid()} unittest {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\n")
        return subprocess.run(command, capture_output=True, text=True, timeout=timeout)
    finally:
        shutil.rmtree(LOCK, ignore_errors=True)


@unittest.skipUnless(PALACE and os.path.isfile(PALACE) and os.path.isfile(V3J) and os.path.isfile(S6_CONFIG),
                     "PALACE_EXECUTABLE / the v3-j mirror / the S6 J1 config of record are needed")
class CxxReaderCrossCheck(unittest.TestCase):
    def preflight(self, tmp, library_models, label):
        with open(V3J) as f:
            library = json.load(f)
        library["Models"] = library_models
        library["Name"] = library["Name"] + "-" + label
        lib_path = tmp / f"library-{label}.json"
        lib_path.write_text(json.dumps(library))
        with open(S6_CONFIG) as f:
            config = json.load(f)
        config["Problem"]["Output"] = str(tmp / label / "postpro")
        config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] = str(lib_path)
        config_path = tmp / f"config-{label}.json"
        config_path.write_text(json.dumps(config))
        self.assertNotIn("lane-B/local/dryrun", config_path.read_text().replace(str(S6_CONFIG), ""))   # re-pointed (838)
        return with_heavy_lock([PALACE, "--surface-response-preflight", str(config_path)], timeout=280)

    def test_every_v3j_pec_model_loads_with_python_pieces_and_a_wrong_piece_is_refused(self):
        with open(V3J) as f:
            library = json.load(f)
        radius = float(library["MatchingRadius"])
        with tempfile.TemporaryDirectory(prefix="upieces-") as tmp:
            tmp = Path(tmp)
            base = self.preflight(tmp, library["Models"], "base")
            self.assertEqual(base.returncode, 0, base.stdout[-2000:] + base.stderr[-2000:])
            models, counts = local_piece_library(library["Models"], radius)
            self.assertEqual(len(counts), 16)                       # the 16 PEC models of v3-j
            self.assertEqual(sorted(set(counts.values())), [2, 3, 4, 5])   # A.3 (c): pieces 2 (3-edge) .. 5 (the 9-edge fb9f050d1e24)
            lane = self.preflight(tmp, models, "pieces")
            self.assertEqual(lane.returncode, 0, lane.stdout[-3000:] + lane.stderr[-3000:])
            # the key is never touched by LocalPiece: the manifests' requirements agree row by row
            def requirements(label):
                with open(tmp / label / "postpro" / "surface-response-requirements.json") as f:
                    m = json.load(f)
                return [(r["Hash"], r["Status"], r.get("TotalEdgeLength")) for r in m["Requirements"]], m["Summary"]["Counts"]
            self.assertEqual(requirements("base"), requirements("pieces"))
            # a deliberately wrong piece on one model: refused by name (the C++'s own check)
            wrong = json.loads(json.dumps(models))
            pec = next(m for m in wrong if m.get("Topology") == "ParallelEdgeCluster")
            pec["Edges"][-1]["LocalPiece"] = pec["Edges"][-1]["LocalPiece"] + 1
            pec["ConductorReferences"].append([pec["ConductorReferences"][-1][0] + 0.1, 0.0, 0.0])
            refused = self.preflight(tmp, wrong, "wrong")
            self.assertNotEqual(refused.returncode, 0)
            self.assertIn("LocalPiecesOfParallelEdges", refused.stdout + refused.stderr)


class DefaultOffIdentity(unittest.TestCase):
    """Without --local-pieces every generator output is the one-reference form (byte-identical)."""

    def test_cluster_signature_csv_and_library_entry(self):
        edges = [{"Offset": o, "GapDirection": d, "Conductor": 1} for o, d in zip([0.0, 2.0, 4.0, 6.0], [1, -1, 1, -1])]
        with tempfile.TemporaryDirectory(prefix="pieces-off-") as tmp:
            tmp = Path(tmp)
            (tmp / "off").mkdir()
            (tmp / "on").mkdir()
            _, csv_off = PLANNER.write_cluster_signature(tmp / "off", edges)
            _, csv_default = PLANNER.write_cluster_signature(tmp / "on", edges, local_pieces=False)
            self.assertEqual(csv_off.read_text(), csv_default.read_text())
            self.assertEqual(csv_off.read_text().splitlines()[1:], [f"{o:.17g},{d},1" for o, d in zip([0.0, 2.0, 4.0, 6.0], [1, -1, 1, -1])])
            _, csv_on = PLANNER.write_cluster_signature(tmp / "on", edges, local_pieces=True)
            self.assertEqual([line.split(",")[2] for line in csv_on.read_text().splitlines()[1:]], ["1", "2", "2", "3"])
            layers = {"SA": (0.002, 4.0), "MS": (0.002, 11.45), "MA": (0.002, 10.0)}
            common = (1.9, 1055.0, 0.1, 0.05, 90.0, 0.0, 0.0, 11.45, layers, "k4", 1.0e-3)
            lib_off = CLUSTER.write_library(tmp / "off", "t", edges, [[3.0, 0.0, 0.0]], [], *common)
            model = json.loads(Path(lib_off).read_text())["Models"][0]
            self.assertNotIn("LocalPiece", json.dumps(model))
            self.assertEqual(model["Reference"], [3.0, 0.0, 0.0])
            refs = [[-0.95, 0.0, 0.0], [3.0, 0.0, 0.0], [6.95, 0.0, 0.0]]
            paths = [{"Indices": [1], "StartConductor": 1, "EndConductor": 2}, {"Indices": [2], "StartConductor": 2, "EndConductor": 3},
                     {"Indices": [3], "StartConductor": 3, "EndConductor": 1}]
            lib_on = CLUSTER.write_library(tmp / "on", "t", edges, refs, paths, *common, pieces=[1, 2, 2, 3])
            model = json.loads(Path(lib_on).read_text())["Models"][0]
            self.assertEqual([e["LocalPiece"] for e in model["Edges"]], [1, 2, 2, 3])
            self.assertEqual([e["Conductor"] for e in model["Edges"]], [1, 1, 1, 1])   # the key's labels stay
            self.assertEqual(model["ConductorReferences"], refs)
            self.assertEqual(model["OpenContourPaths"], paths)
            self.assertIn("LocalPiecesRule", model)
            self.assertEqual(edges[0], {"Offset": 0.0, "GapDirection": 1, "Conductor": 1})   # the input edges untouched


def planner_args(**overrides):
    """The planner options the two cache specs read (a SimpleNamespace stands for parse_args())."""
    base = dict(matching_radius=1.9, orders=[2], mesh_order=2, basis_size=6, samples=12, coupon_depth=1.0,
                edge_offset_tolerance=1.0e-3, straight_lc_fine=0.01, straight_lc_far=0.5, cluster_lc_fine=0.01,
                cluster_lc_far=0.5, min_process_feature_elements=4, cluster_h_factors=[2.0, 1.0])
    base.update(overrides)
    return types.SimpleNamespace(**base)


class PlannerCacheForm(unittest.TestCase):
    """Review d1099 MAJOR-1 / MINOR-3 (decision 1101): the planner's straight / parallel-cluster cache keys carry the
    piece form (no silent reuse of the other form's library from one cache root), the one-reference spec has no such
    key, and the piece form's model Name carries the -pieces<P> suffix. Since decision 1136 (4) the form is PLANNED per
    coupon (Preparation.LocalPieces, test_prepare_surface_response_coupons.PlannerLocalPiecesDefault): the specs read
    the planned form (or an explicit bool), not a CLI flag."""

    PARAMETERS = {"metal_thickness": 0.1, "overetch": 0.05}
    CLUSTER = {"Id": "pec", "Topology": "ParallelEdgeCluster", "BoundaryCondition": "Verified", "Preparation": {},
               "Geometry": {"Edges": [{"Offset": o, "GapDirection": d, "Conductor": 1}
                                      for o, d in zip([0.0, 2.0, 4.0, 6.0], [1, -1, 1, -1])]}}
    STRAIGHT = [{"Id": "gap", "Topology": "SameConductorGap", "BoundaryCondition": "Verified", "Geometry": {"Separation": 2.0},
                 "Preparation": {}},
                {"Id": "iso", "Topology": "IsolatedEdge", "BoundaryCondition": "Verified", "Geometry": {}, "Preparation": {}}]

    def test_cache_specs_key_on_the_form_and_the_one_reference_spec_is_unchanged(self):
        edges = PLANNER.normalize_parallel_edges(self.CLUSTER, 1.9)
        args = planner_args()
        for spec_of in (lambda form: PLANNER.parallel_cluster_spec(self.CLUSTER, edges, args, self.PARAMETERS, form),
                        lambda form: PLANNER.straight_spec(self.STRAIGHT, args, self.PARAMETERS, form)):
            spec_unplanned, spec_off, spec_on = spec_of(None), spec_of(False), spec_of(True)
            self.assertEqual(spec_unplanned, spec_off)                    # a plan without the form = the one-reference form
            self.assertNotIn("LocalPieces", spec_off)
            self.assertEqual(spec_on.pop("LocalPieces"), True)           # the form is the ONLY difference
            self.assertEqual(spec_on, spec_off)
            spec_on["LocalPieces"] = True
            self.assertNotEqual(PLANNER.fingerprint(spec_on), PLANNER.fingerprint(spec_off))   # distinct cache roots
        # the planned form (decision 1136 (4)): the one-label CPW (3 pieces) and the gap are planned in the piece form
        planned = json.loads(json.dumps(self.CLUSTER))
        PLANNER.plan_local_pieces_form(planned, 1.9)
        self.assertTrue(PLANNER.coupon_local_pieces(planned))
        self.assertIn("LocalPieces", PLANNER.parallel_cluster_spec(planned, edges, args, self.PARAMETERS))
        straight = json.loads(json.dumps(self.STRAIGHT))
        for coupon in straight:
            PLANNER.plan_local_pieces_form(coupon, 1.9)
        self.assertEqual([PLANNER.coupon_local_pieces(c) for c in straight], [True, False])
        self.assertTrue(PLANNER.straight_family_local_pieces(straight))
        self.assertIn("LocalPieces", PLANNER.straight_spec(straight, args, self.PARAMETERS))
        key = PLANNER.fingerprint(PLANNER.parallel_cluster_spec(self.CLUSTER, edges, args, self.PARAMETERS, False))
        self.assertEqual(PLANNER.parallel_cluster_model_name(key, edges, False), f"parallel-edge-cluster-{key[:12]}")
        key_on = PLANNER.fingerprint(PLANNER.parallel_cluster_spec(self.CLUSTER, edges, args, self.PARAMETERS, True))
        self.assertEqual(PLANNER.parallel_cluster_model_name(key_on, edges, True), f"parallel-edge-cluster-{key_on[:12]}-pieces3")

    def test_piece_form_model_name(self):
        self.assertEqual(lp.local_pieces_model_name("parallel-edge-cluster-abc", [1, 2, 2, 3]), "parallel-edge-cluster-abc-pieces3")
        self.assertEqual(lp.local_pieces_model_name("n", [1, 1]), "n-pieces1")
        # the generator's default Name (no --model-name): the one-reference form's Name + the suffix under --local-pieces
        spec = {"Edges": self.CLUSTER["Geometry"]["Edges"]}
        with tempfile.TemporaryDirectory(prefix="pieces-name-") as tmp:
            tmp = Path(tmp)
            (tmp / "spec.json").write_text(json.dumps(spec))
            names = {}
            for form, flag in (("off", []), ("on", ["--local-pieces"])):
                subprocess.run([sys.executable, str(CPW2D / "generate_edge_cluster_response.py"), str(tmp / "spec.json"),
                                "--output", str(tmp / form), "--thin-mesh", "thin.msh", "--fabricated-mesh", "fab.msh",
                                "--library-name", "t", *flag], check=True, capture_output=True, cwd=tmp)
                names[form] = json.loads((tmp / form / "process-library.json").read_text())["Models"][0]["Name"]
            self.assertRegex(names["off"], r"^parallel-edge-cluster-[0-9a-f]{12}$")
            self.assertEqual(names["on"], names["off"] + "-pieces3")


if __name__ == "__main__":
    unittest.main()


class StraightFamilyNameRule(unittest.TestCase):
    """Fix round R-10 (OI-STRAIGHT-FAMILY-NAME-COLLISION; decisions 1168 / 1230 / 1233 (5)): a NEW straight pair build's Name
    carries the requirement's 12-hex hash, so the 1168 twins 729ca8684d87 (2.9999993 um) / a11f82ccc415 (3.0000012 um) -- both
    "different-conductor-gap-3um" in the record form -- get DISTINCT Names; a build without a hash keeps the record-form Name
    (existing Names untouched); the planner passes the F-22 representative's hash per (topology, separation); the straight builder
    threads it to the pair generator."""

    def setUp(self):
        import generate_edge_pair_response as PAIR
        import build_surface_process_library as BUILDER
        self.PAIR, self.BUILDER = PAIR, BUILDER

    def test_pair_model_name_record_form_and_hash_form(self):
        self.assertEqual(self.PAIR.pair_model_name("different-conductor-gap", 2.9999993), "different-conductor-gap-3um")
        self.assertEqual(self.PAIR.pair_model_name("different-conductor-gap", 3.0000012), "different-conductor-gap-3um")   # the 1168 collision
        self.assertEqual(self.PAIR.pair_model_name("different-conductor-gap", 2.9999993, "729ca8684d87" + "0" * 52),
                         "different-conductor-gap-3um-729ca8684d87")
        self.assertEqual(self.PAIR.pair_model_name("different-conductor-gap", 3.0000012, "A11F82CCC415" + "0" * 52),
                         "different-conductor-gap-3um-a11f82ccc415")
        self.assertEqual(self.PAIR.pair_model_name("same-conductor-strip", 2.0), "same-conductor-strip-2um")
        with self.assertRaises(ValueError):
            self.PAIR.pair_model_name("same-conductor-gap", 1.0, "xyz")

    def test_write_library_name_with_and_without_the_hash(self):
        layers = {"SA": (0.002, 4.0), "MS": (0.002, 11.45), "MA": (0.002, 10.0)}
        with tempfile.TemporaryDirectory(prefix="r10-") as tmp:
            tmp = Path(tmp)
            names = []
            for sub, digest in (("a", None), ("b", "729ca8684d87" + "0" * 52), ("c", "a11f82ccc415" + "0" * 52)):
                (tmp / sub).mkdir()
                path = self.PAIR.write_library(tmp / sub, "t", 3.0, 1.0e-3, 1.9, 1055.0, True, False, [], 0.1, 11.45, layers,
                                               requirement_hash=digest)
                names.append(json.loads(Path(path).read_text())["Models"][0]["Name"])
            self.assertEqual(names, ["different-conductor-gap-3um", "different-conductor-gap-3um-729ca8684d87", "different-conductor-gap-3um-a11f82ccc415"])

    def test_planner_passes_the_f22_representative_hash_per_separation(self):
        coupons = [
            {"Topology": "DifferentConductorGap", "Geometry": {"Separation": 2.9999993}, "Hash": "729ca8684d87" + "0" * 52, "DeviceEdgeLength": 515.0},
            {"Topology": "DifferentConductorGap", "Geometry": {"Separation": 3.0000012}, "Hash": "a11f82ccc415" + "0" * 52, "DeviceEdgeLength": 6193.0},
            # two requirements of ONE separation float: the representative = length desc, hash asc (F-22)
            {"Topology": "SameConductorStrip", "Geometry": {"Separation": 2.0}, "Hash": "8389" + "0" * 60, "DeviceEdgeLength": 614.0},
            {"Topology": "SameConductorStrip", "Geometry": {"Separation": 2.0}, "Hash": "1111" + "0" * 60, "DeviceEdgeLength": 700.0},
            {"Topology": "SameConductorGap", "Geometry": {"Separation": 1.0}},   # a record-era coupon without a Hash: no option
            {"Topology": "IsolatedEdge", "Geometry": {}},
        ]
        hashes = PLANNER.straight_requirement_hashes(coupons)
        self.assertEqual(hashes[("different-conductor-gaps", 2.9999993)][:12], "729ca8684d87")
        self.assertEqual(hashes[("different-conductor-gaps", 3.0000012)][:12], "a11f82ccc415")
        self.assertEqual(hashes[("same-conductor-strips", 2.0)][:4], "1111")
        self.assertNotIn(("same-conductor-gaps", 1.0), hashes)
        options = PLANNER.straight_requirement_hash_options(coupons)
        self.assertEqual([o.split("=")[0] for o in options],
                         ["different-conductor-gap:2.9999993", "different-conductor-gap:3.0000012", "same-conductor-strip:2.0"])
        # the straight builder reads them back exactly (repr floats) and threads the hash to the pair generator
        args = types.SimpleNamespace(requirement_hash=options, matching_radius=1.9, metal_thickness=0.1, basis_size=6, samples=12, order=2,
                                     coupon_depth=1.0, separation_tolerance=1.0e-3, name="t", overetch=0.0, sidewall_angle=90.0, top_radius=0.0,
                                     bottom_radius=0.0, substrate_permittivity=11.45, sa_thickness=0.002, sa_permittivity=4.0, ms_thickness=0.002,
                                     ms_permittivity=11.45, ma_thickness=0.002, ma_permittivity=10.0, edge_distances=None)
        parsed = self.BUILDER.requirement_hashes(args)
        self.assertEqual(parsed[("different-conductor-gap", 2.9999993)][:12], "729ca8684d87")
        self.assertEqual(parsed[("different-conductor-gap", 3.0000012)][:12], "a11f82ccc415")
        command = self.BUILDER.pair_response_command(args, 3.0000012, "--different-conductors", Path("/tmp/x"),
                                                     parsed[("different-conductor-gap", 3.0000012)])
        self.assertEqual(command[command.index("--requirement-hash") + 1][:12], "a11f82ccc415")
        self.assertNotIn("--requirement-hash", self.BUILDER.pair_response_command(args, 1.0, None, Path("/tmp/x")))
        with self.assertRaises(ValueError):
            self.BUILDER.requirement_hashes(types.SimpleNamespace(requirement_hash=["corner:1.0=abc"]))
