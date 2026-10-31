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
    --local-pieces form (no silent reuse of the other form's library from one cache root), the default spec has no such
    key, and the piece form's model Name carries the -pieces<P> suffix."""

    PARAMETERS = {"metal_thickness": 0.1, "overetch": 0.05}
    CLUSTER = {"Id": "pec", "Topology": "ParallelEdgeCluster", "BoundaryCondition": "Verified",
               "Geometry": {"Edges": [{"Offset": o, "GapDirection": d, "Conductor": 1}
                                      for o, d in zip([0.0, 2.0, 4.0, 6.0], [1, -1, 1, -1])]}}
    STRAIGHT = [{"Id": "gap", "Topology": "SameConductorGap", "BoundaryCondition": "Verified", "Geometry": {"Separation": 2.0}},
                {"Id": "iso", "Topology": "IsolatedEdge", "BoundaryCondition": "Verified", "Geometry": {}}]

    def test_cache_specs_key_on_the_form_and_the_default_is_unchanged(self):
        edges = PLANNER.normalize_parallel_edges(self.CLUSTER, 1.9)
        absent, off, on = planner_args(), planner_args(local_pieces=False), planner_args(local_pieces=True)
        for spec_of in (lambda a: PLANNER.parallel_cluster_spec(self.CLUSTER, edges, a, self.PARAMETERS),
                        lambda a: PLANNER.straight_spec(self.STRAIGHT, a, self.PARAMETERS)):
            spec_absent, spec_off, spec_on = spec_of(absent), spec_of(off), spec_of(on)
            self.assertEqual(spec_absent, spec_off)                       # the default spec: the flag absent == off
            self.assertNotIn("LocalPieces", spec_absent)
            self.assertEqual(spec_on.pop("LocalPieces"), True)           # the form is the ONLY difference
            self.assertEqual(spec_on, spec_off)
            spec_on["LocalPieces"] = True
            self.assertNotEqual(PLANNER.fingerprint(spec_on), PLANNER.fingerprint(spec_off))   # distinct cache roots
        key = PLANNER.fingerprint(PLANNER.parallel_cluster_spec(self.CLUSTER, edges, off, self.PARAMETERS))
        self.assertEqual(PLANNER.parallel_cluster_model_name(key, edges, off), f"parallel-edge-cluster-{key[:12]}")
        key_on = PLANNER.fingerprint(PLANNER.parallel_cluster_spec(self.CLUSTER, edges, on, self.PARAMETERS))
        self.assertEqual(PLANNER.parallel_cluster_model_name(key_on, edges, on), f"parallel-edge-cluster-{key_on[:12]}-pieces3")

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
