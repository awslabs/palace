#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""U-PIECES-spatial (identification fix round R0-b; decisions 1100 / 1136 (2) / 1140 (1)): the SPATIAL Edges-only
LocalPiece form. (1) The arrangement rule by value (cluster_signature_geometry.local_pieces_of_arrangement: one piece
per metal face, the order (Conductor, smallest quantized vertex), the vertex-touch refusal by name); (2) the 18 R0-b
spatial candidates of the v3-j library of record (7 one-label multi-piece + 11 multi-piece-label; the census
probes/v3j_local_pieces_v1.out) rebuilt from their Signature with the pieces per label of the census, the default form
untouched; (3) the generator's --local-pieces (one reference per piece, OpenContourPaths, Edges[].LocalPiece with the
key labels, PlanViewBoundary unchanged) and its default-OFF identity; (4) the C++ library reader (the binary
$PALACE_EXECUTABLE, the S6 J1 preflight of record under the lanes' heavy lock) LOADS every candidate with Python's
Edges[].LocalPiece + one reference per piece, and refuses a wrong piece by name."""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cluster_signature_geometry as G  # noqa: E402
import device_coupons  # noqa: E402
import prepare_surface_response_coupons as planner  # noqa: E402

ASSESSMENT = Path(os.environ.get("COUPON_ASSESSMENT_ROOT", "/Users/simlap/bedrock-tests/coupon-accuracy-assessment-20260913"))
V3J = ASSESSMENT / "stage2-20261004/v3h/stretchfix/local/library/s2-r1p9-v3-j-local/process-library-local.json"
CENSUS = ASSESSMENT / "readiness-20261029/idfix/probes/v3j_local_pieces_v1.out"
S6_CONFIG = ASSESSMENT / "readiness-20261029/idfix/lane-B/local/dryrun-f8007abf2a/sct002-S6-thin-preflight-wi/lane/config.json"
PALACE = os.environ.get("PALACE_EXECUTABLE")
LOCK = Path("/tmp/idfix-heavy.lock")

ONE_LABEL = ("e7f44561bbf0", "d19e77e8fe71", "8e2c1b63713b", "4b0d7bf2ad1b", "17ee51c67d9a", "f6d7ca001864", "66217996479f")
TWO_LABEL = ("a2f53cd8630e", "e04e3eb5891d", "4388ea4d7f2a", "31b6896dc54d", "7d88636452b7", "5bda1bc70e63", "f813a57b5d43",
             "e3190c93a8ea", "09e95b2af0a4", "c983b84ada5c", "a287eb4daef0")
CANDIDATES = ONE_LABEL + TWO_LABEL
CONTEXT_ONLY_CONDUCTOR = ("98d99a3d918a", "bfe706240f6d")   # labels [1] + a context-only conductor 2 (one piece per label)


def segment(p0, p1, gap, conductor=1):
    return {"P0": np.asarray(p0, dtype=float), "P1": np.asarray(p1, dtype=float), "Gap": np.asarray(gap, dtype=float),
            "Conductor": conductor}


class ArrangementRule(unittest.TestCase):
    """local_pieces_of_arrangement by value on hand-built plan views (box (-5, -5) .. (5, 5), R = 1)."""

    BOX = ((-5.0, -5.0), (5.0, 5.0))

    def cpw(self, strip_conductor=1):
        # four chain segments spanning the box in y: ground | gap | strip | gap | ground (metal on the -gap side)
        return [segment((-3.0, -5.0), (-3.0, 5.0), (1.0, 0.0)), segment((-1.0, -5.0), (-1.0, 5.0), (-1.0, 0.0), strip_conductor),
                segment((1.0, -5.0), (1.0, 5.0), (1.0, 0.0), strip_conductor), segment((3.0, -5.0), (3.0, 5.0), (-1.0, 0.0))]

    def test_plan_view_faces_reads_the_arrangement(self):
        segments = self.cpw()
        faces = G.plan_view_arrangement(segments, self.BOX, 1.0)
        def plain(entries):
            return [([list(p) for p in polygon], metal, conductor) for polygon, metal, conductor in entries]
        self.assertEqual(plain(G.plan_view_faces(segments, self.BOX, 1.0)),
                         plain((face["Polygon"], face["Metal"], face["Conductor"]) for face in faces))
        self.assertEqual(sorted(face["Metal"] for face in faces), [False, False, True, True, True])
        for face in faces:
            self.assertEqual(sorted(set(face["Owners"])), face["Owners"])
            self.assertEqual(len(face["Nodes"]), len(face["Polygon"]))

    def test_one_label_cpw_is_three_pieces_in_increasing_x(self):
        faces = G.plan_view_arrangement(self.cpw(), self.BOX, 1.0)
        pieces, table = G.local_pieces_of_arrangement(faces, 4)
        self.assertEqual(pieces, [1, 2, 2, 3])
        self.assertEqual([entry["Piece"] for entry in table], [1, 2, 3])
        self.assertEqual([entry["Conductor"] for entry in table], [1, 1, 1])
        # the order: the smallest quantized vertex (x first) = the left ground, the strip, the right ground
        xs = [min(faces[entry["Face"]]["Polygon"], key=lambda p: (p[0], p[1]))[0] for entry in table]
        self.assertEqual(xs, sorted(xs))

    def test_the_order_is_conductor_first(self):
        faces = G.plan_view_arrangement(self.cpw(strip_conductor=2), self.BOX, 1.0)
        pieces, table = G.local_pieces_of_arrangement(faces, 4)
        self.assertEqual(pieces, [1, 3, 3, 2])                 # grounds (conductor 1) = pieces 1 / 2, the strip (2) = 3
        self.assertEqual([(entry["Conductor"], entry["Piece"]) for entry in table], [(1, 1), (1, 2), (2, 3)])

    def test_context_only_pieces_come_first_and_the_order_is_invariant_to_the_row_order(self):
        """Decision 1142 (amending 1100's order): the pieces bounded by no CLAIM row are numbered first, so the C++
        reader's count (the largest LocalPiece over the entry's claims) is the piece count; the numbering is geometric."""
        segments = self.cpw()
        pieces, table = G.local_pieces_of_arrangement(G.plan_view_arrangement(segments, self.BOX, 1.0), 4, claim_count=2)
        self.assertEqual(pieces, [2, 3, 3, 1])                 # the right ground (rows 3 only: context) first
        self.assertEqual([(entry["Piece"], entry["Claimed"]) for entry in table], [(1, False), (2, True), (3, True)])
        self.assertEqual(max(pieces[:2]), 3)                   # the claims' largest piece = P
        with self.assertRaisesRegex(G.SignatureGeometryError, "0 claims of 4 rows"):
            G.local_pieces_of_arrangement(G.plan_view_arrangement(segments, self.BOX, 1.0), 4, claim_count=0)
        for order in ((3, 2, 1, 0), (1, 3, 0, 2), (2, 0, 3, 1)):
            permuted = [segments[i] for i in order]
            again, _ = G.local_pieces_of_arrangement(G.plan_view_arrangement(permuted, self.BOX, 1.0), 4)
            reference, _ = G.local_pieces_of_arrangement(G.plan_view_arrangement(segments, self.BOX, 1.0), 4)
            self.assertEqual(again, [reference[i] for i in order], order)   # the piece of every row is invariant

    def test_a_metal_face_without_an_edge_row_is_refused_by_name(self):
        # the last two segments are not rows (a bridge, say): the right ground is a piece with no edge
        with self.assertRaisesRegex(G.SignatureGeometryError, "bounded by no Edge row"):
            G.local_pieces_of_arrangement(G.plan_view_arrangement(self.cpw(), self.BOX, 1.0), 2)

    def test_two_metal_pieces_touching_at_a_vertex_are_refused_by_name(self):
        # an hourglass: the left and right metal triangles meet at the apex (gap vectors by sign only)
        def hourglass(apex):
            return [segment((-5.0, -5.0), (0.0, 0.0), (1.0, -1.0)), segment((0.0, 0.0), (-5.0, 5.0), (1.0, 1.0)),
                    segment((5.0, 5.0), apex, (-1.0, 1.0)), segment(apex, (5.0, -5.0), (-1.0, -1.0))]
        touching = hourglass((0.0, 0.0))
        faces = G.plan_view_arrangement(touching, self.BOX, 1.0)
        self.assertEqual(sum(face["Metal"] for face in faces), 2)
        with self.assertRaisesRegex(G.SignatureGeometryError, "metal pieces touch at a vertex.*decision 1100"):
            G.local_pieces_of_arrangement(faces, len(touching))
        apart = hourglass((1.0, 0.0))
        pieces, table = G.local_pieces_of_arrangement(G.plan_view_arrangement(apart, self.BOX, 1.0), len(apart))
        self.assertEqual(pieces, [1, 1, 2, 2])
        self.assertEqual(len(table), 2)


def v3j_spatial_models():
    library = json.loads(V3J.read_text())
    return library, {model["Name"].split("_")[-1]: model for model in library["Models"] if model.get("Topology") == "SpatialEdgeCluster"}


def record_of(model):
    return {"Topology": "SpatialEdgeCluster", "Geometry": {"Signature": model["Signature"]}, "Interfaces": model["Interfaces"],
            "BoundaryCondition": model.get("BoundaryCondition", {"Type": "PEC"})}


def census_pieces():
    """{key: {label: pieces}} of the R0-b candidate census (probes/v3j_local_pieces_v1.out)."""
    out = {}
    for line in CENSUS.read_text().splitlines():
        if line.startswith("SpatialEdgeCluster") and "pieces {" in line:
            key = line.split()[1].split("_")[-1]
            out[key] = {int(k): int(v) for k, v in (item.split(":") for item in
                                                   line[line.index("pieces {") + 8:line.index("}")].replace(" ", "").split(","))}
    return out


def strip_pieces(value):
    if isinstance(value, dict):
        return {k: strip_pieces(v) for k, v in value.items() if k not in ("LocalPiece", "LocalPieces", "LocalPiecesRule")}
    if isinstance(value, list):
        return [strip_pieces(v) for v in value]
    return value


@unittest.skipUnless(V3J.is_file() and CENSUS.is_file(), "the v3-j local mirror and the R0-b census are needed")
class V3jCandidates(unittest.TestCase):
    def test_every_candidate_rebuilds_with_the_census_pieces_and_the_default_form_is_untouched(self):
        library, spatial = v3j_spatial_models()
        radius, fab = float(library["MatchingRadius"]), library["Fabrication"]
        census = census_pieces()
        for key in CANDIDATES:
            model = spatial[key]
            default, default_edges = G.cluster_coupon(record_of(model), radius, fab["MetalThickness"], fab["OveretchDepth"])
            pieces, model_edges = G.cluster_coupon(record_of(model), radius, fab["MetalThickness"], fab["OveretchDepth"],
                                                   local_pieces=True)
            self.assertNotIn("LocalPiece", json.dumps([default, default_edges]), key)
            self.assertEqual(strip_pieces([pieces, model_edges]), [default, default_edges], key)   # the default form, plus the leaf
            per_label = {}
            for row in pieces["Geometry"]["Edges"]:
                per_label.setdefault(row["Conductor"], set()).add(row["LocalPiece"])
            self.assertEqual({label: len(found) for label, found in per_label.items()}, census[key], key)
            count = pieces["Geometry"]["LocalPieces"]
            self.assertEqual(count, sum(census[key].values()), key)
            self.assertEqual({row["LocalPiece"] for row in pieces["Geometry"]["Edges"]}, set(range(1, count + 1)), key)
            self.assertEqual({facet["LocalPiece"] for facet in pieces["Geometry"]["PlanViewFacets"]}, set(range(1, count + 1)), key)
            # the model edges = the claims rows (segments[i] <-> rows[i]) with their pieces
            self.assertEqual([edge["LocalPiece"] for edge in model_edges],
                             [row["LocalPiece"] for row in pieces["Geometry"]["Edges"][:len(model_edges)]], key)
            self.assertEqual([edge["Point"] for edge in model_edges], [edge["Point"] for edge in model["Edges"]], key)
            self.assertEqual(pieces["Geometry"]["PlanViewBoundary"], default["Geometry"]["PlanViewBoundary"], key)
            self.assertEqual(pieces["Geometry"]["LocalPiecesRule"], G.LOCAL_PIECES_RULE)
        self.assertEqual(len(ONE_LABEL), 7)
        self.assertEqual(len(TWO_LABEL), 11)


def generate(coupon, output, radius, fab, name, local_pieces):
    output.mkdir(parents=True, exist_ok=True)
    (output / "coupon.json").write_text(json.dumps(coupon, indent=2) + "\n")
    command = [sys.executable, str(HERE / "generate_spatial_response.py"), str(output / "coupon.json"), "--output", str(output),
               "--radius", str(radius), "--metal-thickness", str(fab["MetalThickness"]), "--overetch-depth", str(fab["OveretchDepth"]),
               "--sidewall-angle", str(fab["SidewallAngleDegrees"]), "--top-rounding", "0", "--trench-rounding", "0",
               "--ring-size", "16", "--cap-triangulation", "delaunay", "--cap-interior-spacing", "1.0", "--order", "1",
               "--model-name", name, "--basis-only", *(["--local-pieces"] if local_pieces else [])]
    return subprocess.run(command, capture_output=True, text=True)


def tree_digests(root):
    out = {}
    for path in sorted(p for p in root.rglob("*") if p.is_file() and p.name != "coupon.json"):
        out[str(path.relative_to(root))] = __import__("hashlib").sha256(path.read_bytes()).hexdigest()
    return out


@unittest.skipUnless(V3J.is_file(), "the v3-j local mirror is needed")
class GeneratorPieceForm(unittest.TestCase):
    """generate_spatial_response.py --local-pieces on the 8-edge one-label candidate e7f44561bbf0 (3 pieces; basis-only):
    the piece form's entry, the default-OFF identity (the same coupon without the leaf = the same files), the refusals."""

    def test_piece_form_entry_default_identity_and_refusals(self):
        library, spatial = v3j_spatial_models()
        radius, fab = float(library["MatchingRadius"]), library["Fabrication"]
        model = spatial["e7f44561bbf0"]
        default, _ = G.cluster_coupon(record_of(model), radius, fab["MetalThickness"], fab["OveretchDepth"])
        pieces, _ = G.cluster_coupon(record_of(model), radius, fab["MetalThickness"], fab["OveretchDepth"], local_pieces=True)
        with tempfile.TemporaryDirectory(prefix="spatial-pieces-") as tmp:
            tmp = Path(tmp)
            off = generate(default, tmp / "off", radius, fab, model["Name"], False)
            self.assertEqual(off.returncode, 0, off.stderr[-2000:])
            # the default form on a coupon carrying the leaf (the flag OFF) = the default form: byte-identical files
            off_leaf = generate(pieces, tmp / "off-leaf", radius, fab, model["Name"], False)
            self.assertEqual(off_leaf.returncode, 0, off_leaf.stderr[-2000:])
            self.assertEqual(tree_digests(tmp / "off"), tree_digests(tmp / "off-leaf"))
            on = generate(pieces, tmp / "on", radius, fab, model["Name"] + "-pieces3", True)
            self.assertEqual(on.returncode, 0, on.stderr[-2000:])
            entry_off = json.loads((tmp / "off" / "process-library.json").read_text())["Models"][0]
            entry_on = json.loads((tmp / "on" / "process-library.json").read_text())["Models"][0]
            self.assertNotIn("LocalPiece", json.dumps(entry_off))
            self.assertEqual(entry_off["Reference"], model["Reference"])                  # the record's one reference
            self.assertEqual(len(entry_on["ConductorReferences"]), 3)
            self.assertNotIn("Reference", entry_on)
            self.assertEqual(entry_on["ConductorReferences"][0], model["Reference"])      # piece 1 = the label's reference
            self.assertTrue(entry_on["OpenContourPaths"])
            self.assertNotIn("ContourGroups", entry_on)
            self.assertEqual([edge["LocalPiece"] for edge in entry_on["Edges"]][:8], [3, 2, 3, 2, 1, 2, 1, 2])
            self.assertEqual([edge["Conductor"] for edge in entry_on["Edges"]], [edge["Conductor"] for edge in entry_off["Edges"]])
            self.assertEqual(entry_on["PlanViewBoundary"], entry_off["PlanViewBoundary"])   # the label form (string-matched by the C++)
            self.assertEqual(entry_on["LocalPiecesRule"], __import__("generate_spatial_response").LOCAL_PIECES_RULE)
            self.assertEqual(entry_on["Name"], model["Name"] + "-pieces3")
            contract = json.loads((tmp / "on" / "basis-contract.json").read_text())
            self.assertEqual(contract["ConductorStates"], 2)
            self.assertEqual(sorted(p.name for p in (tmp / "on").glob("conductor-*.csv")), ["conductor-2.csv", "conductor-3.csv"])
            # the mesher's CSVs are tagged per piece
            signature_labels = {line.split(",")[2] for line in (tmp / "on" / "mesh-signature.csv").read_text().splitlines()[1:]}
            self.assertEqual(signature_labels, {"1", "2", "3"})
            boundary_labels = {line.split(",")[2] for line in (tmp / "on" / "plan-view-boundary.csv").read_text().splitlines()[1:]}
            self.assertEqual(boundary_labels, {"1", "2", "3"})
            # --local-pieces on a coupon without the leaf: refused by name, generation-failure.json written
            refused = generate(default, tmp / "refused", radius, fab, model["Name"], True)
            self.assertNotEqual(refused.returncode, 0)
            self.assertIn("carries no positive LocalPiece", refused.stderr)
            self.assertEqual(json.loads((tmp / "refused" / "generation-failure.json").read_text())["Stage"], "SpatialGeometry")

    def test_device_path_stamps_the_exact_edges_with_their_pieces(self):
        library, spatial = v3j_spatial_models()
        radius, parameters = float(library["MatchingRadius"]), planner.process_parameters(library)
        model = spatial["66217996479f"]                     # 4 claims, 2 pieces
        coupon = {"Id": model["Name"], **record_of(model)}
        with tempfile.TemporaryDirectory(prefix="spatial-device-") as tmp:
            for form in (False, True):
                work = Path(tmp) / ("pieces" if form else "one")
                device_coupons.generate_sources(coupon, work, radius=radius, parameters=parameters, ring_size=16, local_pieces_form=form)
                entry = json.loads((work / "process-library.json").read_text())["Models"][0]
                self.assertEqual(entry["Signature"], model["Signature"])
                self.assertEqual(len(entry["Edges"]), 4)                                     # the exact claims
                self.assertEqual([e["Point"] for e in entry["Edges"]], [e["Point"] for e in model["Edges"]])
                if form:
                    self.assertEqual(entry["Name"], model["Name"] + "-pieces2")
                    self.assertEqual([e["LocalPiece"] for e in entry["Edges"]], [1, 2, 2, 1])
                    self.assertEqual(len(entry["ConductorReferences"]), 2)
                    self.assertIn("LocalPiecesRule", entry)
                else:
                    self.assertEqual(entry["Name"], model["Name"])
                    self.assertNotIn("LocalPiece", json.dumps(entry))
                    self.assertEqual(entry["Reference"], model["Reference"])
            with self.assertRaisesRegex(device_coupons.DeviceAdapterError, "version-2 signature record"):
                device_coupons.coupon_geometry({"Id": "v1", "Topology": "SpatialEdgeCluster",
                                                "Geometry": {"Edges": [{}], "PlanViewFacets": []}}, radius, parameters, local_pieces=True)


def piece_form_models(models, radius, fab, keys):
    """Every v3-j spatial model named in `keys` with Python's Edges[].LocalPiece and one reference per piece (a point on
    the metal side of the piece's first claim edge, in the model's canonical frame) + synthetic OpenContourPaths
    (the reader requires open paths connecting every reference; the preflight reads no matrices); the others
    untouched. Returns (models, {key: piece count})."""
    out, counts = [], {}
    for model in models:
        key = model["Name"].split("_")[-1]
        if model.get("Topology") != "SpatialEdgeCluster" or key not in keys:
            out.append(model)
            continue
        model = json.loads(json.dumps(model))
        _, model_edges = G.cluster_coupon(record_of(model), radius, fab["MetalThickness"], fab["OveretchDepth"], local_pieces=True)
        coupon, model_edges = G.cluster_coupon(record_of(model), radius, fab["MetalThickness"], fab["OveretchDepth"],
                                               local_pieces=True)
        assert [e["Point"] for e in model_edges] == [e["Point"] for e in model["Edges"]]
        for edge, exact in zip(model["Edges"], model_edges):
            edge["LocalPiece"] = exact["LocalPiece"]
        references = {}
        for row in coupon["Geometry"]["Edges"]:          # claims + context rows: every piece has a row (1142)
            if row["LocalPiece"] not in references:
                point = np.asarray(row["Point"]) - 0.1 * radius * np.asarray(row["GapDirection"])
                references[row["LocalPiece"]] = [float(v) for v in point]
        count = max(references)
        assert sorted(references) == list(range(1, count + 1)) and count == coupon["Geometry"]["LocalPieces"]
        assert max(e["LocalPiece"] for e in model_edges) == count   # the reader's count rule (1142)
        model.pop("Reference", None)
        model.pop("ContourGroups", None)
        model.pop("ZeroTraceIndices", None)
        model["ConductorReferences"] = [references[p] for p in range(1, count + 1)]
        model["OpenContourPaths"] = [{"Indices": [k + 1], "StartConductor": k + 1, "EndConductor": k + 2} for k in range(count - 1)]
        counts[key] = count
        out.append(model)
    return out, counts


def with_heavy_lock(command, timeout):
    """One palace preflight under the lanes' heavy-process lock (decision 1062; mkdir atomic; released after)."""
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
        (LOCK / "owner").write_text(f"r0b {os.getpid()} unittest {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\n")
        return subprocess.run(command, capture_output=True, text=True, timeout=timeout)
    finally:
        shutil.rmtree(LOCK, ignore_errors=True)


@unittest.skipUnless(PALACE and os.path.isfile(PALACE) and V3J.is_file() and S6_CONFIG.is_file(),
                     "PALACE_EXECUTABLE / the v3-j mirror / the S6 J1 config of record are needed")
class CxxReaderCrossCheckSpatial(unittest.TestCase):
    """U-PIECES-spatial: the C++ library reader (surfaceresponseoperator.cpp: Edges[].LocalPiece on every spatial edge,
    one ConductorReference per piece) LOADS every one of the 18 candidates in Python's piece form through the S6 J1
    preflight of record, with the requirements row-identical to the base; a wrong piece (a piece index beyond the
    references) and a missing LocalPiece are refused by name."""

    def preflight(self, tmp, library_models, label):
        library = json.loads(V3J.read_text())
        library["Models"] = library_models
        library["Name"] = library["Name"] + "-" + label
        lib_path = tmp / f"library-{label}.json"
        lib_path.write_text(json.dumps(library))
        config = json.loads(S6_CONFIG.read_text())
        config["Problem"]["Output"] = str(tmp / label / "postpro")
        config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] = str(lib_path)
        config_path = tmp / f"config-{label}.json"
        config_path.write_text(json.dumps(config))
        self.assertNotIn("lane-B/local/dryrun", config_path.read_text().replace(str(S6_CONFIG), ""))   # re-pointed (838)
        return with_heavy_lock([PALACE, "--surface-response-preflight", str(config_path)], timeout=280)

    def test_every_candidate_loads_in_the_piece_form_and_wrong_pieces_are_refused(self):
        library, spatial = v3j_spatial_models()
        radius, fab = float(library["MatchingRadius"]), library["Fabrication"]
        with tempfile.TemporaryDirectory(prefix="upieces-spatial-") as tmp:
            tmp = Path(tmp)
            base = self.preflight(tmp, library["Models"], "base")
            self.assertEqual(base.returncode, 0, base.stdout[-2000:] + base.stderr[-2000:])
            # the 18 candidates + the two one-piece-per-label models with a CONTEXT-ONLY conductor (decision 1142)
            models, counts = piece_form_models(library["Models"], radius, fab, set(CANDIDATES) | set(CONTEXT_ONLY_CONDUCTOR))
            self.assertEqual(len(counts), 20)
            self.assertEqual({counts[k] for k in CONTEXT_ONLY_CONDUCTOR}, {2})
            self.assertEqual(counts["e7f44561bbf0"], 3)
            self.assertEqual(counts["4b0d7bf2ad1b"], 5)
            self.assertEqual({counts[k] for k in TWO_LABEL}, {3})
            lane = self.preflight(tmp, models, "pieces")
            self.assertEqual(lane.returncode, 0, lane.stdout[-3000:] + lane.stderr[-3000:])

            def requirements(label):
                manifest = json.loads((tmp / label / "postpro" / "surface-response-requirements.json").read_text())
                return ([(r["Hash"], r["Status"], r.get("TotalEdgeLength")) for r in manifest["Requirements"]],
                        manifest["Summary"]["Counts"])
            self.assertEqual(requirements("base"), requirements("pieces"))
            # a piece index beyond the references on one edge: refused by name (the reader's count rule)
            wrong = json.loads(json.dumps(models))
            target = next(m for m in wrong if m["Name"].endswith("e7f44561bbf0"))
            target["Edges"][0]["LocalPiece"] = counts["e7f44561bbf0"] + 1
            refused = self.preflight(tmp, wrong, "wrong")
            self.assertNotEqual(refused.returncode, 0)
            self.assertIn("one reference per local piece is required", refused.stdout + refused.stderr)
            # a missing LocalPiece on one edge of a piece-form model: refused by name
            partial = json.loads(json.dumps(models))
            target = next(m for m in partial if m["Name"].endswith("d19e77e8fe71"))
            del target["Edges"][-1]["LocalPiece"]
            refused = self.preflight(tmp, partial, "partial")
            self.assertNotEqual(refused.returncode, 0)
            self.assertIn("must be present on every edge", refused.stdout + refused.stderr)


if __name__ == "__main__":
    unittest.main()
