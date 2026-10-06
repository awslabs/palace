# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""nearkey_detection / nearkey_transplant / nearkey_reuse: the structure-key pin (the calibrated pair's key
is admitted; isomorphic-looking 19-edge clusters, edge-count / topology changes and the transmon's 3- / 4-edge
clusters are refused StructureKeyNotCalibrated), every fail-closed refusal of DESIGN v2 section 1.4 on the
fixture pair (B2 d7a22d38cefe <- S1p box 1), the decision-431 donor-override rule, the metrics of the pair
(golden values of the design evidence), a synthetic transplant fixture (identity / shift / orphan), the
record schema of a written reused model, the library assembly, and - when the local evidence mirror is
mounted - the six measured pairs' transplanted matrices byte-identical to the sensitivity lane's libraries."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import device_coupons  # noqa: E402
import nearkey_detection as detection  # noqa: E402
import nearkey_predictor as predictor  # noqa: E402
import nearkey_reuse as reuse  # noqa: E402
import nearkey_transplant as transplant  # noqa: E402
from test_nearkey_predictor import activated  # noqa: E402

FIXTURE = json.loads((HERE / "testdata" / "nearkey-calibrated-pair.json").read_text())
RADIUS = float(FIXTURE["MatchingRadius"])
S1P_BOX1_KEY = "2d27be59ff1b0774e86fe766572cda6cecd58c3c9dd14865b201909b51c32caa"
S2 = Path("/Users/simlap/bedrock-tests/coupon-accuracy-assessment-20260913/stage2-20261004")
LIBRARY_OF_RECORD = S2 / "coupons-bc" / "library" / "s2-r1p9-v3-b1" / "process-library.json"
MODELS_MIRROR = S2 / "nearkey-sensitivity" / "work" / "models"
SENS_LIBRARY = S2 / "nearkey-sensitivity" / "library" / "s2-r1p9-v3-b-sens-c3"
SHELLED = S2 / "nearkey-reuse-impl" / "work" / "shelled"
# the design evidence (pair-geometry.json / transplant-report.json): the six pairs' golden metrics and donors
SIX_PAIRS = {"B2": ("078-spatialedgecluster-edgecount-19-5ae3dbdd2c3d", "spatialedgecluster_edgecount-19_dcfd1e0fe3fa",
                    -0.009527170077628535, 0.0007573071691547667),
             "B3": ("081-spatialedgecluster-edgecount-19-1de0718f1fc8", "spatialedgecluster_edgecount-19_dcfd1e0fe3fa",
                    0.014625330182121218, 0.001268858969695948),
             "B4": ("079-spatialedgecluster-edgecount-19-7d88636452b7", "spatialedgecluster_edgecount-19_a2f53cd8630e",
                    0.003123125975976694, 0.00020145052418115068),
             "B5": ("080-spatialedgecluster-edgecount-19-5bda1bc70e63", "spatialedgecluster_edgecount-19_a2f53cd8630e",
                    0.015825391576225978, 0.00101035359560976),
             "C1": ("076-spatialedgecluster-edgecount-19-6b38cd966bbc", "spatialedgecluster_edgecount-19_130f26d56c23",
                    -0.06064641022620769, 0.0073833166876292955),
             "C2": ("077-spatialedgecluster-edgecount-20-31b6896dc54d", "spatialedgecluster_edgecount-20_e04e3eb5891d",
                    -0.11280045501005607, 0.006391354806027371)}


def rule():
    return predictor.load_rule()


def analyse(exact=None, donor=None, **kwargs):
    return detection.analyse_pair(exact or FIXTURE["Exact"]["Signature"], donor or FIXTURE["Donor"], rule=rule(), radius=RADIUS, **kwargs)


def shifted(signature, dx=0.0, dy=0.0, portions=None, context=None):
    """The signature with the given portions / context pieces translated (R units)."""
    out = copy.deepcopy(signature)
    for i in portions or ():
        p = out["Portions"][i]["P"]
        out["Portions"][i]["P"] = [p[0] + dx, p[1] + dy, p[2] + dx, p[3] + dy]
    for i in context or ():
        p = out["Context"][i]["P"]
        out["Context"][i]["P"] = [p[0] + dx, p[1] + dy, p[2] + dx, p[3] + dy]
    return out


def donor_with(signature=None, **fields):
    donor = copy.deepcopy(FIXTURE["Donor"])
    if signature is not None:
        donor["Signature"] = signature
    donor.update(fields)
    return donor


class StructureKey(unittest.TestCase):
    def test_calibrated_pair_shares_an_admissible_key(self):
        key_e, _, topology_e = detection.structure_key(FIXTURE["Exact"]["Signature"])
        key_d, _, topology_d = detection.structure_key(FIXTURE["Donor"]["Signature"])
        self.assertEqual(key_e, S1P_BOX1_KEY)
        self.assertEqual(key_d, S1P_BOX1_KEY)
        self.assertEqual(topology_e, topology_d)
        self.assertIn(key_e, predictor.admissible_structure_keys(rule()))

    def test_other_structures_are_refused_structure_key_not_calibrated(self):
        admissible = predictor.admissible_structure_keys(rule())
        cases = {"transmon 3-edge": FIXTURE["Transmon3Edge"]["Signature"], "transmon 4-edge": FIXTURE["Transmon4Edge"]["Signature"],
                 "isomorphic-looking 19-edge (14 context pieces)": FIXTURE["Isomorphic19Edge"]["Signature"]}
        eighteen = copy.deepcopy(FIXTURE["Exact"]["Signature"])   # the same cluster with one portion fewer
        eighteen["Portions"] = eighteen["Portions"][:-1]
        eighteen["EdgeCount"] = 18
        cases["same cluster, 18 edges"] = eighteen
        three_conductor = copy.deepcopy(FIXTURE["Exact"]["Signature"])
        three_conductor["Portions"][5]["Conductor"] = 3
        cases["three conductors"] = three_conductor
        for label, signature in cases.items():
            key = detection.structure_key(signature)[0]
            self.assertNotIn(key, admissible, label)
            records = detection.candidate_donors(signature, {"MatchingRadius": RADIUS, "Models": [FIXTURE["Donor"]]}, rule=rule())
            self.assertEqual(list(records), ["*"], label)
            self.assertFalse(records["*"]["Qualifies"])
            self.assertTrue(records["*"]["Refusals"][0].startswith(detection.REFUSAL_STRUCTURE_KEY), (label, records["*"]["Refusals"]))
            record = analyse(exact=signature)
            self.assertTrue(record["Refusals"][0].startswith(detection.REFUSAL_STRUCTURE_KEY), label)

    def test_donor_of_another_key_is_refused(self):
        record = analyse(donor=FIXTURE["OtherFamilyDonor"])
        self.assertFalse(record["Qualifies"])
        self.assertTrue(record["Refusals"][0].startswith("StructureKeyMismatch"))
        record = analyse(donor={"Name": "no-signature", "Topology": "SpatialEdgeCluster"})
        self.assertTrue(record["Refusals"][0].startswith("StructureKeyMismatch"))
        record = analyse(exact={"Type": "IsolatedEdge"})
        self.assertTrue(record["Refusals"][0].startswith(detection.REFUSAL_STRUCTURE_KEY))

    @unittest.skipUnless(LIBRARY_OF_RECORD.is_file(), "the library of record is not mounted")
    def test_every_other_cluster_of_the_library_of_record_is_refused(self):
        library = json.loads(LIBRARY_OF_RECORD.read_text())
        admissible = predictor.admissible_structure_keys(rule())
        calibrated = {"dcfd1e0fe3fa", "5ae3dbdd2c3d", "1de0718f1fc8", "a2f53cd8630e", "7d88636452b7", "5bda1bc70e63",
                      "130f26d56c23", "6b38cd966bbc", "e04e3eb5891d", "31b6896dc54d"}
        keys = {}
        for model in library["Models"]:
            if model.get("Topology") != "SpatialEdgeCluster":
                continue
            key = detection.structure_key(model["Signature"])[0]
            if model["Name"].split("_")[-1] in calibrated:
                self.assertIn(key, admissible, model["Name"])
                keys.setdefault(key, set()).add(model["Name"].split("_")[-1])
            else:
                self.assertNotIn(key, admissible, model["Name"])
        self.assertEqual(len(keys), 4)
        self.assertEqual(keys[S1P_BOX1_KEY], {"dcfd1e0fe3fa", "5ae3dbdd2c3d", "1de0718f1fc8"})


class Detection(unittest.TestCase):
    def test_fixture_pair_metrics_are_the_design_values(self):
        record = analyse()
        self.assertTrue(record["Qualifies"], record["Refusals"])
        self.assertTrue(record["DefaultAdmissible"])   # the decision-292 cut-end override is the admitted kind (decision 431)
        self.assertAlmostEqual(record["W"], -0.009527170077628535, places=15)
        self.assertAlmostEqual(record["S"], 0.0007573071691547667, places=15)
        self.assertAlmostEqual(record["VmaxNm"], 28.447, places=2)
        self.assertAlmostEqual(record["CmaxNm"], 28.447, places=2)
        self.assertAlmostEqual(record["BmaxNm"], 4.0185, places=3)
        self.assertAlmostEqual(record["LeadWidthUm"]["Exact"], 0.210, places=3)
        self.assertAlmostEqual(record["LeadWidthUm"]["Donor"], 0.208, places=3)
        self.assertLess(abs(record["GapScalingResidual"]), 1e-6)
        self.assertEqual(record["PortionCorrespondence"]["Method"], "ByIndex")
        self.assertEqual(record["QuantumDifference"]["MaxQuanta"], 14822.000000000113)
        self.assertEqual(record["DonorBuildGateOverride"]["Kind"], "CutNeighborhoods")
        self.assertEqual(record["DonorBuildGateOverride"]["SupersededBy"], 316)

    def test_quantum_near_match_is_not_a_near_key(self):
        record = analyse(donor=donor_with(signature=FIXTURE["Exact"]["Signature"]))
        self.assertFalse(record["Qualifies"])
        self.assertTrue(record["Refusals"][0].startswith("QuantumNearMatch"))

    def test_turned_direction_arc_and_gap_changes_are_refused(self):
        signature = FIXTURE["Donor"]["Signature"]
        turned = copy.deepcopy(signature)   # the same topology key: reverse one portion's direction (P swapped end for end)
        p = turned["Portions"][0]["P"]
        turned["Portions"][0]["P"] = [p[2], p[3], p[0], p[1]]
        record = analyse(donor=donor_with(signature=turned))
        self.assertFalse(record["Qualifies"])
        self.assertTrue(any(r.startswith("PortionCorrespondence") for r in record["Refusals"]), record["Refusals"])
        # the gap normal is nulled in the structure key (a length-class entry): a turned gap normal is item 2's refusal
        gap_turned = copy.deepcopy(signature)
        gap_turned["Portions"][0]["Gap"] = [-g for g in gap_turned["Portions"][0]["Gap"]]
        record = analyse(donor=donor_with(signature=gap_turned))
        self.assertFalse(record["Qualifies"])
        self.assertTrue(any(r.startswith("PortionCorrespondence") for r in record["Refusals"]), record["Refusals"])
        arc = copy.deepcopy(signature)
        arc["Context"][0]["Arc"] = 0.5
        record = analyse(donor=donor_with(signature=arc))
        self.assertFalse(record["Qualifies"])
        self.assertTrue(any(r.startswith("Arc") or r.startswith("StructureKeyMismatch") for r in record["Refusals"]), record["Refusals"])
        # the facing gap moved without the lead width (the junction gap no longer scales with the lead)
        exact = FIXTURE["Exact"]["Signature"]
        gap = detection.facing_gap(exact["Portions"])
        self.assertIsNotNone(gap)
        moved = copy.deepcopy(exact)
        # shift the conductor-2 portions facing conductor 1 across the gap by 2 % of the gap
        for i, portion in enumerate(moved["Portions"]):
            if portion["Conductor"] == 2:
                dx = -0.02 * gap * portion["Gap"][0]
                dy = -0.02 * gap * portion["Gap"][1]
                moved["Portions"][i]["P"] = [portion["P"][0] + dx, portion["P"][1] + dy, portion["P"][2] + dx, portion["P"][3] + dy]
        record = analyse(donor=donor_with(signature=moved))
        self.assertFalse(record["Qualifies"])
        self.assertTrue(any("does not scale with the lead" in r or "not one width" in r or "narrowest strip" in r
                            for r in record["Refusals"]),
                        record["Refusals"])

    def test_displacements_outside_the_domain_are_refused(self):
        exact = FIXTURE["Exact"]["Signature"]
        far_vertex = shifted(exact, dx=0.06, portions=[0])          # a claim vertex moved 0.06 R > 0.05 R
        record = analyse(donor=donor_with(signature=far_vertex))
        self.assertFalse(record["Qualifies"])
        self.assertTrue(any("claim vertex displacement" in r or "PortionCorrespondence" in r for r in record["Refusals"]),
                        record["Refusals"])
        far_context = shifted(exact, dy=0.06, context=[0])          # a context endpoint moved 0.06 R
        record = analyse(donor=donor_with(signature=far_context))
        self.assertFalse(record["Qualifies"])
        box = copy.deepcopy(exact)
        box["Box"] = [box["Box"][0] - 0.02, box["Box"][1], box["Box"][2], box["Box"][3]]   # a box face moved 0.02 R > 0.015 R
        record = analyse(donor=donor_with(signature=box))
        self.assertFalse(record["Qualifies"])
        self.assertTrue(any("box face displacement" in r for r in record["Refusals"]), record["Refusals"])

    def test_wide_lead_change_is_refused(self):
        # scale the whole exact signature by 15 % about the origin: every strip and gap scale together (W = +15 %), but
        # the vertices move far beyond 0.05 R and |W| > 11.3 % -> refused (recorded reasons)
        exact = FIXTURE["Exact"]["Signature"]
        scaled = copy.deepcopy(exact)
        for entry in scaled["Portions"] + scaled["Context"]:
            entry["P"] = [1.15 * v for v in entry["P"]]
        scaled["Box"] = [1.15 * v for v in scaled["Box"]]
        record = analyse(donor=donor_with(signature=scaled))
        self.assertFalse(record["Qualifies"])

    def test_lead_width_above_wmax_is_refused_directly(self):
        """Decision 438 (7): the |W| > WMax branch asserted on its own - the fixture pair (|W| 0.953 %) against a rule whose
        WMax is lowered below it; every other item of the pair is inside its limit, so the |W| refusal is the only one."""
        narrow = copy.deepcopy(rule())
        narrow["Policy"]["Domain"]["WMax"] = 0.005
        record = detection.analyse_pair(FIXTURE["Exact"]["Signature"], FIXTURE["Donor"], rule=narrow, radius=RADIUS)
        self.assertFalse(record["Qualifies"])
        self.assertEqual(len(record["Refusals"]), 1, record["Refusals"])
        self.assertTrue(record["Refusals"][0].startswith("|W| 0.00953 > 0.005"), record["Refusals"])
        self.assertAlmostEqual(record["W"], -0.009527170077628535, places=15)   # the metric is still recorded
        # the shipped WMax (11.3 %) admits the calibrated maximum (C2, 11.28 %) and refuses 11.31 %
        self.assertGreaterEqual(rule()["Policy"]["Domain"]["WMax"], 0.11280045501005607)
        self.assertLess(rule()["Policy"]["Domain"]["WMax"], 0.1131)

    def test_donor_status_refusals(self):
        for fields, needle in (({"QualificationStatus": "PendingQualification"}, "not Qualified"), ({"StatusProvisional": True},
                                                                                                    "StatusProvisional"),
                               ({"QualificationStatus": detection.STATUS_REUSED},
                                "reused model"), ({"ReusedFrom": {"Donor": "x"}}, "reused model"),
                               ({"TransplantedFrom": {"Donor": "x"}}, "reused model")):
            record = analyse(donor=donor_with(**fields))
            self.assertFalse(record["Qualifies"], fields)
            self.assertTrue(any(needle in r for r in record["Refusals"]), (fields, record["Refusals"]))

    def test_decision_431_override_rule(self):
        # the decision-292 / 311 cut-end kind: default admitted, SupersededBy 316
        record = analyse()
        self.assertTrue(record["DefaultAdmissible"])
        # B2's decision-353 CornerNeighborhoods override as a would-be donor: qualifies, fallback only
        donor = copy.deepcopy(FIXTURE["ExactOverrideBuilt"])
        donor["Signature"] = FIXTURE["Donor"]["Signature"]
        donor["StatusProvisional"] = False
        record = analyse(donor=donor)
        self.assertTrue(record["Qualifies"], record["Refusals"])
        self.assertFalse(record["DefaultAdmissible"])
        self.assertEqual(record["DonorBuildGateOverride"]["Kind"], "CornerNeighborhoods")
        self.assertIsNone(record["DonorBuildGateOverride"]["SupersededBy"])
        self.assertIn("fallback only", record["DefaultRefusal"])
        # a build-override mesh path without a record: fail closed for default
        unrecorded = donor_with(CouponMesh={"Path": "/x/build-override/case/identity.msh", "SHA256": "0" * 64})
        unrecorded.pop("BuildGateOverride")
        record = analyse(donor=unrecorded)
        self.assertTrue(record["Qualifies"])
        self.assertFalse(record["DefaultAdmissible"])
        self.assertEqual(record["DonorBuildGateOverride"]["Kind"], "UnrecordedBuildOverrideMeshPath")
        # another gate: refused for default
        other = donor_with()
        other["BuildGateOverride"]["Gate"] = "trace-diagonal-band"
        record = analyse(donor=other)
        self.assertFalse(record["DefaultAdmissible"])

    def test_vertex_list_difference_is_recorded_not_refused(self):
        donor = copy.deepcopy(FIXTURE["Donor"])
        donor["Signature"]["Vertices"] = donor["Signature"]["Vertices"][:-2]
        record = analyse(donor=donor)
        self.assertTrue(record["Qualifies"], record["Refusals"])
        self.assertTrue(record["VertexListDifference"]["Differs"])
        self.assertEqual(record["VertexListDifference"]["DonorVertices"], len(FIXTURE["Donor"]["Signature"]["Vertices"]) - 2)


def write_synthetic_basis(directory, *, box, nx=3, nz=2, metal_x=None, shift_knots=0.0):
    """A synthetic two-conductor trace basis on the four vertical faces of a box (the caps omitted:
    the transplant locates knots on the faces that carry them): nx knots along each in-plane edge, nz
    levels; conductor 1 occupies the x0 face's middle column, conductor 2 the x1 face's middle column
    (both at every z); every other vertex is a free knot. Returns (n_free, points)."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    x0, y0, x1, y1, z0, z1 = box
    xs = np.linspace(x0, x1, nx)
    ys = np.linspace(y0, y1, nx)
    zs = np.linspace(z0, z1, nz)
    vertices, index = [], {}

    def add(point, conductor):
        key = tuple(np.round(point, 12))
        if key in index:
            return index[key]
        index[key] = len(vertices)
        vertices.append((np.array(point), conductor))
        return index[key]

    triangles = []

    def face(points_2d, embed, conductor_of):
        grid = [[add(embed(u, v), conductor_of(u, v)) for v in zs] for u in points_2d]
        for i in range(len(points_2d) - 1):
            for k in range(nz - 1):
                a, b, c, d = grid[i][k], grid[i + 1][k], grid[i + 1][k + 1], grid[i][k + 1]
                triangles.append((a, b, c))
                triangles.append((a, c, d))
    middle = nx // 2
    face(ys, lambda y, z: [x0, y, z], lambda y, z: 1 if np.isclose(y, ys[middle]) else 0)
    face(ys, lambda y, z: [x1, y, z], lambda y, z: 2 if np.isclose(y, ys[middle]) else 0)
    face(xs, lambda x, z: [x, y0, z], lambda x, z: 0)
    face(xs, lambda x, z: [x, y1, z], lambda x, z: 0)
    points = np.array([p for p, _ in vertices])
    conductors = np.array([c for _, c in vertices])
    if shift_knots:
        # move the free knots of the y0 face's interior column (not on a box edge) along x
        for n, (p, c) in enumerate(vertices):
            if c == 0 and np.isclose(p[1], y0) and not (np.isclose(p[0], x0) or np.isclose(p[0], x1)):
                points[n, 0] += shift_knots
    basis = np.zeros(len(points), dtype=int)
    k = 0
    for n in range(len(points)):
        if conductors[n] == 0:
            k += 1
            basis[n] = k
    with open(directory / "trace-vertices.csv", "w") as out:
        out.write("vertex,x,y,z,basis,conductor\n")
        for n, p in enumerate(points):
            out.write(f"{n + 1},{p[0]:.16e},{p[1]:.16e},{p[2]:.16e},{basis[n]},{conductors[n]}\n")
    with open(directory / "trace-triangles.csv", "w") as out:
        out.write("triangle,vertex_i,vertex_j,vertex_k\n")
        for n, (a, b, c) in enumerate(triangles):
            out.write(f"{n + 1},{a + 1},{b + 1},{c + 1}\n")
    with open(directory / "basis-points.csv", "w") as out:
        out.write("x,y,z\n")
        for n in np.argsort(basis[basis > 0]):
            p = points[basis > 0][n]
            out.write(f"{p[0]:.16e},{p[1]:.16e},{p[2]:.16e}\n")
    return k, points


def write_synthetic_matrices(directory, size, seed=0):
    """SPD matrices in the library CSV forms: the fab / thin domain and surface (3 interfaces x 6 columns)."""
    rng = np.random.default_rng(seed)
    directory = Path(directory)

    def spd(scale):
        a = rng.standard_normal((size, size))
        return scale * (a @ a.T + size * np.eye(size))
    files = {}
    for role, _, name in transplant.MATRIX_ROLES[:4]:
        path = directory / name
        if role.endswith("domain"):
            transplant.write_domain(path, "  basis_i,  basis_j,                   Q_ij (J)", spd(1e-18))
        else:
            per = {k: {c: spd(1e-19 * k) for c in ("Q_ij (J)", "Q_ij normal (J)", "Q_ij tangential (J)")} for k in (1, 2, 3)}
            meta = {"columns": ["Q_ij (J)", "Q_ij normal (J)", "Q_ij tangential (J)"], "R (m)": "+1.900000000000e-06", "edge": 1,
                    "header": "interface,     edge,                      R (m),  basis_i,  basis_j,                   Q_ij (J),"
                              "            Q_ij normal (J),        Q_ij tangential (J)"}
            transplant.write_surface(path, meta, per)
        files[role] = str(path)
    return files


class SyntheticTransplant(unittest.TestCase):
    BOX = (-5.0, -5.0, 5.0, 5.0, -1.9, 2.0)
    INTERFACES = [{"Slot": 0, "Type": "MA", "Coupon": 1}, {"Slot": 0, "Type": "MS", "Coupon": 2}, {"Slot": 0, "Type": "SA", "Coupon": 3}]

    def bases(self, tmp, donor_box=None, shift_knots=0.0, nx_donor=3):
        exact_dir, donor_dir = Path(tmp) / "exact", Path(tmp) / "donor"
        n_e, _ = write_synthetic_basis(exact_dir, box=self.BOX)
        n_d, _ = write_synthetic_basis(donor_dir, box=donor_box or self.BOX, shift_knots=shift_knots, nx=nx_donor)
        exact = transplant.TraceBasis(exact_dir, "exact", interfaces=self.INTERFACES)
        donor = transplant.TraceBasis(donor_dir, "donor", interfaces=self.INTERFACES)
        donor.load_matrices(write_synthetic_matrices(donor_dir, donor.size))
        return exact, donor

    def test_identical_bases_give_the_identity_map_and_identical_matrices(self):
        with tempfile.TemporaryDirectory() as tmp:
            exact, donor = self.bases(tmp)
            result = transplant.transplant(exact, donor, RADIUS)
            self.assertTrue(np.array_equal(result["P"], np.eye(donor.size)))
            self.assertTrue(result["GatesPassed"] and result["T2Passed"])
            self.assertEqual(result["Tests"]["T1"]["MaxError"], 0.0)
            self.assertEqual(result["Tests"]["T4"]["WorstRelative"], 0.0)
            self.assertEqual(result["Tests"]["Knots"]["Matched"], exact.n_free)
            self.assertEqual(result["Tests"]["Knots"]["DonorOrphans"], [])
            for role in ("fab_domain", "thin_domain"):
                self.assertTrue(np.array_equal(result["Reused"][role]["Q_ij (J)"], donor.matrices[role]["Q_ij (J)"]))
            written, map_record = transplant.write_matrices(Path(tmp) / "out", donor, result)
            self.assertEqual(sorted(written), ["fab_domain", "fab_surface", "thin_domain", "thin_surface"])
            for role, value in written.items():
                self.assertLess(value["ReadBackMaxRelativeError"], 1e-12)
                self.assertEqual(Path(value["Path"]).read_bytes(), Path(donor.matrix_files[role]).read_bytes(), role)
            self.assertTrue(Path(map_record["Path"]).is_file())

    def test_affine_box_change_snaps_the_shared_knots(self):
        with tempfile.TemporaryDirectory() as tmp:
            # the donor box 2 nm wider per face: every knot maps onto an exact knot (snap), P = identity
            exact, donor = self.bases(tmp, donor_box=(-5.002, -5.002, 5.002, 5.002, -1.9, 2.0))
            result = transplant.transplant(exact, donor, RADIUS)
            self.assertTrue(np.array_equal(result["P"], np.eye(donor.size)))
            self.assertTrue(all(r["Snapped"] for r in result["Map"]["Rows"]))
            self.assertTrue(result["GatesPassed"])

    def test_displaced_knots_interpolate_rows_sum_to_one_and_the_knot_gate_reads_the_shift(self):
        with tempfile.TemporaryDirectory() as tmp:
            exact, donor = self.bases(tmp, shift_knots=0.05)      # 50 nm = 0.026 R along the y0 face
            result = transplant.transplant(exact, donor, RADIUS)
            P = result["P"]
            self.assertTrue(np.allclose(P[:donor.n_free].sum(axis=1), 1.0))   # every free row partitions unity (incl. the state column)
            self.assertTrue(result["Tests"]["T1"]["Passed"])
            self.assertTrue(result["Tests"]["T4"]["Passed"])
            self.assertTrue(result["Tests"]["T5"]["Passed"])
            knots = result["Tests"]["Knots"]
            self.assertGreater(len(knots["Displaced"]), 0)
            self.assertAlmostEqual(knots["MaxMatchedDisplacementOverR"], 0.05 / RADIUS, places=9)
            self.assertTrue(knots["Passed"])
            self.assertFalse(all(r["Exact"] for r in result["Map"]["Rows"]))
            too_far = transplant.transplant(exact, donor, RADIUS, domain_limits={"OrphanKnotsMax": 2, "KnotShiftMaxOverR": 0.02})
            self.assertFalse(too_far["Tests"]["Knots"]["Passed"])
            self.assertFalse(too_far["GatesPassed"])

    def test_orphan_knots_are_counted_in_both_directions(self):
        with tempfile.TemporaryDirectory() as tmp:
            exact, donor = self.bases(tmp, nx_donor=5)   # the donor carries more knots per edge: donor orphans
            result = transplant.transplant(exact, donor, RADIUS)
            knots = result["Tests"]["Knots"]
            self.assertGreater(len(knots["DonorOrphans"]), 2)
            self.assertFalse(knots["Passed"])
            self.assertFalse(result["GatesPassed"])
            # the extra donor knots on the conductor faces interpolate from the conductor vertices: T1 reads it and names the rows
            self.assertFalse(result["Tests"]["T1"]["Passed"])
            self.assertGreater(len(result["Tests"]["T1"]["RowsTouchingConductor1"]), 0)

    def test_missing_shelled_matrix_fails_closed_outside_measurement_only(self):
        with tempfile.TemporaryDirectory() as tmp:
            exact, donor = self.bases(tmp)
            library = {"Name": "synthetic", "Version": 3, "MatchingRadius": RADIUS,
                       "Models": [{"Name": "donor", "Topology": "SpatialEdgeCluster",
                                   "FabricatedMatrix": "donor/fabricated-domain-response-matrix.csv",
                                   "FabricatedSurfaceMatrix": "donor/fabricated-surface-response-matrix.csv",
                                   "ThinMatrix": "donor/thin-domain-response-matrix.csv",
                                   "ThinSurfaceMatrix": "donor/thin-surface-response-matrix.csv", "FabricatedSurfaceMatrixShelled": None,
                                   "Interfaces": self.INTERFACES}]}
            path = Path(tmp) / "process-library.json"
            path.write_text(json.dumps(library))
            with self.assertRaises(reuse.NearKeyReuseError):
                reuse.load_donor(library["Models"][0], Path(tmp), require_shelled=True)
            donor_loaded, paths, missing = reuse.load_donor(library["Models"][0], Path(tmp), require_shelled=False)
            self.assertTrue(missing)
            self.assertEqual(donor_loaded.size, donor.size)


class ReuseRecords(unittest.TestCase):
    def test_fallback_mode_needs_stop_record_and_approval_and_default_mode_needs_activation(self):
        with tempfile.TemporaryDirectory() as tmp:
            library = {"Name": "x", "Version": 3, "MatchingRadius": RADIUS, "Models": []}
            library_path = Path(tmp) / "process-library.json"
            library_path.write_text(json.dumps(library))
            common = {"exact_signature": FIXTURE["Exact"]["Signature"], "exact_basis_dir": tmp, "exact_model_entry": FIXTURE["Exact"],
                      "library": library, "library_path": library_path, "rule": rule(), "output": Path(tmp) / "out"}
            with self.assertRaises(reuse.NearKeyReuseError):
                reuse.reuse_requirement(mode="default", **common)
            with self.assertRaises(reuse.NearKeyReuseError):
                reuse.reuse_requirement(mode="fallback", approval="x", **common)
            with self.assertRaises(reuse.NearKeyReuseError):
                reuse.reuse_requirement(mode="fallback", stop_record={"Path": "p", "SHA256": "s"}, **common)
            with self.assertRaises(reuse.NearKeyReuseError):
                reuse.reuse_requirement(mode="off", **common)
            # decision 438 (5): default mode refuses a rule-file override even when that file activates the default
            override = Path(tmp) / "rule.json"
            override.write_text(json.dumps({k: v for k, v in activated(rule()).items() if not k.startswith("_")}))
            with self.assertRaises(reuse.NearKeyReuseError) as context:
                reuse.reuse_requirement(mode="default", **{**common, "rule": predictor.load_rule(override)})
            self.assertIn("refuses a rule-file override", str(context.exception))

    def test_stop_record_class_and_key_fields_are_checked(self):
        """Decision 438 (4): the STOP record's Status / StoppedBy class and its binding to the exact key."""
        name, case = "spatialedgecluster_edgecount-19_5ae3dbdd2c3d", "spatial-19-edge-5163ddf143c2"
        good = {"Case": case, "Status": "failed", "StoppedBy": {"Kind": "Stage", "Id": "mesh", "Stage": "mesh", "Message": "family-6"}}
        with tempfile.TemporaryDirectory() as tmp:
            stop = Path(tmp) / "stop.json"
            with self.assertRaises(reuse.NearKeyReuseError):
                reuse.stop_record_from_path(Path(tmp) / "missing.json", name)
            for bad in ("not json", json.dumps([good]), json.dumps({**good, "Status": "registered"}),
                        json.dumps({**good, "Status": "built", "Passed": True}), json.dumps({**good, "StoppedBy": "family-6"}),
                        json.dumps({**good, "StoppedBy": {"Kind": "Weather", "Id": "x"}}), json.dumps({**good,
                                                                                                       "StoppedBy": {"Kind": "Stage"}}),
                        json.dumps({**good, "Case": "other-case"})):
                stop.write_text(bad)
                with self.assertRaises(reuse.NearKeyReuseError, msg=bad):
                    reuse.stop_record_from_path(stop, name, case)
            stop.write_text(json.dumps({**good, "Case": "spatial-19-edge-000000000000"}))   # another case, hash12 absent
            with self.assertRaises(reuse.NearKeyReuseError):
                reuse.stop_record_from_path(stop, name)
            stop.write_text(json.dumps(good))
            record = reuse.stop_record_from_path(stop, name, case)
            self.assertEqual((record["Status"], record["StoppedBy"]["Kind"], record["Case"], len(record["SHA256"])), ("failed",
                                                                                                                      "Stage", case, 64))
            # a build case record: Passed false + StoppedBy of the build matrix
            stop.write_text(json.dumps({"Case": case, "Status": "failed", "Passed": False,
                                        "StoppedBy": {"Kind": "HeadroomGate", "Id": "MaximumElements", "Stage": "headroom"}}))
            self.assertEqual(reuse.stop_record_from_path(stop, name, case)["StoppedBy"]["Id"], "MaximumElements")
            # without a case id the record must name the requirement's hash12 in Case / Requirement / Key
            stop.write_text(json.dumps({**good, "Case": "x", "Requirement": name}))
            self.assertEqual(reuse.stop_record_from_path(stop, name)["Case"], "x")

    def test_requirement_key_must_be_the_signature_hash(self):
        """Decision 438 (3): signature_hash(requirement Signature) == the requirement key, asserted before any transplant."""
        with tempfile.TemporaryDirectory() as tmp:
            library = {"Name": "x", "Version": 3, "MatchingRadius": RADIUS, "Models": [FIXTURE["Donor"]]}
            library_path = Path(tmp) / "process-library.json"
            library_path.write_text(json.dumps(library))
            common = {"exact_signature": FIXTURE["Exact"]["Signature"], "exact_basis_dir": Path(tmp) / "no-basis",
                      "exact_model_entry": FIXTURE["Exact"], "library": library, "library_path": library_path, "rule": rule(),
                      "mode": "measurement-only", "output": Path(tmp) / "out", "log": lambda *_: None}
            with self.assertRaises(reuse.NearKeyReuseError) as context:
                reuse.reuse_requirement(requirement_key="0123456789ab", **common)
            self.assertIn("!= the requirement key", str(context.exception))
            # the right key (the full hash or its >= 12-hex prefix) passes the assertion and fails later on the absent basis
            key = detection.signature_library.signature_hash(FIXTURE["Exact"]["Signature"])
            for given in (key, key[:12]):
                with self.assertRaises((transplant.TransplantError, OSError)):
                    reuse.reuse_requirement(requirement_key=given, **common)

    def test_donor_stored_record_required_and_sha_checked(self):
        """Decision 438 (2): in Default / Fallback the donor's stored (F) record is required and must re-hash to the
        library's SpatialQualification.RecordSHA256."""
        with tempfile.TemporaryDirectory() as tmp:
            record = Path(tmp) / "spatial-qualification.json"
            record.write_text(json.dumps({"Traces": [{"Name": "state-2", "MatrixIdentity": {"SA": {"Predicted": 1.0},
                                                                                            "Domain": {"Predicted": 2.0}}}]}))
            sha = transplant.sha256(record)
            donor = {"Name": "donor", "SpatialQualification": {"Record": "/cluster/root/spatial-qualification.json", "RecordSHA256": sha}}
            with self.assertRaises(reuse.NearKeyReuseError):              # unreadable, required
                reuse.donor_stored_state2(donor, required=True)
            self.assertEqual(reuse.donor_stored_state2(donor, required=False), (None, None))   # measurement-only: recorded, not refused
            predicted, actual = reuse.donor_stored_state2(donor, record_path=record, required=True)
            self.assertEqual((predicted["SA"], actual), (1.0, sha))
            mapped, _ = reuse.donor_stored_state2(donor, record_roots={"/cluster/root/": str(Path(tmp)) + "/"}, required=True)
            self.assertEqual(mapped["Domain"], 2.0)
            with self.assertRaises(reuse.NearKeyReuseError):              # sha mismatch (in every mode)
                reuse.donor_stored_state2({**donor, "SpatialQualification": {"Record": None, "RecordSHA256": "0" * 64}}, record_path=record)
            with self.assertRaises(reuse.NearKeyReuseError):              # no recorded sha to check against
                reuse.donor_stored_state2({"Name": "donor"}, record_path=record)
            record.write_text(json.dumps({"Traces": []}))
            donor["SpatialQualification"]["RecordSHA256"] = transplant.sha256(record)
            with self.assertRaises(reuse.NearKeyReuseError):              # no state-2 identity in the record
                reuse.donor_stored_state2(donor, record_path=record)

    def test_donor_tail_and_structure_key_refusal_record(self):
        self.assertAlmostEqual(reuse.donor_tail(FIXTURE["Donor"]), 0.022262609816823264, places=12)
        with self.assertRaises(reuse.NearKeyReuseError):
            reuse.donor_tail({"Name": "no-ma"})
        with tempfile.TemporaryDirectory() as tmp:
            library = {"Name": "x", "Version": 3, "MatchingRadius": RADIUS, "Models": [FIXTURE["Donor"]]}
            library_path = Path(tmp) / "process-library.json"
            library_path.write_text(json.dumps(library))
            exact = copy.deepcopy(FIXTURE["Transmon4Edge"])
            # a requirement outside the admissible keys: refused before any basis is read (no basis directory needed)
            record = reuse.reuse_requirement(exact_signature=exact["Signature"], exact_basis_dir=Path(tmp) / "no-basis",
                                             exact_model_entry=exact,
                                             library=library, library_path=library_path, rule=rule(), mode="measurement-only",
                                             output=Path(tmp) / "out", log=lambda *_: None)
            self.assertFalse(record["Reused"])
            self.assertTrue(record["Refused"]["Reason"].startswith(detection.REFUSAL_STRUCTURE_KEY))
            self.assertTrue((Path(tmp) / "out" / reuse.REFUSAL_RECORD).is_file())

    @staticmethod
    def reused_model(the_rule, name="spatialedgecluster_edgecount-19_5ae3dbdd2c3d-reused", mode="Fallback", **fields):
        """A complete reused model record (the fields check_reused_model reads) for the assembling rule."""
        model = {"Name": name, "QualificationStatus": detection.STATUS_REUSED, "ReuseMode": mode,
                 "PredictedReuseError": {"RuleVersion": the_rule["RuleVersion"], "RuleFileSHA256": the_rule["_sha256"]},
                 "ReusedFrom": {"Donor": "donor"}, "NearKey": {"W": -0.01}, "TransplantTests": {"GatesPassed": True}}
        if mode == "Fallback":
            model.update({"FallbackStopRecord": {"Path": "stop.json", "SHA256": "a" * 64}, "Approval": "supervisor decision NNN"})
        if mode == "Default":
            model["DefaultActivation"] = the_rule["Policy"]["DefaultActivation"]
        model.update(fields)
        return model

    def test_assemble_adds_models_and_the_header(self):
        library = {"Name": "base", "Version": 3, "MatchingRadius": RADIUS, "Note": "base note", "Models": [copy.deepcopy(FIXTURE["Donor"])]}
        the_rule = rule()
        out = reuse.assemble_library(library, [self.reused_model(the_rule)], rule=the_rule, name="base+reuse")
        self.assertEqual(len(out["Models"]), 2)
        header = out["NearKeyReuse"]
        self.assertEqual(header["RuleVersion"], "nearkey-reuse-rule-v1")
        self.assertEqual(len(header["AdmissibleStructureKeys"]), 4)
        self.assertEqual(header["Policy"]["DefaultBoundPct"], {"SA": 0.5, "MS": 0.5, "MA_sharp": 1.0})
        self.assertFalse(header["Policy"]["DefaultActive"])
        self.assertEqual(header["Counts"]["Fallback"], 1)
        self.assertNotIn("MeasurementOnly", out)
        active = activated(the_rule)
        default = reuse.assemble_library(library, [self.reused_model(active, mode="Default")], rule=active, name="base+default")
        self.assertEqual(default["NearKeyReuse"]["Counts"]["Default"], 1)
        self.assertTrue(default["NearKeyReuse"]["Policy"]["DefaultActive"])
        same_name = self.reused_model(the_rule, name=FIXTURE["Donor"]["Name"], mode="MeasurementOnly")
        with self.assertRaises(reuse.NearKeyReuseError):
            reuse.assemble_library(library, [same_name], rule=the_rule, name="clash")
        demo = reuse.assemble_library(library, [same_name], rule=the_rule, name="demo", measurement_only=True)
        self.assertTrue(demo["MeasurementOnly"])
        self.assertEqual(demo["NearKeyReuse"]["Replaced"], 1)
        self.assertIn("NEVER a library of record", demo["Note"])

    def test_assemble_is_fail_closed_per_model(self):
        """Decision 438 (1) / review MAJOR-1: the four reproduced acceptances are refusals, plus every listed check."""
        library = {"Name": "base", "Version": 3, "MatchingRadius": RADIUS, "Models": []}
        the_rule = rule()
        active = activated(the_rule)

        def refused(model, assembling_rule=the_rule, needle=None, **kwargs):
            with self.assertRaises(reuse.NearKeyReuseError) as context:
                reuse.assemble_library(library, [model], rule=assembling_rule, name="bad", **kwargs)
            if needle:
                self.assertIn(needle, str(context.exception))
        # review case 1: a MeasurementOnly model (exact Name kept) assembled WITHOUT --measurement-only
        refused(self.reused_model(the_rule, name="spatialedgecluster_edgecount-19_5ae3dbdd2c3d", mode="MeasurementOnly"),
                needle="MeasurementOnly model enters a measurement-only assembly only")
        # review case 2: a Default model while the rule's default is inactive
        refused(self.reused_model(the_rule, mode="Default", DefaultActivation={"ValidationPair": "pair-5", "RecordSHA256": "f" * 64}),
                needle="default is inactive")
        # review case 3: a Fallback model without its records
        refused({"Name": "spatialedgecluster_edgecount-19_5ae3dbdd2c3d-reused", "QualificationStatus": detection.STATUS_REUSED,
                 "ReuseMode": "Fallback"}, needle="FallbackStopRecord")
        refused(self.reused_model(the_rule, FallbackStopRecord={"Path": "stop.json"}), needle="FallbackStopRecord")
        refused(self.reused_model(the_rule, Approval="   "), needle="Approval")
        # review case 4: the model's RuleVersion / RuleFileSHA256 vs the assembling rule's
        refused(self.reused_model(the_rule, PredictedReuseError={"RuleVersion": "nearkey-reuse-rule-v2",
                                                                 "RuleFileSHA256": the_rule["_sha256"]}),
                needle="RuleVersion / RuleFileSHA256")
        refused(self.reused_model(the_rule, PredictedReuseError={"RuleVersion": the_rule["RuleVersion"], "RuleFileSHA256": "0" * 64}),
                needle="RuleVersion / RuleFileSHA256")
        # Default: the model's DefaultActivation must equal the rule's
        refused(self.reused_model(active, mode="Default", DefaultActivation={**active["Policy"]["DefaultActivation"],
                                                                             "RecordSHA256": "0" * 64}),
                assembling_rule=active, needle="DefaultActivation")
        # the remaining checks
        refused({**self.reused_model(the_rule), "QualificationStatus": "Qualified"}, needle="not a ReusedResponse")
        refused(self.reused_model(the_rule, mode="Always"), needle="ReuseMode")
        for key in ("ReusedFrom", "NearKey", "TransplantTests", "PredictedReuseError"):
            model = self.reused_model(the_rule)
            del model[key]
            refused(model, needle=key)
        refused(self.reused_model(the_rule, TransplantTests={"GatesPassed": False}), needle="GatesPassed")
        refused(self.reused_model(the_rule, name="spatialedgecluster_edgecount-19_5ae3dbdd2c3d"), needle="-reused")
        # a complete Fallback model passes; the same model is accepted under --measurement-only too
        out = reuse.assemble_library(library, [self.reused_model(the_rule)], rule=the_rule, name="ok")
        self.assertEqual(out["NearKeyReuse"]["Counts"]["Fallback"], 1)


class BuildHook(unittest.TestCase):
    """device_coupons: the --nearkey-reuse options fail closed before any discovery; a reused coupon is not registered."""

    def test_option_validation_fails_closed_before_discovery(self):
        with tempfile.TemporaryDirectory() as tmp:
            common = {"palace": "/nonexistent/palace", "output": Path(tmp) / "out", "manifest_path": Path(tmp) / "manifest.json"}
            for kwargs in ({"nearkey_reuse_mode": "always"},
                           {"nearkey_reuse_mode": "default"},                                   # no DefaultActivation record
                           {"nearkey_reuse_mode": "fallback"},                                  # no approval / stop record
                           {"nearkey_reuse_mode": "fallback", "nearkey_fallback_approval": "x"},
                           {"nearkey_reuse_mode": "fallback", "nearkey_fallback_stop_records": ["5ae3=stop.json"]},
                           {"nearkey_reuse_mode": "fallback", "nearkey_fallback_approval": "x",
                            "nearkey_fallback_stop_records": ["zz=stop.json"]},
                           {"nearkey_reuse_mode": "fallback", "nearkey_fallback_approval": "x",
                            "nearkey_fallback_stop_records": ["5a=a.json", "5a=b.json"]},
                           {"nearkey_reuse_mode": "off", "nearkey_fallback_approval": "x"}):
                with self.assertRaises(device_coupons.DeviceAdapterError, msg=str(kwargs)):
                    device_coupons.prepare_device_sources(Path(tmp) / "device.json", **common, **kwargs)
            self.assertFalse((Path(tmp) / "out").exists())

    def test_reused_coupon_is_not_registered(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "manifest.json"
            manifest.write_text(json.dumps({"Cases": [], "Gates": {}}))
            record = {"Device": {"Config": "d.json", "SHA256": "0" * 64}, "Output": tmp,
                      "Coupons": [{"Case": "spatial-19-edge-000000000000", "Requirement": "r", "ContentHash": "0" * 64, "Directory": tmp,
                                   "NearKeyReuse": {"Reused": True, "Donor": "donor", "ReuseMode": "Fallback"}}]}
            out = device_coupons.register_device_sources(record, manifest_path=manifest, mesh_recipe="recipe", probe=None)
            coupon = out["Coupons"][0]
            self.assertEqual(coupon["Registration"]["Status"], device_coupons.STATUS_NEARKEY_REUSED)
            self.assertIsNone(coupon["ThinCase"])
            self.assertTrue((Path(tmp) / device_coupons.DEVICE_RECORD).is_file())

    def test_version_1_record_and_missing_stop_record_are_recorded_not_offered(self):
        rule_v1 = rule()
        with tempfile.TemporaryDirectory() as tmp:
            coupon = {"Id": "spatialedgecluster_edgecount-19_5ae3dbdd2c3d", "Topology": "SpatialEdgeCluster",
                      "Geometry": {"Edges": [1]}, "Hash": "5ae3"}
            result = device_coupons.nearkey_reuse_for_coupon(coupon, tmp, "case", library={}, library_path=tmp, rule=rule_v1,
                                                             mode="fallback",
                                                             output=Path(tmp) / "reused", stop_record_path=None, approval="x")
            self.assertFalse(result["Reused"])
            self.assertTrue(result["Reason"].startswith("NotApplicable"))
            coupon["Geometry"] = {"Signature": FIXTURE["Exact"]["Signature"]}
            result = device_coupons.nearkey_reuse_for_coupon(coupon, tmp, "case", library={}, library_path=tmp, rule=rule_v1,
                                                             mode="fallback",
                                                             output=Path(tmp) / "reused", stop_record_path=None, approval="x")
            self.assertTrue(result["Reason"].startswith("NoStopRecord"))


@unittest.skipUnless(LIBRARY_OF_RECORD.is_file() and MODELS_MIRROR.is_dir() and SENS_LIBRARY.is_dir(),
                     "the local evidence mirror is not mounted")
class SixPairsIdentity(unittest.TestCase):
    """The production transplant reproduces the sensitivity lane's transplanted matrices BYTE-IDENTICALLY (the four
    coupon matrices and interpolation-map-P.csv of s2-r1p9-v3-b-sens-c3) on the six measured pairs; the written model
    carries every record of DESIGN 4.3. The shelled matrix is transplanted when the donor's file is mirrored."""

    def test_six_pairs(self):
        library = json.loads(LIBRARY_OF_RECORD.read_text())
        by_name = {m["Name"]: m for m in library["Models"]}
        the_rule = rule()
        with tempfile.TemporaryDirectory() as tmp:
            for pair, (exact_dir, donor_name, w, s) in SIX_PAIRS.items():
                exact_hash = exact_dir.split("-")[-1]
                exact_entry = copy.deepcopy(next(m for m in library["Models"] if m["Name"].endswith(exact_hash)))
                shelled = SHELLED / donor_name.split("_")[-1] / "surface-response-matrix.csv"
                record = reuse.reuse_requirement(
                    exact_signature=exact_entry["Signature"], exact_basis_dir=MODELS_MIRROR / exact_dir, exact_model_entry=exact_entry,
                    library=library, library_path=LIBRARY_OF_RECORD, rule=the_rule, mode="measurement-only", output=Path(tmp) / pair,
                    donors=[donor_name], matrices_root=MODELS_MIRROR,
                    shelled_paths={donor_name: str(shelled)} if shelled.is_file() else None,
                    log=lambda *_: None)
                self.assertTrue(record["Reused"], (pair, record["Refused"]))
                model = record["Model"]
                self.assertEqual(model["QualificationStatus"], detection.STATUS_REUSED)
                self.assertEqual(model["ReuseMode"], "MeasurementOnly")
                self.assertAlmostEqual(model["NearKey"]["W"], w, places=15)
                self.assertAlmostEqual(model["NearKey"]["S"], s, places=15)
                self.assertEqual(model["ReusedFrom"]["Donor"], donor_name)
                for key in ("ReusedFrom", "NearKey", "PredictedReuseError", "TransplantTests", "TailFromDonor", "MA", "Note"):
                    self.assertIn(key, model, (pair, key))
                self.assertEqual(model["PredictedReuseError"]["RuleVersion"], "nearkey-reuse-rule-v1")
                self.assertTrue(model["PredictedReuseError"]["MA_sharp"]["TailFromDonor"])
                self.assertEqual(model["TransplantTests"]["T2"]["Passed"], pair != "B3")   # decision 380: B3's orphan knots
                self.assertTrue(model["TransplantTests"]["GatesPassed"], pair)
                model_dir = Path(record["ModelDirectory"])
                sens_dir = SENS_LIBRARY / "models" / f"{exact_dir}-sens"
                for name in ("fabricated-domain-response-matrix.csv", "fabricated-surface-response-matrix.csv",
                             "thin-domain-response-matrix.csv",
                             "thin-surface-response-matrix.csv", "interpolation-map-P.csv"):
                    self.assertEqual(transplant.sha256(model_dir / name), transplant.sha256(sens_dir / name), (pair, name))
                if shelled.is_file():
                    self.assertEqual(model["ShelledMatrix"]["Status"], "Transplanted")
                    self.assertTrue((model_dir / "fabricated-surface-response-matrix-shelled.csv").is_file())
                else:
                    self.assertTrue(model["ShelledMatrix"]["Status"].startswith("NotTransplanted"))
                # the pair's information-only policy verdicts (Option A, as the rule stands: default inactive)
                decisions = model["PredictedReuseError"]["PolicyDecisions"]
                self.assertFalse(decisions["default"]["Allowed"])
                self.assertTrue(decisions["FallbackIfStopRecorded"]["Allowed"], pair)
                for name in model_dir.iterdir():
                    name.unlink()

    def test_build_hook_fallback_reuse_of_b3(self):
        """B3 (T2 failed, SA bound 0.568 > 0.5): fallback-only; through the device_coupons hook with a STOP record + approval."""
        library = json.loads(LIBRARY_OF_RECORD.read_text())
        exact_dir, donor_name, _, _ = SIX_PAIRS["B3"]
        exact_entry = copy.deepcopy(next(m for m in library["Models"] if m["Name"].endswith("1de0718f1fc8")))
        shelled = SHELLED / donor_name.split("_")[-1] / "surface-response-matrix.csv"
        if not shelled.is_file():
            self.skipTest("the donor's shelled matrix is not mirrored")
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp) / "work"
            work.mkdir()
            for name in ("trace-vertices.csv", "trace-triangles.csv", "basis-points.csv"):
                (work / name).write_bytes((MODELS_MIRROR / exact_dir / name).read_bytes())
            generated = {key: exact_entry[key] for key in ("Name", "Topology", "Signature", "Interfaces", "SupportBox", "Edges",
                                                           "ContextEdges")}
            (work / "process-library.json").write_text(json.dumps({"Models": [generated]}))
            stop = Path(tmp) / "stop.json"
            stop.write_text(json.dumps({"Case": "spatial-19-edge-74b842a93a96", "Status": "failed",
                                        "StoppedBy": {"Kind": "Stage", "Id": "mesh", "Stage": "mesh",
                                                      "Message": "family-6 (synthetic test record)"}}))
            # a local library whose donor paths resolve: the models' relative paths against the mirror root
            local = copy.deepcopy(library)
            for model in local["Models"]:
                for key in ("FabricatedMatrix", "FabricatedSurfaceMatrix", "ThinMatrix", "ThinSurfaceMatrix", "BasisPoints"):
                    if key in model:
                        model[key] = str(MODELS_MIRROR / Path(model[key]).parent.name / Path(model[key]).name)
                if model["Name"] == donor_name:
                    model["FabricatedSurfaceMatrixShelled"] = str(shelled)
                    local_donor = model
            library_path = Path(tmp) / "process-library.json"
            library_path.write_text(json.dumps(local))
            coupon = {"Id": exact_entry["Name"], "Topology": "SpatialEdgeCluster", "Geometry": {"Signature": exact_entry["Signature"]},
                      "Hash": "5d3d8ac62dd1", "Interfaces": exact_entry["Interfaces"]}
            roots = {"/data/home/simlap/bedrock-tests/": "/Users/simlap/bedrock-tests/"}   # the donors' (F) records of record, mirrored
            hook = {"library": local, "library_path": library_path, "rule": rule(), "mode": "fallback", "output": Path(tmp) / "reused",
                    "stop_record_path": str(stop), "approval": "supervisor decision NNN (test)", "log": lambda *_: None}
            # decision 438 (2): without the donors' stored (F) records every candidate is NotLoadable -> refused, recorded
            refused = device_coupons.nearkey_reuse_for_coupon(coupon, work, "spatial-19-edge-74b842a93a96", **hook)
            self.assertFalse(refused["Reused"])
            self.assertIn("rejected as NotLoadable", refused["Reason"])
            refusal = json.loads(Path(refused["Record"]).read_text())
            donor_entry = next(c for c in refusal["Candidates"] if c["Model"] == donor_name)
            self.assertTrue(donor_entry["RejectedReason"].startswith("NotLoadable") and "(F) record" in donor_entry["RejectedReason"],
                            donor_entry)
            # decision 438 (3): a requirement key that is not the signature hash stops the reuse
            with self.assertRaises(device_coupons.DeviceAdapterError):
                device_coupons.nearkey_reuse_for_coupon({**coupon, "Hash": "000000000000"}, work, "spatial-19-edge-74b842a93a96",
                                                        record_roots=roots, **hook)
            # decision 438 (4): a STOP record of another case is refused
            other = Path(tmp) / "other-stop.json"
            other.write_text(json.dumps({"Case": "spatial-19-edge-000000000000", "Status": "failed",
                                         "StoppedBy": {"Kind": "Stage", "Id": "mesh", "Stage": "mesh"}}))
            with self.assertRaises(device_coupons.DeviceAdapterError):
                device_coupons.nearkey_reuse_for_coupon(coupon, work, "spatial-19-edge-74b842a93a96", record_roots=roots,
                                                        **{**hook, "stop_record_path": str(other)})
            result = device_coupons.nearkey_reuse_for_coupon(coupon, work, "spatial-19-edge-74b842a93a96", record_roots=roots, **hook)
            self.assertTrue(result["Reused"], result)
            self.assertEqual(result["ReuseMode"], "Fallback")
            self.assertEqual(result["Donor"], donor_name)
            model = json.loads((Path(result["ModelDirectory"]) / reuse.REUSED_MODEL_RECORD).read_text())
            self.assertEqual(model["Name"], exact_entry["Name"] + "-reused")
            self.assertEqual(model["FallbackStopRecord"]["Path"], str(stop))
            self.assertEqual((model["FallbackStopRecord"]["Status"], model["FallbackStopRecord"]["StoppedBy"]["Kind"]), ("failed", "Stage"))
            self.assertEqual(model["Approval"], "supervisor decision NNN (test)")
            self.assertTrue(model["TransplantTests"]["T4"]["StoredRecordUsed"])
            self.assertEqual(model["ReusedFrom"]["DonorRecordSHA256"], local_donor["SpatialQualification"]["RecordSHA256"])
            key_check = model["ReusedFrom"]["GeneratedSignatureHashEqualsKey"]
            self.assertTrue(key_check["Equal"])
            self.assertTrue(key_check["SignatureHash"].startswith("5d3d8ac62dd1"))
            self.assertEqual(key_check["SignatureHash"], key_check["GeneratedSignatureHash"])
            # the written Fallback model passes the library-of-record assembly checks (decision 438 (1))
            assembled = reuse.assemble_library({"Name": "base", "Version": 3, "MatchingRadius": 1.9, "Models": []}, [model], rule=rule(),
                                               name="base+b3")
            self.assertEqual(assembled["NearKeyReuse"]["Counts"]["Fallback"], 1)
            self.assertFalse(model["TransplantTests"]["T2"]["Passed"])
            self.assertTrue(model["TransplantTests"]["GatesPassed"])
            self.assertEqual(model["ShelledMatrix"]["Status"], "Transplanted")
            self.assertTrue(model["PredictedReuseError"]["PolicyDecisions"]["fallback"]["Allowed"])
            self.assertFalse(model["PredictedReuseError"]["PolicyDecisions"]["default"]["Allowed"])
            record = json.loads(Path(result["Record"]).read_text())
            self.assertEqual([c["Model"] for c in record["Candidates"] if c["Chosen"]], [donor_name])
            # the other candidate of the family (B2's exact model, W +2.44 %) is recorded with its bound, not chosen
            self.assertTrue(any(c["Model"].endswith("5ae3dbdd2c3d") for c in record["Candidates"]))


if __name__ == "__main__":
    unittest.main()
