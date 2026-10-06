# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""`coupon-library build --device` on the transmon example (the device the gallery cases
came from): discovery (version-2 identification manifest) -> planner ->
cluster_signature_geometry (the v2 cluster contract: the coupon in the canonical frame of
the record's Signature) -> generate_spatial_response.py --basis-only -> content-hashed
source directories -> register_case (idempotent) -> preflight.  The other families are
recorded out of scope.  Needs the Palace executable (build/bin/palace) for the geometry
preflights."""
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "qualify"))
sys.path.insert(0, str(HERE.parents[1] / "cpw2d"))
import numpy as np  # noqa: E402

import case_inputs  # noqa: E402
import coupon_library  # noqa: E402
import device_coupons  # noqa: E402
import prepare_surface_response_coupons as planner  # noqa: E402
import register_case  # noqa: E402
import trace_basis  # noqa: E402
from derive_semantic_contract import derive as derive_contract  # noqa: E402

REPOSITORY = HERE.parents[2]
PALACE = Path(os.environ.get("PALACE_EXECUTABLE", REPOSITORY / "build" / "bin" / "palace"))
TRANSMON = REPOSITORY / "examples" / "transmon"
DEVICE_CONFIG = TRANSMON / "transmon_surface_coarse.json"
PROCESS_SEED = TRANSMON / "benchmark" / "transmon_surface_process_seed.json"


def mask_areas(path):
    """Facet area per (conductor, plane) of a plan-view-mask.csv (shoelace, facets disjoint)."""
    facets = {}
    with open(path, newline="") as stream:
        for row in csv.DictReader(stream):
            facets.setdefault((int(row["Facet"]), int(row["Conductor"]), float(row["Plane"])), []).append(
                (float(row["X"]), float(row["Y"])))
    areas = {}
    for (_, conductor, plane), points in facets.items():
        twice = sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(points, points[1:] + points[:1]))
        areas[(conductor, plane)] = areas.get((conductor, plane), 0.0) + abs(twice) / 2.0
    return areas


def available():
    return PALACE.is_file() and DEVICE_CONFIG.is_file() and PROCESS_SEED.is_file() and shutil.which("git") is not None


def device_config(tmp):
    """The transmon preflight config (examples/transmon/prepare_surface_response_preflight.py)
    against a process seed whose interface layers equal the config's Dielectric entries -
    the checked-in seed's MS permittivity (11.45) differs from the config's (11.47) and
    Palace refuses the mismatch (ValidateLibraryInterfaceLayers), so the test binds the
    config's values into a copy of the seed (recorded here, not a change of either fixture)."""
    config = json.loads(DEVICE_CONFIG.read_text())
    seed = json.loads(PROCESS_SEED.read_text())
    for entry in config["Boundaries"]["Postprocessing"]["Dielectric"]:
        seed["Fabrication"]["InterfaceLayers"][entry["Type"]] = {"Thickness": entry["Thickness"],
                                                                "Permittivity": entry["Permittivity"]}
    seed_path = tmp / "process-seed.json"
    seed_path.write_text(json.dumps(seed, indent=2) + "\n")
    mesh = Path(config["Model"]["Mesh"])
    config["Model"]["Mesh"] = str(mesh if mesh.is_absolute() else (DEVICE_CONFIG.parent / mesh).resolve())
    config["Problem"]["Output"] = str(tmp / "postpro")
    config["Solver"]["SurfaceResponseCorrection"] = {"Library": str(seed_path), "TargetInterfaces": [1, 2, 3],
                                                     "UnmatchedPolicy": "Warn"}
    path = tmp / "device.json"
    path.write_text(json.dumps(config, indent=2) + "\n")
    return path


def census_probe(labels):
    """A probe standing in for the stages-only build (test_register_case.census_probe)."""
    def probe(probe_manifest, case_id, root, log):
        root.mkdir(parents=True)
        Path(log).write_text("stub probe\n")
        (root / "build-census.json").write_text(json.dumps(
            {"InterfaceAreas": [{"Attribute": label, "Area": 1.0} for label in labels]}))
        return {"Case": case_id, "Commit": "stub", "Root": str(root), "Status": "built", "Stage": "labels-only",
                "ReturnCode": 0, "ScopeGuard": None, "Message": None}
    return probe


class DeviceBasisDefaultTest(unittest.TestCase):
    def test_delaunay_caps_are_the_device_default_and_ear_clipping_an_explicit_option(self):
        """Decision 57: `build --device` triangulates the device basis's box caps with Delaunay
        flips unless `--cap-triangulation ear-clipping` (the gallery producer's) is given."""
        self.assertEqual(device_coupons.DEFAULT_CAP_TRIANGULATION, "delaunay")
        parser = coupon_library.build_parser()
        common = ["build", "--device", "device.json", "--palace", "palace", "--root", "root"]
        self.assertEqual(parser.parse_args(common).cap_triangulation, "delaunay")
        self.assertEqual(parser.parse_args(common + ["--cap-triangulation", "ear-clipping"]).cap_triangulation,
                         "ear-clipping")
        with self.assertRaises(SystemExit):
            parser.parse_args(common + ["--cap-triangulation", "fan"])

    def test_omit_requirement_is_a_repeatable_device_option(self):
        """`build --device --omit-requirement HASH_PREFIX` (repeatable) reaches the discovery; it
        is refused without --device (it only names discovery placeholders)."""
        parser = coupon_library.build_parser()
        common = ["build", "--device", "device.json", "--palace", "palace", "--root", "root"]
        self.assertEqual(parser.parse_args(common).omit_requirement, [])
        self.assertEqual(parser.parse_args(common + ["--omit-requirement", "69ca648cc16f", "--omit-requirement", "abcd"])
                         .omit_requirement, ["69ca648cc16f", "abcd"])
        with self.assertRaises(SystemExit):
            coupon_library.main(["build", "--manifest", "m.json", "--omit-requirement", "69ca648cc16f"])

    def test_per_requirement_span_and_element_caps_are_device_options_with_recorded_reasons(self):
        """`build --device --support-span-cap HASH_PREFIX=OVER_R --support-span-cap-reason TEXT` and
        `--element-cap HASH_PREFIX=ELEMENTS --element-cap-approval TEXT --element-cap-reason TEXT`
        (supervisor decision of 2026-10-02 for the S1p loop end: one closed 38-edge feature of plan
        span 18.24 R > the generator's 16R bound, ~12-13 M elements > the 6 M gate) name one spatial
        coupon each; they are refused without --device, without their texts, with a non-hex prefix,
        a non-positive value, a repeated prefix, or a prefix no spatial coupon carries."""
        parser = coupon_library.build_parser()
        common = ["build", "--device", "device.json", "--palace", "palace", "--root", "root"]
        args = parser.parse_args(common + ["--support-span-cap", "5ed91f8890c0=20", "--support-span-cap-reason", "closed loop",
                                           "--element-cap", "5ed91f8890c0=16000000", "--element-cap-approval", "supervisor",
                                           "--element-cap-reason", "12-13 M elements"])
        self.assertEqual((args.support_span_cap, args.element_cap), (["5ed91f8890c0=20"], ["5ed91f8890c0=16000000"]))
        self.assertEqual(device_coupons.requirement_option_kwargs(args),
                         {"support_span_caps": ["5ed91f8890c0=20"], "support_span_cap_reason": "closed loop",
                          "element_caps": ["5ed91f8890c0=16000000"], "element_cap_approval": "supervisor",
                          "element_cap_reason": "12-13 M elements", "nearkey_reuse_mode": "off", "nearkey_fallback_approval": None,
                          "nearkey_fallback_stop_records": [], "nearkey_rule": None, "nearkey_record_roots": []})
        self.assertEqual(parser.parse_args(common).support_span_cap, [])
        with self.assertRaises(SystemExit):
            coupon_library.main(["build", "--manifest", "m.json", "--element-cap", "5ed91f8890c0=16000000"])
        # the near-key reuse options (decision 428) are --device options too: refused without --device
        with self.assertRaises(SystemExit):
            coupon_library.main(["build", "--manifest", "m.json", "--nearkey-reuse", "fallback"])
        with self.assertRaises(SystemExit):
            coupon_library.main(["build", "--manifest", "m.json", "--nearkey-fallback-approval", "x"])
        nearkey = parser.parse_args(common + ["--nearkey-reuse", "fallback", "--nearkey-fallback-approval", "decision NNN",
                                              "--nearkey-fallback-stop-record", "5d3d8ac62dd1=stop.json"])
        self.assertEqual(device_coupons.requirement_option_kwargs(nearkey)["nearkey_reuse_mode"], "fallback")
        self.assertEqual(device_coupons.requirement_option_kwargs(nearkey)["nearkey_fallback_stop_records"], ["5d3d8ac62dd1=stop.json"])
        self.assertEqual(device_coupons.requirement_options(["5ed91f8890c0=20", "abcd=18.5"], float, "--support-span-cap"),
                         {"5ed91f8890c0": 20.0, "abcd": 18.5})
        self.assertEqual(device_coupons.requirement_options(["5ed91f8890c0=16000000"], int, "--element-cap"),
                         {"5ed91f8890c0": 16000000})
        for bad in ("5ed91f8890c0", "=20", "XYZ=20", "5ed9=0", "5ed9=-3", "5ed9=sixteen"):
            with self.assertRaises(device_coupons.DeviceAdapterError):
                device_coupons.requirement_options([bad], float, "--support-span-cap")
        with self.assertRaisesRegex(device_coupons.DeviceAdapterError, "given twice"):
            device_coupons.requirement_options(["5ed9=20", "5ed9=21"], float, "--support-span-cap")
        options = {"5ed91f8890c0": 20.0}
        full = "5ed91f8890c014ae48ce6958221a2ecc7f0b796a59cc40470fb95a1d1a3690ed"
        self.assertEqual(device_coupons.requirement_option_for(options, full), ("5ed91f8890c0", 20.0))
        self.assertIsNone(device_coupons.requirement_option_for(options, "bd43654a77c6" + "0" * 52))
        with self.assertRaisesRegex(device_coupons.DeviceAdapterError, "several prefixes"):
            device_coupons.requirement_option_for({"5ed9": 1.0, "5ed91f": 2.0}, full)
        device_coupons.check_requirement_options_consumed(options, ["5ed91f8890c0"], "--support-span-cap")
        with self.assertRaisesRegex(device_coupons.DeviceAdapterError, "no spatial coupon of the discovery"):
            device_coupons.check_requirement_options_consumed(options, [], "--support-span-cap")
        # The texts are mandatory with the options (fail closed before any discovery).
        for kwargs in ({"support_span_caps": ["5ed9=20"]},
                       {"element_caps": ["5ed9=16000000"], "element_cap_reason": "r"},
                       {"element_caps": ["5ed9=16000000"], "element_cap_approval": "a"}):
            with self.assertRaisesRegex(device_coupons.DeviceAdapterError, "needs --"):
                device_coupons.prepare_device_sources("device.json", palace="palace", output=Path("/nonexistent/out"),
                                                      manifest_path="m.json", **kwargs)
        # The generator receives the cap for the named coupon only.
        self.assertIn("--support-span-cap", device_coupons.generate_sources.__doc__)


@unittest.skipUnless(available(), "the Palace executable and the transmon fixture are needed")
class DeviceCouponsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="coupon-device-test-"))
        cls.device = device_config(cls.tmp)
        manifest = json.loads((HERE / "geometry-independence-suite.json").read_text())
        manifest["RepositoryRoot"] = str((HERE / manifest["RepositoryRoot"]).resolve())
        cls.manifest_path = cls.tmp / "manifest.json"
        cls.manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        cls.record = device_coupons.prepare_device_sources(cls.device, palace=PALACE, output=cls.tmp / "device",
                                                           manifest_path=cls.manifest_path, log=lambda message: None)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, True)

    def test_discovery_maps_to_content_hashed_spatial_source_directories(self):
        record = self.record
        self.assertEqual(record["Device"]["SHA256"], case_inputs.sha256(self.device))
        self.assertTrue(record["Discovery"]["Complete"] is False)     # against the empty seed every requirement is Missing
        self.assertGreaterEqual(len(record["Coupons"]), 3)
        self.assertTrue(record["OutOfScope"])
        self.assertEqual({item["Method"] for item in record["OutOfScope"]}, {"CornerCoupon", "StraightEdgeBuilder"})
        # The device basis default is the Delaunay cap triangulation (decision 57), recorded in
        # the device record, every basis contract and every provenance.
        self.assertEqual(record["TraceBasis"]["CapTriangulation"], "delaunay")
        for coupon in record["Coupons"]:
            directory = Path(coupon["Directory"])
            self.assertEqual(directory.name, coupon["Case"])
            contract = json.loads((directory / "basis-contract.json").read_text())
            self.assertEqual(contract["CapTriangulation"]["Method"], "delaunay", coupon["Case"])
            provenance = json.loads((directory / "provenance.json").read_text())
            self.assertEqual(provenance["Generator"]["CapTriangulation"], "delaunay")
            self.assertIn("delaunay", provenance["Generator"]["Command"])
            digest, digests = device_coupons.content_hash(directory)
            self.assertEqual(coupon["Case"], f"spatial-{coupon['EdgeCount']}-edge-{digest[:12]}")
            self.assertEqual(set(digests), set(device_coupons.CONTENT_ROLES))
            self.assertEqual(provenance["ContentHash"]["SHA256"], digest)
            self.assertEqual(provenance["EtchFootprint"]["Declared"], "producer-default")
            self.assertFalse((directory / "retained-etch.csv").exists())
            self.assertFalse((directory / "traces").exists())
            # The trace basis loads in the mesh frame; the contract records its own digests.
            basis = trace_basis.load_trace_basis(directory / "basis-contract.json", directory / "trace-vertices.csv",
                                                 directory / "trace-triangles.csv", directory / "process-library.json")
            contract = json.loads((directory / "basis-contract.json").read_text())
            self.assertEqual(int(basis["Basis"].max()) + contract["ConductorStates"], contract["Sources"])
            self.assertEqual(contract["FrameFitResidual"], 0.0)
            self.assertEqual(coupon["Sources"], contract["Sources"])
        # The v2 cluster contract: every device coupon is built in the canonical frame of its
        # signature; the model carries the record's Signature (the matcher's exact key) and,
        # as Edges, the claimed portions exactly; the mask covers every conductor, its
        # physical boundary segments lie on the portions' lines, and the mesh-signature rows
        # (the generator's frame) reproduce the portions' lines up to the frame rotation.
        closure = json.loads(Path(record["Discovery"]["Manifest"]).read_text())
        by_id = {planner.coupon_id(requirement): requirement for requirement in closure["Requirements"]}
        radius = record["ProcessLibrary"]["MatchingRadius"]
        for coupon in record["Coupons"]:
            requirement = by_id[coupon["Requirement"]]
            signature = requirement["Signature"]
            self.assertEqual(signature["Type"], "SpatialEdgeCluster")
            self.assertEqual(coupon["EdgeCount"], len(signature["Portions"]))
            directory = Path(coupon["Directory"])
            model = json.loads((directory / "process-library.json").read_text())["Models"][0]
            self.assertEqual(model["Signature"], signature)
            portions = sorted(tuple(round(v * radius, 6) for v in p["P"]) for p in signature["Portions"])
            reconstructed = []
            for edge in model["Edges"]:
                gap, point = np.asarray(edge["GapDirection"]), np.asarray(edge["Point"])
                tangent = np.cross(gap, np.asarray(edge["ProcessNormal"]))
                a, b = [(point + s * tangent)[:2] for s in edge["Interval"]]
                a, b = (tuple(round(float(v), 6) for v in a), tuple(round(float(v), 6) for v in b))
                reconstructed.append(min(a, b) + max(a, b))
            self.assertEqual(sorted(reconstructed), portions, coupon["Case"])
            # Spatial-support contract v3 (decision 282): a signature with Box + Context is
            # built in its own frame from the claims plus the context pieces (continuation
            # chains, foreign edges — here the fixture's second lead and the port-cut end
            # edges); the mask then covers the context conductors too, the model carries
            # the SupportBox / ContextEdges / ForeignEdges and the mesh-signature rows are
            # the exact claims plus the context rows (Context / Chain columns).
            context = signature.get("Context", [])
            conductors = {int(p["Conductor"]) for p in signature["Portions"] + context}
            areas = mask_areas(directory / "plan-view-mask.csv")
            self.assertEqual({conductor for conductor, _ in areas}, conductors, coupon["Case"])
            self.assertTrue(all(area > 0.0 for area in areas.values()))
            with open(directory / "plan-view-boundary.csv", newline="") as stream:
                loops = list(csv.DictReader(stream))
            with open(directory / "mesh-signature.csv", newline="") as stream:
                rows = list(csv.DictReader(stream))
            if "Box" in signature:
                self.assertEqual([v * radius for v in signature["Box"]], model["SupportBox"])
                self.assertEqual(len(model["ContextEdges"]), len(rows) - coupon["EdgeCount"])
                self.assertEqual(sum(int(r["Context"]) for r in rows), len(rows) - coupon["EdgeCount"])
                self.assertEqual(len(model["ForeignEdges"]), sum(1 for c in context if not c["Chain"]))
                self.assertGreater(len(model["ContextEdges"]), 0)
            else:
                self.assertNotIn("SupportBox", model)
            lines = [((float(r["Px"]), float(r["Py"])), (float(r["Tx"]), float(r["Ty"]))) for r in rows]
            rows = [r for r in rows if not int(r.get("Context", 0))]
            self.assertEqual(len(rows), coupon["EdgeCount"])
            physical = [row for row in loops if row["Class"] == "Physical"]
            self.assertTrue(physical)
            for row in physical:
                x, y = float(row["X"]), float(row["Y"])
                self.assertTrue(any(abs((x - px) * ty - (y - py) * tx) <= 1e-6 for (px, py), (tx, ty) in lines),
                                f"{coupon['Case']}: physical boundary vertex ({x}, {y}) off every edge line")

    def test_omitted_requirement_is_recorded_out_of_scope_and_not_built(self):
        """prepare_device_sources(omit_requirements=[prefix]): the discovery gives the matching
        Missing requirement no placeholder, the adapter builds no source directory for it (whatever
        its family) and records it out of scope with Method OmittedRequirement, its family method
        and the prefix; the device record's Discovery lists it; every other coupon is unchanged."""
        omitted = min(self.record["Coupons"], key=lambda coupon: coupon["Sources"])
        closure = json.loads(Path(self.record["Discovery"]["Manifest"]).read_text())
        by_id = {planner.coupon_id(requirement): requirement for requirement in closure["Requirements"]}
        prefix = by_id[omitted["Requirement"]]["Hash"][:12]
        record = device_coupons.prepare_device_sources(self.device, palace=PALACE, output=self.tmp / "device-omit",
                                                       manifest_path=self.manifest_path, log=lambda message: None,
                                                       omit_requirements=[prefix])
        self.assertEqual([coupon["Case"] for coupon in record["Coupons"]],
                         [coupon["Case"] for coupon in self.record["Coupons"] if coupon["Case"] != omitted["Case"]])
        self.assertFalse((self.tmp / "device-omit" / "sources" / omitted["Case"]).exists())
        out_of_scope = [item for item in record["OutOfScope"] if item["Method"] == device_coupons.OMITTED_METHOD]
        self.assertEqual(len(out_of_scope), 1)
        self.assertEqual(out_of_scope[0]["Id"], omitted["Requirement"])
        self.assertEqual(out_of_scope[0]["FamilyMethod"], device_coupons.SPATIAL_METHOD)
        self.assertEqual(out_of_scope[0]["Hash"][:12], prefix)
        self.assertIn(prefix, out_of_scope[0]["Reason"])
        self.assertEqual([item["Hash"][:12] for item in record["Discovery"]["OmittedRequirements"]], [prefix])
        self.assertEqual(record["Discovery"]["OmittedRequirements"][0]["Topology"], "SpatialEdgeCluster")
        closure_omit = json.loads(Path(record["Discovery"]["Manifest"]).read_text())
        self.assertEqual([item["Hash"] for item in closure_omit["OmittedRequirements"]],
                         [record["Discovery"]["OmittedRequirements"][0]["Hash"]])
        self.assertEqual({requirement["Status"] for requirement in closure_omit["Requirements"]
                          if requirement.get("Hash", "").startswith(prefix)}, {"Missing"})
        history = json.loads((self.tmp / "device-omit" / "discovery" / "closure-history.json").read_text())
        self.assertNotIn(by_id[omitted["Requirement"]]["Hash"][:16],
                         {placeholder["Id"] for entry in history["Passes"] for placeholder in entry["AddedPlaceholders"]})

    def test_rerun_is_idempotent_by_content(self):
        again = device_coupons.prepare_device_sources(self.device, palace=PALACE, output=self.tmp / "device",
                                                      manifest_path=self.manifest_path, log=lambda message: None)
        self.assertEqual([coupon["Case"] for coupon in again["Coupons"]], [coupon["Case"] for coupon in self.record["Coupons"]])
        self.assertTrue(all(coupon["Status"] == "existing" for coupon in again["Coupons"]))

    def test_registration_and_preflight_of_the_two_smallest(self):
        """The two smallest device coupons register (stub census probe from the provisional
        contract's labels: the mesher is the live run's evidence) into a manifest copy,
        idempotently by content, with InventoryStatus DeviceDerived and the shared mesh
        recipe; the manifest copy passes the matrix preflight."""
        record = json.loads(json.dumps(self.record))
        # The two smallest by plan rows (claims + context: the device-plan coupon of the
        # fixture's JJ carries the second lead and exceeds the headroom gate) then sources.
        record["Coupons"] = sorted(record["Coupons"],
                                   key=lambda coupon: (coupon["EdgeCount"] + coupon["ContextEdgeCount"], coupon["Sources"]))[:2]
        labels = {}
        for coupon in record["Coupons"]:
            directory = Path(coupon["Directory"])
            for kind in ("fabricated", "thin"):
                provisional = derive_contract(directory, None, signature=directory / "mesh-signature.csv",
                                              boundary=directory / "plan-view-boundary.csv",
                                              process_library=directory / "process-library.json", kind=kind)
                case_id = coupon["Case"] if kind == "fabricated" else register_case.thin_case_id(coupon["Case"])
                labels[case_id] = [item["Attribute"] for item in provisional["BoundaryLabels"]]

        def probe(probe_manifest, case_id, root, log):
            return census_probe(labels[case_id])(probe_manifest, case_id, root, log)
        registered = device_coupons.register_device_sources(record, manifest_path=self.manifest_path,
                                                            work=self.tmp / "register", log=lambda message: None, probe=probe,
                                                            jobs=3)
        manifest = json.loads(self.manifest_path.read_text())
        by_id = {case["Id"]: case for case in manifest["Cases"]}
        # Decision 62(2): the probes / derivations ran in a pool of 3, the manifest appends
        # serially in the device record's coupon order.
        self.assertEqual(registered["RegistrationPool"]["Jobs"], 3)
        # Decision 66: every fabricated device coupon is followed by its thin pair (Kind
        # thin, FabricatedCase, the same directory, its own thin contract file).
        self.assertEqual([case["Id"] for case in manifest["Cases"] if case["Id"] in by_id and case["Id"].startswith("spatial-")],
                         [coupon["Case"] for coupon in registered["Coupons"]]
                         + [coupon["ThinCase"] for coupon in registered["Coupons"]])
        for coupon in registered["Coupons"]:
            self.assertEqual(coupon["ThinRegistration"]["Status"], register_case.STATUS_REGISTERED, coupon["ThinRegistration"])
            thin = by_id[coupon["ThinCase"]]
            self.assertEqual((thin["Kind"], thin["FabricatedCase"]), ("thin", coupon["Case"]))
            self.assertEqual(thin["Source"]["Directory"], by_id[coupon["Case"]]["Source"]["Directory"])
            self.assertEqual(thin["Source"]["Files"]["SemanticContract"]["Name"], register_case.THIN_SEMANTIC_CONTRACT_FILE)
            self.assertIn("ThinMetal", thin["Features"])
            thin_contract = json.loads((Path(coupon["Directory"]) / register_case.THIN_SEMANTIC_CONTRACT_FILE).read_text())
            self.assertEqual(sorted(item["Attribute"] // 1000 for item in thin_contract["BoundaryLabels"]),
                             sorted([0, 3] + [4] * (len(thin_contract["BoundaryLabels"]) - 2)))
        with self.assertRaises(device_coupons.DeviceAdapterError):
            device_coupons.register_device_sources(record, manifest_path=self.manifest_path, work=self.tmp / "register-0",
                                                   log=lambda message: None, probe=probe, jobs=0)
        parser = coupon_library.build_parser()
        self.assertEqual(parser.parse_args(["build", "--root", "r"]).register_jobs, device_coupons.DEFAULT_REGISTER_JOBS)
        self.assertEqual(parser.parse_args(["build", "--root", "r", "--register-jobs", "3"]).register_jobs, 3)
        for coupon in registered["Coupons"]:
            self.assertEqual(coupon["Registration"]["Status"], register_case.STATUS_REGISTERED, coupon["Registration"])
            case = by_id[coupon["Case"]]
            self.assertEqual(case["InventoryStatus"], device_coupons.INVENTORY_STATUS)
            self.assertEqual(case["Source"]["EtchFootprint"], "producer-default")
            self.assertEqual(case["Source"]["Files"]["MeshRecipe"]["RepositoryPath"], registered["MeshRecipe"])
            self.assertIn("TraceBasis", case["Features"])
            self.assertTrue((Path(coupon["Directory"]) / "semantic-contract.json").is_file())
        again = device_coupons.register_device_sources(record, manifest_path=self.manifest_path, work=self.tmp / "register-2",
                                                       log=lambda message: None, probe=probe)
        self.assertTrue(all(coupon["Registration"]["Status"] == register_case.STATUS_REUSED for coupon in again["Coupons"]))
        self.assertTrue(all(coupon["ThinRegistration"]["Status"] == register_case.STATUS_REUSED for coupon in again["Coupons"]))
        fresh = json.loads(json.dumps(self.record))
        fresh["Coupons"] = sorted(fresh["Coupons"],
                                  key=lambda coupon: (coupon["EdgeCount"] + coupon["ContextEdgeCount"], coupon["Sources"]))[:2]
        without = device_coupons.register_device_sources(fresh, manifest_path=self.manifest_path, work=self.tmp / "register-3",
                                                         log=lambda message: None, probe=probe, thin=False)
        self.assertTrue(all("ThinCase" not in coupon for coupon in without["Coupons"]))
        result = subprocess.run([sys.executable, str(HERE / "run_general_mesh_suite.py"), "--manifest", str(self.manifest_path),
                                 "--root", str(self.tmp / "preflight"), "--preflight-only"], text=True, capture_output=True)
        summary = json.loads((self.tmp / "preflight" / "summary.json").read_text())
        preflight = {case["Id"]: case for case in summary["Cases"]}
        # A device-plan coupon (decision 282) carries the plan inside its box (the fixture's
        # 4-edge: 19.8 R of ground edge, the box grown 0.5 R on two faces): its pre-build
        # estimate may exceed the suite's MaximumElements, a legitimate fail-closed outcome
        # of the headroom gate (recorded, never a silent pass); every other case passes.
        passed_spatial = 0
        for coupon in registered["Coupons"]:
            for case_id in (coupon["Case"], coupon["ThinCase"]):
                case = preflight[case_id]
                if not case["Passed"]:
                    self.assertIn("headroom gate", case.get("Error", ""), case)
                    self.assertFalse(case["BuildCostEstimate"]["Passed"])
                    continue
                self.assertTrue(case["BuildCostEstimate"]["Passed"])
                passed_spatial += 1
            # The thin estimate: one tube per side, fewer prisms than the fabricated pair.
            self.assertLess(preflight[coupon["ThinCase"]]["BuildCostEstimate"]["EstimatedPrisms"],
                            preflight[coupon["Case"]]["BuildCostEstimate"]["EstimatedPrisms"])
        self.assertGreaterEqual(passed_spatial, 2, result.stdout + result.stderr)
        registered_ids = {c for coupon in registered["Coupons"] for c in (coupon["Case"], coupon["ThinCase"])}
        if all(preflight[case_id]["Passed"] for case_id in registered_ids):
            self.assertTrue(summary["PreflightPassed"], result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
