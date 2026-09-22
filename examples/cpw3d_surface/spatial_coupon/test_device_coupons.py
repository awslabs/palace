# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""`coupon-library build --device` on the transmon example (the device the gallery cases
came from): discovery -> planner -> generate_spatial_response.py --basis-only ->
content-hashed source directories -> register_case (idempotent) -> preflight.  The
four-edge / three-edge / two-edge / ten-edge device coupons reproduce the gallery
cases' frozen geometry files byte for byte; the other families are recorded out of
scope.  Needs the Palace executable (build/bin/palace) for the geometry preflights."""
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
import case_inputs  # noqa: E402
import coupon_library  # noqa: E402
import device_coupons  # noqa: E402
import register_case  # noqa: E402
import trace_basis  # noqa: E402
from derive_semantic_contract import derive as derive_contract  # noqa: E402

REPOSITORY = HERE.parents[2]
PALACE = Path(os.environ.get("PALACE_EXECUTABLE", REPOSITORY / "build" / "bin" / "palace"))
TRANSMON = REPOSITORY / "examples" / "transmon"
DEVICE_CONFIG = TRANSMON / "transmon_surface_coarse.json"
PROCESS_SEED = TRANSMON / "benchmark" / "transmon_surface_process_seed.json"
# The gallery cases the transmon's discovery closure names (their producer model digests).
GALLERY_MODELS = {"9d2cb9bbb3fe": "four-edge-9d2cb9bbb3fe", "419576fdab24": "three-edge-419576fdab24",
                  "3f8992613e95": "two-edge-3f8992613e95", "6791f1c84123": "ten-edge-6791f1c84123"}
# Byte-identical to the gallery's frozen files; the plan-view mask's facet tessellation
# follows the device mesh (same footprint: equal facet area per conductor and plane).
GEOMETRY_FILES = ("mesh-signature.csv", "plan-view-boundary.csv", "process.toml")


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

    def coupon_of(self, model_digest):
        return next(coupon for coupon in self.record["Coupons"] if coupon["Requirement"].endswith(model_digest))

    def test_discovery_maps_to_content_hashed_spatial_source_directories(self):
        record = self.record
        self.assertEqual(record["Device"]["SHA256"], case_inputs.sha256(self.device))
        self.assertTrue(record["Discovery"]["Complete"] is False)     # against the empty seed every requirement is Missing
        self.assertGreaterEqual(len(record["Coupons"]), 4)
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
        for digest, case_id in GALLERY_MODELS.items():
            coupon = self.coupon_of(digest)
            for name in GEOMETRY_FILES:
                self.assertEqual(case_inputs.sha256(Path(coupon["Directory"]) / name),
                                 case_inputs.sha256(HERE / "testdata" / case_id / name), f"{case_id}/{name}")
            device_areas = mask_areas(Path(coupon["Directory"]) / "plan-view-mask.csv")
            gallery_areas = mask_areas(HERE / "testdata" / case_id / "plan-view-mask.csv")
            self.assertEqual(set(device_areas), set(gallery_areas), case_id)
            for key, value in gallery_areas.items():
                self.assertAlmostEqual(device_areas[key], value, delta=1e-9 * max(1.0, value), msg=f"{case_id} {key}")

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
        record["Coupons"] = sorted(record["Coupons"], key=lambda coupon: coupon["Sources"])[:2]
        labels = {}
        for coupon in record["Coupons"]:
            directory = Path(coupon["Directory"])
            provisional = derive_contract(directory, None, signature=directory / "mesh-signature.csv",
                                          boundary=directory / "plan-view-boundary.csv",
                                          process_library=directory / "process-library.json")
            labels[coupon["Case"]] = [item["Attribute"] for item in provisional["BoundaryLabels"]]

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
        self.assertEqual([case["Id"] for case in manifest["Cases"] if case["Id"] in by_id and case["Id"].startswith("spatial-")],
                         [coupon["Case"] for coupon in registered["Coupons"]])
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
        result = subprocess.run([sys.executable, str(HERE / "run_general_mesh_suite.py"), "--manifest", str(self.manifest_path),
                                 "--root", str(self.tmp / "preflight"), "--preflight-only"], text=True, capture_output=True)
        summary = json.loads((self.tmp / "preflight" / "summary.json").read_text())
        self.assertTrue(summary["PreflightPassed"], result.stdout + result.stderr)
        preflight = {case["Id"]: case for case in summary["Cases"]}
        for coupon in registered["Coupons"]:
            self.assertTrue(preflight[coupon["Case"]]["Passed"], preflight[coupon["Case"]])
            self.assertTrue(preflight[coupon["Case"]]["BuildCostEstimate"]["Passed"])


if __name__ == "__main__":
    unittest.main()
