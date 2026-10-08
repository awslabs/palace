#!/usr/bin/env python3

import copy
import csv
import importlib.util
import json
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np


CPW2D = Path(__file__).parent
PREPARE_SPEC = importlib.util.spec_from_file_location(
    "prepare_surface_response_coupons",
    CPW2D / "prepare_surface_response_coupons.py",
)
PREPARE = importlib.util.module_from_spec(PREPARE_SPEC)
PREPARE_SPEC.loader.exec_module(PREPARE)

SPATIAL_PATH = (
    CPW2D.parent / "cpw3d_surface" / "spatial_coupon"
    / "generate_spatial_response.py"
)
SPATIAL_SPEC = importlib.util.spec_from_file_location(
    "generate_spatial_response", SPATIAL_PATH
)
SPATIAL = importlib.util.module_from_spec(SPATIAL_SPEC)
SPATIAL_SPEC.loader.exec_module(SPATIAL)

CORNER_PATH = (
    CPW2D.parent / "cpw3d_surface" / "corner_coupon"
    / "generate_corner_response.py"
)
CORNER_SPEC = importlib.util.spec_from_file_location(
    "generate_corner_response", CORNER_PATH
)
CORNER = importlib.util.module_from_spec(CORNER_SPEC)
CORNER_SPEC.loader.exec_module(CORNER)

COMBINER_PATH = (
    CPW2D.parent / "cpw3d_surface" / "corner_coupon"
    / "combine_process_libraries.py"
)
COMBINER_SPEC = importlib.util.spec_from_file_location(
    "combine_process_libraries", COMBINER_PATH
)
COMBINER = importlib.util.module_from_spec(COMBINER_SPEC)
COMBINER_SPEC.loader.exec_module(COMBINER)

# The stored stage-2 library of record and the round-2 concave-60 node (decision 505 Q1: the name repeat
# combine refuses; --supersede is the tool for the next library version); records only, matrices on the cluster.
ASSESSMENT = Path(os.environ.get("COUPON_ASSESSMENT_ROOT", CPW2D.parents[2] / "coupon-accuracy-assessment-20260913"))
V3B1_LIBRARY = ASSESSMENT / "stage2-20261004" / "coupons-bc" / "library" / "s2-r1p9-v3-b1" / "process-library.json"
ROUND2_CONCAVE_60 = (ASSESSMENT / "curved-clusters-20261005" / "family4" / "round2" / "concavecorner-60-r2"
                     / "corner-d6def3107d2c886bd4f9a1edfff823cc93c8b1c8b02be861409d49d7ce3cd11b")

CORNER_COMPARE_PATH = (
    CPW2D.parent / "cpw3d_surface" / "corner_coupon"
    / "compare_probe_convergence.py"
)
CORNER_COMPARE_SPEC = importlib.util.spec_from_file_location(
    "compare_probe_convergence", CORNER_COMPARE_PATH
)
CORNER_COMPARE = importlib.util.module_from_spec(CORNER_COMPARE_SPEC)
CORNER_COMPARE_SPEC.loader.exec_module(CORNER_COMPARE)
sys.modules["compare_probe_convergence"] = CORNER_COMPARE

CORNER_CONVERGENCE_PATH = (
    CPW2D.parent / "cpw3d_surface" / "corner_coupon"
    / "run_probe_convergence.py"
)
CORNER_CONVERGENCE_SPEC = importlib.util.spec_from_file_location(
    "run_probe_convergence", CORNER_CONVERGENCE_PATH
)
CORNER_CONVERGENCE = importlib.util.module_from_spec(CORNER_CONVERGENCE_SPEC)
CORNER_CONVERGENCE_SPEC.loader.exec_module(CORNER_CONVERGENCE)

CORNER_INTERPOLATION_PATH = (
    CPW2D.parent / "cpw3d_surface" / "corner_coupon"
    / "qualify_corner_interpolation.py"
)
CORNER_INTERPOLATION_SPEC = importlib.util.spec_from_file_location(
    "qualify_corner_interpolation", CORNER_INTERPOLATION_PATH
)
CORNER_INTERPOLATION = importlib.util.module_from_spec(
    CORNER_INTERPOLATION_SPEC
)
CORNER_INTERPOLATION_SPEC.loader.exec_module(CORNER_INTERPOLATION)


def pec():
    return {"Type": "PEC"}


def interfaces(slot=0):
    return [
        {"Slot": slot, "Type": interface, "Target": index}
        for index, interface in enumerate(("SA", "MS", "MA"), start=1)
    ]


def endpoint_coupon():
    return {
        "Id": "endpoint",
        "Topology": "Endpoint",
        "Geometry": {
            "SignatureVersion": 1,
            "ArmCount": 1,
            "ArmAnglesDegrees": [],
            "Arms": [
                {
                    "Direction": [1.0, 0.0, 0.0],
                    "GapDirection": [0.0, 1.0, 0.0],
                    "ProcessNormal": [0.0, 0.0, 1.0],
                    "Interval": [0.0, 2.0],
                    "Conductor": 1,
                    "InterfaceSlot": 0,
                    "BoundaryCondition": pec(),
                }
            ],
        },
        "Interfaces": interfaces(),
        "BoundaryCondition": pec(),
        "CoverageStatus": "Missing",
    }


def process_parameters():
    return {
        "metal_thickness": 0.1,
        "overetch": 0.05,
        "sidewall_angle": 80.0,
        "top_radius": 0.01,
        "bottom_radius": 0.01,
        "substrate_permittivity": 11.47,
        "sa_thickness": 0.002,
        "sa_permittivity": 4.0,
        "ms_thickness": 0.0003,
        "ms_permittivity": 11.47,
        "ma_thickness": 0.03,
        "ma_permittivity": 10.0,
    }


def overlapping_spatial_coupon(masked=False):
    geometry = {
        "Edges": [
            {
                "Point": [0.0, 0.0, -2.0],
                "GapDirection": [0.0, 0.0, -1.0],
                "ProcessNormal": [0.0, 1.0, 0.0],
                "Interval": [-2.0, 2.0],
                "Conductor": 1,
                "InterfaceSlot": 0,
                "BoundaryCondition": pec(),
            },
            {
                "Point": [0.0, 0.0, 0.0],
                "GapDirection": [1.0, 0.0, 0.0],
                "ProcessNormal": [0.0, 1.0, 0.0],
                "Interval": [0.0, 2.0],
                "Conductor": 2,
                "InterfaceSlot": 0,
                "BoundaryCondition": pec(),
            },
        ]
    }
    if masked:
        geometry["PlanViewFacets"] = [
            {
                "Conductor": 1,
                "Points": [
                    [-2.0, 0.0, -2.0],
                    [2.0, 0.0, -2.0],
                    [2.0, 0.0, -1.0],
                    [-2.0, 0.0, -1.0],
                ],
            },
            {
                "Conductor": 2,
                "Points": [
                    [-1.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 2.0],
                    [-1.0, 0.0, 2.0],
                ],
            },
        ]
    return {
        "Topology": "SpatialEdgeCluster",
        "Geometry": geometry,
        "Interfaces": interfaces(),
        "BoundaryCondition": pec(),
    }


class PrepareSurfaceResponseCouponsTest(unittest.TestCase):
    def test_combiner_accepts_metadata_only_seed(self):
        with tempfile.TemporaryDirectory() as directory:
            seed = Path(directory) / "process-library.json"
            seed.write_text(
                '{"Version": 3, "MatchingRadius": 2.0, '
                '"Fabrication": {"InterfaceLayers": {}}, "Models": []}\n'
            )
            path, library = COMBINER.load_library(seed)
            self.assertEqual(path, seed.resolve())
            self.assertEqual(library["Models"], [])

    def test_combiner_can_preserve_metadata_when_every_coupon_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            seed = root / "seed.json"
            seed.write_text(
                '{"Version": 3, "MatchingRadius": 2.0, '
                '"Fabrication": {"InterfaceLayers": {}}, "Models": []}\n'
            )
            output = root / "combined"
            with mock.patch.object(
                sys,
                "argv",
                [
                    "combine_process_libraries.py",
                    "--output",
                    str(output),
                    "--allow-empty",
                    str(seed),
                ],
            ):
                COMBINER.main()
            combined = PREPARE.load_json(output / "process-library.json")
            self.assertEqual(combined["Version"], 3)
            self.assertTrue(combined["ExhaustiveSpatialClosure"])
            self.assertEqual(combined["Models"], [])
            self.assertEqual(combined["Fabrication"], {"InterfaceLayers": {}})

    def test_combiner_records_source_libraries_and_refuses_a_null_matrix_path(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            header = '"Version": 3, "MatchingRadius": 2.0, "Fabrication": {"InterfaceLayers": {}}, '
            (root / "a").mkdir()
            (root / "a" / "fab.csv").write_text("basis_i\n")
            (root / "a" / "thin.csv").write_text("basis_i\n")
            first = root / "a" / "process-library.json"
            first.write_text(
                '{' + header + '"Name": "spatial", "Models": [{"Name": "cluster", '
                '"FabricatedMatrix": "fab.csv", "ThinMatrix": "thin.csv"}]}\n'
            )
            second = root / "b.json"
            second.write_text(
                '{' + header + '"Name": "corners", "Models": [{"Name": "corner", '
                '"FabricatedMatrix": "' + str(root / "a" / "fab.csv") + '", "ThinMatrix": null, '
                '"NotLoadable": {"Reason": "no thin matrices"}}]}\n'
            )
            output = root / "combined"
            argv = ["combine_process_libraries.py", "--output", str(output), str(first)]
            with mock.patch.object(sys, "argv", argv):
                COMBINER.main()
            combined = PREPARE.load_json(output / "process-library.json")
            resolved = str(first.resolve())
            self.assertEqual(
                combined["Sources"],
                [{"Path": resolved, "SHA256": COMBINER.sha256(first), "Name": "spatial",
                  "Version": 3, "Models": ["cluster"]}],
            )
            self.assertEqual(combined["Models"][0]["CombinedFrom"],
                             {"Path": resolved, "SHA256": COMBINER.sha256(first)})
            self.assertTrue((output / combined["Models"][0]["ThinMatrix"]).is_file())
            with mock.patch.object(sys, "argv", argv + [str(second)]):
                with self.assertRaisesRegex(ValueError, "corner field ThinMatrix is None.*no thin matrices"):
                    COMBINER.main()

    def test_corner_radius_interpolation_qualification(self):
        metadata = {
            "MatchingRadius": 2.0,
            "Fabrication": {"MetalThickness": 0.1},
            "Topology": "ConvexCorner",
            "Angle": 90.0,
            "ContourGroups": [3],
            "Interfaces": [{"Type": "SA", "Coupon": 1}],
            "BoundaryCondition": "PEC",
        }
        thin_domain = np.diag([2.0, 3.0])
        fabricated_domain = np.diag([3.0, 4.0])

        def case(name, radius, surface):
            return {
                "Name": name,
                "Root": Path(name),
                "ModelName": f"corner-r{radius:g}",
                "Radius": radius,
                "Metadata": metadata,
                "Basis": np.array([[0.0, 0.0], [1.0, 0.0]]),
                "Coefficients": np.array([1.0, 0.5]),
                "Active": np.array([0, 1]),
                "Responses": {
                    "thin": {
                        "domain": thin_domain,
                        "surfaces": {1: surface},
                    },
                    "fabricated": {
                        "domain": fabricated_domain,
                        "surfaces": {1: surface},
                    },
                },
                "InterfaceNames": {1: "SA"},
            }

        lower_surface = np.diag([1.0, 2.0])
        upper_surface = np.diag([3.0, 4.0])
        report = CORNER_INTERPOLATION.qualify_cases(
            case("lower", 0.25, lower_surface),
            case("heldout", 0.5, 0.5 * (lower_surface + upper_surface)),
            case("upper", 0.75, upper_surface),
            5.0,
            10.0,
        )

        self.assertTrue(report["Passed"])
        self.assertEqual(report["Weights"], {"Lower": 0.5, "Upper": 0.5})
        self.assertEqual(
            report["LibraryRecord"]["Qualification"]["Method"],
            "HeldOutCoupon",
        )

    def test_corner_interpolation_fixed_flux_uses_active_subspace(self):
        case = {
            "Basis": np.zeros((2, 2)),
            "Coefficients": np.array([1.0, 7.0]),
            "Active": np.array([0]),
            "Responses": {
                "thin": {"domain": np.array([[4.0, 20.0], [20.0, 6.0]])},
                "fabricated": {
                    "domain": np.array([[2.0, 10.0], [10.0, 3.0]])
                },
            },
        }
        np.testing.assert_allclose(
            CORNER_INTERPOLATION.fixed_flux(case),
            np.array([2.0, 0.0]),
        )
        np.testing.assert_allclose(
            CORNER_INTERPOLATION.fixed_flux(case, np.array([2.0, 9.0])),
            np.array([4.0, 0.0]),
        )

    def test_combiner_embeds_only_passed_corner_interpolation_reports(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fabrication = {"InterfaceLayers": {}}

            def library(path, model):
                path.write_text(
                    json.dumps(
                        {
                            "Version": 3,
                            "MatchingRadius": 2.0,
                            "Fabrication": fabrication,
                            "Models": [{"Name": model}],
                        }
                    )
                    + "\n"
                )

            lower = root / "lower.json"
            upper = root / "upper.json"
            library(lower, "corner-r0.25")
            library(upper, "corner-r0.75")
            report = root / "interpolation.json"
            record = {
                "LowerModel": "corner-r0.25",
                "UpperModel": "corner-r0.75",
                "Qualification": {
                    "Method": "HeldOutCoupon",
                    "Passed": True,
                    "HeldoutRadius": 0.5,
                },
            }
            report.write_text(
                json.dumps(
                    {
                        "Study": "CornerRadiusInterpolation",
                        "Passed": True,
                        "LibraryRecord": record,
                    }
                )
                + "\n"
            )
            output = root / "combined"
            with mock.patch.object(
                sys,
                "argv",
                [
                    "combine_process_libraries.py",
                    "--output",
                    str(output),
                    "--corner-interpolation-qualification",
                    str(report),
                    str(lower),
                    str(upper),
                ],
            ):
                COMBINER.main()
            combined = PREPARE.load_json(output / "process-library.json")
            self.assertEqual(combined["CornerRadiusInterpolation"], [record])

            failed = root / "failed-interpolation.json"
            failed.write_text(
                json.dumps(
                    {
                        "Study": "CornerRadiusInterpolation",
                        "Passed": False,
                        "LibraryRecord": record,
                    }
                )
                + "\n"
            )
            with mock.patch.object(
                sys,
                "argv",
                [
                    "combine_process_libraries.py",
                    "--output",
                    str(root / "rejected"),
                    "--corner-interpolation-qualification",
                    str(failed),
                    str(lower),
                    str(upper),
                ],
            ):
                with self.assertRaisesRegex(
                    ValueError, "not a passed corner-radius interpolation report"
                ):
                    COMBINER.main()

    @staticmethod
    def write_library(path, name, models, *, root=None):
        """A version-3 library at `path` whose models' matrix files are written beside it (content = the file name)."""
        root = root or path.parent
        root.mkdir(parents=True, exist_ok=True)
        entries = []
        for model_name in models:
            entry = {"Name": model_name}
            for field, filename in COMBINER.PATH_NAMES.items():
                relative = Path("src") / model_name / filename
                (root / relative).parent.mkdir(parents=True, exist_ok=True)
                (root / relative).write_text(f"{name}:{model_name}:{filename}\n")
                entry[field] = str(relative)
            entries.append(entry)
        path.write_text(json.dumps({"Version": 3, "MatchingRadius": 2.0, "Name": name,
                                    "Fabrication": {"InterfaceLayers": {}}, "Models": entries}, indent=2) + "\n")
        return path

    def combine(self, output, *argv):
        with mock.patch.object(sys, "argv", ["combine_process_libraries.py", "--output", str(output), *map(str, argv)]):
            COMBINER.main()
        return PREPARE.load_json(output / "process-library.json")

    def test_supersede_replaces_the_named_base_model_in_place_and_records_it(self):
        """Decision 505 Q1: --supersede NAME=<record> replaces the base library's model NAME by the later input's
        model of the same name AT THE BASE MODEL'S INDEX (every other model and its directory byte-identical to
        the flag-less combine of the base alone), records Supersedes {Name, BaseModelSHA, NewModelSHA, Record,
        Decision}; without the flag the repeated name is refused as before."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            base = self.write_library(root / "base" / "process-library.json", "base", ["alpha", "concave-corner-60deg", "omega"])
            round2 = self.write_library(root / "round2" / "process-library.json", "round2", ["concave-corner-60deg"])
            record = root / "concavecorner-60-r2-entry.json"
            record.write_text(json.dumps({"Node": "r2"}) + "\n")
            decision = "328 (3) / 505 Q1: the round-2 concave 60 replaces the round-1 node"
            plain = self.combine(root / "plain", base)
            with self.assertRaisesRegex(ValueError, "repeated"):
                self.combine(root / "refused", base, round2)
            superseded = self.combine(root / "superseded", "--supersede", f"concave-corner-60deg={record}", "--supersede-decision", decision,
                                      base, round2)
            self.assertEqual([model["Name"] for model in superseded["Models"]], ["alpha", "concave-corner-60deg", "omega"])
            self.assertEqual([model["Name"] for model in plain["Models"]], [model["Name"] for model in superseded["Models"]])
            # Every other model entry and file byte-identical to the flag-less combine of the base; the superseded
            # one is the round-2 model at the base slot (index 002), its files the round-2 files.
            for index, (before, after) in enumerate(zip(plain["Models"], superseded["Models"])):
                if before["Name"] != "concave-corner-60deg":
                    self.assertEqual(before, after)
                    for field in COMBINER.PATH_FIELDS:
                        self.assertEqual((root / "plain" / before[field]).read_bytes(), (root / "superseded" / after[field]).read_bytes())
                else:
                    self.assertEqual(after["CombinedFrom"], {"Path": str(round2.resolve()), "SHA256": COMBINER.sha256(round2)})
                    self.assertEqual(after["FabricatedMatrix"], "models/002-concave-corner-60deg/fabricated-domain-response-matrix.csv")
                    self.assertEqual((root / "superseded" / after["FabricatedMatrix"]).read_text(),
                                     "round2:concave-corner-60deg:fabricated-domain-response-matrix.csv\n")
                    self.assertEqual((root / "plain" / before["FabricatedMatrix"]).read_text(),
                                     "base:concave-corner-60deg:fabricated-domain-response-matrix.csv\n")
            self.assertEqual([entry["Name"] for entry in superseded["Supersedes"]], ["concave-corner-60deg"])
            entry = superseded["Supersedes"][0]
            base_model = PREPARE.load_json(base)["Models"][1]
            new_model = PREPARE.load_json(round2)["Models"][0]
            self.assertEqual((entry["BaseModelSHA"], entry["NewModelSHA"]),
                             (COMBINER.model_sha256(base_model, base.parent), COMBINER.model_sha256(new_model, round2.parent)))
            self.assertEqual(COMBINER.model_sha256(base_model, base.parent), COMBINER.model_sha256(dict(base_model), root / "base"))
            self.assertNotEqual(entry["BaseModelSHA"], entry["NewModelSHA"])
            self.assertEqual(entry["Record"], {"Path": str(record.resolve()), "SHA256": COMBINER.sha256(record)})
            self.assertEqual(entry["Decision"], decision)        # the --supersede-decision text, never read from the record
            self.assertEqual((entry["BaseSource"]["SHA256"], entry["NewSource"]["SHA256"]), (COMBINER.sha256(base), COMBINER.sha256(round2)))
            self.assertEqual(entry["Rule"], COMBINER.SUPERSEDE_RULE)
            self.assertEqual([source["Name"] for source in superseded["Sources"]], ["base", "round2"])
            self.assertNotIn("Supersedes", plain)
            # Refusals: a name the base does not hold, a missing record, an unreadable record, a name no later
            # input (or two) supplies, --supersede without a decision text (decision 513 (2): mandatory, non-blank)
            # or a decision text without --supersede; nothing is written before the refusal.
            ruled = ("--supersede-decision", decision)
            for argv, message in (
                    (("--supersede", f"beta={record}", *ruled, base, round2), "has no model 'beta'"),
                    (("--supersede", f"concave-corner-60deg={root / 'absent.json'}", *ruled, base, round2), "is missing"),
                    (("--supersede", f"concave-corner-60deg={root / 'base' / 'src' / 'alpha' / 'basis-points.csv'}", *ruled, base, round2), "is unreadable"),
                    (("--supersede", f"concave-corner-60deg={record}", *ruled, base), "0 later inputs supply"),
                    (("--supersede", f"concave-corner-60deg={record}", *ruled, base, round2, round2), "2 later inputs supply"),
                    (("--supersede", f"alpha={record}", *ruled, base, round2), "0 later inputs supply"),
                    (("--supersede", "concave-corner-60deg", *ruled, base, round2), "not NAME=<record path>"),
                    (("--supersede", f"concave-corner-60deg={record}", base, round2), "requires --supersede-decision TEXT"),
                    (("--supersede", f"concave-corner-60deg={record}", "--supersede-decision", "  ", base, round2), "requires --supersede-decision TEXT"),
                    (("--supersede-decision", decision, base), "without --supersede")):
                with self.assertRaisesRegex(ValueError, message):
                    self.combine(root / "refusal", *argv)
                self.assertFalse((root / "refusal" / "process-library.json").exists(), message)

    @unittest.skipUnless(V3B1_LIBRARY.is_file() and (ROUND2_CONCAVE_60 / "process-library.json").is_file()
                         and (ROUND2_CONCAVE_60 / "heldout-qualification.json").is_file(),
                         "the stored s2-r1p9-v3-b1 library and the round-2 concave-60 node records are needed")
    def test_supersede_round2_concave_60_replaces_the_round1_node_of_v3_b1_on_copies(self):
        """Decision 505 Q1 on COPIES of the stored records (the matrices live on the cluster: stand-in files at
        the recorded relative paths): combining v3-b1 with the round-2 concave-corner-60deg node is REFUSED
        without the flag (the refusal of record), and --supersede concave-corner-60deg=<the node's held-out
        qualification record> with its --supersede-decision replaces the round-1 node (v3-b1's model 75) in place: every other model entry and
        file byte-identical to the flag-less combine of v3-b1 alone, Supersedes recorded against the record."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            def staged(source, destination):
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, destination)
                for model in PREPARE.load_json(source)["Models"]:
                    for relative in [model[f] for f in COMBINER.PATH_FIELDS if isinstance(model.get(f), str)] + list(model.get("TraceMesh", {}).values()):
                        self.assertFalse(Path(relative).is_absolute(), relative)
                        target = destination.parent / relative
                        target.parent.mkdir(parents=True, exist_ok=True)
                        if not target.exists():
                            target.write_text(f"{destination.parent.name}:{relative}\n")
                return destination
            base = staged(V3B1_LIBRARY, root / "v3-b1" / "process-library.json")
            round2 = staged(ROUND2_CONCAVE_60 / "process-library.json", root / "concavecorner-60-r2" / "process-library.json")
            record = ROUND2_CONCAVE_60 / "heldout-qualification.json"
            base_names = [model["Name"] for model in PREPARE.load_json(base)["Models"]]
            self.assertEqual(base_names.index("concave-corner-60deg"), 74)
            self.assertEqual([model["Name"] for model in PREPARE.load_json(round2)["Models"]], ["concave-corner-60deg"])
            with self.assertRaisesRegex(ValueError, "'concave-corner-60deg' is empty or repeated"):
                self.combine(root / "refused", "--name", "probe", base, round2)
            plain = self.combine(root / "plain", "--name", "s2-r1p9-v3-b1-copy", base)
            decision = ("decisions 328 (3) / 505 Q1 / 513 (2): the round-2 concave-corner-60deg node replaces the round-1 node of "
                        "s2-r1p9-v3-b1 in the next library version")
            with self.assertRaisesRegex(ValueError, "requires --supersede-decision TEXT"):
                self.combine(root / "undecided", "--name", "probe", "--supersede", f"concave-corner-60deg={record}", base, round2)
            superseded = self.combine(root / "superseded", "--name", "s2-r1p9-v3-b1-copy", "--supersede", f"concave-corner-60deg={record}",
                                      "--supersede-decision", decision, base, round2)
            self.assertEqual([model["Name"] for model in superseded["Models"]], base_names)
            self.assertEqual(len(superseded["Models"]), 81)
            for before, after in zip(plain["Models"], superseded["Models"]):
                if before["Name"] == "concave-corner-60deg":
                    self.assertEqual(after["CombinedFrom"]["SHA256"], COMBINER.sha256(round2))
                    self.assertEqual(after["FabricatedMatrix"], "models/075-concave-corner-60deg/fabricated-domain-response-matrix.csv")
                    self.assertTrue((root / "superseded" / after["FabricatedMatrix"]).read_text().startswith("concavecorner-60-r2:"))
                    self.assertTrue((root / "plain" / before["FabricatedMatrix"]).read_text().startswith("v3-b1:"))
                    self.assertEqual(after["Topology"], "ConcaveCorner")
                    continue
                self.assertEqual(before, after)
                for field in COMBINER.PATH_FIELDS:
                    if isinstance(before.get(field), str):
                        self.assertEqual((root / "plain" / before[field]).read_bytes(), (root / "superseded" / after[field]).read_bytes())
            entry = superseded["Supersedes"][0]
            self.assertEqual((entry["Name"], entry["Record"]), ("concave-corner-60deg", {"Path": str(record.resolve()), "SHA256": COMBINER.sha256(record)}))
            self.assertEqual(entry["Decision"], decision)
            self.assertNotEqual(entry["BaseModelSHA"], entry["NewModelSHA"])
            self.assertEqual(entry["BaseSource"]["SHA256"], COMBINER.sha256(base))
            self.assertEqual(entry["NewSource"]["SHA256"], COMBINER.sha256(round2))
            self.assertEqual({key: value for key, value in superseded.items() if key not in ("Models", "Sources", "Supersedes")},
                             {key: value for key, value in plain.items() if key not in ("Models", "Sources")})

    def test_material_options_emit_one_value_per_flag(self):
        options = PREPARE.material_options(process_parameters())
        self.assertEqual(
            options,
            [
                "--substrate-permittivity",
                11.47,
                "--sa-thickness",
                0.002,
                "--sa-permittivity",
                4.0,
                "--ms-thickness",
                0.0003,
                "--ms-permittivity",
                11.47,
                "--ma-thickness",
                0.03,
                "--ma-permittivity",
                10.0,
            ],
        )

    def test_plan_routes_every_supported_family(self):
        requirements = []
        for topology, geometry in (
            ("IsolatedEdge", {}),
            (
                "ParallelEdgeCluster",
                {
                    "EdgeCount": 3,
                    "Edges": [
                        {"Offset": 0.0, "GapDirection": 1, "Conductor": 1},
                        {"Offset": 1.0, "GapDirection": -1, "Conductor": 1},
                        {"Offset": 2.0, "GapDirection": 1, "Conductor": 2},
                    ],
                },
            ),
            ("ConvexCorner", {"AngleDegrees": 45.0, "CornerRadius": 0.2}),
            ("Endpoint", endpoint_coupon()["Geometry"]),
        ):
            requirements.append(
                {
                    "Topology": topology,
                    "Geometry": geometry,
                    "Interfaces": interfaces(),
                    "BoundaryCondition": pec(),
                    "Status": "Missing",
                }
            )
        manifest = {
            "Library": {"MatchingRadius": 2.0},
            "Requirements": requirements,
        }
        plan = PREPARE.plan_from_manifest(
            Path("requirements.json"),
            manifest,
            Path("library.json"),
            {"Fabrication": {}},
            False,
        )
        methods = {
            coupon["Topology"]: coupon["Preparation"]["Method"]
            for coupon in plan["Coupons"]
        }
        self.assertEqual(methods["IsolatedEdge"], "StraightEdgeBuilder")
        self.assertEqual(methods["ParallelEdgeCluster"], "ParallelClusterCoupon")
        self.assertEqual(methods["ConvexCorner"], "CornerCoupon")
        self.assertEqual(methods["Endpoint"], "SpatialCoupon")
        self.assertEqual(plan["Summary"]["Unsupported"], 0)

    def test_corner_planner_rejects_only_invalid_angles(self):
        for angle in (30.0, 90.0, 150.0):
            requirement = {
                "Topology": "ConcaveCorner",
                "Geometry": {"AngleDegrees": angle, "CornerRadius": 0.2},
                "BoundaryCondition": pec(),
            }
            self.assertEqual(
                PREPARE.preparation(requirement)["Method"], "CornerCoupon"
            )
        for angle in (0.0, float("nan")):
            requirement = {
                "Topology": "ConvexCorner",
                "Geometry": {"AngleDegrees": angle, "CornerRadius": 0.0},
                "BoundaryCondition": pec(),
            }
            self.assertEqual(
                PREPARE.preparation(requirement)["Method"], "Unsupported"
            )
        # 180 deg is the corner family's straight anchor (USER decision 121 (C)) — sharp only.
        for radius, method in ((0.0, "CornerCoupon"), (0.2, "Unsupported")):
            requirement = {
                "Topology": "ConvexCorner",
                "Geometry": {"AngleDegrees": 180.0, "CornerRadius": radius},
                "BoundaryCondition": pec(),
            }
            self.assertEqual(PREPARE.preparation(requirement)["Method"], method)

    def test_corner_config_edge_lines_are_the_sa_perimeter(self):
        # The fabricated corner coupon is a 3D metal slab whose every edge is a fold between
        # non-coplanar metal faces: the version-2 automatic perimeter extraction retains no
        # one-sided metal edge on it ("No physical metal perimeter was found"). The coupon
        # therefore names its edge lines explicitly as the perimeter of the SA surface minus
        # the matching box, for the thin sheet and the fabricated slab alike.
        for fabricated in (False, True):
            config = CORNER.make_config(
                Path("."),
                "corner",
                Path("mesh.msh"),
                [Path("trace.csv")],
                1.9,
                2,
                fabricated,
                11.47,
                {"SA": (0.002, 4.0), "MS": (0.002, 11.47), "MA": (0.002, 10.0)},
            )
            dielectrics = config["Boundaries"]["Postprocessing"]["Dielectric"]
            self.assertEqual([entry["Type"] for entry in dielectrics], ["SA", "MS", "MA"])
            for entry in dielectrics:
                self.assertNotIn("AutomaticEdges", entry)
                self.assertEqual(entry["EdgeAttributes"], [CORNER.SA_ATTRIBUTE])
                self.assertEqual(
                    entry["EdgeExcludeAttributes"], [CORNER.MATCHING_SURFACE_ATTRIBUTE]
                )
                self.assertEqual(entry["EdgeFrameNormal"], [0.0, 0.0, 1.0])
                self.assertEqual(entry["EdgeDistances"], [1.9])
            self.assertEqual(dielectrics[0]["Attributes"], [CORNER.SA_ATTRIBUTE])
            self.assertEqual(
                config["Boundaries"]["Ground"]["Attributes"], [2, 4] if fabricated else [2]
            )

    def test_corner_trace_mask_and_reference_follow_requested_angle(self):
        center = CORNER.corner_center(60.0, 0.5)
        np.testing.assert_allclose(
            center, [0.5 / np.tan(np.deg2rad(30.0)), 0.5]
        )
        points = np.asarray(
            [
                [1.0, 0.2, 0.0],
                [0.2, 1.0, 0.0],
                [0.0, 0.0, 0.0],
                [*center, 0.0],
            ]
        )
        convex = CORNER.metal_footprint_mask(
            points, 2.0, 60.0, 0.5, 0.0, "convex"
        )
        concave = CORNER.metal_footprint_mask(
            points, 2.0, 60.0, 0.5, 0.0, "concave"
        )
        np.testing.assert_array_equal(convex, [True, False, False, True])
        np.testing.assert_array_equal(concave, [False, True, True, False])

        with tempfile.TemporaryDirectory() as directory:
            library_path = CORNER.write_library(
                Path(directory),
                2.0,
                60.0,
                0.5,
                [8],
                [],
                "convex",
                0.1,
                0.05,
                80.0,
                0.01,
                0.01,
                11.47,
                {
                    "SA": (0.002, 4.0),
                    "MS": (0.0003, 11.47),
                    "MA": (0.03, 10.0),
                },
            )
            model = PREPARE.load_json(library_path)["Models"][0]
            self.assertEqual(model["Angle"], 60.0)
            np.testing.assert_allclose(model["Reference"], [*center, 0.0])

    def test_corner_builder_rejects_fillet_outside_matching_box(self):
        coupon = {
            "Id": "acute-rounded-corner",
            "Topology": "ConvexCorner",
            "Geometry": {"AngleDegrees": 10.0, "CornerRadius": 0.5},
            "BoundaryCondition": pec(),
        }
        args = SimpleNamespace(
            matching_radius=2.0,
            corner_lc_fine=0.02,
            min_process_feature_elements=2.0,
        )
        with self.assertRaisesRegex(
            ValueError, "tangency distance .* must be smaller"
        ):
            PREPARE.build_corner(
                coupon, args, process_parameters(), Path("unused")
            )

    def test_corner_builder_rejects_nonfinite_radius(self):
        coupon = {
            "Id": "invalid-rounded-corner",
            "Topology": "ConvexCorner",
            "Geometry": {
                "AngleDegrees": 90.0,
                "CornerRadius": float("nan"),
            },
            "BoundaryCondition": pec(),
        }
        args = SimpleNamespace(
            matching_radius=2.0,
            corner_lc_fine=0.02,
            min_process_feature_elements=2.0,
        )
        with self.assertRaisesRegex(
            ValueError, "corner radius must be finite"
        ):
            PREPARE.build_corner(
                coupon, args, process_parameters(), Path("unused")
            )

    def test_corner_refined_trace_basis_default_on_sharp_and_rounded_corners(self):
        # The planner's refined default (all-rings-follow-metal) applies to EVERY corner,
        # sharp or rounded (decision 511, fillet-basis design 2026-10-07): the box trace basis
        # is independent of CornerRadius (the fillet lies inside the matching box; a rounded
        # corner's basis files are byte-identical to the sharp corner's at the same angle),
        # and the legacy MetalRingsOnly rule it replaces fails the held-out self-check on
        # every 90-degree corner (S7's rounded corners, decision 310). Until decision 511 this
        # test pinned the opposite policy (a rounded corner forced onto the legacy rule, a
        # refined request on it refused): that guard was a qualification scope, not a
        # geometric constraint, and is lifted. An explicit legacy request stays honoured on
        # either kind of corner.
        def generator_command(coupon, args):
            calls = []

            def record(command, check=True):
                calls.append([str(value) for value in command])
                return 0

            with tempfile.TemporaryDirectory() as directory:
                cache = Path(directory)
                with (
                    mock.patch.object(PREPARE, "run", side_effect=record),
                    mock.patch.object(
                        PREPARE,
                        "run_probe_convergence",
                        return_value=(1, cache / "p.json", {"Passed": False}),
                    ),
                    self.assertRaisesRegex(RuntimeError, "probe convergence failed"),
                ):
                    PREPARE.build_corner(coupon, args, process_parameters(), cache)
                spec = PREPARE.load_json(next(cache.glob("corner-*/coupon-spec.json")))
            generator = next(
                command
                for command in calls
                if command[1].endswith("generate_corner_response.py")
            )
            return (
                generator[generator.index("--trace-basis") + 1],
                generator[generator.index("--ring-size") + 1],
                spec["Response"],
            )

        def corner(identifier, radius):
            return {
                "Id": identifier,
                "Topology": "ConvexCorner",
                "Geometry": {"AngleDegrees": 90.0, "CornerRadius": radius},
                "BoundaryCondition": pec(),
            }

        def arguments(**overrides):
            return SimpleNamespace(
                matching_radius=2.0,
                orders=[2, 3],
                corner_lc_fine=0.02,
                corner_lc_far=0.3,
                mesh_order=1,
                min_process_feature_elements=2.0,
                force=True,
                julia="julia",
                julia_project=None,
                **overrides,
            )

        resolvability = {
            "MeshSizing": "KnotGap",
            "MinimumActiveNodesPerOrderSquared": 5,
        }
        sharp = generator_command(corner("sharp-90", 0.0), arguments())
        self.assertEqual(sharp[:2], ("all-rings-follow-metal", "16"))
        self.assertEqual(
            sharp[2],
            {
                "RingSize": 16,
                "TraceBasis": "all-rings-follow-metal",
                "TraceResolvability": resolvability,
            },
        )
        rounded = generator_command(corner("rounded-90", 0.5), arguments())
        self.assertEqual(rounded[:2], ("all-rings-follow-metal", "16"))
        self.assertEqual(rounded[2], sharp[2])
        # The CLI default (None) resolves the same way; an explicit legacy request is honoured
        # on a sharp and on a rounded corner; an explicit refined request on a rounded corner
        # is the default.
        self.assertEqual(
            generator_command(
                corner("rounded-90", 0.5), arguments(corner_trace_basis=None)
            )[0],
            "all-rings-follow-metal",
        )
        for identifier, radius in (("sharp-90", 0.0), ("rounded-90", 0.5)):
            legacy = generator_command(
                corner(identifier, radius), arguments(corner_trace_basis="legacy")
            )
            self.assertEqual(legacy[:2], ("legacy", "8"))
            self.assertEqual(
                legacy[2],
                {"RingSize": 8, "TraceBasis": "legacy", "TraceResolvability": resolvability},
            )
        self.assertEqual(
            generator_command(
                corner("rounded-90", 0.5),
                arguments(corner_trace_basis="all-rings-follow-metal"),
            )[:2],
            ("all-rings-follow-metal", "16"),
        )
        with self.assertRaisesRegex(ValueError, "unknown corner trace basis rule"):
            PREPARE.corner_trace_basis(
                "rounded-90", SimpleNamespace(corner_trace_basis="fillet-knots")
            )

    def test_corner_mesh_refinement_scales_only_fine_size(self):
        args = SimpleNamespace(
            force=True,
            julia="julia",
            julia_project=None,
            matching_radius=2.0,
            orders=[3, 4],
            corner_lc_fine=0.02,
            corner_lc_far=0.3,
        )
        spec = {"Mesh": {"Order": 2}}
        commands = []

        def record(command, check=True):
            commands.append([str(value) for value in command])
            return 0

        with (
            tempfile.TemporaryDirectory() as directory,
            mock.patch.object(PREPARE, "run", side_effect=record),
        ):
            mesh_root, meshes = PREPARE.generate_corner_meshes(
                Path(directory),
                "convex",
                45.0,
                0.5,
                args,
                process_parameters(),
                spec,
                2.0,
            )

        self.assertEqual(mesh_root.name, "h-2")
        self.assertEqual(set(meshes), {"thin", "fabricated"})
        # Two mesher runs sized at the generator's trace mesh (decision 328), then the trace
        # resolvability gate on both meshes at the final solve order.
        self.assertEqual(len(commands), 3)
        for command in commands[:2]:
            self.assertEqual(command[command.index("--angle") + 1], "45.0")
            self.assertEqual(command[command.index("--lc-fine") + 1], "0.04")
            self.assertEqual(command[command.index("--lc-far") + 1], "0.3")
            self.assertEqual(command[command.index("--mesh-order") + 1], "2")
            self.assertEqual(command[command.index("--trace-mesh") + 1], directory)
        gate = commands[2]
        self.assertTrue(gate[1].endswith("trace_resolvability.py"))
        self.assertEqual(gate[2], directory)
        self.assertEqual(gate[gate.index("--order") + 1], "4")
        self.assertEqual(gate[gate.index("--radius") + 1], "2.0")
        self.assertEqual(gate.count("--mesh"), 2)
        self.assertEqual(
            gate[gate.index("--report") + 1], str(mesh_root / "trace-resolvability.json")
        )

    def test_plan_accepts_verified_finite_impedance(self):
        requirements = [
            {
                "Topology": "IsolatedEdge",
                "Geometry": {},
                "Interfaces": interfaces(),
                "BoundaryCondition": {"Type": "Impedance", "Ls": 1.0e-12},
                "Status": "Missing",
            },
            endpoint_coupon(),
        ]
        requirements[1]["Geometry"]["Arms"][0]["BoundaryCondition"] = {
            "Type": "Impedance",
            "Ls": 1.0e-12,
        }
        requirements[1]["Status"] = "Missing"
        manifest = {
            "Library": {"MatchingRadius": 2.0},
            "Requirements": requirements,
        }
        plan = PREPARE.plan_from_manifest(
            Path("requirements.json"),
            manifest,
            Path("library.json"),
            {"Fabrication": {}},
            False,
        )
        self.assertEqual(plan["Summary"]["Unsupported"], 0)
        methods = {
            coupon["Topology"]: coupon["Preparation"]["Method"]
            for coupon in plan["Coupons"]
        }
        self.assertEqual(methods["IsolatedEdge"], "StraightEdgeBuilder")
        self.assertEqual(methods["Endpoint"], "SpatialCoupon")
        for coupon in plan["Coupons"]:
            self.assertEqual(
                coupon["Preparation"]["BoundaryLawQualification"], "Missing"
            )

    def test_plan_rejects_unverified_finite_impedance(self):
        for condition in (
            "Impedance",
            {
                "Type": "Impedance",
                "Parameters": [0.0, 1.0e-12, 0.0],
                "ParametersVerified": True,
            },
            {"Type": "Impedance"},
        ):
            requirement = {
                "Topology": "IsolatedEdge",
                "Geometry": {},
                "Interfaces": interfaces(),
                "BoundaryCondition": condition,
                "Status": "Missing",
            }
            preparation = PREPARE.preparation(requirement)
            self.assertEqual(preparation["Method"], "Unsupported")
            self.assertIn("boundary-law", preparation["Reason"])

    def test_stamp_library_preserves_finite_impedance_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "process-library.json"
            path.write_text(
                json.dumps(
                    {
                        "Version": 3,
                        "Models": [
                            {
                                "Name": "isolated",
                                "Topology": "IsolatedEdge",
                            },
                            {
                                "Name": "spatial",
                                "Topology": "SpatialEdgeCluster",
                                "Edges": [{}, {}],
                            },
                        ],
                    }
                )
                + "\n"
            )
            isolated = {
                "Id": "isolated",
                "Topology": "IsolatedEdge",
                "Geometry": {},
                "BoundaryCondition": {
                    "Type": "Conductivity",
                    "Conductivity": 5.8e7,
                    "Permeability": 1.2,
                    "Thickness": 2.0e-7,
                    "External": False,
                },
            }
            spatial = overlapping_spatial_coupon()
            spatial["Id"] = "spatial"
            spatial["BoundaryCondition"] = {"Type": "PEC"}
            spatial["Geometry"]["Edges"][1]["BoundaryCondition"] = {
                "Type": "RationalImpedance",
                "Numerator": [1.0e-7, 0.0],
                "Denominator": [1.0e-19, 2.0e-9, 100.0],
            }
            PREPARE.stamp_library_boundary_conditions(
                path, [isolated, spatial]
            )
            models = PREPARE.load_json(path)["Models"]
            self.assertEqual(
                models[0]["BoundaryCondition"],
                isolated["BoundaryCondition"],
            )
            self.assertEqual(
                models[0]["BoundaryLawQualification"]["Status"],
                "Unqualified",
            )
            self.assertFalse(
                models[0]["BoundaryLawQualification"]["FrequencyUniversal"]
            )
            self.assertNotIn("BoundaryCondition", models[1])
            self.assertEqual(
                models[1]["Edges"][1]["BoundaryCondition"],
                spatial["Geometry"]["Edges"][1]["BoundaryCondition"],
            )
            self.assertEqual(
                models[1]["BoundaryLawQualification"]["Calibration"],
                "QuasiElectrostatic",
            )

    def test_content_ids_include_complete_geometry(self):
        first = endpoint_coupon()
        second = endpoint_coupon()
        second["Geometry"]["Arms"][0]["Interval"] = [0.0, 3.0]
        self.assertNotEqual(PREPARE.coupon_id(first), PREPARE.coupon_id(second))
        self.assertEqual(
            PREPARE.fingerprint({"a": 1, "b": 2}),
            PREPARE.fingerprint({"b": 2, "a": 1}),
        )

    def test_parallel_signature_is_canonical_and_one_based(self):
        coupon = {
            "Id": "cluster",
            "Geometry": {
                "Edges": [
                    {"Offset": [0.0, 0.0], "GapDirection": [1.0, 0.0], "Conductor": 0},
                    {"Offset": [1.0, 0.0], "GapDirection": [-1.0, 0.0], "Conductor": 0},
                    {"Offset": [2.0, 0.0], "GapDirection": [1.0, 0.0], "Conductor": 1},
                ]
            },
        }
        edges = PREPARE.normalize_parallel_edges(coupon, 2.0)
        self.assertEqual([edge["Conductor"] for edge in edges], [1, 1, 2])
        self.assertEqual([edge["Offset"] for edge in edges], [0.0, 1.0, 2.0])

    def test_parallel_stack_is_bounded_per_consecutive_link(self):
        """Decision 82(2): a stack's consecutive edges interact below 2R; its span may
        exceed 2R (the 4-edge 2 / 2 / 2 um flux line at R 1.9 spans 3.16 R), a link at or
        beyond 2R is refused."""
        def coupon(offsets):
            return {"Id": "stack", "Geometry": {"Edges": [
                {"Offset": [x, 0.0], "GapDirection": [1.0 if k % 2 == 0 else -1.0, 0.0], "Conductor": (1, 2, 2, 1)[k]}
                for k, x in enumerate(offsets)]}}
        R = 1.9
        edges = PREPARE.normalize_parallel_edges(coupon([0.0, 2.0, 4.0, 6.0]), R)
        self.assertEqual([edge["Offset"] for edge in edges], [0.0, 2.0, 4.0, 6.0])
        with self.assertRaises(ValueError):
            PREPARE.normalize_parallel_edges(coupon([0.0, 2.0, 4.0, 4.0 + 2.0 * R + 1e-6]), R)

    def test_version2_signatures_are_carried_and_stamped(self):
        """A version-2 record's Signature travels with the plan coupon and is stamped on the
        written model that corresponds to it (pairs by Separation, stacks by their Edges,
        corners by Angle / CornerRadius); a model matching two signatures fails closed."""
        signature = {"Type": "SameConductorStrip", "SeparationOverR": 1.0526316,
                     "Edges": [{"OffsetOverR": 0.0, "GapSide": -1, "Conductor": 1}, {"OffsetOverR": 1.0526316, "GapSide": 1, "Conductor": 1}]}
        stack_signature = {"Type": "ParallelEdgeCluster", "Edges": [
            {"OffsetOverR": 0.0, "GapSide": 1, "Conductor": 1}, {"OffsetOverR": 1.0526316, "GapSide": -1, "Conductor": 2},
            {"OffsetOverR": 2.1052632, "GapSide": 1, "Conductor": 2}]}
        manifest = {
            "Version": 2,
            "Library": {"MatchingRadius": 1.9, "Path": "seed.json"},
            "Requirements": [
                {"Topology": "SameConductorStrip", "Status": "Missing", "Count": 171, "TotalEdgeLength": 614.43,
                 "Geometry": {"EdgeCount": 2, "Separation": 2.0000008}, "Interfaces": [], "BoundaryCondition": {"Type": "PEC"},
                 "Signature": signature, "Hash": "ab" * 32, "Instances": 8, "DistinctSignatures": 1, "ParameterSpread": 0.0,
                 "ExactParameters": True},
                {"Topology": "ParallelEdgeCluster", "Status": "Missing", "Count": 10, "TotalEdgeLength": 30.0,
                 "Geometry": {"EdgeCount": 3, "Edges": [{"Offset": [0.0, 0.0], "GapDirection": [1.0, 0.0], "Conductor": 1},
                                                        {"Offset": [2.0, 0.0], "GapDirection": [-1.0, 0.0], "Conductor": 2},
                                                        {"Offset": [4.0, 0.0], "GapDirection": [1.0, 0.0], "Conductor": 2}]},
                 "Interfaces": [], "BoundaryCondition": {"Type": "PEC"}, "Signature": stack_signature, "Hash": "cd" * 32,
                 "Instances": 2, "DistinctSignatures": 1, "ParameterSpread": 0.0, "ExactParameters": True,
                 "NearKeys": ["12" * 32], "SpanCapAllowance": {"Label": "cdcdcdcdcdcd", "SpanCapOverR": 22.0,
                                                               "Reason": "unit", "Approval": "decision 303",
                                                               "MatchedQuanta": 1.0}},
                {"Topology": "ConvexCorner", "Status": "Missing", "Count": 12, "TotalEdgeLength": 45.6,
                 "Geometry": {"AngleDegrees": 90.0, "CornerRadius": 0.0}, "Interfaces": [], "BoundaryCondition": {"Type": "PEC"},
                 "Signature": {"Type": "ConvexCorner", "AngleDegrees": 90.0, "CornerRadiusOverR": 0.0}, "Hash": "ef" * 32,
                 "Instances": 12, "DistinctSignatures": 1, "ParameterSpread": 0.0, "ExactParameters": True},
            ],
        }
        plan = PREPARE.plan_from_manifest(Path("m.json"), manifest, Path("seed.json"), {}, False)
        by_topology = {coupon["Topology"]: coupon for coupon in plan["Coupons"]}
        self.assertEqual(by_topology["SameConductorStrip"]["Signature"], signature)
        self.assertEqual(by_topology["SameConductorStrip"]["Instances"], 8)
        self.assertEqual(by_topology["ParallelEdgeCluster"]["Preparation"]["Method"], "ParallelClusterCoupon")
        # Block (b) steps 1-2: the near-match group keys and the resolved span-cap allowance
        # travel with the plan coupon (device_coupons reads the allowance as the generator's
        # --support-span-cap).
        self.assertEqual(by_topology["ParallelEdgeCluster"]["NearKeys"], ["12" * 32])
        self.assertEqual(by_topology["ParallelEdgeCluster"]["SpanCapAllowance"]["SpanCapOverR"], 22.0)
        self.assertNotIn("SpanCapAllowance", by_topology["SameConductorStrip"])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "process-library.json"
            models = [
                {"Name": "isolated", "Topology": "IsolatedEdge"},
                {"Name": "strip", "Topology": "SameConductorStrip", "Separation": 2.0000008},
                {"Name": "other-strip", "Topology": "SameConductorStrip", "Separation": 3.0},
                {"Name": "stack", "Topology": "ParallelEdgeCluster",
                 "Edges": [{"Offset": 0.0, "GapDirection": 1, "Conductor": 1}, {"Offset": 2.0, "GapDirection": -1, "Conductor": 2},
                           {"Offset": 4.0, "GapDirection": 1, "Conductor": 2}]},
                {"Name": "corner", "Topology": "ConvexCorner", "Angle": 90.0, "CornerRadius": 0.0},
            ]
            PREPARE.write_json(path, {"Version": 3, "Models": models})
            self.assertEqual(PREPARE.stamp_library_signatures(path, plan["Coupons"]), 3)
            stamped = {model["Name"]: model for model in PREPARE.load_json(path)["Models"]}
            self.assertEqual(stamped["strip"]["Signature"], signature)
            self.assertEqual(stamped["strip"]["Instances"], 8)
            self.assertEqual(stamped["stack"]["Signature"], stack_signature)
            self.assertEqual(stamped["corner"]["Signature"]["Type"], "ConvexCorner")
            self.assertNotIn("Signature", stamped["isolated"])
            self.assertNotIn("Signature", stamped["other-strip"])
            duplicate = copy.deepcopy(by_topology["SameConductorStrip"])
            duplicate["Signature"] = {**signature, "SeparationOverR": 1.0526317}
            duplicate["BoundaryCondition"] = {"Type": "Impedance", "Ls": 1.0e-13}
            # The message names every matching coupon with its boundary law (the geometry
            # match ignores the law; a second law on one geometry is the usual cause).
            with self.assertRaisesRegex(ValueError, r"matches 2 version-2 coupon signatures.*Impedance"):
                PREPARE.stamp_library_signatures(path, plan["Coupons"] + [duplicate])

    def test_parallel_cluster_mesh_failure_prevents_full_response_solves(self):
        coupon = {
            "Id": "cluster",
            "Topology": "ParallelEdgeCluster",
            "Geometry": {
                "Edges": [
                    {"Offset": 0.0, "GapDirection": 1, "Conductor": 1},
                    {"Offset": 1.0, "GapDirection": -1, "Conductor": 1},
                    {"Offset": 2.0, "GapDirection": 1, "Conductor": 2},
                ]
            },
            "BoundaryCondition": pec(),
        }
        args = SimpleNamespace(
            matching_radius=2.0,
            orders=[2, 3],
            cluster_lc_fine=0.002,
            cluster_lc_far=0.05,
            cluster_h_factors=[2.0, 1.0],
            mesh_order=1,
            basis_size=16,
            samples=32,
            coupon_depth=10.0,
            edge_offset_tolerance=1.0e-3,
            min_process_feature_elements=2.0,
            force=True,
            name="test",
            palace=Path("palace"),
            julia="julia",
            julia_project=None,
            ranks=1,
            max_fabricated_matrix_change=5.0,
            max_fabricated_energy_change=10.0,
            max_domain_defect_change=5.0,
            max_heldout_error=10.0,
        )
        with tempfile.TemporaryDirectory() as directory:
            cache = Path(directory)
            p_report = cache / "p-convergence.json"
            h_report = cache / "h-convergence.json"
            calls = []

            def record(command, check=True):
                calls.append([str(value) for value in command])
                return 0

            with (
                mock.patch.object(PREPARE, "run", side_effect=record),
                mock.patch.object(
                    PREPARE,
                    "run_probe_convergence",
                    return_value=(0, p_report, {"Passed": True}),
                ),
                mock.patch.object(
                    PREPARE,
                    "prepare_probe_mesh_calibration",
                    return_value=cache / "coarse-calibration",
                ),
                mock.patch.object(
                    PREPARE,
                    "run_mesh_convergence",
                    return_value=(1, h_report, {"Passed": False}),
                ),
            ):
                with self.assertRaisesRegex(
                    RuntimeError, "Parallel-cluster mesh convergence failed"
                ):
                    PREPARE.build_parallel_cluster(
                        coupon, args, process_parameters(), cache
                    )

            self.assertTrue(
                any(
                    "--mesh-order" in command
                    and command[command.index("--mesh-order") + 1] == "2"
                    for command in calls
                )
            )
            self.assertFalse(any(command[0] == "palace" for command in calls))
            qualification = next(
                cache.glob("parallel-cluster-*/qualification.json")
            )
            result = PREPARE.load_json(qualification)
            self.assertEqual(result["MeshConvergenceReport"], str(h_report))
            self.assertFalse(result["Passed"])

    def test_spatial_matching_surface_encloses_process_geometry(self):
        coupon = endpoint_coupon()
        frame, edges, facets = SPATIAL.normalize_geometry(coupon, 2.0)
        self.assertEqual(facets, [])
        lower, upper = SPATIAL.coupon_bounds(edges, 2.0, 0.1, 0.05)
        self.assertAlmostEqual(lower[0], -4.0)
        self.assertAlmostEqual(upper[0], 4.0)
        levels = [lower[2], -0.05, 0.0, 0.1, upper[2]]
        points, triangles, groups = SPATIAL.build_matching_surface(
            np.vstack((lower, upper)),
            levels,
            8,
            edges,
            2.0,
            0.1,
            80.0,
            facets,
        )
        labels = SPATIAL.conductor_at_points(
            points, edges, 2.0, 0.1, 80.0, facets
        )
        self.assertGreaterEqual(groups[0], 8)
        self.assertEqual(set(np.unique(labels)), {0, 1})
        self.assertGreater(np.count_nonzero(labels == 0), 0)
        self.assertEqual(sum(groups), len(points))
        areas = np.linalg.norm(
            np.cross(
                points[triangles[:, 1]] - points[triangles[:, 0]],
                points[triangles[:, 2]] - points[triangles[:, 0]],
            ),
            axis=1,
        )
        self.assertTrue(np.all(areas > 0.0))
        canonical = SPATIAL.canonical_points(points, frame)
        self.assertEqual(canonical.shape, points.shape)

    def test_spatial_conductor_reference_is_order_invariant_and_inside_mask(self):
        edges = [
            {
                "Point": [0.0, 0.0, 0.0],
                "GapDirection": [-1.0, 0.0, 0.0],
                "Interval": [-1.0, 1.0],
                "Conductor": 1,
            },
            {
                "Point": [2.0, 0.0, 0.0],
                "GapDirection": [1.0, 0.0, 0.0],
                "Interval": [-1.0, 1.0],
                "Conductor": 1,
            },
        ]
        facets = [
            {
                "Conductor": 1,
                "Plane": 0.0,
                "Points": [[0.0, -1.0], [2.0, -1.0], [2.0, 1.0], [0.0, 1.0]],
            }
        ]
        first = SPATIAL.reference_points({}, edges, facets, np.eye(3), 2.0)
        second = SPATIAL.reference_points({}, list(reversed(edges)), facets, np.eye(3), 2.0)
        self.assertEqual(first, second)
        self.assertGreater(first[0][0], 0.0)
        self.assertLess(first[0][0], 2.0)

    def test_spatial_mask_labels_continuation_metal_outside_edge_interval(self):
        points = np.asarray([[-3.0, 3.0, 0.0]])
        edges = [
            {
                "Point": [0.0, 0.0, 0.0],
                "GapDirection": [1.0, 0.0, 0.0],
                "Tangent": [0.0, 1.0, 0.0],
                "ProcessNormal": [0.0, 0.0, 1.0],
                "Interval": [-0.5, 0.5],
                "Conductor": 1,
            }
        ]
        facets = [
            {
                "Conductor": 1,
                "Plane": 0.0,
                "Points": [[-4.0, -4.0], [0.0, -4.0], [0.0, 4.0], [-4.0, 4.0]],
            }
        ]
        labels = SPATIAL.conductor_at_points(points, edges, 2.0, 0.1, 90.0, facets)
        np.testing.assert_array_equal(labels, [1])

    def test_spatial_conductor_lift_is_one_only_on_excluded_contact_nodes(self):
        points = np.asarray(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
            ]
        )
        labels = np.asarray([1, 0, 2])
        lifts = SPATIAL.conductor_trace_lifts(points, labels, 2)
        np.testing.assert_array_equal(lifts[2], [0.0, 0.0, 1.0])

    def test_spatial_probe_cutoff_is_compatible_with_thin_and_fabricated_cuts(self):
        points = np.asarray(
            [
                [0.0, 0.0, -1.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.05],
                [0.0, 0.0, 0.1],
                [0.0, 0.0, 1.0],
            ]
        )
        edges = [
            {
                "Point": [0.0, 0.0, 0.0],
                "ProcessNormal": [0.0, 0.0, 1.0],
            }
        ]
        cutoff = SPATIAL.spatial_metal_band_cutoff(points, edges, 2.0, 0.1)
        self.assertEqual(cutoff[1], 0.0)
        self.assertEqual(cutoff[2], 0.0)
        self.assertEqual(cutoff[3], 0.0)
        self.assertGreater(cutoff[0], 0.0)
        self.assertGreater(cutoff[4], 0.0)

    def test_spatial_heldout_cutoff_vanishes_on_pec_knots_only(self):
        """USER decision 149 (6), option (c): the held-out trace's cutoff is the smoothstep
        over R / 3 of the distance to the nearest PEC contact knot, so a free knot of the
        metal-band rings (here the ring knots at z = 0 and z = 0.1 over the gap, which the
        band cutoff zeroed) is excited while every PEC knot stays at zero."""
        points = np.asarray(
            [
                [-2.0, 0.0, 0.0],  # PEC contact knot on the z = 0 ring
                [-2.0, 0.0, 0.1],  # PEC contact knot on the metal-top ring
                [-2.0, 1.0, 0.0],  # free knot of the z = 0 ring, 1.0 from the contact
                [-2.0, 1.0, 0.1],
                [-2.0, 2.0, 0.0],  # a box corner of the z = 0 ring, 2.0 away
                [-2.0, 0.0, 2.0],  # the top ring
                [-2.0, 0.0, -0.05],  # the trench ring right below the contact
            ]
        )
        labels = np.asarray([1, 1, 0, 0, 0, 0, 0])
        edges = [{"Point": [0.0, 0.0, 0.0], "ProcessNormal": [0.0, 0.0, 1.0]}]
        band = SPATIAL.spatial_metal_band_cutoff(points, edges, 2.0, 0.1)
        np.testing.assert_array_equal(band[:5], 0.0)
        cutoff = SPATIAL.spatial_heldout_cutoff(points, labels, 2.0)
        np.testing.assert_array_equal(cutoff[:2], 0.0)
        self.assertTrue(np.all(cutoff[2:] > 0.0))
        np.testing.assert_allclose(cutoff[2:6], 1.0)
        coordinate = 0.05 / (2.0 / 3.0)
        self.assertAlmostEqual(cutoff[6], coordinate * coordinate * (3.0 - 2.0 * coordinate))
        np.testing.assert_array_equal(
            SPATIAL.spatial_heldout_cutoff(points, np.zeros(len(points), dtype=int), 2.0), 1.0
        )

    def test_spatial_strip_extension_stops_at_finite_intervals(self):
        vertex_arm = {"Interval": [0.0, 2.0], "VertexArm": True}
        continuing = {"Interval": [-2.0, 2.0], "VertexArm": False}
        finite = {"Interval": [-1.0, 1.0], "VertexArm": False}
        self.assertEqual(SPATIAL.extended_interval(vertex_arm, 2.0), (0.0, 6.0))
        self.assertEqual(SPATIAL.extended_interval(continuing, 2.0), (-6.0, 6.0))
        self.assertEqual(SPATIAL.extended_interval(finite, 2.0), (-1.0, 1.0))

    def test_spatial_coupon_rejects_overlapping_conductor_half_strips(self):
        coupon = overlapping_spatial_coupon()
        _, edges, facets = SPATIAL.normalize_geometry(coupon, 2.0)
        with self.assertRaisesRegex(
            ValueError, "overlapping plan-view metal"
        ):
            SPATIAL.validate_plan_view_geometry(edges, 2.0, facets)

    def test_spatial_coupon_accepts_explicit_plan_view_mask(self):
        coupon = overlapping_spatial_coupon(masked=True)
        _, edges, facets = SPATIAL.normalize_geometry(coupon, 2.0)
        SPATIAL.validate_plan_view_geometry(edges, 2.0, facets)
        self.assertEqual({facet["Conductor"] for facet in facets}, {1, 2})

    def test_spatial_coupon_preserves_finite_impedance_law(self):
        coupon = endpoint_coupon()
        condition = {"Type": "Impedance", "Rs": 0.01, "Ls": 1.0e-12}
        coupon["BoundaryCondition"] = condition
        coupon["Geometry"]["Arms"][0]["BoundaryCondition"] = condition
        _, edges, _ = SPATIAL.normalize_geometry(coupon, 2.0)
        self.assertEqual(edges[0]["BoundaryCondition"], condition)

    def test_spatial_planner_fails_closed_without_mask_and_passes_mask_to_mesher(self):
        requirements = []
        for masked in (False, True):
            requirement = overlapping_spatial_coupon(masked)
            requirement["Status"] = "Missing"
            requirements.append(requirement)
        plan = PREPARE.plan_from_manifest(
            Path("requirements.json"),
            {
                "Library": {"MatchingRadius": 2.0},
                "Requirements": requirements,
            },
            Path("library.json"),
            {"Fabrication": {}},
            False,
        )
        self.assertEqual(
            [coupon["Preparation"]["Method"] for coupon in plan["Coupons"]],
            ["SpatialCoupon", "SpatialCoupon"],
        )
        unmasked = next(
            coupon
            for coupon in plan["Coupons"]
            if "PlanViewFacets" not in coupon["Geometry"]
        )
        masked = next(
            coupon
            for coupon in plan["Coupons"]
            if "PlanViewFacets" in coupon["Geometry"]
        )
        args = SimpleNamespace(
            matching_radius=2.0,
            orders=[1],
            spatial_lc_fine=0.02,
            spatial_lc_far=0.3,
            mesh_order=1,
            spatial_ring_size=8,
            min_process_feature_elements=2.0,
            force=False,
            palace=Path("palace"),
            julia="julia",
            julia_project=None,
            ranks=1,
            max_fabricated_matrix_change=5.0,
            max_fabricated_energy_change=10.0,
            max_domain_defect_change=5.0,
            max_heldout_error=10.0,
        )

        class MeshingReached(Exception):
            pass

        real_run = PREPARE.run
        meshing_commands = []

        def run_through_signature(command, check=True):
            if "--signature-only" in command:
                return real_run(command, check)
            meshing_commands.append([str(value) for value in command])
            raise MeshingReached

        with tempfile.TemporaryDirectory() as directory:
            cache = Path(directory)
            with mock.patch.object(PREPARE, "run", side_effect=run_through_signature):
                with self.assertRaisesRegex(
                    RuntimeError, "overlapping plan-view metal"
                ):
                    PREPARE.build_spatial(
                        unmasked, args, process_parameters(), cache
                    )
                with self.assertRaises(MeshingReached):
                    PREPARE.build_spatial(
                        masked, args, process_parameters(), cache
                    )
            self.assertEqual(len(meshing_commands), 1)
            self.assertIn("--mask", meshing_commands[0])
            self.assertEqual(
                meshing_commands[0][meshing_commands[0].index("--mesh-order") + 1],
                "2",
            )
            mask = Path(meshing_commands[0][meshing_commands[0].index("--mask") + 1])
            self.assertTrue(mask.is_file())
            generated_coupon = PREPARE.load_json(mask.parent / "coupon.json")
            self.assertIn(
                "PlanViewBoundary", generated_coupon["Geometry"]
            )
            self.assertEqual(
                generated_coupon["Geometry"]["PlanViewBoundary"],
                PREPARE.canonical_plan_view_boundary(
                    masked["Geometry"]["PlanViewFacets"], 2.0
                ),
            )

    def test_cross_layer_geometry_round_trips_through_canonical_frame(self):
        coupon = {
            "Topology": "SpatialEdgeCluster",
            "Geometry": {
                "Edges": [
                    {
                        "Point": [1.0, 2.0, 3.0],
                        "GapDirection": [0.0, 1.0, 0.0],
                        "ProcessNormal": [0.0, 0.0, 1.0],
                        "Interval": [-2.0, 1.0],
                        "Conductor": 1,
                        "InterfaceSlot": 0,
                        "BoundaryCondition": pec(),
                    },
                    {
                        "Point": [-1.0, 0.5, 8.0],
                        "GapDirection": [1.0, 0.0, 0.0],
                        "ProcessNormal": [0.0, 0.0, -1.0],
                        "Interval": [-0.5, 2.0],
                        "Conductor": 2,
                        "InterfaceSlot": 1,
                        "BoundaryCondition": pec(),
                    },
                ]
            },
            "Interfaces": interfaces(0) + interfaces(1),
            "BoundaryCondition": pec(),
        }
        frame, edges, facets = SPATIAL.normalize_geometry(coupon, 2.0)
        self.assertEqual(facets, [])
        original = np.asarray(
            [edge["Point"] for edge in coupon["Geometry"]["Edges"]]
        )
        local = np.asarray([edge["Point"] for edge in edges])
        np.testing.assert_allclose(
            SPATIAL.canonical_points(local, frame), original, atol=1.0e-12
        )
        np.testing.assert_allclose(frame @ frame.T, np.eye(3), atol=1.0e-12)
        self.assertEqual([edge["InterfaceSlot"] for edge in edges], [0, 1])

    def test_open_paths_partition_free_trace(self):
        labels = np.zeros(24, dtype=int)
        active = np.flatnonzero(labels == 0).tolist()
        paths = SPATIAL.open_paths([8, 8, 8], labels, active, 3)
        assigned = sorted(
            index for path in paths for index in path["Indices"]
        )
        self.assertEqual(assigned, list(range(1, len(active) + 1)))
        self.assertEqual(
            {
                conductor
                for path in paths
                for conductor in (
                    path["StartConductor"],
                    path["EndConductor"],
                )
            },
            {1, 2, 3},
        )
        self.assertTrue(
            all(
                path["StartConductor"] != path["EndConductor"]
                for path in paths
            )
        )

    def test_spatial_attributes_separate_conductors_from_interface_slots(self):
        edges = [
            {"Conductor": 1, "InterfaceSlot": 0},
            {"Conductor": 2, "InterfaceSlot": 0},
            {"Conductor": 2, "InterfaceSlot": 1},
        ]
        self.assertEqual(SPATIAL.edge_attributes(edges, False, 0), [3000])
        self.assertEqual(SPATIAL.edge_attributes(edges, True, 0), [3000, 3100])
        self.assertEqual(SPATIAL.edge_attributes(edges, True, 1), [3001, 3101])
        self.assertEqual(
            SPATIAL.edge_attributes(edges, True, 0, {3100}), [3100]
        )
        self.assertEqual(
            SPATIAL.edge_attributes(edges, False, 0, {3001}), [3001]
        )
        self.assertEqual(
            SPATIAL.conductor_attributes(edges, False, 2), [4002, 4102]
        )
        self.assertEqual(
            SPATIAL.conductor_attributes(edges, True, 2),
            [5002, 5102, 6002, 6102],
        )
        self.assertEqual(
            SPATIAL.interface_attributes(edges, False, 0, "MS"),
            [4001, 4002],
        )
        self.assertEqual(
            SPATIAL.interface_attributes(edges, True, 1, "MA"), [6102]
        )
        self.assertEqual(
            SPATIAL.interface_attributes(edges, True, 1, "MS"), [5102]
        )
        self.assertEqual(
            SPATIAL.interface_attributes(
                edges, True, 0, "SA", {3100, 5001, 6001}
            ),
            [3100],
        )
        with self.assertRaisesRegex(ValueError, "SA interface slot 1"):
            SPATIAL.interface_attributes(edges, True, 1, "SA", {3100})
        self.assertEqual(
            SPATIAL.interface_attributes(
                edges, True, 1, "SA", {3100}, allow_empty=True
            ),
            [],
        )

    def test_spatial_multislot_interface_attributes_are_disjoint(self):
        edges = [
            {"Conductor": 1, "InterfaceSlot": 0},
            {"Conductor": 1, "InterfaceSlot": 1},
        ]
        for fabricated in (False, True):
            for interface_type in ("MA", "MS"):
                first = set(
                    SPATIAL.interface_attributes(
                        edges, fabricated, 0, interface_type
                    )
                )
                second = set(
                    SPATIAL.interface_attributes(
                        edges, fabricated, 1, interface_type
                    )
                )
                self.assertTrue(first.isdisjoint(second))

    def test_spatial_multislot_legacy_mesh_fails_closed(self):
        edges = [
            {
                "Conductor": 1,
                "InterfaceSlot": slot,
                "ProcessNormal": [0.0, 1.0, 0.0],
            }
            for slot in (0, 1)
        ]
        with self.assertRaisesRegex(ValueError, "legacy conductor-wide"):
            SPATIAL.make_config(
                Path("."),
                "legacy",
                Path("mesh.msh"),
                [Path("trace.csv")],
                Path("zero.csv"),
                [],
                True,
                1,
                2.0,
                11.45,
                {"SA": (0.002, 4.0), "MS": (0.002, 11.45), "MA": (0.002, 10.0)},
                [{"Slot": 0, "Type": "MA"}, {"Slot": 1, "Type": "MA"}],
                edges,
                available_attributes={1, 3000, 3001, 5001, 6001},
            )
        # Explicit mesher metadata distinguishes a valid disappeared slot from a legacy
        # conductor-wide mesh whose physical names happen to look identical.
        SPATIAL.validate_metal_slot_partitioning(
            edges,
            True,
            {1, 3000, 3001, 5001, 6001},
            slot_partitioned=True,
        )

    def test_spatial_config_omits_physically_absent_interface(self):
        edges = [
            {
                "Conductor": 1,
                "InterfaceSlot": 0,
                "ProcessNormal": [0.0, 1.0, 0.0],
            },
            {
                "Conductor": 1,
                "InterfaceSlot": 1,
                "ProcessNormal": [0.0, 1.0, 0.0],
            },
        ]
        interfaces = [
            {"Slot": 0, "Type": "MS"},
            {"Slot": 1, "Type": "SA"},
        ]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            trace = root / "trace.csv"
            trace.write_text("0,0,0,0\n")
            config = SPATIAL.make_config(
                root,
                "test",
                root / "mesh.msh",
                [trace],
                trace,
                [],
                True,
                1,
                2.0,
                11.47,
                {
                    "SA": (0.002, 4.0),
                    "MS": (0.002, 11.47),
                    "MA": (0.002, 10.0),
                },
                interfaces,
                edges,
                available_attributes={1, 3101, 6001, 6101},
            )
        dielectric = config["Boundaries"]["Postprocessing"]["Dielectric"]
        self.assertEqual(config["Solver"]["Linear"]["EstimatorTol"], 5.0e-1)
        self.assertEqual(config["Solver"]["Linear"]["EstimatorMaxIts"], 5)
        self.assertEqual([entry["Index"] for entry in dielectric], [2])
        self.assertEqual(dielectric[0]["Attributes"], [3101])
        self.assertNotIn("EdgeExcludeAttributes", dielectric[0])

    def test_spatial_matching_support_is_explicit_and_bounded(self):
        frame = np.eye(3)
        points = SPATIAL.matching_support_points(
            np.asarray([-4.0, -3.0, -1.0]),
            np.asarray([4.0, 3.0, 1.0]),
            frame,
            1.0,
        )
        self.assertEqual(len(points), 8)
        with self.assertRaisesRegex(ValueError, "spans 16.20R .* above the cap of 16R"):
            SPATIAL.matching_support_points(
                np.asarray([-8.1, -3.0, -1.0]),
                np.asarray([8.1, 3.0, 1.0]),
                frame,
                1.0,
            )

    def test_spatial_matching_support_span_cap_is_raised_per_coupon_only(self):
        """A single closed feature wider than 16R (the S1p loop end, 18.24 R) builds with an
        explicit --support-span-cap; the default stays 16R, the cap must be positive and
        the span is still judged against the raised cap."""
        frame = np.eye(3)
        lower, upper = np.asarray([-9.12, -7.0, -1.0]), np.asarray([9.12, 7.0, 1.0])
        self.assertEqual(SPATIAL.DEFAULT_SUPPORT_SPAN_CAP_OVER_R, 16.0)
        with self.assertRaisesRegex(ValueError, "spans 18.24R .* above the cap of 16R"):
            SPATIAL.matching_support_points(lower, upper, frame, 1.0)
        self.assertEqual(len(SPATIAL.matching_support_points(lower, upper, frame, 1.0, 20.0)), 8)
        with self.assertRaisesRegex(ValueError, "above the cap of 18R"):
            SPATIAL.matching_support_points(lower, upper, frame, 1.0, 18.0)
        with self.assertRaisesRegex(ValueError, "must be positive"):
            SPATIAL.matching_support_points(lower, upper, frame, 1.0, 0.0)

    def test_spatial_library_uses_reference_for_single_conductor(self):
        coupon = {
            "Topology": "SpatialEdgeCluster",
            "BoundaryCondition": pec(),
            "Geometry": {
                "Edges": [
                    {
                        "Conductor": 1,
                        "InterfaceSlot": 0,
                        "Point": [0.0, 0.0, 0.0],
                        "GapDirection": [1.0, 0.0, 0.0],
                        "ProcessNormal": [0.0, 1.0, 0.0],
                        "Interval": [-1.0, 1.0],
                    }
                ]
            },
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            basis = root / "basis-points.csv"
            basis.write_text("x,y,z\n0,0,0\n")
            path = SPATIAL.write_library(
                root,
                coupon,
                2.0,
                basis,
                [[x, y, z] for x in (-2.0, 2.0) for y in (-2.0, 2.0) for z in (-2.0, 2.0)],
                [1],
                [],
                [],
                [[0.0, 0.0, 0.0]],
                [{"Slot": 0, "Type": "SA"}],
                process_parameters(),
                "single-conductor",
            )
            library = PREPARE.load_json(path)
            model = library["Models"][0]
        self.assertTrue(library["ExhaustiveSpatialClosure"])
        self.assertEqual(model["Reference"], [0.0, 0.0, 0.0])
        self.assertEqual(len(model["SupportPoints"]), 8)
        self.assertNotIn("ConductorReferences", model)

    def test_spatial_mesh_boundary_attributes_reads_named_surface_groups(self):
        contents = (
            b"$MeshFormat\n2.2 1 8\n"
            b"$EndMeshFormat\n$PhysicalNames\n"
            b"7\n3 1 \"substrate\"\n2 1 \"matching_surface\"\n"
            b"2 3100 \"surface_3100\"\n2 5001 \"surface_5001\"\n"
            b"2 4101 \"surface_4101\"\n2 5101 \"surface_5101\"\n"
            b"2 6101 \"surface_6101\"\n"
            b"$EndPhysicalNames\n$Nodes\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            mesh = Path(directory) / "coupon.msh"
            mesh.write_bytes(contents)
            self.assertEqual(
                SPATIAL.mesh_boundary_attributes(mesh),
                {1, 3100, 4101, 5001, 5101, 6101},
            )

    def test_spatial_cache_key_reuses_exact_coupon_and_invalidates_geometry(self):
        coupon = endpoint_coupon()
        parameters = process_parameters()
        args = SimpleNamespace(
            matching_radius=2.0,
            orders=[1, 2],
            spatial_lc_fine=0.02,
            spatial_lc_far=0.3,
            mesh_order=1,
            spatial_ring_size=8,
            min_process_feature_elements=2.0,
            force=False,
        )
        exact = PREPARE.spatial_spec(coupon, args, parameters)
        self.assertEqual(exact["Mesh"]["Order"], 2)
        self.assertEqual(exact["Mesh"]["HRefinementFactors"], [2.0, 1.0])
        changed = endpoint_coupon()
        changed["Geometry"]["Arms"][0]["Interval"] = [0.0, 1.5]
        self.assertNotEqual(
            PREPARE.fingerprint(exact),
            PREPARE.fingerprint(PREPARE.spatial_spec(changed, args, parameters)),
        )

        with tempfile.TemporaryDirectory() as directory:
            cache = Path(directory)
            key = PREPARE.fingerprint(exact)
            root = cache / f"spatial-{key}"
            root.mkdir()
            library = root / "process-library.json"
            qualification = root / "qualification.json"
            library.write_text('{"Version": 3, "Models": []}\n')
            qualification.write_text(
                '{"Fingerprint": "' + key + '", "Passed": true}\n'
            )
            with mock.patch.object(
                PREPARE, "run", side_effect=AssertionError("cache was rebuilt")
            ):
                reused_library, reused_qualification = PREPARE.build_spatial(
                    coupon, args, parameters, cache
                )
            self.assertEqual(reused_library, library)
            self.assertEqual(reused_qualification, qualification)

    def test_spatial_mesher_selection_is_explicit_and_fingerprinted(self):
        coupon = endpoint_coupon()
        parameters = process_parameters()
        common = {
            "matching_radius": 2.0,
            "orders": [1, 2],
            "spatial_lc_fine": 0.02,
            "spatial_lc_tangent": 0.1,
            "spatial_lc_far": 0.3,
            "mesh_order": 2,
            "spatial_ring_size": 8,
            "min_process_feature_elements": 2.0,
        }
        field_args = SimpleNamespace(**common, spatial_mesher="field")
        swept_args = SimpleNamespace(**common, spatial_mesher="swept")
        self.assertEqual(
            PREPARE.spatial_mesh_path(field_args),
            PREPARE.SPATIAL_FIELD_MESH,
        )
        self.assertEqual(
            PREPARE.spatial_mesh_path(swept_args),
            PREPARE.SPATIAL_SWEPT_MESH,
        )
        swept_preparation = PREPARE.preparation(
            coupon,
            PREPARE.SPATIAL_SWEPT_MESH,
        )
        self.assertTrue(
            swept_preparation["MeshGenerator"].endswith(
                "mesh_spatial_coupon_swept.jl"
            )
        )
        field = PREPARE.spatial_spec(coupon, field_args, parameters)
        swept = PREPARE.spatial_spec(coupon, swept_args, parameters)
        self.assertEqual(field["Mesh"]["Mesher"], "field")
        self.assertEqual(swept["Mesh"]["Mesher"], "swept")
        self.assertNotEqual(
            PREPARE.fingerprint(field),
            PREPARE.fingerprint(swept),
        )

    def test_swept_spatial_command_keeps_tangent_fixed_during_normal_refinement(self):
        args = SimpleNamespace(
            spatial_mesher="swept",
            force=True,
            julia="julia",
            julia_project=None,
            matching_radius=2.0,
            spatial_lc_fine=0.02,
            spatial_lc_tangent=0.1,
            spatial_lc_far=0.3,
            spatial_process_core_width=1.0,
            spatial_process_fine_width=0.0,
            spatial_process_grading_power=1.7,
            spatial_max_nodes=500_000,
            spatial_max_elements=2_000_000,
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            signature = root / "mesh-signature.csv"
            signature.write_text("signature\n")
            commands = []
            with mock.patch.object(
                PREPARE,
                "run",
                side_effect=lambda command: commands.append(
                    [str(value) for value in command]
                ),
            ):
                PREPARE.generate_spatial_meshes(
                    root,
                    signature,
                    root / "missing-mask.csv",
                    root / "missing-boundary.csv",
                    args,
                    process_parameters(),
                    {"Mesh": {"Order": 2}},
                    2.0,
                )
        self.assertEqual(len(commands), 2)
        for command in commands:
            self.assertEqual(command[1], str(PREPARE.SPATIAL_SWEPT_MESH))
            self.assertEqual(
                command[command.index("--lc-fine") + 1],
                "0.04",
            )
            self.assertEqual(
                command[command.index("--lc-tangent") + 1],
                "0.1",
            )

    def test_swept_spatial_parser_requires_positive_independent_tangent(self):
        with tempfile.TemporaryDirectory() as directory:
            common = [
                "prepare_surface_response_coupons.py",
                "manifest.json",
                "--output",
                directory,
                "--spatial-mesher",
                "swept",
            ]
            with mock.patch.object(sys, "argv", common):
                with self.assertRaises(SystemExit):
                    PREPARE.parse_args()
            with mock.patch.object(
                sys,
                "argv",
                common + ["--spatial-lc-tangent", "0.01"],
            ):
                args = PREPARE.parse_args()
        self.assertEqual(args.spatial_mesher, "swept")
        self.assertEqual(args.spatial_lc_tangent, 0.01)
        self.assertEqual(args.spatial_max_nodes, 3_000_000)

    def test_field_spatial_parser_checks_coarsest_normal_size(self):
        with tempfile.TemporaryDirectory() as directory:
            command = [
                "prepare_surface_response_coupons.py",
                "manifest.json",
                "--output",
                directory,
                "--spatial-mesher",
                "field",
                "--spatial-lc-fine",
                "0.02",
                "--spatial-lc-tangent",
                "0.02",
                "--spatial-h-factors",
                "2",
                "1",
            ]
            with mock.patch.object(sys, "argv", command):
                with self.assertRaises(SystemExit):
                    PREPARE.parse_args()

    def test_spatial_cache_key_uses_union_boundary_not_facet_triangulation(self):
        args = SimpleNamespace(
            matching_radius=2.0,
            orders=[1, 2],
            spatial_lc_fine=0.02,
            spatial_lc_far=0.3,
            mesh_order=1,
            spatial_ring_size=8,
            min_process_feature_elements=2.0,
        )
        corners = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
        coarse = endpoint_coupon()
        coarse["Geometry"]["PlanViewFacets"] = [
            {"Conductor": 1, "Points": [corners[0], corners[1], corners[2]]},
            {"Conductor": 1, "Points": [corners[0], corners[2], corners[3]]},
        ]
        refined = endpoint_coupon()
        center = [1.0, 0.5, 0.0]
        midpoints = [
            [1.0, 0.0, 0.0],
            [2.0, 0.5, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.5, 0.0],
        ]
        refined["Geometry"]["PlanViewFacets"] = [
            {
                "Conductor": 1,
                "Points": [
                    center,
                    corners[index],
                    midpoints[index],
                ],
            }
            for index in range(4)
        ] + [
            {
                "Conductor": 1,
                "Points": [
                    center,
                    midpoints[index],
                    corners[(index + 1) % 4],
                ],
            }
            for index in range(4)
        ]
        coarse_spec = PREPARE.spatial_spec(
            coarse, args, process_parameters()
        )
        refined_spec = PREPARE.spatial_spec(
            refined, args, process_parameters()
        )
        self.assertEqual(
            coarse_spec["Coupon"]["Geometry"]["PlanViewBoundary"],
            refined_spec["Coupon"]["Geometry"]["PlanViewBoundary"],
        )
        self.assertEqual(
            PREPARE.fingerprint(coarse_spec),
            PREPARE.fingerprint(refined_spec),
        )
        self.assertNotIn(
            "PlanViewFacets", coarse_spec["Coupon"]["Geometry"]
        )

    def test_spatial_plan_aggregates_equivalent_mask_triangulations(self):
        corners = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
        center = [1.0, 0.5, 0.0]
        coarse = endpoint_coupon()
        coarse["Geometry"]["PlanViewFacets"] = [
            {"Conductor": 1, "Points": [corners[0], corners[1], corners[2]]},
            {"Conductor": 1, "Points": [corners[0], corners[2], corners[3]]},
        ]
        refined = endpoint_coupon()
        refined["Geometry"]["PlanViewFacets"] = [
            {
                "Conductor": 1,
                "Points": [center, corners[index], corners[(index + 1) % 4]],
            }
            for index in range(4)
        ]
        boundary = PREPARE.canonical_plan_view_boundary(
            coarse["Geometry"]["PlanViewFacets"], 2.0, 2
        )
        coarse["Geometry"]["PlanViewBoundary"] = boundary
        refined["Geometry"]["PlanViewBoundary"] = boundary
        coarse.update({"Status": "Missing", "Count": 2, "TotalEdgeLength": 3.0})
        refined.update({"Status": "Missing", "Count": 5, "TotalEdgeLength": 7.0})

        plan = PREPARE.plan_from_manifest(
            Path("requirements.json"),
            {
                "Library": {"MatchingRadius": 2.0},
                "Requirements": [coarse, refined],
            },
            Path("library.json"),
            {"Fabrication": {}},
            False,
        )
        self.assertEqual(plan["Summary"]["CouponCount"], 1)
        self.assertEqual(plan["Coupons"][0]["DeviceOccurrences"], 7)
        self.assertEqual(plan["Coupons"][0]["DeviceEdgeLength"], 10.0)

    def test_spatial_cache_boundary_rejects_malformed_facets(self):
        with self.assertRaisesRegex(ValueError, "not on one process plane"):
            PREPARE.canonical_plan_view_boundary(
                [
                    {
                        "Conductor": 1,
                        "Points": [
                            [0.0, 0.0, 0.0],
                            [1.0, 0.0, 0.0],
                            [1.0, 0.0, 1.0],
                            [0.0, 0.1, 1.0],
                        ],
                    }
                ],
                2.0,
            )

        with self.assertRaisesRegex(ValueError, "nonmanifold"):
            PREPARE.canonical_plan_view_boundary(
                [
                    {
                        "Conductor": 1,
                        "Points": [
                            [0.0, 0.0, 0.0],
                            [1.0, 0.0, 0.0],
                            [0.0, 0.0, 1.0],
                        ],
                    },
                    {
                        "Conductor": 1,
                        "Points": [
                            [0.0, 0.0, 0.0],
                            [1.0, 0.0, 0.0],
                            [0.5, 0.0, -1.0],
                        ],
                    },
                    {
                        "Conductor": 1,
                        "Points": [
                            [0.0, 0.0, 0.0],
                            [1.0, 0.0, 0.0],
                            [1.0, 0.0, -1.0],
                        ],
                    },
                ],
                2.0,
            )

    def test_spatial_cache_boundary_uses_cpp_rounding_and_integer_topology(self):
        boundary = PREPARE.canonical_plan_view_boundary(
            [
                {
                    "Conductor": 1,
                    "Points": [
                        [-1.5, 0.0, 0.0],
                        [0.5, 0.0, 0.0],
                        [0.5, 0.0, 1.5],
                        [-1.5, 0.0, 1.5],
                    ],
                }
            ],
            1.0e9,
        )
        self.assertEqual(
            boundary,
            [
                {
                    "Conductor": 1,
                    "Segments": [
                        [[-2, 0, 0], [-2, 0, 2]],
                        [[-2, 0, 0], [1, 0, 0]],
                        [[-2, 0, 2], [1, 0, 2]],
                        [[1, 0, 0], [1, 0, 2]],
                    ],
                }
            ],
        )

    def test_spatial_cache_boundary_merges_near_collinear_joints(self):
        # Block (b) DESIGN A3 (1): a 1e-7-rad quantisation kink between two boundary sides
        # (the oblique claim / context junctions of the O1 / O4 coupons) is one side; a
        # 2e-4-rad kink stays a vertex (the smallest real kink of the census is 3.9e-4 rad).
        def facets(kink):
            return [
                {
                    "Conductor": 1,
                    "Points": [
                        [0.0, 0.0, 0.0],
                        [4.0, 0.0, 3.0],
                        [8.0, 0.0, 6.0 + 5.0 * kink],
                        [8.0, 0.0, 10.0],
                        [0.0, 0.0, 10.0],
                    ],
                }
            ]

        merged = PREPARE.canonical_plan_view_boundary(facets(1.0e-7), 1.0)
        self.assertEqual(len(merged[0]["Segments"]), 4)
        self.assertIn([[0, 0, 0], [8000000000, 0, 6000000500]], merged[0]["Segments"])
        kept = PREPARE.canonical_plan_view_boundary(facets(2.0e-4), 1.0)
        self.assertEqual(len(kept[0]["Segments"]), 5)
        exact = PREPARE.canonical_plan_view_boundary(facets(0.0), 1.0)
        self.assertEqual(len(exact[0]["Segments"]), 4)
        self.assertIn([[0, 0, 0], [8000000000, 0, 6000000000]], exact[0]["Segments"])

    def test_spatial_boundary_loops_merge_near_collinear_joints_of_one_class(self):
        # The mesher's plan-view boundary (generate_spatial_response.plan_view_boundary_loops):
        # the same rule with the vertex classes - a 1e-7-rad Physical / Physical kink merges, a
        # Physical / Continuation joint never does, a 2e-4-rad kink stays.
        def loops(kink):
            apex = [0.0, 10.0]
            chain = [[0.0, 0.0], [4.0, 3.0], [8.0, 6.0 + 5.0 * kink], [8.0, 10.0]]
            facets = [{"Conductor": 1, "Plane": 0.0, "Points": [apex, a, b]} for a, b in zip(chain, chain[1:])]
            return SPATIAL.plan_view_boundary_loops(
                facets, 1.0, np.asarray([-1.0, -1.0, -1.0]), np.asarray([8.0, 11.0, 1.0])
            )

        merged = loops(1.0e-7)
        self.assertEqual(len(merged), 1)
        self.assertEqual(len(merged[0]["Points"]), 4)
        self.assertEqual(merged[0]["Classes"].count("Continuation"), 1)
        kept = loops(2.0e-4)
        self.assertEqual(len(kept[0]["Points"]), 5)
        self.assertEqual(len(loops(0.0)[0]["Points"]), 4)
        self.assertTrue(SPATIAL.near_collinear((1, 0), (1000000, 50)))
        self.assertFalse(SPATIAL.near_collinear((1, 0), (1000000, 250)))
        self.assertFalse(SPATIAL.near_collinear((1, 0), (-1000000, 50)))  # a near-spike, not a joint

    def test_spatial_mask_facets_follow_the_merged_boundary(self):
        # Block (b) step 4.1: the merged near-collinear vertex is moved onto the merged side in
        # the mask facets too (the mesher tests its metal surfaces against the mask at 1e-7 R,
        # below the kink offset), so the mask and the boundary bound the same metal; an exactly
        # collinear vertex and a loop vertex are untouched; a 2e-4-rad kink (kept in the
        # boundary) moves nothing.
        def facets_and_loops(kink):
            apex = [0.0, 10.0]
            chain = [[0.0, 0.0], [4.0, 3.0], [8.0, 6.0 + 5.0 * kink], [8.0, 10.0]]
            facets = [{"Conductor": 1, "Plane": 0.0, "Points": [list(apex), list(a), list(b)]}
                      for a, b in zip(chain, chain[1:])]
            loops = SPATIAL.plan_view_boundary_loops(
                facets, 1.0, np.asarray([-1.0, -1.0, -1.0]), np.asarray([8.0, 11.0, 1.0])
            )
            return facets, loops

        facets, loops = facets_and_loops(1.0e-7)
        before = [list(map(tuple, facet["Points"])) for facet in facets]
        moved = SPATIAL.reconcile_mask_with_boundary(facets, loops, 1.0)
        self.assertEqual(len(moved), 2)                       # the dropped (4, 3) in two facets
        for origin, target in moved:
            self.assertEqual(origin, [4.0, 3.0])
            # Onto the merged side (0, 0) -> (8, 6 + 5e-7): the foot of the perpendicular.
            self.assertAlmostEqual(target[1] - 0.75 * target[0], 0.0, delta=1.0e-6)
            self.assertGreater(abs(target[1] - 3.0) + abs(target[0] - 4.0), 0.0)
        loop_vertices = {tuple(np.round(p, 9)) for p in loops[0]["Points"]}
        for facet in facets:
            for point in facet["Points"]:
                rounded = tuple(np.round(point, 9))
                if rounded not in loop_vertices:
                    self.assertAlmostEqual(point[1] - 0.75 * point[0], 0.0, delta=1.0e-6)
        facets, loops = facets_and_loops(0.0)                 # exactly collinear: byte-identical
        before = [[list(p) for p in facet["Points"]] for facet in facets]
        self.assertEqual(SPATIAL.reconcile_mask_with_boundary(facets, loops, 1.0), [])
        self.assertEqual([[list(p) for p in facet["Points"]] for facet in facets], before)
        facets, loops = facets_and_loops(2.0e-4)              # a kept kink: nothing to move
        self.assertEqual(SPATIAL.reconcile_mask_with_boundary(facets, loops, 1.0), [])

    def test_spatial_boundary_loops_classify_only_coupon_clipping_edges(self):
        facets = [
            {
                "Conductor": 1,
                "Plane": 0.0,
                "Points": [[0.0, 0.0], [2.0, 0.0], [2.0, 1.0]],
            },
            {
                "Conductor": 1,
                "Plane": 0.0,
                "Points": [[0.0, 0.0], [2.0, 1.0], [0.0, 1.0]],
            },
        ]
        loops = SPATIAL.plan_view_boundary_loops(
            facets,
            1.0,
            np.asarray([0.0, -1.0, -1.0]),
            np.asarray([3.0, 2.0, 1.0]),
        )
        self.assertEqual(len(loops), 1)
        self.assertFalse(loops[0]["Hole"])
        self.assertEqual(
            loops[0]["Classes"].count("Continuation"),
            1,
        )
        continuation = loops[0]["Classes"].index("Continuation")
        first = loops[0]["Points"][continuation]
        second = loops[0]["Points"][(continuation + 1) % len(loops[0]["Points"])]
        self.assertEqual(first[0], 0.0)
        self.assertEqual(second[0], 0.0)

    def test_spatial_boundary_loops_preserve_exported_classification(self):
        facets = [
            {
                "Conductor": 1,
                "Plane": 0.0,
                "Points": [[0.0, 0.0], [2.0, 0.0], [2.0, 1.0]],
            },
            {
                "Conductor": 1,
                "Plane": 0.0,
                "Points": [[0.0, 0.0], [2.0, 1.0], [0.0, 1.0]],
            },
        ]
        segments = [
            {
                "Conductor": 1,
                "Plane": 0.0,
                "Points": np.asarray([[2.0, 0.0], [2.0, 1.0]]),
            }
        ]
        loops = SPATIAL.plan_view_boundary_loops(
            facets,
            1.0,
            np.asarray([0.0, -1.0, -1.0]),
            np.asarray([3.0, 2.0, 1.0]),
            segments,
        )
        self.assertEqual(loops[0]["Classes"].count("Continuation"), 1)
        continuation = loops[0]["Classes"].index("Continuation")
        first = loops[0]["Points"][continuation]
        second = loops[0]["Points"][
            (continuation + 1) % len(loops[0]["Points"])
        ]
        self.assertEqual(first[0], 2.0)
        self.assertEqual(second[0], 2.0)

    def test_spatial_classified_boundary_transforms_to_mesher_frame(self):
        geometry = endpoint_coupon()["Geometry"]
        geometry["PlanViewBoundary"] = [
            {
                "Conductor": 1,
                "Segments": [
                    [[0, 0, 0], [1000000000, 0, 0]],
                ],
                "ContinuationSegments": [
                    [[0, 0, 0], [1000000000, 0, 0]],
                ],
            }
        ]
        frame = SPATIAL.frame_from_geometry("Endpoint", geometry)
        segments = SPATIAL.classified_continuation_segments(
            geometry, frame, 1.0
        )
        self.assertEqual(len(segments), 1)
        np.testing.assert_allclose(
            segments[0]["Points"],
            [[0.0, 0.0], [0.0, -1.0]],
            atol=1.0e-15,
        )

    def test_spatial_probe_failure_prevents_full_response_solves(self):
        args = SimpleNamespace(
            matching_radius=2.0,
            orders=[1, 2],
            spatial_lc_fine=0.02,
            spatial_lc_far=0.3,
            mesh_order=1,
            spatial_ring_size=8,
            min_process_feature_elements=2.0,
            force=False,
            palace=Path("palace"),
            julia="julia",
            julia_project=None,
            ranks=1,
            max_fabricated_matrix_change=5.0,
            max_fabricated_energy_change=10.0,
            max_domain_defect_change=5.0,
            max_heldout_error=10.0,
        )
        with tempfile.TemporaryDirectory() as directory:
            cache = Path(directory)
            report = cache / "probe-convergence.json"
            calls = []

            def record(command, check=True):
                calls.append([str(value) for value in command])
                return 0

            with (
                mock.patch.object(PREPARE, "run", side_effect=record),
                mock.patch.object(
                    PREPARE,
                    "run_probe_convergence",
                    return_value=(1, report, {"Passed": False}),
                ),
            ):
                with self.assertRaisesRegex(
                    RuntimeError, "Spatial probe convergence failed"
                ):
                    PREPARE.build_spatial(
                        endpoint_coupon(), args, process_parameters(), cache
                    )

            palace_configs = [
                argument
                for command in calls
                for argument in command
                if argument.endswith(".json")
            ]
            self.assertFalse(
                any(
                    Path(config).name
                    in {
                        "spatial_thin.json",
                        "spatial_fabricated.json",
                        "heldout_spatial_thin.json",
                        "heldout_spatial_fabricated.json",
                    }
                    for config in palace_configs
                )
            )

    def test_corner_trace_resolvability_gate_fails_closed_before_any_solve(self):
        # Decision 328: the generator runs BEFORE the mesher (the mesher sizes the mesh at the
        # trace knots, --trace-mesh = the cache root), the resolvability gate runs on both
        # meshes before any solve, and a failing gate stops the coupon with a recorded
        # qualification (Passed False, the report path, the reason) — no palace call.
        coupon = {
            "Id": "concave-48.75",
            "Topology": "ConcaveCorner",
            "Geometry": {"AngleDegrees": 48.75, "CornerRadius": 0.0},
            "BoundaryCondition": pec(),
        }
        args = SimpleNamespace(
            matching_radius=1.9,
            orders=[3, 4],
            corner_lc_fine=0.02,
            corner_lc_far=0.3,
            mesh_order=2,
            min_process_feature_elements=2.0,
            force=True,
            palace=Path("palace"),
            julia="julia",
            julia_project=None,
            ranks=1,
        )
        with tempfile.TemporaryDirectory() as directory:
            cache = Path(directory)
            calls = []

            def record(command, check=True):
                command = [str(value) for value in command]
                calls.append(command)
                return 1 if command[1].endswith("trace_resolvability.py") else 0

            with (
                mock.patch.object(PREPARE, "run", side_effect=record),
                self.assertRaisesRegex(RuntimeError, "do not resolve the trace basis"),
            ):
                PREPARE.build_corner(coupon, args, process_parameters(), cache)
            root = next(cache.glob("corner-*"))
            scripts = [Path(command[1]).name for command in calls]
            self.assertEqual(
                scripts,
                [
                    "generate_corner_response.py",
                    "mesh_corner_coupon.jl",
                    "mesh_corner_coupon.jl",
                    "trace_resolvability.py",
                ],
            )
            for command in calls[1:3]:
                self.assertEqual(command[command.index("--trace-mesh") + 1], str(root))
            generator = calls[0]
            self.assertEqual(
                generator[generator.index("--thin-mesh") + 1], str(root / "corner_thin.msh")
            )
            self.assertFalse(any(command[0] == "palace" for command in calls))
            qualification = PREPARE.load_json(root / "qualification.json")
            self.assertFalse(qualification["Passed"])
            self.assertEqual(
                qualification["TraceResolvabilityReport"],
                str(root / "trace-resolvability.json"),
            )
            self.assertIn("do not resolve the trace basis", qualification["Reason"])
            spec = PREPARE.load_json(root / "coupon-spec.json")
            self.assertEqual(
                spec["Response"]["TraceResolvability"],
                {"MeshSizing": "KnotGap", "MinimumActiveNodesPerOrderSquared": 5},
            )

    def test_corner_mesh_failure_prevents_full_response_solves(self):
        coupon = {
            "Id": "convex-45",
            "Topology": "ConvexCorner",
            "Geometry": {"AngleDegrees": 45.0, "CornerRadius": 0.5},
            "BoundaryCondition": pec(),
        }
        args = SimpleNamespace(
            matching_radius=2.0,
            orders=[2, 3],
            corner_lc_fine=0.02,
            corner_lc_far=0.3,
            corner_h_factors=[2.0, 1.0],
            mesh_order=1,
            ring_size=8,
            min_process_feature_elements=2.0,
            force=True,
            palace=Path("palace"),
            julia="julia",
            julia_project=None,
            ranks=1,
            max_fabricated_matrix_change=5.0,
            max_fabricated_energy_change=10.0,
            max_domain_defect_change=5.0,
            max_heldout_error=10.0,
        )
        with tempfile.TemporaryDirectory() as directory:
            cache = Path(directory)
            p_report = cache / "p-convergence.json"
            h_report = cache / "h-convergence.json"
            calls = []

            def record(command, check=True):
                calls.append([str(value) for value in command])
                return 0

            with (
                mock.patch.object(PREPARE, "run", side_effect=record),
                mock.patch.object(
                    PREPARE,
                    "run_probe_convergence",
                    return_value=(0, p_report, {"Passed": True}),
                ),
                mock.patch.object(
                    PREPARE,
                    "prepare_probe_mesh_calibration",
                    return_value=cache / "coarse-calibration",
                ),
                mock.patch.object(
                    PREPARE,
                    "run_mesh_convergence",
                    return_value=(1, h_report, {"Passed": False}),
                ),
            ):
                with self.assertRaisesRegex(
                    RuntimeError, "Corner mesh convergence failed"
                ):
                    PREPARE.build_corner(
                        coupon, args, process_parameters(), cache
                    )

            self.assertFalse(any(command[0] == "palace" for command in calls))
            qualification = next(cache.glob("corner-*/qualification.json"))
            result = PREPARE.load_json(qualification)
            self.assertEqual(result["MeshConvergenceReport"], str(h_report))
            self.assertFalse(result["Passed"])

    def test_geometry_coverage_padding_adds_missing_zero_interface(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for kind in ("thin", "fabricated"):
                postpro = root / "postpro" / f"spatial_{kind}"
                postpro.mkdir(parents=True)
                (postpro / "domain-response-matrix.csv").write_text(
                    "basis_i,basis_j,Q_ij (J)\n"
                    "1,1,1.0\n1,2,0.0\n2,2,1.0\n"
                )
                (postpro / "surface-response-matrix.csv").write_text(
                    "interface,edge,R (m),basis_i,basis_j,Q_ij (J),"
                    "Q_ij normal (J),Q_ij tangential (J),Q_total_ij (J),"
                    "Q_total_ij normal (J),Q_total_ij tangential (J)\n"
                    "2,1,2e-6,1,1,1,1,0,1,1,0\n"
                    "2,1,2e-6,1,2,0,0,0,0,0,0\n"
                    "2,1,2e-6,2,2,1,1,0,1,1,0\n"
                )
            PREPARE.pad_missing_surface_response_matrices(root, 2, 2.0)
            for kind in ("thin", "fabricated"):
                path = root / "postpro" / f"spatial_{kind}" / "surface-response-matrix.csv"
                with path.open(newline="") as stream:
                    rows = list(csv.DictReader(stream))
                zero = [row for row in rows if int(row["interface"]) == 1]
                self.assertEqual(len(zero), 3)
                self.assertTrue(
                    all(float(row["Q_total_ij (J)"]) == 0.0 for row in zero)
                )

    def test_geometry_coverage_stamp_is_explicitly_unqualified(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "process-library.json"
            path.write_text(
                '{"Version": 3, "Models": [{"Name": "model"}]}\n'
            )
            PREPARE.stamp_geometry_coverage_only(path)
            qualification = PREPARE.load_json(path)["Models"][0][
                "BoundaryLawQualification"
            ]
            self.assertEqual(qualification["Status"], "Unqualified")
            self.assertEqual(
                qualification["Calibration"], "GeometryCoverageOnly"
            )

    def test_execute_preserves_successes_after_independent_coupon_failure(self):
        failed = endpoint_coupon()
        failed["Id"] = "failed"
        failed["Preparation"] = {"Method": "SpatialCoupon"}
        qualified = endpoint_coupon()
        qualified["Id"] = "qualified"
        qualified["Preparation"] = {"Method": "SpatialCoupon"}
        plan = {
            "SourceManifest": "requirements.json",
            "Coupons": [failed, qualified],
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.json"
            source.write_text('{"Version": 3, "Fabrication": {}}\n')
            args = SimpleNamespace(
                cache=root / "cache",
                output=root / "output",
                name="test-library",
            )

            def build(coupon, args, parameters, cache):
                if coupon["Id"] == "failed":
                    raise RuntimeError("probe did not converge")
                return root / "qualified.json", root / "qualification.json"

            with (
                mock.patch.object(PREPARE, "process_parameters", return_value={}),
                mock.patch.object(PREPARE, "build_spatial", side_effect=build),
                mock.patch.object(PREPARE, "run", return_value=0),
            ):
                complete = PREPARE.execute(plan, source, args)

            self.assertFalse(complete)
            self.assertFalse(plan["Execution"]["Complete"])
            self.assertEqual(len(plan["Execution"]["Failures"]), 1)
            self.assertEqual(plan["Execution"]["Failures"][0]["Ids"], ["failed"])
            manifest = PREPARE.load_json(
                root / "output" / "library" / "qualification-manifest.json"
            )
            self.assertEqual(
                manifest["GeneratedLibraries"], [str(root / "qualified.json")]
            )
            self.assertFalse(manifest["Passed"])

    def test_execute_parallel_coupon_jobs_preserves_plan_order(self):
        coupons = []
        for index in range(3):
            coupon = endpoint_coupon()
            coupon["Id"] = f"coupon-{index}"
            coupon["Preparation"] = {"Method": "SpatialCoupon"}
            coupons.append(coupon)
        plan = {"SourceManifest": "requirements.json", "Coupons": coupons}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.json"
            source.write_text('{"Version": 3, "Fabrication": {}}\n')
            args = SimpleNamespace(
                cache=root / "cache",
                output=root / "output",
                name="parallel-library",
                coupon_jobs=2,
                coverage_only=False,
            )

            def build(coupon, args, parameters, cache):
                index = coupon["Id"].split("-")[-1]
                return root / f"library-{index}.json", root / f"report-{index}.json"

            with (
                mock.patch.object(PREPARE, "process_parameters", return_value={}),
                mock.patch.object(PREPARE, "build_spatial", side_effect=build),
                mock.patch.object(PREPARE, "run", return_value=0),
            ):
                self.assertTrue(PREPARE.execute(plan, source, args))

            manifest = PREPARE.load_json(
                root / "output" / "library" / "qualification-manifest.json"
            )
            self.assertEqual(
                manifest["GeneratedLibraries"],
                [str(root / f"library-{index}.json") for index in range(3)],
            )

    def test_execute_coverage_only_reports_geometry_complete_not_qualified(self):
        coupon = endpoint_coupon()
        coupon["Id"] = "coverage"
        coupon["Preparation"] = {"Method": "SpatialCoupon"}
        plan = {
            "SourceManifest": "requirements.json",
            "Coupons": [coupon],
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.json"
            source.write_text('{"Version": 3, "Fabrication": {}}\n')
            args = SimpleNamespace(
                cache=root / "cache",
                output=root / "output",
                name="coverage-library",
                coverage_only=True,
            )
            with (
                mock.patch.object(PREPARE, "process_parameters", return_value={}),
                mock.patch.object(
                    PREPARE,
                    "build_spatial",
                    return_value=(
                        root / "coverage.json",
                        root / "qualification.json",
                    ),
                ),
                mock.patch.object(PREPARE, "run", return_value=0),
            ):
                complete = PREPARE.execute(plan, source, args)

            self.assertTrue(complete)
            self.assertTrue(plan["Execution"]["GeometryComplete"])
            self.assertFalse(plan["Execution"]["Complete"])
            self.assertFalse(plan["Execution"]["QualificationRequired"])
            manifest = PREPARE.load_json(
                root / "output" / "library" / "qualification-manifest.json"
            )
            self.assertTrue(manifest["GeometryComplete"])
            self.assertFalse(manifest["Passed"])

    def test_execute_probe_study_stops_without_generated_library(self):
        coupon = endpoint_coupon()
        coupon["Id"] = "probe"
        coupon["Preparation"] = {"Method": "SpatialCoupon"}
        plan = {"SourceManifest": "requirements.json", "Coupons": [coupon]}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.json"
            source.write_text('{"Version": 3, "Fabrication": {}}\n')
            args = SimpleNamespace(
                cache=root / "cache",
                output=root / "output",
                name="probe-library",
                coupon_jobs=1,
                coverage_only=False,
                probe_study_only=True,
            )
            with (
                mock.patch.object(PREPARE, "process_parameters", return_value={}),
                mock.patch.object(
                    PREPARE,
                    "build_spatial",
                    return_value=(None, root / "probe-report.json"),
                ),
                mock.patch.object(PREPARE, "run", return_value=0),
            ):
                self.assertTrue(PREPARE.execute(plan, source, args))

            self.assertTrue(plan["Execution"]["ProbeStudyComplete"])
            self.assertFalse(plan["Execution"]["GeometryComplete"])
            self.assertFalse(plan["Execution"]["Complete"])
            manifest = PREPARE.load_json(
                root / "output" / "library" / "qualification-manifest.json"
            )
            self.assertEqual(manifest["GeneratedLibraries"], [])
            self.assertTrue(manifest["ProbeStudyComplete"])
            self.assertTrue(manifest["Passed"])

    def test_execute_reports_unqualified_boundary_law_physics(self):
        coupon = endpoint_coupon()
        coupon["Id"] = "impedance-endpoint"
        coupon["BoundaryCondition"] = {
            "Type": "Impedance",
            "Rs": 0.01,
            "Ls": 1.0e-12,
        }
        coupon["Geometry"]["Arms"][0]["BoundaryCondition"] = coupon[
            "BoundaryCondition"
        ]
        coupon["Preparation"] = {
            "Method": "SpatialCoupon",
            "BoundaryLawQualification": "Missing",
        }
        plan = {
            "SourceManifest": "requirements.json",
            "Coupons": [coupon],
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.json"
            source.write_text('{"Version": 3, "Fabrication": {}}\n')
            args = SimpleNamespace(
                cache=root / "cache",
                output=root / "output",
                name="test-library",
            )

            with (
                mock.patch.object(PREPARE, "process_parameters", return_value={}),
                mock.patch.object(
                    PREPARE,
                    "build_spatial",
                    return_value=(
                        root / "qualified.json",
                        root / "qualification.json",
                    ),
                ),
                mock.patch.object(PREPARE, "run", return_value=0),
            ):
                complete = PREPARE.execute(plan, source, args)

            self.assertTrue(complete)
            manifest = PREPARE.load_json(
                root / "output" / "library" / "qualification-manifest.json"
            )
            self.assertTrue(manifest["Passed"])
            self.assertEqual(
                manifest["BoundaryLawPhysics"],
                {
                    "Complete": False,
                    "UnqualifiedCoupons": ["impedance-endpoint"],
                },
            )

    def test_missing_probe_report_is_a_failed_result(self):
        args = SimpleNamespace(
            palace=Path("palace"),
            orders=[2, 3],
            ranks=1,
            max_fabricated_matrix_change=5.0,
            max_fabricated_energy_change=10.0,
            max_domain_defect_change=5.0,
            force=False,
        )
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "convergence"
            with mock.patch.object(PREPARE, "run", return_value=1):
                code, report, result = PREPARE.run_probe_convergence(
                    Path(directory), output, args
                )
            self.assertEqual(code, 1)
            self.assertEqual(report, output / "probe-convergence.json")
            self.assertFalse(result["Passed"])
            self.assertIn("did not write", result["Failure"])

    def test_missing_mesh_probe_report_is_a_failed_result(self):
        args = SimpleNamespace(
            palace=Path("palace"),
            orders=[2, 3],
            ranks=1,
            max_fabricated_matrix_change=5.0,
            max_fabricated_energy_change=10.0,
            max_domain_defect_change=5.0,
            force=False,
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "mesh-convergence"
            calibrations = [
                ("h-2", root / "coarse"),
                ("h-1", root / "fine"),
            ]
            with mock.patch.object(PREPARE, "run", return_value=1) as runner:
                code, report, result = PREPARE.run_mesh_convergence(
                    calibrations, output, args
                )
            command = runner.call_args.args[0]
            self.assertIn("--fixed-order", command)
            self.assertEqual(command[command.index("--fixed-order") + 1], 3)
            self.assertEqual(code, 1)
            self.assertEqual(report, output / "probe-convergence.json")
            self.assertFalse(result["Passed"])
            self.assertEqual(result["Study"], "MeshResolution")

    def test_domain_defect_reports_full_response_scale(self):
        def response(thin, fabricated):
            return {
                "responses": {
                    "thin": {
                        "domain": np.array([[thin]]),
                        "surfaces": {},
                    },
                    "fabricated": {
                        "domain": np.array([[fabricated]]),
                        "surfaces": {},
                    },
                },
                "interface_names": {},
            }

        args = SimpleNamespace(
            max_fabricated_matrix_change=5.0,
            max_fabricated_energy_change=10.0,
            max_domain_defect_change=6.0,
        )
        passed, results = CORNER_CONVERGENCE.compare_cases(
            response(8.85, 9.8),
            response(9.0, 10.0),
            args,
            True,
        )
        defect = next(row for row in results if row["Kind"] == "defect")

        self.assertTrue(passed)
        self.assertAlmostEqual(defect["MatrixChangePercent"], 5.0)
        self.assertAlmostEqual(
            defect["NormRelativeToFabricatedPercent"], 10.0
        )
        self.assertAlmostEqual(
            defect["ChangeRelativeToFabricatedPercent"], 0.5
        )

    def test_process_resolution_rejects_underresolved_fabrication(self):
        parameters = {
            "metal_thickness": 0.1,
            "overetch": 0.05,
            "top_radius": 0.01,
            "bottom_radius": 0.02,
        }
        report = PREPARE.process_resolution(parameters, 0.02, 2.0)
        self.assertTrue(report["Passed"])
        self.assertEqual(report["LimitingFeature"], "OveretchDepth")
        self.assertEqual(report["RoundingMinimumCirclePoints"], 24)
        with self.assertRaisesRegex(ValueError, "only 1 elements"):
            PREPARE.process_resolution(parameters, 0.05, 2.0)


if __name__ == "__main__":
    unittest.main()
