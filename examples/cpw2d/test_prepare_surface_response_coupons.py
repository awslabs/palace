#!/usr/bin/env python3

import copy
import csv
import importlib.util
import json
import os
import shutil
import subprocess
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

    def test_combiner_refuses_two_models_of_one_signature_before_writing(self):
        """OI-LIBRARY-DUPLICATE-SIGNATURE (decisions 1103 / 1126 (5) / 1128 Q4): two models of one canonical Signature under
        DIFFERENT Names (the one-reference and the -pieces3 form of one key; the Edges differ, the Signature does not) are
        refused by name before any model directory is written; the same Signature under one Name is the name refusal."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            signature = {"Type": "ParallelEdgeCluster", "Edges": [{"Conductor": 1, "GapSide": 1, "OffsetOverR": 0.0},
                                                                  {"Conductor": 1, "GapSide": -1, "OffsetOverR": 1.05}]}
            first = self.write_library(root / "a" / "process-library.json", "one-reference", ["parallel-edge-cluster-9c177156962b"])
            second = self.write_library(root / "b" / "process-library.json", "pieces", ["parallel-edge-cluster-9c177156962b-pieces3"])
            for path, pieces in ((first, False), (second, True)):
                library = PREPARE.load_json(path)
                # the same Signature with its keys in another order (the canonical form sorts keys; the Edges ORDER is the key's)
                library["Models"][0]["Signature"] = {k: signature[k] for k in (("Edges", "Type") if pieces else ("Type", "Edges"))}
                if pieces:
                    library["Models"][0]["LocalPiecesRule"] = "LocalPiecesOfParallelEdges"
                path.write_text(json.dumps(library, indent=2) + "\n")
            output = root / "combined"
            with mock.patch.object(sys, "argv", ["combine_process_libraries.py", "--output", str(output), str(first), str(second)]):
                with self.assertRaisesRegex(ValueError, "OI-LIBRARY-DUPLICATE-SIGNATURE.*parallel-edge-cluster-9c177156962b, "
                                                        "parallel-edge-cluster-9c177156962b-pieces3") as refused:
                    COMBINER.main()
            self.assertEqual(type(refused.exception).__name__, "DuplicateSignatureError")
            self.assertFalse((output / "process-library.json").exists())
            self.assertEqual(sorted(p.name for p in output.iterdir()), [])     # refused before the first model directory
            # the two forms in SEPARATE libraries stay loadable (the 1103 rule of record)
            self.combine(root / "only-first", first)
            self.combine(root / "only-second", second)
            # a Signature-less model never collides on the signature (its Name is the only key)
            third = self.write_library(root / "c" / "process-library.json", "corners", ["convex-corner-45deg", "convex-corner-90deg"])
            self.assertEqual([m["Name"] for m in self.combine(root / "with-corners", first, third)["Models"]],
                             ["parallel-edge-cluster-9c177156962b", "convex-corner-45deg", "convex-corner-90deg"])

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

    def test_replace_by_signature_puts_the_piece_form_in_the_one_reference_models_slot(self):
        """Decision 1103's rule in the assembly (idfix R0-b, decisions 1136 (1) / 1140 (4)): the piece form of a key (another
        Name, the -pieces<P> suffix, the SAME canonical Signature) REPLACES the base's one-reference model of that Signature
        in place (--replace-by-signature NEWNAME=<record> + --supersede-decision), recorded under Supersedes with Kind
        Signature; without the flag the two models of one signature are refused by name (OI-LIBRARY-DUPLICATE-SIGNATURE)
        and nothing is written; a library never holds two models of one signature."""
        pec = {"Type": "ParallelEdgeCluster", "Edges": [{"OffsetOverR": o, "GapSide": d, "Conductor": 1}
                                                       for o, d in zip([0.0, 1.5789474, 3.1578947, 4.7368421], [1, -1, 1, -1])]}
        other = {"Type": "ParallelEdgeCluster", "Edges": [{"OffsetOverR": 0.0, "GapSide": 1, "Conductor": 1},
                                                         {"OffsetOverR": 1.0526316, "GapSide": -1, "Conductor": 1}]}

        def stamp(path, signatures):
            library = PREPARE.load_json(path)
            for model in library["Models"]:
                if model["Name"] in signatures:
                    model["Signature"] = signatures[model["Name"]]
            PREPARE.write_json(path, library)

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            base = self.write_library(root / "base" / "process-library.json", "base",
                                      ["alpha", "parallel-edge-cluster-ff3a5bee4e2d", "omega"])
            stamp(base, {"parallel-edge-cluster-ff3a5bee4e2d": pec, "omega": other})
            pieces = self.write_library(root / "pieces" / "process-library.json", "pieces", ["parallel-edge-cluster-0123456789ab-pieces3"])
            stamp(pieces, {"parallel-edge-cluster-0123456789ab-pieces3": pec})
            record = root / "pieces3-qualification.json"
            record.write_text(json.dumps({"Passed": True}) + "\n")
            decision = "1136 (1): the 2D one-label multi-piece models of v3-j in the piece form -> v3-k-2d"
            plain = self.combine(root / "plain", base)
            with self.assertRaisesRegex(ValueError, "OI-LIBRARY-DUPLICATE-SIGNATURE.*parallel-edge-cluster-0123456789ab-pieces3, parallel-edge-cluster-ff3a5bee4e2d"):
                self.combine(root / "refused", base, pieces)
            self.assertFalse((root / "refused" / "process-library.json").exists())
            replaced = self.combine(root / "replaced", "--replace-by-signature", f"parallel-edge-cluster-0123456789ab-pieces3={record}",
                                    "--supersede-decision", decision, base, pieces)
            self.assertEqual([model["Name"] for model in replaced["Models"]], ["alpha", "parallel-edge-cluster-0123456789ab-pieces3", "omega"])
            for before, after in zip(plain["Models"], replaced["Models"]):
                if before["Name"] != "parallel-edge-cluster-ff3a5bee4e2d":
                    self.assertEqual(before, after)                     # every other model entry untouched
                else:
                    self.assertEqual(after["Signature"], pec)
                    self.assertEqual(after["CombinedFrom"], {"Path": str(pieces.resolve()), "SHA256": COMBINER.sha256(pieces)})
                    self.assertEqual(after["FabricatedMatrix"],
                                     "models/002-parallel-edge-cluster-0123456789ab-pieces3/fabricated-domain-response-matrix.csv")
                    self.assertEqual((root / "replaced" / after["FabricatedMatrix"]).read_text(),
                                     "pieces:parallel-edge-cluster-0123456789ab-pieces3:fabricated-domain-response-matrix.csv\n")
            self.assertEqual(len(replaced["Supersedes"]), 1)
            entry = replaced["Supersedes"][0]
            self.assertEqual((entry["Kind"], entry["Name"], entry["NewName"]),
                             ("Signature", "parallel-edge-cluster-ff3a5bee4e2d", "parallel-edge-cluster-0123456789ab-pieces3"))
            self.assertEqual(entry["SignatureHash"], __import__("signature_library").signature_hash(pec))
            self.assertEqual(entry["Decision"], decision)
            self.assertEqual(entry["Record"], {"Path": str(record.resolve()), "SHA256": COMBINER.sha256(record)})
            self.assertEqual(entry["Rule"], COMBINER.REPLACE_BY_SIGNATURE_RULE)
            self.assertNotEqual(entry["BaseModelSHA"], entry["NewModelSHA"])
            # one model per signature after the replacement, by construction
            self.assertEqual(len({__import__("signature_library").signature_hash(m["Signature"]) for m in replaced["Models"] if "Signature" in m}), 2)
            # Refusals (nothing written): a new name no later input supplies, a new model without a Signature, a Signature the
            # base does not hold, a name the base already holds (--supersede's case), no decision text.
            unsigned = self.write_library(root / "unsigned" / "process-library.json", "unsigned", ["lonely-pieces2"])
            foreign = self.write_library(root / "foreign" / "process-library.json", "foreign", ["foreign-pieces2"])
            stamp(foreign, {"foreign-pieces2": {"Type": "SameConductorGap", "SeparationOverR": 0.5}})
            same_name = self.write_library(root / "same" / "process-library.json", "same", ["omega"])
            stamp(same_name, {"omega": other})
            ruled = ("--supersede-decision", decision)
            for argv, message in (
                    (("--replace-by-signature", f"absent={record}", *ruled, base, pieces), "0 later inputs supply"),
                    (("--replace-by-signature", f"lonely-pieces2={record}", *ruled, base, unsigned), "carries no Signature"),
                    (("--replace-by-signature", f"foreign-pieces2={record}", *ruled, base, foreign), "holds 0 model"),
                    (("--replace-by-signature", f"omega={record}", *ruled, base, same_name), "already holds a model 'omega'"),
                    (("--replace-by-signature", f"parallel-edge-cluster-0123456789ab-pieces3={record}", base, pieces),
                     "requires --supersede-decision TEXT")):
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
            "MeshSizing": "KnotGapAndRingSpacing",
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
        # The metal-graded rule (the concave-MS remedy's decide step, decisions 590 / 629 O-9)
        # is requested by name on a sharp or a rounded corner (RingSize 20); it is NOT the
        # default (the defaults above are unchanged).
        for identifier, radius in (("sharp-90", 0.0), ("rounded-90", 0.5)):
            graded = generator_command(
                corner(identifier, radius),
                arguments(corner_trace_basis="all-rings-follow-metal-metal-graded"),
            )
            self.assertEqual(graded[:2], ("all-rings-follow-metal-metal-graded", "20"))
            self.assertEqual(
                graded[2],
                {
                    "RingSize": 20,
                    "TraceBasis": "all-rings-follow-metal-metal-graded",
                    "TraceResolvability": resolvability,
                },
            )
        self.assertEqual(PREPARE.CORNER_TRACE_BASIS_DEFAULT, "all-rings-follow-metal")
        # The convexity-conditional SELECTOR (decision 739 / S7-DESIGN 3.2 R2): a CONVEX node
        # resolves to the default rule (its coupon spec is the default build's: byte-identical
        # basis files by construction); a CONCAVE node to the metal-graded rule (RingSize 20).
        selector = PREPARE.CORNER_TRACE_BASIS_CONCAVE_GRADED_SELECTOR
        self.assertEqual(selector, "all-rings-follow-metal-concave-graded")
        for identifier, radius in (("sharp-90", 0.0), ("rounded-90", 0.5)):
            convex_selected = generator_command(
                corner(identifier, radius), arguments(corner_trace_basis=selector)
            )
            self.assertEqual(
                convex_selected, generator_command(corner(identifier, radius), arguments())
            )
            concave_coupon = dict(corner(identifier, radius), Topology="ConcaveCorner")
            concave_selected = generator_command(
                concave_coupon, arguments(corner_trace_basis=selector)
            )
            self.assertEqual(
                concave_selected,
                generator_command(
                    concave_coupon,
                    arguments(corner_trace_basis="all-rings-follow-metal-metal-graded"),
                ),
            )
            self.assertEqual(
                concave_selected[:2], ("all-rings-follow-metal-metal-graded", "20")
            )
        self.assertEqual(
            PREPARE.corner_trace_basis(
                "c", SimpleNamespace(corner_trace_basis=selector), "concave"
            ),
            "all-rings-follow-metal-metal-graded",
        )
        self.assertEqual(
            PREPARE.corner_trace_basis(
                "c", SimpleNamespace(corner_trace_basis=selector), "convex"
            ),
            "all-rings-follow-metal",
        )
        with self.assertRaisesRegex(ValueError, "no rule for topology"):
            PREPARE.corner_trace_basis(
                "c", SimpleNamespace(corner_trace_basis=selector), "straight"
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
                                # the generator writes the coupon's Edges on the model: the match is by
                                # them (Point / GapDirection / Conductor), never by topology alone
                                # (OI-PLANNER-SKIPPED-SIGNATURE-STAMP, decisions 1117 / 1137)
                                "Edges": [
                                    {"Point": edge["Point"], "GapDirection": edge["GapDirection"], "Conductor": edge["Conductor"]}
                                    for edge in overlapping_spatial_coupon()["Geometry"]["Edges"]
                                ],
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

    def test_stamp_library_signatures_refuses_two_models_of_one_signature_before_writing(self):
        """OI-LIBRARY-DUPLICATE-SIGNATURE (decisions 1103 / 1126 (5) / 1128 Q4) at the planner's merge: a base model already
        carrying the coupon's Signature and this run's model of the same key (another Name: the -pieces3 form; both match the
        coupon's geometry) would both carry one Signature after the stamping -> refused by name, the file left as written."""
        signature = {"Type": "SameConductorStrip", "SeparationOverR": 1.0526316,
                     "Edges": [{"OffsetOverR": 0.0, "GapSide": -1, "Conductor": 1}, {"OffsetOverR": 1.0526316, "GapSide": 1, "Conductor": 1}]}
        coupon = {"Id": "strip-2um", "Topology": "SameConductorStrip", "Geometry": {"EdgeCount": 2, "Separation": 2.0000008},
                  "Signature": signature, "Instances": 8}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "process-library.json"
            models = [{"Name": "same-conductor-strip-2um", "Topology": "SameConductorStrip", "Separation": 2.0000008, "Signature": signature},
                      {"Name": "same-conductor-strip-2um-pieces2", "Topology": "SameConductorStrip", "Separation": 2.0000008,
                       "LocalPiecesRule": "LocalPiecesOfParallelEdges"}]
            PREPARE.write_json(path, {"Version": 3, "Models": models})
            before = path.read_text()
            with self.assertRaisesRegex(ValueError, "OI-LIBRARY-DUPLICATE-SIGNATURE.*same-conductor-strip-2um, same-conductor-strip-2um-pieces2"):
                PREPARE.stamp_library_signatures(path, [coupon])
            self.assertEqual(path.read_text(), before)       # fail closed before the write
            PREPARE.write_json(path, {"Version": 3, "Models": models[1:]})
            self.assertEqual(PREPARE.stamp_library_signatures(path, [coupon]), 1)      # one model per signature: stamped

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
                {
                    "MeshSizing": "KnotGapAndRingSpacing",
                    "MinimumActiveNodesPerOrderSquared": 5,
                },
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

    def test_execute_refuses_a_piece_form_spatial_coupon(self):
        """Fix round R-3 (OI-PLANNER-EXECUTE-SPATIAL-ONE-REFERENCE; decisions 1157 Q1 / 1230): a SpatialCoupon planned in the
        piece form is REFUSED by --execute before any build (build_spatial has no --local-pieces); --skip-methods SpatialCoupon
        leaves it to device_coupons; a one-reference plan builds as before."""
        pieces = endpoint_coupon()
        pieces["Id"] = "pieces"
        pieces["Preparation"] = {"Method": "SpatialCoupon", "LocalPieces": True, "LocalPiecesReason": "test"}
        plain = endpoint_coupon()
        plain["Id"] = "plain"
        plain["Preparation"] = {"Method": "SpatialCoupon", "LocalPieces": False}
        plan = {"SourceManifest": "requirements.json", "Coupons": [pieces, plain]}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.json"
            source.write_text('{"Version": 3, "Fabrication": {}}\n')
            args = SimpleNamespace(cache=root / "cache", output=root / "output", name="test-library")
            calls = []

            def build(coupon, args, parameters, cache):
                calls.append(coupon["Id"])
                return root / "built.json", root / "qualification.json"

            with (
                mock.patch.object(PREPARE, "process_parameters", return_value={}),
                mock.patch.object(PREPARE, "build_spatial", side_effect=build),
                mock.patch.object(PREPARE, "run", return_value=0),
            ):
                with self.assertRaisesRegex(ValueError, "refuses the SpatialCoupon.*piece form.*\\['pieces'\\].*device_coupons.py") as context:
                    PREPARE.execute(plan, source, args)
                self.assertIn("fix round R-3", str(context.exception))
                self.assertEqual(calls, [])
                # --skip-methods SpatialCoupon: the family is left to the device builder, nothing refused
                skipping = SimpleNamespace(cache=root / "cache", output=root / "output", name="test-library", skip_methods=["SpatialCoupon"])
                PREPARE.execute(plan, source, skipping)
                self.assertEqual(calls, [])
                self.assertEqual(sorted(s["Ids"][0] for s in plan["Execution"]["Skipped"]), ["pieces", "plain"])
                # the one-reference plan alone builds
                plan["Coupons"] = [plain]
                PREPARE.execute(plan, source, args)
                self.assertEqual(calls, ["plain"])

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

    # --- the held-out escalation rule (decisions 1116 / 1126 (6) / 1128 Q3) ---

    def test_heldout_only_failure_and_escalated_args(self):
        passed = {"Passed": True, "MatrixChecks": [{"Passed": True}], "ConvergenceChecks": [], "HeldoutChecks": [{"Passed": True}]}
        heldout = {"Passed": False, "MatrixChecks": [{"Passed": True}], "ConvergenceChecks": [{"Passed": True}],
                   "HeldoutChecks": [{"Passed": True}, {"Quantity": "interface-2", "RelativeErrorPercent": 11.09, "Passed": False}]}
        matrix = {"Passed": False, "MatrixChecks": [{"Passed": False}], "ConvergenceChecks": [], "HeldoutChecks": [{"Passed": False}]}
        self.assertFalse(PREPARE.heldout_only_failure(passed))
        self.assertTrue(PREPARE.heldout_only_failure(heldout))            # the RG-7 L2 case: MS +11.09 % alone
        self.assertFalse(PREPARE.heldout_only_failure(matrix))            # not a held-out-only failure: no escalation
        args = SimpleNamespace(basis_size=96, samples=1200, orders=[2, 3], force=False)
        escalated = PREPARE.escalated_args(args)
        self.assertEqual((escalated.basis_size, escalated.samples), (144, 1800))   # RG-7's escalation of record
        self.assertEqual((args.basis_size, args.samples), (96, 1200))             # the original untouched
        odd = PREPARE.escalated_args(SimpleNamespace(basis_size=50, samples=50))
        self.assertEqual((odd.basis_size, odd.samples), (76, 76))                 # 75 -> even 76; samples >= basis
        self.assertEqual(PREPARE.HELDOUT_ESCALATION_FACTOR, 1.5)

    def escalation_builder(self, root, outcomes):
        """A builder recording the (basis, samples) it was called with; outcomes[i] = 'heldout' / 'fail' / 'ok' per call."""
        calls = []

        def build(args):
            index = len(calls)
            calls.append((args.basis_size, args.samples))
            path = root / f"attempt-{index}" / "qualification.json"
            path.parent.mkdir(parents=True)
            if outcomes[index] == "ok":
                PREPARE.write_json(path, {"Passed": True, "HeldOutFailure": False})
                return root / f"attempt-{index}" / "process-library.json", path
            PREPARE.write_json(path, {"Passed": False, "HeldOutFailure": outcomes[index] == "heldout"})
            if outcomes[index] == "heldout":
                raise PREPARE.HeldOutQualificationFailure("held-out failed", path)
            raise RuntimeError("probe convergence failed")
        return build, calls

    def test_heldout_escalation_retries_once_at_1p5x_and_records_both_attempts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args = SimpleNamespace(basis_size=96, samples=1200, heldout_escalation=True)
            build, calls = self.escalation_builder(root, ["heldout", "ok"])
            library, qualification = PREPARE.with_heldout_escalation(build, args, "coupon-a")
            self.assertEqual(calls, [(96, 1200), (144, 1800)])
            self.assertEqual(library, root / "attempt-1" / "process-library.json")
            second = PREPARE.load_json(qualification)["HeldOutEscalation"]
            self.assertEqual((second["Attempt"], second["Outcome"], second["Settings"]), (2, "Qualified", {"BasisSize": 144, "Samples": 1800}))
            self.assertEqual(second["FirstAttempt"], {"Qualification": str(root / "attempt-0" / "qualification.json"), "BasisSize": 96, "Samples": 1200})
            first = PREPARE.load_json(root / "attempt-0" / "qualification.json")["HeldOutEscalation"]
            self.assertEqual((first["Attempt"], first["Outcome"], first["EscalatedTo"]), (1, "HeldOutFailure", str(qualification)))
            self.assertEqual(first["Escalation"], {"BasisSize": 144, "Samples": 1800})
            self.assertIn("decision 1116", first["Rule"])

    def test_heldout_escalation_second_failure_is_unqualified_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args = SimpleNamespace(basis_size=96, samples=1200, heldout_escalation=True)
            build, calls = self.escalation_builder(root, ["heldout", "heldout"])
            with self.assertRaisesRegex(RuntimeError, r"coupon-b: Unqualified \(fail closed.*basis 96 / samples 1200.*1.5x escalation basis 144 / samples 1800"):
                PREPARE.with_heldout_escalation(build, args, "coupon-b")
            self.assertEqual(calls, [(96, 1200), (144, 1800)])                      # exactly ONE retry
            second = PREPARE.load_json(root / "attempt-1" / "qualification.json")["HeldOutEscalation"]
            self.assertTrue(second["Unqualified"])
            self.assertEqual((second["Attempt"], second["Outcome"]), (2, "HeldOutFailure"))
            first = PREPARE.load_json(root / "attempt-0" / "qualification.json")["HeldOutEscalation"]
            self.assertTrue(first["Unqualified"])
            self.assertEqual(first["SecondAttempt"], str(root / "attempt-1" / "qualification.json"))
            # a non-held-out failure is never escalated; the rule switched off never retries
            build, calls = self.escalation_builder(root / "matrix", ["fail", "ok"])
            with self.assertRaisesRegex(RuntimeError, "probe convergence failed"):
                PREPARE.with_heldout_escalation(build, args, "coupon-c")
            self.assertEqual(calls, [(96, 1200)])
            build, calls = self.escalation_builder(root / "off", ["heldout", "ok"])
            with self.assertRaises(PREPARE.HeldOutQualificationFailure):
                PREPARE.with_heldout_escalation(build, SimpleNamespace(basis_size=96, samples=1200, heldout_escalation=False), "coupon-d")
            self.assertEqual(calls, [(96, 1200)])

    def test_non_heldout_failure_at_the_escalated_attempt_propagates_fail_closed(self):
        # review d1126 MINOR-5 / decision 1133 Q3 (a): a NON-held-out failure at the 1.5x attempt (e.g. mesh convergence) is
        # not escalated again and propagates as the builder's own error; the first record keeps its Attempt-1 block as
        # written (Outcome HeldOutFailure, Escalation 144 / 1800) with no second-attempt cross-reference (no write beyond 1116)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args = SimpleNamespace(basis_size=96, samples=1200, heldout_escalation=True)
            build, calls = self.escalation_builder(root, ["heldout", "fail"])
            with self.assertRaisesRegex(RuntimeError, "^probe convergence failed$"):
                PREPARE.with_heldout_escalation(build, args, "coupon-e")
            self.assertEqual(calls, [(96, 1200), (144, 1800)])                      # no third attempt
            first = PREPARE.load_json(root / "attempt-0" / "qualification.json")
            self.assertTrue(first["HeldOutFailure"])
            block = first["HeldOutEscalation"]
            self.assertEqual((block["Attempt"], block["Outcome"], block["Escalation"]), (1, "HeldOutFailure", {"BasisSize": 144, "Samples": 1800}))
            self.assertEqual(sorted(block), ["Attempt", "Escalation", "Factor", "Outcome", "Rule", "Settings"])
            second = PREPARE.load_json(root / "attempt-1" / "qualification.json")
            self.assertEqual((second["Passed"], second["HeldOutFailure"]), (False, False))
            self.assertNotIn("HeldOutEscalation", second)

    def test_recorded_heldout_failure_is_not_resolved_and_the_rule_is_default_on(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "qualification.json"
            PREPARE.write_json(path, {"Fingerprint": "k", "Passed": False, "HeldOutFailure": True})
            args = SimpleNamespace(force=False)
            with self.assertRaisesRegex(PREPARE.HeldOutQualificationFailure, "recorded, not re-solved"):
                PREPARE.raise_recorded_heldout_failure(PREPARE.load_json(path), "k", path, args, "Parallel-cluster")
            PREPARE.raise_recorded_heldout_failure(PREPARE.load_json(path), "other-key", path, args, "Parallel-cluster")   # another fingerprint
            PREPARE.raise_recorded_heldout_failure({"Fingerprint": "k", "Passed": False, "HeldOutFailure": False}, "k", path, args, "x")
            PREPARE.raise_recorded_heldout_failure(PREPARE.load_json(path), "k", path, SimpleNamespace(force=True), "x")      # --force rebuilds
            # default ON through the planner's execute: the manifest records the rule and the families
            coupon = endpoint_coupon()
            coupon["Id"] = "pec"
            coupon["Topology"] = "ParallelEdgeCluster"
            coupon["Preparation"] = {"Method": "ParallelClusterCoupon"}
            plan = {"SourceManifest": "requirements.json", "Coupons": [coupon]}
            source = root / "source.json"
            source.write_text('{"Version": 3, "Fabrication": {}}\n')
            args = SimpleNamespace(cache=root / "cache", output=root / "output", name="lib", basis_size=96, samples=1200)
            seen = []

            def build(coupon, args, parameters, cache):
                seen.append((args.basis_size, args.samples))
                if len(seen) == 1:
                    first = root / "first-qualification.json"
                    PREPARE.write_json(first, {"Passed": False, "HeldOutFailure": True})
                    raise PREPARE.HeldOutQualificationFailure("held-out", first)
                second = root / "second-qualification.json"
                PREPARE.write_json(second, {"Passed": True})
                return root / "escalated-library.json", second

            with (mock.patch.object(PREPARE, "process_parameters", return_value={}),
                  mock.patch.object(PREPARE, "build_parallel_cluster", side_effect=build),
                  mock.patch.object(PREPARE, "run", return_value=0)):
                self.assertTrue(PREPARE.execute(plan, source, args))
            self.assertEqual(seen, [(96, 1200), (144, 1800)])
            self.assertEqual(plan["Execution"]["Failures"], [])
            manifest = PREPARE.load_json(root / "output" / "library" / "qualification-manifest.json")
            self.assertEqual(manifest["GeneratedLibraries"], [str(root / "escalated-library.json")])
            self.assertEqual(manifest["HeldOutEscalation"]["Enabled"], True)
            self.assertEqual(manifest["HeldOutEscalation"]["Families"], ["StraightEdgeBuilder", "ParallelClusterCoupon"])
            self.assertEqual(manifest["HeldOutEscalation"]["NotApplicable"], ["CornerCoupon", "SpatialCoupon"])
            self.assertEqual(PREPARE.load_json(root / "second-qualification.json")["HeldOutEscalation"]["Outcome"], "Qualified")

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


V3J_MIRROR = (ASSESSMENT / "stage2-20261004" / "v3h" / "stretchfix" / "local" / "library" / "s2-r1p9-v3-j-local"
              / "process-library-local.json")
RG7_MANIFEST = (ASSESSMENT / "readiness-20261029" / "idfix" / "regrade" / "rg7" / "manifests"
                / "item0-stub-driven-requirements-cluster.json")


def pec_cluster_requirement(offsets, directions, conductors, digest):
    return {"Topology": "ParallelEdgeCluster", "Status": "Missing", "Count": 4, "TotalEdgeLength": 40.0,
            "Geometry": {"EdgeCount": len(offsets),
                         "Edges": [{"Offset": [o, 0.0], "GapDirection": [float(d), 0.0], "Conductor": c}
                                   for o, d, c in zip(offsets, directions, conductors)]},
            "Interfaces": [], "BoundaryCondition": {"Type": "PEC"}, "Hash": digest * 32,
            "Signature": {"Type": "ParallelEdgeCluster",
                          "Edges": [{"OffsetOverR": o / 1.9, "GapSide": d, "Conductor": c}
                                    for o, d, c in zip(offsets, directions, conductors)]},
            "Instances": 1, "DistinctSignatures": 1, "ParameterSpread": 0.0, "ExactParameters": True}


class PlannerLocalPiecesDefault(unittest.TestCase):
    """Decision 1136 (4) (the 1134 result): the planner plans the PIECE form for every coupon whose cross-section / plan
    view holds a label with >= 2 local metal pieces (Preparation.LocalPieces), the one-reference form otherwise;
    --one-reference opts every coupon out; the builders read the planned form."""

    def manifest(self):
        return {"Version": 2, "Library": {"MatchingRadius": 1.9, "Path": "seed.json"},
                "Requirements": [
                    pec_cluster_requirement([0.0, 2.0, 4.0, 6.0], [1, -1, 1, -1], [1, 1, 1, 1], "aa"),   # a one-label CPW: 3 pieces
                    pec_cluster_requirement([0.0, 2.0, 4.0], [1, -1, 1], [1, 2, 2], "bb"),                # two labels, one piece each
                    pec_cluster_requirement([0.0, 2.0, 4.0, 6.0], [1, -1, 1, -1], [1, 2, 2, 1], "cc"),   # grounds 1 (2 pieces), strip 2 (1)
                    {"Topology": "SameConductorGap", "Status": "Missing", "Count": 2, "TotalEdgeLength": 3.0,
                     "Geometry": {"EdgeCount": 2, "Separation": 2.0}, "Interfaces": [], "BoundaryCondition": {"Type": "PEC"},
                     "Signature": {"Type": "SameConductorGap", "SeparationOverR": 1.0526316}, "Hash": "dd" * 32},
                    {"Topology": "IsolatedEdge", "Status": "Missing", "Count": 2, "TotalEdgeLength": 3.0, "Geometry": None,
                     "Interfaces": [], "BoundaryCondition": {"Type": "PEC"}, "Signature": {"Type": "IsolatedEdge"}, "Hash": "ee" * 32},
                    {"Topology": "ConvexCorner", "Status": "Missing", "Count": 1, "TotalEdgeLength": 1.0,
                     "Geometry": {"AngleDegrees": 90.0, "CornerRadius": 0.0}, "Interfaces": [], "BoundaryCondition": {"Type": "PEC"},
                     "Signature": {"Type": "ConvexCorner", "AngleDegrees": 90.0, "CornerRadiusOverR": 0.0}, "Hash": "ff" * 32},
                ]}

    def test_the_piece_form_is_planned_for_multi_piece_labels_only(self):
        plan = PREPARE.plan_from_manifest(Path("m.json"), self.manifest(), Path("seed.json"), {}, False)
        forms = {coupon["Hash"][:2]: coupon["Preparation"] for coupon in plan["Coupons"] if coupon.get("Hash")}
        self.assertEqual({key: value["LocalPieces"] for key, value in forms.items()},
                         {"aa": True, "bb": False, "cc": True, "dd": True, "ee": False, "ff": False})
        self.assertEqual(forms["aa"]["PiecesPerLabel"], {"1": 3})
        self.assertEqual(forms["bb"]["PiecesPerLabel"], {"1": 1, "2": 1})
        self.assertEqual(forms["cc"]["PiecesPerLabel"], {"1": 2, "2": 1})
        self.assertEqual(forms["dd"]["PiecesPerLabel"], {"1": 2})
        self.assertNotIn("PiecesPerLabel", forms["ee"])          # a class without a piece form
        self.assertIn("without a piece form", forms["ff"]["LocalPiecesReason"])
        self.assertEqual(plan["Summary"]["LocalPiecesForm"], 3)
        self.assertEqual(len(plan["LocalPiecesForm"]["Coupons"]), 3)
        self.assertFalse(plan["LocalPiecesForm"]["OneReference"])
        # the builders read the planned form: the CPW's spec / Name are the piece form's, the two-label stack's are not
        by_hash = {coupon["Hash"][:2]: coupon for coupon in plan["Coupons"] if coupon.get("Hash")}
        args = SimpleNamespace(matching_radius=1.9, orders=[2], mesh_order=2, basis_size=6, samples=12, coupon_depth=1.0,
                                     edge_offset_tolerance=1.0e-3, straight_lc_fine=0.01, straight_lc_far=0.5,
                                     cluster_lc_fine=0.01, cluster_lc_far=0.5, min_process_feature_elements=4,
                                     cluster_h_factors=[2.0, 1.0])
        parameters = {"metal_thickness": 0.1, "overetch": 0.05}
        for key, expected in (("aa", True), ("bb", False), ("cc", True)):
            edges = PREPARE.normalize_parallel_edges(by_hash[key], 1.9)
            spec = PREPARE.parallel_cluster_spec(by_hash[key], edges, args, parameters)
            self.assertEqual("LocalPieces" in spec, expected, key)
            name = PREPARE.parallel_cluster_model_name(PREPARE.fingerprint(spec), edges, PREPARE.coupon_local_pieces(by_hash[key]))
            self.assertEqual(name.endswith("-pieces3") or name.endswith("-pieces2"), expected, key)
        straight = [by_hash["dd"], by_hash["ee"]]
        self.assertTrue(PREPARE.straight_family_local_pieces(straight))
        self.assertIn("LocalPieces", PREPARE.straight_spec(straight, args, parameters))

    def test_one_reference_opts_every_coupon_out_and_is_the_pre_1136_plan(self):
        plan = PREPARE.plan_from_manifest(Path("m.json"), self.manifest(), Path("seed.json"), {}, False, one_reference=True)
        self.assertTrue(all(coupon["Preparation"]["LocalPieces"] is False for coupon in plan["Coupons"]))
        self.assertTrue(all("--one-reference" in coupon["Preparation"]["LocalPiecesReason"] for coupon in plan["Coupons"]))
        self.assertEqual(plan["Summary"]["LocalPiecesForm"], 0)
        self.assertTrue(plan["LocalPiecesForm"]["OneReference"])
        # a plan written before the rule (no Preparation.LocalPieces) is read as the one-reference form
        self.assertFalse(PREPARE.coupon_local_pieces({"Preparation": {"Method": "ParallelClusterCoupon"}}))
        for coupon in plan["Coupons"]:
            del coupon["Preparation"]["LocalPieces"]
        self.assertFalse(PREPARE.straight_family_local_pieces(plan["Coupons"]))
        self.assertFalse(any(PREPARE.coupon_local_pieces(coupon) for coupon in plan["Coupons"]))

    @unittest.skipUnless(V3J_MIRROR.is_file(), "the v3-j local mirror is needed")
    def test_spatial_coupons_are_planned_by_their_plan_view_pieces(self):
        """The 7 one-label multi-piece and 11 two-label multi-piece-label spatial models of v3-j (the R0-b candidate census,
        probes/v3j_local_pieces_v1.out) are planned in the piece form from their Signature; a one-piece-per-label spatial
        model is not; a library without Fabrication metadata leaves the spatial pieces not derivable = Unsupported (1145)."""
        library = json.loads(V3J_MIRROR.read_text())
        spatial = {model["Name"].split("_")[-1]: model for model in library["Models"] if model.get("Topology") == "SpatialEdgeCluster"}
        multi = {"e7f44561bbf0", "d19e77e8fe71", "8e2c1b63713b", "4b0d7bf2ad1b", "17ee51c67d9a", "f6d7ca001864", "66217996479f",
                 "a2f53cd8630e", "e04e3eb5891d", "4388ea4d7f2a", "31b6896dc54d", "7d88636452b7", "5bda1bc70e63", "f813a57b5d43",
                 "e3190c93a8ea", "09e95b2af0a4", "c983b84ada5c", "a287eb4daef0"}
        single = {"13c7cc02fd9b", "da0179c64591", "ada4e3c1920d"}
        requirements = []
        for key in sorted(multi | single):
            model = spatial[key]
            requirements.append({"Topology": "SpatialEdgeCluster", "Status": "Missing", "Count": 1, "TotalEdgeLength": 1.0,
                                 "Geometry": {"EdgeCount": len(model["Edges"]), "Signature": model["Signature"]},
                                 "Interfaces": model["Interfaces"], "BoundaryCondition": model.get("BoundaryCondition", {"Type": "PEC"}),
                                 "Signature": model["Signature"], "Hash": key * 6 + "0" * 4})
        manifest = {"Version": 2, "Library": {"MatchingRadius": library["MatchingRadius"], "Path": str(V3J_MIRROR)},
                    "Requirements": requirements}
        plan = PREPARE.plan_from_manifest(Path("m.json"), manifest, V3J_MIRROR, library, False)
        forms = {coupon["Hash"][:12]: coupon["Preparation"] for coupon in plan["Coupons"]}
        self.assertEqual({key for key, form in forms.items() if form["LocalPieces"]}, multi)
        self.assertEqual(forms["e7f44561bbf0"]["PiecesPerLabel"], {"1": 3})
        self.assertEqual(forms["4b0d7bf2ad1b"]["PiecesPerLabel"], {"1": 5})
        self.assertEqual(forms["4388ea4d7f2a"]["PiecesPerLabel"], {"1": 2, "2": 1})
        self.assertEqual(forms["13c7cc02fd9b"]["PiecesPerLabel"], {"1": 1, "2": 1})
        self.assertEqual(plan["Summary"]["LocalPiecesForm"], 18)
        bare = PREPARE.plan_from_manifest(Path("m.json"), manifest, V3J_MIRROR, {"MatchingRadius": 1.9}, False)
        for coupon in bare["Coupons"]:
            # decision 1145: no Fabrication = the pieces not derivable = Unsupported by name, never the one-reference form
            self.assertFalse(coupon["Preparation"]["LocalPieces"])
            self.assertIn("not derivable", coupon["Preparation"]["LocalPiecesReason"])
            self.assertEqual(coupon["Preparation"]["Method"], "Unsupported")
            self.assertIn("version-3 Fabrication", coupon["Preparation"]["Reason"])
        self.assertEqual(bare["Summary"]["Unsupported"], len(multi | single))
        self.assertEqual(bare["Summary"]["AutomaticallyPreparatable"], 0)

    FABRICATION = {"MetalThickness": 0.1, "OveretchDepth": 0.05, "SidewallAngleDegrees": 90.0, "TopRoundingRadius": 0.0,
                   "BottomRoundingRadius": 0.0, "SubstratePermittivity": 11.45, "LengthUnit": "um",
                   "InterfaceLayers": {"SA": {"Thickness": 0.002, "Permittivity": 4.0}, "MS": {"Thickness": 0.0003, "Permittivity": 11.45},
                                       "MA": {"Thickness": 0.002, "Permittivity": 10.0}}}

    def hourglass_requirement(self):
        """A version-2 SpatialEdgeCluster record whose two metal triangles touch at the origin (the hourglass of
        test_spatial_local_pieces): the one-reference form builds, the piece form is refused by name (decision 1100)."""
        s = 0.5 ** 0.5
        def portion(p0, p1, gap):
            return {"Conductor": 1, "Gap": list(gap), "Interfaces": ["MA", "MS", "SA"], "Law": "{\"Type\":\"PEC\"}",
                    "P": [*p0, *p1]}
        signature = {"Type": "SpatialEdgeCluster", "EdgeCount": 4,
                     "Portions": [portion((-2.0, -2.0), (0.0, 0.0), (s, -s)), portion((0.0, 0.0), (-2.0, 2.0), (s, s)),
                                  portion((2.0, 2.0), (0.0, 0.0), (-s, s)), portion((0.0, 0.0), (2.0, -2.0), (-s, -s))],
                     "Vertices": [{"P": [0.0, 0.0], "TurnDegrees": 90.0, "Type": "ConvexCorner"}]}
        return {"Topology": "SpatialEdgeCluster", "Status": "Missing", "Count": 1, "TotalEdgeLength": 1.0,
                "Geometry": {"EdgeCount": 4, "Signature": signature}, "Signature": signature, "Hash": "ab" * 32,
                "Interfaces": [{"Slot": 0, "Type": "MA", "Coupon": 1}, {"Slot": 0, "Type": "MS", "Coupon": 2},
                               {"Slot": 1, "Type": "MA", "Coupon": 3}, {"Slot": 1, "Type": "MS", "Coupon": 4},
                               {"Slot": 1, "Type": "SA", "Coupon": 5}],
                "BoundaryCondition": {"Type": "PEC"}}

    def test_a_non_derivable_piece_form_is_unsupported_unless_one_reference(self):
        """Decision 1145 Q1 (a) (review d1140 MAJOR-1): a coupon of a class WITH a piece form whose pieces the piece rule
        refuses by name (the spatial vertex-touch, 1100; lane C's non-cross-section parallel edges) is Unsupported by
        default - never the one-reference form 1134 found wrong - and the one-reference form under --one-reference."""
        manifest = {"Version": 2, "Library": {"MatchingRadius": 1.9, "Path": "seed.json"},
                    "Requirements": [self.hourglass_requirement(),
                                     pec_cluster_requirement([0.0, 2.0, 4.0], [1, 1, -1], [1, 1, 1], "cd")]}   # gaps 0 / 1 on one side
        library = {"Version": 3, "MatchingRadius": 1.9, "Fabrication": self.FABRICATION, "Models": []}
        plan = PREPARE.plan_from_manifest(Path("m.json"), manifest, Path("seed.json"), library, False)
        by_hash = {coupon["Hash"][:2]: coupon["Preparation"] for coupon in plan["Coupons"]}
        self.assertEqual(set(by_hash), {"ab", "cd"})
        for key, refusal in (("ab", "metal pieces touch at a vertex"), ("cd", "both have their gap on the +u side")):
            form = by_hash[key]
            self.assertEqual(form["Method"], "Unsupported", key)
            self.assertIn("not derivable", form["Reason"])
            self.assertIn(refusal, form["Reason"])
            self.assertIn("--one-reference", form["Reason"])
            self.assertFalse(form["LocalPieces"])
            self.assertIn(refusal, form["LocalPiecesReason"])
            self.assertFalse({"Driver", "MeshGenerator", "ResponseGenerator"} & set(form))   # not built
        self.assertEqual(plan["Summary"], {"CouponCount": 2, "LocalPiecesForm": 0, "AutomaticallyPreparatable": 0, "Unsupported": 2})
        explicit = PREPARE.plan_from_manifest(Path("m.json"), manifest, Path("seed.json"), library, False, one_reference=True)
        by_hash = {coupon["Hash"][:2]: coupon["Preparation"] for coupon in explicit["Coupons"]}
        self.assertEqual(by_hash["ab"]["Method"], "SpatialCoupon")
        self.assertEqual(by_hash["cd"]["Method"], "ParallelClusterCoupon")
        self.assertTrue(all(form["LocalPieces"] is False and "--one-reference" in form["LocalPiecesReason"] for form in by_hash.values()))
        self.assertEqual(explicit["Summary"], {"CouponCount": 2, "LocalPiecesForm": 0, "AutomaticallyPreparatable": 2, "Unsupported": 0})

    def test_cli_refuses_the_default_flag_together_with_the_opt_out(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "m.json"
            manifest.write_text(json.dumps(self.manifest()))
            (Path(tmp) / "seed.json").write_text(json.dumps({"Version": 3, "MatchingRadius": 1.9, "Models": []}))
            result = subprocess.run([sys.executable, str(CPW2D / "prepare_surface_response_coupons.py"), str(manifest), "--output", str(Path(tmp) / "out"),
                                     "--local-pieces", "--one-reference"], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("mutually exclusive", result.stderr)
            for flags, expected in (([], 3), (["--one-reference"], 0), (["--local-pieces"], 3)):
                subprocess.run([sys.executable, str(CPW2D / "prepare_surface_response_coupons.py"), str(manifest), "--output", str(Path(tmp) / "out"), *flags],
                               check=True, capture_output=True, text=True)
                plan = json.loads((Path(tmp) / "out" / "coupon-plan.json").read_text())
                self.assertEqual(plan["Summary"]["LocalPiecesForm"], expected, flags)


class SkippedCouponNeverStampsASignature(unittest.TestCase):
    """OI-PLANNER-SKIPPED-SIGNATURE-STAMP (decisions 1117 / 1137; owned by the R0-b lane): RG-7's planner run built ONE
    PEC coupon with the stub's SpatialEdgeCluster requirement under --skip-methods SpatialCoupon, and the merged library's
    spatial model spatialedgecluster_edgecount-10_13c7cc02fd9b (v3-j) came out re-stamped with the SKIPPED stub cluster's
    6-edge Signature (the binary refused the library: VerifySpatialEdgesInSignatureFrame). Two generic rules: a spatial
    model matches a coupon by its KEY (model_matches_coupon: the same Signature, or the same Edges), never by topology
    alone; and only the coupons BUILT by the execution stamp the merged library (built_coupons_for_stamping)."""

    def spatial_model(self, name, signature, edges):
        return {"Name": name, "Topology": "SpatialEdgeCluster", "Signature": signature,
                "Edges": [{"Point": [float(k), 0.0, 0.0], "GapDirection": [0.0, 1.0, 0.0], "Conductor": 1} for k in range(edges)],
                "PlanViewBoundary": [{"Conductor": 1, "Segments": []}]}

    def test_a_spatial_model_matches_by_signature_or_edges_never_by_topology(self):
        own = {"Type": "SpatialEdgeCluster", "Box": [0, 0, 1, 1], "Portions": [{"P": [0, 0, 1, 0]}] * 10}
        stub = {"Type": "SpatialEdgeCluster", "Box": [0, 0, 2, 2], "Portions": [{"P": [0, 0, 2, 0]}] * 6}
        model = self.spatial_model("spatialedgecluster_edgecount-10_13c7cc02fd9b", own, 10)
        skipped = {"Id": "spatialedgecluster_edgecount-6_cd97d9e672ba", "Topology": "SpatialEdgeCluster", "Geometry": {"Signature": stub},
                   "Signature": stub}
        self.assertFalse(PREPARE.model_matches_coupon(model, skipped))
        self.assertTrue(PREPARE.model_matches_coupon(model, {"Id": "own", "Topology": "SpatialEdgeCluster", "Geometry": {}, "Signature": own}))
        # a version-1 model (no Signature) matches a coupon carrying the same Edges only
        v1 = {k: v for k, v in model.items() if k != "Signature"}
        geometry = {"Edges": [{"Point": [float(k), 0.0, 0.0], "GapDirection": [0.0, 1.0, 0.0], "Conductor": 1} for k in range(10)]}
        self.assertTrue(PREPARE.model_matches_coupon(v1, {"Id": "v1", "Topology": "SpatialEdgeCluster", "Geometry": geometry}))
        self.assertFalse(PREPARE.model_matches_coupon(v1, skipped))
        other = copy.deepcopy(geometry)
        other["Edges"][3]["Conductor"] = 2
        self.assertFalse(PREPARE.model_matches_coupon(v1, {"Id": "v1b", "Topology": "SpatialEdgeCluster", "Geometry": other}))

    def test_the_rg7_plan_leaves_the_base_spatial_model_as_written(self):
        own = {"Type": "SpatialEdgeCluster", "Box": [0, 0, 1, 1], "Portions": [{"P": [0, 0, 1, 0]}] * 10}
        stub = {"Type": "SpatialEdgeCluster", "Box": [0, 0, 2, 2], "Portions": [{"P": [0, 0, 2, 0]}] * 6}
        pec = {"Type": "ParallelEdgeCluster", "Edges": [{"OffsetOverR": o, "GapSide": d, "Conductor": 1}
                                                       for o, d in zip([0.0, 1.5789474, 3.1578947, 4.7368421], [1, -1, 1, -1])]}
        built_coupon = {"Id": "parallel-edge-cluster-a670", "Topology": "ParallelEdgeCluster", "Signature": pec, "Instances": 1,
                        "Geometry": {"Edges": [{"Offset": [o, 0.0], "GapDirection": [float(d), 0.0], "Conductor": 1}
                                               for o, d in zip([0.0, 3.0, 6.0, 9.0], [1, -1, 1, -1])]},
                        "Preparation": {"Method": "ParallelClusterCoupon"}}
        skipped_coupon = {"Id": "spatialedgecluster_edgecount-6_cd97d9e672ba", "Topology": "SpatialEdgeCluster", "Signature": stub,
                          "Geometry": {"Signature": stub}, "Preparation": {"Method": "SpatialCoupon"}}
        plan_coupons = [built_coupon, skipped_coupon]
        skipped = [{"Ids": [skipped_coupon["Id"]], "Topologies": ["SpatialEdgeCluster"], "Method": "SpatialCoupon",
                    "Reason": "left to another builder (--skip-methods)"}]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "process-library.json"
            models = [self.spatial_model("spatialedgecluster_edgecount-10_13c7cc02fd9b", own, 10),
                      {"Name": "parallel-edge-cluster-b20f56c8e9fa-pieces3", "Topology": "ParallelEdgeCluster",
                       "Edges": [{"Offset": o, "GapDirection": d, "Conductor": 1, "LocalPiece": p}
                                 for o, d, p in zip([0.0, 3.0, 6.0, 9.0], [1, -1, 1, -1], [1, 2, 2, 3])]}]
            PREPARE.write_json(path, {"Version": 3, "Models": models})
            # (1) the execution's stamping set = the built coupons; a skipped Id listed as built is refused by name
            stamping = PREPARE.built_coupons_for_stamping([built_coupon], skipped, [])
            self.assertEqual([coupon["Id"] for coupon in stamping], [built_coupon["Id"]])
            with self.assertRaisesRegex(ValueError, "OI-PLANNER-SKIPPED-SIGNATURE-STAMP"):
                PREPARE.built_coupons_for_stamping(plan_coupons, skipped, [])
            self.assertEqual(PREPARE.stamp_library_signatures(path, stamping), 1)
            stamped = {model["Name"]: model for model in PREPARE.load_json(path)["Models"]}
            self.assertEqual(stamped["spatialedgecluster_edgecount-10_13c7cc02fd9b"]["Signature"], own)   # as written
            self.assertEqual(stamped["parallel-edge-cluster-b20f56c8e9fa-pieces3"]["Signature"], pec)
            # (2) even the whole plan (the pre-fix call) cannot re-stamp the spatial model: no match by topology alone
            PREPARE.write_json(path, {"Version": 3, "Models": models})
            self.assertEqual(PREPARE.stamp_library_signatures(path, plan_coupons), 1)
            self.assertEqual(PREPARE.load_json(path)["Models"][0]["Signature"], own)

    @unittest.skipUnless(V3J_MIRROR.is_file() and RG7_MANIFEST.is_file(), "the v3-j mirror and the RG-7 stub manifest are needed")
    def test_the_rg7_records_reproduce_the_case(self):
        """The RG-7 records: the stub manifest's two requirements (a670 ParallelEdgeCluster + the stub SpatialEdgeCluster
        cd97d9e672ba... Hash c06c2f00) planned with the stub cluster Skipped, against v3-j's models: 13c7cc02fd9b and every other
        v3-j spatial model keep their own Signature (the pre-fix stamping over plan["Coupons"] rewrote 13c7cc02fd9b's)."""
        manifest = json.loads(RG7_MANIFEST.read_text())
        library = json.loads(V3J_MIRROR.read_text())
        plan = PREPARE.plan_from_manifest(RG7_MANIFEST, manifest, V3J_MIRROR, library, False)
        by_topology = {coupon["Topology"]: coupon for coupon in plan["Coupons"]}
        self.assertEqual(set(by_topology), {"ParallelEdgeCluster", "SpatialEdgeCluster"})
        self.assertEqual(by_topology["SpatialEdgeCluster"]["Hash"][:8], "c06c2f00")
        self.assertTrue(by_topology["ParallelEdgeCluster"]["Preparation"]["LocalPieces"])   # the a670 CPW: the piece form (1136 (4))
        spatial = [model for model in library["Models"] if model.get("Topology") == "SpatialEdgeCluster"]
        for model in spatial:
            self.assertFalse(PREPARE.model_matches_coupon(model, by_topology["SpatialEdgeCluster"]), model["Name"])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "process-library.json"
            before = {model["Name"]: model.get("Signature") for model in spatial}
            PREPARE.write_json(path, {"Version": 3, "Models": copy.deepcopy(spatial)})
            skipped = [{"Ids": [by_topology["SpatialEdgeCluster"]["Id"]], "Method": "SpatialCoupon"}]
            stamping = PREPARE.built_coupons_for_stamping([by_topology["ParallelEdgeCluster"]], skipped, [])
            self.assertEqual(PREPARE.stamp_library_signatures(path, stamping), 0)
            self.assertEqual({model["Name"]: model.get("Signature") for model in PREPARE.load_json(path)["Models"]}, before)
            self.assertEqual(PREPARE.stamp_library_signatures(path, plan["Coupons"]), 0)     # the pre-fix call: no match either
            self.assertEqual({model["Name"]: model.get("Signature") for model in PREPARE.load_json(path)["Models"]}, before)


TWIN_OFFSETS_4 = [0.0, 2.0000008, 3.9999997, 6.0000005]          # v3-j ff3a5bee4e2d / 799fccf7b589 (the 2 / 2 / 2 um CPW)
TWIN_OFFSETS_3 = [0.0, 2.0000008, 3.9999997]                     # v3-j e4c250357351 / bf97438bcd9d


def pec_cluster_model(name, offsets, directions, conductors, **extra):
    return {"Name": name, "Topology": "ParallelEdgeCluster",
            "Edges": [{"Offset": o, "GapDirection": d, "Conductor": c} for o, d, c in zip(offsets, directions, conductors)],
            **extra}


def pec_cluster_coupon(identifier, offsets, directions, conductors, digest, law=None):
    coupon = {"Id": identifier, **pec_cluster_requirement(offsets, directions, conductors, digest)}
    if law is not None:
        coupon["BoundaryCondition"] = law
    return coupon


class LabelAwareStampMatcher(unittest.TestCase):
    """OI-PLANNER-STAMP-LABEL-BLIND (decisions 1146 / 1147 (d) / 1148; the R0-b cluster lane's stamp step): model_matches_coupon
    compared a ParallelEdgeCluster's Offsets + GapDirections only, never Edges[].Conductor, at both stamp sites. v3-j holds
    two geometry twins that differ in the labels only: ff3a5bee4e2d (1,1,1,1) / 799fccf7b589 (1,2,2,1) at 0 / 2.0000008 /
    3.9999997 / 6.0000005 um and e4c250357351 (1,1,2) / bf97438bcd9d (1,1,1) at 0 / 2.0000008 / 3.9999997 um, so a plan holding
    both twin coupons made the stamp ambiguous, and a plan holding one twin's coupon over a base holding the other twin
    re-stamped the base model silently. ONE label-aware matcher serves both sites; >= 2 candidates fail by name."""

    CONDUCTIVITY = {"Type": "Conductivity", "Conductivity": 5.8e7, "Permeability": 1.0, "Thickness": 2.0e-7, "External": False}

    def twin_pairs(self):
        d4 = [1, -1, 1, -1]
        d3 = [-1, 1, -1]
        return [
            (pec_cluster_model("parallel-edge-cluster-ff3a5bee4e2d-pieces3", TWIN_OFFSETS_4, d4, [1, 1, 1, 1]),
             pec_cluster_model("parallel-edge-cluster-799fccf7b589", TWIN_OFFSETS_4, d4, [1, 2, 2, 1]),
             pec_cluster_coupon("pec4_f43855f20324", TWIN_OFFSETS_4, d4, [1, 1, 1, 1], "f4"),
             pec_cluster_coupon("pec4_52bf974324b2", TWIN_OFFSETS_4, d4, [1, 2, 2, 1], "52", self.CONDUCTIVITY)),
            (pec_cluster_model("parallel-edge-cluster-bf97438bcd9d-pieces3", TWIN_OFFSETS_3, d3, [1, 1, 1]),
             pec_cluster_model("parallel-edge-cluster-e4c250357351", TWIN_OFFSETS_3, d3, [1, 1, 2]),
             pec_cluster_coupon("pec3_475cf24e121f", TWIN_OFFSETS_3, d3, [1, 1, 1], "47"),
             pec_cluster_coupon("pec3_e4c250357351", TWIN_OFFSETS_3, d3, [1, 1, 2], "e4", self.CONDUCTIVITY)),
        ]

    def test_both_v3j_twin_pairs_stamp_unambiguously_at_both_sites(self):
        for one_label, two_label, one_coupon, two_coupon in self.twin_pairs():
            with self.subTest(pair=one_label["Name"]):
                self.assertTrue(PREPARE.model_matches_coupon(one_label, one_coupon))
                self.assertTrue(PREPARE.model_matches_coupon(two_label, two_coupon))
                self.assertFalse(PREPARE.model_matches_coupon(one_label, two_coupon))      # the label-blind match said True
                self.assertFalse(PREPARE.model_matches_coupon(two_label, one_coupon))
                with tempfile.TemporaryDirectory() as tmp:
                    path = Path(tmp) / "process-library.json"
                    PREPARE.write_json(path, {"Version": 3, "Models": [one_label, two_label]})
                    self.assertEqual(PREPARE.stamp_library_signatures(path, [two_coupon, one_coupon]), 2)
                    PREPARE.stamp_library_boundary_conditions(path, [two_coupon, one_coupon])
                    stamped = {model["Name"]: model for model in PREPARE.load_json(path)["Models"]}
                    self.assertEqual(stamped[one_label["Name"]]["Signature"], one_coupon["Signature"])
                    self.assertEqual(stamped[two_label["Name"]]["Signature"], two_coupon["Signature"])
                    self.assertEqual(stamped[one_label["Name"]]["BoundaryCondition"], {"Type": "PEC"})
                    self.assertEqual(stamped[two_label["Name"]]["BoundaryCondition"], self.CONDUCTIVITY)
                    self.assertFalse(any("-bc-" in name for name in stamped))

    def test_one_twin_coupon_never_restamps_the_other_twin_base_model(self):
        """Effect (ii) of the record: a plan holding bf97's coupon over a base holding e4c2 (its own Signature) stamps the new
        bf97 piece model only; e4c2 keeps its Signature (pre-fix: re-stamped, then refused as a duplicate signature)."""
        one_label, two_label, one_coupon, two_coupon = self.twin_pairs()[1]
        two_label["Signature"] = copy.deepcopy(two_coupon["Signature"])
        two_label["Instances"] = 7
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "process-library.json"
            PREPARE.write_json(path, {"Version": 3, "Models": [two_label, one_label]})
            self.assertEqual(PREPARE.stamp_library_signatures(path, [one_coupon]), 1)
            stamped = {model["Name"]: model for model in PREPARE.load_json(path)["Models"]}
            self.assertEqual(stamped[two_label["Name"]]["Signature"], two_coupon["Signature"])
            self.assertEqual(stamped[two_label["Name"]]["Instances"], 7)
            self.assertEqual(stamped[one_label["Name"]]["Signature"], one_coupon["Signature"])

    def test_zero_based_coupon_labels_read_as_their_canonical_model_labels(self):
        """normalize_parallel_edges shifts a 0-based coupon row to one-based: the matcher reads (0,1,1,0) as (1,2,2,1)."""
        _, two_label, _, two_coupon = self.twin_pairs()[0]
        zero_based = copy.deepcopy(two_coupon)
        for edge, label in zip(zero_based["Geometry"]["Edges"], [0, 1, 1, 0]):
            edge["Conductor"] = label
        self.assertEqual(PREPARE.canonical_conductor_labels([0, 1, 1, 0]), (1, 2, 2, 1))
        self.assertEqual(PREPARE.canonical_conductor_labels([2, 1, 1, 2]), (1, 2, 2, 1))
        self.assertTrue(PREPARE.model_matches_coupon(two_label, zero_based))

    def test_an_ambiguous_library_fails_by_name_at_both_sites(self):
        """Two coupons of ONE geometry AND ONE label pattern (two Ids / two Signatures) still match one model after the label-aware
        fix: refused by name at both sites, each candidate listed with its labels and law; the file is left as written."""
        one_label, _, one_coupon, _ = self.twin_pairs()[0]
        other = copy.deepcopy(one_coupon)
        other["Id"] = "pec4_duplicate_key"
        other["Signature"] = {**one_coupon["Signature"], "Type": "ParallelEdgeCluster", "Edges": one_coupon["Signature"]["Edges"][::-1]}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "process-library.json"
            PREPARE.write_json(path, {"Version": 3, "Models": [one_label]})
            before = path.read_text()
            pattern = (r"model parallel-edge-cluster-ff3a5bee4e2d-pieces3 matches 2 version-2 coupon signatures \(an ambiguous stamp; "
                       r"OI-PLANNER-STAMP-LABEL-BLIND, decision 1148\): pec4_duplicate_key \(Conductor labels \(1, 1, 1, 1\); "
                       r"BoundaryCondition \{\"Type\": \"PEC\"\}\), pec4_f43855f20324 \(Conductor labels \(1, 1, 1, 1\);")
            with self.assertRaisesRegex(ValueError, pattern):
                PREPARE.stamp_library_signatures(path, [one_coupon, other])
            self.assertEqual(path.read_text(), before)
            with self.assertRaisesRegex(ValueError, r"matches 2 coupons \(an ambiguous stamp.*pec4_duplicate_key.*pec4_f43855f20324"):
                PREPARE.stamp_library_boundary_conditions(path, [one_coupon, other])
            self.assertEqual(path.read_text(), before)

    def test_two_laws_on_one_geometry_are_refused_at_the_boundary_condition_site(self):
        """Decision 1148: a response model is computed for ONE boundary law; two laws on one geometry + labels are two coupons /
        two keys, never '-bc-<digest>' variants of one model (the former branch is gone)."""
        pec = {"Id": "scg-2um-pec", "Topology": "SameConductorGap", "Geometry": {"EdgeCount": 2, "Separation": 2.0000008},
               "BoundaryCondition": {"Type": "PEC"}}
        lossy = {**copy.deepcopy(pec), "Id": "scg-2um-conductivity", "BoundaryCondition": self.CONDUCTIVITY}
        model = {"Name": "same-conductor-gap-2um-pieces2", "Topology": "SameConductorGap", "Separation": 2.0000008, "LocalPieces": 2}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "process-library.json"
            PREPARE.write_json(path, {"Version": 3, "Models": [model]})
            before = path.read_text()
            with self.assertRaisesRegex(ValueError, r"same-conductor-gap-2um-pieces2 matches 2 coupons .*scg-2um-conductivity "
                                                    r"\(Conductor labels \(1, 1\); BoundaryCondition \{\"Conductivity\".*scg-2um-pec \("):
                PREPARE.stamp_library_boundary_conditions(path, [pec, lossy])
            self.assertEqual(path.read_text(), before)
            PREPARE.stamp_library_boundary_conditions(path, [lossy])
            stamped = PREPARE.load_json(path)["Models"]
            self.assertEqual([m["Name"] for m in stamped], ["same-conductor-gap-2um-pieces2"])
            self.assertEqual(stamped[0]["BoundaryCondition"], self.CONDUCTIVITY)
            self.assertEqual(stamped[0]["BoundaryLawQualification"]["Status"], "Unqualified")

    def test_paired_topologies_match_by_separation_and_labels(self):
        """SCG by Separation + its labels: a model without Edges carries the paired topology's pattern (SCG / SCS (1,1); DCG
        (1,2)) or its Signature's; a coupon's labels come from its Signature's Edges when its Geometry has none."""
        signature = {"Type": "SameConductorGap", "SeparationOverR": 1.052632,
                     "Edges": [{"OffsetOverR": 0.0, "GapSide": 1, "Conductor": 1}, {"OffsetOverR": 1.052632, "GapSide": -1, "Conductor": 1}]}
        coupon = {"Id": "scg-2um", "Topology": "SameConductorGap", "Geometry": {"EdgeCount": 2, "Separation": 2.0000008},
                  "Signature": signature}
        bare = {"Name": "same-conductor-gap-2um", "Topology": "SameConductorGap", "Separation": 2.0000008}
        self.assertEqual(PREPARE.conductor_labels(bare, "SameConductorGap", name="m"), (1, 1))
        self.assertEqual(PREPARE.conductor_labels(coupon, "SameConductorGap", name="c"), (1, 1))
        self.assertEqual(PREPARE.conductor_labels({"Topology": "DifferentConductorGap"}, "DifferentConductorGap", name="m"), (1, 2))
        self.assertTrue(PREPARE.model_matches_coupon(bare, coupon))
        self.assertTrue(PREPARE.model_matches_coupon({**bare, "Signature": signature}, coupon))
        self.assertFalse(PREPARE.model_matches_coupon({**bare, "Separation": 3.0}, coupon))
        mislabelled = copy.deepcopy(coupon)
        mislabelled["Signature"]["Edges"][1]["Conductor"] = 2
        self.assertFalse(PREPARE.model_matches_coupon(bare, mislabelled))

    def test_an_edge_without_a_conductor_label_is_refused_by_name(self):
        one_label, _, one_coupon, _ = self.twin_pairs()[0]
        del one_label["Edges"][2]["Conductor"]
        with self.assertRaisesRegex(ValueError, r"model parallel-edge-cluster-ff3a5bee4e2d-pieces3: a ParallelEdgeCluster edge carries no Conductor"):
            PREPARE.model_matches_coupon(one_label, one_coupon)
        one_label, _, one_coupon, _ = self.twin_pairs()[0]
        del one_coupon["Geometry"]["Edges"][0]["Conductor"]
        with self.assertRaisesRegex(ValueError, r"coupon pec4_f43855f20324: a ParallelEdgeCluster edge carries no Conductor"):
            PREPARE.model_matches_coupon(one_label, one_coupon)


if __name__ == "__main__":
    unittest.main()
