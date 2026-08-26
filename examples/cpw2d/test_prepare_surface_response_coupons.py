#!/usr/bin/env python3

import csv
import importlib.util
import json
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
            self.assertEqual(combined["Models"], [])
            self.assertEqual(combined["Fabrication"], {"InterfaceLayers": {}})

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
        for angle in (0.0, 180.0, float("nan")):
            requirement = {
                "Topology": "ConvexCorner",
                "Geometry": {"AngleDegrees": angle, "CornerRadius": 0.0},
                "BoundaryCondition": pec(),
            }
            self.assertEqual(
                PREPARE.preparation(requirement)["Method"], "Unsupported"
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

    def test_corner_mesh_refinement_scales_only_fine_size(self):
        args = SimpleNamespace(
            force=True,
            julia="julia",
            julia_project=None,
            matching_radius=2.0,
            corner_lc_fine=0.02,
            corner_lc_far=0.3,
        )
        spec = {"Mesh": {"Order": 2}}
        commands = []
        with (
            tempfile.TemporaryDirectory() as directory,
            mock.patch.object(
                PREPARE,
                "run",
                side_effect=lambda command: commands.append(
                    [str(value) for value in command]
                ),
            ),
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
        self.assertEqual(len(commands), 2)
        for command in commands:
            self.assertEqual(command[command.index("--angle") + 1], "45.0")
            self.assertEqual(command[command.index("--lc-fine") + 1], "0.04")
            self.assertEqual(command[command.index("--lc-far") + 1], "0.3")
            self.assertEqual(command[command.index("--mesh-order") + 1], "2")

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
        with self.assertRaisesRegex(ValueError, "exceeds 12R"):
            SPATIAL.matching_support_points(
                np.asarray([-6.1, -3.0, -1.0]),
                np.asarray([6.1, 3.0, 1.0]),
                frame,
                1.0,
            )

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
            model = PREPARE.load_json(path)["Models"][0]
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
