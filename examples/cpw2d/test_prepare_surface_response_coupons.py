#!/usr/bin/env python3

import importlib.util
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


class PrepareSurfaceResponseCouponsTest(unittest.TestCase):
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
            ("ConvexCorner", {"AngleDegrees": 90.0, "CornerRadius": 0.2}),
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

    def test_plan_reports_finite_impedance_as_unsupported(self):
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
        self.assertEqual(plan["Summary"]["Unsupported"], 2)
        for coupon in plan["Coupons"]:
            self.assertEqual(coupon["Preparation"]["Method"], "Unsupported")
            self.assertIn("Impedance", coupon["Preparation"]["Reason"])

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

    def test_spatial_matching_surface_encloses_process_geometry(self):
        coupon = endpoint_coupon()
        frame, edges = SPATIAL.normalize_geometry(coupon, 2.0)
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
        )
        labels = SPATIAL.conductor_at_points(
            points, edges, 2.0, 0.1, 80.0
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
        frame, edges = SPATIAL.normalize_geometry(coupon, 2.0)
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

    def test_fabricated_edge_localization_uses_metal_substrate_perimeters(self):
        edges = [
            {"Conductor": 1, "InterfaceSlot": 0},
            {"Conductor": 2, "InterfaceSlot": 0},
            {"Conductor": 2, "InterfaceSlot": 1},
        ]
        self.assertEqual(SPATIAL.edge_attributes(edges, False, 0), [3000])
        self.assertEqual(SPATIAL.edge_attributes(edges, True, 0), [5001, 5002])
        self.assertEqual(SPATIAL.edge_attributes(edges, True, 1), [5102])

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
