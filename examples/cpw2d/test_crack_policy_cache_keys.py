#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""P3-39 (upstream integration P3; decisions 1570 / 1586): the crack policy reaches the library-build orchestrators
and their cache keys. Every planner family spec (straight, parallel cluster, corner, spatial) and the straight
builder's coupon specs carry CrackPolicy and fingerprint crack_policy.py, so a coupon solved under one policy is never
reused for a build under the other; the planner passes --crack-policy to the builder / generators it runs and the
straight builder to the pair / isolated-edge generators; the legacy policy stays the default."""

import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

CPW2D = Path(__file__).resolve().parent
sys.path.insert(0, str(CPW2D))
sys.path.insert(0, str(CPW2D.parent / "surface_response_identification"))
sys.path.insert(0, str(CPW2D.parent / "cpw3d_surface" / "spatial_coupon"))
import build_surface_process_library as BUILDER  # noqa: E402
import crack_policy as C  # noqa: E402
import prepare_surface_response_coupons as PLANNER  # noqa: E402

PARAMETERS = {"metal_thickness": 0.1, "overetch": 0.05, "sidewall_angle": 80.0, "top_radius": 0.01,
              "bottom_radius": 0.01, "substrate_permittivity": 11.47, "sa_thickness": 0.002, "sa_permittivity": 4.0,
              "ms_thickness": 0.0003, "ms_permittivity": 11.47, "ma_thickness": 0.03, "ma_permittivity": 10.0}
STRAIGHT = [{"Id": "gap", "Topology": "SameConductorGap", "BoundaryCondition": {"Type": "PEC"},
             "Geometry": {"Separation": 2.0}, "Preparation": {}},
            {"Id": "iso", "Topology": "IsolatedEdge", "BoundaryCondition": {"Type": "PEC"}, "Geometry": {},
             "Preparation": {}}]
CLUSTER = {"Id": "pec", "Topology": "ParallelEdgeCluster", "BoundaryCondition": {"Type": "PEC"}, "Preparation": {},
           "Geometry": {"Edges": [{"Offset": o, "GapDirection": d, "Conductor": 1}
                                  for o, d in zip([0.0, 2.0, 4.0, 6.0], [1, -1, 1, -1])]}}
CORNER = {"Id": "sharp-90", "Topology": "ConvexCorner", "Geometry": {"AngleDegrees": 90.0, "CornerRadius": 0.0},
          "BoundaryCondition": {"Type": "PEC"}}
ENDPOINT = {"Id": "endpoint", "Topology": "Endpoint",
            "Geometry": {"SignatureVersion": 1, "ArmCount": 1, "ArmAnglesDegrees": [],
                         "Arms": [{"Direction": [1.0, 0.0, 0.0], "GapDirection": [0.0, 1.0, 0.0],
                                   "ProcessNormal": [0.0, 0.0, 1.0], "Interval": [0.0, 2.0], "Conductor": 1,
                                   "InterfaceSlot": 0, "BoundaryCondition": {"Type": "PEC"}}]},
            "Interfaces": [{"Slot": 0, "Type": t, "Target": i} for i, t in enumerate(("SA", "MS", "MA"), start=1)],
            "BoundaryCondition": {"Type": "PEC"}, "CoverageStatus": "Missing"}


class Stop(Exception):
    """Raised by the recorded run() at the first writer command (the test needs the command only)."""


def planner_args(**overrides):
    """The planner options the specs and the builder commands read (a SimpleNamespace stands for parse_args())."""
    base = dict(matching_radius=2.0, orders=[2, 3], mesh_order=1, basis_size=6, samples=12, coupon_depth=1.0,
                edge_offset_tolerance=1.0e-3, straight_lc_fine=0.01, straight_lc_far=0.5, cluster_lc_fine=0.01,
                cluster_lc_far=0.5, corner_lc_fine=0.02, corner_lc_far=0.3, spatial_lc_fine=0.02, spatial_lc_far=0.3,
                spatial_ring_size=8, min_process_feature_elements=2.0, cluster_h_factors=[2.0, 1.0], force=True,
                name="lib", palace=Path("/nowhere/palace"), julia="julia", julia_project=None, ranks=1)
    base.update(overrides)
    return types.SimpleNamespace(**base)


def family_specs(args):
    """The four planner family specs whose fingerprints name the cache roots (the corner spec is built inside
    build_corner: read back from its coupon-spec.json)."""
    edges = PLANNER.normalize_parallel_edges(CLUSTER, args.matching_radius)
    _, corner_spec = corner_generator_command(args)
    return {"straight": PLANNER.straight_spec(STRAIGHT, args, PARAMETERS),
            "parallel-cluster": PLANNER.parallel_cluster_spec(CLUSTER, edges, args, PARAMETERS),
            "corner": corner_spec,
            "spatial": PLANNER.spatial_spec(ENDPOINT, args, PARAMETERS)}


def recorded_run(calls, stop_at):
    def record(command, check=True, cwd=None):
        calls.append([str(value) for value in command])
        if any(str(value).endswith(stop_at) for value in command):
            raise Stop()
        return 0
    return record


def corner_generator_command(args):
    """build_corner's generator command and the spec it wrote to coupon-spec.json (run() recorded, stopped there)."""
    calls = []
    with tempfile.TemporaryDirectory() as directory:
        cache = Path(directory)
        with mock.patch.object(PLANNER, "run", side_effect=recorded_run(calls, "generate_corner_response.py")):
            try:
                PLANNER.build_corner(CORNER, args, PARAMETERS, cache)
            except Stop:
                pass
        spec = PLANNER.load_json(next(cache.glob("corner-*/coupon-spec.json")))
    return calls[-1], spec


def policy_value(command):
    return command[command.index("--crack-policy") + 1] if "--crack-policy" in command else None


class PlannerCacheKeysCarryTheCrackPolicy(unittest.TestCase):
    def test_every_family_key_changes_with_the_policy_and_the_policy_is_the_only_difference(self):
        legacy, default = family_specs(planner_args()), family_specs(planner_args(crack_policy="default"))
        unset = family_specs(planner_args())   # no crack_policy attribute = the shared default
        for family in ("straight", "parallel-cluster", "spatial"):
            self.assertEqual(legacy[family]["CrackPolicy"], "legacy", family)
            self.assertEqual(default[family]["CrackPolicy"], "default", family)
            self.assertEqual(unset[family], legacy[family], family)
            self.assertNotEqual(PLANNER.fingerprint(legacy[family]), PLANNER.fingerprint(default[family]), family)
            stripped = {key: value for key, value in default[family].items() if key != "CrackPolicy"}
            self.assertEqual(stripped, {key: value for key, value in legacy[family].items() if key != "CrackPolicy"}, family)
        # the corner spec (written by build_corner before its generator runs) likewise
        self.assertEqual((legacy["corner"]["CrackPolicy"], default["corner"]["CrackPolicy"]), ("legacy", "default"))
        self.assertNotEqual(PLANNER.fingerprint(legacy["corner"]), PLANNER.fingerprint(default["corner"]))

    def test_tool_fingerprints_include_crack_policy_py(self):
        self.assertEqual(C.CRACK_POLICY_SOURCE.name, "crack_policy.py")
        straight = PLANNER.straight_spec(STRAIGHT, planner_args(), PARAMETERS)
        self.assertEqual(straight["ToolFingerprint"], PLANNER.tool_fingerprint(
            (PLANNER.STRAIGHT_BUILDER, PLANNER.STRAIGHT_QUALIFIER, C.CRACK_POLICY_SOURCE, Path(PLANNER.__file__))))
        self.assertNotEqual(straight["ToolFingerprint"], PLANNER.tool_fingerprint(
            (PLANNER.STRAIGHT_BUILDER, PLANNER.STRAIGHT_QUALIFIER, Path(PLANNER.__file__))))

    def test_an_unknown_policy_is_refused(self):
        with self.assertRaises(ValueError):
            PLANNER.straight_spec(STRAIGHT, planner_args(crack_policy="true"), PARAMETERS)

    def test_planner_passes_the_policy_to_every_writer_it_runs(self):
        for policy in ("legacy", "default"):
            args = planner_args(crack_policy=policy)
            command, _ = corner_generator_command(args)
            self.assertEqual(policy_value(command), policy)
            calls = []
            with tempfile.TemporaryDirectory() as directory:
                with mock.patch.object(PLANNER, "run", side_effect=recorded_run(calls, "build_surface_process_library.py")), \
                        self.assertRaises(Stop):
                    PLANNER.build_straight(STRAIGHT, args, PARAMETERS, Path(directory))
                self.assertEqual(policy_value(calls[-1]), policy)
                calls.clear()
                meshes = {"thin": Path(directory) / "t.msh", "fabricated": Path(directory) / "f.msh"}
                with mock.patch.object(PLANNER, "run", side_effect=recorded_run(calls, "generate_edge_cluster_response.py")), \
                        mock.patch.object(PLANNER, "generate_parallel_cluster_meshes", return_value=(None, meshes)), \
                        self.assertRaises(Stop):
                    PLANNER.build_parallel_cluster(CLUSTER, args, PARAMETERS, Path(directory))
                self.assertEqual(policy_value(calls[-1]), policy)
                calls.clear()
                with mock.patch.object(PLANNER, "run", side_effect=recorded_run(calls, "generate_spatial_response.py")), \
                        self.assertRaises(Stop):
                    PLANNER.build_spatial(ENDPOINT, args, PARAMETERS, Path(directory))
                self.assertEqual(policy_value(calls[-1]), policy)

    def test_cli_flag_defaults_to_the_shared_constant(self):
        with mock.patch.object(sys, "argv", ["prepare", "manifest.json", "--output", "out"]):
            self.assertEqual(PLANNER.parse_args().crack_policy, C.CRACK_POLICY_DEFAULT)
        with mock.patch.object(sys, "argv", ["prepare", "manifest.json", "--output", "out", "--crack-policy", "default"]):
            self.assertEqual(PLANNER.parse_args().crack_policy, "default")


def builder_args(**overrides):
    base = dict(matching_radius=1.9, metal_thickness=0.1, overetch=0.05, sidewall_angle=80.0, top_radius=0.01,
                bottom_radius=0.01, substrate_permittivity=11.47, sa_thickness=0.002, sa_permittivity=4.0,
                ms_thickness=0.002, ms_permittivity=11.47, ma_thickness=0.002, ma_permittivity=10.0, lc_fine=0.002,
                lc_far=0.05, mesh_order=2, basis_size=6, samples=12, order=2, coupon_depth=1.0,
                separation_tolerance=1.0e-3, edge_distances=None, name="t", requirement_hash=[])
    base.update(overrides)
    return types.SimpleNamespace(**base)


class StraightBuilderCarriesTheCrackPolicy(unittest.TestCase):
    def test_generator_commands_pass_the_policy(self):
        for policy in ("legacy", "default"):
            args = builder_args(crack_policy=policy)
            self.assertEqual(policy_value(BUILDER.edge_response_command(args, Path("/tmp/x"))), policy)
            self.assertEqual(policy_value(BUILDER.pair_response_command(args, 2.0, None, Path("/tmp/x"))), policy)
        self.assertEqual(policy_value(BUILDER.edge_response_command(builder_args(), Path("/tmp/x"))), C.CRACK_POLICY_DEFAULT)

    def test_a_solved_coupon_directory_is_not_reused_across_policies(self):
        """build_isolated / build_coupon reuse a directory whose coupon-spec.json equals the spec."""
        for spec_of in (BUILDER.isolated_spec, lambda args: BUILDER.coupon_spec(args, "same-conductor-gap", 2.0)):
            legacy, default = spec_of(builder_args()), spec_of(builder_args(crack_policy="default"))
            self.assertEqual((legacy["CrackPolicy"], default["CrackPolicy"]), ("legacy", "default"))
            self.assertNotEqual(legacy, default)
            self.assertEqual({k: v for k, v in legacy.items() if k != "CrackPolicy"},
                             {k: v for k, v in default.items() if k != "CrackPolicy"})
        self.assertEqual(BUILDER.isolated_spec(builder_args())["GeneratorFingerprint"], BUILDER.generator_fingerprint(
            (BUILDER.EDGE_MESH_GENERATOR, BUILDER.EDGE_RESPONSE_GENERATOR, C.CRACK_POLICY_SOURCE, Path(BUILDER.__file__))))


if __name__ == "__main__":
    unittest.main()
