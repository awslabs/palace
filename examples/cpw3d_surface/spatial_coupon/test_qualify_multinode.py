# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Multi-node stage planning of `coupon-library qualify` (decision 457): the per-stage node
count from the per-node memory model, the node-scaled time estimate, the per-job node
counts in the split and the plans (PBS select = N, the hostfile mpirun, the per-node
admission guard), the MaximumNodesPerJob fail-closed rule, the runner's multi-node helpers,
and the bitwise proof that one-node plans are unchanged (replay_plans on the stored plans
of record when the assessment tree is present).  The loop-end coupon 1b26671c9080 (846
sources, 129.7 M H1 at p4) is the acceptance example."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "qualify"))
import build_plan  # noqa: E402
import estimate_stages  # noqa: E402
import job_split  # noqa: E402
import qualify_library  # noqa: E402
import replay_plans  # noqa: E402
import run_stages  # noqa: E402

ASSESSMENT = Path(os.environ.get("COUPON_ASSESSMENT_ROOT", HERE.parents[3] / "coupon-accuracy-assessment-20260913"))
# The loop-end fab coupon 1b26671c9080 (PBS 57103 build; loopend-first-case/evidence/dryrun).
LOOP_END_COUNTS = {"Vertices": 2162767, "Edges": 12359558, "TriangleFaces": 17373016, "QuadFaces": 1918807,
                   "Tetrahedra": 7802147, "Prisms": 1200537, "Pyramids": 92349}
LOOP_END_SOURCES = 846
# A NodeScaling calibration of the shape fit_node_scaling.py writes (synthetic values: the
# measured block of cost-model.json is the record of truth; these exercise the arithmetic).
# MeasuredNodes 8 so the arithmetic tests exercise every node count of the cap (the credit cap
# of decision 468 (4) has its own test).
NODE_SCALING = {"MeasuredNodes": 8,
                "Memory": {"Worker": {"ReplicatedFraction": 0.20}, "Reducer": {"ReplicatedFraction": 0.10},
                           "LocalEdge": {"ReplicatedFraction": 0.20}},
                "Time": {"WorkerPerSource": {"Exponent": 0.8}, "WorkerNonSource": {"Exponent": 0.3},
                         "ReducerSetup": {"Exponent": 0.3}, "ReducerReduction": {"Exponent": 0.9},
                         "LocalEdge": {"Exponent": 0.0}}}


def model_with_scaling(path=estimate_stages.DEVICE_COST_MODEL):
    return {**estimate_stages.load_cost_model(path), "NodeScaling": json.loads(json.dumps(NODE_SCALING))}


class NodePlanTest(unittest.TestCase):
    """The node planning arithmetic on the decision-64a device model (DEVICE_COST_MODEL: the
    model the loop end's one-node refusal of record was made with, reducer 2,988 GiB node used
    for 846 sources x 129.7 M H1) with a synthetic NodeScaling block."""

    @classmethod
    def setUpClass(cls):
        cls.profile = json.loads((HERE / "qualify" / "cluster-profile.json").read_text())
        cls.model = model_with_scaling()
        cls.plain = estimate_stages.load_cost_model(estimate_stages.DEVICE_COST_MODEL)
        cls.worst = f"{max(cls.plain['PCGFactors']):.1f}"
        cls.layout = qualify_library.stage_layout("le", [4], [3, 5], LOOP_END_SOURCES, 8)
        cls.stages = [(item["EstimateKey"], item["Order"], item["Sources"]) for item in cls.layout if item["Kind"] == "response"]
        local = next(item for item in cls.layout if item["Kind"] == "local-edge")
        cls.local = (local["EstimateKey"], local["Order"], local["Sources"])

    def estimate(self, model=None, profile=None, counts=LOOP_END_COUNTS):
        return estimate_stages.estimate(counts, self.stages, local_edge=self.local, model=model or self.model,
                                        profile=profile or self.profile)

    @property
    def wide(self):
        """The profile with the node cap lifted to 8: under the cost model of record the loop
        end's p4 reducer (2,988 GiB node used) needs 5 r8g nodes at a 0.1 replicated fraction."""
        return {**self.profile, "MaximumNodesPerJob": 8}

    def test_profile_carries_the_multi_node_keys(self):
        self.assertEqual(self.profile["MaximumNodesPerJob"], 4)
        self.assertEqual(self.profile["RanksPerNode"], self.profile["Ranks"])
        self.assertEqual(self.profile["SelectResourcesPerNode"], "ncpus=192:mpiprocs=192")
        self.assertEqual(self.profile["MultiNodeMPIExecArguments"], ["--hostfile", "$PBS_NODEFILE", "--map-by", "ppr:{ranks_per_node}:node"])
        self.assertIn("768 ranks = 4 nodes", self.profile["MaximumNodesPerJobOrigin"])
        self.assertEqual(self.profile["SelectResources"], "select=1:ncpus=192:mpiprocs=192")   # the one-node line of record

    def test_one_node_coupon_is_planned_as_before(self):
        """A coupon that fits one node: every NodesRequired 1, MultiNode false, the Decision
        text of record, no scaled figure (the plan replay below proves the bytes)."""
        counts = self.plain["MeasuredMesh"]["EntityCounts"]
        estimate = self.estimate(counts=counts)
        self.assertEqual(set(estimate["Nodes"]["Required"].values()), {1})
        self.assertFalse(estimate["Nodes"]["MultiNode"])
        self.assertTrue(estimate["Nodes"]["Fits"])
        self.assertEqual(estimate["Nodes"]["Decision"], "every stage fits one node under the one-node rules of record")
        for stage in estimate["Stages"].values():
            self.assertTrue(stage["NodePlan"]["OneNodeFits"])
            self.assertNotIn("Nodes", stage)
            self.assertNotIn("PerNodeUsedGiBEstimateWorker", stage)
        self.assertIsNone(qualify_library.node_assignment(estimate["Nodes"], self.layout)[0])
        # Without a NodeScaling block the one-node estimate is byte-identical in its figures.
        without = self.estimate(model=self.plain, counts=counts)
        for name, stage in estimate["Stages"].items():
            self.assertEqual(stage["ByPCGFactor"], without["Stages"][name]["ByPCGFactor"])
        self.assertEqual(estimate["Decision"], without["Decision"])

    def test_loop_end_fails_closed_at_the_profile_cap_under_the_model_of_record(self):
        """The cost model of record (reducer 2,988 GiB node used for 846 sources x 129.7 M H1)
        puts the loop end's p4 reducer above every guard at 4 nodes: fail closed with the reason."""
        estimate = self.estimate()
        self.assertFalse(estimate["Nodes"]["Fits"])
        self.assertFalse(estimate["Nodes"]["MultiNode"])
        self.assertIsNone(estimate["Nodes"]["Required"]["p4-846"])
        self.assertEqual(estimate["Nodes"]["Required"]["local-edge-p4-8"], 2)
        self.assertIn("MaximumNodesPerJob = 4 nodes", estimate["Stages"]["p4-846"]["NodePlan"]["Reason"])
        self.assertIn("fail closed", estimate["Nodes"]["Decision"])
        self.assertEqual(sorted(estimate["Stages"]["p4-846"]["NodePlan"]["PerNodeUsedGiB"]), ["2", "3", "4"])

    def test_loop_end_stages_need_the_minimum_node_count_that_fits_an_instance(self):
        estimate = self.estimate(profile=self.wide)
        required = estimate["Nodes"]["Required"]
        self.assertTrue(estimate["Nodes"]["MultiNode"])
        self.assertTrue(estimate["Nodes"]["Fits"])
        main = estimate["Stages"]["p4-846"]
        self.assertFalse(main["NodePlan"]["OneNodeFits"])
        self.assertGreater(required["p4-846"], 1)
        # The minimum: N - 1 nodes hold no instance, N nodes hold the chosen one.
        plan = main["NodePlan"]
        nodes = plan["NodesRequired"]
        fraction = self.profile["MemoryFitFraction"]
        margin = estimate["Nodes"]["PerNodeGuardMargin"]
        self.assertEqual(margin, 0.10)   # decision 463: the estimated per-node peak stays under the guard by >= 10 %
        guards = [fraction * item["MemoryGiB"] * (1 - margin) for item in self.profile["Instances"]]
        for smaller in range(2, nodes):
            self.assertGreater(max(plan["PerNodeUsedGiB"][str(smaller)].values()), max(guards))
        self.assertLessEqual(max(plan["PerNodeUsedGiB"][str(nodes)].values()), fraction * plan["Instance"]["MemoryGiB"] * (1 - margin))
        self.assertEqual(plan["Instance"]["PerNodeGuardGiB"], fraction * plan["Instance"]["MemoryGiB"] * (1 - margin))
        self.assertEqual(plan["Instance"]["Nodes"], nodes)
        self.assertEqual(plan["Instance"]["MinimumMemAvailableBytes"], int(fraction * plan["Instance"]["MemoryGiB"] * 1024 ** 3))
        # The per-node formula: one-node figure x (r + (1 - r) / N) per kind.
        r = NODE_SCALING["Memory"]["Reducer"]["ReplicatedFraction"]
        self.assertAlmostEqual(plan["PerNodeUsedGiB"]["2"]["Reducer"], main["NodeUsedGiBEstimateReducer"] * (r + (1 - r) / 2))
        # The p3 control fits one node; the local-edge stage (1,093 GB Palace, 1,122 GiB node
        # used by the model of record) does not: the fixed group's N is its maximum.
        self.assertEqual(required["p3-control-8"], 1)
        self.assertGreater(required["local-edge-p4-8"], 1)
        nodes_, assignment = qualify_library.node_assignment(estimate["Nodes"], self.layout)
        self.assertEqual(nodes_["Main"], required["p4-846"])
        self.assertEqual(nodes_["Fixed"], max(required["p5-control-8"], required["p3-control-8"], required["local-edge-p4-8"]))
        self.assertEqual(assignment["p3-control-8"], nodes_["Fixed"])

    def test_fails_closed_above_maximum_nodes_per_job_and_without_a_calibration(self):
        tight = {**self.profile, "MaximumNodesPerJob": 1}
        estimate = self.estimate(profile=tight)
        self.assertFalse(estimate["Nodes"]["Fits"])
        self.assertIsNone(estimate["Nodes"]["Required"]["p4-846"])
        self.assertIsNone(estimate["Nodes"]["Required"]["local-edge-p4-8"])
        self.assertIn("MaximumNodesPerJob", estimate["Nodes"]["Decision"])
        self.assertIn("fail closed", estimate["Nodes"]["Decision"])
        # No NodeScaling block in the model: a stage that does not fit one node is refused with the reason.
        estimate = self.estimate(model=self.plain)
        self.assertFalse(estimate["Nodes"]["Fits"])
        self.assertIn("NodeScaling", estimate["Stages"]["p4-846"]["NodePlan"]["Reason"])
        with self.assertRaisesRegex(ValueError, "MaximumNodesPerJob"):
            build_plan.multi_node_fields(self.profile, 5, {"MinimumMemAvailableBytes": 1})
        with self.assertRaisesRegex(ValueError, "MaximumNodesPerJob"):
            build_plan.select_resources(self.profile, 5)
        self.assertEqual(build_plan.select_resources(self.profile, 1), self.profile["SelectResources"])
        self.assertEqual(build_plan.select_resources(self.profile, 3), "select=3:ncpus=192:mpiprocs=192")

    def test_scaled_estimate_divides_every_part_by_n_to_the_measured_exponent(self):
        estimate = self.estimate(profile=self.wide)
        nodes, assignment = qualify_library.node_assignment(estimate["Nodes"], self.layout)
        scaled = estimate_stages.scale_to_nodes(estimate, assignment, self.model, self.wide)
        main, main_scaled = estimate["Stages"]["p4-846"], scaled["Stages"]["p4-846"]
        n = nodes["Main"]
        self.assertEqual(main_scaled["Nodes"], n)
        self.assertAlmostEqual(main_scaled["WorkerNonSourceSecondsEstimate"], main["WorkerNonSourceSecondsEstimate"] / n ** 0.3)
        self.assertAlmostEqual(main_scaled["ReducerSecondsEstimateParts"]["Setup"], main["ReducerSecondsEstimateParts"]["Setup"] / n ** 0.3)
        self.assertAlmostEqual(main_scaled["ReducerSecondsEstimateParts"]["Evaluation"], main["ReducerSecondsEstimateParts"]["Evaluation"] / n ** 0.9)
        self.assertAlmostEqual(main_scaled["ReducerSecondsEstimate"], sum(main_scaled["ReducerSecondsEstimateParts"].values()))
        for factor, figures in main["ByPCGFactor"].items():
            per_source = figures["PerSourceSecondsEstimate"] / n ** 0.8
            self.assertAlmostEqual(main_scaled["ByPCGFactor"][factor]["PerSourceSecondsEstimate"], per_source)
            self.assertAlmostEqual(main_scaled["ByPCGFactor"][factor]["WorkerSecondsEstimate"],
                                   main_scaled["WorkerNonSourceSecondsEstimate"] + LOOP_END_SOURCES * per_source)
            self.assertAlmostEqual(main_scaled["ByPCGFactor"][factor]["StageSecondsEstimate"],
                                   main_scaled["ByPCGFactor"][factor]["WorkerSecondsEstimate"] + main_scaled["ReducerSecondsEstimate"])
        self.assertEqual(main_scaled["OneNode"]["ByPCGFactor"], main["ByPCGFactor"])
        self.assertEqual(main_scaled["Instance"]["Nodes"], n)
        self.assertLessEqual(estimate_stages.stage_per_node_used_gib(main_scaled), 0.9 * 0.6 * main_scaled["Instance"]["MemoryGiB"])
        # The local-edge stage scales with the LocalEdge exponent (0: no speedup assumed).
        local, local_scaled = estimate["Stages"]["local-edge-p4-8"], scaled["Stages"]["local-edge-p4-8"]
        self.assertEqual(local_scaled["ByPCGFactor"][self.worst]["StageSecondsEstimate"], local["ByPCGFactor"][self.worst]["StageSecondsEstimate"])
        self.assertLess(local_scaled["PerNodeUsedGiBEstimate"], local["NodeUsedGiBEstimate"])
        # Job totals are recomputed from the scaled stages; the Decision names the nodes.
        self.assertAlmostEqual(scaled["JobSecondsEstimateByPCGFactor"]["1.0"],
                               sum(stage["ByPCGFactor"]["1.0"]["StageSecondsEstimate"] for stage in scaled["Stages"].values()))
        self.assertIn("multi-node coupon (decision 457)", scaled["Decision"])
        self.assertIn(f"p4-846 on {n} node(s)", scaled["Decision"])
        self.assertEqual(scaled["Nodes"]["Assigned"], assignment)
        # An exponent above 1 is capped (no superlinear speedup is budgeted).
        superlinear = json.loads(json.dumps(self.model))
        superlinear["NodeScaling"]["Time"]["WorkerPerSource"]["Exponent"] = 1.4
        self.assertEqual(estimate_stages.time_exponent(superlinear, "WorkerPerSource"), 1.0)
        with self.assertRaisesRegex(ValueError, "needs"):
            estimate_stages.scale_to_nodes(estimate, {**assignment, "p4-846": 1}, self.model, self.wide)
        with self.assertRaisesRegex(ValueError, "MaximumNodesPerJob"):
            estimate_stages.scale_to_nodes(estimate, {**assignment, "p4-846": 9}, self.model, self.wide)

    def split(self, scaled, nodes, mode, max_jobs):
        policy = job_split.normalize_policy(mode, max_jobs=max_jobs, walltime_seconds=self.profile["WalltimeSeconds"],
                                            user_job_cap=self.profile["UserJobCap"])
        return job_split.plan_split(indices=list(range(1, LOOP_END_SOURCES + 1)), layout=self.layout, estimate=scaled, policy=policy,
                                    model=self.model, profile=self.wide, nodes=nodes)

    def test_split_separates_the_fixed_stages_when_their_node_count_differs(self):
        estimate = self.estimate(profile=self.wide)
        nodes, assignment = qualify_library.node_assignment(estimate["Nodes"], self.layout)
        scaled = estimate_stages.scale_to_nodes(estimate, assignment, self.model, self.wide)
        self.assertNotEqual(nodes["Main"], nodes["Fixed"])
        record = self.split(scaled, nodes, "speed", 6)
        self.assertEqual(record["Nodes"], {"Main": nodes["Main"], "Fixed": nodes["Fixed"], "SeparateFixedStages": True,
                                           "Rule": record["Nodes"]["Rule"]})
        self.assertFalse(record["Candidates"][0]["Fits"])   # one job cannot hold two node counts
        self.assertTrue(record["Fits"])
        self.assertEqual(record["ControlsJob"], "separate")
        self.assertEqual(record["Blocks"][0], [])
        self.assertEqual(sum(len(block) for block in record["Blocks"]), LOOP_END_SOURCES)
        jobs = {job["Name"]: job for job in record["Jobs"]}
        self.assertEqual(jobs["worker-1"]["Nodes"], nodes["Fixed"])
        for name, job in jobs.items():
            if name != "worker-1":
                self.assertEqual(job["Nodes"], nodes["Main"])
        worst = record["WorstPCGFactor"]
        self.assertAlmostEqual(record["NodeSecondsEstimate"][worst],
                               sum(job["SecondsEstimateWithPreflightAndMargin"][worst] * job["Nodes"] for job in record["Jobs"]))
        # Equal node counts: the fixed stages share job 1 as in the recorded planner.
        same = {"Main": nodes["Main"], "Fixed": nodes["Main"]}
        same_assignment = {key: nodes["Main"] for key in assignment}
        same_scaled = estimate_stages.scale_to_nodes(estimate, same_assignment, self.model, self.wide)
        record = self.split(same_scaled, same, "frugal", 6)
        self.assertEqual(record["ControlsJob"], "worker-1")
        self.assertTrue(all(job["Nodes"] == nodes["Main"] for job in record["Jobs"]))
        # The one-node planner (nodes None) writes no Nodes key: byte-identical records.
        plain = self.split(estimate, None, "frugal", 6)
        self.assertNotIn("Nodes", plain)
        self.assertTrue(all("Nodes" not in job for job in (plain["Jobs"] or [])))

    def test_multi_node_plans_carry_nodes_ranks_launch_arguments_and_the_guard(self):
        estimate = self.estimate(profile=self.wide)
        nodes, assignment = qualify_library.node_assignment(estimate["Nodes"], self.layout)
        scaled = estimate_stages.scale_to_nodes(estimate, assignment, self.model, self.wide)
        record = self.split(scaled, nodes, "speed", 6)
        digests = {"le-p4": {f"worker-block{k}.json": f"{k}" * 64 for k in range(1, record["N"] + 1)} | {"reducer.json": "r" * 64},
                   "le-p5-control": {"worker.json": "a" * 64, "reducer.json": "b" * 64},
                   "le-p3-control": {"worker.json": "c" * 64, "reducer.json": "d" * 64},
                   "le-p4-local-edge": {"config.json": "e" * 64}}
        common = dict(case_id="le", remote_case_root="/r/case", mesh={"Remote": "/r/case/mesh/m.msh", "SHA256": "m" * 64, "Local": "/l/m.msh"},
                      stage_layout=self.layout, estimate=scaled, config_digests=digests, trace_pins={}, profile=self.wide,
                      binary="/r/b.bin", binary_sha256="0" * 64, mpiexec="/r/mpiexec_bound.sh", purpose="test",
                      factors=[f"{factor:.1f}" for factor in self.plain["PCGFactors"]])
        plans = {job["Name"]: build_plan.build_job_plan(job_name=job["Name"], split_job=job, **common) for job in record["Jobs"]}
        reducer = plans["reducer"]
        self.assertEqual(reducer["Nodes"], nodes["Main"])
        self.assertEqual(reducer["RanksPerNode"], 192)
        self.assertEqual(reducer["Ranks"], 192 * nodes["Main"])
        self.assertEqual(reducer["MPIExecArguments"], ["--hostfile", "$PBS_NODEFILE", "--map-by", "ppr:192:node"])
        self.assertEqual(reducer["NodeGuard"]["MinimumMemAvailableBytes"], reducer["MinimumMemAvailableBytes"])
        self.assertEqual(reducer["Instance"]["Nodes"], nodes["Main"])
        self.assertLessEqual(reducer["Instance"]["PerNodeUsedGiBEstimate"], 0.9 * 0.6 * reducer["Instance"]["MemoryGiB"])
        self.assertEqual(reducer["Instance"]["PerNodeUsedGiBEstimate"], scaled["Stages"]["p4-846"]["PerNodeUsedGiBEstimateReducer"])
        self.assertEqual(plans["worker-1"]["Nodes"], nodes["Fixed"])
        self.assertEqual(plans["worker-1"]["StageNames"], ["le-p5-control-worker", "le-p5-control-reducer", "le-p3-control-worker",
                                                           "le-p3-control-reducer", "le-p4-local-edge"])
        self.assertEqual(plans["worker-2"]["Nodes"], nodes["Main"])
        self.assertEqual(plans["worker-2"]["Ranks"], reducer["Ranks"])   # the archive of every block is reducible
        for plan in plans.values():
            self.assertEqual(plan["Version"], build_plan.PLAN_VERSION)
            self.assertEqual(plan["PBSDsh"], self.profile["PBSDsh"])
            self.assertTrue(plan["Purpose"].endswith(build_plan.MULTI_NODE_REDUCTION_NOTE))   # decision 466
            self.assertEqual(plan["MultiNodeReductionNote"], build_plan.MULTI_NODE_REDUCTION_NOTE)
        script = build_plan.render_job_script(profile=self.wide, remote_root="/r", remote_case_root="/r/case", runner="/r/run/run_stages.py",
                                              job_name="j", walltime_seconds=21600, instance_type=reducer["Instance"]["Type"],
                                              job_directory="/r/case/main/jobs/reducer", nodes=reducer["Nodes"])
        self.assertIn(f"#PBS -l select={nodes['Main']}:ncpus=192:mpiprocs=192\n", script)
        self.assertIn(f"#PBS -l instance_type={reducer['Instance']['Type']}\n", script)
        self.assertIn("#PBS -l efa_support=True,subnet_id=subnet-0c98d793bbcebb39a\n", script)
        # A one-node job plan has none of the multi-node keys and the one-node select line.
        one = build_plan.build_job_plan(job_name="w", split_job={**record["Jobs"][1], "Nodes": 1}, **{**common, "estimate": estimate})
        for key in ("Nodes", "RanksPerNode", "MPIExecArguments", "NodeGuard", "PBSDsh", "MultiNodeRule", "MultiNodeReductionNote"):
            self.assertNotIn(key, one)
        self.assertEqual(one["Purpose"], "test")
        self.assertEqual(one["Ranks"], 192)
        one_script = build_plan.render_job_script(profile=self.profile, remote_root="/r", remote_case_root="/r/case", runner="/r/run/run_stages.py",
                                                  job_name="j", walltime_seconds=21600, instance_type="m8g.48xlarge")
        self.assertIn("#PBS -l select=1:ncpus=192:mpiprocs=192\n", one_script)
        # The whole-coupon plan at one node count (the single job of a coupon whose groups agree).
        same_assignment = {key: nodes["Main"] for key in assignment}
        same_scaled = estimate_stages.scale_to_nodes(estimate, same_assignment, self.model, self.wide)
        whole = build_plan.build_plan(**{**common, "estimate": same_scaled}, nodes=nodes["Main"])
        self.assertEqual((whole["Nodes"], whole["Ranks"]), (nodes["Main"], 192 * nodes["Main"]))
        self.assertEqual(whole["Purpose"], f"test; {build_plan.MULTI_NODE_REDUCTION_NOTE}")   # decision 468 (2): the single-job path too
        self.assertIn("<= 8e-6 per interface of the column maximum, per-Type <= 6e-6", build_plan.MULTI_NODE_REDUCTION_NOTE)
        self.assertIn("up to 4 nodes; decisions 466 / 468 / 482", build_plan.MULTI_NODE_REDUCTION_NOTE)
        self.assertEqual(whole["Instance"]["PerNodeUsedGiBEstimate"],
                         max(estimate_stages.stage_per_node_used_gib(stage) for stage in same_scaled["Stages"].values()))

    def test_multi_node_reduction_note_is_the_measured_floor(self):
        """Decision 490: the note carries the MEASURED cross-partition floor of record - the pair-5 2-node
        figures (decisions 466 / 468) and the loop end's MS 5.3e-6 per Type at 4 vs 1 nodes (decision 482,
        cross-partition-floor.json) - as a bound above every measurement, with its provenance; the former
        "per-Type sums <= 3e-6" sat below the 4-node measurement."""
        measurements = build_plan.MULTI_NODE_REDUCTION_MEASUREMENTS
        self.assertEqual([item["Nodes"] for item in measurements], [2, 4])
        self.assertEqual(build_plan.MULTI_NODE_REDUCTION_MAX_NODES, 4)
        per_type = max(item["PerType"] for item in measurements)
        per_interface = max(item["PerInterfaceOfColumnMaximum"] for item in measurements if item["PerInterfaceOfColumnMaximum"] is not None)
        self.assertAlmostEqual(per_type, 5.3e-6)
        self.assertAlmostEqual(per_interface, 7.1e-6)
        self.assertLessEqual(per_type, 6e-6)
        self.assertLessEqual(per_interface, 8e-6)
        self.assertGreater(per_type, 3e-6)   # the superseded per-Type figure
        self.assertNotIn("3e-6", build_plan.MULTI_NODE_REDUCTION_NOTE)
        for item in measurements:
            self.assertTrue(item["Coupon"] and item["Quantity"] and item["Decisions"])
        floor = LOOP_END_F / "cross-partition-floor.json"
        if floor.is_file():
            recorded = json.loads(floor.read_text())["MaxAbsPerType"]
            self.assertAlmostEqual(max(recorded.values()), measurements[1]["PerType"], delta=5e-8)
            self.assertLessEqual(max(recorded.values()), 6e-6)


class SpeedupCreditTest(unittest.TestCase):
    """Decision 468 (4): beyond the largest MEASURED node count the time estimate takes no further
    speedup (the caps stay conservative); the memory model keeps the replicated-fraction law."""

    def test_time_credit_stops_at_the_measured_node_count(self):
        model = {**estimate_stages.load_cost_model(estimate_stages.DEVICE_COST_MODEL),
                 "NodeScaling": {**json.loads(json.dumps(NODE_SCALING)), "MeasuredNodes": 2}}
        self.assertEqual([estimate_stages.speedup_nodes(model, n) for n in (1, 2, 3, 4, 8)], [1, 2, 2, 2, 2])
        two = estimate_stages.scaled_seconds(model, "WorkerPerSource", 100.0, 2)
        self.assertAlmostEqual(two, 100.0 / 2 ** 0.8)
        self.assertEqual(estimate_stages.scaled_seconds(model, "WorkerPerSource", 100.0, 4), two)
        # Memory still divides by the real node count.
        self.assertAlmostEqual(estimate_stages.per_node_used_gib(model, "Reducer", 1000.0, 4), 1000.0 * (0.1 + 0.9 / 4))
        self.assertLess(estimate_stages.per_node_used_gib(model, "Reducer", 1000.0, 4), estimate_stages.per_node_used_gib(model, "Reducer", 1000.0, 2))
        stage = {"Order": 4, "Sources": 10, "NodeUsedGiBEstimateWorker": 100.0, "NodeUsedGiBEstimateReducer": 100.0, "WorkerPalacePeakGBEstimate": 90.0,
                 "ReducerPalacePeakGBEstimate": 90.0, "WorkerNonSourceSecondsEstimate": 100.0, "ReducerSecondsEstimate": 300.0,
                 "ReducerSecondsEstimateParts": {"Setup": 100.0, "Evaluation": 100.0, "Gram": 100.0},
                 "ByPCGFactor": {"1.0": {"MeanPCGIterations": 10, "PerSourceSecondsEstimate": 50.0, "WorkerSecondsEstimate": 600.0, "StageSecondsEstimate": 900.0}}}
        at_four = estimate_stages.scale_stage_to_nodes(model, stage, 4)
        at_two = estimate_stages.scale_stage_to_nodes(model, stage, 2)
        self.assertEqual(at_four["SpeedupNodes"], 2)
        self.assertEqual(at_four["ByPCGFactor"]["1.0"]["StageSecondsEstimate"], at_two["ByPCGFactor"]["1.0"]["StageSecondsEstimate"])
        self.assertLess(at_four["PerNodeUsedGiBEstimateReducer"], at_two["PerNodeUsedGiBEstimateReducer"])
        # The committed model measured 2 nodes: a 4-node stage is credited 2 nodes of speedup.
        committed = estimate_stages.load_cost_model()
        self.assertEqual(committed["NodeScaling"]["MeasuredNodes"], 2)
        self.assertEqual(estimate_stages.speedup_nodes(committed, 4), 2)


class JobNodeHoursTest(unittest.TestCase):
    """Decision 468 (3): the coupon's node-hours sum every job at ITS node count; every stage is
    charged the nodes of the job that ran it."""

    def test_per_job_node_hours_and_stage_nodes(self):
        import summarize_cost
        jobs = [{"Name": "worker-1", "Nodes": 2, "StageNames": ["c-p5-control-worker", "c-p5-control-reducer", "c-p4-local-edge"]},
                {"Name": "worker-2", "Nodes": 4, "StageNames": ["c-p4-worker-block2"]},
                {"Name": "reducer", "Nodes": 4, "StageNames": ["c-p4-reducer"]}]
        statuses = {"worker-1": {"TotalSeconds": 3600.0}, "worker-2": {"TotalSeconds": 1800.0}, "reducer": {"TotalSeconds": 900.0}}
        per_job, total = summarize_cost.job_node_hours(statuses, jobs, 1)
        self.assertEqual(per_job, {"worker-1": 2.0, "worker-2": 2.0, "reducer": 1.0})
        self.assertEqual(total, 5.0)   # not 4 x (3600 + 1800 + 900) / 3600 = 7
        nodes = summarize_cost.stage_nodes(jobs, 1)
        self.assertEqual(nodes["c-p4-worker"], 4)   # the merged block worker takes its blocks' job count
        self.assertEqual((nodes["c-p5-control-worker"], nodes["c-p4-local-edge"], nodes["c-p4-reducer"]), (2, 2, 4))
        # One-node jobs without a Nodes key take the default.
        per_job, total = summarize_cost.job_node_hours({"single": {"TotalSeconds": 7200.0}}, [{"Name": "single", "StageNames": []}], 1)
        self.assertEqual((per_job, total), ({"single": 2.0}, 2.0))
        # summarize charges the stages per their job's nodes.
        status = {"TotalSeconds": 6300.0, "Stages": [{"Name": "c-p4-local-edge", "State": "complete", "WallSeconds": 360.0, "Parsed": {}}]}
        summary = summarize_cost.summarize(status, full_sources=10, nodes=4, nodes_by_stage=nodes)
        self.assertEqual(summary["Stages"]["c-p4-local-edge"]["NodeHours"], 0.2)   # 360 s x 2 nodes


class RunnerMultiNodeTest(unittest.TestCase):
    def test_hosts_conflicts_launch_arguments_and_per_node_peaks(self):
        nodefile = "\n".join(["ip-10-0-0-1"] * 3 + ["ip-10-0-0-2"] * 3 + ["ip-10-0-0-1"]) + "\n"
        self.assertEqual(run_stages.unique_hosts(nodefile), ["ip-10-0-0-1", "ip-10-0-0-2"])
        self.assertEqual(run_stages.conflicting_processes("1 bash\n2 palace-x\n3 prterun\n4 python3\n"), ["2 palace-x", "3 prterun"])
        with tempfile.TemporaryDirectory() as tmp:
            # No pbsdsh anywhere: the node shell is ssh in batch mode and the command's output
            # is read back from the shared job directory.
            saved = os.environ.get("PATH")
            os.environ["PATH"] = "/nonexistent-dir"
            nodefile_text = "\n".join(["h1"] * 3 + ["h2"] * 3) + "\n"
            try:
                shell = run_stages.NodeShell({"PBSDsh": "/nonexistent/pbsdsh"}, tmp, ["h2"], nodefile_text)
            finally:
                os.environ["PATH"] = saved
            self.assertIsNone(shell.pbsdsh)
            self.assertEqual(shell.prefix("h2")[:3], ["ssh", "-o", "BatchMode=yes"])
            self.assertEqual(run_stages.first_slot_index(nodefile_text), {"h1": 0, "h2": 3})
            with self.assertRaisesRegex(SystemExit, "no slot"):
                run_stages.NodeShell({"PBSDsh": "/nonexistent/pbsdsh"}, tmp, ["h3"], nodefile_text)
            # A fake pbsdsh (PBS 23: `-n <task slot index>`, the slot = the host's first line in
            # PBS_NODEFILE; the task's output routed to the job's stdout): a command's output is
            # read from the file it writes on the shared job directory, the hostname verified.
            pbsdsh = Path(tmp) / "pbsdsh"
            pbsdsh.write_text("#!/bin/sh\n"
                              "# fake pbsdsh -n IDX -- program: slot 3 is h2, every other slot h1\n"
                              'idx=$2; shift 3; if [ "$idx" = 3 ]; then HOSTNAME_FAKE=h2 "$@"; else HOSTNAME_FAKE=h1 "$@"; fi\n')
            pbsdsh.chmod(0o755)
            fake_bin = Path(tmp) / "bin"
            fake_bin.mkdir()
            (fake_bin / "hostname").write_text("#!/bin/sh\necho ${HOSTNAME_FAKE:-h2}\n")
            (fake_bin / "hostname").chmod(0o755)
            os.environ["PATH"] = f"{fake_bin}:{saved}"
            try:
                shell = run_stages.NodeShell({"PBSDsh": str(pbsdsh)}, tmp, ["h2"], nodefile_text)
                self.assertEqual(shell.index, {"h1": 0, "h2": 3})
                self.assertEqual(shell.prefix("h2"), [str(pbsdsh), "-n", "3", "--"])
                text = shell.run("h2", "hostname; echo ---MEM---; echo 'MemTotal: 10 kB'; echo ---PS---; echo '1 bash'", Path(tmp) / "out.txt")
                self.assertIn("---PS---", text)
                with self.assertRaisesRegex(SystemExit, "landed on"):
                    run_stages.remote_node_preflight(run_stages.NodeShell({"PBSDsh": str(pbsdsh)}, tmp, ["h1", "h2"], "\n".join(["h2"] * 3 + ["h1"] * 3)), "h1", tmp)
                if Path("/proc/meminfo").exists():
                    values, processes = run_stages.remote_node_preflight(shell, "h2", tmp)
                    self.assertGreater(values["MemAvailableBytes" if "MemAvailableBytes" in values else "MemAvailable"], 0)
                    self.assertIn("bash", processes + "bash")
                else:
                    with self.assertRaisesRegex(SystemExit, "no MemTotal"):
                        run_stages.remote_node_preflight(shell, "h2", tmp)
            finally:
                os.environ["PATH"] = saved
            samples = Path(tmp) / "memory-samples-h2.csv"
            samples.write_text("unix,host,used_bytes\n100,h2,5\n110,h2,9\n130,h2,7\n")
            self.assertEqual(run_stages.per_node_peaks([samples], 105, 120), {"h2": 9})
            self.assertEqual(run_stages.per_node_peaks([samples], 125, 140), {"h2": 7})
            self.assertEqual(run_stages.per_node_peaks([Path(tmp) / "absent.csv"], 0, 200), {})
            os.environ["PBS_NODEFILE"] = str(Path(tmp) / "nodes")
            try:
                arguments = [os.path.expandvars(a) for a in ["--hostfile", "$PBS_NODEFILE", "--map-by", "ppr:192:node"]]
            finally:
                del os.environ["PBS_NODEFILE"]
            self.assertEqual(arguments, ["--hostfile", str(Path(tmp) / "nodes"), "--map-by", "ppr:192:node"])

    @unittest.skipUnless(Path("/proc/meminfo").exists(), "the node sampler reads /proc/meminfo (Linux)")
    def test_node_sampler_mode_writes_until_the_stop_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            csv_path, stop = Path(tmp) / "samples.csv", Path(tmp) / "stop"
            process = subprocess.Popen([sys.executable, str(HERE / "qualify" / "run_stages.py"), "--node-sampler", str(csv_path), str(stop)])
            time.sleep(1.5)
            stop.write_text("now\n")
            self.assertEqual(process.wait(timeout=30), 0)
            lines = csv_path.read_text().splitlines()
            self.assertEqual(lines[0], "unix,host,used_bytes")
            self.assertGreaterEqual(len(lines), 2)


LOOP_END_RECORDS = ASSESSMENT / "curved-clusters-20261005" / "loopend-first-case" / "records"
LOOP_END_F = ASSESSMENT / "curved-clusters-20261005" / "loopend-first-case" / "f-qualification" / "spatial-38-edge-f0461584cccf"
STORED_RUNS = [ASSESSMENT / "stage2-20261004" / "coupons-a" / "qualify" / "sct002-S2p" / "fab-spatial-3-edge-2c54db92028f" / "library-qualification.json",
               ASSESSMENT / "stage2-20261004" / "coupons-bc" / "qualify" / "ctx003-C3" / "fab-spatial-19-edge-12b9d5c5c1bf" / "library-qualification.json"]


class ReplayStoredPlansTest(unittest.TestCase):
    """Decision 457 (1) / 458: the plans of record of stored one-node coupons (a single-job
    coupon and a two-worker split) are byte-identical under the kept cost model they planned
    with (the digest the recorded estimate names)."""

    @unittest.skipUnless(all(path.is_file() for path in STORED_RUNS), "the stage-2 qualify records are not present")
    def test_stored_one_node_plans_are_byte_identical(self):
        for path in STORED_RUNS:
            case = json.loads(path.read_text())["Cases"][0]
            estimate = json.loads((path.parent / Path(case["Root"]).name / "preflight" / "stage-estimate.json").read_text())
            digest = estimate["CostModel"]["SHA256"]
            model_path = next((candidate for candidate in (HERE / "qualify").glob("cost-model*.json")
                               if replay_plans.sha256(candidate) == digest), None)
            self.assertIsNotNone(model_path, f"no kept cost model with digest {digest[:12]} in qualify/")
            record = replay_plans.replay(path, cost_model_path=model_path)
            self.assertTrue(record["AllIdentical"], json.dumps(record, indent=1)[:3000])
            self.assertEqual(record["PlansIdentical"], record["Plans"])
            for case_record in record["Cases"]:
                self.assertFalse(case_record["Estimate"]["MultiNode"])


if __name__ == "__main__":
    unittest.main()


class DenseTwinPlanTest(unittest.TestCase):
    """Decision 457 (4): the (F) dense twins as one run_stages plan on the minimum node count
    (dense_twin_plan on the committed dense-twin model; the pair-5 twins reproduce themselves)."""

    @classmethod
    def setUpClass(cls):
        import dense_twin_plan
        cls.dense_twin_plan = dense_twin_plan
        cls.model = model_with_scaling(estimate_stages.COST_MODEL)
        cls.dense_model = json.loads(dense_twin_plan.DENSE_TWIN_MODEL.read_text())
        cls.profile = json.loads((HERE / "qualify" / "cluster-profile.json").read_text())

    def test_dense_twin_model_reproduces_the_measured_pair5_twins(self):
        self.assertEqual(sorted(self.dense_model["Orders"]), ["p4", "p5"])
        for name, run in self.dense_model["Runs"].items():
            rates = self.dense_model["Orders"][f"p{run['Order']}"]
            self.assertGreaterEqual(rates["PalaceGBPerMillionH1"] * run["H1"] / 1e6, run["PalacePeakGB"] - 1e-6, name)
            self.assertGreaterEqual((rates["NonSolveSecondsPerMillionH1"] + rates["SolveSecondsPerMillionH1PerTrace"] * run["Traces"]) * run["H1"] / 1e6,
                                    run["PalaceTotalSeconds"] - 1e-6, name)
            self.assertEqual(run["Traces"], 5)
        self.assertEqual(self.dense_model["Executable"], "cc7c4091fa47bde8739fe57678d1b72a22a4aee51b0de8da18ccb9912a5aec88")

    def test_plan_puts_the_runs_on_the_minimum_node_count_and_pins_the_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mesh = root / "identity.msh"
            mesh.write_text("mesh\n")
            configs = {}
            for name, order in (("fabricated-p4", 4), ("fabricated-p5", 5)):
                traces = []
                for k in range(5):
                    trace = root / f"trace-{k}.csv"
                    trace.write_text(f"{k}\n")
                    traces.append({"Index": 100 + k, "DataFile": str(trace)})
                config = {"Problem": {"Type": "Electrostatic", "Output": str(root / name / "postpro")}, "Model": {"Mesh": str(mesh)},
                          "Solver": {"Order": order}, "Boundaries": {"PrescribedPotential": traces}}
                (root / name).mkdir()
                (root / name / "config.json").write_text(json.dumps(config))
                configs[name] = str(root / name / "config.json")
            (root / "dense-traces.json").write_text(json.dumps({"Version": 1, "Model": "le", "Configs": configs, "Traces": []}))
            profile = {**self.profile, "MaximumNodesPerJob": 8}
            plan_kwargs = dict(dense_dir=root, counts=LOOP_END_COUNTS, binary_sha256="0" * 64, remote_root="/r", job_dir="/r/dense",
                               model=self.model, dense_model=self.dense_model, profile=profile, case_id="le")
            # Decisions 482 / 487 (d): a job holding the identity twin (fabricated-p4) takes the reducer's
            # partition: without it the planner fails closed; the p5 twin alone plans as before.
            with self.assertRaisesRegex(ValueError, "identity twin"):
                self.dense_twin_plan.plan_dense_twins(runs=["fabricated-p4", "fabricated-p5"], **plan_kwargs)
            reducer = {"Nodes": 4, "Ranks": 768, "Job": "reducer", "PBSJobID": "57910.fake", "Case": "le", "Record": "/r/lq.json"}
            plan, script = self.dense_twin_plan.plan_dense_twins(runs=["fabricated-p4", "fabricated-p5"], reducer_partition=reducer,
                                                                 **plan_kwargs)
            with self.assertRaisesRegex(ValueError, "below"):
                self.dense_twin_plan.plan_dense_twins(runs=["fabricated-p5"], nodes=1, **plan_kwargs)
            alone, _ = self.dense_twin_plan.plan_dense_twins(runs=["fabricated-p5"], **plan_kwargs)
            # The identity twin on the reducer's partition exactly: a differing --nodes, a reducer whose
            # rank count this profile cannot launch, or a twin needing more nodes than the reducer fail closed.
            with self.assertRaisesRegex(ValueError, "differs from the reducer"):
                self.dense_twin_plan.plan_dense_twins(runs=["fabricated-p4"], reducer_partition=reducer, nodes=2, **plan_kwargs)
            with self.assertRaisesRegex(ValueError, "cannot take the reducer's partition here"):
                self.dense_twin_plan.plan_dense_twins(runs=["fabricated-p4"], reducer_partition={**reducer, "Ranks": 384}, **plan_kwargs)
            with self.assertRaisesRegex(ValueError, "no tolerance is widened"):
                self.dense_twin_plan.plan_dense_twins(runs=["fabricated-p4", "fabricated-p5"], reducer_partition={**reducer, "Nodes": 1, "Ranks": 192},
                                                      **plan_kwargs)
            one_node, _ = self.dense_twin_plan.plan_dense_twins(runs=["fabricated-p4"], reducer_partition={**reducer, "Nodes": 1, "Ranks": 192},
                                                                **plan_kwargs)
        runs = plan["Estimate"]["Runs"]
        # The loop-end fab p4 twin (887 GB Palace) fits one node; the p5 twin (1,702 GB) does not: alone
        # it runs on its own node count; with the identity twin the job takes the reducer's 4 nodes.
        self.assertEqual(runs["fabricated-p4"]["NodePlan"]["NodesRequired"], 1)
        self.assertGreater(runs["fabricated-p5"]["NodePlan"]["NodesRequired"], 1)
        self.assertEqual(alone["Nodes"], alone["Estimate"]["Runs"]["fabricated-p5"]["NodePlan"]["NodesRequired"])
        self.assertIsNone(alone["IdentityPartition"])
        self.assertEqual((plan["Nodes"], plan["Ranks"]), (4, 768))
        self.assertEqual(plan["IdentityPartition"]["Run"], "fabricated-p4")
        self.assertEqual((plan["IdentityPartition"]["Nodes"], plan["IdentityPartition"]["Ranks"]), (4, 768))
        self.assertEqual(plan["IdentityPartition"]["Reducer"], reducer)
        self.assertTrue(plan["IdentityPartition"]["MatchesReducer"])
        self.assertEqual(plan["IdentityPartition"]["Rule"], self.dense_twin_plan.IDENTITY_PARTITION_RULE)
        self.assertIn("57910.fake", plan["Purpose"])
        self.assertNotIn("Nodes", one_node)   # a one-node reducer: the one-node plan, byte-identical keys
        self.assertEqual((one_node["Ranks"], one_node["IdentityPartition"]["Nodes"]), (192, 1))
        self.assertEqual(plan["Ranks"], 192 * plan["Nodes"])
        self.assertEqual(plan["Instance"]["Nodes"], plan["Nodes"])
        self.assertLessEqual(plan["Estimate"]["PerNodeUsedGiB"], 0.9 * 0.6 * plan["Instance"]["MemoryGiB"])
        self.assertEqual([stage["Name"] for stage in plan["Stages"]], ["fabricated-p4-dense", "fabricated-p5-dense"])
        for stage in plan["Stages"]:
            self.assertEqual(stage["Environment"], {})   # an ordinary Palace run
            self.assertGreaterEqual(stage["CapSeconds"], stage["MinimumSeconds"])
        self.assertEqual(len(plan["PinnedSHA256"]), 2 + 1 + 5)   # two configs, the mesh, the five shared traces
        self.assertEqual(plan["UnpinnedInputs"], [])
        self.assertIn(f"#PBS -l select={plan['Nodes']}:ncpus=192:mpiprocs=192\n", script)
        self.assertIn("D=/r/dense\n", script)
        with self.assertRaisesRegex(ValueError, "no measured order"):
            self.dense_twin_plan.estimate_run(self.model, self.dense_model, LOOP_END_COUNTS, 3, 5, profile=self.profile)

    def test_reducer_partition_is_read_from_the_qualify_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            record = Path(tmp) / "library-qualification.json"
            split = {"Case": "le", "Jobs": [{"Name": "worker-1", "Kind": "worker", "Nodes": 2, "Ranks": 384, "Submission": {"Job": "1.h"}},
                                            {"Name": "reducer", "Kind": "reducer", "Nodes": 4, "Ranks": 768, "Submission": {"Job": "57910.h"}}]}
            single = {"Case": "s", "Jobs": [{"Name": "single", "Kind": "single", "Submission": {"Job": "2.h"}}]}
            record.write_text(json.dumps({"Cases": [split, single]}))
            partition = self.dense_twin_plan.reducer_partition_of(record, "le", profile=self.profile)
            self.assertEqual({key: partition[key] for key in ("Nodes", "Ranks", "Job", "PBSJobID", "Case")},
                             {"Nodes": 4, "Ranks": 768, "Job": "reducer", "PBSJobID": "57910.h", "Case": "le"})
            partition = self.dense_twin_plan.reducer_partition_of(record, "s", profile=self.profile)
            self.assertEqual((partition["Nodes"], partition["Ranks"], partition["Job"]), (1, 192, "single"))
            with self.assertRaisesRegex(ValueError, "name the coupon"):
                self.dense_twin_plan.reducer_partition_of(record, None, profile=self.profile)
            with self.assertRaisesRegex(ValueError, "no case"):
                self.dense_twin_plan.reducer_partition_of(record, "x", profile=self.profile)
            record.write_text(json.dumps({"Cases": [{"Case": "p", "Jobs": []}]}))
            with self.assertRaisesRegex(ValueError, "no reducer"):
                self.dense_twin_plan.reducer_partition_of(record, None, profile=self.profile)

    @unittest.skipUnless((LOOP_END_F / "dense" / "fabricated-p4-4n-plan.json").is_file() and (LOOP_END_RECORDS / "qualify-fab" / "library-qualification.json").is_file(),
                         "the loop-end (F) records are not present")
    def test_loop_end_identity_twin_reproduces_the_decision_482_four_node_plan(self):
        """The acceptance example of decisions 482 / 487 (d): the loop-end fab p4 identity twin, first
        planned on 1 node (PBS 58012, identity MS 1.6-5.3e-6 against the 4-node Gram) and re-solved by
        the ruling on the reducer's 4 x r8g partition (PBS 58043, `--nodes 4`, identity <= 6.0e-8), is
        now planned on that partition from the qualify record itself: the same Nodes / Ranks / instance /
        stage caps as the stored 4-node plan, with IdentityPartition recorded."""
        stored = json.loads((LOOP_END_F / "dense" / "fabricated-p4-4n-plan.json").read_text())
        qualify_record = LOOP_END_RECORDS / "qualify-fab" / "library-qualification.json"
        case_id = json.loads(qualify_record.read_text())["Cases"][0]["Case"]
        partition = self.dense_twin_plan.reducer_partition_of(qualify_record, case_id, profile=self.profile)
        self.assertEqual((partition["Nodes"], partition["Ranks"], partition["Job"]), (4, 768, "reducer"))
        build = json.loads((LOOP_END_RECORDS / "library-build.json").read_text())
        counts = next(case["H1"]["EntityCounts"] for case in build["Cases"] if case["Case"] == case_id)
        dense = json.loads((LOOP_END_F / "dense" / "dense-traces.json").read_text())
        with tempfile.TemporaryDirectory() as tmp:
            # The stored dense directory relocated to this mirror (its configs are the recorded bytes).
            relocated = Path(tmp) / "dense"
            relocated.mkdir()
            configs = {name: str(LOOP_END_F / "dense" / name / "config.json") for name in dense["Configs"]}
            (relocated / "dense-traces.json").write_text(json.dumps({**dense, "Configs": configs}))
            model = estimate_stages.load_cost_model()
            plan, script = self.dense_twin_plan.plan_dense_twins(
                dense_dir=relocated, runs=["fabricated-p4"], counts=counts, binary_sha256=stored["BinarySHA256"],
                remote_root=stored["MPIExec"].rsplit("/", 1)[0], job_dir=stored["Stages"][0]["Config"].rsplit("/", 2)[0], model=model,
                dense_model=self.dense_model, profile=self.profile, case_id=case_id, reducer_partition=partition)
        self.assertEqual((plan["Nodes"], plan["Ranks"], plan["Instance"]["Type"]), (stored["Nodes"], stored["Ranks"], stored["Instance"]["Type"]))
        self.assertEqual((plan["Nodes"], plan["Ranks"]), (4, 768))
        self.assertEqual([(stage["Name"], stage["CapSeconds"], stage["MinimumSeconds"]) for stage in plan["Stages"]],
                         [(stage["Name"], stage["CapSeconds"], stage["MinimumSeconds"]) for stage in stored["Stages"]])
        self.assertEqual(plan["Estimate"]["NodesRequired"], stored["Estimate"]["NodesRequired"])
        self.assertAlmostEqual(plan["Estimate"]["PerNodeUsedGiB"], stored["Estimate"]["PerNodeUsedGiB"], places=6)
        self.assertEqual(plan["MinimumMemAvailableBytes"], stored["MinimumMemAvailableBytes"])
        self.assertEqual(plan["PinnedSHA256"][configs["fabricated-p4"]], stored["PinnedSHA256"][dense["Configs"]["fabricated-p4"]])
        self.assertEqual(plan["IdentityPartition"]["Reducer"]["PBSJobID"], partition["PBSJobID"])
        self.assertEqual((plan["IdentityPartition"]["Nodes"], plan["IdentityPartition"]["Ranks"]), (4, 768))
        self.assertIn("#PBS -l select=4:ncpus=192:mpiprocs=192\n", script)


class FitNodeScalingTest(unittest.TestCase):
    """fit_node_scaling: the replicated fraction and the time exponents from a one-node status
    set and a two-node runner status (synthetic figures with a known answer)."""

    def test_fractions_and_exponents(self):
        import fit_node_scaling
        gib = 2**30
        timing = [{"Index": i, "Iterations": 20, "SolveSeconds": 30.0, "TotalSeconds": 40.0} for i in range(10)]
        report = ("Elapsed Time Report (s)           Min.        Max.        Avg.\n"
                  "==============================================================\n"
                  "Initialization                   1.0       1.0       {init}\n"
                  "  Archive Reduction              1.0       1.0       {red}\n"
                  "--------------------------------------------------------------\n"
                  "Total                           {tot}     {tot}     {tot}\n")
        one = {"PBSJobID": "1.h", "Stages": [
            {"Name": "c-p4-worker-block1", "State": "complete", "WallSeconds": 500.0, "NodePeakUsedBytesSampled": 200 * gib,
             "Parsed": {"SourceTiming": timing, "PalacePeakMemory": {"Total": "150G", "Max": "150G"}}},
            {"Name": "c-p4-reducer", "State": "complete", "WallSeconds": 100.0, "NodePeakUsedBytesSampled": 300 * gib,
             "Parsed": {"ElapsedTimeReport": report.format(init=40.0, red=60.0, tot=100.0), "PalacePeakMemory": {"Total": "250G", "Max": "250G"}}}]}
        two_timing = [{**t, "TotalSeconds": 25.0, "SolveSeconds": 18.0} for t in timing]
        two = {"PBSJobID": "2.h", "Nodes": ["h1", "h2"], "Stages": [
            {"Name": "c-p4-worker-block1", "State": "complete", "WallSeconds": 330.0,
             "NodePeakUsedBytesSampled": 120 * gib, "NodePeakUsedBytesSampledPerNode": {"h1": 120 * gib, "h2": 110 * gib},
             "Parsed": {"SourceTiming": two_timing, "PalacePeakMemory": {"Total": "160G", "Max": "80G"}}},
            {"Name": "c-p4-reducer", "State": "complete", "WallSeconds": 70.0,
             "NodePeakUsedBytesSampled": 165 * gib, "NodePeakUsedBytesSampledPerNode": {"h1": 160 * gib, "h2": 165 * gib},
             "Parsed": {"ElapsedTimeReport": report.format(init=35.0, red=35.0, tot=70.0), "PalacePeakMemory": {"Total": "260G", "Max": "130G"}}}]}
        block = fit_node_scaling.fit([one], two)
        # Worker: 120 of 200 GiB per node -> r = 2 x 0.6 - 1 = 0.2; reducer 165 of 300 -> 0.1.
        self.assertAlmostEqual(block["Memory"]["Worker"]["ReplicatedFraction"], 0.2)
        self.assertAlmostEqual(block["Memory"]["Reducer"]["ReplicatedFraction"], 0.1)
        self.assertEqual(block["Memory"]["LocalEdge"]["ReplicatedFraction"], 0.2)   # carried from Worker without an ordinary stage
        self.assertIn("carried", block["Memory"]["LocalEdge"]["Rule"])
        # Per-source 40 -> 25 s: speedup 1.6; non-source 100 -> 80 s: 1.25; reducer setup 40 -> 35, reduction 60 -> 35.
        self.assertAlmostEqual(block["Time"]["WorkerPerSource"]["Speedup"], 1.6)
        self.assertAlmostEqual(block["Time"]["WorkerPerSource"]["Exponent"], 0.678, places=3)
        self.assertAlmostEqual(block["Time"]["WorkerNonSource"]["Speedup"], 1.25)
        self.assertAlmostEqual(block["Time"]["ReducerSetup"]["Speedup"], 40.0 / 35.0, places=3)
        self.assertAlmostEqual(block["Time"]["ReducerReduction"]["Speedup"], 60.0 / 35.0, places=3)
        self.assertEqual(block["Time"]["LocalEdge"], {"Speedup": 1.0, "Exponent": 0.0, "Rule": "not measured: no speedup assumed"})
        self.assertEqual(block["MeasuredNodes"], 2)
        # Better than even division is planned as even (r clamped at 0); the raw value is kept.
        clamped = fit_node_scaling.replicated_fraction(200 * gib, 90 * gib)
        self.assertEqual(clamped["ReplicatedFraction"], 0.0)
        self.assertLess(clamped["Raw"], 0.0)
        # The block drives estimate_stages: per-node figures and scaled times follow the rules.
        model = {**estimate_stages.load_cost_model(estimate_stages.DEVICE_COST_MODEL), "NodeScaling": block}
        self.assertAlmostEqual(estimate_stages.per_node_used_gib(model, "Reducer", 1000.0, 4), 1000.0 * (0.1 + 0.9 / 4))
        self.assertAlmostEqual(estimate_stages.scaled_seconds(model, "WorkerPerSource", 100.0, 2), 100.0 / 1.6, places=1)   # the exponent is recorded to 4 decimals


class LoopEndUnderMeasuredModelTest(unittest.TestCase):
    """The acceptance example of decision 457: the loop-end fab coupon 1b26671c9080 under the committed
    Version-3 model with the measured NodeScaling (the cluster dry run of record
    `stage2-20261004/qualify-multinode/evidence/loopend-dryrun/v3-nodescaling-speed10/` planned the same)."""

    @classmethod
    def setUpClass(cls):
        cls.model = estimate_stages.load_cost_model()
        cls.profile = json.loads((HERE / "qualify" / "cluster-profile.json").read_text())
        cls.layout = qualify_library.stage_layout("le", [4], [3, 5], LOOP_END_SOURCES, 8)
        stages = [(item["EstimateKey"], item["Order"], item["Sources"]) for item in cls.layout if item["Kind"] == "response"]
        local = next(item for item in cls.layout if item["Kind"] == "local-edge")
        cls.estimate = estimate_stages.estimate(LOOP_END_COUNTS, stages, local_edge=(local["EstimateKey"], local["Order"], local["Sources"]),
                                                model=cls.model, profile=cls.profile)

    def test_loop_end_plans_the_main_stage_on_four_r8g_nodes_and_the_fixed_group_on_two(self):
        self.assertIn("NodeScaling", self.model)
        required = self.estimate["Nodes"]["Required"]
        self.assertEqual(required, {"p4-846": 4, "p5-control-8": 1, "p3-control-8": 1, "local-edge-p4-8": 2})
        main = self.estimate["Stages"]["p4-846"]
        self.assertEqual(main["NodePlan"]["Instance"]["Type"], "r8g.48xlarge")
        self.assertLessEqual(main["NodePlan"]["Instance"]["PerNodeUsedGiBEstimate"], 0.9 * 0.6 * 1536)
        self.assertGreater(max(main["NodePlan"]["PerNodeUsedGiB"]["3"].values()), 0.9 * 0.6 * 1536)   # 3 nodes do not hold the reducer
        self.assertAlmostEqual(main["NodeUsedGiBEstimateReducer"], 1769.0, delta=5.0)   # the refit reducer line x 1.164 (was 2,988 GiB)
        nodes, assignment = qualify_library.node_assignment(self.estimate["Nodes"], self.layout)
        self.assertEqual(nodes, {"Main": 4, "Fixed": 2})
        scaled = estimate_stages.scale_to_nodes(self.estimate, assignment, self.model, self.profile)
        policy = job_split.normalize_policy("speed", max_jobs=10, walltime_seconds=self.profile["WalltimeSeconds"], user_job_cap=40)
        split = job_split.plan_split(indices=list(range(1, LOOP_END_SOURCES + 1)), layout=self.layout, estimate=scaled, policy=policy,
                                     model=self.model, profile=self.profile, nodes=nodes)
        self.assertTrue(split["Fits"])
        # Decision 468 (4): the 4-node stages are credited the measured 2-node speedup only (the
        # 4-node runs of the loop end are the measurement that extends the model).
        self.assertEqual(scaled["Stages"]["p4-846"]["SpeedupNodes"], 2)
        self.assertEqual((split["N"], split["ControlsJob"]), (10, "separate"))
        self.assertEqual([len(block) for block in split["Blocks"]], [0] + [94] * 9)
        self.assertEqual([job["Nodes"] for job in split["Jobs"]], [2] + [4] * 10)
        worst = split["WorstPCGFactor"]
        self.assertAlmostEqual(split["CriticalPathEstimateSeconds"][worst] / 60, 192.0, delta=1.0)
        self.assertAlmostEqual(split["NodeSecondsEstimate"][worst] / 3600, 78.9, delta=0.5)
        self.assertAlmostEqual(split["NodeSecondsEstimate"]["1.0"] / 3600, 57.0, delta=1.5)
        # The profile cap binds: at MaximumNodesPerJob 3 the main stage fails closed.
        tight = estimate_stages.estimate(LOOP_END_COUNTS, [("p4-846", 4, LOOP_END_SOURCES)], model=self.model,
                                         profile={**self.profile, "MaximumNodesPerJob": 3})
        self.assertFalse(tight["Nodes"]["Fits"])
