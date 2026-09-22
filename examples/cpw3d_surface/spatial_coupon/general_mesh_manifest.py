#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Fail-closed manifest and evidence gate for coupon mesh generality."""
import csv
import hashlib
import json
import math
import re
from pathlib import Path
import tomllib

from audit_edge_metric_mesh import (ANISOTROPY_GATE_APPLIED, ANISOTROPY_GATE_NOT_APPLICABLE,
                                    BAND_GATE_CUTOFFS, LAYER_ADJACENT_BAND_RULE)
from edge_volume_metric import EDGE_LAYER_QUALITY_RULE
from canonical_mesh_build import same_canonical_build, validate_build_record
from mesh_stage_contract import (GMSH_ONLY_PIPELINE, LEGACY_MMG_PIPELINE, PIPELINE_BUILD_STAGE,
                                 PLACEMENT_STAGE_ORDER, STAGE_TOOLS, TRACE_BASIS_RATIO_OPTION,
                                 canonical_stage_order, pipeline_of, scope_classes_of_case_inputs,
                                 stage_order, unsupported_scope_classes, validate_stage_dag)
from semantic_mesh_contract import (REQUIRED_ROLES, load_semantic_contract,
                                    validate_feature_topology)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
                                     separators=(",", ":")).encode()).hexdigest()


def _artifact_path(base, item):
    path = Path(item["Path"])
    return path if path.is_absolute() else (base / path).resolve()


def _check_artifact(base, item, description):
    if not isinstance(item, dict) or not item.get("Path") or not item.get("SHA256"):
        raise ValueError(f"{description} binding is incomplete")
    path = _artifact_path(base, item)
    if not path.is_file() or sha256(path) != item["SHA256"]:
        raise ValueError(f"{description} artifact hash mismatch")
    return path


def parent_labeled_mesh_sha256(evidence):
    """The digest the exact source-seed covariance judges: the reference variant's
    published mesh before the label-only radial MA shell relabel (decision 61a;
    ParentLabeledMeshSHA256 of its transform receipt) - the canonical mesh itself for an
    identity placement - or the mesh digest of evidence without that field."""
    return evidence.get("ParentLabeledMeshSHA256") or evidence["Mesh"]["SHA256"]


def _validate_mesh(path, contract, *, mesh=None):
    """The audited mesh reads and analyzes against the contract; `mesh` is the already
    loaded audit mesh (general_mesh_audit_producer.read_audit_mesh) when the caller reads
    it once for every check of the variant (decision 62 step 3, proposal 2): the memoized
    view and analysis are then shared with the recomputation."""
    from general_mesh_audit_producer import analyze as memoized_analyze, read_audit_mesh, simplicial_view as memoized_view
    try:
        memoized_analyze(memoized_view(read_audit_mesh(path) if mesh is None else mesh), contract,
                         require_material_names=True)
    except SystemExit as error:
        raise ValueError("audited mesh is not a readable Gmsh mesh") from error


def _transform_points(points, transform):
    result = []
    for x, y, z in points:
        value = (x, y, z, 1.0)
        result.append([sum(transform[4 * row + column] * value[column]
                           for column in range(4)) for row in range(3)])
    return result


def _same_points(expected, actual, tolerance):
    if not expected or not actual or len(expected) != len(actual):
        return False
    unused = [tuple(float(x) for x in point) for point in actual]
    for point in expected:
        point = tuple(float(x) for x in point)
        match = next((i for i, other in enumerate(unused)
                      if len(point) == len(other) and math.dist(point, other) <= tolerance), None)
        if match is None:
            return False
        unused.pop(match)
    return True


def _finite_number(value, *, nonnegative=False, positive=False):
    good = isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
    return good and (not nonnegative or value >= 0) and (not positive or value > 0)


def _is_rigid_transform(transform, tolerance=1e-12):
    if transform[12:] != [0, 0, 0, 1] and transform[12:] != [0.0, 0.0, 0.0, 1.0]:
        return False
    rotation = [[transform[4 * row + column] for column in range(3)]
                for row in range(3)]
    gram = [[sum(rotation[k][i] * rotation[k][j] for k in range(3))
             for j in range(3)] for i in range(3)]
    determinant = (rotation[0][0] * (rotation[1][1] * rotation[2][2] -
                                     rotation[1][2] * rotation[2][1]) -
                   rotation[0][1] * (rotation[1][0] * rotation[2][2] -
                                     rotation[1][2] * rotation[2][0]) +
                   rotation[0][2] * (rotation[1][0] * rotation[2][1] -
                                     rotation[1][1] * rotation[2][0]))
    return (all(abs(gram[i][j] - float(i == j)) <= tolerance
                for i in range(3) for j in range(3)) and
            abs(determinant - 1.0) <= tolerance)


PRODUCER_DEFAULT_ETCH_FOOTPRINT = "producer-default"
# The production recipe: the options every production build executes, recorded in
# the production manifest and bound to the recorded stage commands of every
# production case (the canonical cache key does not encode recipe options; the
# census / recipe records are bound by the stage contract).  Under the Gmsh-only
# pipeline (supervisor decision 38) the recipe is the build command's tube, corner,
# band and far-field options (BuildCommandOptions); under the legacy MMG pipeline
# (decisions 34B/35) the seed, metric and adaptation options.  A labeled calibration
# manifest never carries it: its cases declare their own options against the
# recorded pre-34B production values.
PRODUCTION_RECIPE_KEY = "ProductionRecipe"
PIPELINE_PRODUCTION_RECIPE_STAGE_OPTIONS = {
    GMSH_ONLY_PIPELINE: {"gmsh-build": "BuildCommandOptions"},
    LEGACY_MMG_PIPELINE: {"seed-generation": "SeedCommandOptions",
                          "metric-preparation": "MetricCommandOptions",
                          "native-adaptation-mmg": "AdaptationCommandOptions"},
}
PRODUCTION_RECIPE_STAGE_OPTIONS = PIPELINE_PRODUCTION_RECIPE_STAGE_OPTIONS[LEGACY_MMG_PIPELINE]
MANIFEST_PIPELINE_KEY = "Pipeline"


def manifest_pipeline(manifest):
    """The canonical pipeline of a manifest: the one whose exact stage set its
    StageToolSHA256 freezes; an explicit `Pipeline` key must name the same one."""
    stage_tools = manifest.get("StageToolSHA256")
    if not isinstance(stage_tools, dict):
        raise ValueError("Manifest must freeze every stage tool digest")
    pipeline = pipeline_of(stage_tools)
    declared = manifest.get(MANIFEST_PIPELINE_KEY)
    if declared is not None and declared != pipeline:
        raise ValueError(f"Manifest declares pipeline {declared} but freezes the {pipeline} stages")
    return pipeline


def option_values(command, option):
    """Every value the recorded argv passes to `option` (separate tokens), as floats."""
    values = []
    for index, token in enumerate(command):
        if token == option:
            if index + 1 >= len(command):
                raise ValueError(f"recorded command ends with option {option}")
            values.append(float(command[index + 1]))
    return values


def validate_production_recipe(manifest):
    """ProductionRecipe (if present) names, per stage of PRODUCTION_RECIPE_STAGE_OPTIONS,
    a non-empty dict of option -> finite number, in a manifest without a Calibration
    block.  Returns the block or None."""
    recipe = manifest.get(PRODUCTION_RECIPE_KEY)
    if recipe is None:
        return None
    if "Calibration" in manifest:
        raise ValueError("A calibration manifest cannot carry a production recipe")
    keys = PIPELINE_PRODUCTION_RECIPE_STAGE_OPTIONS[manifest_pipeline(manifest)].values()
    if (not isinstance(recipe, dict) or
            any(not isinstance(recipe.get(key), dict) or not recipe[key] or
                any(not isinstance(option, str) or not option.startswith("--") or
                    not _finite_number(value)
                    for option, value in recipe[key].items())
                for key in keys)):
        raise ValueError("Production recipe must bind finite option values for the pipeline's "
                         "recipe stages: " + ", ".join(keys))
    if manifest_pipeline(manifest) == GMSH_ONLY_PIPELINE:
        validate_build_cost_estimate_model(recipe)
        validate_physics_run(recipe)
    return recipe


BUILD_COST_ESTIMATE_KEY = "BuildCostEstimate"
PHYSICS_RUN_KEY = "PhysicsRun"


JOB_POLICY_KEY = "JobPolicy"
JOB_POLICY_MODES = ("speed", "frugal", "fixed")
REDUCER_BLOCK_SIZE_KEY = "ReducerBlockSize"
FROZEN_EXECUTABLE_KEY = "FrozenExecutable"
SHA256_HEX = re.compile(r"^[0-9a-f]{64}$")


def validate_physics_run(recipe):
    """The Gmsh-only production recipe binds the Palace run parameters of every coupon
    (supervisor decision 52): PhysicsRun.Order (integer >= 1), PhysicsRun.LinearTol (0 <
    float < 1) with a Provenance and a Binding text; qualify derives every run config
    from the case's sources and these values (qualify/case_inputs.py)."""
    block = recipe.get(PHYSICS_RUN_KEY)
    if (not isinstance(block, dict) or not isinstance(block.get("Order"), int) or isinstance(block.get("Order"), bool) or
            block["Order"] < 1 or not isinstance(block.get("LinearTol"), float) or not 0.0 < block["LinearTol"] < 1.0 or
            any(not isinstance(block.get(key), str) or not block[key] for key in ("Rule", "Provenance", "Binding"))):
        raise ValueError("Gmsh-only production recipe must carry the PhysicsRun block "
                         "(Order int >= 1, LinearTol in (0, 1), Rule, Provenance, Binding)")
    policy = block.get(JOB_POLICY_KEY)
    if policy is not None:
        # The recorded default of qualify's per-coupon source split (decision 61b): Mode
        # speed / frugal / fixed, FixedJobs (int >= 1) exactly for fixed, a Rule text.
        mode = policy.get("Mode") if isinstance(policy, dict) else None
        fixed = policy.get("FixedJobs") if isinstance(policy, dict) else None
        if (mode not in JOB_POLICY_MODES or
                (mode == "fixed") != (isinstance(fixed, int) and not isinstance(fixed, bool) and fixed >= 1) or
                not isinstance(policy.get("Rule"), str) or not policy["Rule"]):
            raise ValueError("PhysicsRun.JobPolicy must carry Mode speed | frugal | fixed, FixedJobs (int >= 1) "
                             "exactly for fixed, and a Rule")
    reducer_block_size = block.get(REDUCER_BLOCK_SIZE_KEY)
    if reducer_block_size is not None:
        # The recorded default of the reducer's PALACE_RESPONSE_BLOCK_SIZE (decision 62(1)):
        # Value (int >= 1) with a Rule text (the memory rationale).
        value = reducer_block_size.get("Value") if isinstance(reducer_block_size, dict) else None
        if (not isinstance(value, int) or isinstance(value, bool) or value < 1 or
                not isinstance(reducer_block_size.get("Rule"), str) or not reducer_block_size["Rule"]):
            raise ValueError("PhysicsRun.ReducerBlockSize must carry Value (int >= 1) and a Rule")
    frozen = block.get(FROZEN_EXECUTABLE_KEY)
    if frozen is not None:
        # The recorded default frozen Palace executable of every qualify stage (decision 63):
        # SHA256 (64 lowercase hex), the PreviousSHA256 it replaced (or None), Rule, Provenance.
        digest = frozen.get("SHA256") if isinstance(frozen, dict) else None
        previous = frozen.get("PreviousSHA256") if isinstance(frozen, dict) else None
        if (not isinstance(digest, str) or not SHA256_HEX.match(digest) or
                (previous is not None and (not isinstance(previous, str) or not SHA256_HEX.match(previous))) or
                any(not isinstance(frozen.get(key), str) or not frozen[key] for key in ("Rule", "Provenance"))):
            raise ValueError("PhysicsRun.FrozenExecutable must carry SHA256 (64 hex), PreviousSHA256 (64 hex or null), "
                             "Rule and Provenance")
    return block
SCOPE_KEY = "Scope"
UNSUPPORTED_CLASS_KEY = "UnsupportedClass"


class UnsupportedClassError(ValueError):
    """A case exhibits a class the Gmsh-only recipe guards (supervisor decision 48):
    recorded as "unsupported class <id>" distinctly from any other preflight failure."""

    def __init__(self, guard):
        super().__init__(f"unsupported class {guard}")
        self.guard = guard


def preflight_recipe_scope(manifest, paths, case):
    """The recipe scope classification of a Gmsh-only case from its frozen inputs
    (signature, boundary, process, footprint declaration, trace basis): the exhibited
    classes and the guarded ones among them; None for a legacy manifest."""
    if manifest.get("Pipeline") != GMSH_ONLY_PIPELINE:
        return None
    process = tomllib.loads(paths["Process"].read_text())
    classes = scope_classes_of_case_inputs(paths["Signature"], paths["Boundary"], process,
                                          device_footprint="RetainedEtch" in paths,
                                          trace_basis="BasisContract" in paths)
    return {"ExhibitedClasses": classes, "UnsupportedClasses": unsupported_scope_classes(classes)}


def validate_build_cost_estimate_model(recipe):
    """The Gmsh-only production recipe carries the pre-build element estimate model
    (estimate_build_cost.py): a positive finite TetrahedraPerCubicSize with its rule,
    calibration record and binding text.  The estimate is gated by the manifest's
    MaximumElements before any build (preflight and run_gmsh_only_case.py)."""
    model = recipe.get(BUILD_COST_ESTIMATE_KEY)
    if (not isinstance(model, dict) or not _finite_number(model.get("TetrahedraPerCubicSize"), positive=True) or
            any(not isinstance(model.get(key), str) or not model[key] for key in ("Rule", "Binding")) or
            not isinstance(model.get("Calibration"), dict) or not model["Calibration"]):
        raise ValueError("Gmsh-only production recipe must carry the BuildCostEstimate model "
                         "(TetrahedraPerCubicSize, Rule, Calibration, Binding)")
    return model


def preflight_build_cost(manifest, manifest_path, case):
    """The pre-build element estimate of a case against MaximumElements (fail closed): a
    production case with its recipe's options and model; a case of a labeled Gmsh-only
    calibration manifest with its own options and the production manifest's model
    (estimate_build_cost.build_options_and_model); None for a legacy manifest."""
    if manifest.get("Pipeline") != GMSH_ONLY_PIPELINE:
        return None
    recipe = manifest.get(PRODUCTION_RECIPE_KEY)
    if "Calibration" not in manifest and (recipe is None or BUILD_COST_ESTIMATE_KEY not in recipe):
        return None
    from estimate_build_cost import gate as estimate_gate
    return estimate_gate(manifest, manifest_path, case)


def validate_production_recipe_commands(manifest, case, bounded_stages):
    """A production case (no Calibration block) of a manifest carrying ProductionRecipe
    executed every recorded option exactly once at its recorded value in the seed,
    metric and adaptation commands.  The dimensionless trace-basis size ratio
    (`--trace-basis-size-ratio`, supervisor decision 42) is passed with the bound
    trace basis only: a case freezing the trace basis executes it exactly once at
    the recipe value, a case without one never executes it.  Raises ValueError
    otherwise."""
    recipe = manifest.get(PRODUCTION_RECIPE_KEY)
    if recipe is None or case.get("Calibration") is not None:
        return
    files = case.get("Source", {}).get("Files", {})
    trace_basis_bound = all(role in files for role in TRACE_BASIS_ROLES)
    for stage, key in PIPELINE_PRODUCTION_RECIPE_STAGE_OPTIONS[manifest_pipeline(manifest)].items():
        command = bounded_stages[stage]["Command"]
        for option, value in recipe[key].items():
            executed = option_values(command, option)
            expected = [float(value)]
            if option == TRACE_BASIS_RATIO_OPTION and not trace_basis_bound:
                expected = []
            if executed != expected:
                raise ValueError(f"{stage} command does not execute the production recipe "
                                 f"option {option}={value} exactly once (executed {executed})"
                                 if expected else
                                 f"{stage} command executes the production recipe option "
                                 f"{option} without a bound trace basis (executed {executed})")

# A Gmsh-only calibration case (labeled calibration manifest of the Gmsh-only
# pipeline, supervisor decision 41) declares the build options that differ from
# production (Calibration.BuildCommandOptions, any subset) against the production
# values it was declared at (Calibration.ProductionValues); the recorded gmsh-build
# command is bound to them by verify_canonical_case_entries.validate_calibration_commands.
CALIBRATION_BUILD_OPTIONS_KEY = "BuildCommandOptions"
# The production values a calibration case's options are declared against: under the
# legacy MMG pipeline those of the production recipe before decision 34B (seed
# --lc-tangent 0.1, metric --far-growth 1.0, no edge layer, adapter --hmin
# NormalSize, no corner grading); under the Gmsh-only pipeline the production
# BuildCommandOptions at the time of the decision-41 study (before decision 42:
# --trace-basis-size-ratio 1.0).  Both are historical baselines: an adopted case
# (EL4c under 34B, V-a under 42) stays declared against them and is labeled
# AdoptedAsProductionRecipe.
PIPELINE_CALIBRATION_PRODUCTION_VALUES_KEY = {LEGACY_MMG_PIPELINE: "ProductionValuesBefore34B",
                                              GMSH_ONLY_PIPELINE: "ProductionValues"}
# A label-only calibration case (supervisor decision 56): its mesh is a relabel of a
# built base case's identity mesh (relabel_radial_ma_shells.py), declared under
# Calibration.Relabel instead of BuildCommandOptions; the mesher never builds it.
RELABEL_KEY = "Relabel"
RELABEL_KINDS = ("radial-ma-shells",)


def validate_calibration_relabel(case, calibration):
    """Calibration.Relabel (if present) binds Kind (RELABEL_KINDS), BaseCase equal to
    Calibration.BaseCase, the ParentMeshSHA256 of the relabeled identity mesh, the tube
    RingRadii (strictly increasing positive) and the Tool name; a relabel case declares
    no BuildCommandOptions.  Returns the block or None."""
    relabel = calibration.get(RELABEL_KEY)
    if relabel is None:
        return None
    radii = relabel.get("RingRadii") if isinstance(relabel, dict) else None
    digest = relabel.get("ParentMeshSHA256") if isinstance(relabel, dict) else None
    if (not isinstance(relabel, dict) or relabel.get("Kind") not in RELABEL_KINDS or
            relabel.get("BaseCase") != calibration.get("BaseCase") or not isinstance(relabel.get("BaseCase"), str) or
            not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest) or
            not isinstance(radii, list) or not radii or any(not _finite_number(r, positive=True) for r in radii) or
            any(b <= a for a, b in zip(radii, radii[1:])) or
            not isinstance(relabel.get("Tool"), str) or not relabel["Tool"] or
            CALIBRATION_BUILD_OPTIONS_KEY in calibration):
        raise ValueError(f"{case.get('Id')} declares a Relabel block that is not a labeled label-only relabel of its "
                         f"base case (Kind, BaseCase, ParentMeshSHA256, increasing RingRadii, Tool; no BuildCommandOptions)")
    return relabel


def validate_calibration_case_options(manifest, case):
    """Under a Gmsh-only calibration manifest every case declares a Calibration block
    with a Label, a non-empty BuildCommandOptions dict of option -> finite number and a
    ProductionValues dict holding every declared option at a different finite value.
    Legacy-pipeline calibration cases are declared per stage (verified by the
    per-case verifier) and are not judged here.  Raises ValueError otherwise."""
    if "Calibration" not in manifest or manifest_pipeline(manifest) != GMSH_ONLY_PIPELINE:
        return
    calibration = case.get("Calibration")
    if isinstance(calibration, dict) and validate_calibration_relabel(case, calibration) is not None:
        production = calibration.get(PIPELINE_CALIBRATION_PRODUCTION_VALUES_KEY[GMSH_ONLY_PIPELINE])
        if (not isinstance(calibration.get("Label"), str) or not calibration["Label"] or not isinstance(production, dict) or
                not production or any(not isinstance(option, str) or not option.startswith("--") or not _finite_number(value)
                                      for option, value in production.items())):
            raise ValueError(f"{case.get('Id')} must declare its Gmsh-only calibration label and the production values "
                             f"its relabeled base mesh was built at")
        return
    options = calibration.get(CALIBRATION_BUILD_OPTIONS_KEY) if isinstance(calibration, dict) else None
    production = (calibration.get(PIPELINE_CALIBRATION_PRODUCTION_VALUES_KEY[GMSH_ONLY_PIPELINE])
                  if isinstance(calibration, dict) else None)
    if (not isinstance(calibration.get("Label") if isinstance(calibration, dict) else None, str) or
            not calibration["Label"] or
            not isinstance(options, dict) or not options or not isinstance(production, dict) or
            any(not isinstance(option, str) or not option.startswith("--") or not _finite_number(value)
                for option, value in list(options.items()) + list(production.items())) or
            any(option not in production or float(production[option]) == float(value)
                for option, value in options.items())):
        raise ValueError(f"{case.get('Id')} must declare its Gmsh-only calibration label and "
                         f"build options against differing production values")


# The trace basis is bound (seed cut-surface sizing and metric record) exactly when
# a case freezes all four roles; a partial set fails preflight.
TRACE_BASIS_ROLES = ("BasisContract", "TraceVertices", "TraceTriangles", "ProcessLibrary")


# The seeded edge layer of a calibration case is declared under its Calibration
# block (Calibration.EdgeLayer); the layer-local quality rule (decision 32) is the
# manifest gate Gates.EdgeLayerQualityRule, allowed only in a calibration manifest
# that labels it as a gate deviation.
EDGE_LAYER_CASE_KEY = "EdgeLayer"
EDGE_LAYER_QUALITY_RULE_GATE = "EdgeLayerQualityRule"
EDGE_LAYER_QUALITY_RULE_OPTION = "--edge-layer-maximum-aspect"


def validate_edge_layer_quality_rule_gate(manifest):
    """Gates.EdgeLayerQualityRule (if present) is a labeled calibration-only rule:
    a finite MaximumEdgeAspect > 1, a ScaledJacobianRoundoffFloor in (0,
    MinimumScaledJacobian), and the manifest's Calibration.GateDeviations names it
    with ProductionUse FORBIDDEN.  A manifest without a Calibration block cannot
    carry it.  Returns the rule or None."""
    gates = manifest.get("Gates", {})
    rule = gates.get(EDGE_LAYER_QUALITY_RULE_GATE)
    if rule is None:
        return None
    calibration = manifest.get("Calibration")
    if not isinstance(calibration, dict):
        raise ValueError("A production manifest cannot carry an edge-layer quality rule")
    deviation = calibration.get("GateDeviations", {}).get(EDGE_LAYER_QUALITY_RULE_GATE)
    if (not isinstance(rule, dict) or
            not _finite_number(rule.get("MaximumEdgeAspect"), positive=True) or
            rule["MaximumEdgeAspect"] <= 1.0 or
            not _finite_number(rule.get("ScaledJacobianRoundoffFloor"), positive=True) or
            not rule["ScaledJacobianRoundoffFloor"] < gates.get("MinimumScaledJacobian", 0.0) or
            not isinstance(deviation, dict) or
            "FORBIDDEN" not in str(deviation.get("ProductionUse", "")) or
            deviation.get("Calibration") != rule["MaximumEdgeAspect"] or
            deviation.get("Production") is not None):
        raise ValueError("Edge-layer quality rule gate is invalid or not labeled as a "
                         "calibration-only deviation")
    return rule


def case_gates(manifest, case):
    """The gates that judge one case's evidence.  The calibration-only edge-layer
    quality rule (Gates.EdgeLayerQualityRule) judges only a case declaring
    Calibration.EdgeLayerQualityRule, whose seed and label restorer executed its
    bound (verify_canonical_case_entries.validate_edge_layer_quality_rule_binding);
    every other case - a layer case without the declaration included (the 4 nm
    layer case, judged by MinimumScaledJacobian on its whole mesh as recorded) - is
    judged by the manifest gates without the rule.  A case declaring
    Calibration.MaximumElements (a labeled calibration-only element cap,
    validate_case_element_cap) is judged by that cap.  The manifest's Gates stay
    the canonical cache key of the build."""
    gates = dict(manifest["Gates"])
    calibration = case.get("Calibration")
    calibration = calibration if isinstance(calibration, dict) else {}
    if EDGE_LAYER_QUALITY_RULE_GATE in gates and calibration.get(EDGE_LAYER_QUALITY_RULE_GATE) is None:
        del gates[EDGE_LAYER_QUALITY_RULE_GATE]
    if calibration.get(ELEMENT_CAP_GATE) is not None:
        gates[ELEMENT_CAP_GATE] = validate_case_element_cap(manifest, case)
    if calibration.get(JACOBIAN_CONDITION_GATE) is not None:
        gates[JACOBIAN_CONDITION_GATE] = validate_case_jacobian_condition(manifest, case)
    return gates


# A calibration case may carry its own Jacobian condition bound
# (Calibration.MaximumJacobianCondition, supervisor decision 53: 1200 for the 0.125 nm
# tube ring set, whose innermost prisms have twice the production ring / layer aspect -
# 2 x 586.30 = 1172.6 by construction, not a quality loss), labeled in the manifest's
# Calibration.GateDeviations.MaximumJacobianCondition naming the case; every other case
# and every production case keeps the manifest gate.
JACOBIAN_CONDITION_GATE = "MaximumJacobianCondition"


def validate_case_jacobian_condition(manifest, case):
    """Calibration.MaximumJacobianCondition of a case (if present) is a labeled
    calibration-only deviation: a finite number above the manifest gate, equal to the
    deviation's Calibration value, the case named in the deviation's Cases, Production
    equal to the manifest gate and ProductionUse FORBIDDEN.  Returns the bound or None."""
    calibration = case.get("Calibration")
    bound = calibration.get(JACOBIAN_CONDITION_GATE) if isinstance(calibration, dict) else None
    if bound is None:
        return None
    manifest_bound = manifest.get("Gates", {}).get(JACOBIAN_CONDITION_GATE)
    deviation = manifest.get("Calibration", {}).get("GateDeviations", {}).get(JACOBIAN_CONDITION_GATE)
    if (not _finite_number(bound, positive=True) or not _finite_number(manifest_bound, positive=True) or
            bound <= manifest_bound or not isinstance(deviation, dict) or
            deviation.get("Production") != manifest_bound or deviation.get("Calibration") != bound or
            not isinstance(deviation.get("Cases"), list) or case.get("Id") not in deviation["Cases"] or
            "FORBIDDEN" not in str(deviation.get("ProductionUse", ""))):
        raise ValueError(f"{case.get('Id')} declares a Jacobian condition bound that is not labeled as a "
                         f"calibration-only deviation")
    return float(bound)


# A calibration case may carry its own element cap (Calibration.MaximumElements,
# supervisor decision 33: 5,000,000 for the 1 nm aspect-4 layer), labeled in the
# manifest's Calibration.GateDeviations.MaximumElements naming the case; every
# other case and every production case keeps the manifest gate.
ELEMENT_CAP_GATE = "MaximumElements"


# A calibration case may run its physics at its own linear tolerance
# (Calibration.PhysicsRun.LinearTol, supervisor decision 56: 1e-8 for the radial MA
# shell relabel of the two-edge production mesh, so its summed MA reproduces the
# production acceptance run PBS 46023 - which ran the producer default 1e-8 - at the
# deterministic-solve level), labeled in the manifest's Calibration.GateDeviations.
# LinearTol naming the case; every other case runs at the recipe's PhysicsRun value.
LINEAR_TOL_GATE = "LinearTol"


def validate_case_linear_tol(manifest, case):
    """Calibration.PhysicsRun.LinearTol of a case (if present) is a labeled calibration-only
    deviation: a float in (0, 1) equal to the deviation's Calibration value, the case
    named in the deviation's Cases, a Production value in (0, 1) differing from it (bound
    to the recipe by qualify/case_inputs.physics_run_parameters) and ProductionUse
    FORBIDDEN.  Returns the tolerance or None."""
    calibration = case.get("Calibration")
    physics = calibration.get("PhysicsRun") if isinstance(calibration, dict) else None
    if physics is None:
        return None
    tolerance = physics.get(LINEAR_TOL_GATE) if isinstance(physics, dict) else None
    deviation = manifest.get("Calibration", {}).get("GateDeviations", {}).get(LINEAR_TOL_GATE)
    if (not isinstance(physics, dict) or set(physics) != {LINEAR_TOL_GATE} or not isinstance(tolerance, float) or
            not 0.0 < tolerance < 1.0 or not isinstance(deviation, dict) or
            not isinstance(deviation.get("Production"), float) or not 0.0 < deviation["Production"] < 1.0 or
            deviation["Production"] == tolerance or deviation.get("Calibration") != tolerance or
            not isinstance(deviation.get("Cases"), list) or case.get("Id") not in deviation["Cases"] or
            "FORBIDDEN" not in str(deviation.get("ProductionUse", ""))):
        raise ValueError(f"{case.get('Id')} declares a physics-run linear tolerance that is not labeled as a "
                         f"calibration-only deviation")
    return tolerance


def validate_case_element_cap(manifest, case):
    """Calibration.MaximumElements of a case (if present) is a labeled calibration-only
    deviation: an integer above the manifest gate, equal to the deviation's
    Calibration value, the case named in the deviation's Cases, Production equal to
    the manifest gate and ProductionUse FORBIDDEN.  Returns the cap or None."""
    calibration = case.get("Calibration")
    cap = calibration.get(ELEMENT_CAP_GATE) if isinstance(calibration, dict) else None
    if cap is None:
        return None
    manifest_cap = manifest.get("Gates", {}).get(ELEMENT_CAP_GATE)
    deviation = manifest.get("Calibration", {}).get("GateDeviations", {}).get(ELEMENT_CAP_GATE)
    if (isinstance(cap, bool) or not isinstance(cap, int) or not isinstance(manifest_cap, int) or
            cap <= manifest_cap or not isinstance(deviation, dict) or
            deviation.get("Production") != manifest_cap or deviation.get("Calibration") != cap or
            not isinstance(deviation.get("Cases"), list) or case.get("Id") not in deviation["Cases"] or
            "FORBIDDEN" not in str(deviation.get("ProductionUse", ""))):
        raise ValueError(f"{case.get('Id')} declares an element cap that is not labeled as a "
                         f"calibration-only deviation")
    return cap


def validate_manifest(manifest, manifest_path, *, check_available_files=True):
    if manifest.get("Version") != 2 or not isinstance(manifest.get("Cases"), list):
        raise ValueError("Unsupported generality-suite manifest")
    if not manifest["Cases"]:
        raise ValueError("Manifest must declare cases")
    identifiers = [case.get("Id") for case in manifest["Cases"]]
    if any(not value for value in identifiers) or len(set(identifiers)) != len(identifiers):
        raise ValueError("Manifest case identifiers must be nonempty and unique")
    gates = manifest.get("Gates", {})
    required_gates = ("CornerTolerance", "MaximumNormalFactor", "MinimumAchievedAspect",
                      "MaximumCornerAspect", "MinimumNoncornerAspect",
                      "MaximumProtectedMeasureError", "MinimumScaledJacobian",
                      "MaximumJacobianCondition", "MaximumSeconds", "MaximumRSSGiB",
                      "MaximumElements")
    if any(not _finite_number(gates.get(name), nonnegative=True) for name in required_gates):
        raise ValueError("Manifest has missing or invalid mesh gates")
    validate_edge_layer_quality_rule_gate(manifest)
    validate_production_recipe(manifest)
    repository = (manifest_path.parent / manifest["RepositoryRoot"]).resolve()
    tools = manifest.get("Tools")
    if not isinstance(tools, list) or not tools:
        raise ValueError("Manifest must freeze at least one evidence tool")
    tool_hashes = {}
    for tool in tools:
        name, digest = tool.get("Name"), tool.get("SHA256")
        if not name or name in tool_hashes or not digest:
            raise ValueError("Tool names and hashes must be nonempty and unique")
        tool_hashes[name] = digest
        if check_available_files:
            path = repository / tool["Path"]
            if not path.is_file() or sha256(path) != digest:
                raise ValueError(f"frozen tool hash mismatch: {name}")
    stage_tools = manifest.get("StageToolSHA256")
    pipeline = manifest_pipeline(manifest)
    if (any(not isinstance(stage_tools[stage], dict) or
            set(stage_tools[stage]) != STAGE_TOOLS[stage] or
            any(not isinstance(value, str) or len(value) != 64
                for value in stage_tools[stage].values())
            for stage in stage_order(pipeline))):
        raise ValueError("Manifest must freeze every stage tool digest")
    matrix = set()
    comparison_kinds = set()
    calibration_manifest = "Calibration" in manifest
    for case in manifest["Cases"]:
        # A seeded edge layer and its calibration-only gate rule exist only in a
        # labeled calibration manifest: a production case never declares either.
        if not calibration_manifest and any(
                key in case for key in ("Calibration", EDGE_LAYER_CASE_KEY)):
            raise ValueError(f"{case['Id']} declares a calibration or edge-layer block in a "
                             f"production manifest")
        if EDGE_LAYER_CASE_KEY in case:
            raise ValueError(f"{case['Id']} must declare its edge layer under Calibration")
        validate_case_element_cap(manifest, case)
        validate_case_jacobian_condition(manifest, case)
        validate_case_linear_tol(manifest, case)
        validate_calibration_case_options(manifest, case)
        source = case.get("Source", {})
        files = source.get("Files")
        if not isinstance(files, dict) or any(role not in files for role in REQUIRED_ROLES):
            raise ValueError(f"{case['Id']} lacks a required immutable role")
        if "MeshRecipe" not in files:
            raise ValueError(f"{case['Id']} lacks a frozen mesh recipe")
        # The etched footprint is a recorded choice, never an omission: either a
        # bound RetainedEtch file or the explicit producer default.
        declared_default = "EtchFootprint" in source
        if ("RetainedEtch" in files) == declared_default or (
                declared_default and source["EtchFootprint"] != PRODUCER_DEFAULT_ETCH_FOOTPRINT):
            raise ValueError(f"{case['Id']} must declare exactly one etch footprint: "
                             f"a RetainedEtch file or EtchFootprint "
                             f"\"{PRODUCER_DEFAULT_ETCH_FOOTPRINT}\"")
        declared_basis = {role for role in TRACE_BASIS_ROLES if role in files}
        if declared_basis and declared_basis != set(TRACE_BASIS_ROLES):
            raise ValueError(f"{case['Id']} must freeze all trace basis roles "
                             f"{list(TRACE_BASIS_ROLES)} or none")
        variants = case.get("Variants")
        if not isinstance(variants, list) or not variants:
            raise ValueError(f"{case['Id']} has no variants")
        variant_ids = []
        for variant in variants:
            if (not isinstance(variant, dict) or not variant.get("Id") or
                    not isinstance(variant.get("Transform"), list) or
                    len(variant["Transform"]) != 16 or
                    not all(_finite_number(x) for x in variant["Transform"]) or
                    not _is_rigid_transform(variant["Transform"])):
                raise ValueError(f"{case['Id']} has an invalid variant transform")
            variant_ids.append(variant["Id"])
            matrix.add((case["Id"], variant["Id"]))
        if len(set(variant_ids)) != len(variant_ids) or "identity" not in variant_ids:
            raise ValueError(f"{case['Id']} variants must be unique and include identity")
        by_id = {variant["Id"]: variant["Transform"] for variant in variants}
        identity = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
        if by_id["identity"] != identity:
            raise ValueError(f"{case['Id']} identity variant is not the identity transform")
        comparison = case.get("TransformComparison")
        physical_tolerances = (
            "MaximumRelativeVolumeError", "MaximumRelativeSurfaceMeasureError",
            "MaximumProtectedSupportHausdorff", "MaximumProtectedMeasureError",
            "MaximumQualityDistributionRelativeError", "MaximumAnisotropyRelativeError",
            "MaximumComplexityRatio")
        if (not isinstance(comparison, dict) or comparison.get("Reference") not in variant_ids or
                comparison.get("Transformed") not in variant_ids or
                comparison.get("Reference") == comparison.get("Transformed") or
                by_id.get(comparison.get("Reference")) ==
                by_id.get(comparison.get("Transformed")) or
                any(not _finite_number(comparison.get(name), nonnegative=True)
                    for name in physical_tolerances) or
                comparison.get("MaximumComplexityRatio", 0) < 1.0):
            raise ValueError(f"{case['Id']} has no frozen physical covariance comparison")
    comparisons = manifest.get("ScalingComparisons")
    if not isinstance(comparisons, list) or not comparisons:
        raise ValueError("Manifest must declare scaling comparisons")
    for comparison in comparisons:
        kind = comparison.get("Kind")
        refs = (tuple(comparison.get("Reference", [])), tuple(comparison.get("Compared", [])))
        if kind not in ("feature-scaling", "cad-subdivision-sensitivity") or any(ref not in matrix for ref in refs):
            raise ValueError("Invalid scaling comparison")
        if not _finite_number(comparison.get("MaximumNormalizedDOFRatio"), positive=True):
            raise ValueError("Scaling comparison has no positive bound")
        comparison_kinds.add(kind)
    # A labeled element-cap deviation names exactly the cases declaring the cap.
    cap_deviation = manifest.get("Calibration", {}).get("GateDeviations", {}).get(ELEMENT_CAP_GATE)
    if cap_deviation is not None:
        declaring = sorted(case["Id"] for case in manifest["Cases"]
                           if isinstance(case.get("Calibration"), dict) and
                           case["Calibration"].get(ELEMENT_CAP_GATE) is not None)
        if not isinstance(cap_deviation, dict) or sorted(cap_deviation.get("Cases") or []) != declaring:
            raise ValueError("Element cap deviation must name exactly the cases declaring the cap")
    condition_deviation = manifest.get("Calibration", {}).get("GateDeviations", {}).get(JACOBIAN_CONDITION_GATE)
    if condition_deviation is not None:
        declaring = sorted(case["Id"] for case in manifest["Cases"]
                           if isinstance(case.get("Calibration"), dict) and
                           case["Calibration"].get(JACOBIAN_CONDITION_GATE) is not None)
        if not isinstance(condition_deviation, dict) or sorted(condition_deviation.get("Cases") or []) != declaring:
            raise ValueError("Jacobian condition deviation must name exactly the cases declaring the bound")
    tolerance_deviation = manifest.get("Calibration", {}).get("GateDeviations", {}).get(LINEAR_TOL_GATE)
    if tolerance_deviation is not None:
        declaring = sorted(case["Id"] for case in manifest["Cases"]
                           if isinstance(case.get("Calibration"), dict) and
                           isinstance(case["Calibration"].get("PhysicsRun"), dict) and
                           case["Calibration"]["PhysicsRun"].get(LINEAR_TOL_GATE) is not None)
        if not isinstance(tolerance_deviation, dict) or sorted(tolerance_deviation.get("Cases") or []) != declaring:
            raise ValueError("Linear tolerance deviation must name exactly the cases declaring the tolerance")
    if comparison_kinds != {"feature-scaling", "cad-subdivision-sensitivity"}:
        raise ValueError("Both feature and CAD-subdivision scaling controls are required")
    return repository, tool_hashes, matrix


def _recompute_mesh_measurements(mesh_path, binding, source_paths, bounded, *, mesh=None):
    """Rerun the frozen producer implementation instead of trusting record JSON; `mesh`
    is the loaded audit mesh when the caller reads it once (the three producers then
    share it and its memoized kernels; decision 62 step 3, proposal 2)."""
    from general_mesh_audit_producer import (bounded_record, complexity_record,
                                             invariants_record, topology_record)
    base = {"Transform": binding["Transform"]}
    # Recover reference and ownership paths from the independently validated
    # embedded stage bindings.
    reference = build_volume_binding(bounded)["Path"]
    restoration_recipe, build_census = feature_record_bindings(bounded)
    publication = bounded["proper-rigid-publication"]["Artifacts"]
    ownership = publication["ownership-partition"]["Path"]
    quadrature = publication["ownership-quadrature-partition"]["Path"]
    topology = topology_record(dict(base), mesh_path, source_paths["SemanticContract"],
        source_paths["MeshRecipe"], source_paths["Process"], source_paths["Signature"],
        reference, ownership, quadrature,
        None if restoration_recipe is None else restoration_recipe["Path"],
        build_census_path=None if build_census is None else build_census["Path"], mesh=mesh)["Measurements"]
    complexity = complexity_record(dict(base), mesh_path, source_paths["SemanticContract"],
                                   source_paths["MeshRecipe"], mesh=mesh)["Measurements"]
    invariants = invariants_record(dict(base), mesh_path, mesh=mesh)["Measurements"]
    return {**topology, **complexity, **invariants}


def bounded_pipeline(bounded):
    """The pipeline of a bounded-stage dict (keyed by stage name)."""
    return pipeline_of(bounded)


def build_volume_binding(bounded):
    """The source-local volume build artifact binding: the Gmsh-only build mesh or
    the legacy seed mesh."""
    pipeline = bounded_pipeline(bounded)
    role = "gmsh-mesh" if pipeline == GMSH_ONLY_PIPELINE else "seed-mesh"
    return bounded[PIPELINE_BUILD_STAGE[pipeline]]["Artifacts"][role]


def feature_record_bindings(bounded):
    """(restoration recipe binding, build census binding): exactly one per pipeline."""
    if bounded_pipeline(bounded) == GMSH_ONLY_PIPELINE:
        return None, bounded["gmsh-build"]["Artifacts"]["build-census"]
    return bounded["metric-preparation"]["Artifacts"]["restoration-recipe"], None


def _validate_source_transformation(reports, binding, source_paths):
    """Independently validate source-transform inputs, outputs, and metric linkage."""
    pipeline = bounded_pipeline(reports)
    from transform_coupon_source_contract import (
        transform_semantic_contract, transformed_supports, validate_rigid_transform)
    stage = reports["canonical-source-validation"]
    role_names = {"source-semantic-contract": "SemanticContract",
                  "source-signature": "Signature", "source-boundary": "Boundary",
                  "source-mask": "Mask"}
    for stage_name, source_role in role_names.items():
        if stage["Inputs"][stage_name]["SHA256"] != binding["InputSHA256"][source_role]:
            raise ValueError("source-transform input differs from immutable source")
    transform_path = Path(stage["Inputs"]["canonical-transform"]["Path"])
    transform = json.loads(transform_path.read_text())
    if isinstance(transform, dict):
        transform = transform.get("Transform")
    identity = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    if transform != identity:
        raise ValueError("canonical source validation is not source-local identity")
    matrix = validate_rigid_transform(transform)
    semantic_path = Path(stage["Artifacts"]["canonical-semantic-contract"]["Path"])
    supports_path = Path(stage["Artifacts"]["canonical-supports"]["Path"])
    source_semantic = json.loads(Path(source_paths["SemanticContract"]).read_text())
    transform_sha256 = sha256(transform_path)
    semantic_sha256 = sha256(source_paths["SemanticContract"])
    expected_semantic = transform_semantic_contract(source_semantic, matrix)
    expected_semantic["SourceSemanticContractSHA256"] = semantic_sha256
    expected_semantic["CanonicalTransformSHA256"] = transform_sha256
    expected_supports = transformed_supports(
        Path(source_paths["Signature"]).parent, matrix,
        signature=source_paths["Signature"], boundary=source_paths["Boundary"],
        mask=source_paths["Mask"], transform_sha256=transform_sha256,
        semantic_sha256=semantic_sha256)
    actual_semantic = json.loads(semantic_path.read_text())
    actual_supports = json.loads(supports_path.read_text())
    if actual_semantic != expected_semantic or actual_supports != expected_supports:
        raise ValueError("transformed semantic/support artifact differs from source transform")

    build_stage = reports[PIPELINE_BUILD_STAGE[pipeline]]
    seed_inputs = build_stage["Inputs"]
    retained_etch = binding["InputSHA256"].get("RetainedEtch")
    if retained_etch is None:
        if "source-retained-etch" in seed_inputs:
            raise ValueError("seed bound an etch footprint the case does not declare")
    elif seed_inputs.get("source-retained-etch", {}).get("SHA256") != retained_etch:
        raise ValueError("seed did not consume the immutable retained etch footprint")

    from mesh_stage_contract import TRACE_BASIS_INPUTS
    declared_basis = {role: binding["InputSHA256"].get(role) for role in TRACE_BASIS_ROLES}
    basis_stages = [("seed", seed_inputs)]
    if pipeline == LEGACY_MMG_PIPELINE:
        basis_stages.append(("metric", reports["metric-preparation"]["Inputs"]))
    for stage_name, stage_inputs in basis_stages:
        for name, role in TRACE_BASIS_INPUTS.items():
            if declared_basis[role] is None:
                if name in stage_inputs:
                    raise ValueError(f"{stage_name} bound a trace basis the case does not declare")
            elif stage_inputs.get(name, {}).get("SHA256") != declared_basis[role]:
                raise ValueError(f"{stage_name} did not consume the immutable trace basis {role}")
    seed_semantic = Path(build_stage["Inputs"]["canonical-semantic-contract"]["Path"])
    if seed_semantic.resolve() != semantic_path.resolve():
        raise ValueError("volume build did not consume the reconstructed canonical semantic contract")
    if pipeline == LEGACY_MMG_PIPELINE:
        metric = reports["metric-preparation"]
        recipe_path = Path(metric["Artifacts"]["restoration-recipe"]["Path"])
        recipe = json.loads(recipe_path.read_text())
        metric_semantic = Path(metric["Inputs"]["canonical-semantic-contract"]["Path"])
        metric_supports = Path(metric["Inputs"]["canonical-supports"]["Path"])
        if (metric_semantic.resolve() != semantic_path.resolve() or
                metric_supports.resolve() != supports_path.resolve() or
                recipe.get("TransformedSupportsArtifact") != str(supports_path.resolve()) or
                recipe.get("TransformedSupportsSHA256") != sha256(supports_path) or
                recipe.get("SemanticContract") != actual_semantic or
                recipe.get("TransformedSupports") != actual_supports):
            raise ValueError("metric recipe did not consume the reconstructed canonical supports")
    else:
        census = json.loads(Path(build_stage["Artifacts"]["build-census"]["Path"]).read_text())
        if census.get("SemanticContractSHA256") != sha256(semantic_path):
            raise ValueError("build census did not consume the reconstructed canonical semantic contract")

    placement = reports["proper-rigid-publication"]
    placement_roles = {"source-semantic-contract": "SemanticContract",
                       "source-signature": "Signature", "source-boundary": "Boundary",
                       "source-mask": "Mask", "source-process": "Process"}
    if any(placement["Inputs"][name]["SHA256"] != binding["InputSHA256"][role]
           for name, role in placement_roles.items()):
        raise ValueError("placement publication source differs from immutable source")
    placement_transform = Path(placement["Inputs"]["placement-transform"]["Path"])
    placement_value = json.loads(placement_transform.read_text())
    if isinstance(placement_value, dict):
        placement_value = placement_value.get("Transform")
    if placement_value != binding["Transform"]:
        raise ValueError("placement transform file differs from the bound variant")
    placement_matrix = validate_rigid_transform(binding["Transform"])
    expected_semantic = transform_semantic_contract(source_semantic, placement_matrix)
    expected_semantic["SourceSemanticContractSHA256"] = semantic_sha256
    expected_semantic["CanonicalTransformSHA256"] = sha256(placement_transform)
    expected_supports = transformed_supports(
        Path(source_paths["Signature"]).parent, placement_matrix,
        signature=source_paths["Signature"], boundary=source_paths["Boundary"],
        mask=source_paths["Mask"], transform_sha256=sha256(placement_transform),
        semantic_sha256=semantic_sha256)
    final_semantic = json.loads(Path(placement["Artifacts"]
        ["transformed-semantic-contract"]["Path"]).read_text())
    final_supports = json.loads(Path(placement["Artifacts"]["transformed-supports"]["Path"]).read_text())
    if final_semantic != expected_semantic or final_supports != expected_supports:
        raise ValueError("published semantic/support objects differ from independent reconstruction")


def topology_binds_feature_record(topology, bounded):
    """The topology record binds the pipeline's feature record: the legacy
    restoration recipe or the Gmsh-only build census, never both."""
    recipe, census = feature_record_bindings(bounded)
    if recipe is not None:
        return (topology.get("RestorationRecipeSHA256") == recipe["SHA256"] and
                "BuildCensusSHA256" not in topology)
    return (topology.get("BuildCensusSHA256") == census["SHA256"] and
            "RestorationRecipeSHA256" not in topology)


def _validate_bound_records(evidence_path, evidence, binding, source_paths, *, mesh=None, identity_meshes=None):
    """`mesh` is the loaded audit mesh of the variant (read once by the caller);
    `identity_meshes` a dict the caller keeps across the variants of one case, caching
    the loaded canonical (identity) mesh by its SHA-256 so both placements' covariance
    comparisons share it and its memoized kernels (decision 62 step 3, proposal 2)."""
    records = evidence.get("AuditRecords")
    required = {"bounded-run", "mesh-topology-quality", "mesh-complexity",
                "mesh-invariants", "variant-transform"}
    if (not isinstance(records, list) or len(records) != len(required) or
            {item.get("Kind") for item in records} != required):
        raise ValueError("exactly one record of every required audit kind is required")
    measurements = {}
    digests = {evidence["Mesh"]["SHA256"]}
    expected_dependencies = {
        "bounded-run": {},
        "mesh-topology-quality": {
            role: binding["InputSHA256"][role]
            for role in ("SemanticContract", "MeshRecipe", "Process", "Signature")},
        "mesh-complexity": {
            role: binding["InputSHA256"][role]
            for role in ("MeshRecipe", "SemanticContract")},
        "mesh-invariants": {},
        "variant-transform": {},
    }
    records_by_kind = {}
    canonical_stage_digests, placement_stage_digests = set(), set()
    for item in records:
        path = _check_artifact(evidence_path.parent, item, "audit record")
        digest = item["SHA256"]
        if digest in digests:
            raise ValueError("audit records must be content-distinct")
        digests.add(digest)
        record = json.loads(path.read_text())
        records_by_kind[item["Kind"]] = record
        expected = {"Version": 1, "Kind": item["Kind"],
                    "CaseId": binding["CaseId"], "Variant": binding["Variant"],
                    "MeshSHA256": evidence["Mesh"]["SHA256"],
                    "InputSHA256": binding["InputSHA256"],
                    "Transform": binding["Transform"],
                    "TransformSHA256": binding["TransformSHA256"]}
        if any(record.get(key) != value for key, value in expected.items()):
            raise ValueError("audit record has stale or mismatched bindings")
        producer = record.get("Producer", {})
        if (producer.get("Name") not in binding["ToolSHA256"] or
                binding["ToolSHA256"][producer["Name"]] != producer.get("SHA256") or
                not record.get("Command") or not isinstance(record.get("Environment"), dict)):
            raise ValueError("audit record producer/command/environment is not frozen")
        for section, value in record.get("Measurements", {}).items():
            if section in measurements:
                raise ValueError("measurement section has multiple producers")
            measurements[section] = value
        dependencies = record.get("Dependencies", {})
        if dependencies != expected_dependencies[item["Kind"]]:
            raise ValueError("audit record dependencies differ from frozen inputs")
        if item["Kind"] == "variant-transform":
            if (record.get("TransformVerified") is not True or
                    evidence.get("IdentityMeshSHA256") != record.get("IdentityMeshSHA256") or
                    evidence.get("IdentitySeedMeshSHA256") !=
                    record.get("IdentitySeedMeshSHA256") or
                    evidence.get("ParentLabeledMeshSHA256") !=
                    record.get("ParentLabeledMeshSHA256") or
                    evidence.get("RadialShells") != record.get("RadialShells") or
                    evidence.get("TransformMaximumCoordinateError") !=
                    record.get("TransformMaximumCoordinateError")):
                raise ValueError("variant transform result differs from its bound audit")
        if item["Kind"] == "bounded-run":
            stage_items = record.get("BoundedStageRecords")
            if (not isinstance(stage_items, list) or
                    {stage.get("Stage") for stage in stage_items} != set(binding["StageToolSHA256"]) or
                    len(stage_items) != len(binding["StageToolSHA256"])):
                raise ValueError("bounded stage records are incomplete")
            reports, stage_digests = validate_stage_dag(
                {stage["Stage"]: Path(stage["Path"]) for stage in stage_items},
                _artifact_path(evidence_path.parent, evidence["Mesh"]),
                "run_bounded_mesher.py", binding["ToolSHA256"]["run_bounded_mesher.py"],
                binding["StageToolSHA256"])
            if (reports != record.get("BoundedStages") or
                    sorted(stage_digests) != record.get("StageRecordSHA256")):
                raise ValueError("bounded stage DAG differs from its bound producer output")
            if digests & stage_digests:
                raise ValueError("bounded stage and audit artifacts must be content-distinct")
            canonical_stage_digests = {
                sha256(stage["Path"]) for stage in stage_items
                if stage["Stage"] not in PLACEMENT_STAGE_ORDER}
            placement_stage_digests = stage_digests - canonical_stage_digests
            digests.update(stage_digests)
    bounded = records_by_kind["bounded-run"]["BoundedStages"]
    pipeline = bounded_pipeline(bounded)
    placement = bounded["proper-rigid-publication"]
    canonical_record_path = Path(placement["Inputs"]["canonical-build-record"]["Path"])
    canonical_record = json.loads(canonical_record_path.read_text())
    canonical_tools = {
        f"{stage}/{role}": digest
        for stage in canonical_stage_order(pipeline)
        for role, digest in binding["StageToolSHA256"][stage].items()}
    validate_build_record(canonical_record, binding["InputSHA256"], binding["Gates"],
                          canonical_tools)
    if (evidence.get("CanonicalBuildId") != canonical_record["CanonicalBuildId"] or
            evidence.get("CanonicalBuildSHA256") != canonical_record["CanonicalBuildSHA256"] or
            evidence.get("CanonicalArtifactSHA256") != {
                name: item["SHA256"] for name, item in
                canonical_record["CanonicalArtifacts"].items()}):
        raise ValueError("variant canonical-build reference differs from bound build")
    _validate_source_transformation(bounded, binding, source_paths)
    topology = records_by_kind["mesh-topology-quality"]
    publication = bounded["proper-rigid-publication"]["Artifacts"]
    if (topology.get("ReferenceMeshSHA256") != build_volume_binding(bounded)["SHA256"] or
            not topology_binds_feature_record(topology, bounded) or
            topology.get("OwnershipReportSHA256") !=
            publication["ownership-partition"]["SHA256"] or
            topology.get("OwnershipQuadratureSHA256") !=
            publication["ownership-quadrature-partition"]["SHA256"]):
        raise ValueError("topology audit does not bind staged reference/ownership artifacts")
    mesh_path = _artifact_path(evidence_path.parent, evidence["Mesh"])
    from general_mesh_audit_producer import read_audit_mesh
    if mesh is None:
        mesh = read_audit_mesh(mesh_path)
    recomputed = _recompute_mesh_measurements(mesh_path, binding, source_paths, bounded, mesh=mesh)
    recorded_mesh_measurements = {
        key: value for key, value in measurements.items()
        if key not in {"Resources", "PhysicalCovariance"}
    }
    # Variant-transform has no measurement section. Resources are independently
    # derived from the validated stage reports and parsed final tetrahedra below.
    if recomputed != recorded_mesh_measurements:
        raise ValueError("audit measurements differ from an independent producer rerun")
    variant = records_by_kind["variant-transform"]
    identity_path = Path(variant.get("IdentityMeshPath", ""))
    identity_seed_path = Path(variant.get("IdentitySeedMeshPath", ""))
    if (not identity_path.is_file() or sha256(identity_path) !=
            variant.get("IdentityMeshSHA256") or not identity_seed_path.is_file() or
            sha256(identity_seed_path) != variant.get("IdentitySeedMeshSHA256")):
        raise ValueError("variant identity mesh bindings changed")
    from general_mesh_audit_producer import _physical_covariance_report
    matrix = __import__("numpy").asarray(binding["Transform"], dtype=float).reshape(4, 4)
    identity_digest = variant["IdentityMeshSHA256"]
    if identity_meshes is None or identity_digest not in identity_meshes:
        identity_mesh = read_audit_mesh(identity_path)
        if identity_meshes is not None:
            identity_meshes[identity_digest] = identity_mesh
    else:
        identity_mesh = identity_meshes[identity_digest]
    recomputed_physical = _physical_covariance_report(
        identity_mesh, mesh, load_semantic_contract(source_paths["SemanticContract"]), matrix)
    if measurements.get("PhysicalCovariance") != recomputed_physical:
        raise ValueError("physical covariance differs from independent normalization")
    resources = measurements.get("Resources", {})
    canonical_reports = [bounded[name] for name in canonical_stage_order(pipeline)]
    placement_reports = [bounded[name] for name in PLACEMENT_STAGE_ORDER]
    expected_canonical = {
        "Seconds": sum(report["Seconds"] for report in canonical_reports),
        "PeakRSSGiB": max(report["PeakProcessTreeRSSBytes"]
                          for report in canonical_reports) / 2**30}
    expected_placement = {
        "Seconds": sum(report["Seconds"] for report in placement_reports),
        "PeakRSSGiB": max(report["PeakProcessTreeRSSBytes"]
                          for report in placement_reports) / 2**30}
    expected_resources = {
        "ExitCode": 0,
        "Seconds": expected_canonical["Seconds"] + expected_placement["Seconds"],
        "PeakRSSGiB": max(expected_canonical["PeakRSSGiB"],
                           expected_placement["PeakRSSGiB"]),
        "CanonicalBuild": expected_canonical,
        "PlacementPublication": expected_placement,
    }
    if any(resources.get(key) != value for key, value in expected_resources.items()):
        raise ValueError("canonical/placement resource measurements differ from stage reports")
    if (pipeline == GMSH_ONLY_PIPELINE and
            not tube_design_statement(measurements.get("AchievedAnisotropy"))):
        raise ValueError("Gmsh-only evidence lacks the prism tube design statement")

    required_sections = {"Resources", "ActualVolumeMaterials", "ActualBoundaryAttributes",
                         "ActualAdjacency", "OwnershipClosure", "ActualSemanticCorners",
                         "CornerNeighborhoods", "SubdivisionNeighborhoods", "CutNeighborhoods",
                         "ProtectedSurfaces", "AchievedAnisotropy", "TraceDiagonal",
                         "MeshQuality", "Complexity", "ComparisonInvariants",
                         "PhysicalCovariance"}
    if set(measurements) != required_sections:
        raise ValueError("bound records do not supply the exact measurement schema")
    for section, value in measurements.items():
        if evidence.get(section) != value:
            raise ValueError("normalized measurements differ from bound producer records")
    variant_digests = (digests - canonical_stage_digests -
                       {evidence["Mesh"]["SHA256"]})
    return variant_digests, canonical_stage_digests, canonical_record


def _bounded_stages(evidence_path, evidence):
    """The recorded stage reports of the evidence's bounded-run record."""
    item = next(record for record in evidence["AuditRecords"] if record.get("Kind") == "bounded-run")
    record = json.loads(_check_artifact(evidence_path.parent, item, "audit record").read_text())
    return record["BoundedStages"]


def audit_manifest_evidence(evidence, gates, contract, binding):
    """Judge actual measurements against a separately frozen semantic contract."""
    failures = []
    required_binding = {
        "CaseId": binding["CaseId"], "Variant": binding["Variant"],
        "TransformSHA256": binding["TransformSHA256"],
        "InputSHA256": binding["InputSHA256"],
        "ProcessSHA256": binding["InputSHA256"]["Process"],
        "SemanticContractSHA256": binding["InputSHA256"]["SemanticContract"],
        "RecipeSHA256": binding["InputSHA256"]["MeshRecipe"],
        "ToolSHA256": binding["ToolSHA256"],
        "StageToolSHA256": binding["StageToolSHA256"],
    }
    if evidence.get("Version") != 3 or any(evidence.get(key) != value
                                            for key, value in required_binding.items()):
        failures.append("provenance-binding")

    expected_materials = sorted(contract["VolumeMaterials"], key=lambda item: item["Attribute"])
    actual_materials = evidence.get("ActualVolumeMaterials")
    expected_labels = sorted(item["Attribute"] for item in contract["BoundaryLabels"])
    actual_labels = evidence.get("ActualBoundaryAttributes")
    if (not actual_materials or not actual_labels or
            sorted(actual_materials, key=lambda item: item.get("Attribute", -1)) != expected_materials or
            sorted(actual_labels) != expected_labels):
        failures.append("exact-labels-materials")
    expected_adjacency = {str(item["Attribute"]): sorted(item["AdjacentMaterials"])
                          for item in contract["BoundaryLabels"]}
    actual_adjacency = evidence.get("ActualAdjacency")
    if (not actual_adjacency or
            {str(key): sorted(value) for key, value in actual_adjacency.items()} !=
            expected_adjacency):
        failures.append("material-adjacency")

    ownership = evidence.get("OwnershipClosure", {})
    physical_coverage = ownership.get("PhysicalSurfaceCoverage", {})
    response_ownership = ownership.get("ResponseOwnership", {})
    diagnostics = ownership.get("WholeElementAmbiguityDiagnostics", {})
    if (physical_coverage.get("Complete") is not True or
            response_ownership.get("UnmatchedPolicy") != contract["UnmatchedPolicy"] or
            response_ownership.get("PositiveWeights") is not True or
            response_ownership.get("UnmatchedPoints") != 0 or
            response_ownership.get("OverlappingPoints") != 0 or
            response_ownership.get("Exhaustive") is not True or
            response_ownership.get("RelativeClosureError", float("inf")) >
            response_ownership.get("ClosureTolerance", -1.0) or
            response_ownership.get("OwnerPartitionRelativeClosure", float("inf")) > 1e-12 or
            response_ownership.get("NoDuplicateOrMissingOwners") is not True or
            response_ownership.get("ExpectedOwnerAttributes") !=
            response_ownership.get("OwnerAttributes") or
            diagnostics.get("AuthoritativeForResponseOwnership") is not False):
        failures.append("ownership-exhaustive-closure")
    expected_corners = _transform_points(contract["SemanticCorners"], binding["Transform"])
    if not _same_points(expected_corners, evidence.get("ActualSemanticCorners"),
                        float(gates["CornerTolerance"])):
        failures.append("semantic-corners")
    def neighborhood_failure(name, expected, maximum=None, minimum=None):
        values = evidence.get(name)
        if not isinstance(values, list) or len(values) != len(expected):
            return True
        if not expected:
            return False
        if not _same_points(expected, [item.get("Point") for item in values],
                            float(gates["CornerTolerance"])):
            return True
        aspects = [item.get("MaximumAspect") for item in values]
        return (any(not _finite_number(value, positive=True) for value in aspects) or
                (maximum is not None and any(value > maximum for value in aspects)) or
                (minimum is not None and any(value < minimum for value in aspects)))
    topology = contract["FeatureTopology"]
    subdivisions = _transform_points(topology["CADSubdivisionEndpoints"], binding["Transform"])
    cuts = _transform_points(topology["CutEndpoints"], binding["Transform"])
    if (neighborhood_failure("CornerNeighborhoods", expected_corners,
                             maximum=gates["MaximumCornerAspect"]) or
            neighborhood_failure("SubdivisionNeighborhoods", subdivisions,
                                 minimum=gates["MinimumNoncornerAspect"]) or
            neighborhood_failure("CutNeighborhoods", cuts,
                                 minimum=gates["MinimumNoncornerAspect"])):
        failures.append("semantic-corner-and-endpoint-anisotropy")
    protected = evidence.get("ProtectedSurfaces", {})
    if (not protected.get("Actual") or protected.get("PlaneSupportsMatch") is not True or
            protected.get("TopologyMatches") is not True or
            sorted(protected.get("Actual", [])) != sorted(contract["ProtectedSupports"]) or
            not _finite_number(protected.get("MaximumRelativeMeasureError"),
                               nonnegative=True) or
            protected.get("MaximumRelativeMeasureError", math.inf) >
            gates["MaximumProtectedMeasureError"] or
            not _finite_number(protected.get("MaximumSupportVertexDistance"),
                               nonnegative=True) or
            protected.get("MaximumSupportVertexDistance", math.inf) >
            gates["CornerTolerance"]):
        failures.append("protected-surfaces")

    widths = evidence.get("AchievedAnisotropy", {})
    if not layer_covered_band(widths) and not tube_design_statement(widths):
        values = [widths.get(name) for name in ("Transverse1P90", "Transverse2P90",
                                                 "NormalTarget", "TangentialP50")]
        if (widths.get("Gate") != ANISOTROPY_GATE_APPLIED or
                not isinstance(widths.get("Samples"), int) or widths.get("Samples", 0) <= 0 or
                any(not _finite_number(x, positive=True) for x in values) or
                max(values[:2]) > gates["MaximumNormalFactor"] * values[2] or
                values[3] < gates["MinimumAchievedAspect"] * max(values[:2])):
            failures.append("achieved-anisotropy")
    if evidence.get("TraceDiagonal", {}).get("GlobalDiagonalBands") != 0:
        failures.append("trace-diagonal-overrefinement")

    quality = evidence.get("MeshQuality", {})
    judged = quality
    rule = gates.get(EDGE_LAYER_QUALITY_RULE_GATE)
    layer = quality.get("EdgeLayer") if isinstance(quality, dict) else None
    if rule is not None and isinstance(layer, dict) and layer.get("Cells", 0) > 0:
        # Decision 32 (calibration only): the recorded layer's cells are judged by
        # orientation above the roundoff floor and the edge-aspect bound; every
        # other quality gate judges the cells outside the layer.  Without the rule
        # the whole mesh (layer included) is judged as before.
        judged = quality.get("OutsideEdgeLayer", {})
        if (not isinstance(judged.get("Samples"), int) or
                judged["Samples"] + layer["Cells"] != quality.get("Samples") or
                layer.get("PositiveOrientation") is not True or
                not _finite_number(layer.get("MinimumScaledJacobian")) or
                layer["MinimumScaledJacobian"] <= rule["ScaledJacobianRoundoffFloor"] or
                not _finite_number(layer.get("MaximumEdgeAspect"), positive=True) or
                layer["MaximumEdgeAspect"] > rule["MaximumEdgeAspect"] or
                not _finite_number(layer.get("Reach"), positive=True) or
                layer.get("Rule") != EDGE_LAYER_QUALITY_RULE):
            failures.append("edge-layer-quality")
    if (not isinstance(judged.get("Samples"), int) or judged.get("Samples", 0) <= 0 or
            judged.get("PositiveOrientation") is not True or
            not _finite_number(judged.get("MinimumScaledJacobian"), nonnegative=True) or
            judged.get("MinimumScaledJacobian", -1) < gates["MinimumScaledJacobian"] or
            not _finite_number(judged.get("MaximumJacobianCondition"), positive=True) or
            judged.get("MaximumJacobianCondition", math.inf) > gates["MaximumJacobianCondition"] or
            not _per_type_quality_passes(quality, gates)):
        failures.append("mesh-quality-jacobian")
    resources = evidence.get("Resources", {})
    resource_names = ("Seconds", "PeakRSSGiB", "Elements")
    canonical_resources = resources.get("CanonicalBuild", {})
    placement_resources = resources.get("PlacementPublication", {})
    separated = all(_finite_number(item.get(name), nonnegative=True)
                    for item in (canonical_resources, placement_resources)
                    for name in ("Seconds", "PeakRSSGiB"))
    if (resources.get("ExitCode") != 0 or not separated or
            any(not _finite_number(resources.get(name), nonnegative=True) for name in resource_names) or
            canonical_resources.get("Seconds", math.inf) > gates["MaximumSeconds"] or
            canonical_resources.get("PeakRSSGiB", math.inf) > gates["MaximumRSSGiB"] or
            placement_resources.get("Seconds", math.inf) > gates["MaximumSeconds"] or
            placement_resources.get("PeakRSSGiB", math.inf) > gates["MaximumRSSGiB"] or
            resources.get("Elements", math.inf) > gates["MaximumElements"]):
        failures.append("bounded-resources")
    complexity = evidence.get("Complexity", {})
    if (not isinstance(complexity.get("H1DOFs"), int) or complexity.get("H1DOFs", 0) <= 0 or
            not isinstance(complexity.get("FeatureCount"), int) or
            complexity.get("FeatureCount", 0) <= 0 or
            not isinstance(complexity.get("CADSubdivisionCount"), int) or
            complexity.get("CADSubdivisionCount", -1) < 0):
        failures.append("complexity-counts")
    invariants = evidence.get("ComparisonInvariants")
    if (not isinstance(invariants, dict) or not invariants or
            any(not _finite_number(value) for value in invariants.values())):
        failures.append("comparison-invariants")
    physical = evidence.get("PhysicalCovariance", {})
    if (physical.get("ComparisonFrame") != "SourceLocal" or
            physical.get("LabelsMaterialsAdjacencyMatch") is not True):
        failures.append("physical-covariance-contract")
    return failures


def _per_type_quality_passes(quality, gates):
    """Mixed-element quality (MeshQuality.ByType, Gmsh-only pipeline): every volume
    element type is positively oriented within the Jacobian condition gate and the
    tetrahedra are within the scaled-Jacobian gate; the top-level record must be
    the aggregate of the types.  A record without ByType (legacy tetrahedral
    evidence) passes this check and is judged by the top-level values alone."""
    by_type = quality.get("ByType") if isinstance(quality, dict) else None
    if by_type is None:
        return True
    if not isinstance(by_type, dict) or "Tetrahedron" not in by_type or not by_type:
        return False
    for name, record in by_type.items():
        if (not isinstance(record, dict) or not isinstance(record.get("Samples"), int) or
                record["Samples"] <= 0 or record.get("PositiveOrientation") is not True or
                record.get("NonpositiveCells") != 0 or
                not _finite_number(record.get("MaximumJacobianCondition"), positive=True) or
                record["MaximumJacobianCondition"] > gates["MaximumJacobianCondition"] or
                not _finite_number(record.get("MinimumScaledJacobian"), nonnegative=True)):
            return False
    tetrahedra = by_type["Tetrahedron"]
    if tetrahedra["MinimumScaledJacobian"] < gates["MinimumScaledJacobian"]:
        return False
    return (quality.get("Samples") == sum(record["Samples"] for record in by_type.values()) and
            quality.get("MaximumJacobianCondition") ==
            max(record["MaximumJacobianCondition"] for record in by_type.values()) and
            quality.get("MinimumScaledJacobian") == tetrahedra["MinimumScaledJacobian"])


def _relative_error(a, b):
    scale = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / scale


ANISOTROPY_STATISTICS = ("TangentialP50", "Transverse1P90", "Transverse2P90")
# The prism tube design statement compared between placements of one canonical
# build (all recorded by the build; equal for every placement).
TUBE_DESIGN_STATISTICS = ("InnerSize", "GrowthRatio", "TangentialSize", "SpacingMinimum",
                          "SpacingMaximum", "Prisms", "Pyramids", "MeshPrisms", "MeshPyramids")


def tube_design_statement(widths):
    """True when the achieved-anisotropy design gate is replaced by the prism tube
    design statement (Gmsh-only pipeline, supervisor decision 38): the record
    declares TUBE_DESIGN_GATE with Samples 0, the recorded tube design (finite
    positive inner size / ratio > 1 / spacing within the tangential size, rings >= 1,
    positive prism and pyramid counts, the decision-40 layer record: the layer rule
    and axis size law named, LayerGrowthCap equal to the growth ratio, the layer
    thickness Minimum <= P50 <= Maximum <= TangentialSize equal to the spacing
    extremes, the neighbour ratio within the cap and the count of layers below
    TangentialSize / GrowthRatio within the layer count) and the mesh's prism and pyramid
    counts equal to the census (CensusMatchesMesh).  Anything less is judged by the
    anisotropy gate as an ordinary band sample (and fails on Samples 0)."""
    from general_mesh_audit_producer import TUBE_DESIGN_GATE, TUBE_DESIGN_RULE
    if not isinstance(widths, dict) or widths.get("Gate") != TUBE_DESIGN_GATE:
        return False
    if (widths.get("Samples") != 0 or widths.get("Rule") != TUBE_DESIGN_RULE or
            widths.get("CensusMatchesMesh") is not True or
            not _finite_number(widths.get("InnerSize"), positive=True) or
            not _finite_number(widths.get("GrowthRatio"), positive=True) or
            widths["GrowthRatio"] <= 1.0 or
            not _finite_number(widths.get("TangentialSize"), positive=True) or
            not _finite_number(widths.get("NormalTarget"), positive=True) or
            not _finite_number(widths.get("SpacingMinimum"), positive=True) or
            not _finite_number(widths.get("SpacingMaximum"), positive=True) or
            not widths["SpacingMinimum"] <= widths["SpacingMaximum"] <= widths["TangentialSize"] or
            not isinstance(widths.get("Rings"), int) or widths["Rings"] < 1 or
            not isinstance(widths.get("RingSizes"), list) or len(widths["RingSizes"]) != widths["Rings"] or
            any(not isinstance(widths.get(name), int) or widths[name] <= 0
                for name in ("TubeCount", "Layers", "Prisms", "Pyramids", "MeshPrisms", "MeshPyramids")) or
            widths["MeshPrisms"] != widths["Prisms"] or widths["MeshPyramids"] != widths["Pyramids"]):
        return False
    layers = widths.get("LayerThickness")
    if (not all(isinstance(widths.get(name), str) and widths[name]
                for name in ("LayerRule", "TubeAxisSizeLaw")) or
            not _finite_number(widths.get("LayerGrowthCap"), positive=True) or
            widths["LayerGrowthCap"] != widths["GrowthRatio"] or
            not isinstance(layers, dict) or
            any(not _finite_number(layers.get(name), positive=True)
                for name in ("Minimum", "P50", "Maximum", "MaximumNeighbourRatio")) or
            not layers["Minimum"] <= layers["P50"] <= layers["Maximum"] <= widths["TangentialSize"] or
            layers["Minimum"] != widths["SpacingMinimum"] or layers["Maximum"] != widths["SpacingMaximum"] or
            layers["MaximumNeighbourRatio"] > widths["LayerGrowthCap"] or
            isinstance(layers.get("LayersBelowTangentialSizeOverGrowthRatio"), bool) or
            not isinstance(layers.get("LayersBelowTangentialSizeOverGrowthRatio"), int) or
            not 0 <= layers["LayersBelowTangentialSizeOverGrowthRatio"] <= widths["Layers"]):
        return False
    return True


def layer_covered_band(widths):
    """True when the achieved-anisotropy design gate is not applicable by
    construction (supervisor decision 35): the record declares
    ANISOTROPY_GATE_NOT_APPLICABLE with no band cell outside the recorded edge layer
    within one NormalSize (Samples 0, ExcludedEdgeLayerCells > 0), a layer record
    with cells and finite statistics, and a layer-adjacent band record with cells,
    finite statistics and its transverse P90 over NormalSize (the
    MaximumNormalFactor-equivalent, reported and never gated).  Anything less is
    judged by the gate as an ordinary band sample (and fails on Samples 0)."""
    if not isinstance(widths, dict) or widths.get("Gate") != ANISOTROPY_GATE_NOT_APPLICABLE:
        return False
    layer, adjacent = widths.get("EdgeLayer"), widths.get("LayerAdjacentBand")
    normal = widths.get("NormalTarget")
    if (widths.get("Samples") != 0 or
            not isinstance(widths.get("ExcludedEdgeLayerCells"), int) or
            widths["ExcludedEdgeLayerCells"] <= 0 or
            not _finite_number(normal, positive=True) or
            widths.get("DistanceCutoff") != normal * BAND_GATE_CUTOFFS[0] or
            not isinstance(layer, dict) or not isinstance(adjacent, dict) or
            any(not isinstance(item.get("Cells"), int) or item["Cells"] <= 0 or
                any(not _finite_number(item.get(name), positive=True)
                    for name in ANISOTROPY_STATISTICS)
                for item in (layer, adjacent)) or
            adjacent.get("DistanceCutoff") != normal * BAND_GATE_CUTOFFS[-1] or
            adjacent.get("Rule") != LAYER_ADJACENT_BAND_RULE or
            not _finite_number(adjacent.get("TransverseP90OverNormalSize"), positive=True) or
            adjacent["TransverseP90OverNormalSize"] !=
            max(adjacent["Transverse1P90"], adjacent["Transverse2P90"]) / normal or
            not isinstance(adjacent.get("NearestSpanVertexDistance"), dict) or
            any(not _finite_number(adjacent["NearestSpanVertexDistance"].get(name), nonnegative=True)
                for name in ("Minimum", "Maximum"))):
        return False
    return True


def _compared_anisotropy(evidence):
    """The statistics the covariance comparison uses: the gated band sample, the
    layer-adjacent band when the gate is not applicable, or the tube design
    statement (mapped onto the comparison keys: the recorded tube statistics)."""
    widths = evidence.get("AchievedAnisotropy", {})
    if layer_covered_band(widths):
        return dict(widths["LayerAdjacentBand"], Gate=widths["Gate"])
    if tube_design_statement(widths):
        return {"Gate": widths["Gate"],
                "Statistics": [float(widths[name]) for name in TUBE_DESIGN_STATISTICS]}
    return dict(widths, Gate=widths.get("Gate") if isinstance(widths, dict) else None)


def _physical_comparison_failures(reference_evidence, transformed_evidence, comparison):
    """Compare final meshes physically; deterministic topology is diagnostic only."""
    failures = []
    physical = transformed_evidence.get("PhysicalCovariance", {})
    if physical.get("LabelsMaterialsAdjacencyMatch") is not True:
        failures.append("physical labels/material adjacency")
        return failures
    left, right = (physical.get(name, {}) for name in
                   ("ReferenceInvariants", "TransformedInvariants"))
    if left.keys() != right.keys() or not left:
        return ["physical measure keys differ"]
    volume_error = max((_relative_error(left[name], right[name]) for name in left
                        if name.startswith("Volume:")), default=math.inf)
    area_error = max((_relative_error(left[name], right[name]) for name in left
                      if name.startswith("Area:")), default=math.inf)
    if volume_error > comparison["MaximumRelativeVolumeError"]:
        failures.append("material volumes")
    if area_error > comparison["MaximumRelativeSurfaceMeasureError"]:
        failures.append("boundary surface measures")
    protected = physical.get("ProtectedSurfaces", {})
    if (protected.get("PlaneSupportsMatch") is not True or
            protected.get("TopologyMatches") is not True or
            protected.get("MaximumSupportVertexDistance", math.inf) >
            comparison["MaximumProtectedSupportHausdorff"] or
            protected.get("MaximumRelativeMeasureError", math.inf) >
            comparison["MaximumProtectedMeasureError"]):
        failures.append("protected support geometry/topology")
    qualities = [physical.get(name, {}) for name in
                 ("ReferenceQuality", "TransformedQuality")]
    quality_values = []
    for name in ("ScaledJacobianQuantiles", "JacobianConditionQuantiles"):
        if (not all(isinstance(item.get(name), list) for item in qualities) or
                len(qualities[0].get(name, [])) != len(qualities[1].get(name, []))):
            quality_values = [math.inf]
            break
        quality_values.extend(_relative_error(a, b)
                              for a, b in zip(qualities[0][name], qualities[1][name]))
    if (not all(item.get("PositiveOrientation") is True for item in qualities) or
            max(quality_values, default=math.inf) >
            comparison["MaximumQualityDistributionRelativeError"]):
        failures.append("orientation/quality distribution")
    anisotropy = [_compared_anisotropy(item)
                  for item in (reference_evidence, transformed_evidence)]
    statistics = [item["Statistics"] if "Statistics" in item else
                  [item.get(name) for name in ANISOTROPY_STATISTICS] for item in anisotropy]
    if (anisotropy[0]["Gate"] != anisotropy[1]["Gate"] or
            any(not _finite_number(value) for values in statistics for value in values) or
            max(_relative_error(a, b) for a, b in zip(*statistics)) >
            comparison["MaximumAnisotropyRelativeError"]):
        failures.append("local-frame anisotropy")
    complexity_values = [
        [item["Complexity"]["H1DOFs"], item["Resources"]["Elements"]]
        for item in (reference_evidence, transformed_evidence)]
    complexity_ratio = max(max(a, b) / max(min(a, b), 1)
                           for a, b in zip(*complexity_values))
    if complexity_ratio > comparison["MaximumComplexityRatio"]:
        failures.append("complexity ratio")
    return failures


def run_manifest(args):
    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text())
    repository, tool_hashes, _ = validate_manifest(manifest, manifest_path)
    cases_by_id = {case["Id"]: case for case in manifest["Cases"]}
    overrides = {}
    for item in args.input:
        if "=" not in item:
            raise ValueError("--input must be CASE=DIRECTORY")
        key, value = item.split("=", 1)
        if key in overrides or key not in cases_by_id:
            raise ValueError("Duplicate or unknown input override: " + key)
        overrides[key] = Path(value).resolve()

    records, sources, preflight_ok = [], {}, True
    for case in manifest["Cases"]:
        record = {"Id": case["Id"], "Passed": False,
                  "Variants": [variant["Id"] for variant in case["Variants"]]}
        try:
            source = case["Source"]
            directory = overrides.get(case["Id"])
            if directory is None and source.get("Directory"):
                candidate = Path(source["Directory"])
                directory = candidate if candidate.is_absolute() else repository / candidate
            if directory is None or not directory.is_dir():
                raise ValueError("required immutable input directory is unavailable")
            hashes, paths = {}, {}
            for role, entry in source["Files"].items():
                expected, name = entry.get("SHA256"), entry.get("Name")
                repository_name = entry.get("RepositoryPath")
                if not expected or (not name and not repository_name) or (name and repository_name):
                    raise ValueError(f"{role} has no unambiguous frozen path and SHA256")
                candidate = Path(repository_name or name)
                if candidate.is_absolute():
                    path = candidate
                elif repository_name:
                    path = repository / candidate
                else:
                    path = directory / candidate
                if not path.is_file() or sha256(path) != expected:
                    raise ValueError(f"immutable {role} hash mismatch")
                hashes[role], paths[role] = expected, path
            contract = load_semantic_contract(paths["SemanticContract"])
            validate_feature_topology(contract, paths["Signature"], paths["Boundary"])
            signature_role = source["SignatureRole"]
            with paths[signature_role].open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            required_columns = set(source["SignatureColumns"])
            if not rows or not required_columns.issubset(rows[0]):
                raise ValueError("empty or malformed edge signature")
            record.update({"InputDirectory": str(directory), "InputSHA256": hashes,
                           "DiscoveredEdgeCount": len(rows),
                           "DiscoveredSlots": sorted({int(row["Slot"]) for row in rows}),
                           "DiscoveredConductors": sorted({int(row["Conductor"]) for row in rows})})
            scope = preflight_recipe_scope(manifest, paths, case)
            if scope is not None:
                record[SCOPE_KEY] = scope
                if scope["UnsupportedClasses"]:
                    raise UnsupportedClassError(scope["UnsupportedClasses"][0])
            cost = preflight_build_cost(manifest, manifest_path, case)
            if cost is not None:
                record[BUILD_COST_ESTIMATE_KEY] = {
                    key: cost[key] for key in ("EstimatedTetrahedra", "EstimatedPrisms", "EstimatedPyramids",
                                               "EstimatedElements", "MaximumElements", "EstimateOverCap",
                                               "Integrals", "TetrahedraPerCubicSize", "Passed")}
                if not cost["Passed"]:
                    raise ValueError(f"pre-build element estimate {cost['EstimatedElements']:.0f} exceeds "
                                     f"MaximumElements {cost['MaximumElements']} (headroom gate, fail closed)")
            sources[case["Id"]] = (hashes, paths, contract)
            record["Passed"] = True
        except UnsupportedClassError as error:
            # An unsupported class is not a preflight failure of the matrix: the case is
            # recorded as such and never built or judged (decision 48).
            record["Error"] = str(error)
            record[UNSUPPORTED_CLASS_KEY] = error.guard
        except (KeyError, OSError, TypeError, ValueError) as error:
            record["Error"] = str(error)
            preflight_ok = False
        records.append(record)

    unsupported = {record["Id"]: record[UNSUPPORTED_CLASS_KEY] for record in records
                   if UNSUPPORTED_CLASS_KEY in record}
    summary = {"Version": 2, "Scope": "Mesh-only geometry-independence gates",
               "Manifest": str(manifest_path), "PreflightPassed": preflight_ok,
               "UnsupportedClassCases": unsupported,
               "Cases": records, "Passed": False}
    args.root.mkdir(parents=True, exist_ok=False)
    if not preflight_ok or args.preflight_only:
        summary["Passed"] = preflight_ok and args.preflight_only
        (args.root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        return summary["Passed"]
    if args.audit_root is None:
        raise ValueError("--audit-root is required unless --preflight-only is used")

    evidence_by_key, used_variant_audits, used_meshes = {}, set(), set()
    canonical_reuse = {}
    for case, record in zip(manifest["Cases"], records):
        record["VariantResults"] = []
        if UNSUPPORTED_CLASS_KEY in record:
            continue
        hashes, _, contract = sources[case["Id"]]
        for variant in case["Variants"]:
            variant_id = variant["Id"]
            result = {"Variant": variant_id, "Passed": False}
            path = args.audit_root / f"{case['Id']}--{variant_id}.json"
            try:
                evidence = json.loads(path.read_text())
                binding = {"CaseId": case["Id"], "Variant": variant_id,
                           "Transform": variant["Transform"],
                           "TransformSHA256": canonical_sha256(variant["Transform"]),
                           "InputSHA256": hashes, "ToolSHA256": tool_hashes,
                           "StageToolSHA256": manifest["StageToolSHA256"],
                           "Gates": manifest["Gates"]}
                failures = audit_manifest_evidence(evidence, case_gates(manifest, case), contract,
                                                   binding)
                mesh_path = _check_artifact(path.parent, evidence.get("Mesh"), "audited mesh")
                _validate_mesh(mesh_path, contract)
                mesh_digest = evidence["Mesh"]["SHA256"]
                if mesh_digest in used_meshes:
                    raise ValueError("audited meshes must be content-distinct per matrix entry")
                used_meshes.add(mesh_digest)
                variant_digests, canonical_digests, canonical_record = _validate_bound_records(
                    path, evidence, binding, sources[case["Id"]][1])
                validate_production_recipe_commands(manifest, case, _bounded_stages(path, evidence))
                if used_variant_audits & variant_digests:
                    raise ValueError("variant audit/placement records must be content-distinct")
                used_variant_audits.update(variant_digests)
                build_id = canonical_record["CanonicalBuildId"]
                previous = canonical_reuse.get(build_id)
                if previous is None:
                    if any(canonical_digests & item[1] for item in canonical_reuse.values()):
                        raise ValueError("canonical stages were reused under a different cache key")
                    canonical_reuse[build_id] = (canonical_record, canonical_digests)
                elif (not same_canonical_build(previous[0], canonical_record) or
                      previous[1] != canonical_digests):
                    raise ValueError("shared canonical stages require exact cache key and hashes")
                result.update({"AuditEvidence": str(path), "GateFailures": failures,
                               "Passed": not failures})
                evidence_by_key[(case["Id"], variant_id)] = evidence
            except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError) as error:
                result["Error"] = str(error)
            record["VariantResults"].append(result)
        record["Passed"] = all(item["Passed"] for item in record["VariantResults"])

    comparison_failures = []
    for case in manifest["Cases"]:
        comparison = case["TransformComparison"]
        keys = [(case["Id"], comparison[name]) for name in ("Reference", "Transformed")]
        if any(key not in evidence_by_key for key in keys):
            comparison_failures.append(case["Id"] + ": missing transform evidence")
            continue
        reference_evidence, transformed_evidence = (evidence_by_key[key] for key in keys)
        identity_digest = parent_labeled_mesh_sha256(reference_evidence)
        identity_seed_digest = reference_evidence.get("IdentitySeedMeshSHA256")
        coordinate_error = transformed_evidence.get("TransformMaximumCoordinateError")
        if (reference_evidence.get("IdentityMeshSHA256") != identity_digest or
                transformed_evidence.get("IdentityMeshSHA256") != identity_digest or
                transformed_evidence.get("IdentitySeedMeshSHA256") != identity_seed_digest or
                not _finite_number(coordinate_error, nonnegative=True) or
                coordinate_error > manifest["Gates"]["CornerTolerance"]):
            comparison_failures.append(case["Id"] + ": exact source-seed covariance")
            continue
        comparison_failures.extend(
            case["Id"] + ": " + failure for failure in _physical_comparison_failures(
                reference_evidence, transformed_evidence, comparison))
    scaling_failures = []
    for comparison in manifest["ScalingComparisons"]:
        keys = [tuple(comparison[name]) for name in ("Reference", "Compared")]
        if any(key not in evidence_by_key for key in keys):
            scaling_failures.append(comparison["Id"] + ": missing evidence")
            continue
        complexity = [evidence_by_key[key]["Complexity"] for key in keys]
        normalized = [item["H1DOFs"] / item["FeatureCount"] for item in complexity]
        ratio = max(normalized) / min(normalized)
        if (comparison["Kind"] == "cad-subdivision-sensitivity" and
                (complexity[0]["FeatureCount"] != complexity[1]["FeatureCount"] or
                 complexity[0]["CADSubdivisionCount"] == complexity[1]["CADSubdivisionCount"])):
            scaling_failures.append(comparison["Id"] + ": invalid CAD subdivision control")
        elif ratio > comparison["MaximumNormalizedDOFRatio"]:
            scaling_failures.append(comparison["Id"] + ": H1 DOF scaling")
    summary["TransformComparisonFailures"] = comparison_failures
    summary["ScalingComparisonFailures"] = scaling_failures
    summary["Passed"] = (all(record["Passed"] for record in records) and
                         not comparison_failures and not scaling_failures)
    (args.root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary["Passed"]
