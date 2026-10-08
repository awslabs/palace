#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Validation for reusable canonical-build and per-placement publication DAGs."""
import csv
import hashlib
import json
import math
from pathlib import Path
import re

from edge_volume_metric import (COPLANAR_TOLERANCE, EDGE_LAYER_ORIENTATION_FLOOR,
                                EDGE_LAYER_QUALITY_RULE)
from trace_basis import NEEDLE_ALTITUDE_OVER_SHORTEST_EDGE
from semantic_mesh_contract import (boundary_attributes, contract_is_corner_free, cut_surface_attributes,
                                    invariant_corners, material_interface_attributes)


# Two canonical-build pipelines share the source validation, the Gmsh publication
# and the rigid placement.  The production pipeline (supervisor decision 38) is
# Gmsh-only: the mesher's prism-tube build is the canonical volume mesh.  The
# legacy MMG pipeline (seed -> metric -> MMG -> label restoration) is retired from
# production and remains available under the labeled calibration manifest.  A
# pipeline is identified by its exact stage set (`pipeline_of`); the module
# constants below name the legacy order for the tools that predate decision 38.
GMSH_ONLY_PIPELINE = "gmsh-only"
LEGACY_MMG_PIPELINE = "legacy-mmg"
PIPELINE_CANONICAL_STAGES = {
    GMSH_ONLY_PIPELINE: ("canonical-source-validation", "gmsh-build",
                         "canonical-gmsh-publication"),
    LEGACY_MMG_PIPELINE: ("canonical-source-validation", "seed-generation", "metric-preparation",
                          "native-adaptation-mmg", "label-restoration",
                          "canonical-gmsh-publication"),
}
PLACEMENT_STAGE_ORDER = ("proper-rigid-publication",)
CANONICAL_STAGE_ORDER = PIPELINE_CANONICAL_STAGES[LEGACY_MMG_PIPELINE]
STAGE_ORDER = CANONICAL_STAGE_ORDER + PLACEMENT_STAGE_ORDER
# The volume mesh stage of every pipeline and the volume artifact the publication
# consumes from it.
PIPELINE_BUILD_STAGE = {GMSH_ONLY_PIPELINE: "gmsh-build", LEGACY_MMG_PIPELINE: "seed-generation"}
PIPELINE_PUBLISHED_VOLUME = {GMSH_ONLY_PIPELINE: "gmsh-mesh", LEGACY_MMG_PIPELINE: "restored-mesh"}


def canonical_stage_order(pipeline):
    if pipeline not in PIPELINE_CANONICAL_STAGES:
        raise ValueError(f"Unknown canonical pipeline: {pipeline}")
    return PIPELINE_CANONICAL_STAGES[pipeline]


def stage_order(pipeline):
    return canonical_stage_order(pipeline) + PLACEMENT_STAGE_ORDER


def pipeline_of(stages, *, canonical_only=False):
    """The pipeline whose exact canonical (+ placement) stage set is `stages`; fails
    closed on any other set."""
    names = set(stages)
    for pipeline in PIPELINE_CANONICAL_STAGES:
        expected = set(canonical_stage_order(pipeline) if canonical_only else stage_order(pipeline))
        if names == expected:
            return pipeline
    raise ValueError("Stage set differs from every canonical pipeline: " + ", ".join(sorted(names)))


def pipeline_of_tool_roles(roles):
    """The pipeline whose canonical `stage/role` tool names are exactly `roles`."""
    names = set(roles)
    for pipeline in PIPELINE_CANONICAL_STAGES:
        expected = {f"{stage}/{role}" for stage in canonical_stage_order(pipeline)
                    for role in STAGE_TOOLS[stage]}
        if names == expected:
            return pipeline
    raise ValueError("Canonical tool roles differ from the exact canonical schema of every pipeline")


def canonical_artifact_roles(pipeline):
    return frozenset(role for stage in canonical_stage_order(pipeline)
                     for role in stage_outputs(stage, pipeline))


def canonical_tool_roles(pipeline):
    return frozenset(f"{stage}/{role}" for stage in canonical_stage_order(pipeline)
                     for role in STAGE_TOOLS[stage])


STAGE_TOOLS = {
    "canonical-source-validation": {"runtime", "source-validator"},
    "seed-generation": {"runtime", "mesher"},
    "gmsh-build": {"runtime", "mesher"},
    "metric-preparation": {"runtime", "metric-preparer"},
    "native-adaptation-mmg": {"runtime", "adaptation-wrapper", "adapter-mmg", "mmg-library"},
    "label-restoration": {"runtime", "label-restorer"},
    "canonical-gmsh-publication": {"runtime", "publisher"},
    "proper-rigid-publication": {"runtime", "rigid-publisher", "ownership-runtime",
                                 "ownership-auditor"},
}
STAGE_INPUTS = {
    "canonical-source-validation": {"source-semantic-contract", "source-signature",
                                    "source-boundary", "source-mask", "canonical-transform"},
    "seed-generation": {"source-signature", "source-boundary", "source-mask",
                        "canonical-semantic-contract"},
    "gmsh-build": {"source-signature", "source-boundary", "source-mask",
                   "canonical-semantic-contract"},
    "metric-preparation": {"seed-mesh", "seed-corner-census", "canonical-semantic-contract",
                           "canonical-supports"},
    "native-adaptation-mmg": {"mmg-seed", "metric", "pins", "fixed-triangles",
                              "required-tetrahedra", "restoration-recipe"},
    "label-restoration": {"adapted-mesh", "restoration-recipe"},
    # The publication consumes the pipeline's volume mesh (PIPELINE_PUBLISHED_VOLUME);
    # stage_inputs(stage, pipeline) adds it.
    "canonical-gmsh-publication": {"source-process", "source-signature", "source-boundary"},
    "proper-rigid-publication": {
        "canonical-candidate-mesh", "canonical-build-record", "placement-transform",
        "source-semantic-contract", "source-signature", "source-boundary", "source-mask",
        "source-process"},
}


def stage_inputs(stage, pipeline):
    inputs = set(STAGE_INPUTS[stage])
    if stage == "canonical-gmsh-publication":
        inputs.add(PIPELINE_PUBLISHED_VOLUME[pipeline])
    return inputs


def stage_outputs(stage, pipeline):
    return set(STAGE_OUTPUTS[stage])
# Inputs a stage binds only when the case declares them. The device etch footprint
# (retained-etch.csv) is bound for seed generation exactly when the case declares
# a RetainedEtch source; otherwise the case records the producer default. The
# trace basis (basis contract, trace vertices/triangles, process library for its
# frame) is bound for seed generation and metric preparation exactly when the case
# declares all four roles; the seed sizes the frozen cut surface with it and the
# metric records the rule.
TRACE_BASIS_INPUTS = {"source-basis-contract": "BasisContract",
                      "source-trace-vertices": "TraceVertices",
                      "source-trace-triangles": "TraceTriangles",
                      "source-process-library": "ProcessLibrary"}
TRACE_BASIS_OPTIONS = {"--trace-basis-contract": "source-basis-contract",
                       "--trace-vertices": "source-trace-vertices",
                       "--trace-triangles": "source-trace-triangles",
                       "--process-library": "source-process-library"}
STAGE_OPTIONAL_INPUTS = {"seed-generation": {"source-retained-etch", *TRACE_BASIS_INPUTS},
                         "gmsh-build": {"source-retained-etch", *TRACE_BASIS_INPUTS},
                         "metric-preparation": set(TRACE_BASIS_INPUTS)}
STAGE_OUTPUTS = {
    "canonical-source-validation": {"canonical-semantic-contract", "canonical-supports"},
    "seed-generation": {"seed-mesh", "seed-corner-census"},
    # The Gmsh-only build: the labeled mixed-element volume mesh and its census
    # (the build report: tubes, per-type quality, size laws, footprint, junctions).
    "gmsh-build": {"gmsh-mesh", "build-census"},
    "metric-preparation": {"mmg-seed", "metric", "pins", "fixed-triangles",
                           "required-tetrahedra", "restoration-recipe"},
    "native-adaptation-mmg": {"adapted-mesh", "adaptation-receipt"},
    "label-restoration": {"source-local-restored-mesh", "restored-mesh"},
    "canonical-gmsh-publication": {"canonical-candidate-mesh",
                                   "canonical-ownership-partition",
                                   "canonical-ownership-quadrature-partition"},
    "proper-rigid-publication": {"candidate-mesh", "transformed-semantic-contract",
                                 "transformed-supports", "ownership-partition",
                                 "ownership-quadrature-partition", "transform-receipt"},
}
STAGE_PRIMARY_TOOL = {
    "canonical-source-validation": "source-validator",
    "seed-generation": "mesher",
    "gmsh-build": "mesher",
    "metric-preparation": "metric-preparer",
    "native-adaptation-mmg": "adaptation-wrapper",
    "label-restoration": "label-restorer",
    "canonical-gmsh-publication": "publisher",
    "proper-rigid-publication": "rigid-publisher",
}
# Named command options that must equal a bound tool, input, or output exactly.
STAGE_TOOL_OPTIONS = {
    "native-adaptation-mmg": {"--adapter": "adapter-mmg", "--mmg-library": "mmg-library"},
    "proper-rigid-publication": {"--ownership-runtime": "ownership-runtime",
                                 "--ownership-auditor": "ownership-auditor"},
}
STAGE_BINDING_OPTIONS = {
    "canonical-source-validation": {
        "--semantic-input": ("Inputs", "source-semantic-contract"),
        "--signature": ("Inputs", "source-signature"),
        "--boundary": ("Inputs", "source-boundary"),
        "--mask": ("Inputs", "source-mask")},
    "seed-generation": {"--mask": ("Inputs", "source-mask"),
                        "--boundary": ("Inputs", "source-boundary"),
                        "--semantic-contract": ("Inputs", "canonical-semantic-contract"),
                        "--corner-census": ("Artifacts", "seed-corner-census")},
    "gmsh-build": {"--mask": ("Inputs", "source-mask"),
                   "--boundary": ("Inputs", "source-boundary"),
                   "--semantic-contract": ("Inputs", "canonical-semantic-contract"),
                   "--corner-census": ("Artifacts", "build-census")},
    "metric-preparation": {
        "--semantic-contract": ("Inputs", "canonical-semantic-contract"),
        "--transformed-supports": ("Inputs", "canonical-supports"),
        "--seed-census": ("Inputs", "seed-corner-census")},
    "native-adaptation-mmg": {"--fixed-triangles": ("Inputs", "fixed-triangles"),
                              "--required-tetrahedra": ("Inputs", "required-tetrahedra")},
    "label-restoration": {
        "--source-local-output": ("Artifacts", "source-local-restored-mesh")},
    "canonical-gmsh-publication": {
        "--process": ("Inputs", "source-process"),
        "--signature": ("Inputs", "source-signature"),
        "--boundary": ("Inputs", "source-boundary")},
    "proper-rigid-publication": {
        "--semantic-input": ("Inputs", "source-semantic-contract"),
        "--signature": ("Inputs", "source-signature"),
        "--boundary": ("Inputs", "source-boundary"),
        "--mask": ("Inputs", "source-mask"),
        "--process": ("Inputs", "source-process"),
        "--canonical-build-record": ("Inputs", "canonical-build-record"),
        "--transformed-semantic": ("Artifacts", "transformed-semantic-contract"),
        "--transformed-supports": ("Artifacts", "transformed-supports"),
        "--ownership": ("Artifacts", "ownership-partition"),
        "--ownership-quadrature": ("Artifacts", "ownership-quadrature-partition")},
}
# Options bound to an optional input: required with the bound path when the input
# is bound, forbidden when it is not (an undeclared footprint is fail-closed).
STAGE_OPTIONAL_BINDING_OPTIONS = {
    "seed-generation": {"--etch-boundary": ("Inputs", "source-retained-etch"),
                        **{option: ("Inputs", name) for option, name in TRACE_BASIS_OPTIONS.items()}},
    "gmsh-build": {"--etch-boundary": ("Inputs", "source-retained-etch"),
                   **{option: ("Inputs", name) for option, name in TRACE_BASIS_OPTIONS.items()}},
    "metric-preparation": {option: ("Inputs", name) for option, name in TRACE_BASIS_OPTIONS.items()},
}
# Bound inputs/outputs that are consumed positionally and must occur in argv.
STAGE_BINDING_ARGUMENTS = {
    "canonical-source-validation": {
        ("Inputs", "canonical-transform"), ("Artifacts", "canonical-semantic-contract"),
        ("Artifacts", "canonical-supports")},
    "seed-generation": {("Inputs", "source-signature"), ("Artifacts", "seed-mesh")},
    "gmsh-build": {("Inputs", "source-signature"), ("Artifacts", "gmsh-mesh")},
    "metric-preparation": {("Inputs", "seed-mesh")},
    "native-adaptation-mmg": {
        ("Inputs", "mmg-seed"), ("Inputs", "metric"), ("Inputs", "pins"),
        ("Inputs", "restoration-recipe"), ("Artifacts", "adapted-mesh"),
        ("Artifacts", "adaptation-receipt")},
    "label-restoration": {("Inputs", "adapted-mesh"), ("Inputs", "restoration-recipe"),
                          ("Artifacts", "restored-mesh")},
    # The publication's consumed volume mesh is added per pipeline by
    # stage_binding_arguments.
    "canonical-gmsh-publication": {("Artifacts", "canonical-candidate-mesh")},
    "proper-rigid-publication": {
        ("Inputs", "canonical-candidate-mesh"), ("Inputs", "placement-transform"),
        ("Artifacts", "candidate-mesh"), ("Artifacts", "transform-receipt")},
}


def stage_binding_arguments(stage, pipeline):
    arguments = set(STAGE_BINDING_ARGUMENTS.get(stage, set()))
    if stage == "canonical-gmsh-publication":
        arguments.add(("Inputs", PIPELINE_PUBLISHED_VOLUME[pipeline]))
    return arguments


OWNERSHIP_AUDIT_OPTIONS = {"--process": "source-process", "--signature": "source-signature",
                           "--boundary": "source-boundary"}
# Seed options whose values must equal the metric recipe's corner-isotropy
# prescription, so the seed and the metric honor one ball around one corner set.
SEED_RECIPE_VALUE_OPTIONS = {"--lc-fine": "NormalSize",
                             "--corner-isotropy-radius": "CornerIsotropyRadius"}

_INTERPRETER_OPTIONS_WITH_VALUE = {
    "-H", "--home", "-J", "--sysimage", "-C", "--cpu-target", "-t", "--threads",
    "-p", "--procs", "--machine-file", "--gcthreads", "-W", "-X",
    "--check-hash-based-pycs",
}
_INTERPRETER_OPTIONS_WITH_ATTACHED_VALUE = (
    "--home=", "--sysimage=", "--cpu-target=", "--threads=", "--procs=",
    "--machine-file=", "--gcthreads=", "--project=", "--startup-file=",
    "--history-file=", "--compiled-modules=", "--pkgimages=", "--banner=",
    "--check-hash-based-pycs=", "-W", "-X",
)
_INTERPRETER_OPTIONS_WITHOUT_VALUE = {
    "-b", "-bb", "-B", "-d", "-h", "--help", "-i", "-I", "-O", "-OO", "-P",
    "-q", "--quiet", "-s", "-S", "-u", "-v", "-V", "--version", "-x", "--project",
}
_INTERPRETER_EXECUTION_OPTIONS = {"-c", "-m", "-e", "--eval", "--print", "-L", "--load"}


def _first_interpreter_script(command, primary):
    index = 1
    while index < len(command):
        argument = command[index]
        if argument == "--":
            return index + 1 if index + 1 < len(command) else None
        if (argument in _INTERPRETER_EXECUTION_OPTIONS or
                (argument == "-E" and Path(primary).suffix == ".jl")):
            return None
        if argument in _INTERPRETER_OPTIONS_WITH_VALUE:
            index += 2; continue
        if (argument in _INTERPRETER_OPTIONS_WITHOUT_VALUE or
                (argument == "-E" and Path(primary).suffix == ".py") or
                argument.startswith(_INTERPRETER_OPTIONS_WITH_ATTACHED_VALUE)):
            index += 1; continue
        return index
    return None


def _resolved_argument(argument, base):
    path = Path(argument)
    return (path if path.is_absolute() else base / path).resolve()


def _require_path_option(command, option, expected, working_directory=None):
    positions = [index for index, value in enumerate(command) if value == option]
    if len(positions) != 1 or positions[0] + 1 >= len(command):
        raise ValueError(f"Stage command must provide exactly one {option}")
    base = Path(working_directory or Path.cwd())
    if _resolved_argument(command[positions[0] + 1], base) != Path(expected).resolve():
        raise ValueError(f"Stage command {option} differs from its bound input or tool")


def _require_path_argument(command, expected, description, working_directory=None):
    base = Path(working_directory or Path.cwd())
    if str(Path(expected).resolve()) not in [str(_resolved_argument(value, base))
                                             for value in command]:
        raise ValueError(f"Stage command did not consume its bound {description}")


def _require_runtime_and_script(command, runtime, script, description, working_directory=None):
    base = Path(working_directory or Path.cwd())
    runtime, script = str(Path(runtime).resolve()), str(Path(script).resolve())
    def matches(argument, expected):
        return str(_resolved_argument(argument, base)) == expected
    if not isinstance(command, list) or not command or not matches(command[0], runtime):
        raise ValueError(f"{description} runtime is not argv[0]")
    script_position = _first_interpreter_script(command, script)
    if script_position is None or not matches(command[script_position], script):
        raise ValueError(f"{description} tool is not in the executed script position")


def validate_tool_invocation(stage, command, tools, working_directory=None):
    try:
        _require_runtime_and_script(command, tools["runtime"], tools[STAGE_PRIMARY_TOOL[stage]],
                                    "Stage", working_directory)
    except ValueError as error:
        raise ValueError(f"{error}: {stage}") from error
    for option, role in STAGE_TOOL_OPTIONS.get(stage, {}).items():
        _require_path_option(command, option, tools[role], working_directory)


def validate_command_bindings(stage, command, inputs, artifacts, working_directory=None,
                              pipeline=LEGACY_MMG_PIPELINE):
    """Require every consumed source/input/output path to be the bound one."""
    bound = {"Inputs": inputs, "Artifacts": artifacts}
    for option, (section, name) in STAGE_BINDING_OPTIONS.get(stage, {}).items():
        _require_path_option(command, option, bound[section][name], working_directory)
    for option, (section, name) in STAGE_OPTIONAL_BINDING_OPTIONS.get(stage, {}).items():
        if name in bound[section]:
            _require_path_option(command, option, bound[section][name], working_directory)
        elif option in command:
            raise ValueError(f"Stage command passes {option} without a bound {name} input")
    for section, name in stage_binding_arguments(stage, pipeline):
        _require_path_argument(command, bound[section][name], f"{section.lower()} {name}",
                               working_directory)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def binding(path):
    path = Path(path).resolve()
    if not path.is_file():
        raise ValueError(f"Stage artifact does not exist: {path}")
    return {"Path": str(path), "SHA256": sha256(path)}


def _validate_bindings(items, expected, description, optional=frozenset()):
    if (not isinstance(items, dict) or not expected <= set(items) or
            not set(items) <= expected | set(optional)):
        raise ValueError(f"{description} names do not match the frozen stage contract")
    for name, item in items.items():
        if (not isinstance(item, dict) or not item.get("Path") or not item.get("SHA256") or
                not Path(item["Path"]).is_file() or sha256(item["Path"]) != item["SHA256"]):
            raise ValueError(f"{description} binding changed: {name}")


def validate_stage_report(report, stage, launcher_name=None, launcher_sha256=None,
                          expected_tool_sha256=None, pipeline=None):
    """Validate one bounded stage report.  `pipeline` names the canonical pipeline
    the stage belongs to; when None the publication stage is judged by the volume
    input it bound (exactly one pipeline's), every other stage is pipeline-free."""
    if pipeline is None:
        pipeline = _report_pipeline(report, stage)
    if (report.get("Version") != 3 or report.get("Stage") != stage or
            report.get("ReturnCode") != 0 or report.get("StopReason") is not None or
            not isinstance(report.get("Command"), list) or not report["Command"] or
            not isinstance(report.get("Environment"), dict) or
            not isinstance(report.get("WorkingDirectory"), str) or
            not Path(report["WorkingDirectory"]).is_absolute()):
        raise ValueError(f"Invalid or unsuccessful bounded stage: {stage}")
    producer = report.get("Producer", {})
    if launcher_name is not None and (producer.get("Name") != launcher_name or
                                      producer.get("SHA256") != launcher_sha256):
        raise ValueError("Stage launcher identity differs from the frozen launcher")
    _validate_bindings(report.get("Inputs"), stage_inputs(stage, pipeline), "stage input",
                       STAGE_OPTIONAL_INPUTS.get(stage, frozenset()))
    _validate_bindings(report.get("Artifacts"), stage_outputs(stage, pipeline), "stage output")
    tools = report.get("Tools")
    if not isinstance(tools, dict) or set(tools) != STAGE_TOOLS[stage]:
        raise ValueError(f"Stage tools are incomplete: {stage}")
    for role, item in tools.items():
        if (not isinstance(item, dict) or not item.get("Path") or not item.get("SHA256") or
                not Path(item["Path"]).is_file() or sha256(item["Path"]) != item["SHA256"] or
                (expected_tool_sha256 is not None and
                 item["SHA256"] != expected_tool_sha256.get(role))):
            raise ValueError(f"Stage tool binding changed: {stage}/{role}")
    tool_paths = {role: item["Path"] for role, item in tools.items()}
    validate_tool_invocation(stage, report["Command"], tool_paths,
                             report.get("WorkingDirectory"))
    validate_command_bindings(
        stage, report["Command"],
        {name: item["Path"] for name, item in report["Inputs"].items()},
        {name: item["Path"] for name, item in report["Artifacts"].items()},
        report.get("WorkingDirectory"), pipeline)
    if stage == "canonical-source-validation":
        transform = json.loads(Path(report["Inputs"]["canonical-transform"]["Path"]).read_text())
        if isinstance(transform, dict): transform = transform.get("Transform")
        if transform != [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]:
            raise ValueError("Canonical source validation must use the identity transform")
    elif stage == "native-adaptation-mmg":
        recipe = report["Inputs"]["restoration-recipe"]["Path"]
        receipt = json.loads(Path(report["Artifacts"]["adaptation-receipt"]["Path"]).read_text())
        recipe_data = json.loads(Path(recipe).read_text())
        hmax = recipe_data.get("FarFieldBudgetPolicy", {}).get("EffectiveFarSize")
        if (receipt.get("RecipeSHA256") != sha256(recipe) or
                receipt.get("EffectiveFarSize") != hmax or
                float(receipt.get("HmaxArgument", "nan")) != hmax or
                receipt.get("OutputSHA256") != report["Artifacts"]["adapted-mesh"]["SHA256"]):
            raise ValueError("Native MMG hmax receipt differs from the bound metric policy")
        required = report["Inputs"]["required-tetrahedra"]
        if (receipt.get("RequiredTetrahedraSHA256") != required["SHA256"] or
                receipt.get("RequiredTetrahedra") !=
                recipe_data.get("RequiredTetrahedra", {}).get("Count")):
            raise ValueError("Native MMG receipt required tetrahedra differ from the bound list")
        if (receipt.get("AdapterSHA256") != tools["adapter-mmg"]["SHA256"] or
                receipt.get("MMGLibrarySHA256") != tools["mmg-library"]["SHA256"] or
                str(Path(receipt.get("MMGLibraryPath", "")).resolve()) !=
                str(Path(tools["mmg-library"]["Path"]).resolve())):
            raise ValueError("Native MMG receipt adapter/library differ from the bound stage tools")
    elif stage == "proper-rigid-publication":
        _validate_ownership_audit_binding(report)
    return report


def stage_pipeline_of_inputs(stage, inputs):
    """The pipeline a stage belongs to when the stage is pipeline-specific
    (gmsh-build, the legacy stages, or the publication by its bound volume input
    name); the legacy pipeline for the shared stages (their contract is identical)."""
    if stage == "gmsh-build":
        return GMSH_ONLY_PIPELINE
    if stage == "canonical-gmsh-publication":
        bound = [pipeline for pipeline, name in PIPELINE_PUBLISHED_VOLUME.items() if name in inputs]
        if len(bound) != 1:
            raise ValueError("Publication stage must bind exactly one pipeline's volume mesh")
        return bound[0]
    return LEGACY_MMG_PIPELINE


def _report_pipeline(report, stage):
    inputs = report.get("Inputs") if isinstance(report.get("Inputs"), dict) else {}
    return stage_pipeline_of_inputs(stage, inputs)


def _validate_ownership_audit_binding(report):
    """The transformed-ownership audit must consume only bound source and outputs."""
    inputs, artifacts, tools = report["Inputs"], report["Artifacts"], report["Tools"]
    receipt = json.loads(Path(artifacts["transform-receipt"]["Path"]).read_text())
    command = receipt.get("OwnershipCommand")
    working_directory = report.get("WorkingDirectory")
    try:
        _require_runtime_and_script(command, tools["ownership-runtime"]["Path"],
                                    tools["ownership-auditor"]["Path"], "Ownership audit",
                                    working_directory)
        for option, name in OWNERSHIP_AUDIT_OPTIONS.items():
            _require_path_option(command, option, inputs[name]["Path"], working_directory)
        for name in ("candidate-mesh", "ownership-partition"):
            _require_path_argument(command, artifacts[name]["Path"], f"artifacts {name}",
                                   working_directory)
    except ValueError as error:
        raise ValueError(f"Ownership audit command is not bound: {error}") from error
    expected_hashes = {
        "OwnershipRuntimeSHA256": tools["ownership-runtime"]["SHA256"],
        "OwnershipAuditorSHA256": tools["ownership-auditor"]["SHA256"],
        "OwnershipSHA256": artifacts["ownership-partition"]["SHA256"],
        "OwnershipQuadratureSHA256": artifacts["ownership-quadrature-partition"]["SHA256"],
        "TransformedSemanticSHA256": artifacts["transformed-semantic-contract"]["SHA256"],
        "TransformedSupportsSHA256": artifacts["transformed-supports"]["SHA256"],
        "SourceInputSHA256": {
            name: inputs[name]["SHA256"]
            for name in ("source-semantic-contract", "source-signature", "source-boundary",
                         "source-mask", "source-process")},
    }
    if any(receipt.get(key) != value for key, value in expected_hashes.items()):
        raise ValueError("Rigid publication receipt source/ownership hashes differ from bindings")


def _validate_reports(report_paths, order, launcher_name, launcher_sha256,
                      expected_tool_sha256, pipeline):
    if set(report_paths) != set(order):
        raise ValueError("Stage report set differs from the required DAG role")
    reports, report_digests = {}, set()
    for stage in order:
        path = Path(report_paths[stage])
        expected = None if expected_tool_sha256 is None else expected_tool_sha256[stage]
        report = validate_stage_report(json.loads(path.read_text()), stage, launcher_name,
                                       launcher_sha256, expected, pipeline)
        digest = sha256(path)
        if digest in report_digests:
            raise ValueError("Bounded stage reports must be content-distinct")
        report_digests.add(digest); reports[stage] = report
    return reports, report_digests


def _option_value(command, option):
    positions = [index for index, value in enumerate(command) if value == option]
    if len(positions) != 1 or positions[0] + 1 >= len(command):
        raise ValueError(f"Stage command must provide exactly one {option}")
    return command[positions[0] + 1]


def _recipe_number(recipe, name):
    value = recipe.get(name)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"Restoration recipe lacks a numeric {name}")
    return float(value)


def _count(value, description):
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{description} must be a non-negative integer count")
    return value


def validate_protected_corner_balls(recipe):
    """The metric stage must freeze the seed's corner balls: the recipe records the
    rule's radius (the recipe's CornerIsotropyRadius) and the frozen-triangle count
    of every semantic corner. Counts are reported, not gated."""
    balls = recipe.get("ProtectedCornerBalls")
    corners = recipe.get("TruePhysicalCorners")
    if (not isinstance(balls, dict) or not isinstance(corners, list) or
            (not corners and not contract_is_corner_free(recipe.get("SemanticContract"))) or
            not isinstance(balls.get("PerCorner"), list)):
        raise ValueError("Restoration recipe lacks the protected corner balls")
    if _recipe_number(balls, "Radius") != _recipe_number(recipe, "CornerIsotropyRadius"):
        raise ValueError("Protected corner ball radius differs from the recipe CornerIsotropyRadius")
    frozen = _count(balls.get("FrozenTriangles"), "Protected corner ball frozen triangles")
    if frozen > _count(recipe.get("FixedSurfaceTriangles"), "Fixed surface triangles"):
        raise ValueError("Protected corner balls exceed the fixed surface triangles")
    per_corner = balls["PerCorner"]
    if (len(per_corner) != len(corners) or
            any(not isinstance(row, dict) for row in per_corner) or
            sorted(row.get("Point") for row in per_corner) != sorted(corners)):
        raise ValueError("Protected corner balls differ from the recipe semantic corners")
    if sum(_count(row.get("FrozenTriangles"), "Protected corner ball frozen triangles")
           for row in per_corner) < frozen:
        raise ValueError("Protected corner ball counts do not cover the frozen triangles")
    return balls


PRODUCER_DEFAULT_FOOTPRINT_PROVENANCE = "producer-default"


def footprint_provenance(census):
    """Provenance of the census footprint: the bound retained-etch SHA-256 or the
    producer default."""
    digest = census.get("EtchBoundarySHA256")
    if digest is None:
        if census.get("EtchBoundary") != PRODUCER_DEFAULT_FOOTPRINT_PROVENANCE:
            raise ValueError("Seed census footprint provenance is neither a bound file nor the default")
        return PRODUCER_DEFAULT_FOOTPRINT_PROVENANCE
    if (not isinstance(digest, str) or len(digest) != 64 or
            any(character not in "0123456789abcdef" for character in digest)):
        raise ValueError("Seed census footprint hash is invalid")
    return digest


def census_coupon_kind(census):
    """The coupon kind a Gmsh-only build census records (PrismTubes.Section.Kind, decision
    66); fabricated for a census without the record (the legacy seed census)."""
    section = (census.get("PrismTubes") or {}).get("Section") if isinstance(census, dict) else None
    kind = section.get("Kind", "fabricated") if isinstance(section, dict) else "fabricated"
    if kind not in COUPON_KINDS:
        raise ValueError(f"the census records an unknown coupon kind {kind!r}")
    return kind


def footprint_segments(census):
    """Every edge of every simplified footprint polygon as a 3D segment
    [x0, y0, z, x1, y1, z] on the polygon's process plane, in census order (none for a
    thin coupon)."""
    segments = []
    for polygon in validate_footprint_polygons(census, census_coupon_kind(census)):
        points, plane = polygon["Points"], float(polygon["Plane"])
        for index, point in enumerate(points):
            following = points[(index + 1) % len(points)]
            segments.append([float(point[0]), float(point[1]), plane,
                             float(following[0]), float(following[1]), plane])
    return segments


def validate_footprint_segments(recipe, census):
    """The metric recipe's FootprintSegments must be exactly the bound census's
    simplified footprint edges with their provenance and tolerance."""
    record = recipe.get("FootprintSegments")
    if (not isinstance(record, dict) or
            record.get("Provenance") != footprint_provenance(census) or
            record.get("EtchBoundary") != census.get("EtchBoundary") or
            record.get("Tolerance") != census.get("FootprintCollinearTolerance") or
            record.get("Segments") != footprint_segments(census) or
            record.get("Polygons") != len(census["FootprintPolygons"])):
        raise ValueError("Restoration recipe FootprintSegments differ from the bound seed census")
    return record


def validate_junction_segments(recipe, census):
    """The metric recipe's JunctionSegments are the cut-surface/material-interface
    junction lines of the seed: finite nondegenerate segments whose count and total
    length are self-consistent, whose label sets are the recipe contract's cut and
    two-material interface labels, and whose total length equals the census's CAD
    junction curves within the shared dimensionless tolerance."""
    record = recipe.get("JunctionSegments")
    semantic = recipe.get("SemanticContract")
    if (not isinstance(record, dict) or not isinstance(record.get("Segments"), list) or
            not record["Segments"] or not isinstance(semantic, dict)):
        raise ValueError("Restoration recipe lacks the junction segments")
    total = 0.0
    for segment in record["Segments"]:
        if (not isinstance(segment, list) or len(segment) != 6 or
                any(isinstance(value, bool) or not isinstance(value, (int, float)) or
                    not math.isfinite(value) for value in segment)):
            raise ValueError("Restoration recipe junction segment is invalid")
        length = math.dist(segment[:3], segment[3:])
        if length <= 0:
            raise ValueError("Restoration recipe junction segment is degenerate")
        total += length
    if (record.get("Count") != len(record["Segments"]) or
            abs(_recipe_number(record, "TotalLength") - total) > 1e-12 * total or
            record.get("CutSurfaceAttributes") != sorted(cut_surface_attributes(semantic)) or
            record.get("MaterialInterfaceAttributes") !=
            sorted(material_interface_attributes(semantic))):
        raise ValueError("Restoration recipe junction segments are inconsistent with their contract")
    curves = census.get("JunctionCurves")
    if (not isinstance(curves, dict) or
            abs(_recipe_number(curves, "TotalLength") - total) > COPLANAR_TOLERANCE * total or
            _count(curves.get("Count"), "Seed census junction curves") <= 0):
        raise ValueError("Seed census junction curves differ from the recipe junction segments")
    return record


def validate_footprint_polygons(census, kind="fabricated"):
    """The seed simplified every etch footprint polygon (device or producer default)
    with the shared collinearity tolerance before CAD face creation and recorded the
    result: removed vertices and a maximum deviation within tolerance x local scale.
    A thin coupon (decision 66) etches nothing: it records no polygon."""
    tolerance = census.get("FootprintCollinearTolerance")
    if isinstance(tolerance, bool) or tolerance != COPLANAR_TOLERANCE:
        raise ValueError("Seed census footprint collinearity tolerance differs from COPLANAR_TOLERANCE")
    polygons = census.get("FootprintPolygons")
    summary = census.get("FootprintSimplification")
    if not isinstance(polygons, list) or (bool(polygons) != (kind == "fabricated")) or not isinstance(summary, dict):
        raise ValueError("Seed census lacks the simplified footprint polygons" if kind == "fabricated" else
                         "Thin coupon census records etch footprint polygons")
    removed, worst = 0, 0.0
    for polygon in polygons:
        record = polygon.get("Simplification") if isinstance(polygon, dict) else None
        points = polygon.get("Points") if isinstance(polygon, dict) else None
        if (not isinstance(record, dict) or not isinstance(points, list) or len(points) < 3 or
                any(not isinstance(point, list) or len(point) != 2 or
                    any(isinstance(value, bool) or not isinstance(value, (int, float))
                        for value in point) for point in points) or
                not isinstance(polygon.get("Hole"), bool) or
                isinstance(polygon.get("Plane"), bool) or
                not isinstance(polygon.get("Plane"), (int, float)) or
                _count(polygon.get("Conductor"), "Footprint conductor") <= 0):
            raise ValueError("Seed census footprint polygon is incomplete")
        indices = record.get("RemovedVertexIndices")
        original = _count(record.get("OriginalVertices"), "Footprint original vertices")
        count = _count(record.get("RemovedVertexCount"), "Footprint removed vertices")
        if (not isinstance(indices, list) or len(indices) != count or
                any(_count(index, "Footprint removed vertex index") <= 0 or index > original
                    for index in indices) or len(set(indices)) != count or
                _count(record.get("Vertices"), "Footprint vertices") != len(points) or
                original - count != len(points) or record.get("Tolerance") != tolerance):
            raise ValueError("Seed census footprint simplification record is inconsistent")
        deviation = _recipe_number(record, "MaximumDeviation")
        scale = _recipe_number(record, "MaximumDeviationLocalScale")
        relative = _recipe_number(record, "MaximumRelativeDeviation")
        # The bound is deviation <= tolerance x local scale; the relative value is
        # the recorded quotient (1e-12 covers its floating-point roundoff only).
        if (deviation < 0 or scale < 0 or relative < 0 or relative > tolerance or
                deviation > tolerance * scale * (1 + 1e-12)):
            raise ValueError("Seed census footprint simplification exceeds the collinearity tolerance")
        removed += count; worst = max(worst, relative)
    if (summary.get("Polygons") != len(polygons) or summary.get("RemovedVertices") != removed or
            _recipe_number(summary, "MaximumRelativeDeviation") != worst):
        raise ValueError("Seed census footprint simplification summary differs from its polygons")
    return polygons


def validate_seed_corner_isotropy(seed_report, recipe_path):
    """The seed's corner ball must be the recipe's: same size, radius and corners,
    and the metric stage must have protected the balls it received."""
    recipe = json.loads(Path(recipe_path).read_text())
    validate_protected_corner_balls(recipe)
    command = seed_report["Command"]
    for option, name in SEED_RECIPE_VALUE_OPTIONS.items():
        try:
            value = float(_option_value(command, option))
        except ValueError as error:
            raise ValueError(f"Seed command {option} is not a bound number: {error}") from error
        if value != _recipe_number(recipe, name):
            raise ValueError(f"Seed command {option} differs from the recipe {name}")
    census = json.loads(Path(seed_report["Artifacts"]["seed-corner-census"]["Path"]).read_text())
    if (not isinstance(census, dict) or census.get("Version") != 1 or
            not isinstance(census.get("Corners"), list) or
            not isinstance(census.get("LongitudinalFaces"), list)):
        raise ValueError("Seed corner census has an unsupported schema")
    # The etched footprint the seed used must be the bound one (or the recorded
    # producer default) and its per-label interface areas are recorded.
    etch = seed_report["Inputs"].get("source-retained-etch")
    if etch is None:
        if census.get("EtchBoundary") != "producer-default":
            raise ValueError("Seed census names an etch footprint the stage did not bind")
    elif census.get("EtchBoundarySHA256") != etch["SHA256"]:
        raise ValueError("Seed census etch footprint differs from the bound retained etch")
    areas = census.get("InterfaceAreas")
    if (not isinstance(areas, list) or not areas or
            any(not isinstance(row, dict) or isinstance(row.get("Attribute"), bool) or
                not isinstance(row.get("Attribute"), int) or
                not isinstance(row.get("Area"), (int, float)) or
                isinstance(row.get("Area"), bool) or not row["Area"] > 0
                for row in areas)):
        raise ValueError("Seed census lacks positive per-label interface areas")
    # The census must measure the labels the seed is written with (after any
    # slot/conductor relabeling), which are exactly the contract's boundary labels.
    semantic = recipe.get("SemanticContract")
    if not isinstance(semantic, dict) or not isinstance(semantic.get("BoundaryLabels"), list):
        raise ValueError("Restoration recipe lacks the semantic contract boundary labels")
    labels = [row["Attribute"] for row in areas]
    if len(labels) != len(set(labels)) or set(labels) != boundary_attributes(semantic):
        raise ValueError("Seed census interface-area labels differ from the semantic contract")
    validate_footprint_polygons(census)
    validate_footprint_segments(recipe, census)
    validate_junction_segments(recipe, census)
    # The ridge-to-ridge face census (interior nodes and full-height triangles) is
    # recorded for every seed; its values are reported, not gated.
    for row in census["LongitudinalFaces"]:
        if (not isinstance(row, dict) or
                any(not isinstance(row.get(key), int) or isinstance(row.get(key), bool)
                    for key in ("Surface", "Triangles", "InteriorNodes",
                                "InteriorNodesAwayFromCorners", "FullHeightTriangles",
                                "FullHeightTrianglesAwayFromCorners"))):
            raise ValueError("Seed corner census longitudinal face rows are incomplete")
    if (_recipe_number(census, "IsotropicSize") != _recipe_number(recipe, "NormalSize") or
            _recipe_number(census, "CornerIsotropyRadius") !=
            _recipe_number(recipe, "CornerIsotropyRadius")):
        raise ValueError("Seed corner census size or radius differs from the recipe")
    corners = recipe.get("TruePhysicalCorners")
    if (not isinstance(corners, list) or (not corners and not contract_is_corner_free(semantic)) or
            sorted(census.get("SemanticCorners", [])) != sorted(corners) or
            len(census["Corners"]) != len(corners)):
        raise ValueError("Seed corner census corners differ from the recipe semantic corners")
    return census


# Edge layer options of the seed and metric commands bound to the recipe's
# EdgeLayer record (seed and metric prescribe one layer: EdgeSize, GrowthRatio,
# EdgeLayerAspect); the adapter hmin must be the layer's EdgeSize.
EDGE_LAYER_SEED_OPTIONS = {"--edge-size": "EdgeSize", "--edge-growth-ratio": "GrowthRatio",
                           "--edge-layer-aspect": "Aspect"}
EDGE_LAYER_METRIC_OPTIONS = {"--edge-size": "EdgeSize", "--edge-growth-ratio": "GrowthRatio",
                             "--edge-layer-aspect": "Aspect"}
EDGE_LAYER_DEFAULTS = {"GrowthRatio": 2.0, "Aspect": 4.0}
ADAPTATION_MINIMUM_SIZE_OPTION = "--hmin"


def _option_or_default(command, option, default):
    positions = [index for index, value in enumerate(command) if value == option]
    if not positions:
        if default is None:
            raise ValueError(f"Stage command must provide {option}")
        return default
    try:
        return float(_option_value(command, option))
    except ValueError as error:
        raise ValueError(f"Stage command {option} is not a bound number: {error}") from error


def validate_edge_layer(seed_report, metric_report, adaptation_report, recipe, census):
    """The seed, the metric stage and the adapter honor one edge layer or none.

    With a recipe EdgeLayer record the seed and metric commands pass EdgeSize,
    GrowthRatio and EdgeLayerAspect equal to the record (ratio/aspect may be the
    documented defaults), the census records the same layer (sizes, rows, span
    length, row nodes), the recipe reach and thickness follow from the sizes,
    the layer lies inside the frozen band, and the adapter hmin is EdgeSize.
    Without the record, neither command asks for a layer and the census records
    none; the adapter hmin is then the recipe NormalSize.
    """
    layer = recipe.get("EdgeLayer")
    seed_command = seed_report["Command"]
    metric_command = metric_report["Command"]
    census_layer = census.get("EdgeLayer") if isinstance(census, dict) else None
    normal = _recipe_number(recipe, "NormalSize")
    if layer is None:
        if _option_or_default(seed_command, "--edge-size", 0.0) != 0.0:
            raise ValueError("Seed command seeds an edge layer the recipe does not record")
        if any(option in metric_command for option in EDGE_LAYER_METRIC_OPTIONS):
            raise ValueError("Metric command binds an edge layer the recipe does not record")
        if isinstance(census_layer, dict) and census_layer.get("EdgeSize", 0.0) != 0.0:
            raise ValueError("Seed census records an edge layer the recipe does not record")
        minimum = normal
    else:
        if not isinstance(layer, dict) or not isinstance(census_layer, dict):
            raise ValueError("Edge layer recipe or census record is missing")
        edge_size = _recipe_number(layer, "EdgeSize")
        ratio = _recipe_number(layer, "GrowthRatio")
        aspect = _recipe_number(layer, "Aspect")
        if not 0 < edge_size < normal or ratio <= 1 or aspect < 1:
            raise ValueError("Edge layer sizes are invalid")
        for command, options, stage in ((seed_command, EDGE_LAYER_SEED_OPTIONS, "Seed"),
                                        (metric_command, EDGE_LAYER_METRIC_OPTIONS, "Metric")):
            for option, name in options.items():
                value = _option_or_default(command, option, EDGE_LAYER_DEFAULTS.get(name))
                if value != _recipe_number(layer, name):
                    raise ValueError(f"{stage} command {option} differs from the recipe edge layer {name}")
        for name in ("EdgeSize", "GrowthRatio", "Aspect", "Layers", "LayerThickness",
                     "TotalSpanLength", "Rows"):
            if census_layer.get(name) != layer.get(name if name != "Rows" else "SeedRows"):
                raise ValueError(f"Seed census edge layer {name} differs from the recipe")
        if (census_layer.get("RowOffsets") != layer.get("RowOffsets") or
                census_layer.get("RowNodes") != layer.get("SeedRowNodes") or
                census_layer.get("NormalSize") != normal):
            raise ValueError("Seed census edge layer rows differ from the recipe")
        offsets = layer["RowOffsets"]
        expected = []
        size = edge_size
        while size < normal:
            expected.append(expected[-1] + size if expected else size)
            size *= ratio
        if (len(offsets) != len(expected) or
                any(abs(a - b) > 1e-12 * normal for a, b in zip(offsets, expected)) or
                _recipe_number(layer, "Reach") != (normal - edge_size) / (ratio - 1.0) or
                _recipe_number(layer, "LayerThickness") != offsets[-1]):
            raise ValueError("Edge layer rows do not follow EdgeSize and GrowthRatio")
        if offsets[-1] > _recipe_number(recipe, "SurfaceProtectionRadius"):
            raise ValueError("Edge layer is thicker than the frozen band")
        spans = layer.get("Spans")
        if (not isinstance(spans, list) or len(spans) != len(census_layer.get("Curves", [])) or
                layer.get("SpanCount") != len(spans)):
            raise ValueError("Edge layer spans differ from the seed census curves")
        minimum = edge_size
    if _option_or_default(adaptation_report["Command"], ADAPTATION_MINIMUM_SIZE_OPTION, None) != minimum:
        raise ValueError("Adaptation command --hmin differs from the recipe minimum size")
    return layer


# Corner grading (supervisor decision 33): the seed and the metric stage prescribe
# one CornerSize inside the corner balls or none (0 / absent = uniform NormalSize).
CORNER_SIZE_OPTION = "--corner-size"


def validate_corner_grading(seed_report, metric_report, recipe, census):
    """The seed, the metric stage, the census and the recipe agree on the corner
    grading.  With a recipe CornerGrading record: the seed and metric commands pass
    --corner-size equal to its CornerSize (0 < CornerSize < NormalSize), its
    GrowthRatio is the edge-layer growth ratio of both commands (or the default),
    its radius is the recipe CornerIsotropyRadius, its Reach = (NormalSize -
    CornerSize) / (GrowthRatio - 1) lies inside the radius, its ShellRadii are the
    cumulative geometric offsets below NormalSize followed by the radius, and the
    census CornerGrading records the same values; with a recorded edge layer, the
    layer reaches the ball (census LayerReachesCornerBall, no taper) and records
    UnlayeredEdgeLengthPerCorner.  Without the record neither command passes a
    corner size and the census records none.  Returns the record or None."""
    grading = recipe.get("CornerGrading")
    seed_value = _option_or_default(seed_report["Command"], CORNER_SIZE_OPTION, 0.0)
    metric_value = _option_or_default(metric_report["Command"], CORNER_SIZE_OPTION, 0.0)
    census_grading = census.get("CornerGrading") if isinstance(census, dict) else None
    census_size = (float(census_grading.get("CornerSize", 0.0))
                   if isinstance(census_grading, dict) else 0.0)
    if grading is None:
        if seed_value != 0.0 or metric_value != 0.0 or census_size != 0.0:
            raise ValueError("Corner grading requested without a recipe CornerGrading record")
        return None
    normal = _recipe_number(recipe, "NormalSize")
    radius = _recipe_number(recipe, "CornerIsotropyRadius")
    corner_size = _recipe_number(grading, "CornerSize")
    ratio = _recipe_number(grading, "GrowthRatio")
    if not 0.0 < corner_size < normal or ratio <= 1.0:
        raise ValueError("Corner grading sizes are invalid")
    if seed_value != corner_size or metric_value != corner_size or census_size != corner_size:
        raise ValueError("Seed command, metric command and census differ from the recipe CornerSize")
    for command, stage in ((seed_report["Command"], "Seed"), (metric_report["Command"], "Metric")):
        if _option_or_default(command, "--edge-growth-ratio", EDGE_LAYER_DEFAULTS["GrowthRatio"]) != ratio:
            raise ValueError(f"{stage} command --edge-growth-ratio differs from the corner grading ratio")
    expected = []
    size = corner_size
    while size < normal:
        expected.append(expected[-1] + size if expected else size)
        size *= ratio
    reach = (normal - corner_size) / (ratio - 1.0)
    sizes = [corner_size * ratio**k for k in range(len(expected))] + [normal]
    if (_recipe_number(grading, "NormalSize") != normal or _recipe_number(grading, "Radius") != radius or
            _recipe_number(grading, "Reach") != reach or not reach <= radius or
            grading.get("ShellRadii") != expected + [radius] or grading.get("ShellSizes") != sizes):
        raise ValueError("Corner grading shells do not follow CornerSize, GrowthRatio and the radius")
    if (not isinstance(census_grading, dict) or
            any(census_grading.get(name) != grading.get(name)
                for name in ("CornerSize", "GrowthRatio", "NormalSize", "Radius", "Reach", "ShellRadii",
                             "ShellSizes"))):
        raise ValueError("Seed census corner grading differs from the recipe")
    layer = census.get("EdgeLayer") if isinstance(census, dict) else None
    if isinstance(layer, dict) and layer.get("EdgeSize", 0.0) != 0.0:
        rows = layer.get("UnlayeredEdgeLengthPerCorner")
        if (layer.get("LayerReachesCornerBall") is not True or layer.get("TaperSubdivisions") != [] or
                not isinstance(rows, list) or len(rows) != len(recipe.get("TruePhysicalCorners", [])) or
                any(not isinstance(row, dict) or not isinstance(row.get("UnlayeredLengths"), list)
                    for row in rows)):
            raise ValueError("Seed census edge layer does not record the layer reaching the corner balls")
    return grading


# Quality gates the seed stage enforces on the MMG required region (decision 30;
# the Jacobian condition bound, the manifest MaximumJacobianCondition, decision 34)
# and the label restorer enforces on the adapted mesh: both commands carry them,
# with equal values.
REQUIRED_REGION_GATE_OPTIONS = ("--maximum-corner-aspect", "--minimum-scaled-jacobian",
                                "--maximum-jacobian-condition",
                                "--maximum-quality-displacement-over-normal")
# The edge-layer quality rule (supervisor decision 32, calibration manifests only):
# both commands carry the same bound, or neither; it needs a recipe edge layer.
EDGE_LAYER_QUALITY_RULE_OPTION = "--edge-layer-maximum-aspect"


def _required_indices(path, tetrahedra):
    values = Path(path).read_text().split()
    if not values or any(not value.isdigit() for value in values):
        raise ValueError("Required tetrahedron list must be non-empty positive integers")
    indices = [int(value) for value in values]
    if (indices != sorted(set(indices)) or indices[0] < 1 or indices[-1] > tetrahedra):
        raise ValueError("Required tetrahedron list must be sorted, unique and in range")
    return indices


def _restoration_quality_report(restoration_report):
    """The label restorer's own report (next to its restored-mesh artifact)."""
    restored = Path(restoration_report["Artifacts"]["restored-mesh"]["Path"])
    report = json.loads(restored.with_suffix(".projection.json").read_text())
    if not isinstance(report, dict):
        raise ValueError("Label restoration report is not a record")
    return report


def validate_required_region(seed_report, restoration_report, recipe, census, required_path):
    """The metric stage's required tetrahedra are the recipe's corner balls and edge
    layer, the list is the recorded one, and the seed stage optimized and gated the
    same region with the restorer's gate values.

    The recipe RequiredTetrahedra record names the rule, CornerIsotropyRadius, one
    count per semantic corner (each positive), the layer reach LayerRequiredReach =
    LayerThickness x (1 + RowZigzag) + EdgeSize (also the recipe EdgeLayer
    RequiredReach) and one count per span when a layer is recorded (no reach and no
    spans otherwise); the list has exactly Count sorted unique seed indices.  The
    seed command passes the four gate options with the restorer's values; the
    census SeedQualityOptimization record carries them, gated exactly Count cells
    (the set recomputed on the final seed positions), no required cell below the
    scaled-Jacobian gate or above the Jacobian condition gate and every corner
    aspect within the corner gate; the label restorer found exactly Count required
    tetrahedra in the adapted mesh (MMG kept them all) and reports their maximum
    Jacobian condition within the gate.
    """
    record = recipe.get("RequiredTetrahedra")
    corners = recipe.get("TruePhysicalCorners")
    if (not isinstance(record, dict) or not isinstance(corners, list) or
            (not corners and not contract_is_corner_free(recipe.get("SemanticContract"))) or
            not isinstance(record.get("PerCorner"), list) or
            not isinstance(record.get("PerSpan"), list) or record.get("IndexBase") != 1):
        raise ValueError("Restoration recipe lacks the required-tetrahedra record")
    tetrahedra = _count(recipe.get("Tetrahedra"), "Recipe tetrahedra")
    indices = _required_indices(required_path, tetrahedra)
    if _count(record.get("Count"), "Required tetrahedra") != len(indices) or not indices:
        raise ValueError("Required tetrahedron list differs from the recipe record")
    if _recipe_number(record, "CornerRadius") != _recipe_number(recipe, "CornerIsotropyRadius"):
        raise ValueError("Required corner-ball radius differs from the recipe CornerIsotropyRadius")
    per_corner = record["PerCorner"]
    if (len(per_corner) != len(corners) or any(not isinstance(row, dict) for row in per_corner) or
            sorted(row.get("Point") for row in per_corner) != sorted(corners) or
            any(_count(row.get("Tetrahedra"), "Required corner tetrahedra") <= 0
                for row in per_corner)):
        raise ValueError("Required corner balls differ from the recipe semantic corners")
    layer = recipe.get("EdgeLayer")
    if layer is None:
        if record.get("LayerRequiredReach") is not None or record["PerSpan"]:
            raise ValueError("Required edge layer recorded without a recipe edge layer")
        regions = sum(row["Tetrahedra"] for row in per_corner)
    else:
        reach = (_recipe_number(layer, "LayerThickness") *
                 (1.0 + _recipe_number(layer, "RowZigzag")) + _recipe_number(layer, "EdgeSize"))
        if (_recipe_number(record, "LayerRequiredReach") != reach or
                _recipe_number(layer, "RequiredReach") != reach or
                len(record["PerSpan"]) != layer.get("SpanCount") or
                any(not isinstance(row, dict) or
                    _count(row.get("Tetrahedra"), "Required span tetrahedra") <= 0
                    for row in record["PerSpan"]) or
                [row.get("Span") for row in record["PerSpan"]] != layer.get("Spans")):
            raise ValueError("Required edge-layer record differs from the recipe edge layer")
        regions = sum(row["Tetrahedra"] for row in per_corner + record["PerSpan"])
    if not max(row["Tetrahedra"] for row in per_corner) <= len(indices) <= regions:
        raise ValueError("Required tetrahedron counts do not cover the list")
    gates = {}
    for option in REQUIRED_REGION_GATE_OPTIONS:
        seed_value = _option_or_default(seed_report["Command"], option, None)
        restorer_value = _option_or_default(restoration_report["Command"], option, None)
        if seed_value != restorer_value or not math.isfinite(seed_value) or seed_value <= 0:
            raise ValueError(f"Seed command {option} differs from the label-restoration command")
        gates[option] = seed_value
    quality = census.get("SeedQualityOptimization") if isinstance(census, dict) else None
    if (not isinstance(quality, dict) or
            _recipe_number(quality, "MaximumCornerAspect") != gates["--maximum-corner-aspect"] or
            _recipe_number(quality, "MinimumScaledJacobian") != gates["--minimum-scaled-jacobian"] or
            _recipe_number(quality, "MaximumJacobianCondition") !=
            gates["--maximum-jacobian-condition"] or
            _recipe_number(quality, "DisplacementBoundOverNormal") !=
            gates["--maximum-quality-displacement-over-normal"] or
            _count(quality.get("RequiredCellsBelowGateAfter"), "Required cells below the gate") != 0 or
            _count(quality.get("RequiredCellsAboveConditionAfter"),
                   "Required cells above the condition gate") != 0 or
            _count(quality.get("RequiredTetrahedra"), "Seed required tetrahedra") != len(indices) or
            not isinstance(quality.get("CornerAspectsAfter"), list) or
            len(quality["CornerAspectsAfter"]) != len(corners) or
            any(isinstance(value, bool) or not isinstance(value, (int, float)) or
                not math.isfinite(value) or value > gates["--maximum-corner-aspect"]
                for value in quality["CornerAspectsAfter"])):
        raise ValueError("Seed census does not record a gated required-region optimization")
    # Round 3 class (1): the seed records null when no cell is scaled-Jacobian-gated (a
    # corner-free coupon under the layer rule); accepted only with ScaledJacobianGateCells 0.
    if quality.get("RequiredMinimumScaledJacobianAfter") is None:
        if _count(quality.get("ScaledJacobianGateCells"), "Scaled-Jacobian-gated cells") != 0:
            raise ValueError("Seed required region lacks RequiredMinimumScaledJacobianAfter with gated cells")
    elif (_recipe_number(quality, "RequiredMinimumScaledJacobianAfter") <
            gates["--minimum-scaled-jacobian"]):
        raise ValueError("Seed required region is below the scaled-Jacobian gate")
    if (_recipe_number(quality, "RequiredMaximumJacobianConditionAfter") >
            gates["--maximum-jacobian-condition"]):
        raise ValueError("Seed required region is above the Jacobian condition gate")
    restoration = _restoration_quality_report(restoration_report)
    if _count(restoration.get("RequiredTetrahedra"), "Restored required tetrahedra") != len(indices):
        raise ValueError("Label restoration found a different number of required tetrahedra "
                         "than the recipe record")
    # None: no required cell is judged by the condition (all are layer cells under
    # the layer quality rule).
    if (restoration.get("RequiredMaximumJacobianCondition") is not None and
            _recipe_number(restoration, "RequiredMaximumJacobianCondition") >
            gates["--maximum-jacobian-condition"]):
        raise ValueError("Label restoration found required tetrahedra above the Jacobian "
                         "condition gate")
    validate_edge_layer_quality_rule(seed_report, restoration_report, layer, quality, restoration)
    return record


def validate_edge_layer_quality_rule(seed_report, restoration_report, layer, quality, restoration):
    """The seed and the label restorer apply the same edge-layer quality rule bound
    (EDGE_LAYER_QUALITY_RULE_OPTION) or none.  With a bound: the recipe records an
    edge layer, the census SeedQualityOptimization.EdgeLayerQualityRule records that
    bound with no layer cell above it or below the roundoff floor, and the restorer's
    EdgeLayerQuality passes with the same bound.  Without: neither record carries a
    rule.  Returns the bound or None."""
    seed_value = _option_or_default(seed_report["Command"], EDGE_LAYER_QUALITY_RULE_OPTION, 0.0)
    restorer_value = _option_or_default(restoration_report["Command"],
                                        EDGE_LAYER_QUALITY_RULE_OPTION, 0.0)
    if seed_value != restorer_value:
        raise ValueError(f"Seed command {EDGE_LAYER_QUALITY_RULE_OPTION} differs from the "
                         "label-restoration command")
    seed_rule = quality.get("EdgeLayerQualityRule")
    restorer_rule = restoration.get("EdgeLayerQuality")
    if seed_value == 0.0:
        if seed_rule is not None or restorer_rule is not None:
            raise ValueError("Edge layer quality rule recorded without the bound option")
        return None
    if layer is None:
        raise ValueError("Edge layer quality rule without a recipe edge layer")
    if not math.isfinite(seed_value) or seed_value <= 1.0:
        raise ValueError("Edge layer maximum edge aspect must be a finite bound above 1")
    if (not isinstance(seed_rule, dict) or
            _recipe_number(seed_rule, "MaximumEdgeAspect") != seed_value or
            _recipe_number(seed_rule, "ScaledJacobianRoundoffFloor") != EDGE_LAYER_ORIENTATION_FLOOR or
            _count(seed_rule.get("LayerCells"), "Seed layer cells") <= 0 or
            _count(seed_rule.get("CellsAboveBoundAfter"), "Layer cells above the bound") != 0 or
            _count(seed_rule.get("CellsBelowRoundoffFloorAfter"), "Flat layer cells") != 0 or
            _recipe_number(seed_rule, "MaximumEdgeAspectAfter") > seed_value):
        raise ValueError("Seed census does not record a gated edge-layer quality rule")
    if (not isinstance(restorer_rule, dict) or restorer_rule.get("Passes") is not True or
            _recipe_number(restorer_rule, "MaximumEdgeAspectBound") != seed_value or
            _recipe_number(restorer_rule, "ScaledJacobianRoundoffFloor") != EDGE_LAYER_ORIENTATION_FLOOR or
            _recipe_number(restorer_rule, "MaximumEdgeAspect") > seed_value or
            restorer_rule.get("PositiveOrientation") is not True or
            _count(restorer_rule.get("Cells"), "Restored layer cells") <= 0 or
            restorer_rule.get("Rule") != EDGE_LAYER_QUALITY_RULE):
        raise ValueError("Label restoration does not record a passing edge-layer quality rule")
    return seed_value


TRACE_BASIS_RATIO_OPTION = "--trace-basis-size-ratio"
# Report-only needle classification fraction recorded by the census (decision 43):
# defined once in trace_basis.py (origin documented there); the census value must equal it.
TRACE_BASIS_NEEDLE_ALTITUDE_OVER_SHORTEST_EDGE = NEEDLE_ALTITUDE_OVER_SHORTEST_EDGE


def bound_trace_basis(report):
    """Digests of the trace basis inputs a stage bound, keyed by source role, or None
    when it bound none; a partial binding fails closed."""
    inputs = report["Inputs"]
    present = {name for name in TRACE_BASIS_INPUTS if name in inputs}
    if not present:
        return None
    if present != set(TRACE_BASIS_INPUTS):
        raise ValueError("Stage bound an incomplete trace basis")
    return {role: inputs[name]["SHA256"] for name, role in TRACE_BASIS_INPUTS.items()}


def _optional_ratio(command, bound):
    positions = [index for index, value in enumerate(command) if value == TRACE_BASIS_RATIO_OPTION]
    if not bound:
        if positions:
            raise ValueError("Stage command passes a trace basis size ratio without a bound trace basis")
        return None
    try:
        ratio = float(_option_value(command, TRACE_BASIS_RATIO_OPTION))
    except ValueError as error:
        raise ValueError(f"Stage command {TRACE_BASIS_RATIO_OPTION} is not a bound number: {error}") from error
    if not math.isfinite(ratio) or ratio <= 0:
        raise ValueError("Trace basis size ratio must be a positive finite number")
    return ratio


def validate_trace_basis_sizing(seed_report, metric_report, recipe, census):
    """The seed and the metric stage bind the same trace basis (or none); when bound,
    both commands pass the same dimensionless TraceBasisSizeRatio, the census and the
    recipe record that ratio, the bound input digests and the same mesh-frame basis
    triangles, and the recipe records the cut-surface size statistics."""
    seed_basis = bound_trace_basis(seed_report)
    metric_basis = bound_trace_basis(metric_report)
    if seed_basis != metric_basis:
        raise ValueError("Seed and metric stages bound different trace bases")
    seed_ratio = _optional_ratio(seed_report["Command"], seed_basis is not None)
    metric_ratio = _optional_ratio(metric_report["Command"], metric_basis is not None)
    record = recipe.get("TraceBasisSizing")
    census_record = census.get("TraceBasisSizing")
    if seed_basis is None:
        if record is not None or census_record is not None:
            raise ValueError("Trace basis sizing recorded without a bound trace basis")
        return None
    if (not isinstance(record, dict) or not isinstance(census_record, dict) or
            seed_ratio != metric_ratio or _recipe_number(record, "Ratio") != seed_ratio or
            _recipe_number(census_record, "Ratio") != seed_ratio or
            record.get("RatioIsDimensionless") is not True or
            record.get("InputSHA256") != seed_basis or census_record.get("InputSHA256") != seed_basis):
        raise ValueError("Trace basis sizing records differ from the bound trace basis and ratio")
    recipe_triangles = record.get("MeshFrameTriangles")
    census_triangles = census_record.get("MeshFrameTriangles")
    if (not isinstance(recipe_triangles, list) or not recipe_triangles or
            not isinstance(census_triangles, list) or
            len(recipe_triangles) != len(census_triangles) or
            record.get("Triangles") != len(recipe_triangles)):
        raise ValueError("Trace basis sizing triangles are missing or inconsistent")
    scale = 0.0
    for name in ("Lower", "Upper"):
        values = record.get(name)
        if (not isinstance(values, list) or len(values) != 3 or
                any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v)
                    for v in values)):
            raise ValueError("Trace basis sizing box is invalid")
    scale = max(u - l for l, u in zip(record["Lower"], record["Upper"]))
    if scale <= 0:
        raise ValueError("Trace basis sizing box is invalid")
    for first, second in zip(recipe_triangles, census_triangles):
        if (not isinstance(first, list) or not isinstance(second, list) or len(first) != 3 or
                len(second) != 3 or any(
                    not isinstance(a, list) or not isinstance(b, list) or len(a) != 3 or len(b) != 3 or
                    any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x)
                        for x in a + b) or
                    any(abs(x - y) > 1e-8 * scale for x, y in zip(a, b))
                    for a, b in zip(first, second))):
            raise ValueError("Trace basis sizing triangles differ between the seed census and the recipe")
    sizes = record.get("CutSurfaceSize")
    if (not isinstance(sizes, dict) or
            any(_recipe_number(sizes, name) <= 0 for name in ("Minimum", "Median", "Maximum")) or
            _count(sizes.get("CutTriangles"), "Cut-surface triangles") <= 0 or
            _count(record.get("BasisEdgesBelowFarSize"), "Basis edges below the far size") < 0):
        raise ValueError("Trace basis sizing lacks the cut-surface size statistics")
    _validate_trace_basis_edges(recipe, recipe_triangles, seed_basis, scale)
    return record


def _validate_trace_basis_edges(recipe, triangles, digests, scale):
    """The recipe's TraceBasisEdges (the audit's source-driven band lines) are exactly
    the unique edges of the recorded basis triangles placed by the recipe's contract
    transform, with the bound input digests.  Matching is one-to-one: every recorded
    segment consumes one expected edge, so a duplicated segment cannot stand in for a
    missing one."""
    record = recipe.get("TraceBasisEdges")
    if (not isinstance(record, dict) or record.get("InputSHA256") != digests or
            not isinstance(record.get("Segments"), list)):
        raise ValueError("Restoration recipe lacks the bound trace basis edges")
    placement = recipe.get("SemanticContract", {}).get(
        "RigidTransform", [1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 1.])
    matrix = [[float(placement[4 * row + column]) for column in range(4)] for row in range(4)]
    def place(point):
        return [sum(matrix[row][column] * point[column] for column in range(3)) + matrix[row][3]
                for row in range(3)]
    expected = {}
    for triangle in triangles:
        placed = [tuple(place(point)) for point in triangle]
        for a, b in ((0, 1), (1, 2), (2, 0)):
            key = tuple(sorted((placed[a], placed[b])))
            expected[key] = None
    segments = record["Segments"]
    if record.get("Count") != len(segments) or len(segments) != len(expected):
        raise ValueError("Restoration recipe trace basis edges differ from the basis triangles")
    unmatched = list(expected)
    for segment in segments:
        if (not isinstance(segment, list) or len(segment) != 6 or
                any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v)
                    for v in segment)):
            raise ValueError("Restoration recipe trace basis edge is invalid")
        ends = tuple(sorted((tuple(segment[:3]), tuple(segment[3:]))))
        match = next((index for index, key in enumerate(unmatched)
                      if all(abs(x - y) <= 1e-8 * scale
                             for p, q in zip(ends, key) for x, y in zip(p, q))), None)
        if match is None:
            raise ValueError("Restoration recipe trace basis edges differ from the basis triangles")
        del unmatched[match]
    if unmatched:
        raise ValueError("Restoration recipe trace basis edges differ from the basis triangles")
    return record


# Options of the Gmsh-only build command (mesh_spatial_coupon.jl --prism-tubes) bound
# to the census records: the tube recipe (inner ring size = corner size, ring ratio,
# extrusion spacing, band sizes and growth), the corner ball, the element cap and
# the fail-closed quality gates the mesher applies per element type.
GMSH_BUILD_TUBE_OPTION = "--prism-tubes"
GMSH_BUILD_RECIPE_OPTIONS = {"--edge-size": "InnerSize", "--edge-growth-ratio": "GrowthRatio",
                             "--lc-fine": "NormalSize", "--lc-far": "FarSize",
                             "--far-growth": "FarGrowth"}
# Coupon-scale size bound (mesh_spatial_coupon.jl SIZE_BOUND_RULE): the census
# TangentialSize is min(--lc-tangent, --lc-far) - FarSize = FarSizeOverRadius x Radius is
# the coarsest size the coupon admits, so the along-edge tube spacing never exceeds it;
# the request and the bound are recorded in census SizeBounds.
GMSH_BUILD_TANGENTIAL_OPTION = "--lc-tangent"
# Tube cross-section option with the mesher's default (mesh_spatial_coupon.jl
# --tube-sector-degrees) bound to the census Section.SectorDegrees.
GMSH_BUILD_SECTOR_OPTION = ("--tube-sector-degrees", 30.0)
# Mesher constants of the tube recipe bound to the census records: the pyramid
# height over the outermost ring size (TUBE_PYRAMID_HEIGHT_OVER_OUTER_RING), the
# band law's radial growth inside the protected distance (BAND_RADIAL_GROWTH) and
# the protected distance over NormalSize (BAND_PROTECTED_DISTANCE_OVER_NORMAL).
GMSH_BUILD_PYRAMID_HEIGHT_OVER_OUTER_RING = 0.5
GMSH_BUILD_BAND_RADIAL_GROWTH = 1.0
GMSH_BUILD_BAND_PROTECTED_DISTANCE_OVER_NORMAL = 2.0
GMSH_BUILD_GATE_OPTIONS = {"--maximum-corner-aspect": "MaximumCornerAspect",
                           "--minimum-scaled-jacobian": "MinimumScaledJacobian",
                           "--maximum-jacobian-condition": "MaximumJacobianCondition",
                           "--maximum-quality-displacement-over-normal": "DisplacementBoundOverNormal"}
GMSH_BUILD_ELEMENT_CAP_OPTION = "--max-elements"
GMSH_BUILD_VOLUME_TYPES = ("Tetrahedron", "Prism", "Pyramid")
# Mesher design round 2 F5-A (supervisor decisions 351 / 358 / 363): the verdict bound of
# the INVARIANT (non-perpendicular) semantic corners, kappa_reg <= CornerShapeGate =
# min(E_pop, 5.0) (manifest Gates.CornerShapeGate -> the mesher's --corner-shape-gate);
# legacy corners keep MaximumCornerAspect. A build with an invariant corner needs it; a
# rectilinear build does not carry it.
GMSH_BUILD_CORNER_SHAPE_GATE_OPTION = "--corner-shape-gate"
# Mesher design round 2 F6 (decisions 437 / 443): the smallest ring count validated by (F) for
# the coupon kind, the mesher's qualified ring-count range (Gates.MinimumQualifiedRings of the
# manifest per kind -> run_gmsh_only_case -> this option).
GMSH_BUILD_MINIMUM_QUALIFIED_RINGS_OPTION = "--minimum-qualified-rings"
RINGS_PER_SIDE_RULE = ("mesher design round 2 F6 (decisions 347 / 349 / 437 / 443): every tubed side's ring count K_side is "
                       "the largest K with r_K + h_K <= FacingBound = min(TransverseBound, FacingWidth / 2), FacingWidth its "
                       "smallest facing width across the metal; a side at the coupon's Rings records nothing, every other "
                       "Tubes[].Rings / FacingWidth / FacingBound; K_side below MinimumQualifiedRings fails closed at "
                       "ScopeGuard[UnqualifiedRingCount]")


def ring_count_within(inner_size, ratio, bound):
    """The largest K with r_K + h_K <= bound (the mesher's tube_ring_count; 0 when none fits)."""
    rings = 0
    while True:
        k = rings + 1
        radius = inner_size * (ratio**k - 1.0) / (ratio - 1.0)
        if radius + inner_size * ratio**(k - 1) > bound:
            return rings
        rings = k


def validate_tube_rings_per_side(tubes, rows, command):
    """Design round 2 F6: the per-side ring records. Section.Rings is the coupon's (process)
    count; a row carrying Rings is a side reduced by its facing width - Rings in [1,
    Section.Rings), equal to the largest K with r_K + h_K <= FacingBound, FacingBound ==
    min(Section.TransverseBound, FacingWidth / 2) < the bound of Section.Rings, FacingWidth >=
    Section.MetalFacingWidth, and Rings >= the command's --minimum-qualified-rings (required
    when any row is reduced); Section.MinimumRings / ReducedSides / FacingBound /
    MinimumQualifiedRings follow the rows and the command. A census without the records
    (before F6) validates when no row carries Rings."""
    section = tubes.get("Section")
    if not isinstance(section, dict):
        raise ValueError("Prism tube section is missing")
    coupon_rings = _count(section.get("Rings"), "Tube rings")
    inner = tubes["InnerSize"]
    ratio = tubes["GrowthRatio"]
    # The measurement-only ring cap (--maximum-rings; F6 2.3): the coupon's count is min(law, cap).
    cap_option = (int(_option_or_default(command, "--maximum-rings", None)) if "--maximum-rings" in command else None)
    recorded_cap = section.get("RingsCap")
    if (recorded_cap is None) != (cap_option is None) or (recorded_cap is not None and recorded_cap != cap_option):
        raise ValueError("Prism tube section RingsCap differs from the build command --maximum-rings")
    if "TransverseBound" in section:
        law = ring_count_within(inner, ratio, _census_number(section, "TransverseBound", "Tube section"))
        if coupon_rings != (min(law, cap_option) if cap_option is not None else law):
            raise ValueError("Prism tube Rings do not follow the transverse bound (and the measurement-only cap)")
    reduced = [row for row in rows if "Rings" in row]
    minimum_option = (int(_option_or_default(command, GMSH_BUILD_MINIMUM_QUALIFIED_RINGS_OPTION, None))
                      if GMSH_BUILD_MINIMUM_QUALIFIED_RINGS_OPTION in command else None)
    if minimum_option is not None and minimum_option < 1:
        raise ValueError(f"Gmsh-only build command carries an invalid {GMSH_BUILD_MINIMUM_QUALIFIED_RINGS_OPTION}")
    if "MinimumRings" not in section:
        if reduced:
            raise ValueError("Prism tube rows record per-side rings without the section's per-side record")
        return
    transverse = _census_number(section, "TransverseBound", "Tube section")
    facing_bound = _census_number(section, "FacingBound", "Tube section")
    minimum_rings = _count(section.get("MinimumRings"), "Tube section MinimumRings")
    reduced_count = _count(section.get("ReducedSides"), "Tube section ReducedSides")
    recorded_minimum = section.get("MinimumQualifiedRings")
    if (recorded_minimum is None) != (minimum_option is None) or (
            recorded_minimum is not None and recorded_minimum != minimum_option):
        raise ValueError("Prism tube section MinimumQualifiedRings differs from the build command")
    if not isinstance(section.get("RingsPerSideRule"), str) or "UnqualifiedRingCount" not in section["RingsPerSideRule"]:
        raise ValueError("Prism tube section lacks the per-side ring rule")
    if reduced and minimum_option is None:
        raise ValueError("Prism tube rows are reduced by their facing width without a qualified ring-count range")
    smallest_bound = transverse
    smallest_rings = coupon_rings
    for row in reduced:
        rings = _count(row.get("Rings"), "Tube row rings")
        width = _census_number(row, "FacingWidth", "Tube row")
        bound = _census_number(row, "FacingBound", "Tube row")
        if not 1 <= rings < coupon_rings or rings < minimum_option:
            raise ValueError("Prism tube row rings lie outside [MinimumQualifiedRings, Section.Rings)")
        if abs(bound - min(transverse, 0.5 * width)) > 1e-12 * transverse or ring_count_within(inner, ratio, bound) != rings:
            raise ValueError("Prism tube row rings do not follow the largest K with r_K + h_K <= min(TransverseBound, FacingWidth / 2)")
        smallest_bound = min(smallest_bound, bound)
        smallest_rings = min(smallest_rings, rings)
    tubes_per_side = _count(section.get("TubesPerSide"), "Tubes per side")
    width = section.get("MetalFacingWidth")
    # Section.FacingBound is the coupon's smallest per-side bound min(TransverseBound, w / 2) over
    # EVERY tubed side (a bound below the transverse bound need not reduce a side: the ring law
    # steps, and a measurement cap may already hold the count).
    expected_bound = transverse if width is None else min(transverse, 0.5 * _census_number(section, "MetalFacingWidth", "Tube section"))
    if (minimum_rings != smallest_rings or reduced_count * tubes_per_side != len(reduced) or
            abs(facing_bound - expected_bound) > 1e-12 * transverse or smallest_bound < expected_bound * (1.0 - 1e-12)):
        raise ValueError("Prism tube section MinimumRings / ReducedSides / FacingBound do not follow the rows")
    if width is not None and any(row["FacingWidth"] < width * (1.0 - 1e-12) for row in reduced):
        raise ValueError("Prism tube row facing width lies below the section's MetalFacingWidth")
    if reduced and width is None:
        raise ValueError("Prism tube rows are reduced by their facing width but the section records none")

CORNER_SHAPE_GATE = "CornerShapeGate"
CORNER_MEASURES = {"Legacy": "VertexFrameCondition", "Invariant": "RegularCondition"}
INVARIANT_CORNER_TARGET = 3.8
# Round 3 class (1) (decision 510): a corner-free build (the contract's Derivation.CornerFree,
# SemanticCorners []) has no corner ball and no corner to judge - MaximumCornerAspect /
# CornerShapeGate are NOT applicable and the census says so explicitly (the mesher's
# CORNER_GATES_NOT_APPLICABLE, one spelling); 0 CornerMeasures rows, 0 seed corner census
# rows and 0 tube cap regions are accepted only under it.
CORNER_GATES_RECORD = "CornerGates"
CORNER_GATES_NOT_APPLICABLE = "NotApplicable (0 semantic corners, Derivation.CornerFree)"
# The round-2 census records (F5-A / F5-B / SEAM; decisions 363 / 365 / 368) and the rule
# round of a census (supervisor decision 392 MINOR-6 ruling, the 4.1 / 4.2 convention): a
# census declaring NONE of them is a PRE-RULE census and validates only as one - judged by
# the pre-rule corner rule (every CornerAspectsAfter <= MaximumCornerAspect; no
# --corner-shape-gate, no contract invariant corner, no seam census) and bound to its
# declared tool version (its build report's mesher digest is not the round-2 mesher's beside
# this module: ROUND2_MESHER); a census declaring any of them is a ROUND-2 census and every
# round-2 record is required (a missing one fails closed).
ROUND2_CENSUS_RECORDS = ("SemanticCornerKinds", "InvariantCorners", "TipBisectors", "ThinSheetSeams")
ROUND2_OPTIMIZATION_RECORDS = ("CornerMeasures", "InvariantCorners", "CornerShapeGate",
                               "InvariantCornerTarget", "InvariantCornerRule", "CornerReconnectionRule")
CENSUS_RULE_ROUNDS = ("pre-rule", "round-2")
ROUND2_MESHER = Path(__file__).resolve().parent / "mesh_spatial_coupon.jl"


def census_rule_round(census, optimization):
    """'round-2' when the census or its SeedQualityOptimization declares any round-2 record,
    'pre-rule' when it declares none (a census written before mesher design round 2)."""
    declared = ([name for name in ROUND2_CENSUS_RECORDS if name in census] +
                [name for name in ROUND2_OPTIMIZATION_RECORDS
                 if isinstance(optimization, dict) and name in optimization])
    return "round-2" if declared else "pre-rule"


def _mesher_digest(build_report):
    """The build report's recorded mesher digest (None when absent)."""
    tools = build_report.get("Tools")
    mesher = tools.get("mesher") if isinstance(tools, dict) else None
    digest = mesher.get("SHA256") if isinstance(mesher, dict) else None
    return digest if isinstance(digest, str) and digest else None


def validate_pre_rule_census_tool(build_report):
    """A pre-rule census is bound to its declared tool: the build report's mesher digest must
    be recorded and must differ from the round-2 mesher beside this module (which writes
    every round-2 record, so a census from it lacking them is not a pre-rule census)."""
    tools = build_report.get("Tools")
    mesher = tools.get("mesher") if isinstance(tools, dict) else None
    digest = mesher.get("SHA256") if isinstance(mesher, dict) else None
    if not isinstance(digest, str) or not digest:
        raise ValueError("A pre-rule build census (no round-2 record) needs the build report's mesher digest")
    if digest == sha256(ROUND2_MESHER):
        raise ValueError("Build census from the round-2 mesher lacks the round-2 records "
                         "(CornerMeasures / SemanticCornerKinds / InvariantCorners / TipBisectors / ThinSheetSeams)")
    return digest


def validate_corner_measures(optimization, corners, maximum_corner_aspect, corner_shape_gate,
                             semantic_invariant, rule_round="round-2"):
    """The seed optimization's per-corner verdict records (CornerMeasures: one per semantic
    corner, Point among the corners) judged by kind: a Legacy corner's VertexFrameCondition
    After <= MaximumCornerAspect (target 0.95 x the bound); an Invariant corner's
    RegularCondition After <= CornerShapeGate (the command's --corner-shape-gate, required
    when any corner is invariant), Target 3.8, its BridgingSlivers counts recorded (supervisor
    decision 365: the measure is the verdict, a candidate above the gate the reconnection
    trigger, so none may remain above the gate); the set of Invariant corners equals the
    contract's Derivation.InvariantCorners (fail closed on a disagreement); CornerAspectsAfter
    mirrors the per-corner After values. A PRE-RULE census (rule_round 'pre-rule', decision
    392 MINOR-6: no round-2 record declared) is judged by the pre-rule rule alone - every
    CornerAspectsAfter <= MaximumCornerAspect, no --corner-shape-gate, no contract invariant
    corner (such a contract needs a round-2 build). Returns the number of invariant corners."""
    if rule_round not in CENSUS_RULE_ROUNDS:
        raise ValueError(f"Unknown census rule round {rule_round!r}")
    after_values = optimization.get("CornerAspectsAfter")
    if rule_round == "pre-rule":
        if corner_shape_gate is not None:
            raise ValueError("A pre-rule build census (no CornerMeasures) cannot be judged by --corner-shape-gate")
        if semantic_invariant:
            raise ValueError("A pre-rule build census (no CornerMeasures) cannot judge the contract's invariant "
                             "corners: rebuild under the round-2 mesher")
        if (not isinstance(after_values, list) or len(after_values) != len(corners) or
                any(value > maximum_corner_aspect for value in after_values)):
            raise ValueError("Pre-rule build census corner aspects exceed MaximumCornerAspect")
        return 0
    measures = optimization.get("CornerMeasures")
    if (not isinstance(measures, list) or len(measures) != len(corners) or
            any(not isinstance(row, dict) for row in measures)):
        raise ValueError("Seed optimization does not record one CornerMeasures row per semantic corner")
    points = [row.get("Point") for row in measures]
    if sorted(points) != sorted(corners):
        raise ValueError("Seed CornerMeasures points differ from the semantic corners")
    if not isinstance(after_values, list) or len(after_values) != len(measures):
        raise ValueError("Seed CornerAspectsAfter differs from the CornerMeasures rows")
    invariant_points = []
    for row, after in zip(measures, after_values):
        kind = row.get("Kind")
        if kind not in CORNER_MEASURES or row.get("Measure") != CORNER_MEASURES[kind]:
            raise ValueError(f"Seed corner measure row has an unknown kind or measure: {kind!r}")
        value = _census_number(row, "After", "Seed corner measure")
        if value != after or _census_number(row, "Before", "Seed corner measure") <= 0.0:
            raise ValueError("Seed corner measure After differs from CornerAspectsAfter")
        slivers = row.get("BridgingSlivers")
        if (not isinstance(slivers, dict) or
                _count(slivers.get("Before"), "Bridging slivers before") < 0 or
                _count(slivers.get("After"), "Bridging slivers after") < 0 or
                _count(slivers.get("AboveGateAfter"), "Bridging slivers above the gate") < 0):
            raise ValueError("Seed corner measure row lacks the BridgingSlivers counts")
        if kind == "Legacy":
            expected_gate, expected_target = maximum_corner_aspect, 0.95 * maximum_corner_aspect
        else:
            if corner_shape_gate is None:
                raise ValueError("An invariant semantic corner was judged without --corner-shape-gate")
            expected_gate, expected_target = corner_shape_gate, INVARIANT_CORNER_TARGET
            invariant_points.append(row["Point"])
        if (_census_number(row, "Gate", "Seed corner measure") != expected_gate or
                _census_number(row, "Target", "Seed corner measure") != expected_target):
            raise ValueError(f"Seed corner measure gate / target differ from the {kind} rule")
        if value > expected_gate or slivers["AboveGateAfter"] != 0 or row.get("Passed") is not True:
            raise ValueError(f"Seed {kind} corner {row['Point']} fails its gate: {value} > {expected_gate} "
                             f"or {slivers['AboveGateAfter']} bridging-sliver cells above the gate remain")
    if sorted(invariant_points) != sorted(semantic_invariant):
        raise ValueError("Seed invariant corners differ from the contract's Derivation.InvariantCorners")
    if optimization.get("InvariantCorners") != len(invariant_points):
        raise ValueError("Seed optimization InvariantCorners count differs from its rows")
    recorded_gate = optimization.get(CORNER_SHAPE_GATE)
    if (corner_shape_gate is None) != (recorded_gate is None) or (
            corner_shape_gate is not None and
            _census_number(optimization, CORNER_SHAPE_GATE, "Seed optimization") != corner_shape_gate):
        raise ValueError("Seed optimization CornerShapeGate differs from the build command")
    return len(invariant_points)

# Recipe scope (supervisor decision 48; mesh_spatial_coupon.jl RECIPE_SCOPE_*): the
# input classes the prism-tube recipe builds, the guarded classes it fails closed on
# (a guard's error message carries "ScopeGuard[<id>]"; DetectedFrom "inputs" when the
# class is visible in the frozen inputs, "build" when only a derived quantity shows it)
# and the classification of an input from its frozen inputs.  The census Scope block
# records the same lists and is bound here; the build drivers record an unsupported
# class distinctly from any other failure with these ids.  UntubedShortEdges (block (b)
# design A9 family 3, supervisor decision 304) is the one supported class known only at
# the build (a side's tube interval against the derived corner clearances): it is
# exhibited by the census Scope.UntubedEdges records, not by the input classification.
RECIPE_SCOPE_RECIPE = "prism-tubes"
RECIPE_SCOPE_SUPPORTED_CLASSES = ("ArcSides", "ContinuationVertices", "DeviceFootprint", "DownwardLayers",
                                  "ExteriorLoops", "HoleLoops", "MultipleConductors",
                                  "MultipleLayers", "MultipleSlots", "ThinMetal", "TraceBasis",
                                  "UntubedShortEdges")
RECIPE_SCOPE_GUARDS = {
    "ArcTubeRadiusVsCurvature": "build",
    # Decision 391 MAJOR-2 (ii) / decision 437 (3) / round 3 decision 556: arc face ends and arc
    # joints beyond the BUILT ranges (ARC_SMOOTH_JOINT_TURN_BOUND_RADIANS, the per-kind
    # ARC_CORNER_JOINT_TURN_RANGE_RADIANS, ARC_FACE_END_TILT_RANGE_DEGREES; the loop end's
    # ARC_JOINT_TURN_BOUND_RADIANS inside) fail closed (mesh_spatial_coupon.jl spells the same list).
    "ArcFaceEnds": "build", "ArcJointTilt": "build",
    # Mesher design round 3 class (11) interim (decisions 497 / 500 / 510; part M 1.3): a smooth
    # joint of two DISTINCT arc runs has no joint owner (the 32b0083dad90 build failure) - refused
    # by name before any CAD until fix 11 lands.
    "ArcArcJoint": "build",
    "TopRounding": "inputs", "TrenchRounding": "inputs", "SlopedSidewalls": "inputs",
    "NoTrench": "inputs", "ShallowTrench": "build",
    "NarrowTransverseBound": "build", "NarrowHoles": "build", "NarrowLayerGap": "build",
    # Mesher design round 2 F6 (decisions 347 / 349 / 437 / 443): a per-side ring count below the
    # smallest ring count validated by (F) for the coupon kind (NarrowMetal retired into the
    # per-side bound).
    "UnqualifiedRingCount": "build",
    # Mesher design round 2 F2b (decisions 358 / 363 / 437): a face end beyond the validity
    # ceiling of the capped end block (2 h_pyr_end x slope >= lc_cap); round 3 class (8) 8B (part M
    # 5.2) / decision 566: an end above the largest built tilt of its kind
    # (THIN_FACE_END_TILT_BOUND_DEGREES / FABRICATED_FACE_END_TILT_BOUND_DEGREES) fails closed at
    # the same guard.
    "SteepFaceCrossing": "build",
    "FreeEdgeEnds": "build", "FootprintWithoutEdge": "build",
    "FootprintTopology": "build",
    # Mesher design round 3 class (5) interim 5B (part M 2.2): a concave Physical arc of a
    # fabricated coupon leaving the box by less than the trench-collar width (its shrunk collar
    # circle misses the box-face line; the 32dc558f4810 refusal) and the untested line-arc corner
    # kink whose offsets do not meet.
    "CollarFaceEnd": "build"}
# The largest arc-joint turn any full build has exercised: the loop end 1b26671c9080's
# smooth arc / line joints (3.0e-8 .. 1.6e-6 rad on the signature); the mesher's
# ScopeGuard[ArcJointTilt] bound (decision 391 MAJOR-2 (ii)).
ARC_JOINT_TURN_BOUND_RADIANS = 1.6e-6
# Mesher design round 2b (decision 437 (3)): the guards are LIFTED to the ranges the four synthetic
# full builds tested (the mesher's ARC_SMOOTH_JOINT_TURN_BOUND / ARC_CORNER_JOINT_TURN_RANGE /
# ARC_FACE_END_TILT_RANGE, spelled identically): a smooth arc joint turning by <= 5e-5 rad, a corner
# arc joint turning within the coupon KIND's built range, an arc box-face cut end whose tilt lies in
# the BUILT range build; everything beyond fails closed at the same guards.
ARC_SMOOTH_JOINT_TURN_BOUND_RADIANS = 5.0e-5
# Round 3 class (6) (decisions 491 / 556): the corner-joint range is the BUILT range PER KIND. THIN
# 2e-4 rad .. 30 degrees = the round-2b corner2e-4 / kink30 production-size builds (PBS 57706 /
# 57892); FABRICATED 2e-4 rad .. 15 degrees = the round-3 B4 production-size builds of the 2e-4-rad
# corner and the 0.1 / 5.33 / 8.53 / 15-degree kinks, once the E4 root cause (the tube cap-ray node
# order on a periodic trench-wall face) was fixed in the mesher; the fabricated 20 / 25 / 30-degree
# kinks fail the quality gates (a named follow-up, not admitted).
ARC_CORNER_JOINT_TURN_RANGE_RADIANS = {"fabricated": (2.0e-4, math.radians(15.0)),
                                       "thin": (2.0e-4, math.radians(30.0))}
# Mesher design round 3 (decisions 497 ERRATUM / 510 MAJOR-1 / 510 O6; DESIGN R9): the admitted arc
# face-end tilt range is the BUILT range (lowest built, largest built) in degrees - ONE constant, the
# mesher's ARC_FACE_END_TILT_RANGE in radians. Provenance: 0.1 degrees = fe0p1 of the round-3 B2
# record run (fe0p1 / fe0p5 / fe2p1 / fe8 fab + thin at the production sizes); 70 degrees = fe70 fab
# + thin of the round-2b record run (PBS 57706 / 57892; the thin top is the largest built THIN
# angle). Round 2b had built 15 / 45 / 70 only and admitted (0, 15) untested (the 497 erratum).
ARC_FACE_END_TILT_RANGE_DEGREES = (0.1, 70.0)
# Round 3 class (8) interim 8B (part M 5.2; E3, decision 475 (7)): the largest BUILT thin face-end
# tilt (the V10 thin 70-degree production-size build); a thin end above it fails closed at
# ScopeGuard[SteepFaceCrossing] until the class-8 fix builds it.
THIN_FACE_END_TILT_BOUND_DEGREES = 70.0
# Round 3 B4 (decision 566): the fabricated tested-range bound - the largest fabricated straight face-end
# tilt BUILT under the code that runs (the V10 fabricated 70-degree production build; under 8A the derived
# SteepFaceCrossing ceiling no longer binds for the fabricated kind, so this bound is the admission); the
# B4 part-2 record run raises it to the largest fabricated tilt it builds and passes under 8A.
FABRICATED_FACE_END_TILT_BOUND_DEGREES = 70.0
# Round 3 class (8), 8A-bitwise (part M 5.2; decisions 510 O8 / 563): above the kind's largest BUILT
# face-end tilt the end block's pyramids take h_pyr_end = min(h_pyr, TangentialSize / (4 slope)) (the
# mesher's FACE_END_DERIVED_APEX_ABOVE, spelled identically: 70 degrees fabricated = the V10 fabricated
# 70-degree production build, 70 degrees thin = the V10 thin 70-degree production build = the thin
# admission top THIN_FACE_END_TILT_BOUND_DEGREES; part M 5.2's "45 thin" predates B2's 70); at or
# below it h_pyr_end = h_pyr (the A2 (4) block, bitwise). FaceEnds[].PyramidHeight carries h_pyr_end.
FACE_END_DERIVED_APEX_ABOVE_DEGREES = {"fabricated": 70.0, "thin": THIN_FACE_END_TILT_BOUND_DEGREES}
# Metal thickness option of the mesher command with its default; the top tube of a
# process layer with normal Nz lies at plane + Nz x MetalThickness (decision 48).
GMSH_BUILD_THICKNESS_OPTION = ("--metal-thickness", 0.1)
# Coupon kinds (the mesher command's second positional token): a fabricated coupon
# carries a top and a bottom tube per straight metal side, a thin coupon (decision 66:
# ThinMetal is a supported class) one sheet tube per side at the process plane, whose
# inner ring size is the recorded cutoff of every thin surface participation.
COUPON_KINDS = ("fabricated", "thin")
TUBES_PER_SIDE = {"fabricated": 2, "thin": 1}
TUBE_EDGES = {"fabricated": ("top", "bottom"), "thin": ("sheet",)}


def build_coupon_kind(command):
    """The coupon kind of a mesher / publisher command: the token of COUPON_KINDS it
    carries (fabricated when it carries none - the production mesher always names its
    kind; fail closed on both)."""
    kinds = [kind for kind in COUPON_KINDS if kind in command]
    if len(kinds) > 1:
        raise ValueError(f"the command names both coupon kinds {kinds}")
    return kinds[0] if kinds else "fabricated"
SCOPE_GUARD_PATTERN = re.compile(r"ScopeGuard\[([A-Za-z]+)\]")
# Process options of the mesher command with the mesher's defaults (a command without
# the option builds the default) that classify an input.
SCOPE_PROCESS_OPTIONS = {"--sidewall-angle": 80.0, "--top-radius": 0.01, "--bottom-radius": 0.01,
                         "--overetch": 0.05}


def scope_guard_in_text(text):
    """The ScopeGuard id a mesher log or message carries (the first one), or None; an
    id outside RECIPE_SCOPE_GUARDS fails closed (the mesher and this contract must
    spell one list)."""
    match = SCOPE_GUARD_PATTERN.search(text)
    if match is None:
        return None
    if match.group(1) not in RECIPE_SCOPE_GUARDS:
        raise ValueError(f"unknown scope guard {match.group(1)}")
    return match.group(1)


def scope_classes(signature_rows, boundary_rows, *, fabricated=True, sidewall_angle=90.0,
                  top_rounding=0.0, trench_rounding=0.0, overetch, device_footprint=False,
                  trace_basis=False):
    """The sorted classes an input exhibits (mesh_spatial_coupon.jl
    exhibited_scope_classes): from the signature rows (Slot, Conductor, Pz, Nz), the
    plan-view boundary rows (Hole, Class) and the process values."""
    classes = []
    holes = {int(row["Loop"]): int(float(row["Hole"])) for row in boundary_rows}
    if any(hole == 0 for hole in holes.values()):
        classes.append("ExteriorLoops")
    if any(hole != 0 for hole in holes.values()):
        classes.append("HoleLoops")
    if any(row["Class"] == "Continuation" for row in boundary_rows):
        classes.append("ContinuationVertices")
    if len({int(row["Slot"]) for row in signature_rows}) > 1:
        classes.append("MultipleSlots")
    if len({int(row["Conductor"]) for row in signature_rows}) > 1:
        classes.append("MultipleConductors")
    layers = {(float(row["Pz"]), int(float(row.get("Nz", 1) or 1))) for row in signature_rows}
    if len(layers) > 1:
        classes.append("MultipleLayers")
    if any(sign < 0 for _, sign in layers):
        classes.append("DownwardLayers")
    if device_footprint:
        classes.append("DeviceFootprint")
    if trace_basis:
        classes.append("TraceBasis")
    if top_rounding > 0.0:
        classes.append("TopRounding")
    if trench_rounding > 0.0:
        classes.append("TrenchRounding")
    if sidewall_angle != 90.0:
        classes.append("SlopedSidewalls")
    if not fabricated:
        classes.append("ThinMetal")
    if overetch == 0.0:
        classes.append("NoTrench")
    # Block (b) design A1 (4): a boundary carrying arc tags (one tagged chord side at least).
    if boundary_rows and "ArcId" in boundary_rows[0] and any(row.get("ArcId") not in (None, "") for row in boundary_rows):
        classes.append("ArcSides")
    return sorted(classes)


# The parts an arc side of signed sweep `sweep` (radians) is split into (mesh_spatial_coupon.jl
# arc_part_count): <= 90-degree parts of equal angular fraction.
def arc_part_count(sweep):
    return max(1, math.ceil(abs(sweep) / (0.5 * math.pi) - 1.0e-9))


def boundary_arc_runs(rows):
    """Per loop (Loop order) the tagged arc runs of plan-view boundary rows (block (b) design
    A1 (4)): [{ArcId, Centre, Radius, Sign, Chords, Sweep (signed, radians: the angular travel
    from the run's first to its last vertex about the centre), Parts, EdgeIndices}]; empty
    without the arc columns.  The sweep is a function of the two end vertices alone
    (mesh_spatial_coupon.jl run_sweep), so the part count does not depend on the chord count."""
    loops = {}
    for row in rows:
        loops.setdefault(int(row["Loop"]), []).append(row)
    result = []
    for _, loop_rows in sorted(loops.items()):
        loop_rows.sort(key=lambda row: int(row["Vertex"]))
        points = [(float(row["X"]), float(row["Y"])) for row in loop_rows]
        n = len(points)
        ids = [int(row["ArcId"]) if row.get("ArcId") not in (None, "") else 0 for row in loop_rows]
        runs = []
        if any(ids):
            start = next((i for i in range(n) if ids[i] and ids[i - 1] != ids[i]), None)
            if start is None:
                raise ValueError("a plan-view loop made of one closed arc is not supported by the tube recipe")
            seen = set()
            for step in range(n):
                index = (start + step) % n
                if not ids[index] or ids[index - 1] == ids[index]:
                    continue
                if ids[index] in seen:
                    raise ValueError(f"arc {ids[index]} of the plan-view boundary is not one consecutive run of chords")
                seen.add(ids[index])
                edges = []
                k = index
                while ids[k] == ids[index] and len(edges) < n:
                    edges.append(k)
                    k = (k + 1) % n
                row = loop_rows[index]
                centre = (float(row["ArcCx"]), float(row["ArcCy"]))
                first, last = points[edges[0]], points[(edges[-1] + 1) % n]
                steps = []
                for e in edges:
                    a, b = points[e], points[(e + 1) % n]
                    ra, rb = (a[0] - centre[0], a[1] - centre[1]), (b[0] - centre[0], b[1] - centre[1])
                    steps.append(math.atan2(ra[0] * rb[1] - ra[1] * rb[0], ra[0] * rb[0] + ra[1] * rb[1]))
                orientation = 1.0 if sum(steps) >= 0.0 else -1.0
                theta_first = math.atan2(first[1] - centre[1], first[0] - centre[0])
                theta_last = math.atan2(last[1] - centre[1], last[0] - centre[0])
                travel = (orientation * (theta_last - theta_first)) % (2.0 * math.pi)
                sweep = orientation * (2.0 * math.pi if first == last else travel)
                runs.append({"ArcId": ids[index], "Centre": centre, "Radius": float(row["ArcR"]),
                             "Sign": int(row["ArcSign"]), "Chords": len(edges), "Sweep": sweep,
                             "Parts": arc_part_count(sweep), "EdgeIndices": edges})
        result.append(runs)
    return result


def unsupported_scope_classes(classes):
    """The guarded classes among `classes`, in guard order."""
    return [guard for guard in RECIPE_SCOPE_GUARDS if guard in classes]


def read_csv_rows(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def scope_classes_of_case_inputs(signature, boundary, process, *, device_footprint, trace_basis,
                                 fabricated=True):
    """The classes of a manifest case from its frozen inputs (signature, boundary,
    process.toml values) and its kind (a thin case exhibits ThinMetal; decision 66)."""
    return scope_classes(read_csv_rows(signature), read_csv_rows(boundary), fabricated=fabricated,
                         sidewall_angle=float(process["SidewallAngle"]),
                         top_rounding=float(process["TopRounding"]),
                         trench_rounding=float(process["TrenchRounding"]),
                         overetch=float(process["Overetch"]),
                         device_footprint=device_footprint, trace_basis=trace_basis)


def scope_classes_of_build(build_report):
    """The classes of a gmsh-build stage's input from its bound inputs and command."""
    command = build_report["Command"]
    inputs = build_report["Inputs"]
    process = {option: _option_or_default(command, option, default)
               for option, default in SCOPE_PROCESS_OPTIONS.items()}
    return scope_classes(read_csv_rows(inputs["source-signature"]["Path"]),
                         read_csv_rows(inputs["source-boundary"]["Path"]),
                         fabricated=build_coupon_kind(command) == "fabricated",
                         sidewall_angle=process["--sidewall-angle"],
                         top_rounding=process["--top-radius"],
                         trench_rounding=process["--bottom-radius"],
                         overetch=process["--overetch"],
                         device_footprint="source-retained-etch" in inputs,
                         trace_basis=bound_trace_basis(build_report) is not None)


def metal_loop_sides(boundary_rows, lower, upper, tolerance):
    """Per plan-view loop (in Loop order) the straight sides not lying on an outer box
    face (mesh_spatial_coupon.jl metal_loop_records); each carries TubesPerSide tubes
    unless it is an untubed short side (validate_untubed_edges)."""
    return [len(sides) for sides in metal_loop_side_points(boundary_rows, lower, upper, tolerance)]


def metal_loop_side_points(boundary_rows, lower, upper, tolerance):
    """Per plan-view loop (in Loop order) the (start, stop) points of the straight sides
    not lying on an outer box face, in vertex order; the chords of a tagged arc are not
    straight sides (metal_loop_arc_parts counts the arc's parts)."""
    loops = {}
    for row in boundary_rows:
        loops.setdefault(int(row["Loop"]), []).append((int(row["Vertex"]), float(row["X"]), float(row["Y"])))
    arc_runs = boundary_arc_runs(boundary_rows)
    sides = []
    for loop_index, (_, vertices) in enumerate(sorted(loops.items())):
        points = [(x, y) for _, x, y in sorted(vertices)]
        arc_edges = {e for run in arc_runs[loop_index] for e in run["EdgeIndices"]}
        loop_sides = []
        for i, p in enumerate(points):
            if i in arc_edges:
                continue
            q = points[(i + 1) % len(points)]
            on_face = any((abs(p[d] - lower[d]) <= tolerance and abs(q[d] - lower[d]) <= tolerance) or
                          (abs(p[d] - upper[d]) <= tolerance and abs(q[d] - upper[d]) <= tolerance)
                          for d in range(2))
            if not on_face:
                loop_sides.append((p, q))
        sides.append(loop_sides)
    return sides


def metal_loop_arc_parts(boundary_rows):
    """Per plan-view loop (Loop order) the number of arc PARTS its tagged arcs carry (block (b)
    design 1.2 (2): every part is a side with TubesPerSide tubes); 0 without arcs."""
    return [sum(run["Parts"] for run in runs) for runs in boundary_arc_runs(boundary_rows)]


def validate_untubed_edges(scope, census, side_points, tolerance):
    """The untubed short sides of the census Scope (block (b) design A9 family 3, decision
    304): every record names a straight metal side of the bound plan-view boundary (its
    Start / Stop are consecutive loop vertices not on a box face), a positive Span equal
    to the side's length, two non-negative Clearances, Interval = Span - their sum below
    the tube's inner ring size (PrismTubes.InnerSize: an interval of at least one inner
    ring carries a tube), its two Corners flags and CoveredByBalls exactly when the
    corner balls (the census CornerIsotropyRadius about each corner end) reach over the
    span; no side is recorded twice.  A census without the key (a pre-rule census)
    records no untubed side.  Returns the record count."""
    rows = scope.get("UntubedEdges", [])
    if not isinstance(rows, list):
        raise ValueError("Build census Scope UntubedEdges is not a list")
    if rows and not (isinstance(scope.get("UntubedShortEdgeRule"), str) and scope["UntubedShortEdgeRule"]):
        raise ValueError("Build census Scope records untubed edges without their rule")
    inner = _census_number(census["PrismTubes"], "InnerSize", "Prism tube record")
    radius = _census_number(census, "CornerIsotropyRadius", "Build census")
    seen = set()
    for row in rows:
        side = row.get("Side") if isinstance(row, dict) else None
        if (not isinstance(side, dict) or
                any(not isinstance(side.get(name), list) or len(side[name]) != 2 or
                    any(isinstance(v, bool) or not isinstance(v, (int, float)) for v in side[name])
                    for name in ("Start", "Stop"))):
            raise ValueError("Build census untubed edge lacks its side")
        start, stop = tuple(side["Start"]), tuple(side["Stop"])
        matches = [(p, q) for loop in side_points for p, q in loop
                   if ((math.dist(p, start) <= tolerance and math.dist(q, stop) <= tolerance) or
                       (math.dist(p, stop) <= tolerance and math.dist(q, start) <= tolerance))]
        if len(matches) != 1 or matches[0] in seen:
            raise ValueError("Build census untubed edge is not a straight metal side of the bound plan-view boundary")
        seen.add(matches[0])
        span = _census_number(row, "Span", "Untubed edge")
        clearances = row.get("Clearances")
        corners = row.get("Corners")
        if (span <= 0.0 or abs(span - math.dist(start, stop)) > tolerance or
                not isinstance(clearances, list) or len(clearances) != 2 or
                any(isinstance(c, bool) or not isinstance(c, (int, float)) or c < 0.0 for c in clearances) or
                not isinstance(corners, list) or len(corners) != 2 or
                any(not isinstance(c, bool) for c in corners) or
                not isinstance(row.get("CoveredByBalls"), bool)):
            raise ValueError("Build census untubed edge lacks its span, clearances, corners or coverage")
        interval = span - clearances[0] - clearances[1]
        if abs(_census_number(row, "Interval", "Untubed edge") - interval) > tolerance or not interval < inner:
            raise ValueError("Build census untubed edge interval is not below the tube inner size")
        if row["CoveredByBalls"] != (radius * sum(corners) >= span):
            raise ValueError("Build census untubed edge coverage disagrees with its corner balls")
    return len(rows)


def validate_recipe_scope(build_report, census):
    """The census Scope block (decision 48) names the recipe, spells the supported and
    guarded class lists of this contract with a statement per guard, exhibits exactly
    the classes recomputed from the bound inputs and command (none of them guarded: a
    guarded class never reaches the census), and counts per loop the sides not on the
    census CouponBox so that TubeCount = TubesPerSide x their sum minus the untubed short
    sides (2 for a fabricated coupon, 1 for a thin one; decision 66; the untubed sides
    per validate_untubed_edges)."""
    scope = census.get("Scope")
    if not isinstance(scope, dict) or not isinstance(scope.get("Rule"), str) or not scope["Rule"]:
        raise ValueError("Build census lacks the recipe Scope record")
    if scope.get("Recipe") != RECIPE_SCOPE_RECIPE:
        raise ValueError("Build census Scope names another recipe")
    if (scope.get("SupportedClasses") != list(RECIPE_SCOPE_SUPPORTED_CLASSES) or
            scope.get("GuardedClasses") != list(RECIPE_SCOPE_GUARDS)):
        raise ValueError("Build census Scope classes differ from the recipe scope of this contract")
    guards = scope.get("Guards")
    if (not isinstance(guards, list) or
            [guard.get("Id") if isinstance(guard, dict) else None for guard in guards] != list(RECIPE_SCOPE_GUARDS) or
            any(guard.get("DetectedFrom") != RECIPE_SCOPE_GUARDS[guard["Id"]] or
                not isinstance(guard.get("Statement"), str) or not guard["Statement"] for guard in guards)):
        raise ValueError("Build census Scope guards lack their id, detection origin or statement")
    exhibited = scope.get("ExhibitedClasses")
    if (not isinstance(exhibited, list) or exhibited != sorted(set(exhibited)) or
            exhibited != scope_classes_of_build(build_report)):
        raise ValueError("Build census Scope exhibited classes differ from the bound inputs")
    if unsupported_scope_classes(exhibited):
        raise ValueError("Build census exhibits a class the recipe guards")
    loops = scope.get("MetalLoops")
    box = census.get("CouponBox") if isinstance(census.get("CouponBox"), dict) else {}
    radius = _census_number(box, "Radius", "Coupon box")
    boundary_rows = read_csv_rows(build_report["Inputs"]["source-boundary"]["Path"])
    side_points = metal_loop_side_points(boundary_rows, box.get("Lower"), box.get("Upper"), 1e-7 * radius)
    arc_parts = metal_loop_arc_parts(boundary_rows)
    expected = [len(sides) + parts for sides, parts in zip(side_points, arc_parts)]
    if (not isinstance(loops, list) or len(loops) != len(expected) or
            any(not isinstance(loop, dict) or not isinstance(loop.get("Hole"), bool) or
                _count(loop.get("Sides"), "Metal loop sides") != sides or
                _count(loop.get("Vertices"), "Metal loop vertices") < 3 or
                _count(loop.get("Loop"), "Metal loop index") != index + 1
                for index, (loop, sides) in enumerate(zip(loops, expected)))):
        raise ValueError("Build census Scope metal loops differ from the bound plan-view boundary")
    # Arc loops (design 1.2 (2)): the straight sides, the arc parts and the arcs themselves.
    for loop, sides, parts, runs in zip(loops, side_points, arc_parts, boundary_arc_runs(boundary_rows)):
        if not runs:
            if any(key in loop for key in ("StraightSides", "ArcParts", "Arcs")):
                raise ValueError("Build census Scope metal loop records arcs the bound boundary does not carry")
            continue
        arcs = loop.get("Arcs")
        if (_count(loop.get("StraightSides"), "Metal loop straight sides") != len(sides) or
                _count(loop.get("ArcParts"), "Metal loop arc parts") != parts or
                not isinstance(arcs, list) or len(arcs) != len(runs) or
                any(not isinstance(arc, dict) or arc.get("ArcId") != run["ArcId"] or
                    _count(arc.get("Parts"), "Metal loop arc parts") != run["Parts"] or
                    _count(arc.get("Chords"), "Metal loop arc chords") != run["Chords"] or
                    abs(_census_number(arc, "Radius", "Metal loop arc") - run["Radius"]) > 1e-9 * run["Radius"] or
                    abs(_census_number(arc, "SweepDegrees", "Metal loop arc") - math.degrees(run["Sweep"])) > 1e-9
                    for arc, run in zip(arcs, runs))):
            raise ValueError("Build census Scope metal loop arcs differ from the bound plan-view boundary")
    if ("ArcSides" in exhibited) != any(arc_parts):
        raise ValueError("Build census Scope arc sides differ from the exhibited classes")
    if ("HoleLoops" in exhibited) != any(loop["Hole"] for loop in loops):
        raise ValueError("Build census Scope hole loops differ from the exhibited classes")
    tubes = census.get("PrismTubes")
    kind = build_coupon_kind(build_report["Command"])
    section = tubes.get("Section") if isinstance(tubes, dict) else None
    if (not isinstance(section, dict) or section.get("Kind") != kind or
            section.get("TubesPerSide") != TUBES_PER_SIDE[kind] or
            ("ThinMetal" in exhibited) != (kind == "thin")):
        raise ValueError("Prism tube section does not name the command's coupon kind and its tubes per side")
    untubed = validate_untubed_edges(scope, census, side_points, 1e-7 * radius)
    tubed_sides = sum(expected) - untubed
    if tubes.get("TubeCount") != TUBES_PER_SIDE[kind] * tubed_sides or tubed_sides <= 0:
        raise ValueError("Prism tube count is not TubesPerSide x the straight sides of every loop minus the untubed short sides")
    cutoff = section.get("ThinCutoff")
    if kind == "thin":
        # Decision 66: the thin coupon is recorded at its cutoff (the inner ring size).
        if (isinstance(cutoff, bool) or not isinstance(cutoff, (int, float)) or
                cutoff != _option_or_default(build_report["Command"], "--edge-size", None)):
            raise ValueError("Thin coupon census does not record its cutoff as the tube inner size")
    elif cutoff is not None:
        raise ValueError("Fabricated coupon census records a thin cutoff")
    validate_tube_layers(build_report, tubes)
    return scope


def validate_tube_layers(build_report, tubes):
    """Every census tube row names its process layer sign (Tubes[].Layer = the signature
    Nz of the rows on its Plane) and lies on that layer's edge: the top tube at
    Plane + Layer x MetalThickness (the command's --metal-thickness), the bottom tube on
    the Plane (decision 48: b = (0, 0, Nz) is the tube frame); a thin coupon's sheet
    tube on the Plane (decision 66)."""
    command = build_report["Command"]
    thickness = _option_or_default(command, *GMSH_BUILD_THICKNESS_OPTION)
    edges = TUBE_EDGES[build_coupon_kind(command)]
    signs = {}
    for row in read_csv_rows(build_report["Inputs"]["source-signature"]["Path"]):
        signs.setdefault(float(row["Pz"]), set()).add(int(float(row.get("Nz", 1) or 1)))
    if any(len(values) != 1 for values in signs.values()):
        raise ValueError("Signature rows on one plane carry both process normals")
    rows = tubes.get("Tubes")
    if not isinstance(rows, list):
        raise ValueError("Prism tube record lacks its rows")
    for row in rows:
        plane = _census_number(row, "Plane", "Tube row")
        origin = row.get("Origin")
        layer = row.get("Layer")
        matching = [sign for z, values in signs.items() if abs(z - plane) <= 1e-9 * max(1.0, abs(plane))
                    for sign in values]
        if (isinstance(layer, bool) or layer not in (-1, 1) or matching != [layer] or
                not isinstance(origin, list) or len(origin) != 3 or
                any(isinstance(x, bool) or not isinstance(x, (int, float)) for x in origin) or
                row.get("Edge") not in edges):
            raise ValueError("Prism tube row lacks its process layer sign or lies on another layer")
        expected = plane + layer * thickness if row["Edge"] == "top" else plane
        if abs(origin[2] - expected) > 1e-9 * max(1.0, abs(expected)):
            raise ValueError("Prism tube row does not lie on its layer's metal edge (plane + Nz x thickness / plane)")
    return rows


def _census_number(record, name, description):
    value = record.get(name) if isinstance(record, dict) else None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{description} lacks a finite {name}")
    return float(value)


def validate_size_bounds(census, command):
    """The census SizeBounds record binds the coupon-scale size bound to the build
    command: RequestedTangentialSize = --lc-tangent, FarSize = --lc-far, TangentialSize
    = min(request, FarSize) with the flag saying whether the bound acted."""
    bounds = census.get("SizeBounds")
    if not isinstance(bounds, dict) or not isinstance(bounds.get("Rule"), str) or not bounds["Rule"]:
        raise ValueError("Build census lacks the SizeBounds record")
    requested = _option_or_default(command, GMSH_BUILD_TANGENTIAL_OPTION, None)
    far = _option_or_default(command, "--lc-far", None)
    if (_census_number(bounds, "RequestedTangentialSize", "Size bounds") != requested or
            _census_number(bounds, "FarSize", "Size bounds") != far or
            _census_number(bounds, "TangentialSize", "Size bounds") != min(requested, far) or
            bounds.get("TangentialSizeBoundByFarSize") is not (far < requested)):
        raise ValueError("Build census SizeBounds do not follow TangentialSize = min(--lc-tangent, FarSize)")
    return bounds


def validate_coupon_box(census):
    """The census CouponBox record (decision 47): the box rule text, the radius, a
    consistent box, and every CAD-subdivision chain with its rows, union length and the
    extension verdict (a chain is extended iff its union reaches Radius from its
    midpoint); ChainedRows / ExtendedChains are the recomputed counts."""
    box = census.get("CouponBox")
    if (not isinstance(box, dict) or not isinstance(box.get("Rule"), str) or
            "subdivision" not in box["Rule"] or not isinstance(box.get("EdgeChains"), list)):
        raise ValueError("Build census lacks the CouponBox rule record")
    radius = _census_number(box, "Radius", "Coupon box")
    lower, upper = box.get("Lower"), box.get("Upper")
    if (radius <= 0.0 or not isinstance(lower, list) or not isinstance(upper, list) or
            len(lower) != 3 or len(upper) != 3 or
            any(not isinstance(a, (int, float)) or not isinstance(b, (int, float)) or not a < b
                for a, b in zip(lower, upper))):
        raise ValueError("Build census CouponBox radius or bounds are inconsistent")
    rows, extended = 0, 0
    for chain in box["EdgeChains"]:
        chain_rows = chain.get("Rows") if isinstance(chain, dict) else None
        if (not isinstance(chain_rows, list) or len(chain_rows) < 2 or
                any(isinstance(row, bool) or not isinstance(row, int) or row <= 0 for row in chain_rows) or
                len(set(chain_rows)) != len(chain_rows)):
            raise ValueError("Build census CouponBox chain lacks two or more distinct rows")
        length = _census_number(chain, "UnionLength", "Coupon box chain")
        rows += len(chain_rows)
        extended += length / 2.0 >= radius - 1e-10 * radius
    if (_count(box.get("ChainedRows"), "Chained rows") != rows or
            _count(box.get("ExtendedChains"), "Extended chains") != extended):
        raise ValueError("Build census CouponBox chain counts do not follow the recorded chains")
    return box


FACE_END_FACES = ("x0", "x1", "y0", "y1")
FACE_END_RULE = ("block (b) design A2 (supervisor decisions 302 / 320): a tube end on a box face with a single metal side "
                 "at the vertex and a tilt theta > 0 ends ON the face with m = ceil(2 (R + h_pyr) |tan theta| / lc_end) "
                 "sheared layers of spacing lc_end = max(TangentialSize, 4 h_pyr |tan theta|) (regime I); mesher design "
                 "round 2 F2b 3.2 (decisions 358 / 363 / 437): above 4 h_pyr |tan theta| > lc_cap the spacing is capped at "
                 "lc_cap = FACE_END_CONDITION_MARGIN x MaximumJacobianCondition x the smallest planar singular value of "
                 "the section's prism corner frames and m = max(ceil((R + h_pyr) |tan theta| / (lc_cap - 2 h_pyr |tan theta|)), "
                 "the regime-I count at lc_cap) (regime II), so the thinnest layer keeps t_min >= 2 h_pyr |tan theta|; "
                 "2 h_pyr |tan theta| >= lc_cap fails closed at ScopeGuard[SteepFaceCrossing]; mesher design round 3 "
                 "9H (part M 3.3): an ARC tube's formulas read its crossing slope max_u |s'(u)| = rho h / (r sqrt(r^2 - "
                 "h^2)) at the inner node circle r = rho - (R + h_pyr) in place of |tan theta| (FaceEnds[].CrossingSlope; "
                 "tan theta at the axis), a derived ceiling INSIDE the admitted tilt range; round 3 8A-bitwise (part M "
                 "5.2): above the kind's largest built tilt FACE_END_DERIVED_APEX_ABOVE_DEGREES the block's pyramids take "
                 "h_pyr_end = min(h_pyr, TangentialSize / (4 s)) (FaceEnds[].PyramidHeight) in every formula, so the "
                 "sheared layers stay regular")
FACE_END_CONDITION_MARGIN = 0.95
FACE_END_REGIMES = ("I", "II")
GMSH_BUILD_CONDITION_OPTION = "--maximum-jacobian-condition"


def section_frame_singular_values(ring_radii, angles):
    """The planar corner frames of a tube section's prisms (mesher design round 2 F2b; the
    mesher's section_frame_singular_values): for every sector triangle of the section - the
    inner triangle (edge point, ring 1 at both rays) and, per outer ring k, the two triangles
    (r_{k-1} j, r_k j, r_k j+1) and (r_{k-1} j, r_k j+1, r_{k-1} j+1) - the two edges leaving
    each vertex as the columns of a 2 x 2 matrix; returns (sigma_1, sigma_2) per frame in
    closed form.  A prism layer of axial spacing lc adds the orthogonal column lc e, so its
    corner frame has the singular values {lc, sigma_1, sigma_2}."""
    radii = [0.0] + [float(r) for r in ring_radii]
    rays = [float(a) for a in angles]
    if len(radii) < 2 or len(rays) < 2:
        raise ValueError("Tube section frames need at least one ring and one sector")

    def node(k, j):
        angle = math.radians(rays[j])
        return (radii[k] * math.cos(angle), radii[k] * math.sin(angle))

    frames = []
    for j in range(len(rays) - 1):
        triangles = [(node(0, j), node(1, j), node(1, j + 1))]
        for k in range(2, len(radii)):
            a, b, c, d = node(k - 1, j), node(k, j), node(k, j + 1), node(k - 1, j + 1)
            triangles.append((a, b, c))
            triangles.append((a, c, d))
        for triangle in triangles:
            for v in range(3):
                p = triangle[v]
                e1 = (triangle[(v + 1) % 3][0] - p[0], triangle[(v + 1) % 3][1] - p[1])
                e2 = (triangle[(v + 2) % 3][0] - p[0], triangle[(v + 2) % 3][1] - p[1])
                n1 = e1[0] ** 2 + e1[1] ** 2
                n2 = e2[0] ** 2 + e2[1] ** 2
                dot = e1[0] * e2[0] + e1[1] * e2[1]
                half = 0.5 * (n1 + n2)
                root = 0.5 * math.sqrt((n1 - n2) ** 2 + 4.0 * dot ** 2)
                frames.append((math.sqrt(half + root), math.sqrt(max(half - root, 0.0))))
    return frames


def section_prism_condition(ring_radii, angles, spacing):
    """The largest Jacobian condition over the prism corner frames of one regular layer of
    axial spacing `spacing` (the mesher's section_prism_condition)."""
    if not spacing > 0.0:
        raise ValueError("A prism layer needs a positive spacing")
    return max(max(spacing, s1) / min(spacing, s2) for s1, s2 in section_frame_singular_values(ring_radii, angles))


def face_end_spacing_cap(ring_radii, angles, maximum_jacobian_condition):
    """The end-spacing cap lc_cap of a face end (design round 2 F2b 3.2): FACE_END_CONDITION_MARGIN
    x the Jacobian-condition ceiling x the smallest sigma_2 over the section's prism frames,
    valid only in the spacing-dominated regime (every sigma_1 at or below it)."""
    if (isinstance(maximum_jacobian_condition, bool) or not isinstance(maximum_jacobian_condition, (int, float)) or
            not math.isfinite(maximum_jacobian_condition) or not maximum_jacobian_condition > 1.0):
        raise ValueError("A face end's end-spacing cap needs a finite Jacobian-condition ceiling > 1")
    ceiling = FACE_END_CONDITION_MARGIN * maximum_jacobian_condition
    frames = section_frame_singular_values(ring_radii, angles)
    cap = ceiling * min(s2 for _, s2 in frames)
    if any(s1 > cap for s1, _ in frames):
        raise ValueError("The face-end spacing cap is not in the spacing-dominated regime of the prism frames")
    return cap


def section_frame_angles(section):
    """The rays of the section whose prism frames the face-end cap reads: the fabricated top
    section or the thin sheet section (both carry every ring of the coupon's tubes)."""
    for name in ("Top", "Sheet"):
        rays = section.get(name)
        if isinstance(rays, dict) and isinstance(rays.get("Angles"), list) and len(rays["Angles"]) >= 2:
            return rays["Angles"]
    raise ValueError("Tube section lacks the Top / Sheet rays")


def arc_crossing_slope(rho, envelope_radius, face_distance):
    """Mesher design round 3 class (9) 9H (part M 3.3 Fact 2): the largest crossing slope of an
    arc tube's end block over its section, |s'(u)| = rho h / (r sqrt(r^2 - h^2)) at the inner
    node circle r = rho - (Radius + PyramidHeight), h the distance from the arc centre to the
    face plane (tan theta at the axis; the mesher's `arc_crossing_slope`, spelled identically:
    sqrt only, no transcendental, so no deterministic_math route is needed)."""
    if not (rho > 0.0 and 0.0 <= envelope_radius < rho):
        raise ValueError("Tube face end crossing slope needs 0 <= envelope < rho")
    r = rho - envelope_radius
    if not r > face_distance >= 0.0:
        raise ValueError("Tube face end of an arc row: the inner node circle does not reach the face plane")
    return rho * face_distance / (r * math.sqrt((r - face_distance) * (r + face_distance)))


def validate_tube_face_ends(row, tangential_size, section, condition_ceiling=None, round2b=False,
                            inner_size=None, growth_ratio=None, box=None, fabricated=None):
    """The face-end records of one census tube row (Tubes[].FaceEnds, absent on a plain
    tube): every record names a box face and an end, its tilt lies in (0, 90) degrees,
    its spacing and layer count follow FACE_END_RULE from the section's radius and pyramid
    height, and at most one record per end.  A record carrying Regime (mesher design round
    2 F2b) is judged by its regime against the cap recomputed from the section's rings and
    rays and the command's Jacobian-condition ceiling (`condition_ceiling`): regime I below
    4 h_pyr s <= lc_cap with the A2 (4) formulas, regime II above with lc_end == lc_cap and
    the apex-rule layer count; every record binds EndSpacing <= lc_cap, EndSpacingCap ==
    lc_cap, ApexThickness == 2 h_pyr s and the apex inequality LayerThicknessRange[0] >=
    ApexThickness.  The slope s (round 3 9H) is |tan theta| for a straight row; for an ARC
    row it is the arc's crossing slope recomputed from the row's Arc (Centre, Radius), the
    row's envelope and the face plane of the census CouponBox (`box` = (Lower, Upper)), and
    the record's CrossingSlope must carry it (a straight row's CrossingSlope, when recorded,
    is |tan theta|).  The block's pyramid height (round 3 8A-bitwise) is h_pyr at or below the
    kind's FACE_END_DERIVED_APEX_ABOVE_DEGREES (`fabricated` names the kind) and h_pyr_end =
    min(h_pyr, TangentialSize / (4 s)) above it; a record above the threshold must carry
    PyramidHeight = h_pyr_end, and every recorded PyramidHeight is bound.  A record without
    Regime is a pre-F2b record (the regime-I formulas alone); from the round-2b mesher
    (`round2b`) it fails closed.  Returns the largest EndSpacing (0 without face ends)."""
    records = row.get("FaceEnds", [])
    if not isinstance(records, list):
        raise ValueError("Tube row FaceEnds is not a list")
    if not records:
        return 0.0
    if not isinstance(section, dict):
        raise ValueError("Prism tube section is missing")
    radius = _census_number(section, "Radius", "Tube section")
    pyramid_height = _census_number(section, "PyramidHeight", "Tube section")
    ring_radii = section.get("RingRadii")
    if "Rings" in row:
        # A side reduced by its facing width (design round 2 F6): its own tube's radius,
        # pyramid height and rings bound its face ends.
        rings = _count(row.get("Rings"), "Tube row rings")
        if inner_size is None or growth_ratio is None:
            raise ValueError("Tube row records per-side rings without the recipe's inner size and ratio")
        ring_radii = [inner_size * (growth_ratio**k - 1.0) / (growth_ratio - 1.0) for k in range(1, rings + 1)]
        radius = ring_radii[-1]
        pyramid_height = GMSH_BUILD_PYRAMID_HEIGHT_OVER_OUTER_RING * inner_size * growth_ratio**(rings - 1)
    bound = 0.0
    ends = []
    cap = None
    for record in records:
        if (not isinstance(record, dict) or record.get("Face") not in FACE_END_FACES or
                record.get("End") not in ("start", "end")):
            raise ValueError("Tube face end lacks its face or end")
        theta = _census_number(record, "ThetaDegrees", "Tube face end")
        if not 0.0 < theta < 90.0:
            raise ValueError("Tube face end tilt is outside (0, 90) degrees")
        slope = abs(math.tan(math.radians(theta)))
        if "Arc" in row:
            # Round 3 9H: the arc's own crossing slope, from the row's circle and the face plane.
            if box is None:
                raise ValueError("Tube face end of an arc row needs the census CouponBox to recompute its crossing slope")
            arc = row["Arc"]
            axis = 0 if record["Face"][0] == "x" else 1
            face_value = (box[0] if record["Face"][1] == "0" else box[1])[axis]
            slope = arc_crossing_slope(_census_number(arc, "Radius", "Tube arc"), radius + pyramid_height,
                                       abs(face_value - _census_number({"c": arc["Centre"][axis]}, "c", "Tube arc centre")))
            if "CrossingSlope" not in record:
                raise ValueError("Tube face end of an arc row lacks its CrossingSlope")
        if "CrossingSlope" in record and \
                abs(_census_number(record, "CrossingSlope", "Tube face end") - slope) > 1e-12 * slope:
            raise ValueError("Tube face end CrossingSlope does not follow the row's crossing geometry")
        # Round 3 8A-bitwise: the block's own pyramid height above the kind's largest built tilt.
        apex_height = pyramid_height
        if fabricated is not None:
            threshold = FACE_END_DERIVED_APEX_ABOVE_DEGREES["fabricated" if fabricated else "thin"]
            if theta > threshold * (1.0 + 1e-9):
                apex_height = min(pyramid_height, tangential_size / (4.0 * slope))
        if "PyramidHeight" in record:
            if abs(_census_number(record, "PyramidHeight", "Tube face end") - apex_height) > 1e-12 * apex_height:
                raise ValueError("Tube face end PyramidHeight does not follow the kind's derived apex rule")
        elif apex_height != pyramid_height:
            raise ValueError("Tube face end above the kind's largest built tilt lacks its PyramidHeight")
        spacing = _census_number(record, "EndSpacing", "Tube face end")
        layers = _count(record.get("Layers"), "Tube face end layers")
        regime_one_spacing = max(tangential_size, 4.0 * apex_height * slope)

        def regime_one_layers(lc):
            return max(1, math.ceil(2.0 * (radius + pyramid_height) * slope / lc * (1.0 - 1e-9)))

        if "Regime" in record:
            if cap is None:
                if condition_ceiling is None:
                    raise ValueError("Tube face end records a regime without the command's Jacobian-condition ceiling")
                if not isinstance(ring_radii, list) or not ring_radii:
                    raise ValueError("Tube section lacks its ring radii")
                cap = face_end_spacing_cap(ring_radii, section_frame_angles(section), condition_ceiling)
                # The section's cap is the coupon section's (its RingRadii); a reduced row's own
                # cap is bound per record below.
                recorded_cap = _census_number(section, "FaceEndSpacingCap", "Tube section")
                coupon_cap = face_end_spacing_cap(section.get("RingRadii"), section_frame_angles(section), condition_ceiling)
                if abs(recorded_cap - coupon_cap) > 1e-12 * coupon_cap or not (
                        isinstance(section.get("FaceEndSpacingCapRule"), str) and section["FaceEndSpacingCapRule"]):
                    raise ValueError("Tube section FaceEndSpacingCap does not follow the section's prism frames")
            regime = record["Regime"]
            apex = 2.0 * apex_height * slope
            if regime not in FACE_END_REGIMES:
                raise ValueError("Tube face end regime is unknown")
            if abs(_census_number(record, "EndSpacingCap", "Tube face end") - cap) > 1e-12 * cap or \
                    abs(_census_number(record, "ApexThickness", "Tube face end") - apex) > 1e-12 * apex:
                raise ValueError("Tube face end cap or apex thickness does not follow the section")
            if 4.0 * apex_height * slope <= cap:
                expected_regime, expected_spacing = "I", regime_one_spacing
                expected_layers = regime_one_layers(expected_spacing)
            else:
                if not apex < cap:
                    raise ValueError("Tube face end lies beyond the validity ceiling of the capped end block")
                expected_regime, expected_spacing = "II", cap
                expected_layers = max(math.ceil((radius + pyramid_height) * slope / (cap - apex)),
                                      regime_one_layers(cap))
            if regime != expected_regime:
                raise ValueError("Tube face end regime does not follow 4 h_pyr |tan theta| against the cap")
            # lc_end <= lc_cap (design 3.2); a regime-I block never exceeds the regular layers
            # either, so TangentialSize above the cap (never a production size: 50 nm against
            # 80.6 nm fabricated) stays the bitwise A2 (4) value.
            if spacing > max(cap, tangential_size) * (1.0 + 1e-9):
                raise ValueError("Tube face end spacing exceeds the end-spacing cap")
            thickness = record.get("LayerThicknessRange")
            if (not isinstance(thickness, list) or len(thickness) != 2 or
                    not thickness[0] >= apex * (1.0 - 1e-9)):
                raise ValueError("Tube face end thinnest layer violates the apex rule")
        else:
            if round2b:
                raise ValueError("Tube face end from the round-2b mesher lacks its Regime record")
            expected_spacing = regime_one_spacing
            expected_layers = regime_one_layers(expected_spacing)
        if abs(spacing - expected_spacing) > 1e-9 * expected_spacing or layers != expected_layers:
            raise ValueError("Tube face end spacing or layer count does not follow the face-end rule")
        thickness = record.get("LayerThicknessRange")
        if (not isinstance(thickness, list) or len(thickness) != 2 or
                not 0.5 * spacing * (1.0 - 1e-9) <= thickness[0] <= thickness[1] <= 1.5 * spacing * (1.0 + 1e-9)):
            raise ValueError("Tube face end layer thickness range is outside [lc_end / 2, 3 lc_end / 2]")
        ends.append(record["End"])
        bound = max(bound, spacing)
    if len(set(ends)) != len(ends):
        raise ValueError("Tube row has two face ends at one end")
    return bound


def validate_tube_face_end_summary(tubes, rows):
    """PrismTubes.FaceEnds: Count = the face-end records over the rows, EndBlockLayers their
    layers, LegacyBoxVertexCorners the recorded theta-0 Physical-class box vertices (a
    non-negative count), the rules named."""
    summary = tubes.get("FaceEnds")
    records = [record for row in rows for record in row.get("FaceEnds", [])]
    legacy = sum(len(row.get("LegacyBoxVertexCorner", [])) for row in rows)
    if summary is None and not records and legacy == 0:
        return          # a census recorded before the face-end rule (no face end, no legacy corner)
    spacing = max([r["EndSpacing"] for r in records], default=0.0)
    if _census_number(tubes, "FaceEndSpacingMaximum", "Prism tube record") != spacing:
        raise ValueError("Prism tube face-end summary does not match the tube rows")
    if (not isinstance(summary, dict) or
            _count(summary.get("Count"), "Face-end count") != len(records) or
            _count(summary.get("EndBlockLayers"), "Face-end block layers") != sum(r["Layers"] for r in records) or
            _count(summary.get("LegacyBoxVertexCorners"), "Legacy box-vertex corners") != legacy or
            not all(isinstance(summary.get(name), str) and summary[name] for name in ("Rule", "BoxVertexRule"))):
        raise ValueError("Prism tube face-end summary does not match the tube rows")


ARC_CURVATURE_BOUND = 0.25


def validate_arc_tubes(tubes, rows, boundary_rows):
    """The arc tubes of a census (block (b) design 1.2 (3) / A3 (2)): PrismTubes.ArcTubes
    {Count, JointEnds, PartSplits, SharedSections = JointEnds + PartSplits, TotalArcLength, Rule,
    SmoothJointRule} against the rows
    carrying an Arc record (ArcId / Centre / Radius / Sign of a tagged run of the bound boundary,
    Part in 1..Parts = the run's parts, SweepDegrees = the run's sweep over its parts, Length =
    Radius x sweep of the tube's interval <= the part's), the section's envelope within
    ARC_CURVATURE_BOUND x Radius, and the rows' Joints records (TiltRadians >= 0, PlaneCut iff
    tilt > 0); a census without arc rows carries no ArcTubes record (None).  Returns the arc row count."""
    arc_rows = [row for row in rows if isinstance(row, dict) and "Arc" in row]
    summary = tubes.get("ArcTubes")
    if not arc_rows:
        if summary is not None:
            raise ValueError("Prism tube record names arc tubes without arc rows")
        if any("Joints" in row for row in rows if isinstance(row, dict)):
            raise ValueError("Prism tube rows record joints without arc tubes")
        return 0
    runs = {run["ArcId"]: run for loop_runs in boundary_arc_runs(boundary_rows) for run in loop_runs}
    section = tubes.get("Section") if isinstance(tubes, dict) else None
    envelope = (_census_number(section, "Radius", "Tube section") +
                _census_number(section, "PyramidHeight", "Tube section"))
    joint_ends = 0
    for row in rows:
        for joint in row.get("Joints", []):
            tilt = _census_number(joint, "TiltRadians", "Tube joint")
            if (tilt < 0.0 or joint.get("End") not in ("start", "end") or
                    joint.get("PlaneCut") is not (tilt > 0.0)):
                raise ValueError("Prism tube joint record lacks its end, tilt or plane-cut flag")
            joint_ends += 1
    for row in arc_rows:
        arc = row["Arc"]
        run = runs.get(arc.get("ArcId")) if isinstance(arc, dict) else None
        if run is None:
            raise ValueError("Prism tube arc row names an arc the bound boundary does not carry")
        radius = _census_number(arc, "Radius", "Tube arc")
        parts = _count(arc.get("Parts"), "Tube arc parts")
        part = _count(arc.get("Part"), "Tube arc part")
        centre = arc.get("Centre")
        if (not isinstance(centre, list) or len(centre) != 2 or
                any(abs(float(c) - r) > 1e-9 * max(1.0, radius) for c, r in zip(centre, run["Centre"])) or
                abs(radius - run["Radius"]) > 1e-9 * run["Radius"] or arc.get("Sign") != run["Sign"] or
                parts != run["Parts"] or not 1 <= part <= parts or
                abs(_census_number(arc, "SweepDegrees", "Tube arc") - math.degrees(abs(run["Sweep"])) / parts) > 1e-9 or
                envelope > ARC_CURVATURE_BOUND * radius * (1.0 + 1e-12) or
                _census_number(row, "Length", "Tube row") > radius * abs(run["Sweep"]) / parts * (1.0 + 1e-9)):
            raise ValueError("Prism tube arc row does not follow its tagged arc (centre, radius, sign, parts, sweep)")
    part_splits = sum(1 for row in arc_rows if _count(row["Arc"].get("Part"), "Tube arc part") <
                      _count(row["Arc"].get("Parts"), "Tube arc parts"))
    if (not isinstance(summary, dict) or
            _count(summary.get("Count"), "Arc tube count") != len(arc_rows) or
            _count(summary.get("JointEnds"), "Arc tube joint ends") != joint_ends or
            _count(summary.get("PartSplits"), "Arc tube part splits") != part_splits or
            _count(summary.get("SharedSections"), "Arc tube shared sections") != joint_ends + part_splits or
            abs(_census_number(summary, "TotalArcLength", "Arc tubes") -
                sum(float(row["Length"]) for row in arc_rows)) > 1e-9 * max(1.0, sum(float(row["Length"]) for row in arc_rows)) or
            not all(isinstance(summary.get(name), str) and summary[name] for name in ("Rule", "SmoothJointRule"))):
        raise ValueError("Prism tube arc summary does not match the arc rows")
    return len(arc_rows)


def validate_gmsh_build_census(build_report, census, semantic):
    """The Gmsh-only build census (the build report of decision 38) is bound to the
    build command and the canonical semantic contract: the corner ball and its
    grading to the tube inner size, the etch footprint provenance and simplified
    polygons, the junction segments, the trace basis sizing (iff bound), positive
    per-label interface areas over exactly the contract labels, the prism tube
    record (sizes equal to the command options, one tube per recorded row, spacing
    within the tangential size, the layers following the size field on the axis
    with their thickness statistics within the growth ratio, geometric rings, the
    sector angle equal to
    --tube-sector-degrees or its default, the pyramid height 0.5 x the outermost
    ring, the band law's RadialGrowth 1 and ProtectedDistance 2 x NormalSize, the
    band curves 1D-meshed at NormalSize with the achieved junction first-layer
    transverse sizes recorded, the decision-39 volume laws growing with FarGrowth
    with their achieved shell statistics), the per-type quality computed by
    the mesher within the command's gates for every volume type (orientation and
    Jacobian condition; scaled Jacobian for the tetrahedra), the cap regions within
    the scaled-Jacobian gate, the corner aspects within the corner gate and the
    element count within the command's cap.  No tetrahedral edge layer is recorded."""
    command = build_report["Command"]
    if _option_value(command, GMSH_BUILD_TUBE_OPTION).lower() != "true":
        raise ValueError(f"Gmsh-only build command does not pass {GMSH_BUILD_TUBE_OPTION} true")
    if (not isinstance(census, dict) or census.get("Version") != 1 or
            not isinstance(census.get("Corners"), list) or
            not isinstance(census.get("LongitudinalFaces"), list)):
        raise ValueError("Build census has an unsupported schema")
    corners = census.get("SemanticCorners")
    corner_free = contract_is_corner_free(semantic)
    if (not isinstance(corners, list) or (not corners and not corner_free) or len(census["Corners"]) != len(corners) or
            sorted(corners) != sorted(semantic.get("SemanticCorners", []))):
        raise ValueError("Build census corners differ from the canonical semantic contract")
    corner_gates = census.get(CORNER_GATES_RECORD)
    if corner_gates is not None and corner_gates != CORNER_GATES_NOT_APPLICABLE:
        raise ValueError(f"Build census {CORNER_GATES_RECORD} record is not {CORNER_GATES_NOT_APPLICABLE!r}")
    if corner_free != (corner_gates is not None):
        raise ValueError(f"Build census {CORNER_GATES_RECORD} record differs from the contract's Derivation.CornerFree")
    normal = _option_or_default(command, "--lc-fine", None)
    radius = _option_or_default(command, "--corner-isotropy-radius", None)
    if (_census_number(census, "IsotropicSize", "Build census") != normal or
            _census_number(census, "CornerIsotropyRadius", "Build census") != radius):
        raise ValueError("Build census size or corner radius differs from the build command")
    etch = build_report["Inputs"].get("source-retained-etch")
    if etch is None:
        if census.get("EtchBoundary") != PRODUCER_DEFAULT_FOOTPRINT_PROVENANCE:
            raise ValueError("Build census names an etch footprint the stage did not bind")
    elif census.get("EtchBoundarySHA256") != etch["SHA256"]:
        raise ValueError("Build census etch footprint differs from the bound retained etch")
    areas = census.get("InterfaceAreas")
    if (not isinstance(areas, list) or not areas or
            any(not isinstance(row, dict) or isinstance(row.get("Attribute"), bool) or
                not isinstance(row.get("Attribute"), int) or
                _census_number(row, "Area", "Interface area row") <= 0 for row in areas)):
        raise ValueError("Build census lacks positive per-label interface areas")
    labels = [row["Attribute"] for row in areas]
    if len(labels) != len(set(labels)) or set(labels) != boundary_attributes(semantic):
        raise ValueError("Build census interface-area labels differ from the semantic contract")
    validate_footprint_polygons(census, build_coupon_kind(command))
    gmsh_build_junction_segments(census)
    validate_build_trace_basis(build_report, census)
    if census.get("EdgeLayer") is not None:
        raise ValueError("Gmsh-only build census records a tetrahedral edge layer")
    grading = census.get("CornerGrading")
    corner_size = _option_or_default(command, CORNER_SIZE_OPTION, 0.0)
    ratio = _option_or_default(command, "--edge-growth-ratio", EDGE_LAYER_DEFAULTS["GrowthRatio"])
    if (not isinstance(grading, dict) or corner_size <= 0.0 or not corner_size < normal or
            ratio <= 1.0 or _census_number(grading, "CornerSize", "Corner grading") != corner_size or
            _census_number(grading, "GrowthRatio", "Corner grading") != ratio or
            _census_number(grading, "NormalSize", "Corner grading") != normal or
            _census_number(grading, "Radius", "Corner grading") != radius):
        raise ValueError("Build census corner grading differs from the build command")
    expected = []
    size = corner_size
    while size < normal:
        expected.append(expected[-1] + size if expected else size)
        size *= ratio
    reach = (normal - corner_size) / (ratio - 1.0)
    if (grading.get("ShellRadii") != expected + [radius] or
            grading.get("ShellSizes") != [corner_size * ratio**k for k in range(len(expected))] + [normal] or
            _census_number(grading, "Reach", "Corner grading") != reach or not reach <= radius):
        raise ValueError("Build census corner grading shells do not follow CornerSize and GrowthRatio")
    tubes = census.get("PrismTubes")
    if not isinstance(tubes, dict):
        raise ValueError("Build census lacks the prism tube record")
    for option, name in GMSH_BUILD_RECIPE_OPTIONS.items():
        if _census_number(tubes, name, "Prism tube record") != _option_or_default(command, option, None):
            raise ValueError(f"Prism tube record {name} differs from the build command {option}")
    validate_size_bounds(census, command)
    validate_coupon_box(census)
    validate_recipe_scope(build_report, census)
    if tubes["TangentialSize"] != _census_number(census["SizeBounds"], "TangentialSize", "Size bounds"):
        raise ValueError("Prism tube record TangentialSize differs from the bound tangential size")
    if tubes["InnerSize"] != corner_size:
        raise ValueError("Prism tube inner size differs from the corner size: one graded law is required")
    rows = tubes.get("Tubes")
    if not isinstance(rows, list) or not rows or tubes.get("TubeCount") != len(rows):
        raise ValueError("Prism tube rows are missing or exceed the tangential spacing")
    # Face ends (block (b) design A2, supervisor decisions 302 / 320): a tube ending on a
    # box face at a tilt theta > 0 carries an end block of m sheared layers of spacing
    # lc_end = max(TangentialSize, 4 h_pyr |tan theta|), m = ceil(2 (R + h_pyr) |tan theta| /
    # lc_end); its Spacing may exceed TangentialSize by exactly that block.
    section = tubes.get("Section")
    face_end_bound = {}
    condition_ceiling = (_option_or_default(command, GMSH_BUILD_CONDITION_OPTION, None)
                         if GMSH_BUILD_CONDITION_OPTION in command else None)
    round2b = _mesher_digest(build_report) == sha256(ROUND2_MESHER)
    coupon_box = census.get("CouponBox")
    face_box = ((coupon_box["Lower"], coupon_box["Upper"])
                if isinstance(coupon_box, dict) and "Lower" in coupon_box and "Upper" in coupon_box else None)
    kind_fabricated = build_coupon_kind(command) == "fabricated"
    if isinstance(section, dict) and "FaceEndDerivedApexAboveDegrees" in section and \
            _census_number(section, "FaceEndDerivedApexAboveDegrees", "Tube section") != \
            FACE_END_DERIVED_APEX_ABOVE_DEGREES["fabricated" if kind_fabricated else "thin"]:
        raise ValueError("Tube section FaceEndDerivedApexAboveDegrees differs from the kind's largest built tilt")
    for index, row in enumerate(rows):
        face_end_bound[index] = validate_tube_face_ends(row, tubes["TangentialSize"], section,
                                                        condition_ceiling, round2b,
                                                        tubes["InnerSize"], tubes["GrowthRatio"],
                                                        box=face_box, fabricated=kind_fabricated)
    if (isinstance(section, dict) and section.get("FaceEndSpacingCap") is not None and
            not any(row.get("FaceEnds") for row in rows)):
        raise ValueError("Tube section records a face-end spacing cap without a face end")
    # The interior layers are capped at TangentialSize x (1 - 1e-9) by the mesher; a face-end
    # block layer is lc_end up to the rounding of its stations (hence the 1e-9 slack on lc_end).
    if any(not isinstance(row, dict) or
           not 0.0 < _census_number(row, "Spacing", "Tube row") <=
           max(tubes["TangentialSize"], face_end_bound[index] * (1.0 + 1e-9)) or
           _count(row.get("Layers"), "Tube layers") <= 0 or
           _census_number(row, "Length", "Tube row") <= 0.0 for index, row in enumerate(rows)):
        raise ValueError("Prism tube rows are missing or exceed the tangential spacing")
    spacing_bound = max([tubes["TangentialSize"]] + [bound * (1.0 + 1e-9) for bound in face_end_bound.values()])
    validate_tube_face_end_summary(tubes, rows)
    validate_tube_rings_per_side(tubes, rows, command)
    validate_arc_tubes(tubes, rows, read_csv_rows(build_report["Inputs"]["source-boundary"]["Path"]))
    # Decision 40: the layers follow the composed size field on the tube axis. Every
    # tube records its layer thickness statistics (the largest layer is its Spacing,
    # the neighbour ratio within the growth ratio), the record names the layer rule
    # and the axis size law, and the growth cap is the ring growth ratio; the count
    # of layers below TangentialSize / GrowthRatio is a count within the layers.
    # Decision 41: at a tube end on the outer box (EndsOnBox) the end layer is the
    # surface value of the field, at most the prescribed size there.
    layer_statistics = tubes.get("LayerThickness")
    if (not all(isinstance(tubes.get(name), str) and tubes[name]
                for name in ("LayerRule", "TubeAxisSizeLaw")) or
            _census_number(tubes, "LayerGrowthCap", "Prism tube record") != ratio or
            not isinstance(layer_statistics, dict) or
            _census_number(layer_statistics, "MaximumNeighbourRatio", "Tube layer thickness") > ratio or
            not 0.0 < _census_number(layer_statistics, "Minimum", "Tube layer thickness") <=
            _census_number(layer_statistics, "P50", "Tube layer thickness") <=
            _census_number(layer_statistics, "Maximum", "Tube layer thickness") <= spacing_bound or
            _census_number(tubes, "SpacingMinimum", "Prism tube record") != layer_statistics["Minimum"] or
            _census_number(tubes, "SpacingMaximum", "Prism tube record") != layer_statistics["Maximum"] or
            _count(layer_statistics.get("LayersBelowTangentialSizeOverGrowthRatio"),
                   "Tube layers below TangentialSize / GrowthRatio") > _count(tubes.get("Layers"), "Tube layers")):
        raise ValueError("Prism tube layers do not record the size-field layer rule within the growth ratio")
    for row in rows:
        layers = row.get("LayerThickness")
        ends_on_box = row.get("EndsOnBox")
        # A face end replaces the decision-41 surface layer of its end by its block: the
        # end layer IS the block spacing (the block / interior ratio is recorded, not gated).
        face_ends = {record["End"]: record for record in row.get("FaceEnds", [])}
        blocks = layers.get("FaceEndBlocks", {}) if isinstance(layers, dict) else {}
        if (not isinstance(layers, dict) or not isinstance(blocks, dict) or
                set(blocks) != {end.capitalize() for end in face_ends} or
                any(_census_number(layers, name, "Tube row layer thickness") <= 0.0
                    for name in ("Minimum", "P50", "Maximum", "AtStart", "AtEnd",
                                 "PrescribedAtStart", "PrescribedAtEnd")) or
                layers["Maximum"] != row["Spacing"] or
                _census_number(layers, "MaximumNeighbourRatio", "Tube row layer thickness") > ratio or
                not isinstance(layers.get("AchievedOverPrescribed"), dict) or
                any(_census_number(layers["AchievedOverPrescribed"], name, "Tube row achieved layers") <= 0.0
                    for name in ("Minimum", "P50", "Maximum")) or
                not isinstance(ends_on_box, list) or len(ends_on_box) != 2 or
                not all(isinstance(flag, bool) for flag in ends_on_box) or
                ("start" in face_ends and not ends_on_box[0]) or ("end" in face_ends and not ends_on_box[1]) or
                ("start" not in face_ends and ends_on_box[0] and
                 layers["AtStart"] > layers["PrescribedAtStart"] * (1.0 + 1e-9)) or
                ("end" not in face_ends and ends_on_box[1] and
                 layers["AtEnd"] > layers["PrescribedAtEnd"] * (1.0 + 1e-9)) or
                ("start" in face_ends and (abs(layers["AtStart"] - face_ends["start"]["EndSpacing"]) >
                                           1e-9 * face_ends["start"]["EndSpacing"] or
                                           _count(blocks.get("Start", {}).get("Layers"), "Face-end block layers") !=
                                           face_ends["start"]["Layers"] or
                                           _census_number(blocks.get("Start", {}), "NeighbourRatio",
                                                          "Face-end block") <= 0.0)) or
                ("end" in face_ends and (abs(layers["AtEnd"] - face_ends["end"]["EndSpacing"]) >
                                         1e-9 * face_ends["end"]["EndSpacing"] or
                                         _count(blocks.get("End", {}).get("Layers"), "Face-end block layers") !=
                                         face_ends["end"]["Layers"] or
                                         _census_number(blocks.get("End", {}), "NeighbourRatio",
                                                        "Face-end block") <= 0.0))):
            raise ValueError("Prism tube row lacks its layer thickness record within the growth ratio "
                             "with the surface layer at an end on the box (or its face-end block)")
    section = tubes.get("Section")
    rings = _count(section.get("Rings") if isinstance(section, dict) else None, "Tube rings")
    sizes = section.get("RingSizes") if isinstance(section, dict) else None
    if (rings <= 0 or not isinstance(sizes, list) or len(sizes) != rings or
            any(abs(size - tubes["InnerSize"] * ratio**k) > 1e-12 * tubes["InnerSize"] * ratio**k
                for k, size in enumerate(sizes))):
        raise ValueError("Prism tube rings do not follow the inner size and growth ratio")
    if _census_number(section, "SectorDegrees", "Tube section") != _option_or_default(
            command, *GMSH_BUILD_SECTOR_OPTION):
        raise ValueError(f"Prism tube sector angle differs from the build command {GMSH_BUILD_SECTOR_OPTION[0]}")
    if (_census_number(section, "PyramidHeightOverOuterRing", "Tube section") !=
            GMSH_BUILD_PYRAMID_HEIGHT_OVER_OUTER_RING or
            abs(_census_number(section, "PyramidHeight", "Tube section") -
                GMSH_BUILD_PYRAMID_HEIGHT_OVER_OUTER_RING * sizes[-1]) > 1e-12 * sizes[-1]):
        raise ValueError("Prism tube pyramid height is not 0.5 x the outermost ring size")
    if (_count(tubes.get("Prisms"), "Tube prisms") <= 0 or
            _count(tubes.get("Pyramids"), "Tube pyramids") <= 0):
        raise ValueError("Prism tube record has no prisms or pyramids")
    laws = tubes.get("SizeLaws")
    if (not isinstance(laws, dict) or
            any(_census_number(laws, name, "Size laws") != tubes[name]
                for name in ("NormalSize", "FarSize", "FarGrowth")) or
            not all(isinstance(laws.get(name), str) and laws[name]
                    for name in ("TubeRule", "BandRule", "Composition"))):
        raise ValueError("Prism tube size laws are not recorded")
    if (_census_number(laws, "RadialGrowth", "Size laws") != GMSH_BUILD_BAND_RADIAL_GROWTH or
            _census_number(laws, "ProtectedDistance", "Size laws") !=
            GMSH_BUILD_BAND_PROTECTED_DISTANCE_OVER_NORMAL * tubes["NormalSize"]):
        raise ValueError("Prism tube band law growth or protected distance differs from the recipe")
    # Decision 39 volume laws: the trace rule and the corner-ball exterior grow with
    # FarGrowth (the trace record's GradingSlope is that growth), the corner radius
    # is the command's, the junction lines carry the band law in the volume, and the
    # achieved shell statistics are recorded (trace apexes iff a basis is bound).
    trace_record = census.get("TraceBasisSizing")
    achieved = laws.get("Achieved")
    if (_census_number(laws, "CornerExteriorGrowth", "Size laws") != tubes["FarGrowth"] or
            _census_number(laws, "TraceBasisVolumeGrowth", "Size laws") != tubes["FarGrowth"] or
            _census_number(laws, "CornerIsotropyRadius", "Size laws") != radius or
            not all(isinstance(laws.get(name), str) and laws[name]
                    for name in ("JunctionVolumeRule", "CornerExteriorRule", "TraceBasisVolumeRule")) or
            (trace_record is not None and
             _census_number(trace_record, "GradingSlope", "Trace basis sizing") != tubes["FarGrowth"])):
        raise ValueError("Prism tube volume size laws (trace, corner exterior, junction) are not bound to FarGrowth")
    if (not isinstance(achieved, dict) or
            any(not isinstance(achieved.get(name), dict) or
                not isinstance(achieved[name].get("Shells"), list) or not achieved[name]["Shells"]
                for name in ("CornerExterior", "JunctionLines")) or
            (trace_record is None) != (achieved.get("TraceApexes") is None) or
            (trace_record is not None and
             (not isinstance(achieved["TraceApexes"], dict) or not achieved["TraceApexes"].get("Shells")))):
        raise ValueError("Build census does not record the achieved volume sizes of the size laws")
    band_curves = tubes.get("BandCurves")
    if (not isinstance(band_curves, dict) or
            _census_number(band_curves, "Spacing", "Band curves") != tubes["NormalSize"] or
            _count(band_curves.get("Count"), "Band curves") < 0):
        raise ValueError("Band curves are not 1D-meshed at NormalSize")
    validate_curve_spacing(census, tubes["NormalSize"], tubes["TangentialSize"], ratio,
                           band_curves["Count"])
    bands = tubes.get("Bands")
    first_layer = bands.get("JunctionFirstLayer") if isinstance(bands, dict) else None
    if (not isinstance(first_layer, dict) or
            _census_number(bands, "Prescribed", "Bands") != tubes["NormalSize"] or
            any(not isinstance(first_layer.get(name), dict) or
                _count(first_layer[name].get("Elements"), f"Junction first layer {name}") <= 0 or
                any(_census_number(first_layer[name], statistic, f"Junction first layer {name}") <= 0
                    for statistic in ("TransverseP50", "TransverseP90", "AchievedOverPrescribedP50"))
                for name in ("CutSurface", "Tetrahedra"))):
        raise ValueError("Build census does not record the achieved junction first-layer sizes")
    gates = {name: _option_or_default(command, option, None)
             for option, name in GMSH_BUILD_GATE_OPTIONS.items()}
    if any(not math.isfinite(value) or value <= 0 for value in gates.values()):
        raise ValueError("Gmsh-only build command lacks the quality gates")
    quality = tubes.get("Quality")
    if not isinstance(quality, dict):
        raise ValueError("Prism tube record lacks the per-type quality")
    total = 0
    for name in GMSH_BUILD_VOLUME_TYPES:
        record = quality.get(name)
        if not isinstance(record, dict):
            raise ValueError(f"Build census quality lacks the {name} record")
        total += _count(record.get("Count"), f"{name} count")
        if (record.get("PositiveOrientation") is not True or
                _count(record.get("NonpositiveCells"), f"{name} nonpositive cells") != 0 or
                _census_number(record, "MaximumJacobianCondition", name) >
                gates["MaximumJacobianCondition"]):
            raise ValueError(f"Build census {name} cells fail orientation or the condition gate")
    if _census_number(quality["Tetrahedron"], "MinimumScaledJacobian", "Tetrahedra") < gates["MinimumScaledJacobian"]:
        raise ValueError("Build census tetrahedra fall below the scaled-Jacobian gate")
    if quality.get("Total") != total or total <= 0:
        raise ValueError("Build census quality total differs from its per-type counts")
    cap = _option_or_default(command, GMSH_BUILD_ELEMENT_CAP_OPTION, None)
    budget = tubes.get("FarFieldBudgetPolicy")
    if (not isinstance(budget, dict) or budget.get("Elements") != total or
            budget.get("MaximumElements") != cap or total > cap or
            _census_number(budget, "EffectiveFarSize", "Far-field budget") != tubes["FarSize"]):
        raise ValueError("Build census element budget differs from the build command cap")
    caps = tubes.get("CapRegions")
    if not isinstance(caps, dict) or (_count(caps.get("Caps"), "Cap regions") == 0) != corner_free:
        raise ValueError("Build census tube cap regions differ from the contract's corners (a cap centre before "
                         "every corner; none on a corner-free coupon)")
    if corner_free:
        if (caps.get("MinimumScaledJacobian") is not None or caps.get("MaximumJacobianCondition") is not None or
                caps.get("Regions") != []):
            raise ValueError("Build census records tube cap region statistics on a corner-free coupon")
    elif (_census_number(caps, "MinimumScaledJacobian", "Cap regions") < gates["MinimumScaledJacobian"] or
            _census_number(caps, "MaximumJacobianCondition", "Cap regions") > gates["MaximumJacobianCondition"]):
        raise ValueError("Build census tube cap regions fail the tetrahedral gates")
    optimization = census.get("SeedQualityOptimization")
    if (not isinstance(optimization, dict) or
            any(_census_number(optimization, name, "Seed quality optimization") != value
                for name, value in gates.items()) or
            _count(optimization.get("RequiredCellsBelowGateAfter"), "Corner cells below the gate") != 0 or
            _count(optimization.get("RequiredCellsAboveConditionAfter"), "Corner cells above the condition gate") != 0 or
            not isinstance(optimization.get("CornerAspectsAfter"), list) or
            len(optimization["CornerAspectsAfter"]) != len(corners) or
            any(isinstance(value, bool) or not isinstance(value, (int, float)) or
                not math.isfinite(value) for value in optimization["CornerAspectsAfter"])):
        raise ValueError("Build census does not record gated corner balls")
    # Design round 2 F5-A: the per-corner verdict by kind (legacy MaximumCornerAspect,
    # invariant CornerShapeGate with no bridging-sliver candidate above it), bound to the
    # contract's record; a pre-rule census (decision 392 MINOR-6) by the pre-rule rule, bound
    # to its declared tool.
    corner_shape_gate = (_option_or_default(command, GMSH_BUILD_CORNER_SHAPE_GATE_OPTION, None)
                         if GMSH_BUILD_CORNER_SHAPE_GATE_OPTION in command else None)
    if corner_shape_gate is not None and (not math.isfinite(corner_shape_gate) or corner_shape_gate <= 1.0):
        raise ValueError("Gmsh-only build command carries an invalid --corner-shape-gate")
    rule_round = census_rule_round(census, optimization)
    validate_corner_measures(optimization, corners, gates["MaximumCornerAspect"], corner_shape_gate,
                             invariant_corners(semantic), rule_round)
    if rule_round == "pre-rule":
        validate_pre_rule_census_tool(build_report)
    else:
        kinds = census.get("SemanticCornerKinds")
        if (not isinstance(kinds, list) or len(kinds) != len(corners) or
                [point for point, kind in zip(census["SemanticCorners"], kinds) if kind == "Invariant"] !=
                [row["Point"] for row in optimization["CornerMeasures"] if row["Kind"] == "Invariant"] or
                any(kind not in CORNER_MEASURES for kind in kinds)):
            raise ValueError("Build census SemanticCornerKinds differ from the seed corner measures")
    validate_thin_sheet_seams(census, build_coupon_kind(command), rule_round)
    return census


def validate_thin_sheet_seams(census, kind, rule_round="round-2"):
    """Design round 2 SEAM (supervisor decisions 363 / 368): a thin build records the
    pinched-seam census ThinSheetSeams {Count, Edges, ...}; its Count is 0 unless the record
    carries UnrefinedCrackSeams (the seams of the convex tips sharper than the bisector's
    minimum opening, every seam attributed: Count equal, RefineCrackElements false) - the
    solve-config writer then sets Model.RefineCrackElements false (case_inputs); a
    fabricated build carries no seam. A PRE-RULE census (rule_round 'pre-rule', decision 392
    MINOR-6: no round-2 record declared) carries no seam census and is accepted as such; a
    round-2 thin census lacking the record fails closed. Returns the UnrefinedCrackSeams
    record or None."""
    if rule_round not in CENSUS_RULE_ROUNDS:
        raise ValueError(f"Unknown census rule round {rule_round!r}")
    seams = census.get("ThinSheetSeams")
    if rule_round == "pre-rule":
        return None         # no round-2 record declared: no seam census to bind
    if kind != "thin":
        if seams is not None and _count(seams.get("Count"), "Thin sheet seams") != 0:
            raise ValueError("A fabricated build census records thin-sheet seams")
        return None
    if (not isinstance(seams, dict) or _count(seams.get("Count"), "Thin sheet seams") < 0 or
            not isinstance(seams.get("Edges"), list) or len(seams["Edges"]) != seams["Count"]):
        raise ValueError("Thin build census lacks the ThinSheetSeams record")
    unrefined = seams.get("UnrefinedCrackSeams")
    if unrefined is None:
        if seams["Count"] != 0:
            raise ValueError(f"Thin build census records {seams['Count']} pinched seams without UnrefinedCrackSeams")
        return None
    if (not isinstance(unrefined, dict) or unrefined.get("RefineCrackElements") is not False or
            _count(unrefined.get("Count"), "Unrefined crack seams") != seams["Count"] or
            not isinstance(unrefined.get("Tips"), list) or not unrefined["Tips"] or
            not isinstance(unrefined.get("Rule"), str) or "RefineCrackElements" not in unrefined["Rule"]):
        raise ValueError("Thin build census UnrefinedCrackSeams record is inconsistent with ThinSheetSeams")
    return unrefined


def validate_curve_spacing(census, normal, tangential, growth, band_count):
    """Decision 43: every explicitly 1D-meshed longitudinal curve records its spacing
    against the composed size field: Count rows (the band curves among them, at
    Spacing = NormalSize; metal parts at TangentialSize), each with positive
    node-spacing statistics ordered Minimum <= P50 <= Maximum <= Spacing, a positive
    PrescribedMinimum <= Spacing, AchievedOverPrescribed statistics (each node interval
    over the limited law at its midpoint) with the maximum within the growth cap - the
    equidistribution bounds an interval by the law over the interval, which steps down
    by at most GrowthRatio inside it (the corner-ball shells; every other law is
    Lipschitz) -, a kept-grid count within the grid intervals, and GradedCurves = the
    rows marked Graded; the growth cap is the recipe's GrowthRatio."""
    record = census.get("CurveSpacing")
    if (not isinstance(record, dict) or not isinstance(record.get("Rule"), str) or
            "composed" not in record["Rule"] or "gradient-limited" not in record["Rule"] or
            _census_number(record, "GrowthRatio", "Curve spacing") != growth or
            not isinstance(record.get("Curves"), list) or
            _count(record.get("Count"), "Curve spacing rows") != len(record["Curves"]) or
            _count(record.get("GradedCurves"), "Graded curves") !=
            sum(1 for row in record["Curves"] if isinstance(row, dict) and row.get("Graded") is True)):
        raise ValueError("Build census does not record the composed curve spacing rule")
    kinds = {"junction": 0, "band": 0, "metal": 0}
    for row in record["Curves"]:
        if not isinstance(row, dict) or row.get("Kind") not in kinds or not isinstance(row.get("Graded"), bool):
            raise ValueError("Curve spacing row is invalid")
        kinds[row["Kind"]] += 1
        spacing = _census_number(row, "Spacing", "Curve spacing row")
        expected = normal if row["Kind"] in ("junction", "band") else tangential
        nodes = row.get("NodeSpacing")
        achieved = row.get("AchievedOverPrescribed")
        if (spacing != expected or _census_number(row, "Length", "Curve spacing row") <= 0.0 or
                not isinstance(nodes, dict) or not isinstance(achieved, dict) or
                not 0.0 < _census_number(nodes, "Minimum", "Curve node spacing") <=
                _census_number(nodes, "P50", "Curve node spacing") <=
                _census_number(nodes, "Maximum", "Curve node spacing") <= spacing * (1.0 + 1e-9) or
                not 0.0 < _census_number(row, "PrescribedMinimum", "Curve spacing row") <= spacing * (1.0 + 1e-9) or
                not 0.0 < _census_number(achieved, "Minimum", "Curve achieved spacing") <=
                _census_number(achieved, "P50", "Curve achieved spacing") <=
                _census_number(achieved, "Maximum", "Curve achieved spacing") <= growth * (1.0 + 1e-9) or
                _count(row.get("GridIntervalsKept"), "Kept grid intervals") >
                _count(row.get("GridIntervals"), "Grid intervals") or
                _count(row.get("InteriorNodes"), "Curve interior nodes") < 0 or
                (not row["Graded"] and row["GridIntervalsKept"] != row["GridIntervals"])):
            raise ValueError("Curve spacing row does not follow the composed size field within its spacing")
    if kinds["junction"] + kinds["band"] != band_count:
        raise ValueError("Curve spacing rows do not cover the band curves")
    return record


def gmsh_build_junction_segments(census):
    """The census junction curves as straight segments [x0 y0 z0 x1 y1 z1]: Count
    finite nondegenerate segments whose total length is the recorded one and none
    curved (the Gmsh-only build requires sharp vertical geometry)."""
    curves = census.get("JunctionCurves")
    if (not isinstance(curves, dict) or not isinstance(curves.get("Segments"), list) or
            _count(curves.get("Count"), "Junction curves") <= 0 or
            len(curves["Segments"]) != curves["Count"] or
            _count(curves.get("CurvedCurves"), "Curved junction curves") != 0):
        raise ValueError("Build census lacks straight junction segments")
    total = 0.0
    for segment in curves["Segments"]:
        if (not isinstance(segment, list) or len(segment) != 6 or
                any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v)
                    for v in segment)):
            raise ValueError("Build census junction segment is invalid")
        length = math.dist(segment[:3], segment[3:])
        if length <= 0:
            raise ValueError("Build census junction segment is degenerate")
        total += length
    if abs(_census_number(curves, "TotalLength", "Junction curves") - total) > COPLANAR_TOLERANCE * total:
        raise ValueError("Build census junction length differs from its segments")
    return curves["Segments"]


def validate_build_trace_basis(build_report, census):
    """The build bound a trace basis (or none); when bound the command passes the
    dimensionless ratio and the census TraceBasisSizing records it with the bound
    digests and nonempty mesh-frame triangles."""
    basis = bound_trace_basis(build_report)
    ratio = _optional_ratio(build_report["Command"], basis is not None)
    record = census.get("TraceBasisSizing")
    if basis is None:
        if record is not None:
            raise ValueError("Build census records a trace basis the stage did not bind")
        return None
    triangles = record.get("MeshFrameTriangles") if isinstance(record, dict) else None
    if (not isinstance(record, dict) or _census_number(record, "Ratio", "Trace basis sizing") != ratio or
            record.get("RatioIsDimensionless") is not True or record.get("InputSHA256") != basis or
            not isinstance(triangles, list) or not triangles or record.get("Triangles") != len(triangles) or
            any(not isinstance(t, list) or len(t) != 3 or
                any(not isinstance(v, list) or len(v) != 3 or
                    any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x)
                        for x in v) for v in t) for t in triangles)):
        raise ValueError("Build census trace basis sizing differs from the bound trace basis")
    validate_trace_basis_size_measure(record, ratio)
    return record


def trace_basis_altitudes_and_shortest_edges(triangles):
    """(minimum altitude, shortest edge) of every mesh-frame basis triangle: the
    altitude is 2 x area / longest edge, the hat-gradient scale of decision 43."""
    rows = []
    for triangle in triangles:
        a, b, c = (tuple(float(x) for x in v) for v in triangle)
        lengths = [math.dist(a, b), math.dist(b, c), math.dist(c, a)]
        ab = [b[i] - a[i] for i in range(3)]
        ac = [c[i] - a[i] for i in range(3)]
        cross = [ab[1] * ac[2] - ab[2] * ac[1], ab[2] * ac[0] - ab[0] * ac[2], ab[0] * ac[1] - ab[1] * ac[0]]
        doubled_area = math.sqrt(sum(x * x for x in cross))
        if doubled_area <= 0.0 or min(lengths) <= 0.0:
            raise ValueError("Build census trace basis triangle is degenerate")
        rows.append((doubled_area / max(lengths), min(lengths)))
    return rows


def validate_trace_basis_size_measure(record, ratio):
    """Decision 43: the trace rule's size measure is the basis triangle's minimum
    altitude (named in SizeMeasure); MinimumBasisAltitude is the minimum over the
    recorded mesh-frame triangles, MinimumRequestedSize = Ratio x it, and the
    report-only needle counts (altitude < NeedleAltitudeOverShortestEdge x shortest
    edge; all, and those below FarSize) are the recomputed ones."""
    rows = trace_basis_altitudes_and_shortest_edges(record["MeshFrameTriangles"])
    minimum_altitude = min(altitude for altitude, _ in rows)
    fraction = _census_number(record, "NeedleAltitudeOverShortestEdge", "Trace basis sizing")
    far = _census_number(record, "FarSize", "Trace basis sizing")
    needles = [altitude < fraction * shortest for altitude, shortest in rows]
    narrow_needles = [needle and ratio * altitude < far for needle, (altitude, _) in zip(needles, rows)]
    if (fraction != TRACE_BASIS_NEEDLE_ALTITUDE_OVER_SHORTEST_EDGE or
            not isinstance(record.get("SizeMeasure"), str) or
            "minimum altitude" not in record["SizeMeasure"] or
            not isinstance(record.get("NeedleRule"), str) or "report-only" not in record["NeedleRule"] or
            abs(_census_number(record, "MinimumBasisAltitude", "Trace basis sizing") - minimum_altitude) >
            1e-12 * minimum_altitude or
            abs(_census_number(record, "MinimumRequestedSize", "Trace basis sizing") - ratio * minimum_altitude) >
            1e-12 * ratio * minimum_altitude or
            _count(record.get("NeedleTriangles"), "Needle triangles") != sum(needles) or
            _count(record.get("NeedleTrianglesBelowFarSize"), "Needle triangles below far") != sum(narrow_needles)):
        raise ValueError("Build census trace basis size measure is not the minimum altitude of the "
                         "recorded basis triangles (decision 43)")
    return record


def trace_basis_edges_of_census(census):
    """Unique edges of the census mesh-frame basis triangles (source-local), or []."""
    record = census.get("TraceBasisSizing")
    if not isinstance(record, dict):
        return []
    edges = {}
    for triangle in record["MeshFrameTriangles"]:
        for a, b in ((0, 1), (1, 2), (2, 0)):
            key = tuple(sorted((tuple(float(x) for x in triangle[a]), tuple(float(x) for x in triangle[b]))))
            edges[key] = None
    return [[*key[0], *key[1]] for key in edges]


def validate_gmsh_only_dag(reports, canonical_mesh):
    source = reports["canonical-source-validation"]
    build = reports["gmsh-build"]
    publication = reports["canonical-gmsh-publication"]
    links = (
        (build["Inputs"]["canonical-semantic-contract"]["SHA256"],
         source["Artifacts"]["canonical-semantic-contract"]["SHA256"]),
        (publication["Inputs"]["gmsh-mesh"]["SHA256"], build["Artifacts"]["gmsh-mesh"]["SHA256"]),
        (publication["Artifacts"]["canonical-candidate-mesh"]["SHA256"], sha256(canonical_mesh)),
    )
    if any(actual != expected for actual, expected in links):
        raise ValueError("Canonical mesh stage input/output digest chain is broken")
    census = json.loads(Path(build["Artifacts"]["build-census"]["Path"]).read_text())
    semantic = json.loads(Path(build["Inputs"]["canonical-semantic-contract"]["Path"]).read_text())
    validate_gmsh_build_census(build, census, semantic)
    return census


def validate_canonical_dag(report_paths, canonical_mesh, launcher_name=None,
                           launcher_sha256=None, expected_tool_sha256=None):
    pipeline = pipeline_of(report_paths, canonical_only=True)
    reports, digests = _validate_reports(report_paths, canonical_stage_order(pipeline),
                                         launcher_name, launcher_sha256,
                                         expected_tool_sha256, pipeline)
    if pipeline == GMSH_ONLY_PIPELINE:
        validate_gmsh_only_dag(reports, canonical_mesh)
        return reports, digests
    source = reports["canonical-source-validation"]
    seed_stage = reports["seed-generation"]
    seed = seed_stage["Artifacts"]["seed-mesh"]["SHA256"]
    metric = reports["metric-preparation"]
    adaptation = reports["native-adaptation-mmg"]
    restoration = reports["label-restoration"]
    publication = reports["canonical-gmsh-publication"]
    links = (
        (seed_stage["Inputs"]["canonical-semantic-contract"]["SHA256"],
         source["Artifacts"]["canonical-semantic-contract"]["SHA256"]),
        (metric["Inputs"]["seed-mesh"]["SHA256"], seed),
        (metric["Inputs"]["seed-corner-census"]["SHA256"],
         seed_stage["Artifacts"]["seed-corner-census"]["SHA256"]),
        (metric["Inputs"]["canonical-semantic-contract"]["SHA256"],
         source["Artifacts"]["canonical-semantic-contract"]["SHA256"]),
        (metric["Inputs"]["canonical-supports"]["SHA256"],
         source["Artifacts"]["canonical-supports"]["SHA256"]),
        *[(adaptation["Inputs"][name]["SHA256"], metric["Artifacts"][name]["SHA256"])
          for name in ("mmg-seed", "metric", "pins", "fixed-triangles",
                       "required-tetrahedra", "restoration-recipe")],
        (restoration["Inputs"]["adapted-mesh"]["SHA256"],
         adaptation["Artifacts"]["adapted-mesh"]["SHA256"]),
        (restoration["Inputs"]["restoration-recipe"]["SHA256"],
         metric["Artifacts"]["restoration-recipe"]["SHA256"]),
        (publication["Inputs"]["restored-mesh"]["SHA256"],
         restoration["Artifacts"]["restored-mesh"]["SHA256"]),
        (publication["Artifacts"]["canonical-candidate-mesh"]["SHA256"],
         sha256(canonical_mesh)),
    )
    if any(actual != expected for actual, expected in links):
        raise ValueError("Canonical mesh stage input/output digest chain is broken")
    validate_seed_corner_isotropy(seed_stage, metric["Artifacts"]["restoration-recipe"]["Path"])
    recipe = json.loads(Path(metric["Artifacts"]["restoration-recipe"]["Path"]).read_text())
    census = json.loads(Path(seed_stage["Artifacts"]["seed-corner-census"]["Path"]).read_text())
    validate_trace_basis_sizing(seed_stage, metric, recipe, census)
    validate_edge_layer(seed_stage, metric, adaptation, recipe, census)
    validate_corner_grading(seed_stage, metric, recipe, census)
    validate_required_region(seed_stage, restoration, recipe, census,
                             metric["Artifacts"]["required-tetrahedra"]["Path"])
    return reports, digests


def validate_placement_dag(report_paths, final_mesh, canonical_record,
                           launcher_name=None, launcher_sha256=None,
                           expected_tool_sha256=None):
    reports, digests = _validate_reports(report_paths, PLACEMENT_STAGE_ORDER,
                                         launcher_name, launcher_sha256,
                                         expected_tool_sha256, LEGACY_MMG_PIPELINE)
    publication = reports["proper-rigid-publication"]
    record = json.loads(Path(canonical_record).read_text())
    if (publication["Inputs"]["canonical-build-record"]["SHA256"] !=
            sha256(canonical_record) or
            publication["Inputs"]["canonical-candidate-mesh"]["SHA256"] !=
            record.get("CanonicalArtifacts", {}).get("canonical-candidate-mesh", {}).get("SHA256") or
            publication["Artifacts"]["candidate-mesh"]["SHA256"] != sha256(final_mesh)):
        raise ValueError("Placement publication is not bound to its canonical build and final mesh")
    receipt = json.loads(Path(publication["Artifacts"]["transform-receipt"]["Path"]).read_text())
    if (receipt.get("CanonicalBuildId") != record.get("CanonicalBuildId") or
            receipt.get("CanonicalBuildSHA256") != record.get("CanonicalBuildSHA256") or
            receipt.get("OutputMeshSHA256") != sha256(final_mesh) or
            receipt.get("TransformSHA256") !=
            publication["Inputs"]["placement-transform"]["SHA256"]):
        raise ValueError("Rigid publication receipt differs from placement-stage bindings")
    return reports, digests


def validate_canonical_record_binding(record, canonical_reports, canonical_report_paths):
    """The canonical build record must bind exactly the pipeline's stage outputs and reports."""
    pipeline = pipeline_of(canonical_report_paths, canonical_only=True)
    order = canonical_stage_order(pipeline)
    artifacts = record.get("CanonicalArtifacts")
    report_hashes = record.get("CanonicalStageReportSHA256")
    if (not isinstance(artifacts, dict) or not isinstance(report_hashes, dict) or
            set(report_hashes) != set(order) or
            set(artifacts) != canonical_artifact_roles(pipeline)):
        raise ValueError("Canonical build record artifact/report roles differ from the stage DAG")
    for stage in order:
        if report_hashes[stage] != sha256(canonical_report_paths[stage]):
            raise ValueError(f"Canonical build record binds a different stage report: {stage}")
        for role, item in canonical_reports[stage]["Artifacts"].items():
            if (artifacts[role].get("SHA256") != item["SHA256"] or
                    str(Path(artifacts[role].get("Path", "")).resolve()) !=
                    str(Path(item["Path"]).resolve())):
                raise ValueError(f"Canonical build record artifact differs from stage output: {role}")


def validate_stage_dag(report_paths, final_mesh, launcher_name=None, launcher_sha256=None,
                       expected_tool_sha256=None, canonical_record=None):
    """Validate the complete canonical + mandatory placement DAG of either pipeline
    (identified by the exact stage set of `report_paths`)."""
    pipeline = pipeline_of(report_paths)
    order = canonical_stage_order(pipeline)
    canonical_paths = {name: report_paths[name] for name in order}
    placement_paths = {name: report_paths[name] for name in PLACEMENT_STAGE_ORDER}
    placement_report = json.loads(Path(placement_paths["proper-rigid-publication"]).read_text())
    record_path = canonical_record or placement_report["Inputs"]["canonical-build-record"]["Path"]
    record = json.loads(Path(record_path).read_text())
    canonical_mesh = record["CanonicalArtifacts"]["canonical-candidate-mesh"]["Path"]
    canonical_expected = None if expected_tool_sha256 is None else {
        name: expected_tool_sha256[name] for name in order}
    placement_expected = None if expected_tool_sha256 is None else {
        name: expected_tool_sha256[name] for name in PLACEMENT_STAGE_ORDER}
    canonical, canonical_digests = validate_canonical_dag(
        canonical_paths, canonical_mesh, launcher_name, launcher_sha256, canonical_expected)
    validate_canonical_record_binding(record, canonical, canonical_paths)
    placement, placement_digests = validate_placement_dag(
        placement_paths, final_mesh, record_path, launcher_name, launcher_sha256,
        placement_expected)
    return {**canonical, **placement}, canonical_digests | placement_digests


def dag_pipeline(reports):
    """The pipeline of a validated report dict (canonical + placement stages)."""
    return pipeline_of(reports)


def build_stage_report(reports):
    """The volume-build stage report of a validated DAG (gmsh-build or seed-generation)."""
    return reports[PIPELINE_BUILD_STAGE[dag_pipeline(reports)]]
