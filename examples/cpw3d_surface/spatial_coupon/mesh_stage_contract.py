#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Validation for reusable canonical-build and per-placement publication DAGs."""
import hashlib
import json
import math
from pathlib import Path

from edge_volume_metric import COPLANAR_TOLERANCE
from semantic_mesh_contract import (boundary_attributes, cut_surface_attributes,
                                    material_interface_attributes)


CANONICAL_STAGE_ORDER = (
    "canonical-source-validation", "seed-generation", "metric-preparation",
    "native-adaptation-mmg", "label-restoration", "canonical-gmsh-publication")
PLACEMENT_STAGE_ORDER = ("proper-rigid-publication",)
STAGE_ORDER = CANONICAL_STAGE_ORDER + PLACEMENT_STAGE_ORDER
STAGE_TOOLS = {
    "canonical-source-validation": {"runtime", "source-validator"},
    "seed-generation": {"runtime", "mesher"},
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
    "metric-preparation": {"seed-mesh", "seed-corner-census", "canonical-semantic-contract",
                           "canonical-supports"},
    "native-adaptation-mmg": {"mmg-seed", "metric", "pins", "fixed-triangles",
                              "required-tetrahedra", "restoration-recipe"},
    "label-restoration": {"adapted-mesh", "restoration-recipe"},
    "canonical-gmsh-publication": {"restored-mesh", "source-process", "source-signature",
                                   "source-boundary"},
    "proper-rigid-publication": {
        "canonical-candidate-mesh", "canonical-build-record", "placement-transform",
        "source-semantic-contract", "source-signature", "source-boundary", "source-mask",
        "source-process"},
}
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
                         "metric-preparation": set(TRACE_BASIS_INPUTS)}
STAGE_OUTPUTS = {
    "canonical-source-validation": {"canonical-semantic-contract", "canonical-supports"},
    "seed-generation": {"seed-mesh", "seed-corner-census"},
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
    "metric-preparation": {option: ("Inputs", name) for option, name in TRACE_BASIS_OPTIONS.items()},
}
# Bound inputs/outputs that are consumed positionally and must occur in argv.
STAGE_BINDING_ARGUMENTS = {
    "canonical-source-validation": {
        ("Inputs", "canonical-transform"), ("Artifacts", "canonical-semantic-contract"),
        ("Artifacts", "canonical-supports")},
    "seed-generation": {("Inputs", "source-signature"), ("Artifacts", "seed-mesh")},
    "metric-preparation": {("Inputs", "seed-mesh")},
    "native-adaptation-mmg": {
        ("Inputs", "mmg-seed"), ("Inputs", "metric"), ("Inputs", "pins"),
        ("Inputs", "restoration-recipe"), ("Artifacts", "adapted-mesh"),
        ("Artifacts", "adaptation-receipt")},
    "label-restoration": {("Inputs", "adapted-mesh"), ("Inputs", "restoration-recipe"),
                          ("Artifacts", "restored-mesh")},
    "canonical-gmsh-publication": {("Inputs", "restored-mesh"),
                                   ("Artifacts", "canonical-candidate-mesh")},
    "proper-rigid-publication": {
        ("Inputs", "canonical-candidate-mesh"), ("Inputs", "placement-transform"),
        ("Artifacts", "candidate-mesh"), ("Artifacts", "transform-receipt")},
}
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


def validate_command_bindings(stage, command, inputs, artifacts, working_directory=None):
    """Require every consumed source/input/output path to be the bound one."""
    bound = {"Inputs": inputs, "Artifacts": artifacts}
    for option, (section, name) in STAGE_BINDING_OPTIONS.get(stage, {}).items():
        _require_path_option(command, option, bound[section][name], working_directory)
    for option, (section, name) in STAGE_OPTIONAL_BINDING_OPTIONS.get(stage, {}).items():
        if name in bound[section]:
            _require_path_option(command, option, bound[section][name], working_directory)
        elif option in command:
            raise ValueError(f"Stage command passes {option} without a bound {name} input")
    for section, name in STAGE_BINDING_ARGUMENTS.get(stage, set()):
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
                          expected_tool_sha256=None):
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
    _validate_bindings(report.get("Inputs"), STAGE_INPUTS[stage], "stage input",
                       STAGE_OPTIONAL_INPUTS.get(stage, frozenset()))
    _validate_bindings(report.get("Artifacts"), STAGE_OUTPUTS[stage], "stage output")
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
        report.get("WorkingDirectory"))
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
                      expected_tool_sha256):
    if set(report_paths) != set(order):
        raise ValueError("Stage report set differs from the required DAG role")
    reports, report_digests = {}, set()
    for stage in order:
        path = Path(report_paths[stage])
        expected = None if expected_tool_sha256 is None else expected_tool_sha256[stage]
        report = validate_stage_report(json.loads(path.read_text()), stage, launcher_name,
                                       launcher_sha256, expected)
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
    if (not isinstance(balls, dict) or not isinstance(corners, list) or not corners or
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


def footprint_segments(census):
    """Every edge of every simplified footprint polygon as a 3D segment
    [x0, y0, z, x1, y1, z] on the polygon's process plane, in census order."""
    segments = []
    for polygon in validate_footprint_polygons(census):
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


def validate_footprint_polygons(census):
    """The seed simplified every etch footprint polygon (device or producer default)
    with the shared collinearity tolerance before CAD face creation and recorded the
    result: removed vertices and a maximum deviation within tolerance x local scale."""
    tolerance = census.get("FootprintCollinearTolerance")
    if isinstance(tolerance, bool) or tolerance != COPLANAR_TOLERANCE:
        raise ValueError("Seed census footprint collinearity tolerance differs from COPLANAR_TOLERANCE")
    polygons = census.get("FootprintPolygons")
    summary = census.get("FootprintSimplification")
    if not isinstance(polygons, list) or not polygons or not isinstance(summary, dict):
        raise ValueError("Seed census lacks the simplified footprint polygons")
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
    if (not isinstance(corners, list) or not corners or
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


# Quality gates the seed stage enforces on the MMG required region (decision 30)
# and the label restorer enforces on the adapted mesh: both commands carry them,
# with equal values.
REQUIRED_REGION_GATE_OPTIONS = ("--maximum-corner-aspect", "--minimum-scaled-jacobian",
                                "--maximum-quality-displacement-over-normal")


def _required_indices(path, tetrahedra):
    values = Path(path).read_text().split()
    if not values or any(not value.isdigit() for value in values):
        raise ValueError("Required tetrahedron list must be non-empty positive integers")
    indices = [int(value) for value in values]
    if (indices != sorted(set(indices)) or indices[0] < 1 or indices[-1] > tetrahedra):
        raise ValueError("Required tetrahedron list must be sorted, unique and in range")
    return indices


def validate_required_region(seed_report, restoration_report, recipe, census, required_path):
    """The metric stage's required tetrahedra are the recipe's corner balls and edge
    layer, the list is the recorded one, and the seed stage optimized and gated the
    same region with the restorer's gate values.

    The recipe RequiredTetrahedra record names the rule, CornerIsotropyRadius, one
    count per semantic corner (each positive), the layer reach LayerThickness x
    (1 + RowZigzag) + EdgeSize and one count per span when a layer is recorded (no
    reach and no spans otherwise); the list has exactly Count sorted unique seed
    indices.  The seed command passes the three gate options with the restorer's
    values; the census SeedQualityOptimization record carries them, no required cell
    below the scaled-Jacobian gate and every corner aspect within the corner gate.
    """
    record = recipe.get("RequiredTetrahedra")
    corners = recipe.get("TruePhysicalCorners")
    if (not isinstance(record, dict) or not isinstance(corners, list) or not corners or
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
            _recipe_number(quality, "DisplacementBoundOverNormal") !=
            gates["--maximum-quality-displacement-over-normal"] or
            _count(quality.get("RequiredCellsBelowGateAfter"), "Required cells below the gate") != 0 or
            _count(quality.get("RequiredTetrahedra"), "Seed required tetrahedra") <= 0 or
            not isinstance(quality.get("CornerAspectsAfter"), list) or
            len(quality["CornerAspectsAfter"]) != len(corners) or
            any(isinstance(value, bool) or not isinstance(value, (int, float)) or
                not math.isfinite(value) or value > gates["--maximum-corner-aspect"]
                for value in quality["CornerAspectsAfter"])):
        raise ValueError("Seed census does not record a gated required-region optimization")
    if (_recipe_number(quality, "RequiredMinimumScaledJacobianAfter") <
            gates["--minimum-scaled-jacobian"]):
        raise ValueError("Seed required region is below the scaled-Jacobian gate")
    return record


TRACE_BASIS_RATIO_OPTION = "--trace-basis-size-ratio"


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


def validate_canonical_dag(report_paths, canonical_mesh, launcher_name=None,
                           launcher_sha256=None, expected_tool_sha256=None):
    reports, digests = _validate_reports(report_paths, CANONICAL_STAGE_ORDER,
                                         launcher_name, launcher_sha256,
                                         expected_tool_sha256)
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
    validate_required_region(seed_stage, restoration, recipe, census,
                             metric["Artifacts"]["required-tetrahedra"]["Path"])
    return reports, digests


def validate_placement_dag(report_paths, final_mesh, canonical_record,
                           launcher_name=None, launcher_sha256=None,
                           expected_tool_sha256=None):
    reports, digests = _validate_reports(report_paths, PLACEMENT_STAGE_ORDER,
                                         launcher_name, launcher_sha256,
                                         expected_tool_sha256)
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
    """The canonical build record must bind exactly the six stages' outputs and reports."""
    artifacts = record.get("CanonicalArtifacts")
    report_hashes = record.get("CanonicalStageReportSHA256")
    if (not isinstance(artifacts, dict) or not isinstance(report_hashes, dict) or
            set(report_hashes) != set(CANONICAL_STAGE_ORDER) or
            set(artifacts) != {role for stage in CANONICAL_STAGE_ORDER
                               for role in STAGE_OUTPUTS[stage]}):
        raise ValueError("Canonical build record artifact/report roles differ from the stage DAG")
    for stage in CANONICAL_STAGE_ORDER:
        if report_hashes[stage] != sha256(canonical_report_paths[stage]):
            raise ValueError(f"Canonical build record binds a different stage report: {stage}")
        for role, item in canonical_reports[stage]["Artifacts"].items():
            if (artifacts[role].get("SHA256") != item["SHA256"] or
                    str(Path(artifacts[role].get("Path", "")).resolve()) !=
                    str(Path(item["Path"]).resolve())):
                raise ValueError(f"Canonical build record artifact differs from stage output: {role}")


def validate_stage_dag(report_paths, final_mesh, launcher_name=None, launcher_sha256=None,
                       expected_tool_sha256=None, canonical_record=None):
    """Validate the complete canonical + mandatory placement DAG."""
    if set(report_paths) != set(STAGE_ORDER):
        raise ValueError("Exactly one report for every canonical and placement stage is required")
    canonical_paths = {name: report_paths[name] for name in CANONICAL_STAGE_ORDER}
    placement_paths = {name: report_paths[name] for name in PLACEMENT_STAGE_ORDER}
    placement_report = json.loads(Path(placement_paths["proper-rigid-publication"]).read_text())
    record_path = canonical_record or placement_report["Inputs"]["canonical-build-record"]["Path"]
    record = json.loads(Path(record_path).read_text())
    canonical_mesh = record["CanonicalArtifacts"]["canonical-candidate-mesh"]["Path"]
    canonical_expected = None if expected_tool_sha256 is None else {
        name: expected_tool_sha256[name] for name in CANONICAL_STAGE_ORDER}
    placement_expected = None if expected_tool_sha256 is None else {
        name: expected_tool_sha256[name] for name in PLACEMENT_STAGE_ORDER}
    canonical, canonical_digests = validate_canonical_dag(
        canonical_paths, canonical_mesh, launcher_name, launcher_sha256, canonical_expected)
    validate_canonical_record_binding(record, canonical, canonical_paths)
    placement, placement_digests = validate_placement_dag(
        placement_paths, final_mesh, record_path, launcher_name, launcher_sha256,
        placement_expected)
    return {**canonical, **placement}, canonical_digests | placement_digests
