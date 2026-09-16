#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Validation for reusable canonical-build and per-placement publication DAGs."""
import hashlib
import json
from pathlib import Path


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
    "metric-preparation": {"seed-mesh", "canonical-semantic-contract",
                           "canonical-supports"},
    "native-adaptation-mmg": {"mmg-seed", "metric", "pins", "fixed-triangles",
                              "restoration-recipe"},
    "label-restoration": {"adapted-mesh", "restoration-recipe"},
    "canonical-gmsh-publication": {"restored-mesh", "source-process", "source-signature",
                                   "source-boundary"},
    "proper-rigid-publication": {
        "canonical-candidate-mesh", "canonical-build-record", "placement-transform",
        "source-semantic-contract", "source-signature", "source-boundary", "source-mask",
        "source-process"},
}
STAGE_OUTPUTS = {
    "canonical-source-validation": {"canonical-semantic-contract", "canonical-supports"},
    "seed-generation": {"seed-mesh", "seed-corner-census"},
    "metric-preparation": {"mmg-seed", "metric", "pins", "fixed-triangles",
                           "restoration-recipe"},
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
        "--transformed-supports": ("Inputs", "canonical-supports")},
    "native-adaptation-mmg": {"--fixed-triangles": ("Inputs", "fixed-triangles")},
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


def _validate_bindings(items, expected, description):
    if not isinstance(items, dict) or set(items) != expected:
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
    _validate_bindings(report.get("Inputs"), STAGE_INPUTS[stage], "stage input")
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


def validate_seed_corner_isotropy(seed_report, recipe_path):
    """The seed's corner ball must be the recipe's: same size, radius and corners."""
    recipe = json.loads(Path(recipe_path).read_text())
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
        (metric["Inputs"]["canonical-semantic-contract"]["SHA256"],
         source["Artifacts"]["canonical-semantic-contract"]["SHA256"]),
        (metric["Inputs"]["canonical-supports"]["SHA256"],
         source["Artifacts"]["canonical-supports"]["SHA256"]),
        *[(adaptation["Inputs"][name]["SHA256"], metric["Artifacts"][name]["SHA256"])
          for name in ("mmg-seed", "metric", "pins", "fixed-triangles",
                       "restoration-recipe")],
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
