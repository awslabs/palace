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
    "native-adaptation-mmg": {"runtime", "adaptation-wrapper", "adapter-mmg"},
    "label-restoration": {"runtime", "label-restorer"},
    "canonical-gmsh-publication": {"runtime", "publisher"},
    "proper-rigid-publication": {"runtime", "rigid-publisher", "ownership-runtime",
                                 "ownership-auditor"},
}
STAGE_INPUTS = {
    "canonical-source-validation": {"source-semantic-contract", "source-signature",
                                    "source-boundary", "source-mask", "canonical-transform"},
    "seed-generation": set(),
    "metric-preparation": {"seed-mesh", "canonical-semantic-contract",
                           "canonical-supports"},
    "native-adaptation-mmg": {"mmg-seed", "metric", "pins", "fixed-triangles",
                              "restoration-recipe"},
    "label-restoration": {"adapted-mesh", "restoration-recipe"},
    "canonical-gmsh-publication": {"restored-mesh"},
    "proper-rigid-publication": {
        "canonical-candidate-mesh", "canonical-build-record", "placement-transform",
        "source-semantic-contract", "source-signature", "source-boundary", "source-mask",
        "source-process"},
}
STAGE_OUTPUTS = {
    "canonical-source-validation": {"canonical-semantic-contract", "canonical-supports"},
    "seed-generation": {"seed-mesh"},
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


def validate_tool_invocation(stage, command, tools, working_directory=None):
    base = Path(working_directory or Path.cwd())
    runtime = str(Path(tools["runtime"]).resolve())
    primary = str(Path(tools[STAGE_PRIMARY_TOOL[stage]]).resolve())
    def matches(argument, expected):
        path = Path(argument); resolved = path if path.is_absolute() else base / path
        return (str(resolved.resolve()) == expected or
                (not path.is_absolute() and expected.endswith("/" + argument)))
    if not matches(command[0], runtime):
        raise ValueError(f"Stage runtime is not argv[0]: {stage}")
    script_position = _first_interpreter_script(command, primary)
    if script_position is None or not matches(command[script_position], primary):
        raise ValueError(f"Stage tool is not in the executed script position: {stage}")
    if stage == "native-adaptation-mmg":
        _require_path_option(command, "--adapter", tools["adapter-mmg"], working_directory)
    if stage == "proper-rigid-publication":
        _require_path_option(command, "--ownership-runtime", tools["ownership-runtime"],
                             working_directory)
        _require_path_option(command, "--ownership-auditor", tools["ownership-auditor"],
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
            not isinstance(report.get("Environment"), dict)):
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
    validate_tool_invocation(stage, report["Command"],
                             {role: item["Path"] for role, item in tools.items()},
                             report.get("WorkingDirectory"))
    if stage == "canonical-source-validation":
        transform = json.loads(Path(report["Inputs"]["canonical-transform"]["Path"]).read_text())
        if isinstance(transform, dict): transform = transform.get("Transform")
        if transform != [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]:
            raise ValueError("Canonical source validation must use the identity transform")
    elif stage == "metric-preparation":
        _require_path_option(report["Command"], "--semantic-contract",
                             report["Inputs"]["canonical-semantic-contract"]["Path"],
                             report.get("WorkingDirectory"))
        _require_path_option(report["Command"], "--transformed-supports",
                             report["Inputs"]["canonical-supports"]["Path"],
                             report.get("WorkingDirectory"))
    elif stage == "native-adaptation-mmg":
        command = report["Command"]
        recipe = report["Inputs"]["restoration-recipe"]["Path"]
        if str(Path(recipe).resolve()) not in [str(_resolved_argument(value,
                Path(report.get("WorkingDirectory") or Path.cwd()))) for value in command]:
            raise ValueError("Native adaptation command did not consume the metric recipe")
        receipt = json.loads(Path(report["Artifacts"]["adaptation-receipt"]["Path"]).read_text())
        recipe_data = json.loads(Path(recipe).read_text())
        hmax = recipe_data.get("FarFieldBudgetPolicy", {}).get("EffectiveFarSize")
        if (receipt.get("RecipeSHA256") != sha256(recipe) or
                receipt.get("EffectiveFarSize") != hmax or
                float(receipt.get("HmaxArgument", "nan")) != hmax or
                receipt.get("OutputSHA256") != report["Artifacts"]["adapted-mesh"]["SHA256"]):
            raise ValueError("Native MMG hmax receipt differs from the bound metric policy")
    return report


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


def validate_canonical_dag(report_paths, canonical_mesh, launcher_name=None,
                           launcher_sha256=None, expected_tool_sha256=None):
    reports, digests = _validate_reports(report_paths, CANONICAL_STAGE_ORDER,
                                         launcher_name, launcher_sha256,
                                         expected_tool_sha256)
    source = reports["canonical-source-validation"]
    seed = reports["seed-generation"]["Artifacts"]["seed-mesh"]["SHA256"]
    metric = reports["metric-preparation"]
    adaptation = reports["native-adaptation-mmg"]
    restoration = reports["label-restoration"]
    publication = reports["canonical-gmsh-publication"]
    links = (
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
            record.get("CanonicalArtifacts", {}).get("candidate-mesh", {}).get("SHA256") or
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
    canonical_mesh = record["CanonicalArtifacts"]["candidate-mesh"]["Path"]
    canonical_expected = None if expected_tool_sha256 is None else {
        name: expected_tool_sha256[name] for name in CANONICAL_STAGE_ORDER}
    placement_expected = None if expected_tool_sha256 is None else {
        name: expected_tool_sha256[name] for name in PLACEMENT_STAGE_ORDER}
    canonical, canonical_digests = validate_canonical_dag(
        canonical_paths, canonical_mesh, launcher_name, launcher_sha256, canonical_expected)
    placement, placement_digests = validate_placement_dag(
        placement_paths, final_mesh, record_path, launcher_name, launcher_sha256,
        placement_expected)
    return {**canonical, **placement}, canonical_digests | placement_digests
