#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Validation for the bounded mesh-production stage DAG."""
import hashlib
import json
from pathlib import Path


STAGE_ORDER = ("source-transformation", "seed-generation", "metric-preparation",
               "native-adaptation-mmg", "label-restoration", "final-gmsh-publication")
STAGE_TOOLS = {
    "source-transformation": {"runtime", "source-transformer"},
    "seed-generation": {"runtime", "mesher"},
    "metric-preparation": {"runtime", "metric-preparer"},
    "native-adaptation-mmg": {"runtime", "adapter-mmg"},
    "label-restoration": {"runtime", "label-restorer"},
    "final-gmsh-publication": {"runtime", "publisher"},
}
STAGE_INPUTS = {
    "source-transformation": {"source-semantic-contract", "source-signature",
                              "source-boundary", "source-mask", "canonical-transform"},
    "seed-generation": set(),
    "metric-preparation": {"seed-mesh", "transformed-semantic-contract",
                           "transformed-supports"},
    "native-adaptation-mmg": {"mmg-seed", "metric", "pins", "fixed-triangles"},
    "label-restoration": {"adapted-mesh", "restoration-recipe"},
    "final-gmsh-publication": {"restored-mesh"},
}
STAGE_OUTPUTS = {
    "source-transformation": {"transformed-semantic-contract", "transformed-supports"},
    "seed-generation": {"seed-mesh"},
    "metric-preparation": {"mmg-seed", "metric", "pins", "fixed-triangles",
                           "restoration-recipe"},
    "native-adaptation-mmg": {"adapted-mesh"},
    "label-restoration": {"source-local-restored-mesh", "restored-mesh"},
    "final-gmsh-publication": {"candidate-mesh", "ownership-partition"},
}
STAGE_PRIMARY_TOOL = {
    "source-transformation": "source-transformer",
    "seed-generation": "mesher",
    "metric-preparation": "metric-preparer",
    "native-adaptation-mmg": "adapter-mmg",
    "label-restoration": "label-restorer",
    "final-gmsh-publication": "publisher",
}


# Interpreter options whose following argument is option data, not a script.  Unknown
# options fail closed by making their following argument the apparent script.
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
_INTERPRETER_EXECUTION_OPTIONS = {"-c", "-m", "-e", "--eval", "--print", "-L",
                                  "--load"}


def _first_interpreter_script(command, primary):
    """Return the first script position after recognized interpreter options."""
    index = 1
    while index < len(command):
        argument = command[index]
        if argument == "--":
            return index + 1 if index + 1 < len(command) else None
        if (argument in _INTERPRETER_EXECUTION_OPTIONS or
                (argument == "-E" and Path(primary).suffix == ".jl")):
            # These execute code, a module, or a loaded file before any positional script.
            return None
        if argument in _INTERPRETER_OPTIONS_WITH_VALUE:
            index += 2
            continue
        if (argument in _INTERPRETER_OPTIONS_WITHOUT_VALUE or
                (argument == "-E" and Path(primary).suffix == ".py") or
                argument.startswith(_INTERPRETER_OPTIONS_WITH_ATTACHED_VALUE)):
            index += 1
            continue
        return index
    return None


def validate_tool_invocation(stage, command, tools, working_directory=None):
    """Require the declared runtime and stage tool in executable/script positions."""
    base = Path(working_directory or Path.cwd())
    runtime = str(Path(tools["runtime"]).resolve())
    primary = str(Path(tools[STAGE_PRIMARY_TOOL[stage]]).resolve())

    def matches(argument, expected):
        path = Path(argument)
        resolved = path if path.is_absolute() else base / path
        return (str(resolved.resolve()) == expected or
                (not path.is_absolute() and expected.endswith("/" + argument)))

    if not matches(command[0], runtime):
        raise ValueError(f"Stage runtime is not argv[0]: {stage}")
    if stage == "native-adaptation-mmg":
        if primary != runtime:
            raise ValueError("Native adapter must be both runtime and adapter-mmg")
        return
    script_position = _first_interpreter_script(command, primary)
    if script_position is None or not matches(command[script_position], primary):
        raise ValueError(f"Stage tool is not in the executed script position: {stage}")


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
        if (not isinstance(item, dict) or not item.get("Path") or
                not item.get("SHA256") or not Path(item["Path"]).is_file() or
                sha256(item["Path"]) != item["SHA256"]):
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
        if (not isinstance(item, dict) or not item.get("Path") or
                not item.get("SHA256") or not Path(item["Path"]).is_file() or
                sha256(item["Path"]) != item["SHA256"] or
                (expected_tool_sha256 is not None and
                 item["SHA256"] != expected_tool_sha256.get(role))):
            raise ValueError(f"Stage tool binding changed: {stage}/{role}")
    validate_tool_invocation(stage, report["Command"],
                             {role: item["Path"] for role, item in tools.items()},
                             report.get("WorkingDirectory"))
    return report


def validate_stage_dag(report_paths, final_mesh, launcher_name=None, launcher_sha256=None,
                       expected_tool_sha256=None):
    if set(report_paths) != set(STAGE_ORDER):
        raise ValueError("Exactly one report for every bounded mesh stage is required")
    reports = {}
    report_digests = set()
    for stage in STAGE_ORDER:
        path = Path(report_paths[stage])
        report = validate_stage_report(
            json.loads(path.read_text()), stage, launcher_name, launcher_sha256,
            None if expected_tool_sha256 is None else expected_tool_sha256[stage])
        digest = sha256(path)
        if digest in report_digests:
            raise ValueError("Bounded stage reports must be content-distinct")
        report_digests.add(digest)
        reports[stage] = report

    source = reports["source-transformation"]
    seed = reports["seed-generation"]["Artifacts"]["seed-mesh"]["SHA256"]
    metric = reports["metric-preparation"]
    adaptation = reports["native-adaptation-mmg"]
    restoration = reports["label-restoration"]
    publication = reports["final-gmsh-publication"]
    links = (
        (metric["Inputs"]["seed-mesh"]["SHA256"], seed),
        (metric["Inputs"]["transformed-semantic-contract"]["SHA256"],
         source["Artifacts"]["transformed-semantic-contract"]["SHA256"]),
        (metric["Inputs"]["transformed-supports"]["SHA256"],
         source["Artifacts"]["transformed-supports"]["SHA256"]),
        (adaptation["Inputs"]["mmg-seed"]["SHA256"],
         metric["Artifacts"]["mmg-seed"]["SHA256"]),
        (adaptation["Inputs"]["metric"]["SHA256"],
         metric["Artifacts"]["metric"]["SHA256"]),
        (adaptation["Inputs"]["pins"]["SHA256"],
         metric["Artifacts"]["pins"]["SHA256"]),
        (adaptation["Inputs"]["fixed-triangles"]["SHA256"],
         metric["Artifacts"]["fixed-triangles"]["SHA256"]),
        (restoration["Inputs"]["adapted-mesh"]["SHA256"],
         adaptation["Artifacts"]["adapted-mesh"]["SHA256"]),
        (restoration["Inputs"]["restoration-recipe"]["SHA256"],
         metric["Artifacts"]["restoration-recipe"]["SHA256"]),
        (publication["Inputs"]["restored-mesh"]["SHA256"],
         restoration["Artifacts"]["restored-mesh"]["SHA256"]),
        (publication["Artifacts"]["candidate-mesh"]["SHA256"], sha256(final_mesh)),
    )
    if any(actual != expected for actual, expected in links):
        raise ValueError("Bounded mesh stage input/output digest chain is broken")
    return reports, report_digests
