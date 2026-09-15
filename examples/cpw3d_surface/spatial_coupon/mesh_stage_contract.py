#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Validation for the bounded mesh-production stage DAG."""
import hashlib
import json
from pathlib import Path


STAGE_ORDER = ("seed-generation", "metric-preparation", "native-adaptation-mmg",
               "label-restoration", "final-gmsh-publication")
STAGE_TOOLS = {
    "seed-generation": {"runtime", "mesher"},
    "metric-preparation": {"runtime", "metric-preparer"},
    "native-adaptation-mmg": {"runtime", "adapter-mmg"},
    "label-restoration": {"runtime", "label-restorer"},
    "final-gmsh-publication": {"runtime", "publisher"},
}
STAGE_INPUTS = {
    "seed-generation": set(),
    "metric-preparation": {"seed-mesh"},
    "native-adaptation-mmg": {"mmg-seed", "metric", "pins", "fixed-triangles"},
    "label-restoration": {"adapted-mesh", "restoration-recipe"},
    "final-gmsh-publication": {"restored-mesh"},
}
STAGE_OUTPUTS = {
    "seed-generation": {"seed-mesh"},
    "metric-preparation": {"mmg-seed", "metric", "pins", "fixed-triangles",
                           "restoration-recipe"},
    "native-adaptation-mmg": {"adapted-mesh"},
    "label-restoration": {"restored-mesh"},
    "final-gmsh-publication": {"candidate-mesh", "ownership-partition"},
}


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


def validate_stage_report(report, stage, launcher_name=None, launcher_sha256=None):
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
    command_paths = {str(Path(value).resolve()) for value in report["Command"]
                     if isinstance(value, str) and Path(value).is_file()}
    for role, item in tools.items():
        if (not isinstance(item, dict) or not item.get("Path") or
                not item.get("SHA256") or not Path(item["Path"]).is_file() or
                sha256(item["Path"]) != item["SHA256"] or
                str(Path(item["Path"]).resolve()) not in command_paths):
            raise ValueError(f"Stage tool was not invoked: {stage}/{role}")
    return report


def validate_stage_dag(report_paths, final_mesh, launcher_name=None, launcher_sha256=None):
    if set(report_paths) != set(STAGE_ORDER):
        raise ValueError("Exactly one report for every bounded mesh stage is required")
    reports = {}
    report_digests = set()
    for stage in STAGE_ORDER:
        path = Path(report_paths[stage])
        report = validate_stage_report(json.loads(path.read_text()), stage,
                                       launcher_name, launcher_sha256)
        digest = sha256(path)
        if digest in report_digests:
            raise ValueError("Bounded stage reports must be content-distinct")
        report_digests.add(digest)
        reports[stage] = report

    seed = reports["seed-generation"]["Artifacts"]["seed-mesh"]["SHA256"]
    metric = reports["metric-preparation"]
    adaptation = reports["native-adaptation-mmg"]
    restoration = reports["label-restoration"]
    publication = reports["final-gmsh-publication"]
    links = (
        (metric["Inputs"]["seed-mesh"]["SHA256"], seed),
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
