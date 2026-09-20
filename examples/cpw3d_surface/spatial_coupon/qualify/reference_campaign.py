# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The reference campaign of a coupon: the graded_v2 inputs (the producer's Palace
config, the trace files, the basis contract) and, when the reference ran to
completion, its reducer response matrices.

A reference campaign directory has the layout of the per-run `reference/` trees of
coupon-accuracy-assessment-20260913: `inputs-<key>/` (spatial_fabricated.json,
traces/basis-NNNN.csv and any other PrescribedPotential DataFile, basis-contract.json,
plan-view-boundary.csv, optionally retained-etch.csv) and `case-<key>-fabricated/`
(worker.json = the config the reference actually ran; reducer/domain- and
surface-response-matrix.csv).  The coupon is matched to its inputs BY CONTENT: the
SHA-256 of the case's frozen basis-contract.json (the manifest binding) equals the
inputs' basis-contract.json digest.  A case without matching inputs has no
configuration to run; a case with inputs but no reducer matrices is run and marked
PendingQualification.
"""
import hashlib
import json
from pathlib import Path

INPUTS_PREFIX = "inputs-"
CASE_SUFFIX = "-fabricated"
REFERENCE_CONFIG = "spatial_fabricated.json"
RAN_CONFIG = "worker.json"
MATRICES = ("domain-response-matrix.csv", "surface-response-matrix.csv")
BASIS_CONTRACT = "basis-contract.json"
TRACES_DIRECTORY = "traces"


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


class ReferenceError(ValueError):
    """A fail-closed reference resolution (recorded per case, never a crash of the run)."""


def find_inputs(campaign, basis_contract_sha256):
    """The `inputs-<key>` directory whose basis-contract.json has the given digest."""
    campaign = Path(campaign)
    matches = []
    for directory in sorted(campaign.glob(f"{INPUTS_PREFIX}*")):
        contract = directory / BASIS_CONTRACT
        if directory.is_dir() and contract.is_file() and sha256(contract) == basis_contract_sha256:
            matches.append(directory)
    if len(matches) > 1:
        raise ReferenceError(f"{len(matches)} inputs directories bind the same basis contract: {matches}")
    return matches[0] if matches else None


def resolve(campaign, basis_contract_sha256):
    """Inputs, reference config, sources and reducer matrices of the coupon under the
    campaign directory (None when the campaign is None or binds no such contract)."""
    if campaign is None:
        return None
    inputs = find_inputs(campaign, basis_contract_sha256)
    if inputs is None:
        return None
    key = inputs.name[len(INPUTS_PREFIX):]
    case_directory = Path(campaign) / f"case-{key}{CASE_SUFFIX}"
    ran = case_directory / RAN_CONFIG
    producer = inputs / REFERENCE_CONFIG
    if ran.is_file():
        config_path, config_role = ran, "the config the reference ran (case worker.json)"
    elif producer.is_file():
        config_path, config_role = producer, "the producer's config (inputs spatial_fabricated.json)"
    else:
        raise ReferenceError(f"{inputs} carries no {REFERENCE_CONFIG} and {case_directory} no {RAN_CONFIG}")
    config = json.loads(config_path.read_text())
    reducer = case_directory / "reducer"
    results = reducer if all((reducer / name).is_file() for name in MATRICES) else None
    contract = json.loads((inputs / BASIS_CONTRACT).read_text())
    sources = source_files(config, inputs)
    return {"Key": key, "Inputs": str(inputs), "Config": str(config_path), "ConfigRole": config_role,
            "ConfigSHA256": sha256(config_path), "ReferenceOrder": int(config["Solver"]["Order"]),
            "LinearTol": config["Solver"]["Linear"]["Tol"],
            "ZeroTraceIndices": [int(i) for i in contract.get("ZeroTraceIndices", [])],
            "Sources": sources, "Results": str(results) if results else None,
            "ResultsSHA256": ({name: sha256(results / name) for name in MATRICES} if results else None),
            "PlanViewBoundary": str(inputs / "plan-view-boundary.csv") if (inputs / "plan-view-boundary.csv").is_file() else None,
            "RetainedEtch": str(inputs / "retained-etch.csv") if (inputs / "retained-etch.csv").is_file() else None}


def source_files(config, inputs):
    """Index -> local DataFile path and digest of every PrescribedPotential source; the
    file is looked up by name under inputs/traces then inputs (the conductor terminal
    traces of multi-conductor cases live next to the basis traces)."""
    inputs = Path(inputs)
    sources = []
    for entry in config["Boundaries"]["PrescribedPotential"]:
        name = Path(entry["DataFile"]).name
        candidates = [inputs / TRACES_DIRECTORY / name, inputs / name]
        path = next((candidate for candidate in candidates if candidate.is_file()), None)
        if path is None:
            raise ReferenceError(f"source {entry['Index']} DataFile {name} is not under {inputs / TRACES_DIRECTORY} or {inputs}")
        sources.append({"Index": int(entry["Index"]), "Name": name, "Path": str(path), "SHA256": sha256(path),
                        "Terminal": bool(entry.get("TerminalAttributes"))})
    indices = [source["Index"] for source in sources]
    if len(set(indices)) != len(indices):
        raise ReferenceError(f"repeated PrescribedPotential indices {indices}")
    return sources
