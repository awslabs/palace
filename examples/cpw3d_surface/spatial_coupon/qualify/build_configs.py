#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Derive the worker / reducer / local-edge Palace configs of a coupon from its
reference config (four-edge-physics-11 / gallery-physics-06b build_configs.py,
parametrized).

Only Model.Mesh, Problem.Output, the PrescribedPotential DataFile directory, the
source subset (controls) and Solver.Order (control orders) change; the local-edge
stage additionally sets SaveLocalEdgeEnergy true on every Dielectric interface
entry.  Every other entry - materials, interfaces, Linear.Tol, MaxIts, the
electrostatic response-matrix options - is byte-for-byte the reference value.

usage: build_configs.py REFERENCE_CONFIG --mesh REMOTE_MESH --traces REMOTE_TRACES_DIR
       --output-root REMOTE_STAGE_DIR --out DIR [--source I ...] [--order P]
       [--save-local-edge-energy]
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def derive(reference, mesh, output_root, sources, traces_directory, *, order=None, save_local_edge_energy=None):
    """(worker, reducer) configs on `mesh` for the given source indices (in reference
    order) whose DataFiles live under `traces_directory` (remote path)."""
    config = copy.deepcopy(reference)
    config["Model"]["Mesh"] = mesh
    if order is not None:
        config["Solver"]["Order"] = int(order)
    if save_local_edge_energy is not None:
        for entry in config["Boundaries"]["Postprocessing"]["Dielectric"]:
            if not (entry.get("LocalizeEdgeEnergy") is True and entry.get("SaveLocalEdgeEnergy") is False):
                raise ValueError("the reference interfaces must localize edge energy without saving it "
                                 f"(entry {entry.get('Index')}): the local-edge stage flips SaveLocalEdgeEnergy only")
            entry["SaveLocalEdgeEnergy"] = save_local_edge_energy
    wanted = set(sources)
    keep = []
    for entry in config["Boundaries"]["PrescribedPotential"]:
        if entry["Index"] in wanted:
            entry = dict(entry)
            entry["DataFile"] = f"{traces_directory}/{Path(entry['DataFile']).name}"
            keep.append(entry)
    if [entry["Index"] for entry in keep] != list(sources):
        raise ValueError(f"sources {list(sources)} are not a subsequence of the reference PrescribedPotential indices")
    config["Boundaries"]["PrescribedPotential"] = keep
    worker = copy.deepcopy(config)
    worker["Problem"]["Output"] = f"{output_root}/worker"
    reducer = copy.deepcopy(config)
    reducer["Problem"]["Output"] = f"{output_root}/reducer"
    return worker, reducer


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")
    return sha256(path)


def write_stage(directory, worker, reducer):
    """worker.json and reducer.json under `directory`; returns their digests."""
    return {"worker.json": write_json(Path(directory) / "worker.json", worker),
            "reducer.json": write_json(Path(directory) / "reducer.json", reducer)}


def write_local_edge_stage(directory, config, output):
    """The ordinary-path config.json (its own output directory, no worker / reducer)."""
    config = copy.deepcopy(config)
    config["Problem"]["Output"] = output
    return {"config.json": write_json(Path(directory) / "config.json", config)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reference_config", type=Path)
    parser.add_argument("--mesh", required=True, help="remote mesh path written into Model.Mesh")
    parser.add_argument("--traces", required=True, help="remote directory of the trace files")
    parser.add_argument("--output-root", required=True, help="remote stage directory (Problem.Output = ROOT/worker, ROOT/reducer)")
    parser.add_argument("--out", type=Path, required=True, help="local directory receiving worker.json / reducer.json")
    parser.add_argument("--source", type=int, action="append", help="source index (repeatable; default: every reference source)")
    parser.add_argument("--order", type=int)
    parser.add_argument("--save-local-edge-energy", action="store_true")
    args = parser.parse_args(argv)
    reference = json.loads(args.reference_config.read_text())
    sources = args.source or [entry["Index"] for entry in reference["Boundaries"]["PrescribedPotential"]]
    worker, reducer = derive(reference, args.mesh, args.output_root, sources, args.traces, order=args.order,
                             save_local_edge_energy=True if args.save_local_edge_energy else None)
    print(json.dumps(write_stage(args.out, worker, reducer), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
