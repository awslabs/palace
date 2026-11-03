#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""device_gate_configs.py -- the (F) chain's device gate solve configs (the "c0" solves) from a device run config.

The (F) qualification of a spatial coupon (coupon_library.py spatial-qualify evaluate --gate palace.json ...) reads the
resolution-aware conductor-consistency gate (spatial_qualification.resolution_aware_gate: Count 0, MaxRatio <= 1e-4 on every
supplied device solve of the registration mesh, decisions 277 / 299) from the device's CHEAP thin solves at Orders 3 / 4 / 5 on
the registration mesh with the interim library -- the c0 solves -- preceded by the geometry-only PREFLIGHT. This tool writes
those configs from the window's device run config in one of two FORMS (fix round R-4, OI-QUALIFY-DEVICE-GATE-ES-ONLY,
decisions 1174 / 1238 / 1239):

  * `--form es` (the DEFAULT; the form of record): the ELECTROSTATIC device config (Problem.Type Electrostatic; the terminal
    excitations) -> `preflight.json` (Problem.Output <out>/preflight/postpro; Solver.Electrostatic.ResponseCorrection.Library =
    the interim library; GeometryFailurePolicy "RecordRaw" on the preflight only, decisions 1029 (2) / 1032) and `c0-p<N>.json`
    per order (Solver.Order N; Electrostatic.Save 0; ResponseCorrection.TraceCoupling "SurfaceMortar"; Output <out>/c0-p<N>/
    postpro); Model.Refinement.MaxIts must already be 0 (the c0 is the registration mesh, no AMR). This is the recipe of the
    (F) chains of record (write_c0_configs v2: rerun-3 R3 -> ... -> spatial-rebuild) reproduced BYTE-IDENTICALLY
    (test_device_gate_configs.py): the chain scripts themselves are untouched.
  * `--form driven`: the DRIVEN device config (Problem.Type Driven; the F-39 LumpedPort excitation; decision 872 / 874) at ONE
    frequency `--frequency F` (GHz; required: no silent default) -> the same two files per the Maxwell path's leaves:
    Solver.SurfaceResponseCorrection.Library = the interim library; TraceCoupling "Collocated" (the Maxwell path refuses
    SurfaceMortar: surfaceresponseoperator.cpp "Maxwell SurfaceMortar response trace coupling is not implemented yet"; the
    Maxwell gate probes the explicit TraceMesh on the collocated path, idfix R0-a); Solver.Driven.Samples = one Point sample
    [F] with SaveStep 0; Model.Refinement.MaxIts set to 0 (a device Driven config of record carries the AMR ladder; the c0
    is the registration mesh). The gate record is the SAME form on both paths (FinishConductorConsistencyGate writes
    SurfaceResponse.Diagnostics.ConductorConsistency with Maxwell / Frequency per record), so spatial_qualification's gate
    reader reads it unchanged against the same bars. UNVERIFIED: no Maxwell-form c0 has been solved yet; the measured
    acceptance test is P3's c0 on the FE6 9f8d interim library (~0.3 node-h, before any F use). A Driven window without an
    electrostatic terminal (FE6: electrostaticsolver.cpp "No terminal or prescribed potential boundaries", 1174) is the
    form's reason.

A base config of the wrong Problem.Type for the form, a Driven base without an excited LumpedPort, or a non-zero
Refinement.MaxIts under `--form es` is REFUSED (ValueError, the reason named). Existing output files are never overwritten
(475). Beside the configs `device-gate-configs.record.json` names the form, the base (sha256), the library, the orders, the
frequency and the leaves set. Pure Python.

usage:
  device_gate_configs.py --base CONFIG --library LIB --out OUTDIR [--form es|driven] [--orders 3 4 5] [--frequency F]
"""

import argparse
import copy
import hashlib
import json
import sys
from pathlib import Path

FORMS = ("es", "driven")
DEFAULT_ORDERS = (3, 4, 5)
PREFLIGHT_NAME = "preflight.json"
RECORD_NAME = "device-gate-configs.record.json"
PROBLEM_TYPE = {"es": "Electrostatic", "driven": "Driven"}
# The device-side leaf of the response correction per form (configfile.cpp: Solver.Electrostatic.ResponseCorrection on the
# electrostatic path; Solver.SurfaceResponseCorrection, parsed by the same reader, on the driven / eigenmode paths).
CORRECTION_PATH = {"es": ("Solver", "Electrostatic", "ResponseCorrection"), "driven": ("Solver", "SurfaceResponseCorrection")}
TRACE_COUPLING = {"es": "SurfaceMortar", "driven": "Collocated"}
PREFLIGHT_GEOMETRY_POLICY = "RecordRaw"


class DeviceGateConfigError(ValueError):
    """The base config does not fit the requested form (the reason is named)."""


def c0_name(order):
    return f"c0-p{int(order)}.json"


def dump_config(config):
    """The chain recipe's serialisation: json.dumps(indent=1) + newline."""
    return json.dumps(config, indent=1) + "\n"


def _leaf(config, path):
    node = config
    for key in path:
        if not isinstance(node, dict) or key not in node:
            raise DeviceGateConfigError(f"the base config carries no {'.'.join(path)} block")
        node = node[key]
    return node


def check_base(config, form):
    """The base config fits the form: Problem.Type, the correction block, (es) Refinement.MaxIts 0, (driven) an excited port."""
    problem_type = config.get("Problem", {}).get("Type")
    if problem_type != PROBLEM_TYPE[form]:
        raise DeviceGateConfigError(f"--form {form} needs a Problem.Type {PROBLEM_TYPE[form]} base config, not {problem_type!r}")
    _leaf(config, CORRECTION_PATH[form])
    if form == "es":
        max_its = config.get("Model", {}).get("Refinement", {}).get("MaxIts")
        if max_its != 0:
            raise DeviceGateConfigError(f"the electrostatic c0 base must carry Model.Refinement.MaxIts 0 (the registration mesh, "
                                        f"no AMR), not {max_its!r}")
    else:
        ports = config.get("Boundaries", {}).get("LumpedPort") or []
        if not any(port.get("Excitation") for port in ports):
            raise DeviceGateConfigError("--form driven needs a Boundaries.LumpedPort entry with Excitation true (the F-39 port form)")
        _leaf(config, ("Solver", "Driven"))


def common_config(base, form, library, frequency):
    """The leaves shared by the preflight and every c0 solve: the library (and, driven, the one-frequency sample + MaxIts 0)."""
    config = copy.deepcopy(base)
    _leaf(config, CORRECTION_PATH[form])["Library"] = library
    if form == "driven":
        config["Solver"]["Driven"]["Samples"] = [{"Type": "Point", "Freq": [float(frequency)], "SaveStep": 0}]
        for legacy in ("MinFreq", "MaxFreq", "FreqStep", "SaveStep"):
            config["Solver"]["Driven"].pop(legacy, None)
        config.setdefault("Model", {}).setdefault("Refinement", {})["MaxIts"] = 0
    return config


def preflight_config(common, form, out):
    config = copy.deepcopy(common)
    config["Problem"]["Output"] = f"{out}/preflight/postpro"
    _leaf(config, CORRECTION_PATH[form])["GeometryFailurePolicy"] = PREFLIGHT_GEOMETRY_POLICY
    return config


def c0_config(common, form, out, order):
    config = copy.deepcopy(common)
    config["Solver"]["Order"] = int(order)
    if form == "es":
        config["Solver"]["Electrostatic"]["Save"] = 0
    _leaf(config, CORRECTION_PATH[form])["TraceCoupling"] = TRACE_COUPLING[form]
    config["Problem"]["Output"] = f"{out}/c0-p{int(order)}/postpro"
    return config


def device_gate_configs(base, *, form, library, out, orders=DEFAULT_ORDERS, frequency=None):
    """{file name -> config dict} for the form; `out` is the output root AS WRITTEN into Problem.Output (a string)."""
    if form not in FORMS:
        raise DeviceGateConfigError(f"unknown form {form!r}; one of {FORMS}")
    if form == "driven" and frequency is None:
        raise DeviceGateConfigError("--form driven needs --frequency F (GHz): the one frequency of the c0 solves, no default")
    check_base(base, form)
    common = common_config(base, form, library, frequency)
    configs = {PREFLIGHT_NAME: preflight_config(common, form, out)}
    for order in orders:
        configs[c0_name(order)] = c0_config(common, form, out, order)
    return configs


def write_configs(base_path, *, form, library, out, orders=DEFAULT_ORDERS, frequency=None):
    """Writes the configs and the record under `out`; refuses to overwrite (475). Returns the record."""
    base_path = Path(base_path)
    base = json.loads(base_path.read_text())
    out_dir = Path(out)
    # The chain recipe wrote f"{Path(outdir)}/...": the Path-normalised string names the Output root.
    configs = device_gate_configs(base, form=form, library=library, out=str(out_dir), orders=orders, frequency=frequency)
    existing = [name for name in list(configs) + [RECORD_NAME] if (out_dir / name).exists()]
    if existing:
        raise DeviceGateConfigError(f"refusing to overwrite {existing} in {out_dir} (decision 475): use a fresh --out")
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, config in configs.items():
        (out_dir / name).write_text(dump_config(config))
    record = {
        "Tool": "device_gate_configs.py (fix round R-4; decisions 1174 / 1238 / 1239)",
        "Form": form,
        "ProblemType": PROBLEM_TYPE[form],
        "Base": str(base_path),
        "BaseSha256": hashlib.sha256(base_path.read_bytes()).hexdigest(),
        "Library": library,
        "Orders": [int(o) for o in orders],
        "FrequencyGHz": None if frequency is None else float(frequency),
        "Files": sorted(configs),
        "Leaves": {"Correction": ".".join(CORRECTION_PATH[form]), "TraceCoupling": TRACE_COUPLING[form],
                   "PreflightGeometryFailurePolicy": PREFLIGHT_GEOMETRY_POLICY,
                   "RefinementMaxIts": 0, "Save": 0 if form == "es" else "Driven.Samples[].SaveStep 0"},
        "Verification": ("the form of record (the (F) chains' write_c0_configs v2 recipe, byte-identical)" if form == "es" else
                         "UNVERIFIED: no Maxwell-form c0 solved yet; the acceptance test = P3's c0 on the FE6 9f8d interim "
                         "library (~0.3 node-h) before any F use (decision 1239)"),
    }
    (out_dir / RECORD_NAME).write_text(json.dumps(record, indent=1) + "\n")
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--base", required=True, help="the window's device run config (Electrostatic for --form es, Driven for --form driven)")
    parser.add_argument("--library", required=True, help="the interim process library (process-library-abs.json) the c0 solves apply")
    parser.add_argument("--out", required=True, help="the output root: the configs' directory AND the Problem.Output root written into them")
    parser.add_argument("--form", choices=FORMS, default="es", help="es (default; the form of record) | driven (UNVERIFIED until P3's 9f8d c0)")
    parser.add_argument("--orders", type=int, nargs="+", default=list(DEFAULT_ORDERS), help="the c0 orders (default 3 4 5)")
    parser.add_argument("--frequency", type=float, help="--form driven: the one frequency (GHz) of the c0 solves (required, no default)")
    args = parser.parse_args(argv)
    try:
        record = write_configs(args.base, form=args.form, library=args.library, out=args.out, orders=args.orders,
                               frequency=args.frequency)
    except DeviceGateConfigError as error:
        print(f"device_gate_configs.py: {error}", file=sys.stderr)
        return 2
    print(f"configs written: {record['Files']} form {record['Form']} base {record['Base']} out {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
