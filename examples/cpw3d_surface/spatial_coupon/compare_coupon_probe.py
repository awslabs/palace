#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Compare coupon probes without gating on singular raw thin-metal SPR.

Thin cases: domain energy only. Fabricated cases: domain and interface energies.
A passing probe is not a qualification of the entire response matrix.
"""
import argparse
import csv
import json
import math
from pathlib import Path


def rows(path):
    with path.open(newline="") as stream:
        return [{key.strip(): float(value) for key, value in row.items()}
                for row in csv.DictReader(stream, skipinitialspace=True)]


def probe_energies(directory, kind):
    """Read ordinary probe outputs or diagonals of an archived response matrix.

    Off-diagonal matrix entries are deliberately outside this probe-only gate.
    Never infer localized interface energy from a compact/ambiguous matrix schema.
    """
    domain, interfaces = {}, {}
    ordinary = (directory / "domain-E.csv").exists()
    seen = set()
    data = rows(directory / ("domain-E.csv" if ordinary else "domain-response-matrix.csv"))
    for row in data:
        source = row["i"] if ordinary else row["basis_i"]
        key = source if ordinary else (source, row["basis_j"])
        if key in seen:
            raise ValueError("duplicate source or matrix entry")
        seen.add(key)
        if ordinary or source == row["basis_j"]:
            domain[source] = row["E_elec (J)" if ordinary else "Q_ij (J)"]
    if kind == "fabricated":
        if ordinary:
            for row in rows(directory / "surface-Q.csv"):
                source = row["i"]
                if source not in domain or source in interfaces:
                    raise ValueError("surface/domain source indices differ or repeat")
                interfaces[source] = {key: value * domain[source] for key, value in row.items()
                                      if key.startswith("p_surf[")}
        else:
            seen = set()
            for row in rows(directory / "surface-response-matrix.csv"):
                if "Q_ij (J)" not in row or "Q_total_ij (J)" not in row:
                    raise ValueError("ambiguous compact surface matrix schema")
                key = (row["basis_i"], row["basis_j"], row["interface"])
                if key in seen:
                    raise ValueError("duplicate surface matrix entry")
                seen.add(key)
                if row["basis_i"] == row["basis_j"]:
                    source = row["basis_i"]
                    if source not in domain:
                        raise ValueError("surface/domain source indices differ")
                    index = row["interface"]
                    if index != int(index):
                        raise ValueError("noninteger interface index")
                    interfaces.setdefault(source, {})[f"p_surf[{int(index)}]"] = row["Q_total_ij (J)"]
        if domain.keys() != interfaces.keys():
            raise ValueError("surface/domain source indices differ")
    return domain, interfaces


# Provisional engineering targets, not automatic library qualification. The strict
# profile is retained for diagnostics when reference convergence warrants it.
PROFILES = {"engineering": (5e-4, 5e-3), "strict": (1e-4, 1e-3)}


def compare(reference, candidate, kind, domain_tol=5e-4, surface_tol=5e-3):
    if kind not in ("thin", "fabricated"):
        raise ValueError("kind must be thin or fabricated")
    if not all(math.isfinite(x) and x > 0 for x in (domain_tol, surface_tol)):
        raise ValueError("tolerances must be finite and positive")
    a, sa = probe_energies(reference, kind)
    b, sb = probe_energies(candidate, kind)
    if not a or not b or a.keys() != b.keys():
        raise ValueError("source indices differ or are empty")
    checks = []

    def check(name, source, x, y, tolerance):
        if not math.isfinite(x) or not math.isfinite(y):
            raise ValueError(f"non-finite {name}")
        relative = abs(y-x)/abs(x) if x else (0.0 if y == 0 else None)
        checks.append({"Quantity": name, "Source": source, "Reference": x,
                       "Candidate": y, "RelativeError": relative,
                       "Tolerance": tolerance,
                       "Passed": relative is not None and relative <= tolerance})

    for source in a:
        check("DomainEnergy", source, a[source], b[source], domain_tol)
    if kind == "fabricated":
        for source in a:
            columns = {key for key in sa[source] if key.startswith("p_surf[")}
            if not columns or columns != {key for key in sb[source] if key.startswith("p_surf[")}:
                raise ValueError("surface interface indices differ or are empty")
            for key in sorted(columns):
                check("InterfaceEnergy:"+key, source,
                      sa[source][key], sb[source][key], surface_tol)
    return {"Version": 1, "Kind": kind,
            "Passed": all(check["Passed"] for check in checks), "Checks": checks,
            "ThinRawSPRGated": False,
            "Tolerances": {"DomainRelative": domain_tol, "InterfaceRelative": surface_tol},
            "LibraryQualified": False,
            "Scope": "Probe energies only; full operator, boundary-trace and geometry qualification remain separate."}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--kind", required=True, choices=("thin", "fabricated"))
    parser.add_argument("--profile", choices=tuple(PROFILES), default="engineering",
                        help="engineering: 0.05%% domain / 0.5%% interface; strict: 0.01%% / 0.1%%")
    parser.add_argument("--domain-tol", type=float, default=None)
    parser.add_argument("--surface-tol", type=float, default=None)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    default_domain, default_surface = PROFILES[args.profile]
    args.domain_tol = default_domain if args.domain_tol is None else args.domain_tol
    args.surface_tol = default_surface if args.surface_tol is None else args.surface_tol
    if not all(math.isfinite(x) and x > 0 for x in (args.domain_tol, args.surface_tol)):
        parser.error("tolerances must be finite and positive")
    report = compare(args.reference, args.candidate, args.kind,
                     args.domain_tol, args.surface_tol)
    report["Profile"] = args.profile
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if report["Passed"] else 1)


if __name__ == "__main__":
    main()
