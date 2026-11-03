#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Attach the magnetic (inductive) leaves of the vacuum twin to a Version-3 fabrication-process library.

mx-sc-magnetic DESIGN-v2 section 3 (USER 1221 D3 / D7; SUPERVISOR AMENDMENT 1202; ERRATUM 1227; decision 1229): in a
translational cross-section the 2D magnetostatic problem is the EXACT dual of the electrostatic one at eps = 1, so a
model's magnetic energy within R as a quadratic form of the flux trace is S_m = Q_dom,vac: the SAME electrostatic coupon
pair (thin / fabricated, the same meshes, hats and piece states) re-run with --substrate-permittivity 1.0. This tool
assembles those vacuum domain-response matrices onto a COPY of a library as the OPTIONAL per-model leaves
`MagneticFabricatedMatrix` / `MagneticThinMatrix` / `MagneticSource` (a library without them loads and behaves bitwise as
today) and writes the top-level `MagneticRules` (the r3 same-mesh tolerance and the two-level r4 PSD ratios). It
REFUSES a twin whose basis file is not byte-identical to the model's (r1), records the regenerated electric matrices'
maximum relative deviation from the model of record (the r3 statistic the loader gates) and the PSD statistic of
dS_m = S_fab - S_thin (the r4 statistic; information here, the loader recomputes it). Nothing of the input library is
rewritten (475 / 501): the output is a NEW library file whose electric entries are the input's.

Attachment manifest (--attach, JSON): a list of objects
  {"Model": <library model Name>,
   "VacuumThin": <domain-response-matrix.csv of the eps = 1 thin coupon>,
   "VacuumFabricated": <... of the eps = 1 fabricated coupon>,
   "VacuumBasisPoints": <the twin generator's basis_points.csv>,
   "RegeneratedElectricThin": <domain-response-matrix.csv of the regenerated eps = eps_r thin coupon> (optional),
   "RegeneratedElectricFabricated": <...> (optional; both or neither; absent = the deviation is NOT recorded: the
                                        loader refuses the leaf by r3 unless the library's SameMeshTolerance admits -1),
   "VacuumPermittivity": 1.0 (optional, default 1.0; the twin's Domains.Materials permittivity),
   "MeshSha256": {"Thin": ..., "Fabricated": ...} (optional provenance),
   "MeshOfRecordSha256Match": {"Thin": true/false, "Fabricated": ...} (optional provenance),
   "Toolchain": "...", "Command": "..." (optional provenance)}
"""

import argparse
import hashlib
import json
import os
import shutil
import sys

import numpy as np

MAGNETIC_SOURCE_VERSION = 1
FABRICATION_KEYS = ("MetalThickness", "OveretchDepth", "SidewallAngle", "TopRoundingRadius", "BottomRoundingRadius",
                    "SubstratePermittivity")
DEFAULT_RULES = {"SameMeshTolerance": 1.0e-12, "PsdFindingRatio": 1.0e-2, "PsdRefusalRatio": 0.1}


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_domain_response_matrix(path):
    """A Palace domain-response-matrix.csv (basis_i, basis_j, Q_ij (J); the upper triangle or the full matrix) as a
    dense symmetric numpy array. The CSV stores the upper triangle (Phase 0 step 9): both triangles are filled."""
    entries = []
    size = 0
    with open(path) as handle:
        header = handle.readline()
        if "basis_i" not in header:
            raise ValueError("%s is not a domain-response-matrix.csv (header %r)" % (path, header.strip()))
        for line in handle:
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 3 or not parts[0]:
                continue
            i, j, value = int(float(parts[0])), int(float(parts[1])), float(parts[2])
            entries.append((i, j, value))
            size = max(size, i, j)
    matrix = np.zeros((size, size))
    seen = np.zeros((size, size), dtype=bool)
    for i, j, value in entries:
        matrix[i - 1, j - 1] = value
        seen[i - 1, j - 1] = True
        if not seen[j - 1, i - 1]:
            matrix[j - 1, i - 1] = value
    return matrix


def relative_deviation(a, b):
    """max_ij |a - b| / max_ij |b| (the Phase-0 D_same / D_diff statistic; 0 for identical files)."""
    scale = np.max(np.abs(b))
    if scale == 0.0:
        return 0.0 if np.max(np.abs(a - b)) == 0.0 else float("inf")
    return float(np.max(np.abs(a - b)) / scale)


def psd_statistic(defect):
    """(-lambda_min / ||dS||_2, the negative eigenvalue count) of the symmetrised defect."""
    symmetric = 0.5 * (defect + defect.T)
    eigenvalues = np.linalg.eigvalsh(symmetric)
    norm = max(float(np.max(np.abs(eigenvalues))), 1.0e-300)
    return max(0.0, -float(eigenvalues.min())) / norm, int(np.count_nonzero(eigenvalues < 0.0))


def normalised_fabrication(fabrication):
    """The library's Fabrication block reduced to the six process keys the loader's r2 rule compares (the coupon
    generators write SidewallAngleDegrees; the library of record carries both spellings)."""
    out = {}
    for key in FABRICATION_KEYS:
        value = fabrication.get(key)
        if value is None and key == "SidewallAngle":
            value = fabrication.get("SidewallAngleDegrees")
        if value is None:
            raise ValueError("the library's Fabrication block has no %s: the twin's process cannot be recorded" % key)
        out[key] = float(value)
    return out


def resolve(path, base_dir):
    return path if os.path.isabs(path) else os.path.normpath(os.path.join(base_dir, path))


def attach(library, library_dir, attachments, output_dir, rules, now=None):
    """Return (the new library dict, the per-model records). Files are copied under output_dir/models/<name>/."""
    models_by_name = {model["Name"]: model for model in library.get("Models", [])}
    if library.get("Version") != 3:
        raise ValueError("magnetic leaves require a Version-3 library (got %r)" % library.get("Version"))
    fabrication = normalised_fabrication(library.get("Fabrication") or {})
    out = json.loads(json.dumps(library))
    out_models = {model["Name"]: model for model in out["Models"]}
    records = []
    for entry in attachments:
        name = entry["Model"]
        if name not in models_by_name:
            raise ValueError("model %r is not in the library" % name)
        model = models_by_name[name]
        if "Blend" in model:
            raise ValueError("model %r is an interpolated (blend) entry: attach the leaves to its source coupons" % name)
        record = {"Model": name}
        # r1: the twin's basis file byte-identical to the model's BasisPoints; the counts.
        model_points = resolve(model["BasisPoints"], library_dir)
        twin_points = entry["VacuumBasisPoints"]
        identical = sha256_of(model_points) == sha256_of(twin_points)
        with open(model_points) as handle:
            basis_count = sum(1 for line in handle if line.strip() and not line.lower().startswith("x"))
        references = model.get("ConductorReferences") or [[0.0, 0.0, 0.0]]
        conductor_states = len(references) - 1
        record["BasisPointsIdentical"] = identical
        record["BasisPointCount"] = basis_count
        record["ConductorStateCount"] = conductor_states
        if not identical:
            raise ValueError("model %r: the vacuum twin's basis_points.csv is not byte-identical to the model's "
                             "BasisPoints (r1): the twin was not generated on the model's hats" % name)
        # The matrices: the twin's, and (r3) the regenerated electric pair vs the model of record.
        vacuum_thin = read_domain_response_matrix(entry["VacuumThin"])
        vacuum_fabricated = read_domain_response_matrix(entry["VacuumFabricated"])
        expected = basis_count + conductor_states
        if vacuum_thin.shape != (expected, expected) or vacuum_fabricated.shape != (expected, expected):
            raise ValueError("model %r: the vacuum matrices are %s / %s, the basis is %d" %
                             (name, vacuum_thin.shape, vacuum_fabricated.shape, expected))
        deviation = -1.0
        if entry.get("RegeneratedElectricThin") or entry.get("RegeneratedElectricFabricated"):
            if not (entry.get("RegeneratedElectricThin") and entry.get("RegeneratedElectricFabricated")):
                raise ValueError("model %r: RegeneratedElectricThin / Fabricated must be given together" % name)
            record_thin = read_domain_response_matrix(resolve(model["ThinMatrix"], library_dir))
            record_fab = read_domain_response_matrix(resolve(model["FabricatedMatrix"], library_dir))
            deviation = max(relative_deviation(read_domain_response_matrix(entry["RegeneratedElectricThin"]), record_thin),
                            relative_deviation(read_domain_response_matrix(entry["RegeneratedElectricFabricated"]),
                                               record_fab))
        record["ElectricRegenerationMaxRelativeDeviation"] = deviation
        psd_ratio, negative = psd_statistic(vacuum_fabricated - vacuum_thin)
        record["PsdRatio"] = psd_ratio
        record["NegativeEigenvalues"] = negative
        record["Verdict"] = ("Refused r3" if deviation < 0.0 or deviation > rules["SameMeshTolerance"] else
                             "Refused r4" if psd_ratio > rules["PsdRefusalRatio"] else
                             "Finding" if psd_ratio > rules["PsdFindingRatio"] else "Applied")
        # Copy the twin files under the new library.
        model_dir = os.path.join(output_dir, "models", name)
        os.makedirs(model_dir, exist_ok=True)
        thin_out = os.path.join(model_dir, "magnetic-thin-domain-response-matrix.csv")
        fab_out = os.path.join(model_dir, "magnetic-fabricated-domain-response-matrix.csv")
        shutil.copyfile(entry["VacuumThin"], thin_out)
        shutil.copyfile(entry["VacuumFabricated"], fab_out)
        source = {
            "Version": MAGNETIC_SOURCE_VERSION,
            "Generator": "attach_magnetic_leaves.py",
            "Command": entry.get("Command", ""),
            "Toolchain": entry.get("Toolchain", ""),
            "VacuumPermittivity": float(entry.get("VacuumPermittivity", 1.0)),
            "Fabrication": fabrication,
            "ElectricRegenerationMaxRelativeDeviation": deviation,
            "BasisPointCount": basis_count,
            "ConductorStateCount": conductor_states,
            "BasisPointsIdentical": identical,
            "PsdRatio": psd_ratio,
            "NegativeEigenvalues": negative,
            "MeshSha256": entry.get("MeshSha256", {}),
            "MeshOfRecordSha256Match": entry.get("MeshOfRecordSha256Match", {}),
            "VacuumThinSha256": sha256_of(thin_out),
            "VacuumFabricatedSha256": sha256_of(fab_out),
        }
        if now:
            source["Date"] = now
        target = out_models[name]
        target["MagneticFabricatedMatrix"] = os.path.abspath(fab_out)
        target["MagneticThinMatrix"] = os.path.abspath(thin_out)
        target["MagneticSource"] = source
        records.append(record)
    out["MagneticRules"] = dict(rules)
    out["Name"] = library.get("Name", "library") + "-mag"
    out["DerivedFrom"] = {"Library": library.get("Name"), "Tool": "attach_magnetic_leaves.py",
                          "Models": [r["Model"] for r in records]}
    return out, records


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--library", required=True, help="the input Version-3 process library JSON (read only)")
    parser.add_argument("--attach", required=True, help="the attachment manifest JSON (see the module docstring)")
    parser.add_argument("--output", required=True, help="the output library JSON (a NEW file; its directory holds models/)")
    parser.add_argument("--same-mesh-tolerance", type=float, default=DEFAULT_RULES["SameMeshTolerance"],
                        help="MagneticRules.SameMeshTolerance (r3; the Phase-0 calibration 1e-12)")
    parser.add_argument("--psd-finding-ratio", type=float, default=DEFAULT_RULES["PsdFindingRatio"],
                        help="MagneticRules.PsdFindingRatio (r4 FINDING above it; 1e-2)")
    parser.add_argument("--psd-refusal-ratio", type=float, default=DEFAULT_RULES["PsdRefusalRatio"],
                        help="MagneticRules.PsdRefusalRatio (r4 REFUSAL above it; provisionally 0.1, ERRATUM 1227)")
    parser.add_argument("--records", help="write the per-model records JSON here (default: beside --output)")
    args = parser.parse_args(argv)
    if os.path.exists(args.output):
        parser.error("--output %s exists: a new library is a NEW file" % args.output)
    with open(args.library) as handle:
        library = json.load(handle)
    with open(args.attach) as handle:
        attachments = json.load(handle)
    rules = {"SameMeshTolerance": args.same_mesh_tolerance, "PsdFindingRatio": args.psd_finding_ratio,
             "PsdRefusalRatio": args.psd_refusal_ratio}
    if not rules["PsdRefusalRatio"] >= rules["PsdFindingRatio"] >= 0.0 or rules["SameMeshTolerance"] < 0.0:
        parser.error("the rules must be nonnegative with PsdRefusalRatio >= PsdFindingRatio")
    output_dir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(output_dir, exist_ok=True)
    out, records = attach(library, os.path.dirname(os.path.abspath(args.library)), attachments, output_dir, rules)
    with open(args.output, "w") as handle:
        json.dump(out, handle, indent=2)
        handle.write("\n")
    records_path = args.records or os.path.splitext(args.output)[0] + "-magnetic-records.json"
    with open(records_path, "w") as handle:
        json.dump({"Rules": rules, "Models": records}, handle, indent=2)
        handle.write("\n")
    for record in records:
        print("%s: %s (psd ratio %.3e, %d negative, same-mesh deviation %.3e)" %
              (record["Model"], record["Verdict"], record["PsdRatio"], record["NegativeEigenvalues"],
               record["ElectricRegenerationMaxRelativeDeviation"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
