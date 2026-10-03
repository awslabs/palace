#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""(F) The spatial-coupon qualification upgrade (USER decision 281, decision 282 with the
ruling on the R0 review's MAJOR-3, decision 285 (5)): the dense held-out traces, their
closure and p-stability criteria, the matrix identity, the conductor-consistency gate as the
rebuild acceptance, the reference box integral hook for windows, and the library statuses
PendingQualification -> Qualified -> WindowValidated (or Failed).

(a) Dense held-out traces. A spatial coupon is solved, fabricated and thin twins, at the
library order p4 and the control order p5 under >= 1 DENSE trace t and must satisfy, per
interface class k in {SA, MS, MA (raw; sharp when the run carries radial shells), Domain}
WITHIN R:

    closure      | E_thin,p4(t) + dE_model(t) - E_fab,p5(t) | / E_fab,p5(t) <= 0.02
    p-stability  | dE_p5(t) - dE_p4(t) | / E_fab,p5(t) <= 0.02   (dE_p = E_fab,p - E_thin,p)
    identity     | E_fab,p4(t) - t^T Q_fab t | / E_fab,p4(t) <= 1e-6

with dE_model(t) = t^T (Q_fab - Q_thin) t from the p4 library matrices (`matrices.py`
format: domain-response-matrix.csv / surface-response-matrix.csv). The dense traces: T1 (a
device-derived coupon) = the thin DEVICE trace at the coupon's knots from one coarse thin
solve of the registration device / window (surface-response-traces.csv, per placed patch);
T2 (every coupon) = the conductor states and the trace of a unit line charge 5 R outside each
of the four in-plane faces at the process plane, superposed with the conductor states so that
every conductor cross-section carries its potential exactly (a dense trace consistent by
construction). The decision-66 caveat holds for MA raw: the thin twin's MA is read at its
cutoff; the sharp reading (radial shells) is reported alongside when present.

(b) The reference box integral (D2a) as the per-class closure of window validation: for a
window reference solve sub-tagged with the placed boxes, `ft_model / REF_A` within the
decision-218 markers (|.| - 1 <= 0.05 validated class / 0.10 new class) with the bracket
[E_in, E_in + E_straddle]; recorded per class and interface.

(3) The conductor-consistency gate (decision 277) is the acceptance check of every rebuilt
coupon: on the registration device the model's plane conductor vertices lie on real device
metal of the same conductor, so Count must be 0 and MaxRatio <= 1e-6.

Statuses: PendingQualification (today's p-sequence controls) -> Qualified (a) for T1 or T2
AND the identity AND gate Count 0) -> WindowValidated ((b) places the class inside its marker
on >= 1 window reference); (a) failing -> Failed (never placed).
"""
import argparse
import csv
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
TOOLS = HERE.parent
for path in (str(HERE), str(TOOLS)):
    if path not in sys.path:
        sys.path.insert(0, path)
import generate_spatial_response as producer  # noqa: E402
import trace_basis  # noqa: E402

CLOSURE_TOLERANCE = 0.02
P_STABILITY_TOLERANCE = 0.02
MATRIX_IDENTITY_TOLERANCE = 1.0e-6
GATE_MAX_RATIO = 1.0e-6
REFERENCE_MARKER_VALIDATED = 0.05
REFERENCE_MARKER_NEW_CLASS = 0.10
LINE_CHARGE_DISTANCE_OVER_R = 5.0

STATUS_PENDING = "PendingQualification"
STATUS_QUALIFIED = "Qualified"
STATUS_WINDOW_VALIDATED = "WindowValidated"
STATUS_FAILED = "Failed"
STATUSES = (STATUS_PENDING, STATUS_QUALIFIED, STATUS_WINDOW_VALIDATED, STATUS_FAILED)

DOMAIN_CLASS = "Domain"
TRACE_FAMILIES = ("T1", "T2")


class SpatialQualificationError(ValueError):
    """The qualification inputs are inconsistent (the reason is named)."""


# --------------------------------------------------------------------------------------
# Criteria (pure functions on per-class energies; every energy a float, classes the keys)
# --------------------------------------------------------------------------------------


def _relative(numerator, denominator):
    if denominator == 0.0:
        return math.inf if numerator != 0.0 else 0.0
    return abs(numerator) / abs(denominator)


def dense_closure(e_thin_p4, de_model, e_fab_p5, tolerance=CLOSURE_TOLERANCE):
    """Per class: residual |E_thin,p4 + dE_model - E_fab,p5| / E_fab,p5 and its verdict."""
    result = {}
    for cls in e_fab_p5:
        residual = _relative(e_thin_p4[cls] + de_model[cls] - e_fab_p5[cls], e_fab_p5[cls])
        result[cls] = {"Residual": residual, "Tolerance": tolerance, "Passed": residual <= tolerance}
    return result


def p_stability(de_p5, de_p4, e_fab_p5, tolerance=P_STABILITY_TOLERANCE):
    """Per class: |dE_p5 - dE_p4| / E_fab,p5 (the correction's p-stability on the dense
    trace; the ruling on the R0 review's MAJOR-3) and its verdict."""
    result = {}
    for cls in e_fab_p5:
        residual = _relative(de_p5[cls] - de_p4[cls], e_fab_p5[cls])
        result[cls] = {"Residual": residual, "Tolerance": tolerance, "Passed": residual <= tolerance}
    return result


def quadratic_form(matrix, coefficients):
    """t^T Q t for a symmetric response matrix (dense numpy or {(i, j): Q_ij} with i <= j,
    one-based indices) and the trace coefficient vector t."""
    t = np.asarray(coefficients, dtype=float)
    if isinstance(matrix, dict):
        total = 0.0
        for (i, j), value in matrix.items():
            if i < 1 or j < 1 or i > len(t) or j > len(t):
                raise SpatialQualificationError(f"matrix entry ({i}, {j}) outside the {len(t)} trace coefficients")
            total += value * t[i - 1] * t[j - 1] * (1.0 if i == j else 2.0)
        return float(total)
    q = np.asarray(matrix, dtype=float)
    if q.shape != (len(t), len(t)):
        raise SpatialQualificationError(f"matrix shape {q.shape} does not match {len(t)} trace coefficients")
    return float(t @ q @ t)


def matrix_identity(e_fab_p4, q_fab, coefficients, tolerance=MATRIX_IDENTITY_TOLERANCE):
    """Per class: |E_fab,p4(t) - t^T Q_fab t| / E_fab,p4(t) (the export / lift conventions,
    incl. the sqrt(Z0) coefficient scale of D2a section 5; D3-C read <= 8.7e-8)."""
    result = {}
    for cls, energy in e_fab_p4.items():
        if cls not in q_fab:
            raise SpatialQualificationError(f"no fabricated matrix for class {cls}")
        predicted = quadratic_form(q_fab[cls], coefficients)
        residual = _relative(energy - predicted, energy)
        result[cls] = {"Energy": energy, "Predicted": predicted, "Residual": residual, "Tolerance": tolerance,
                       "Passed": residual <= tolerance}
    return result


def conductor_consistency_gate(diagnostics, model_name, tolerance=GATE_MAX_RATIO):
    """The decision-277 gate as the rebuild acceptance: from a solve's
    SurfaceResponse.Diagnostics.ConductorConsistency, the model's records must all be probed
    (not Excluded) with MaxRatio <= tolerance and the model must contribute no count."""
    records = [r for r in diagnostics.get("Records", []) if r.get("Model") == model_name]
    if not records:
        return {"Passed": False, "Probed": 0, "Reason": f"no conductor-consistency record for model {model_name} "
                                                       "(the model was not placed on the registration device, or the "
                                                       "solve carried no probes: untestable, never qualified)"}
    max_ratio = max(float(r.get("MaxRatio", math.inf)) for r in records)
    excluded = sum(1 for r in records if r.get("Excluded"))
    count = sum(1 for r in records if float(r.get("MaxRatio", math.inf)) > float(diagnostics.get("Tolerance", tolerance)))
    passed = excluded == 0 and max_ratio <= tolerance and count == 0 and int(diagnostics.get("Count", count)) >= count
    return {"Passed": passed, "Probed": len(records), "Excluded": excluded, "Count": count, "MaxRatio": max_ratio,
            "Tolerance": tolerance}


def reference_box_closure(ft_model, e_in, e_straddle, validated_class):
    """(b): ft_model / REF_A with REF_A the central reading E_in + E_straddle / 2 and the
    bracket [E_in, E_in + E_straddle]; the decision-218 marker (0.05 validated / 0.10 new)."""
    reference = e_in + 0.5 * e_straddle
    marker = REFERENCE_MARKER_VALIDATED if validated_class else REFERENCE_MARKER_NEW_CLASS
    ratio = ft_model / reference if reference != 0.0 else math.inf
    bracket = sorted([ft_model / e_in if e_in != 0.0 else math.inf,
                      ft_model / (e_in + e_straddle) if e_in + e_straddle != 0.0 else math.inf])
    return {"Ratio": ratio, "Bracket": bracket, "Marker": marker, "Reference": reference,
            "Passed": abs(ratio - 1.0) <= marker}


def qualification_status(current, *, dense_passed, identity_passed, gate_passed, window_validated=False):
    """The library status transition: PendingQualification stays until (a) + identity + gate
    pass (Qualified); (b) on >= 1 window reference lifts Qualified to WindowValidated; a failing
    (a) or identity is Failed (never placed: the thin-run guard treats it as Missing); a
    failing gate on a (B) coupon is a mis-keyed library or a placement defect: Failed."""
    if current not in STATUSES:
        raise SpatialQualificationError(f"unknown library status {current!r}")
    if not dense_passed or not identity_passed or not gate_passed:
        return STATUS_FAILED
    if window_validated:
        return STATUS_WINDOW_VALIDATED
    return STATUS_QUALIFIED


# --------------------------------------------------------------------------------------
# Palace outputs: within-R energies per class of a dense-trace run; the p4 matrices
# --------------------------------------------------------------------------------------


def _rows(path):
    with open(path, newline="") as stream:
        return [{key.strip(): value.strip() for key, value in row.items()} for row in csv.DictReader(stream)]


def interface_classes(config):
    """{interface index: class} from a run config's Dielectric entries: the type (MA, MS, SA);
    an MA entry of a radial-shell expansion (index >= its parent stride) reads "MA sharp"."""
    classes = {}
    for entry in config["Boundaries"]["Postprocessing"]["Dielectric"]:
        index = int(entry["Index"])
        cls = str(entry.get("Type", "")).upper()
        if entry.get("RadialShell") is not None or entry.get("Shell") is not None:
            cls = "MA sharp"
        classes[index] = cls
    return classes


def within_r_energies(postpro, classes):
    """Per source i: {class: E within R (J)} from surface-Q.csv (participations), domain-E.csv
    (E_elec) and surface-Q-edge.csv (E_out at the largest R): E_in = p_surf E_elec - E_out,
    summed over the interfaces of a class; Domain = E_elec."""
    postpro = Path(postpro)
    domain = {int(float(row["i"])): float(row["E_elec (J)"]) for row in _rows(postpro / "domain-E.csv")}
    surface = {}
    for row in _rows(postpro / "surface-Q.csv"):
        i = int(float(row["i"]))
        for key, value in row.items():
            if key.startswith("p_surf["):
                interface = int(key[len("p_surf["):-1])
                surface[(i, interface)] = float(value) * domain[i]
    outside = {}
    for row in _rows(postpro / "surface-Q-edge.csv"):
        key = (int(float(row["i"])), int(float(row["interface"])))
        radius = float(row["R (m)"])
        if key not in outside or radius > outside[key][0]:
            outside[key] = (radius, float(row["E_out (J)"]))
    energies = {}
    for i, e_elec in domain.items():
        per_class = {DOMAIN_CLASS: e_elec}
        for (source, interface), e_surf in surface.items():
            if source != i:
                continue
            cls = classes.get(interface)
            if cls is None:
                raise SpatialQualificationError(f"surface-Q.csv interface {interface} has no class in the config")
            if (i, interface) not in outside:
                raise SpatialQualificationError(f"surface-Q-edge.csv carries no E_out for source {i}, interface {interface}")
            per_class[cls] = per_class.get(cls, 0.0) + e_surf - outside[(i, interface)][1]
        energies[i] = per_class
    return energies


def response_matrices(postpro, classes):
    """{class: {(i, j): Q_ij}} within R from the basis run's domain-response-matrix.csv and
    surface-response-matrix.csv (the whole-interface group at the largest radius), summed
    over the interfaces of a class."""
    postpro = Path(postpro)
    matrices = {DOMAIN_CLASS: {}}
    for row in _rows(postpro / "domain-response-matrix.csv"):
        i, j = int(float(row["basis_i"])), int(float(row["basis_j"]))
        matrices[DOMAIN_CLASS][(min(i, j), max(i, j))] = float(row["Q_ij (J)"])
    groups = {}
    for row in _rows(postpro / "surface-response-matrix.csv"):
        interface = int(float(row["interface"]))
        edge, radius = int(float(row["edge"])), float(row["R (m)"])
        i, j = int(float(row["basis_i"])), int(float(row["basis_j"]))
        current = groups.get(interface)
        if current is None or (edge, -radius) < (current[0], -current[1]):
            groups[interface] = (edge, radius)
    for row in _rows(postpro / "surface-response-matrix.csv"):
        interface = int(float(row["interface"]))
        if (int(float(row["edge"])), float(row["R (m)"])) != groups[interface]:
            continue
        cls = classes.get(interface)
        if cls is None:
            raise SpatialQualificationError(f"surface-response-matrix.csv interface {interface} has no class in the config")
        i, j = int(float(row["basis_i"])), int(float(row["basis_j"]))
        key = (min(i, j), max(i, j))
        matrices.setdefault(cls, {})
        matrices[cls][key] = matrices[cls].get(key, 0.0) + float(row["Q_ij (J)"])
    return matrices


def matrix_difference(fabricated, thin):
    """{class: {(i, j): Q_fab - Q_thin}} over the classes of both."""
    result = {}
    for cls in fabricated:
        if cls not in thin:
            raise SpatialQualificationError(f"the thin matrices carry no class {cls}")
        keys = set(fabricated[cls]) | set(thin[cls])
        result[cls] = {key: fabricated[cls].get(key, 0.0) - thin[cls].get(key, 0.0) for key in keys}
    return result


# --------------------------------------------------------------------------------------
# Dense traces: T1 (the device trace) and T2 (the synthetic family); their trace files
# --------------------------------------------------------------------------------------


def device_traces(traces_csv, model_index):
    """T1: {patch: coefficient vector} of a device run's surface-response-traces.csv for the
    model (columns i, patch, model, coefficient, conductor state, value): the basis
    coefficients in order followed by the conductor states, one vector per placed patch of
    the first excitation."""
    rows = _rows(traces_csv)
    per_patch = {}
    for row in rows:
        if int(float(row["i"])) != 1 or int(float(row["model"])) != int(model_index):
            continue
        patch = int(float(row["patch"]))
        coefficient = int(float(row["coefficient"]))
        state = int(float(row["conductor state"]))
        per_patch.setdefault(patch, ([], []))
        (per_patch[patch][1] if state > 0 else per_patch[patch][0]).append(
            (state if state > 0 else coefficient, float(row["value (V)"])))
    traces = {}
    for patch, (coefficients, states) in per_patch.items():
        coefficients.sort()
        states.sort()
        if [k for k, _ in coefficients] != list(range(1, len(coefficients) + 1)):
            raise SpatialQualificationError(f"patch {patch}: the trace coefficients are not 1..N")
        traces[patch] = [v for _, v in coefficients] + [v for _, v in states]
    return traces


def synthetic_traces(basis, labels, radius, conductor_references, distance_over_r=LINE_CHARGE_DISTANCE_OVER_R):
    """T2: the stated family as coefficient vectors [c_1..c_N, s_2..s_C] over the basis
    knots (active basis vertices) and the conductor states: the C - 1 conductor states
    (unit potential on conductor c, zero elsewhere) and, for each of the four in-plane faces,
    the potential of a unit line charge (phi = -ln r, normalised to unit range over the
    knots) placed `distance_over_r` R outside the face at the process plane with every
    conductor GROUNDED (states 0: the conductor cross-sections carry their potential exactly
    — a consistent dense trace with the strong edge field of an exterior influence on
    grounded metal). The DESIGN's variant "superposed with the states read at the conductor
    references" needs conductor potentials other than 0 / 1 V in one excitation, which
    Palace's PrescribedPotential TerminalAttributes (held at one volt) cannot impose for two
    or more conductor states at once: recorded as a limitation (a TerminalPotential field of
    Palace would lift it); `conductor_references` is kept for that extension."""
    points, index_of = np.asarray(basis["Points"]), np.asarray(basis["Basis"])
    lower, upper = np.asarray(basis["Lower"], dtype=float), np.asarray(basis["Upper"], dtype=float)
    count = int(index_of.max())
    conductors = sorted(set(int(label) for label in labels) - {0})
    states = [c for c in conductors if c != 1]
    knots = np.zeros((count, 3))
    for k in range(1, count + 1):
        knots[k - 1] = points[index_of == k][0]
    traces = []
    for state in states:
        traces.append({"Name": f"state-{state}", "Family": "T2",
                       "Coefficients": [0.0] * count + [1.0 if c == state else 0.0 for c in states]})
    plane_z = 0.5 * (lower[2] + upper[2])
    faces = (("x0", np.array([lower[0] - distance_over_r * radius, 0.5 * (lower[1] + upper[1]), plane_z])),
             ("x1", np.array([upper[0] + distance_over_r * radius, 0.5 * (lower[1] + upper[1]), plane_z])),
             ("y0", np.array([0.5 * (lower[0] + upper[0]), lower[1] - distance_over_r * radius, plane_z])),
             ("y1", np.array([0.5 * (lower[0] + upper[0]), upper[1] + distance_over_r * radius, plane_z])))
    references = {int(c): np.asarray(r, dtype=float) for c, r in conductor_references.items()}
    for face, source in faces:
        phi = -np.log(np.linalg.norm((knots - source)[:, :2], axis=1))
        scale = float(np.ptp(phi)) or 1.0
        values = (phi - phi.min()) / scale
        traces.append({"Name": f"line-charge-{face}", "Family": "T2",
                       "Coefficients": values.tolist() + [0.0] * len(states)})
    return traces


def representable_trace(coefficients, states):
    """A dense trace as one Palace excitation: PrescribedPotential imposes the DataFile on
    the matching surface and holds TerminalAttributes at ONE volt, so a trace may carry at
    most one nonzero conductor state s_c: the excitation is the trace scaled by 1 / s_c with
    conductor c as the terminal, its energies scale back by s_c^2 (quadratic). Returns
    (scaled coefficients, terminal conductor or None, energy scale); a trace with two or
    more distinct nonzero states fails closed (recorded, never approximated)."""
    count = len(coefficients) - len(states)
    nonzero = [(c, coefficients[count + n]) for n, c in enumerate(states) if coefficients[count + n] != 0.0]
    if not nonzero:
        return list(coefficients), None, 1.0
    if len(nonzero) > 1:
        raise SpatialQualificationError(f"a dense trace with the conductor states {dict(nonzero)} cannot be imposed in one "
                                        "excitation (TerminalAttributes hold one volt): unsupported, recorded")
    conductor, value = nonzero[0]
    return [v / value for v in coefficients], conductor, value * value


def write_dense_trace(path, basis, labels, coefficients):
    """A dense trace as a PrescribedPotential DataFile: sum_k c_k hat_k + sum_c s_c lift_c on
    the trace mesh vertices (the hats and lifts of case_inputs.regenerate_traces), scaled to
    its representable excitation (representable_trace); returns (path, terminal conductor,
    energy scale)."""
    points, triangles, index_of = basis["Points"], basis["Triangles"], np.asarray(basis["Basis"])
    count = int(index_of.max())
    states = sorted(set(int(label) for label in labels) - {0, 1})
    if len(coefficients) != count + len(states):
        raise SpatialQualificationError(f"a dense trace needs {count} coefficients + {len(states)} states, got "
                                        f"{len(coefficients)}")
    scaled, terminal, energy_scale = representable_trace(coefficients, states)
    values = np.zeros(len(points))
    for k in range(1, count + 1):
        values[index_of == k] = scaled[k - 1]
    for n, conductor in enumerate(states):
        values[np.asarray(labels) == conductor] = scaled[count + n]
    producer.write_surface_trace(Path(path), points, triangles, values)
    return path, terminal, energy_scale


def dense_trace_config(reference_config, mesh, output, traces, order):
    """A solve config of the dense traces on `mesh` at `order`: the case's run config with
    the PrescribedPotential sources replaced by the traces (DataFile per trace; the terminal
    conductor's attributes from the run config's own terminal source), the response matrix
    off (energies per source are the output). `traces`: [{Path, Terminal}]."""
    config = json.loads(json.dumps(reference_config))
    config["Model"]["Mesh"] = str(mesh)
    config["Problem"]["Output"] = str(output)
    config["Solver"]["Order"] = int(order)
    entries = config["Boundaries"]["PrescribedPotential"]
    template = entries[0]
    terminals = [e["TerminalAttributes"] for e in entries if e.get("TerminalAttributes")]
    sources = []
    for index, trace in enumerate(traces, start=1):
        entry = {key: value for key, value in template.items() if key not in ("TerminalAttributes",)}
        entry["Index"] = index
        entry["DataFile"] = str(trace["Path"])
        if trace.get("Terminal") is not None:
            position = int(trace["Terminal"]) - 2  # conductor 2 is the first terminal source
            if position < 0 or position >= len(terminals):
                raise SpatialQualificationError(f"the run config has no terminal source for conductor {trace['Terminal']}")
            entry["TerminalAttributes"] = terminals[position]
        sources.append(entry)
    config["Boundaries"]["PrescribedPotential"] = sources
    electrostatic = config["Solver"].get("Electrostatic", {})
    for key in ("ResponseMatrix", "AggregateResponseMatrix", "ResponseMatrixInterfaces", "ResponseMatrixEdges"):
        electrostatic.pop(key, None)
    return config


# --------------------------------------------------------------------------------------
# The evaluation record
# --------------------------------------------------------------------------------------


def evaluate_trace(name, family, coefficients, *, fab_p4, thin_p4, fab_p5, thin_p5, q_fab, q_thin):
    """The three criteria of one dense trace from its per-class energies at p4 / p5 and the
    p4 matrices."""
    de_p4 = {cls: fab_p4[cls] - thin_p4[cls] for cls in fab_p4}
    de_p5 = {cls: fab_p5[cls] - thin_p5[cls] for cls in fab_p5}
    difference = matrix_difference(q_fab, q_thin)
    de_model = {cls: quadratic_form(difference[cls], coefficients) for cls in fab_p4}
    closure = dense_closure(thin_p4, de_model, fab_p5)
    stability = p_stability(de_p5, de_p4, fab_p5)
    identity = matrix_identity(fab_p4, q_fab, coefficients)
    return {"Name": name, "Family": family,
            "Energies": {"FabricatedP4": fab_p4, "ThinP4": thin_p4, "FabricatedP5": fab_p5, "ThinP5": thin_p5,
                         "ModelCorrection": de_model},
            "Closure": closure, "PStability": stability, "MatrixIdentity": identity,
            "Passed": all(v["Passed"] for v in closure.values()) and all(v["Passed"] for v in stability.values())
            and all(v["Passed"] for v in identity.values())}


def evaluate(traces, *, gate, current_status=STATUS_PENDING, reference_boxes=None, tolerances=None):
    """The qualification record of a coupon: every dense trace's criteria (evaluate_trace
    results), the gate, the optional window references ((b) per class) and the status."""
    if not traces:
        raise SpatialQualificationError("(a) needs at least one dense trace")
    families = {trace["Family"] for trace in traces}
    dense_passed = all(trace["Passed"] for trace in traces)
    identity_passed = all(v["Passed"] for trace in traces for v in trace["MatrixIdentity"].values())
    window = None
    if reference_boxes:
        window = {"Passed": all(entry["Passed"] for entry in reference_boxes), "Entries": reference_boxes}
    status = qualification_status(current_status, dense_passed=dense_passed, identity_passed=identity_passed,
                                  gate_passed=bool(gate["Passed"]), window_validated=bool(window and window["Passed"]))
    return {"Version": 1, "Rule": "decision 282 (F): closure |E_thin,p4 + dE_model - E_fab,p5| / E_fab,p5 <= 0.02 and "
                                  "p-stability |dE_p5 - dE_p4| / E_fab,p5 <= 0.02 per interface class within R on every "
                                  "dense trace, the matrix identity <= 1e-6, the conductor-consistency gate Count 0 "
                                  "(decision 277); (b) the reference box integral within the decision-218 marker lifts "
                                  "Qualified to WindowValidated",
            "Tolerances": {"Closure": CLOSURE_TOLERANCE, "PStability": P_STABILITY_TOLERANCE,
                           "MatrixIdentity": MATRIX_IDENTITY_TOLERANCE, "GateMaxRatio": GATE_MAX_RATIO,
                           **(tolerances or {})},
            "Families": sorted(families), "Traces": traces, "Gate": gate, "WindowReference": window,
            "DensePassed": dense_passed, "IdentityPassed": identity_passed, "GatePassed": bool(gate["Passed"]),
            "PreviousStatus": current_status, "Status": status}


def stamp_library_status(library_path, model_name, record):
    """Write the model's QualificationStatus and the (F) record into process-library.json."""
    library_path = Path(library_path)
    library = json.loads(library_path.read_text())
    models = [m for m in library.get("Models", []) if m.get("Name") == model_name]
    if len(models) != 1:
        raise SpatialQualificationError(f"{library_path}: model {model_name} not found exactly once")
    models[0]["QualificationStatus"] = record["Status"]
    models[0]["SpatialQualification"] = {key: record[key] for key in
                                         ("Version", "Families", "DensePassed", "IdentityPassed", "GatePassed", "Status")}
    models[0]["LibraryQualified"] = record["Status"] in (STATUS_QUALIFIED, STATUS_WINDOW_VALIDATED)
    library_path.write_text(json.dumps(library, indent=2) + "\n")
    return library


# --------------------------------------------------------------------------------------
# Command line
# --------------------------------------------------------------------------------------


def load_basis(source_directory):
    """The trace basis (mesh frame), its conductor labels and the model of a coupon source
    directory (process-library.json, basis-contract.json, trace-vertices.csv,
    trace-triangles.csv)."""
    source_directory = Path(source_directory)
    library = json.loads((source_directory / "process-library.json").read_text())
    model = library["Models"][0]
    basis = trace_basis.load_trace_basis(source_directory / "basis-contract.json", source_directory / "trace-vertices.csv",
                                         source_directory / "trace-triangles.csv", source_directory / "process-library.json")
    with open(source_directory / "trace-vertices.csv", newline="") as stream:
        labels = np.array([int(row["conductor"]) for row in csv.DictReader(stream)])
    references = {}
    if "ConductorReferences" in model:
        for n, reference in enumerate(model["ConductorReferences"], start=1):
            references[n] = basis["Frame"] @ np.asarray(reference, dtype=float)
    elif "Reference" in model:
        references[1] = basis["Frame"] @ np.asarray(model["Reference"], dtype=float)
    return basis, labels, model, float(library["MatchingRadius"]), references


def command_traces(args):
    """Write the dense traces (T2 always; T1 from --device-traces) of a coupon source
    directory and the fabricated / thin solve configs at the control orders."""
    basis, labels, model, radius, references = load_basis(args.source)
    out = Path(args.output)
    (out / "traces").mkdir(parents=True, exist_ok=True)
    traces = synthetic_traces(basis, labels, radius, references)
    if args.device_traces:
        for patch, coefficients in sorted(device_traces(args.device_traces, args.model_index).items()):
            traces.append({"Name": f"device-patch-{patch}", "Family": "T1", "Coefficients": coefficients})
    written = []
    for trace in traces:
        path = out / "traces" / f"{trace['Name']}.csv"
        try:
            _, terminal, energy_scale = write_dense_trace(path, basis, labels, trace["Coefficients"])
        except SpatialQualificationError as error:
            trace["Unsupported"] = str(error)
            continue
        trace.update({"Path": str(path), "Terminal": terminal, "EnergyScale": energy_scale})
        written.append(trace)
    unsupported = [t for t in traces if "Unsupported" in t]
    traces = written
    configs = {}
    for kind, mesh, run_config in (("fabricated", args.fabricated_mesh, args.run_config),
                                   ("thin", args.thin_mesh, args.thin_run_config)):
        if mesh is None:
            continue
        if run_config is None:
            raise SpatialQualificationError(f"the {kind} mesh needs its run config (case_inputs.derive of the {kind} case)")
        reference = json.loads(Path(run_config).read_text())
        for order in args.orders:
            name = f"{kind}-p{order}"
            config = dense_trace_config(reference, mesh, str(out / name / "postpro"), traces, order)
            (out / name).mkdir(parents=True, exist_ok=True)
            (out / name / "config.json").write_text(json.dumps(config, indent=2) + "\n")
            configs[name] = str(out / name / "config.json")
    (out / "dense-traces.json").write_text(json.dumps({"Version": 1, "Model": model["Name"], "Traces": traces,
                                                       "Unsupported": unsupported, "Configs": configs,
                                                       "Orders": list(args.orders)}, indent=2) + "\n")
    print(f"{len(traces)} dense traces ({', '.join(sorted({t['Family'] for t in traces}))}; {len(unsupported)} "
          f"unsupported), {len(configs)} configs -> {out}")
    return 0


def command_evaluate(args):
    """Evaluate (a) + identity + gate [+ (b)] from the run outputs and stamp the status."""
    manifest = json.loads(Path(args.dense_traces).read_text())
    configs = {name: json.loads(Path(path).read_text()) for name, path in manifest["Configs"].items()}
    classes = interface_classes(configs[f"fabricated-p{args.library_order}"])
    runs = {}
    for kind in ("fabricated", "thin"):
        for order in (args.library_order, args.control_order):
            name = f"{kind}-p{order}"
            runs[name] = within_r_energies(Path(configs[name]["Problem"]["Output"]), classes)
    q_fab = response_matrices(args.fabricated_matrices, classes)
    q_thin = response_matrices(args.thin_matrices, classes)
    traces = []
    for index, trace in enumerate(manifest["Traces"], start=1):
        scale = float(trace.get("EnergyScale", 1.0))

        def scaled(run):
            return {cls: scale * value for cls, value in runs[run][index].items()}
        traces.append(evaluate_trace(trace["Name"], trace["Family"], trace["Coefficients"],
                                     fab_p4=scaled(f"fabricated-p{args.library_order}"),
                                     thin_p4=scaled(f"thin-p{args.library_order}"),
                                     fab_p5=scaled(f"fabricated-p{args.control_order}"),
                                     thin_p5=scaled(f"thin-p{args.control_order}"),
                                     q_fab=q_fab, q_thin=q_thin))
    gate_source = json.loads(Path(args.gate).read_text())
    diagnostics = gate_source.get("SurfaceResponse", {}).get("Diagnostics", {}).get("ConductorConsistency")
    if diagnostics is None:
        diagnostics = gate_source.get("Summary", {}).get("ConductorConsistency")
    if diagnostics is None:
        raise SpatialQualificationError(f"{args.gate} carries no ConductorConsistency record")
    gate = conductor_consistency_gate(diagnostics, manifest["Model"])
    reference_boxes = None
    if args.reference_box:
        reference_boxes = [reference_box_closure(e["FtModel"], e["EIn"], e["EStraddle"], e.get("ValidatedClass", False))
                           | {"Class": e.get("Class"), "Interface": e.get("Interface"), "Window": e.get("Window")}
                           for e in json.loads(Path(args.reference_box).read_text())["Entries"]]
    library = json.loads(Path(args.library).read_text())
    model = [m for m in library["Models"] if m["Name"] == manifest["Model"]][0]
    record = evaluate(traces, gate=gate, current_status=model.get("QualificationStatus", STATUS_PENDING),
                      reference_boxes=reference_boxes)
    Path(args.output).write_text(json.dumps(record, indent=2) + "\n")
    if not args.dry_run:
        stamp_library_status(args.library, manifest["Model"], record)
    print(f"{manifest['Model']}: dense {record['DensePassed']} identity {record['IdentityPassed']} gate "
          f"{record['GatePassed']} -> {record['Status']}")
    return 0 if record["Status"] in (STATUS_QUALIFIED, STATUS_WINDOW_VALIDATED) else 1


def add_arguments(parser):
    commands = parser.add_subparsers(dest="spatial_command", required=True)
    traces = commands.add_parser("traces", help="write the dense held-out traces and the twin configs")
    traces.add_argument("--source", required=True, help="the coupon source directory (process-library.json, trace basis)")
    traces.add_argument("--run-config", required=True, help="the fabricated case's run config (case_inputs.derive)")
    traces.add_argument("--thin-run-config", help="the thin case's run config (case_inputs.derive of the thin case)")
    traces.add_argument("--fabricated-mesh", help="the fabricated identity mesh")
    traces.add_argument("--thin-mesh", help="the thin identity mesh")
    traces.add_argument("--orders", type=int, nargs="+", default=[4, 5], help="library and control orders (default 4 5)")
    traces.add_argument("--device-traces", help="a device run's surface-response-traces.csv (T1)")
    traces.add_argument("--model-index", type=int, default=1, help="the model index of the coupon in that run")
    traces.add_argument("--output", required=True)
    traces.set_defaults(func=command_traces)
    evaluate_parser = commands.add_parser("evaluate", help="evaluate (a), the identity, the gate and (b); stamp the status")
    evaluate_parser.add_argument("--dense-traces", required=True, help="dense-traces.json written by `traces`")
    evaluate_parser.add_argument("--library-order", type=int, default=4)
    evaluate_parser.add_argument("--control-order", type=int, default=5)
    evaluate_parser.add_argument("--fabricated-matrices", required=True, help="postpro of the fabricated basis run (p4)")
    evaluate_parser.add_argument("--thin-matrices", required=True, help="postpro of the thin basis run (p4)")
    evaluate_parser.add_argument("--gate", required=True, help="palace.json (or the preflight manifest) with the "
                                                               "ConductorConsistency record of the registration device")
    evaluate_parser.add_argument("--reference-box", help="JSON {Entries: [{Class, Interface, Window, FtModel, EIn, "
                                                         "EStraddle, ValidatedClass}]} from the sub-tagged window reference")
    evaluate_parser.add_argument("--library", required=True, help="process-library.json to stamp")
    evaluate_parser.add_argument("--output", required=True, help="spatial-qualification.json")
    evaluate_parser.add_argument("--dry-run", action="store_true")
    evaluate_parser.set_defaults(func=command_evaluate)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(parser)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
