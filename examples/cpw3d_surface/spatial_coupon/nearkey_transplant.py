#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The near-key TRANSPLANT (DESIGN v2 section 4.1 / 4.2; the nearkey-sensitivity lane's method of
STAGE2-PLAN AMENDMENT 7 A7.5 and decision 380 made production code): a donor model's response
matrices expressed on the exact feature's own trace basis (the `--basis-only` generator output),
with the self-tests T1-T5, the energy-form fidelity T2e and the knot accounting.

  P ((N_D + S) x (N_E + S)): every FREE donor knot is mapped onto the exact box by the face-wise affine
  map D-box -> E-box (per in-plane axis x' = c_E + (x - c_D) L_E / L_D; z unchanged), located in the
  exact trace triangulation and given its barycentric weights on the exact free knots (an exact
  conductor vertex contributes its conductor's state column; conductor 1 = 0 by the (F) convention);
  a mapped knot within 4.5 quanta of 1e-6 R (8.55 pm at R = 1.9 um; decision 380 (3)(A)) of an exact
  vertex IS that vertex; the S state rows are the identity.  Q_reused = P^T Q_D P for the fabricated
  domain, fabricated surface (every interface, every value column), thin domain, thin surface AND the
  fabricated SHELLED surface matrix (decision 422 / 431), written in the exact CSV form.

Gates (DESIGN 4.2; the tolerances from the rule file): T1 constant trace (<= 1e-12); T2 the five (F)
synthetic traces at the exact knots mapped by P vs the analytic values at the donor knots (<= 1e-3 of
the range; default reuse refuses a failure, fallback admits it with T2e <= 2e-3); T2e the same traces in
ENERGY form on the donor matrices (<= 2e-3; 2 x max T2e = the predictor's transplant term); T3 the
round trip E -> D -> E on the synthetic traces (information); T4 state-2 -> state-2 (<= 1e-12 relative to
the donor's own state-2 energy and, when its (F) record is readable, to the stored Predicted value); T5
symmetry / PSD of every transplanted matrix; knots: mutual-nearest matching of the mapped donor knots
and the exact knots - orphans (<= 2 in each direction) and the matched displacement (<= 0.055 R).
"""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

LINE_CHARGE_DISTANCE_OVER_R = 5.0   # qualify/spatial_qualification.LINE_CHARGE_DISTANCE_OVER_R
SNAP_QUANTA, QUANTUM_OVER_R = 4.5, 1e-6
SYNTHETIC_TRACES = ("state-2", "line-charge-x0", "line-charge-x1", "line-charge-y0", "line-charge-y1")
LINE_CHARGES = SYNTHETIC_TRACES[1:]
# the matrix roles of a model entry: (role, library key, file name)
MATRIX_ROLES = (("fab_domain", "FabricatedMatrix", "fabricated-domain-response-matrix.csv"),
                ("fab_surface", "FabricatedSurfaceMatrix", "fabricated-surface-response-matrix.csv"),
                ("thin_domain", "ThinMatrix", "thin-domain-response-matrix.csv"),
                ("thin_surface", "ThinSurfaceMatrix", "thin-surface-response-matrix.csv"),
                ("fab_surface_shelled", "FabricatedSurfaceMatrixShelled", "fabricated-surface-response-matrix-shelled.csv"))
SURFACE_ROLES = ("fab_surface", "thin_surface", "fab_surface_shelled")
# the plain surface matrices carry the three coupon interfaces (Interfaces[].Coupon -> Type)
CLASSES = {1: "MA", 2: "MS", 3: "SA"}
DEFAULT_GATES = {"T1": 1e-12, "T2Default": 1e-3, "T2eMax": 2e-3, "T3Information": 1e-3, "T4": 1e-12, "T5PSD": 1e-9,
                 "T5Symmetry": 1e-12}


class TransplantError(ValueError):
    """A basis / matrix the transplant cannot use (fail closed)."""


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path):
    with open(path, newline="") as handle:
        return [{k.strip(): v.strip() for k, v in r.items()} for r in csv.DictReader(handle)]


def quad(q, t):
    return float(t @ q @ t)


class TraceBasis:
    """A model's trace mesh and basis (trace-vertices.csv, trace-triangles.csv, basis-points.csv) in
    one-based basis order; the matrices are attached by `load_matrices` (a donor) or absent (the exact
    feature's generated basis)."""

    def __init__(self, directory, name, *, support_box=None, interfaces=None):
        self.directory = Path(directory)
        self.name = name
        vertices = rows(self.directory / "trace-vertices.csv")
        self.points = np.array([[float(r[k]) for k in "xyz"] for r in vertices])
        self.basis = np.array([int(r["basis"]) for r in vertices])
        self.conductor = np.array([int(r["conductor"]) for r in vertices])
        self.triangles = np.array([[int(r[k]) - 1 for k in ("vertex_i", "vertex_j", "vertex_k")]
                                   for r in rows(self.directory / "trace-triangles.csv")])
        basis_points = np.array([[float(r[k]) for k in "xyz"] for r in rows(self.directory / "basis-points.csv")])
        self.n_free = int(self.basis.max())
        if (self.basis > 0).sum() != self.n_free or len(basis_points) != self.n_free:
            raise TransplantError(f"{name}: basis columns {(self.basis > 0).sum()} / max {self.n_free} / basis-points "
                                  f"{len(basis_points)} disagree")
        free = self.points[self.basis > 0][np.argsort(self.basis[self.basis > 0])]
        if not np.array_equal(free, basis_points):
            raise TransplantError(f"{name}: basis-points.csv != the free trace vertices in basis order")
        if np.any((self.conductor > 0) & (self.basis > 0)):
            raise TransplantError(f"{name}: a conductor vertex carries a basis column")
        self.conductors = sorted(set(int(c) for c in self.conductor) - {0})
        self.states = [c for c in self.conductors if c != 1]
        if self.states != [2]:
            raise TransplantError(f"{name}: the rule's class has states [2], this basis has {self.states}")
        self.size = self.n_free + len(self.states)
        self.lower, self.upper = self.points.min(0), self.points.max(0)
        if support_box is not None:
            if not (np.allclose(self.lower[:2], support_box[:2]) and np.allclose(self.upper[:2], support_box[2:])):
                raise TransplantError(f"{name}: the trace box {self.lower[:2]} {self.upper[:2]} != SupportBox {support_box}")
        if interfaces is not None:
            coupon_classes = {int(i["Coupon"]): i["Type"] for i in interfaces}
            if coupon_classes != CLASSES:
                raise TransplantError(f"{name}: Interfaces {coupon_classes} are not the coupon classes {CLASSES}")
        self.face_of_triangle = self._triangle_faces()
        self.matrices, self.surface_meta, self.matrix_files, self.domain_header = {}, {}, {}, None

    def free_knots(self):
        order = np.argsort(self.basis[self.basis > 0])
        return self.points[self.basis > 0][order], self.conductor[self.basis > 0][order]

    def load_matrices(self, paths, required=("fab_domain", "fab_surface", "thin_domain", "thin_surface")):
        """`paths` = {role: file path}; the required roles must be present (a missing shelled matrix is
        reported by the caller's policy: absent here, refused by the reuse in Default / Fallback)."""
        for role in required:
            if role not in paths or paths[role] is None or not Path(paths[role]).is_file():
                raise TransplantError(f"{self.name}: the {role} matrix file is not readable ({paths.get(role)})")
        for role, path in paths.items():
            if path is None or not Path(path).is_file():
                continue
            self.matrix_files[role] = str(path)
            if role.endswith("domain"):
                header, matrix = self.read_domain(path)
                self.matrices[role] = {"Q_ij (J)": matrix}
                self.domain_header = self.domain_header or header
            else:
                self.matrices[role], self.surface_meta[role] = self.read_surface(path)

    def read_domain(self, path):
        n = self.size
        q = np.full((n, n), np.nan)
        for r in rows(path):
            i, j = int(float(r["basis_i"])) - 1, int(float(r["basis_j"])) - 1
            if j < i:
                raise TransplantError(f"{path}: not an upper triangle")
            q[i, j] = q[j, i] = float(r["Q_ij (J)"])
        if not np.all(np.isfinite(q)):
            raise TransplantError(f"{path}: incomplete upper triangle for a basis of size {n}")
        with open(path) as handle:
            header = handle.readline().rstrip("\n")
        return header, q

    def read_surface(self, path):
        """{interface: {column: matrix}} for every interface group of the file (3 for the coupon surface
        matrices, every radial shell for the shelled one) and the writer's metadata."""
        n = self.size
        data = rows(path)
        columns = [c for c in data[0] if c.startswith("Q_")]
        groups = sorted(set((int(float(r["interface"])), int(float(r["edge"])), r["R (m)"]) for r in data))
        if any(g[1] != 1 for g in groups) or len(set(g[2] for g in groups)) != 1:
            raise TransplantError(f"{path}: the surface groups {groups} are not one edge / one R")
        out = {g[0]: {c: np.full((n, n), np.nan) for c in columns} for g in groups}
        for r in data:
            k = int(float(r["interface"]))
            i, j = int(float(r["basis_i"])) - 1, int(float(r["basis_j"])) - 1
            if j < i:
                raise TransplantError(f"{path}: not an upper triangle")
            for c in columns:
                out[k][c][i, j] = out[k][c][j, i] = float(r[c])
        for k in out:
            for c in columns:
                if not np.all(np.isfinite(out[k][c])):
                    raise TransplantError(f"{path}: interface {k} column {c} incomplete for a basis of size {n}")
        with open(path) as handle:
            header = handle.readline().rstrip("\n")
        return out, {"columns": columns, "R (m)": groups[0][2], "edge": 1, "header": header, "interfaces": sorted(out)}

    def _triangle_faces(self):
        scale = float(np.max(self.upper - self.lower))
        faces = []
        for tri in self.triangles:
            xyz = self.points[tri]
            face = None
            for axis in range(3):
                for side, bound in ((0, self.lower[axis]), (1, self.upper[axis])):
                    if np.all(np.abs(xyz[:, axis] - bound) <= 1e-9 * scale):
                        face = (axis, side)
            if face is None:
                raise TransplantError(f"{self.name}: a trace triangle lies in no box face plane")
            faces.append(face)
        return faces

    def locate(self, point):
        """(triangle index, barycentric weights) of a point on the box surface: the triangle with the
        largest minimum weight (None when the point is on no face)."""
        scale = float(np.max(self.upper - self.lower))
        best = None
        for t, (axis, side) in enumerate(self.face_of_triangle):
            bound = self.lower[axis] if side == 0 else self.upper[axis]
            if abs(point[axis] - bound) > 1e-9 * scale:
                continue
            keep = [a for a in range(3) if a != axis]
            p = point[keep]
            a, b, c = (self.points[v][keep] for v in self.triangles[t])
            det = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])
            w1 = ((p[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (p[1] - a[1])) / det
            w2 = ((b[0] - a[0]) * (p[1] - a[1]) - (p[0] - a[0]) * (b[1] - a[1])) / det
            w = np.array([1.0 - w1 - w2, w1, w2])
            if w.min() >= -1e-9 and (best is None or w.min() > best[1].min()):
                best = (t, w)
        return best


def model_matrix_paths(entry, root=None):
    """{role: path} of a library model entry; relative paths resolved against `root` (the library's
    directory); a missing role maps to None."""
    paths = {}
    for role, key, _ in MATRIX_ROLES:
        value = entry.get(key)
        if value is None:
            paths[role] = None
            continue
        path = Path(value)
        if not path.is_absolute() and root is not None:
            path = Path(root) / path
        paths[role] = str(path)
    return paths


def affine_map(donor, exact):
    c_d, c_e = 0.5 * (donor.lower + donor.upper), 0.5 * (exact.lower + exact.upper)
    l_d, l_e = donor.upper - donor.lower, exact.upper - exact.lower
    mapped = donor.points.copy()
    for axis in range(2):
        mapped[:, axis] = c_e[axis] + (donor.points[:, axis] - c_d[axis]) * (l_e[axis] / l_d[axis])
    if not (np.allclose(donor.lower[2], exact.lower[2]) and np.allclose(donor.upper[2], exact.upper[2])):
        raise TransplantError("the donor and exact boxes differ in their z extents")
    return mapped, {"DonorCentre": c_d.tolist(), "ExactCentre": c_e.tolist(), "LengthRatio": (l_e / l_d).tolist(),
                    "DonorBox": [donor.lower.tolist(), donor.upper.tolist()], "ExactBox": [exact.lower.tolist(), exact.upper.tolist()]}


def inverse_affine_points(points, affine):
    inverse = points.copy()
    for axis in range(2):
        inverse[:, axis] = affine["DonorCentre"][axis] + (points[:, axis] - affine["ExactCentre"][axis]) / affine["LengthRatio"][axis]
    return inverse


def build_map(source, target, mapped_points, radius):
    """P ((N_src + S) x (N_tgt + S)): the source free knots (at mapped_points) interpolated on the
    target triangulation; returns (P, per-row record)."""
    snap_um = SNAP_QUANTA * QUANTUM_OVER_R * radius
    P = np.zeros((source.size, target.size))
    record = []
    state_column = {2: target.n_free}
    for v in np.where(source.basis > 0)[0]:
        a = source.basis[v] - 1
        located = target.locate(mapped_points[v])
        if located is None:
            raise TransplantError(f"{source.name}: mapped knot {int(a + 1)} at {mapped_points[v].tolist()} lies on no target face")
        t, w = located
        distances = np.linalg.norm(target.points[target.triangles[t]] - mapped_points[v], axis=1)
        if distances.min() <= snap_um:
            w = np.where(distances == distances.min(), 1.0, 0.0)
        touched_conductor = []
        for vertex, weight in zip(target.triangles[t], w):
            if weight == 0.0:
                continue
            if target.basis[vertex] > 0:
                P[a, target.basis[vertex] - 1] += weight
            else:
                touched_conductor.append((int(target.conductor[vertex]), float(weight)))
                if target.conductor[vertex] in state_column:
                    P[a, state_column[target.conductor[vertex]]] += weight
        nearest = int(np.argmin(np.linalg.norm(target.points - mapped_points[v], axis=1)))
        record.append({"SourceKnot": int(a + 1), "Triangle": int(t + 1), "Weights": [float(x) for x in w],
                       "Snapped": bool(distances.min() <= snap_um), "Vertices": [int(x + 1) for x in target.triangles[t]],
                       "NearestTargetVertexDistanceUm": float(np.linalg.norm(target.points[nearest] - mapped_points[v])),
                       "NearestTargetBasis": int(target.basis[nearest]), "Exact": bool(np.max(w) == 1.0),
                       "ConductorVertexWeights": touched_conductor})
    for s, column in state_column.items():
        P[source.n_free + source.states.index(s), column] = 1.0
    return P, record


def congruence(matrices, P):
    """P^T Q P for every (role, interface, column)."""
    out = {}
    for role, value in matrices.items():
        if role.endswith("domain"):
            out[role] = {c: P.T @ q @ P for c, q in value.items()}
        else:
            out[role] = {k: {c: P.T @ q @ P for c, q in columns.items()} for k, columns in value.items()}
    return out


def analytic_traces(exact, radius, points, labels, reference_knots):
    """The (F) T2 family (spatial_qualification.synthetic_traces) at `points`: state-2 and the four
    line charges 5 R outside the exact box's in-plane faces, unit-range normalised over the exact knots."""
    lower, upper = exact.lower, exact.upper
    plane_z = 0.5 * (lower[2] + upper[2])
    faces = (("x0", np.array([lower[0] - LINE_CHARGE_DISTANCE_OVER_R * radius, 0.5 * (lower[1] + upper[1]), plane_z])),
             ("x1", np.array([upper[0] + LINE_CHARGE_DISTANCE_OVER_R * radius, 0.5 * (lower[1] + upper[1]), plane_z])),
             ("y0", np.array([0.5 * (lower[0] + upper[0]), lower[1] - LINE_CHARGE_DISTANCE_OVER_R * radius, plane_z])),
             ("y1", np.array([0.5 * (lower[0] + upper[0]), upper[1] + LINE_CHARGE_DISTANCE_OVER_R * radius, plane_z])))
    out = {"state-2": np.where(labels == 2, 1.0, 0.0)}
    for face, source in faces:
        phi_ref = -np.log(np.linalg.norm((reference_knots - source)[:, :2], axis=1))
        scale = float(np.ptp(phi_ref)) or 1.0
        phi = -np.log(np.linalg.norm((points - source)[:, :2], axis=1))
        values = (phi - phi_ref.min()) / scale
        values[labels > 0] = 0.0
        out[f"line-charge-{face}"] = values
    return out


def with_state(values, name):
    return np.concatenate([values, [1.0 if name == "state-2" else 0.0]])


def knot_matching(mapped_donor_knots, exact_knots, radius):
    """Mutual-nearest matching of the mapped donor free knots and the exact free knots: matched pairs
    with their displacement, donor orphans (no exact counterpart), exact orphans (no donor counterpart)."""
    d = np.linalg.norm(mapped_donor_knots[:, None, :] - exact_knots[None, :, :], axis=2)
    nearest_exact = np.argmin(d, axis=1)
    nearest_donor = np.argmin(d, axis=0)
    matched = [(int(a), int(nearest_exact[a]), float(d[a, nearest_exact[a]])) for a in range(len(mapped_donor_knots))
               if nearest_donor[nearest_exact[a]] == a]
    matched_donor = {a for a, _, _ in matched}
    matched_exact = {e for _, e, _ in matched}
    donor_orphans = [{"DonorKnot": int(a + 1), "NearestExactDistanceUm": float(d[a].min())} for a in range(len(mapped_donor_knots))
                     if a not in matched_donor]
    exact_orphans = [{"ExactKnot": int(e + 1), "NearestDonorDistanceUm": float(d[:, e].min())} for e in range(len(exact_knots))
                     if e not in matched_exact]
    displacement = max((x for _, _, x in matched), default=0.0)
    displaced = [{"DonorKnot": int(a + 1), "ExactKnot": int(e + 1), "DisplacementUm": x} for a, e,
                 x in matched if x > SNAP_QUANTA * QUANTUM_OVER_R * radius]
    return {"Matched": len(matched), "Displaced": displaced, "DonorOrphans": donor_orphans, "ExactOrphans": exact_orphans,
            "MaxMatchedDisplacementUm": displacement, "MaxMatchedDisplacementOverR": displacement / radius}


def transplant(exact, donor, radius, *, gates=DEFAULT_GATES, domain_limits=None, donor_stored_state2=None):
    """The full transplant of `donor` (matrices loaded) onto `exact` (basis only): returns {"P", "Mapped",
    "Affine", "Map", "Reused", "Tests", "Knots", "T2eMax", "GatesPassed", "T2Passed"}. `domain_limits` =
    the rule's Policy.Domain (orphan cap, knot displacement); `donor_stored_state2` = the donor's (F)
    record MatrixIdentity Predicted energies per class ({"SA": .., "MS": .., "MA": .., "Domain": ..}), optional."""
    limits = domain_limits or {"OrphanKnotsMax": 2, "KnotShiftMaxOverR": 0.055}
    mapped, affine = affine_map(donor, exact)
    P, map_record = build_map(donor, exact, mapped, radius)
    reused = congruence(donor.matrices, P)
    tests = {}
    # T1: the constant trace
    ones = np.ones(exact.size)
    t1 = float(np.max(np.abs(P @ ones - 1.0)))
    tests["T1"] = {"MaxError": t1, "Tolerance": gates["T1"], "Passed": bool(t1 <= gates["T1"]),
                   "RowsTouchingConductor1": [r["SourceKnot"] for r in map_record if any(c == 1 for c, _ in r["ConductorVertexWeights"])]}
    # T2 (trace form) and T2e (energy form) on the five synthetic traces
    exact_knots, exact_labels = exact.free_knots()
    donor_order = np.argsort(donor.basis[donor.basis > 0])
    donor_knots, donor_labels = mapped[donor.basis > 0][donor_order], donor.conductor[donor.basis > 0][donor_order]
    analytic_exact = analytic_traces(exact, radius, exact_knots, exact_labels, exact_knots)
    analytic_donor = analytic_traces(exact, radius, donor_knots, donor_labels, exact_knots)
    t2, t2e_per_trace = {}, {}
    for name in SYNTHETIC_TRACES:
        t_e, t_d = with_state(analytic_exact[name], name), with_state(analytic_donor[name], name)
        image = P @ t_e
        err = float(np.max(np.abs(image - t_d)))
        rng = float(np.ptp(t_e[:exact.n_free])) or 1.0
        t2[name] = {"MaxError": err, "Range": rng, "MaxErrorOverRange": err / rng, "WorstDonorRow": int(np.argmax(np.abs(image - t_d))) + 1,
                    "Passed": bool(err / rng <= gates["T2Default"])}
        energies = {}
        for k, cls in CLASSES.items():
            q = donor.matrices["fab_surface"][k]["Q_ij (J)"]
            e_d = quad(q, t_d)
            energies[cls] = abs(quad(q, image) - e_d) / e_d
        q = donor.matrices["fab_domain"]["Q_ij (J)"]
        e_d = quad(q, t_d)
        energies["Domain"] = abs(quad(q, image) - e_d) / e_d
        t2e_per_trace[name] = energies
    t2e_classes = {cls: max(t2e_per_trace[name][cls] for name in LINE_CHARGES) for cls in list(CLASSES.values()) + ["Domain"]}
    t2e_max = max(t2e_classes[cls] for cls in CLASSES.values())
    tests["T2"] = {"Traces": t2, "Tolerance": gates["T2Default"], "MaxErrorOverRange": max(x["MaxErrorOverRange"] for x in t2.values()),
                   "Passed": all(x["Passed"] for x in t2.values())}
    tests["T2e"] = {"PerTrace": t2e_per_trace, "PerClass": t2e_classes, "Max": t2e_max, "Tolerance": gates["T2eMax"],
                    "TransplantTerm": 2.0 * t2e_max, "Passed": bool(t2e_max <= gates["T2eMax"]),
                    "Rule": "max over the four line-charge traces and the classes SA / MS / MA of |E_D(P t_E) - E_D(t_D)| / E_D(t_D) "
                            "on the "
                            "donor's fabricated matrices; the predictor's transplant term = 2 x Max"}
    # T3: the round trip E -> D -> E on the synthetic line charges (information)
    P_back, _ = build_map(exact, donor, inverse_affine_points(exact.points, affine), radius)
    t3 = {}
    for name in LINE_CHARGES:
        t_e = with_state(analytic_exact[name], name)
        t3[name] = float(np.max(np.abs(P_back @ (P @ t_e) - t_e)) / (np.ptp(t_e[:exact.n_free]) or 1.0))
    tests["T3"] = {"Traces": t3, "MaxErrorOverRange": max(t3.values()), "Tolerance": gates["T3Information"],
                   "Passed": bool(max(t3.values()) <= gates["T3Information"]), "Information": True}
    # T4: state-2 -> state-2
    e_s = np.zeros(exact.size)
    e_s[exact.n_free:] = 1.0
    donor_state = np.zeros(donor.size)
    donor_state[donor.n_free:] = 1.0
    t4 = {"ImageMinusDonorState2Max": float(np.max(np.abs(P @ e_s - donor_state))), "Classes": {},
          "StoredRecordUsed": donor_stored_state2 is not None}
    worst = t4["ImageMinusDonorState2Max"]
    for label, q_reused, q_donor, stored_key in (
            [(f"fab {cls}", reused["fab_surface"][k]["Q_ij (J)"], donor.matrices["fab_surface"][k]["Q_ij (J)"], cls) for k,
             cls in CLASSES.items()]
            + [("fab Domain", reused["fab_domain"]["Q_ij (J)"], donor.matrices["fab_domain"]["Q_ij (J)"], "Domain")]):
        e_reused, e_donor = quad(q_reused, e_s), quad(q_donor, donor_state)
        entry = {"E_reused": e_reused, "E_donor_matrix": e_donor, "Relative": abs(e_reused - e_donor) / abs(e_donor)}
        worst = max(worst, entry["Relative"])
        if donor_stored_state2 is not None and stored_key in donor_stored_state2:
            stored = float(donor_stored_state2[stored_key])
            entry["E_donor_stored_predicted"] = stored
            entry["RelativeToStored"] = abs(e_reused - stored) / abs(stored)
            worst = max(worst, entry["RelativeToStored"])
        t4["Classes"][label] = entry
    t4.update({"WorstRelative": worst, "Tolerance": gates["T4"], "Passed": bool(worst <= gates["T4"])})
    tests["T4"] = t4
    # T5: symmetry and PSD of every transplanted matrix (the Q_ij (J) column; the shelled matrix per shell)
    t5 = {}
    for role, value in reused.items():
        items = [("Q_ij (J)", value["Q_ij (J)"])] if role.endswith("domain") else \
            [(f"interface {k} {CLASSES.get(k, 'shell')} Q_ij (J)", value[k]["Q_ij (J)"]) for k in value]
        for label, q in items:
            scale = float(np.max(np.abs(q)))
            sym = float(np.max(np.abs(q - q.T)) / scale) if scale > 0 else 0.0
            eig = np.linalg.eigvalsh(0.5 * (q + q.T))
            t5[f"{role} {label}"] = {"SymmetryRelative": sym, "MinEigRelative": float(eig.min() / eig.max()) if eig.max() > 0 else 0.0,
                                     "MaxEig": float(eig.max()),
                                     "Passed": bool(sym <= gates["T5Symmetry"] and eig.min() >= -gates["T5PSD"] * eig.max())}
    tests["T5"] = {"Matrices": t5, "Passed": all(x["Passed"] for x in t5.values()),
                   "Rule": "symmetric to 1e-12; min eigenvalue >= -T5PSD x max (the domain matrices PSD; the surface matrices PSD "
                           "to roundoff)"}
    # knots
    knots = knot_matching(donor_knots, exact_knots, radius)
    knots.update({"Exact": exact.n_free, "Donor": donor.n_free,
                  "KnotsInterpolated": [r["SourceKnot"] for r in map_record if not r["Exact"]],
                  "MaxNearestVertexDistanceUm": max(r["NearestTargetVertexDistanceUm"] for r in map_record),
                  "OrphanCap": limits["OrphanKnotsMax"], "KnotShiftMaxOverR": limits["KnotShiftMaxOverR"]})
    knots["Passed"] = bool(len(knots["DonorOrphans"]) <= limits["OrphanKnotsMax"] and len(knots["ExactOrphans"]) <= limits["OrphanKnotsMax"]
                           and knots["MaxMatchedDisplacementOverR"] <= limits["KnotShiftMaxOverR"])
    tests["Knots"] = knots
    gates_passed = all(tests[k]["Passed"] for k in ("T1", "T2e", "T4", "T5", "Knots"))
    return {"P": P, "Mapped": mapped, "Affine": affine,
            "Map": {"Rows": map_record, "SnapUm": SNAP_QUANTA * QUANTUM_OVER_R * radius,
                    "SnapRule": "within 4.5 quanta of 1e-6 R of an exact vertex = that vertex (decisions 303 / 317 near-match "
                                "tolerance; 380 (3)(A))"},
            "Reused": reused, "Tests": tests, "T2eMax": t2e_max, "GatesPassed": gates_passed, "T2Passed": tests["T2"]["Passed"],
            "Sizes": {"ExactFree": exact.n_free, "DonorFree": donor.n_free, "States": exact.states}}


def fmt_index(i):
    return f"{i:.2e}"


def write_domain(path, header, q):
    n = q.shape[0]
    with open(path, "w") as out:
        out.write(header + "\n")
        for i in range(n):
            for j in range(i, n):
                out.write(f" {fmt_index(i + 1)}, {fmt_index(j + 1)},        {q[i, j]:+.12e}\n")


def write_surface(path, meta, per_interface):
    n = next(iter(next(iter(per_interface.values())).values())).shape[0]
    with open(path, "w") as out:
        out.write(meta["header"] + "\n")
        for k in sorted(per_interface):
            for i in range(n):
                for j in range(i, n):
                    values = ", ".join(f"{per_interface[k][c][i, j]:+.12e}" for c in meta["columns"])
                    out.write(f" {fmt_index(k)}, {fmt_index(meta['edge'])},        {meta['R (m)']}, {fmt_index(i + 1)}, "
                              f"{fmt_index(j + 1)}, {values}\n")


def check_written(path, role, reference):
    """Re-read a written matrix file and compare with the in-memory Q_reused (the writer's round trip)."""
    back = {}
    for r in rows(path):
        i, j = int(float(r["basis_i"])) - 1, int(float(r["basis_j"])) - 1
        k = int(float(r["interface"])) if "interface" in r else 0
        back.setdefault(k, {})[(i, j)] = float(r["Q_ij (J)"])
    worst = 0.0
    for k, entries in back.items():
        q = reference["Q_ij (J)"] if role.endswith("domain") else reference[k]["Q_ij (J)"]
        scale = float(np.max(np.abs(q)))
        for (i, j), v in entries.items():
            worst = max(worst, abs(v - q[i, j]) / scale)
    return worst


def write_matrices(model_dir, donor, result):
    """The transplanted matrices (every loaded role) + interpolation-map-P.csv into `model_dir`;
    returns {role: {Path, SHA256, ReadBackMaxRelativeError}} and the map CSV record."""
    model_dir = Path(model_dir)
    model_dir.mkdir(parents=True, exist_ok=True)
    written = {}
    names = {role: name for role, _, name in MATRIX_ROLES}
    for role, value in result["Reused"].items():
        path = model_dir / names[role]
        if role.endswith("domain"):
            write_domain(path, donor.domain_header, value["Q_ij (J)"])
        else:
            write_surface(path, donor.surface_meta[role], value)
        written[role] = {"Path": str(path), "SHA256": sha256(path), "ReadBackMaxRelativeError": check_written(path, role, value)}
    map_path = model_dir / "interpolation-map-P.csv"
    np.savetxt(map_path, result["P"], delimiter=",", fmt="%.17g")
    return written, {"Path": str(map_path), "SHA256": sha256(map_path)}


def test_summary(result):
    """The TransplantTests record of DESIGN 4.3 (every gate with its values; the per-row map excluded)."""
    tests = json.loads(json.dumps(result["Tests"]))
    summary = {key: tests[key] for key in ("T1", "T2", "T2e", "T3", "T4", "T5")}
    knots = dict(tests["Knots"])
    knots["KnotsInterpolated"] = knots["KnotsInterpolated"]
    summary["Knots"] = knots
    summary["Passed"] = result["GatesPassed"] and result["T2Passed"]
    summary["GatesPassed"] = result["GatesPassed"]
    summary["T2Passed"] = result["T2Passed"]
    return summary
