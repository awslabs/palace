#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Near-key DETECTION (DESIGN v2 section 1; USER decision 428): which library model is a candidate
DONOR for a Missing SpatialEdgeCluster requirement, and the build-time metrics of the pair.

Structure key (decision 426 MAJOR-3): the runtime's cluster topology key
(signature_library.split_cluster_parameters - every Portions / Context / Box / Vertices number
nulled, entry order preserved) with the `Vertices` member REMOVED, sha256 of its canonical JSON.
RuleVersion v1 admits exactly the four calibrated families; any other key is refused
`StructureKeyNotCalibrated` (fail closed), and a donor of another key `StructureKeyMismatch`.

A library model D qualifies as a candidate donor for the requirement E iff every item of DESIGN 1.2
holds (items 0-9); the metrics W (relative lead-width change of the narrowest claim strip, donor vs
exact), S (claim-length change, kappa = 1), Vmax / Cmax / Bmax (claim-vertex / context-endpoint /
box-face displacements over R), the facing-gap scaling, the portion / context correspondence and
the Vertices-list difference are recorded whether the pair qualifies or not. Units: the signature's
(R = 1); the record quotes um / nm with the library's MatchingRadius.
"""
import hashlib
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
IDENTIFICATION = HERE.parents[1] / "surface_response_identification"
if str(IDENTIFICATION) not in sys.path:
    sys.path.insert(0, str(IDENTIFICATION))
import signature_library  # noqa: E402

STATUS_REUSED = "ReusedResponse"
REFUSAL_STRUCTURE_KEY = "StructureKeyNotCalibrated"
# one signature quantum (1e-6 R): two strip widths within it are one width
WIDTH_QUANTUM_OVER_R = signature_library.LENGTH_QUANTUM_OVER_R
# A pair whose numbers agree within the runtime's quantum near-match is ONE geometry at the grid (the
# matcher resolves it; decisions 303 / 317), not a near key of this rule.
QUANTUM_NEAR_MATCH_MAX = signature_library.CLUSTER_QUANTUM_NEAR_MATCH_MAX_QUANTA + signature_library.CLUSTER_QUANTUM_INCLUSIVE_MARGIN
PORTION_ATTRIBUTES = ("Conductor", "Gap", "Interfaces", "Law")
CONTEXT_ATTRIBUTES = ("Conductor", "Gap", "Chain", "Interfaces", "Law")


class NearKeyDetectionError(ValueError):
    """A signature the detection cannot read (fail closed; distinct from a recorded refusal)."""


def structure_key(signature):
    """(sha256 structure key, canonical JSON, sha256 of the topology key WITH Vertices)."""
    if not signature_library.is_cluster_signature(signature):
        raise NearKeyDetectionError("the structure key is defined for SpatialEdgeCluster signatures only")
    topology_json = signature_library.split_cluster_parameters(signature)[0]
    topology = json.loads(topology_json)
    topology.pop("Vertices", None)
    canonical = json.dumps(topology, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest(), canonical, hashlib.sha256(topology_json.encode()).hexdigest()


def _length(p):
    return math.dist(p[:2], p[2:])


def _direction(p):
    length = _length(p)
    if length <= 0.0:
        raise NearKeyDetectionError(f"zero-length segment {p}")
    return ((p[2] - p[0]) / length, (p[3] - p[1]) / length)


def _midpoint(p):
    return (0.5 * (p[0] + p[2]), 0.5 * (p[1] + p[3]))


def _same_direction(a, b):
    da, db = _direction(a), _direction(b)
    return abs(da[0] * db[0] + da[1] * db[1] - 1.0) <= 1e-9


def _endpoint_displacements(a, b):
    return [math.hypot(b[k] - a[k], b[k + 1] - a[k + 1]) for k in (0, 2)]


def strip_widths(portions):
    """Per portion: the metal width behind it = the separation to the nearest PARALLEL portion of the
    same conductor with the opposite gap normal and an overlapping projection (None without one).
    DESIGN 1.2 item 5 (the rule of the design lane's nearkey_pair_geometry.py)."""
    widths = []
    for i, a in enumerate(portions):
        da = _direction(a["P"])
        best = None
        for j, b in enumerate(portions):
            if i == j or b["Conductor"] != a["Conductor"]:
                continue
            if abs(b["Gap"][0] + a["Gap"][0]) > 1e-9 or abs(b["Gap"][1] + a["Gap"][1]) > 1e-9:
                continue
            db = _direction(b["P"])
            if abs(abs(da[0] * db[0] + da[1] * db[1]) - 1.0) > 1e-6:
                continue
            n = a["Gap"]
            separation = -((b["P"][0] - a["P"][0]) * n[0] + (b["P"][1] - a["P"][1]) * n[1])
            if separation <= 1e-9:
                continue
            ta = sorted((a["P"][0] * da[0] + a["P"][1] * da[1], a["P"][2] * da[0] + a["P"][3] * da[1]))
            tb = sorted((b["P"][0] * da[0] + b["P"][1] * da[1], b["P"][2] * da[0] + b["P"][3] * da[1]))
            if min(ta[1], tb[1]) - max(ta[0], tb[0]) <= 1e-9:
                continue
            if best is None or separation < best:
                best = separation
        widths.append(best)
    return widths


def facing_gap(portions):
    """The smallest separation between a portion and a parallel portion of the OTHER conductor facing
    it across the gap (DESIGN 1.2 item 6); None without a facing pair."""
    best = None
    for a in portions:
        for b in portions:
            if a["Conductor"] == b["Conductor"]:
                continue
            if abs(a["Gap"][0] + b["Gap"][0]) > 1e-9 or abs(a["Gap"][1] + b["Gap"][1]) > 1e-9:
                continue
            da, db = _direction(a["P"]), _direction(b["P"])
            if abs(abs(da[0] * db[0] + da[1] * db[1]) - 1.0) > 1e-6:
                continue
            n = a["Gap"]
            separation = (b["P"][0] - a["P"][0]) * n[0] + (b["P"][1] - a["P"][1]) * n[1]
            if separation <= 1e-9:
                continue
            ta = sorted((a["P"][0] * da[0] + a["P"][1] * da[1], a["P"][2] * da[0] + a["P"][3] * da[1]))
            tb = sorted((b["P"][0] * da[0] + b["P"][1] * da[1], b["P"][2] * da[0] + b["P"][3] * da[1]))
            if min(ta[1], tb[1]) - max(ta[0], tb[0]) <= 1e-9:
                continue
            if best is None or separation < best:
                best = separation
    return best


def correspondence(exact_entries, donor_entries, attributes, tolerance):
    """The donor index of every exact entry (DESIGN 1.2 items 2 / 3): by index when every pair has equal
    attributes, equal direction and endpoints within `tolerance`; else a bijective nearest-midpoint
    matching among the entries of equal attributes (every match within the tolerance). Returns
    (indices, method) or (None, reason)."""
    if len(exact_entries) != len(donor_entries):
        return None, f"entry counts differ ({len(exact_entries)} vs {len(donor_entries)})"

    def compatible(a, b):
        if any(a.get(key) != b.get(key) for key in attributes):
            return False
        if "Arc" in a or "Arc" in b:
            return False
        return _same_direction(a["P"], b["P"]) and max(_endpoint_displacements(a["P"], b["P"])) <= tolerance

    if all(compatible(a, b) for a, b in zip(exact_entries, donor_entries)):
        return list(range(len(exact_entries))), "ByIndex"
    taken = set()
    indices = []
    for a in exact_entries:
        candidates = [(math.dist(_midpoint(a["P"]), _midpoint(b["P"])), j) for j, b in enumerate(donor_entries)
                      if j not in taken and compatible(a, b)]
        if not candidates:
            return None, "no compatible donor entry (attributes / direction / displacement) for an exact entry"
        _, j = min(candidates)
        taken.add(j)
        indices.append(j)
    return indices, "GeometricNearestMidpoint"


def donor_status(model):
    """DESIGN 1.2 item 8: (admissible for default, admissible for fallback, reasons)."""
    reasons = []
    if model.get("QualificationStatus") != "Qualified":
        reasons.append(f"donor QualificationStatus {model.get('QualificationStatus')!r} is not Qualified")
    if model.get("StatusProvisional", False):
        reasons.append("donor is StatusProvisional")
    if model.get("QualificationStatus") == STATUS_REUSED or model.get("ReusedFrom") or model.get("TransplantedFrom"):
        reasons.append("donor is itself a reused model (no chaining)")
    if model.get("Topology") not in (None, "SpatialEdgeCluster"):
        reasons.append(f"donor Topology {model.get('Topology')!r}")
    override = bool(model.get("BuildGateOverride")) or "build-override" in str((model.get("CouponMesh") or {}).get("Path", ""))
    return not reasons and not override, not reasons, reasons, override


def analyse_pair(exact_signature, donor_model, *, rule, radius, exact_interfaces=None, exact_boundary_condition=None):
    """The near-key record of (requirement E, library model D): {"Qualifies", "DefaultAdmissible",
    "Refusals", "StructureKey", "W", "S", ...}. Every refusal of DESIGN 1.4 is a recorded reason; the
    metrics are computed as far as the pair allows."""
    domain = rule["Policy"]["Domain"]
    admissible = {entry["StructureKey"]: entry for entry in rule["AdmissibleStructureKeys"]}
    donor_signature = donor_model.get("Signature")
    record = {"Donor": donor_model.get("Name"), "Qualifies": False, "DefaultAdmissible": False, "Refusals": []}
    refusals = record["Refusals"]
    if not signature_library.is_cluster_signature(exact_signature):
        refusals.append(f"{REFUSAL_STRUCTURE_KEY}: the requirement is not a SpatialEdgeCluster signature")
        return record
    key_e, _, topology_e = structure_key(exact_signature)
    record["StructureKey"] = key_e
    record["TopologyKeyWithVertices"] = topology_e
    if key_e not in admissible:
        refusals.append(f"{REFUSAL_STRUCTURE_KEY}: {key_e[:16]}… is not one of the {len(admissible)} keys RuleVersion "
                        f"{rule['RuleVersion']} admits")
        return record
    record["StructureKeyFamily"] = admissible[key_e]["Family"]
    if not signature_library.is_cluster_signature(donor_signature):
        refusals.append("StructureKeyMismatch: the donor carries no SpatialEdgeCluster signature")
        return record
    key_d, _, topology_d = structure_key(donor_signature)
    record["DonorStructureKey"] = key_d
    if key_d != key_e:
        refusals.append(f"StructureKeyMismatch: donor {key_d[:16]}… vs requirement {key_e[:16]}…")
        return record
    # item 1: counts, Interfaces, Law, BoundaryCondition, Contract 3 (Box + Context on both)
    for name in ("EdgeCount",):
        if exact_signature.get(name) != donor_signature.get(name):
            refusals.append(f"{name} differs ({exact_signature.get(name)} vs {donor_signature.get(name)})")
    pe, pd = exact_signature.get("Portions", []), donor_signature.get("Portions", [])
    ce, cd = exact_signature.get("Context", []), donor_signature.get("Context", [])
    if len(pe) != len(pd):
        refusals.append(f"portion counts differ ({len(pe)} vs {len(pd)})")
    if len(ce) != len(cd):
        refusals.append(f"context counts differ ({len(ce)} vs {len(cd)})")
    if "Box" not in exact_signature or "Box" not in donor_signature or "Context" not in exact_signature or "Context" not in donor_signature:
        refusals.append("Contract: Box + Context (contract 3) must be present on both")
    if signature_library.conductor_count({"Type": "SpatialEdgeCluster", "Signature": exact_signature}) != \
            signature_library.conductor_count({"Type": "SpatialEdgeCluster", "Signature": donor_signature}):
        refusals.append("conductor counts differ")
    if any("Arc" in entry for entry in pe + pd + ce + cd):
        refusals.append("Arc: a curved portion or context piece (not a near key of this rule)")
    laws_e = sorted({entry.get("Law") for entry in pe + ce})
    laws_d = sorted({entry.get("Law") for entry in pd + cd})
    if laws_e != laws_d:
        refusals.append(f"Law sets differ ({laws_e} vs {laws_d})")
    interfaces_e = sorted({i for entry in pe for i in entry.get("Interfaces", [])})
    interfaces_d = sorted({i for entry in pd for i in entry.get("Interfaces", [])})
    if interfaces_e != interfaces_d:
        refusals.append(f"Interfaces differ ({interfaces_e} vs {interfaces_d})")
    if exact_interfaces is not None and donor_model.get("Interfaces") is not None:
        model_types = sorted(entry["Type"] for entry in donor_model["Interfaces"])
        required = sorted(entry["Type"] if isinstance(entry, dict) else str(entry) for entry in exact_interfaces)
        if model_types != required:
            refusals.append(f"model Interfaces {model_types} vs requirement {required}")
    if exact_boundary_condition is not None and donor_model.get("BoundaryCondition") is not None and \
            donor_model["BoundaryCondition"] != exact_boundary_condition:
        refusals.append("BoundaryCondition differs")
    if refusals:
        return record
    # quantum difference (information; a quantum near match is the runtime's exact match, not a near key)
    quantum = signature_library.cluster_quantum_difference(exact_signature, donor_signature)
    record["RuntimeTopologyKeyEqual"] = topology_e == topology_d
    record["QuantumDifference"] = None if quantum is None else {"MaxQuanta": quantum[0], "DifferingNumbers": len(quantum[1])}
    if quantum is not None and quantum[0] <= QUANTUM_NEAR_MATCH_MAX:
        refusals.append(f"QuantumNearMatch: {quantum[0]:.2f} quanta <= {QUANTUM_NEAR_MATCH_MAX} (one geometry at the grid; the "
                        f"runtime matches it exactly)")
        return record
    shift_max = float(domain["ShiftMaxOverR"])
    # items 2 / 3: correspondences
    portion_map, portion_method = correspondence(pe, pd, PORTION_ATTRIBUTES, shift_max)
    context_map, context_method = correspondence(ce, cd, CONTEXT_ATTRIBUTES, shift_max)
    if portion_map is None:
        refusals.append(f"PortionCorrespondence: {portion_method}")
    if context_map is None:
        refusals.append(f"ContextCorrespondence: {context_method}")
    if refusals:
        return record
    record["PortionCorrespondence"] = {"DonorIndexOfExactPortion": portion_map, "Method": portion_method}
    record["ContextCorrespondence"] = {"DonorIndexOfExactPiece": context_map, "Method": context_method}
    pd_ordered = [pd[j] for j in portion_map]
    cd_ordered = [cd[j] for j in context_map]
    # item 4: displacements
    v_max = max(max(_endpoint_displacements(a["P"], b["P"])) for a, b in zip(pe, pd_ordered))
    c_max = max((max(_endpoint_displacements(a["P"], b["P"])) for a, b in zip(ce, cd_ordered)), default=0.0)
    b_max = max(abs(b - a) for a, b in zip(exact_signature["Box"], donor_signature["Box"]))
    record.update({"VmaxOverR": v_max, "CmaxOverR": c_max, "BmaxOverR": b_max,
                   "VmaxNm": v_max * radius * 1e3, "CmaxNm": c_max * radius * 1e3, "BmaxNm": b_max * radius * 1e3})
    if v_max > shift_max:
        refusals.append(f"claim vertex displacement {v_max:.4f} R > {shift_max} R")
    if c_max > shift_max:
        refusals.append(f"context endpoint displacement {c_max:.4f} R > {shift_max} R")
    if b_max > float(domain["BoxShiftMaxOverR"]):
        refusals.append(f"box face displacement {b_max:.4f} R > {domain['BoxShiftMaxOverR']} R")
    # item 5: the narrowest claim strip and W
    we, wd = strip_widths(pe), strip_widths(pd_ordered)
    strips = [i for i, w in enumerate(we) if w is not None]
    if not strips:
        refusals.append("no claim strip (two parallel same-conductor portions with opposite gap normals) on the requirement")
        return record
    w_e = min(we[i] for i in strips)
    lead = [i for i in strips if abs(we[i] - w_e) <= WIDTH_QUANTUM_OVER_R]
    if any(wd[i] is None for i in lead):
        refusals.append("the narrowest strip's portions carry no strip on the donor")
        return record
    w_d = min(wd[i] for i in lead)
    if max(wd[i] for i in lead) - w_d > WIDTH_QUANTUM_OVER_R:
        refusals.append("the donor's corresponding strips are not one width")
    donor_strips = [i for i, w in enumerate(wd) if w is not None]
    if min(wd[i] for i in donor_strips) < w_d - WIDTH_QUANTUM_OVER_R:
        refusals.append("the narrowest strip is not the same set of portions on the donor")
    lengths_e = [_length(p["P"]) for p in pe]
    lengths_d = [_length(p["P"]) for p in pd_ordered]
    width_portions = sorted({i for i in range(len(pe)) if abs(lengths_e[i] - w_e) <= WIDTH_QUANTUM_OVER_R} | set(lead))
    W = (w_d - w_e) / w_e
    record["LeadWidthUm"] = {"Exact": w_e * radius, "Donor": w_d * radius, "Portions": width_portions,
                             "Rule": "the narrowest claim strip (two parallel same-conductor portions with opposite gap normals "
                                     "and overlapping projections); the donor width read on the same portions"}
    record["W"] = W
    record["W_log"] = math.log(w_d / w_e)
    if abs(W) > float(domain["WMax"]):
        refusals.append(f"|W| {abs(W):.5f} > {domain['WMax']} (the calibrated maximum)")
    # item 6: the facing gap scales with the lead
    g_e, g_d = facing_gap(pe), facing_gap(pd_ordered)
    record["JunctionGapUm"] = {"Exact": None if g_e is None else g_e * radius, "Donor": None if g_d is None else g_d * radius}
    if g_e is None or g_d is None:
        refusals.append("no facing gap between the two conductors (the g = w scaling cannot be checked)")
    else:
        gap_change = (g_d - g_e) / g_e
        record["GapScalingResidual"] = gap_change - W
        if abs(gap_change - W) > float(domain["GapScalingTolerance"]):
            refusals.append(f"the facing gap does not scale with the lead: (g_D - g_E) / g_E - W = {gap_change - W:+.5f}")
    # item 7: the claim-length change
    total = sum(lengths_e)
    dl = [b - a for a, b in zip(lengths_e, lengths_d)]
    S = sum(abs(x) for x in dl) / total
    record["S"] = S
    record["S_nonwidth"] = sum(abs(dl[i]) for i in range(len(pe)) if i not in width_portions) / total
    record["S_signed"] = sum(dl) / total
    record["ClaimsLengthUm"] = total * radius
    if S > float(domain["SMax"]):
        refusals.append(f"S {S:.5f} > {domain['SMax']}")
    # item 8: the donor's status
    default_ok, fallback_ok, status_reasons, override = donor_status(donor_model)
    record["DonorBuildGateOverride"] = override
    refusals.extend(status_reasons)
    # item 9: Vertices (derived; may differ)
    ve, vd = exact_signature.get("Vertices", []), donor_signature.get("Vertices", [])
    record["VertexListDifference"] = {"ExactVertices": len(ve), "DonorVertices": len(vd), "Differs": ve != vd}
    record["Qualifies"] = not refusals
    record["DefaultAdmissible"] = record["Qualifies"] and default_ok
    if record["Qualifies"] and not default_ok:
        record["DefaultRefusal"] = "donor built under a BuildGateOverride: fallback only (DESIGN 1.2 item 8)"
    return record


def candidate_donors(exact_signature, library, *, rule, exact_interfaces=None, exact_boundary_condition=None):
    """Every SpatialEdgeCluster model of the process library analysed against the requirement:
    {model name: record}; the qualifying ones carry Qualifies true. A requirement outside the
    admissible structure keys yields one StructureKeyNotCalibrated record under the key "*"."""
    radius = float(library["MatchingRadius"])
    if signature_library.is_cluster_signature(exact_signature):
        key, _, _ = structure_key(exact_signature)
        if key not in {entry["StructureKey"] for entry in rule["AdmissibleStructureKeys"]}:
            return {"*": {"Donor": None, "Qualifies": False, "DefaultAdmissible": False, "StructureKey": key,
                          "Refusals": [f"{REFUSAL_STRUCTURE_KEY}: {key[:16]}… is not one of the "
                                       f"{len(rule['AdmissibleStructureKeys'])} keys RuleVersion {rule['RuleVersion']} admits"]}}
    records = {}
    for model in library.get("Models", []):
        if model.get("Topology") != "SpatialEdgeCluster" or "Signature" not in model:
            continue
        records[model["Name"]] = analyse_pair(exact_signature, model, rule=rule, radius=radius, exact_interfaces=exact_interfaces,
                                              exact_boundary_condition=exact_boundary_condition)
    return records
