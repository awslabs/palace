# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""surface-response-requirements.json: canonical digest, set diff, aggregate reconciliation.

The preflight manifest aggregates requirements by their canonical geometry (Topology,
Geometry, Interfaces, BoundaryCondition, Dimension) with a Count and a TotalEdgeLength; it
carries no per-segment assignment and no feature positions. The canonical form used here
drops everything that depends on the library or on the run (Status, SelectedModels,
NormalizedLibraryDistance, Reason, Library.Path/Name, Statistics, Summary, Complete) so
that two manifests of the same geometry hash identically whatever library matched them
(invariant A3), and re-aggregates records that differ only in those fields.
"""

import hashlib
import json
import re

# Longitudinal (per-length) feature classes; the Curved* classes are the curved-edge chain
# rule's counterparts (bend radius below the recorded straight threshold).
CURVED_PAIRS = ("CurvedSameConductorGap", "CurvedDifferentConductorGap", "CurvedSameConductorStrip")
TRANSLATIONAL = ("IsolatedEdge", "SameConductorGap", "DifferentConductorGap", "SameConductorStrip", "ParallelEdgeCluster", "CurvedEdge") + CURVED_PAIRS
VERTEX = ("ConvexCorner", "ConcaveCorner", "Endpoint", "Junction")
CLUSTER = ("SpatialEdgeCluster",)
RUN_DEPENDENT_FIELDS = ("Status", "SelectedModels", "NormalizedLibraryDistance", "Reason", "Count", "TotalEdgeLength")


def load_manifest(path):
    with open(path) as source:
        return json.load(source)


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def requirement_key(requirement):
    """Geometry-only identity of a requirement record."""
    return canonical_json({k: v for k, v in requirement.items() if k not in RUN_DEPENDENT_FIELDS})


def canonical_requirements(manifest, geometry_only_counts=False):
    """key -> {Count, TotalEdgeLength, Topology} after dropping the run-dependent fields.

    With geometry_only_counts the Count of translational topologies is dropped (it is the
    number of mesh segments, which changes under refinement while the geometry does not),
    so the digest can be compared across meshes of the same layout (invariant A5)."""
    result = {}
    for requirement in manifest["Requirements"]:
        key = requirement_key(requirement)
        entry = result.setdefault(key, {"Count": 0, "TotalEdgeLength": 0.0, "Topology": requirement["Topology"]})
        entry["Count"] += int(requirement["Count"])
        entry["TotalEdgeLength"] += float(requirement["TotalEdgeLength"])
    if geometry_only_counts:
        for entry in result.values():
            if entry["Topology"] in TRANSLATIONAL:
                entry["Count"] = None
    return result


def canonical_digest(manifest, geometry_only_counts=False, length_decimals=6):
    entries = []
    for key, entry in sorted(canonical_requirements(manifest, geometry_only_counts).items()):
        entries.append([key, entry["Count"], round(entry["TotalEdgeLength"], length_decimals)])
    payload = canonical_json(
        {
            "MatchingRadius": manifest.get("Library", {}).get("MatchingRadius"),
            "DecisionQuantization": manifest.get("Library", {}).get("DecisionQuantization"),
            "MeshDimension": manifest.get("MeshDimension"),
            "Maxwell": manifest.get("Maxwell"),
            "Requirements": entries,
        }
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def diff_manifests(a, b, geometry_only_counts=False, length_tolerance=1.0e-6):
    """Requirement set difference between two manifests by geometry key."""
    ca = canonical_requirements(a, geometry_only_counts)
    cb = canonical_requirements(b, geometry_only_counts)
    added = [{"Key": json.loads(k), **cb[k]} for k in sorted(set(cb) - set(ca))]
    removed = [{"Key": json.loads(k), **ca[k]} for k in sorted(set(ca) - set(cb))]
    changed = []
    for k in sorted(set(ca) & set(cb)):
        if ca[k]["Count"] != cb[k]["Count"] or abs(ca[k]["TotalEdgeLength"] - cb[k]["TotalEdgeLength"]) > length_tolerance * max(1.0, abs(ca[k]["TotalEdgeLength"])):
            changed.append({"Key": json.loads(k), "A": ca[k], "B": cb[k]})
    return {"Added": added, "Removed": removed, "Changed": changed, "Identical": not (added or removed or changed)}


def summarize(manifest):
    """Counts and lengths by topology class, weight sums, matching-radius, statistics."""
    by_topology = {}
    weight_defects = []
    for requirement in manifest["Requirements"]:
        topology = requirement["Topology"]
        entry = by_topology.setdefault(topology, {"Records": 0, "Count": 0, "TotalEdgeLength": 0.0, "Missing": 0, "MissingLength": 0.0})
        entry["Records"] += 1
        entry["Count"] += int(requirement["Count"])
        entry["TotalEdgeLength"] += float(requirement["TotalEdgeLength"])
        if requirement.get("Status") == "Missing":
            entry["Missing"] += int(requirement["Count"])
            entry["MissingLength"] += float(requirement["TotalEdgeLength"])
        models = requirement.get("SelectedModels")
        if models:
            total = sum(float(m["Weight"]) for m in models)
            if abs(total - 1.0) > 1.0e-9:
                weight_defects.append({"Topology": topology, "WeightSum": total, "Models": [m["Name"] for m in models]})
    translational_length = sum(v["TotalEdgeLength"] for k, v in by_topology.items() if k in TRANSLATIONAL)
    translational_count = sum(v["Count"] for k, v in by_topology.items() if k in TRANSLATIONAL)
    vertex_count = sum(v["Count"] for k, v in by_topology.items() if k in VERTEX)
    cluster_records = [r for r in manifest["Requirements"] if r["Topology"] in CLUSTER]
    clusters = [
        {
            "EdgeCount": r["Geometry"].get("EdgeCount"),
            "Status": r.get("Status"),
            "Count": r["Count"],
            "Reason": r.get("Reason", ""),
            "Conductors": sorted({e.get("Conductor") for e in r["Geometry"].get("Edges", []) if isinstance(e, dict)}),
        }
        for r in cluster_records
    ]
    return {
        "ByTopology": dict(sorted(by_topology.items())),
        "TranslationalLength": translational_length,
        "TranslationalCount": translational_count,
        "VertexFeatureCount": vertex_count,
        "Clusters": clusters,
        "WeightDefects": weight_defects,
        "MatchingRadius": manifest.get("Library", {}).get("MatchingRadius"),
        "Statistics": manifest.get("Statistics", {}).get("Geometry", {}),
        "Summary": manifest.get("Summary", {}),
    }


LOG_PATTERNS = {
    "OmittedCrossInterface": re.compile(r"Omitting (\d+) of (\d+) three-dimensional target edge segments which are within 2R"),
    "OmittedUnsupported": re.compile(
        r"Omitting (\d+) of (\d+) three-dimensional target edge segments in unsupported local interaction neighborhoods "
        r"\(nonparallel: (\d+), incompatible process normal: (\d+), process-normal offset: (\d+), unclassified topology: (\d+), "
        r"missing library model: (\d+), multi-edge: (\d+)\)"
    ),
    "MatchedSegments": re.compile(r"Matched physical edge segments: (\d+)"),
    "CornerPatches": re.compile(r"Matched corner patches: (\d+)"),
    "EndpointPatches": re.compile(r"Matched endpoint patches: (\d+)"),
    "JunctionPatches": re.compile(r"Matched junction patches: (\d+)"),
    "ClusterPatches": re.compile(r"Matched spatial edge-cluster patches: (\d+)"),
    "UnmatchedVertices": re.compile(r"has (\d+) unmatched corner, endpoint, or junction vertices"),
    "Bisection": re.compile(r"Added (\d+) elements in (\d+) iterations of local bisection"),
    "Ranks": re.compile(r"Running with (\d+) MPI process"),
}


def parse_palace_log(path):
    """Counts the preflight prints but does not write to the manifest."""
    with open(path, errors="replace") as source:
        text = re.sub(r"\x1b\[[0-9;]*m", "", source.read())
    result = {"OmittedCrossInterface": [], "OmittedUnsupported": []}
    for name, pattern in LOG_PATTERNS.items():
        for match in pattern.finditer(text):
            values = [int(v) for v in match.groups()]
            if name == "OmittedCrossInterface":
                result[name].append({"Omitted": values[0], "Of": values[1]})
            elif name == "OmittedUnsupported":
                result[name].append(
                    dict(
                        zip(
                            ("Omitted", "Of", "Nonparallel", "IncompatibleProcessNormal", "ProcessNormalOffset", "UnclassifiedTopology", "MissingLibraryModel", "MultiEdge"),
                            values,
                        )
                    )
                )
            elif name == "Bisection":
                result[name] = {"AddedElements": values[0], "Iterations": values[1]}
            else:
                result[name] = values[0]
    # Segments the classifier dropped from the requirement set: the cross-interface omissions
    # and the unsupported neighborhoods except the missing-library-model ones, which stay in
    # the manifest as Missing translational requirements.
    result["OmittedSegments"] = sum(e["Omitted"] for e in result["OmittedCrossInterface"]) + sum(
        e["Omitted"] - e["MissingLibraryModel"] for e in result["OmittedUnsupported"]
    )
    return result
