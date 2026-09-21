#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Calibration-only, label-only postprocess of a built coupon's identity mesh: the MA
surfaces (6000 + 100 slot + conductor: metal top face and sidewalls) are split into
radial shells following the prism tube rings (supervisor decision 56).

Shell of a surface element = the ring interval [r_{k-1}, r_k) of the tube ring set that
contains the distance of its centroid to the nearest metal edge line (the Physical
plan-view boundary segments, at the metal top z = plane + Nz MetalThickness - the
"top" edges - and at the plane - the "bottom" edges; arcs by their chords); every
element farther than the tube radius from every edge line is the "far" shell.  Each
(shell, parent MA label) pair becomes a distinct boundary label of the same interface
type (MA):

    label = SHELL_LABEL_STRIDE x ordinal + parent,  ordinal 1 = far,
    1 + k = ring k of a top edge, 1 + K + k = ring k of a bottom edge (K rings).

The relabel is byte-level on the MSH 2.2 binary file: only the physical / elementary tag
integers of the MA elements and the $PhysicalNames block change; the $Nodes block and
every other element byte are identical (asserted by re-reading both files).  The shell
areas sum to the parent MA area to roundoff (asserted against the parent's own
ownership partition record), the parent ownership is reproduced element by element
from the parent's partition certificate (asserted), and the radial partition is audited
point-wise (degree-4 positive rules: the straddling measure of elements whose quadrature
points fall in another shell is recorded per shell; the rules' weights are closed against
the elements' cross-product areas at 1e-12 - an independent area, not the same sum).

Outputs under ROOT/<case>/: identity.msh, identity.msh.radial-shells.json (the census the
qualify command consumes: label -> parent / kind / ring / radii / area) and ROOT/
library-build.json (a build record of the relabeled case: Variants.identity, the parent's
entity counts, Relabel binding), consumable by `coupon_library.py qualify`.

usage: relabel_radial_ma_shells.py CASE_ID --manifest CALIBRATION_MANIFEST
       --parent-build-record library-build.json --root DIR
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess
import sys
import time
import tomllib

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import general_mesh_manifest  # noqa: E402

SHELL_LABEL_STRIDE = 10000
MA_FAMILY = 6
FAR_ORDINAL = 1
KINDS = ("far", "top", "bottom")
CENSUS_SUFFIX = ".radial-shells.json"
STATUS_RELABELED = "relabeled"
LABEL_RULE = (f"label = {SHELL_LABEL_STRIDE} x ordinal + parent MA attribute; ordinal {FAR_ORDINAL} = farther than the "
              "tube radius from every metal edge line, 1 + k = ring k of the nearest top edge (metal top face level), "
              "1 + K + k = ring k of the nearest bottom edge (process plane), K = number of tube rings; ring k = "
              "centroid distance in [r_(k-1), r_k) with r_0 = 0")
NODE_COUNT = {1: 2, 2: 3, 3: 4, 4: 4, 5: 8, 6: 6, 7: 5, 15: 1}
SURFACE_TYPES = {2: "Triangle", 3: "Quadrangle"}
# Dunavant degree-4 triangle rule (6 positive points, weights summing to one).
TRIANGLE_RULE = [(0.445948490915965, 0.445948490915965, 0.223381589678011),
                 (0.445948490915965, 0.108103018168070, 0.223381589678011),
                 (0.108103018168070, 0.445948490915965, 0.223381589678011),
                 (0.091576213509771, 0.091576213509771, 0.109951743655322),
                 (0.091576213509771, 0.816847572980459, 0.109951743655322),
                 (0.816847572980459, 0.091576213509771, 0.109951743655322)]
GAUSS3 = [(-math.sqrt(3.0 / 5.0), 5.0 / 9.0), (0.0, 8.0 / 9.0), (math.sqrt(3.0 / 5.0), 5.0 / 9.0)]


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


class RelabelError(ValueError):
    """A fail-closed stop of the relabel (nothing is written)."""


def read_msh22_binary(data):
    """The sections of a MSH 2.2 binary file: physical names [(dim, tag, name)] with
    their byte span, the $Nodes block span and node coordinates {tag: xyz}, and the
    elements in file order: (type, tag byte offset, tags, nodes)."""
    if not data.startswith(b"$MeshFormat\n2.2 1 8\n"):
        raise RelabelError("the mesh is not a MSH 2.2 binary file (little-endian, 8-byte reals)")
    names_start = data.find(b"$PhysicalNames\n")
    names_end = data.find(b"$EndPhysicalNames\n")
    if names_start < 0 or names_end < 0:
        raise RelabelError("the mesh carries no $PhysicalNames block")
    lines = data[names_start + len(b"$PhysicalNames\n"):names_end].decode().splitlines()
    if int(lines[0]) != len(lines) - 1:
        raise RelabelError("the $PhysicalNames count differs from its entries")
    names = []
    for line in lines[1:]:
        dimension, tag, name = line.split(maxsplit=2)
        names.append((int(dimension), int(tag), name.strip('"')))
    nodes_start = data.find(b"$Nodes\n")
    header_end = data.index(b"\n", nodes_start + len(b"$Nodes\n"))
    node_count = int(data[nodes_start + len(b"$Nodes\n"):header_end])
    block = header_end + 1
    nodes_end = block + 28 * node_count
    if data[nodes_end:nodes_end + len(b"\n$EndNodes\n")] != b"\n$EndNodes\n":
        raise RelabelError("the $Nodes block does not end where its count says")
    node_records = np.frombuffer(data[block:nodes_end], dtype=np.dtype([("tag", "<i4"), ("xyz", "<f8", (3,))]))
    coordinates = {int(tag): tuple(float(v) for v in xyz) for tag, xyz in zip(node_records["tag"], node_records["xyz"])}
    elements_start = data.find(b"$Elements\n", nodes_end)
    header_end = data.index(b"\n", elements_start + len(b"$Elements\n"))
    element_count = int(data[elements_start + len(b"$Elements\n"):header_end])
    position = header_end + 1
    elements = []
    while len(elements) < element_count:
        element_type, count, tag_count = struct.unpack_from("<3i", data, position)
        position += 12
        if element_type not in NODE_COUNT:
            raise RelabelError(f"unsupported element type {element_type}")
        width = 1 + tag_count + NODE_COUNT[element_type]
        for _ in range(count):
            values = struct.unpack_from(f"<{width}i", data, position)
            elements.append((element_type, position + 4, values[1:1 + tag_count], values[1 + tag_count:]))
            position += 4 * width
    if data[position:position + len(b"\n$EndElements\n")] != b"\n$EndElements\n":
        raise RelabelError("the $Elements block does not end where its count says")
    return {"PhysicalNames": names, "PhysicalNamesSpan": (names_start, names_end + len(b"$EndPhysicalNames\n")),
            "NodesSpan": (nodes_start, nodes_end + len(b"\n$EndNodes\n")), "Coordinates": coordinates,
            "Elements": elements, "ElementsSpan": (elements_start, position + len(b"\n$EndElements\n"))}


def read_csv_rows(path):
    import csv
    with Path(path).open(newline="") as stream:
        return [{key.strip(): value.strip() for key, value in row.items()} for row in csv.DictReader(stream)]


def metal_edge_lines(boundary_rows, signature_rows, thickness):
    """The metal edge lines as plan-view segments with their z level: for every Physical
    segment of every plan-view loop, a top line at plane + Nz MetalThickness and a
    bottom line at the plane (Nz from the signature edges of the loop's conductor)."""
    normals = {}
    for row in signature_rows:
        key = (int(row["Conductor"]), float(row["Pz"]))
        normals.setdefault(key, set()).add(float(row["Nz"]))
    loops = {}
    for row in boundary_rows:
        loops.setdefault(int(row["Loop"]), []).append(row)
    lines = []
    for loop_index, rows in sorted(loops.items()):
        rows.sort(key=lambda row: int(row["Vertex"]))
        conductor = int(rows[0]["Conductor"])
        plane = float(rows[0]["Plane"])
        candidates = normals.get((conductor, plane))
        if not candidates:
            candidates = {value for (c, _), values in normals.items() if c == conductor for value in values}
        if len(candidates) != 1:
            raise RelabelError(f"loop {loop_index} (conductor {conductor}, plane {plane}) has no unique process normal "
                               f"sign in the signature: {sorted(candidates)}")
        nz = candidates.pop()
        points = [(float(row["X"]), float(row["Y"])) for row in rows]
        for index, row in enumerate(rows):
            if row["Class"] != "Physical":
                continue
            first, last = points[index], points[(index + 1) % len(points)]
            for kind, z in (("top", plane + nz * thickness), ("bottom", plane)):
                lines.append({"Kind": kind, "Loop": loop_index, "Conductor": conductor, "First": list(first),
                              "Last": list(last), "Z": z})
    if not lines:
        raise RelabelError("the plan-view boundary has no Physical segment: no metal edge line")
    return lines


def edge_distances(points, lines):
    """(distance, kind index) of every point (n x 3) to the nearest metal edge line."""
    best = np.full(len(points), np.inf)
    kind = np.zeros(len(points), dtype=int)
    for line in lines:
        a = np.asarray(line["First"])
        b = np.asarray(line["Last"])
        ab = b - a
        length2 = float(ab @ ab)
        if length2 <= 0.0:
            raise RelabelError(f"degenerate metal edge line {line}")
        t = np.clip(((points[:, :2] - a) @ ab) / length2, 0.0, 1.0)
        foot = a + t[:, None] * ab
        dxy2 = np.sum((points[:, :2] - foot) ** 2, axis=1)
        distance = np.sqrt(dxy2 + (points[:, 2] - line["Z"]) ** 2)
        closer = distance < best
        best[closer] = distance[closer]
        kind[closer] = KINDS.index(line["Kind"])
    return best, kind


def shell_ordinal(distance, kind, radii):
    """The shell ordinal of a centroid distance: far beyond the tube radius, else the
    ring of its kind."""
    rings = len(radii)
    k = int(np.searchsorted(radii, distance, side="right")) + 1
    if k > rings:
        return FAR_ORDINAL
    return 1 + k if KINDS[kind] == "top" else 1 + rings + k


def ordinal_description(ordinal, radii):
    rings = len(radii)
    if ordinal == FAR_ORDINAL:
        return "far", 0, radii[-1], math.inf
    k = ordinal - 1 if ordinal <= 1 + rings else ordinal - 1 - rings
    kind = "top" if ordinal <= 1 + rings else "bottom"
    return kind, k, (0.0 if k == 1 else radii[k - 2]), radii[k - 1]


def element_area(element_type, corners):
    """The area of a linear triangle (half the cross product) or quadrangle (the two
    triangles of its diagonal 0-2), independent of the quadrature rule."""
    corners = np.asarray(corners)
    area = 0.5 * np.linalg.norm(np.cross(corners[1] - corners[0], corners[2] - corners[0]))
    if element_type == 3:
        area += 0.5 * np.linalg.norm(np.cross(corners[2] - corners[0], corners[3] - corners[0]))
    return float(area)


def quadrature_points(element_type, corners):
    """(points, measures) of a positive degree-4 rule on a linear triangle / quadrangle."""
    corners = np.asarray(corners)
    if element_type == 2:
        area = 0.5 * np.linalg.norm(np.cross(corners[1] - corners[0], corners[2] - corners[0]))
        points = [l1 * corners[0] + l2 * corners[1] + (1.0 - l1 - l2) * corners[2] for l1, l2, _ in TRIANGLE_RULE]
        return np.asarray(points), np.asarray([w * area for _, _, w in TRIANGLE_RULE])
    points, measures = [], []
    for u, wu in GAUSS3:
        for v, wv in GAUSS3:
            shape = np.asarray([(1 - u) * (1 - v), (1 + u) * (1 - v), (1 + u) * (1 + v), (1 - u) * (1 + v)]) / 4.0
            du = np.asarray([-(1 - v), (1 - v), (1 + v), -(1 + v)]) / 4.0
            dv = np.asarray([-(1 - u), -(1 + u), (1 + u), (1 - u)]) / 4.0
            points.append(shape @ corners)
            measures.append(wu * wv * np.linalg.norm(np.cross(du @ corners, dv @ corners)))
    return np.asarray(points), np.asarray(measures)


def relabel(data, *, lines, radii):
    """The relabeled bytes and the census of the shells of every MA parent label."""
    mesh = read_msh22_binary(data)
    coordinates = mesh["Coordinates"]
    ma_labels = sorted({tag for dimension, tag, _ in mesh["PhysicalNames"] if dimension == 2 and tag // 1000 == MA_FAMILY})
    if not ma_labels:
        raise RelabelError("the mesh defines no MA boundary label (family 6000)")
    max_elementary = max(element[2][1] for element in mesh["Elements"] if len(element[2]) >= 2)
    ma_indices = [index for index, element in enumerate(mesh["Elements"])
                  if element[0] in SURFACE_TYPES and len(element[2]) >= 2 and element[2][0] in ma_labels]
    if not ma_indices:
        raise RelabelError("no surface element carries an MA label")
    corners = [np.asarray([coordinates[node] for node in mesh["Elements"][index][3]]) for index in ma_indices]
    centroids = np.asarray([c.mean(axis=0) for c in corners])
    distances, kinds = edge_distances(centroids, lines)
    ordinals = np.asarray([shell_ordinal(d, k, radii) for d, k in zip(distances, kinds)])
    groups = sorted({(int(ordinal), mesh["Elements"][index][2][0]) for ordinal, index in zip(ordinals, ma_indices)})
    elementary = {group: max_elementary + 1 + position for position, group in enumerate(groups)}
    shells = {group: {"Label": SHELL_LABEL_STRIDE * group[0] + group[1], "Parent": group[1], "Ordinal": group[0],
                      "ElementaryTag": elementary[group], "Elements": 0, "Triangles": 0, "Quadrangles": 0,
                      "Area": [], "QuadratureMeasure": [], "StraddlingMeasure": [],
                      "CentroidDistanceMinimum": math.inf, "CentroidDistanceMaximum": 0.0}
              for group in groups}
    out = bytearray(data)
    relabeled = {}
    for position, index in enumerate(ma_indices):
        element_type, tag_offset, tags, _ = mesh["Elements"][index]
        group = (int(ordinals[position]), tags[0])
        shell = shells[group]
        points, measures = quadrature_points(element_type, corners[position])
        point_distances, point_kinds = edge_distances(points, lines)
        straddling = sum(float(m) for m, d, k in zip(measures, point_distances, point_kinds)
                         if shell_ordinal(d, k, radii) != group[0])
        shell["Elements"] += 1
        shell["Triangles" if element_type == 2 else "Quadrangles"] += 1
        shell["Area"].append(element_area(element_type, corners[position]))
        shell["QuadratureMeasure"].append(float(measures.sum()))
        shell["StraddlingMeasure"].append(straddling)
        shell["CentroidDistanceMinimum"] = min(shell["CentroidDistanceMinimum"], float(distances[position]))
        shell["CentroidDistanceMaximum"] = max(shell["CentroidDistanceMaximum"], float(distances[position]))
        struct.pack_into("<2i", out, tag_offset, shell["Label"], shell["ElementaryTag"])
        relabeled[index] = (shell["Label"], shell["ElementaryTag"], tags[0])
    for shell in shells.values():
        for key in ("Area", "QuadratureMeasure", "StraddlingMeasure"):
            shell[key] = math.fsum(shell[key])
    # The parent MA names leave $PhysicalNames with their last element: a boundary
    # attribute a config may name is one the mesh carries (producer.mesh_boundary_attributes).
    names = [entry for entry in mesh["PhysicalNames"] if not (entry[0] == 2 and entry[1] in ma_labels)]
    for group, shell in shells.items():
        kind, ring, _, _ = ordinal_description(group[0], radii)
        names.append((2, shell["Label"], f"ma_shell_{kind}_{ring}_{shell['Parent']}"))
    names.sort(key=lambda entry: (entry[0], entry[1]))
    block = "$PhysicalNames\n" + str(len(names)) + "\n" + "".join(f'{d} {t} "{n}"\n' for d, t, n in names) + "$EndPhysicalNames\n"
    start, end = mesh["PhysicalNamesSpan"]
    out = bytes(out[:start]) + block.encode() + bytes(out[end:])
    return out, mesh, shells, relabeled, ma_labels


def assert_label_only(original, relabeled_bytes, mesh, relabeled):
    """Nodes byte-identical; every element identical in type / nodes / tag count, tags
    identical except the (physical, elementary) pair of the relabeled MA elements."""
    after = read_msh22_binary(relabeled_bytes)
    before_nodes = original[mesh["NodesSpan"][0]:mesh["NodesSpan"][1]]
    after_nodes = relabeled_bytes[after["NodesSpan"][0]:after["NodesSpan"][1]]
    if before_nodes != after_nodes:
        raise RelabelError("the $Nodes block changed")
    if len(after["Elements"]) != len(mesh["Elements"]):
        raise RelabelError("the element count changed")
    changed_tags = 0
    for index, (before, current) in enumerate(zip(mesh["Elements"], after["Elements"])):
        if before[0] != current[0] or before[3] != current[3] or len(before[2]) != len(current[2]):
            raise RelabelError(f"element {index} changed in type, node list or tag count")
        if index in relabeled:
            label, elementary, parent = relabeled[index]
            if current[2][:2] != (label, elementary) or before[2][0] != parent or before[2][2:] != current[2][2:]:
                raise RelabelError(f"relabeled element {index} does not carry its shell label / elementary tag only")
            changed_tags += 2
        elif before[2] != current[2]:
            raise RelabelError(f"element {index} outside the MA surfaces changed its tags")
    # Every byte outside the $PhysicalNames block and the changed tag integers is identical.
    prefix = mesh["PhysicalNamesSpan"][0]
    if original[:prefix] != relabeled_bytes[:prefix]:
        raise RelabelError("bytes before $PhysicalNames changed")
    shift = after["PhysicalNamesSpan"][1] - mesh["PhysicalNamesSpan"][1]
    tail_before = np.frombuffer(original[mesh["PhysicalNamesSpan"][1]:], dtype=np.uint8)
    tail_after = np.frombuffer(relabeled_bytes[after["PhysicalNamesSpan"][1]:], dtype=np.uint8)
    if len(tail_before) != len(tail_after):
        raise RelabelError("the byte length after $PhysicalNames changed")
    differing = int(np.count_nonzero(tail_before != tail_after))
    if differing > 8 * len(relabeled):
        raise RelabelError(f"{differing} bytes differ after $PhysicalNames, more than the {8 * len(relabeled)} tag bytes")
    return {"NodesBlockIdentical": True, "ElementsIdenticalApartFromLabels": True, "RelabeledElements": len(relabeled),
            "ChangedTagIntegers": changed_tags, "DifferingBytesAfterPhysicalNames": differing,
            "DifferingBytesBound": 8 * len(relabeled), "PhysicalNamesByteShift": shift}


def parent_area_check(shells, ma_labels, parent_partition_path, tolerance=1e-9):
    """Shell areas sum to the parent MA area (the parent's ownership partition record)."""
    rows = read_csv_rows(parent_partition_path)
    recorded = {int(row["attribute"]): float(row["area"]) for row in rows}
    parents = {}
    for label in ma_labels:
        own = [shell for shell in shells.values() if shell["Parent"] == label]
        total = math.fsum(shell["Area"] for shell in own)
        if label not in recorded:
            raise RelabelError(f"the parent partition record {parent_partition_path} has no area for {label}")
        difference = abs(total - recorded[label]) / recorded[label]
        if difference > tolerance:
            raise RelabelError(f"the shells of {label} sum to {total}, the parent records {recorded[label]} "
                               f"(relative {difference:.3e} > {tolerance:.0e})")
        parents[label] = {"ShellSum": total, "RecordedArea": recorded[label], "RelativeDifference": difference,
                          "Shells": len(own), "Elements": sum(shell["Elements"] for shell in own)}
    return parents


def certificate_parents(data, mesh, relabeled, certificate_path):
    """Every relabeled element's parent equals the parent ownership certificate's
    attribute for its element tag (the tag integer precedes the element's tags)."""
    rows = read_csv_rows(certificate_path)
    by_element = {int(row["element"]): int(row["attribute"]) for row in rows}
    mismatches = 0
    for index, (_, _, parent) in relabeled.items():
        element_tag = struct.unpack_from("<i", data, mesh["Elements"][index][1] - 4)[0]
        if by_element.get(element_tag) != parent:
            mismatches += 1
    if mismatches:
        raise RelabelError(f"{mismatches} relabeled elements disagree with the parent ownership certificate")
    return {"CertificateRows": len(by_element), "RelabeledElementsChecked": len(relabeled), "Mismatches": 0}


def closure_check(shells, tolerance=1e-12):
    """The quadrature measure of every shell (the positive rule's weights) against the
    shells' areas from the independent cross-product formula (element_area)."""
    whole = math.fsum(shell["Area"] for shell in shells.values())
    owned = math.fsum(shell["QuadratureMeasure"] for shell in shells.values())
    closure = abs(owned - whole) / whole
    if closure > tolerance:
        raise RelabelError(f"the radial shell quadrature does not close: {closure}")
    straddling = math.fsum(shell["StraddlingMeasure"] for shell in shells.values())
    return {"WholeMeasure": whole, "OwnedMeasure": owned, "RelativeClosure": closure, "Tolerance": tolerance,
            "StraddlingMeasure": straddling, "StraddlingFraction": straddling / whole,
            "Rule": "degree-4 positive rules (Dunavant 6-point triangles, 3 x 3 Gauss quadrangles); a point straddles when "
                    "its own radial shell differs from its element's centroid shell (ring-aligned tube quadrangles never "
                    "straddle; the graded corner-ball and cap triangles may)"}


def tool_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, text=True,
                                       stderr=subprocess.DEVNULL).strip()
    except (subprocess.CalledProcessError, OSError):
        return None


def load_inputs(case_id, manifest_path, parent_record_path):
    manifest = json.loads(manifest_path.read_text())
    general_mesh_manifest.validate_manifest(manifest, manifest_path)
    case = next((item for item in manifest["Cases"] if item["Id"] == case_id), None)
    if case is None:
        raise RelabelError(f"{case_id} is not a case of {manifest_path}")
    relabel_block = (case.get("Calibration") or {}).get(general_mesh_manifest.RELABEL_KEY)
    if not relabel_block:
        raise RelabelError(f"{case_id} declares no Calibration.{general_mesh_manifest.RELABEL_KEY} block")
    parent_record = json.loads(parent_record_path.read_text())
    parent = next((item for item in parent_record["Cases"] if item["Case"] == relabel_block["BaseCase"]), None)
    if parent is None or not parent.get("Passed"):
        raise RelabelError(f"the parent build record {parent_record_path} has no passed case {relabel_block['BaseCase']}")
    identity = parent["Variants"]["identity"]
    if not Path(identity["Path"]).is_file() or sha256(identity["Path"]) != identity["SHA256"]:
        raise RelabelError(f"the parent identity mesh {identity['Path']} is missing or its SHA256 differs from the record")
    if identity["SHA256"] != relabel_block["ParentMeshSHA256"]:
        raise RelabelError(f"the parent identity mesh SHA256 {identity['SHA256']} differs from the case's declared "
                           f"ParentMeshSHA256 {relabel_block['ParentMeshSHA256']}")
    variant = next(item for item in case["Variants"] if item["Id"] == "identity")
    if [float(v) for v in variant["Transform"]] != [1.0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 1.0]:
        raise RelabelError("the relabel classifies in source-local coordinates: the identity variant must have the "
                           "identity transform")
    parent_root = Path(parent["Root"])
    census = json.loads((parent_root / "build-census.json").read_text())
    radii = [float(v) for v in census["PrismTubes"]["Section"]["RingRadii"]]
    declared = [float(v) for v in relabel_block["RingRadii"]]
    if len(radii) != len(declared) or any(abs(a - b) > 1e-12 for a, b in zip(radii, declared)):
        raise RelabelError(f"the parent build census ring radii {radii} differ from the declared RingRadii {declared}")
    directory = (manifest_path.parent / manifest["RepositoryRoot"]).resolve() / case["Source"]["Directory"]
    files = case["Source"]["Files"]
    for role in ("Signature", "Boundary", "Process"):
        path = directory / files[role]["Name"]
        if sha256(path) != files[role]["SHA256"]:
            raise RelabelError(f"{path} differs from the manifest binding")
    process = tomllib.loads((directory / files["Process"]["Name"]).read_text())
    return {"Manifest": manifest, "ManifestPath": manifest_path, "Case": case, "Relabel": relabel_block,
            "ParentRecord": parent_record, "ParentRecordPath": parent_record_path, "Parent": parent, "ParentRoot": parent_root,
            "Identity": identity, "Radii": radii, "Thickness": float(process["MetalThickness"]),
            "Signature": directory / files["Signature"]["Name"], "Boundary": directory / files["Boundary"]["Name"]}


def run(case_id, *, manifest_path, parent_record_path, root):
    start = time.time()
    inputs = load_inputs(case_id, Path(manifest_path).resolve(), Path(parent_record_path).resolve())
    lines = metal_edge_lines(read_csv_rows(inputs["Boundary"]), read_csv_rows(inputs["Signature"]), inputs["Thickness"])
    original = Path(inputs["Identity"]["Path"]).read_bytes()
    out, mesh, shells, relabeled, ma_labels = relabel(original, lines=lines, radii=inputs["Radii"])
    label_only = assert_label_only(original, out, mesh, relabeled)
    parent_root = inputs["ParentRoot"]
    parents = parent_area_check(shells, ma_labels, parent_root / "identity.msh.interface-partition.csv")
    certificate = certificate_parents(original, mesh, relabeled, parent_root / "identity.msh.interface-partition.csv.elements.csv")
    closure = closure_check(shells)
    case_root = Path(root) / case_id
    case_root.mkdir(parents=True, exist_ok=True)
    mesh_path = case_root / "identity.msh"
    if mesh_path.exists():
        raise RelabelError(f"refusing to overwrite {mesh_path}")
    mesh_path.write_bytes(out)
    mesh_sha = sha256(mesh_path)
    if hashlib.sha256(out).hexdigest() != mesh_sha:
        raise RelabelError("the written mesh differs from the relabeled bytes")
    radii = inputs["Radii"]
    shell_records = []
    for group in sorted(shells):
        shell = shells[group]
        kind, ring, inner, outer = ordinal_description(group[0], radii)
        shell_records.append({"Label": shell["Label"], "Parent": shell["Parent"], "Ordinal": shell["Ordinal"], "Kind": kind,
                              "Ring": ring, "InnerRadius": inner, "OuterRadius": (None if math.isinf(outer) else outer),
                              "ElementaryTag": shell["ElementaryTag"], "Elements": shell["Elements"],
                              "Triangles": shell["Triangles"], "Quadrangles": shell["Quadrangles"], "Area": shell["Area"],
                              "QuadratureMeasure": shell["QuadratureMeasure"], "StraddlingMeasure": shell["StraddlingMeasure"],
                              "StraddlingFraction": shell["StraddlingMeasure"] / shell["Area"],
                              "CentroidDistanceMinimum": shell["CentroidDistanceMinimum"],
                              "CentroidDistanceMaximum": shell["CentroidDistanceMaximum"]})
    census = {"Version": 1, "Kind": inputs["Relabel"]["Kind"], "Case": case_id, "BaseCase": inputs["Relabel"]["BaseCase"],
              "Tool": Path(__file__).name, "ToolSHA256": sha256(__file__), "ToolCommit": tool_commit(),
              "ParentMesh": {"Path": inputs["Identity"]["Path"], "SHA256": inputs["Identity"]["SHA256"]},
              "ParentBuildRecord": {"Path": str(inputs["ParentRecordPath"]), "SHA256": sha256(inputs["ParentRecordPath"])},
              "Mesh": {"Path": str(mesh_path), "SHA256": mesh_sha, "Bytes": len(out)},
              "RingRadii": radii, "TubeRadius": radii[-1], "MetalThickness": inputs["Thickness"], "LengthUnit": "um",
              "LabelRule": LABEL_RULE, "LabelStride": SHELL_LABEL_STRIDE, "FarOrdinal": FAR_ORDINAL,
              "EdgeLines": lines, "Parents": {str(k): v for k, v in parents.items()}, "Shells": shell_records,
              "LabelOnly": label_only, "ParentCertificate": certificate, "RadialQuadrature": closure,
              "Seconds": time.time() - start}
    census_path = Path(str(mesh_path) + CENSUS_SUFFIX)
    census_path.write_text(json.dumps(census, indent=2) + "\n")
    manifest = inputs["Manifest"]
    manifest_path = inputs["ManifestPath"]
    parent = inputs["Parent"]
    calibration = inputs["Case"]["Calibration"]
    record = {"Version": 1, "Command": f"{Path(__file__).name} (label-only relabel of a built case)",
              "Root": str(root),
              "Library": {"CasesAttempted": 1, "CasesBuilt": 0, "CasesRelabeled": 1, "CasesPassed": 1, "CasesUnsupported": 0,
                          "CasesFailed": 0, "CasesUnbuilt": 0, "FlaggedCases": [], "WallClockSeconds": time.time() - start,
                          "Jobs": 1, "Commit": tool_commit(),
                          "Bounds": inputs["ParentRecord"]["Library"]["Bounds"],
                          "Manifest": {"Path": str(manifest_path), "SHA256": sha256(manifest_path), "Kind": "calibration",
                                       "ProductionManifest": str((manifest_path.parent / manifest["Calibration"]["ProductionManifest"]).resolve())}},
              "Cases": [{"Case": case_id, "InventoryStatus": inputs["Case"]["InventoryStatus"],
                         "FixtureVersion": inputs["Case"].get("FixtureVersion"),
                         "Calibration": {"Label": calibration["Label"], "BaseCase": calibration.get("BaseCase"),
                                         "BuildCommandOptions": calibration.get("BuildCommandOptions"),
                                         "ProductionValues": calibration["ProductionValues"],
                                         general_mesh_manifest.RELABEL_KEY: inputs["Relabel"],
                                         "PhysicsRun": calibration.get("PhysicsRun")},
                         "Scope": parent.get("Scope"), "Status": STATUS_RELABELED, "Passed": True, "StoppedBy": None,
                         "CanonicalBuildId": parent.get("CanonicalBuildId"),
                         "Variants": {"identity": {"Path": str(mesh_path), "SHA256": mesh_sha}},
                         "Elements": parent.get("Elements"), "H1": parent.get("H1"), "Estimate": parent.get("Estimate"),
                         "Stages": {}, "Verification": parent.get("Verification"), "HeadroomFlags": [],
                         "Relabel": {"Kind": inputs["Relabel"]["Kind"], "Parent": parent["Case"],
                                     "ParentMesh": census["ParentMesh"], "ParentBuildRecord": census["ParentBuildRecord"],
                                     "ParentRoot": str(parent_root), "Shells": {"Path": str(census_path), "SHA256": sha256(census_path)},
                                     "ShellCount": len(shell_records), "MAParents": ma_labels,
                                     "Rule": "geometry, connectivity, element order and every non-MA label of the parent identity "
                                             "mesh are byte-identical (LabelOnly); the parent's audits (mixed-element quality, "
                                             "protected surfaces, ownership quadrature) hold by identity; the shell partition "
                                             "is audited in the census (Parents, ParentCertificate, RadialQuadrature)"},
                         "Root": str(case_root), "DriverReturnCode": 0, "WallSeconds": time.time() - start}]}
    record_path = Path(root) / "library-build.json"
    record_path.write_text(json.dumps(record, indent=2) + "\n")
    return record, census


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case_id")
    parser.add_argument("--manifest", type=Path, required=True, help="the labeled calibration manifest declaring the case")
    parser.add_argument("--parent-build-record", type=Path, required=True,
                        help="library-build.json holding the passed base case whose identity mesh is relabeled")
    parser.add_argument("--root", type=Path, required=True, help="output root (ROOT/<case>/identity.msh, ROOT/library-build.json)")
    args = parser.parse_args(argv)
    try:
        record, census = run(args.case_id, manifest_path=args.manifest, parent_record_path=args.parent_build_record, root=args.root)
    except RelabelError as error:
        print(f"relabel failed closed: {error}", file=sys.stderr)
        return 1
    print(json.dumps({"Mesh": census["Mesh"], "Shells": len(census["Shells"]), "Parents": census["Parents"],
                      "LabelOnly": census["LabelOnly"], "RadialQuadrature": census["RadialQuadrature"],
                      "BuildRecord": str(Path(args.root) / "library-build.json")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
