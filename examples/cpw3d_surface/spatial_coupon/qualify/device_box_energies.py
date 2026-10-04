#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The device-side reading of the (F) DOMAIN twin-consistency check (decision 299, (1a)):
the registration device's thin mesh with the placed coupon boxes sub-tagged in the VOLUME,
so that one thin device solve at two orders on the same (c0) mesh reports the device's box
domain energy bracket [E_in, E_in + E_straddle] and its p-step next to the thin twin's domain
energy under the device trace.

Why the domain only: in the production correction the fixed-trace interface energy is
E_out(R) + sum_patches t^T Q_fab,surf t (the device's own within-R raw surface energy is
dropped, electrostaticsolver.cpp ApplyResponse); the thin twin enters through its DOMAIN
matrix alone (the domain defect Q_fab,dom - Q_thin,dom of the fixed-trace correction and the
self-consistent operator, the fixed-flux transform). The thin twin's surface matrices never
enter ft / ff / sc (surfaceresponseoperator.cpp: the surface defects are built, their only
consumer GetEnergyCorrection has no caller).

Why information only: the device's c0 mesh resolves the sheet edges at micrometres, the thin
twin at its recorded 2-nm cutoff; the domain energy of a sheet-edge field converges like h^1,
so the coarse device reading brackets the twin's but is not a pass / fail quantity. The
criterion of (1a) is the twin's own domain p-stability (spatial_qualification).

Sub-commands:
  relabel  the thin device mesh (Gmsh 2.2, read by meshio) + the library models' SupportPoints
           placed by the device run's patch frames (surface-response-patches.csv, palace.json
           ModelCatalog) -> a copy whose tetrahedra inside / straddling each box carry new
           attributes (one per box, status and original material) + device-boxes.json;
  config   the device run config at the requested orders on the relabelled mesh: the
           materials and the domain energy totals extended by the sub-tags (the totals are
           unchanged), one Domains.Postprocessing.Energy entry per box and status;
  read     the runs' domain-E.csv -> device-box-energies.json (per box and order the bracket,
           the p-step and, when the twin's dense-trace energies are given, the comparison).
"""
import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

STATUS_INSIDE, STATUS_STRADDLE = "Inside", "Straddle"
ATTRIBUTE_BASE = 9000
ENERGY_INDEX_BASE = 1000
ENERGY_INDEX_STATUS = {STATUS_INSIDE: 1, STATUS_STRADDLE: 2}


class DeviceBoxError(ValueError):
    """The device-box inputs are inconsistent (the reason is named)."""


def _rows(path):
    with open(path, newline="") as stream:
        return [{key.strip(): value.strip() for key, value in row.items()} for row in csv.DictReader(stream)]


def placed_boxes(library_path, patches_csv, catalog_json, model_names, length_scale=1.0e-6):
    """[{Model, ModelIndex, Patch, Origin, Axes, SupportPoints, Corners}] of the named models:
    the library's canonical SupportPoints placed by the device run's patch frame (origin and
    axes u, v, w in metres, converted to mesh units by `length_scale`)."""
    library = json.loads(Path(library_path).read_text())
    catalog = json.loads(Path(catalog_json).read_text())["SurfaceResponse"]["ModelCatalog"]
    index_of = {entry["Name"]: int(entry["Index"]) for entry in catalog}
    patches = [{key: float(value) for key, value in row.items()} for row in _rows(patches_csv)]
    boxes = []
    for name in model_names:
        models = [m for m in library["Models"] if m.get("Name") == name]
        if len(models) != 1 or name not in index_of:
            raise DeviceBoxError(f"model {name}: found {len(models)} times in the library, catalog index "
                                 f"{index_of.get(name)}")
        model = models[0]
        if "SupportPoints" not in model:
            raise DeviceBoxError(f"model {name} carries no SupportPoints (not a spatial box model)")
        rows = [row for row in patches if int(row["model"]) == index_of[name]]
        if len(rows) != 1:
            raise DeviceBoxError(f"model {name} (index {index_of[name]}) has {len(rows)} patches in {patches_csv}; "
                                 "the box check needs exactly one placed patch")
        row = rows[0]
        origin = np.array([row["origin x (m)"], row["origin y (m)"], row["origin z (m)"]]) / length_scale
        axes = np.array([[row[f"axis {axis} x"], row[f"axis {axis} y"], row[f"axis {axis} z"]] for axis in "uvw"])
        support = np.asarray(model["SupportPoints"], dtype=float)
        corners = origin + support @ axes
        boxes.append({"Model": name, "ModelIndex": index_of[name], "Patch": int(row.get("patch", -1)) if "patch" in row else None,
                      "Origin": origin.tolist(), "Axes": axes.tolist(), "SupportPoints": support.tolist(),
                      "Corners": corners.tolist()})
    return boxes


def box_membership(points, box, tolerance=1.0e-9):
    """Per point: inside the placed box (the canonical coordinates of the point within the
    SupportPoints' bounds, each side widened by `tolerance` mesh units)."""
    origin, axes = np.asarray(box["Origin"]), np.asarray(box["Axes"])
    local = (np.asarray(points, dtype=float) - origin) @ axes.T
    support = np.asarray(box["SupportPoints"])
    lower, upper = support.min(axis=0) - tolerance, support.max(axis=0) + tolerance
    return np.all((local >= lower) & (local <= upper), axis=1)


def relabel_tetrahedra(points, tetrahedra, materials, boxes):
    """New attributes of the tetrahedra inside / straddling any box; the others keep their
    material. One attribute (ATTRIBUTE_BASE + a running index) per distinct combination of
    the per-box status (none / inside / straddle) and the original material, so that boxes
    overlapping in the volume (the S1p 19-edge pair's margins) are each read exactly from the
    attributes carrying their status. Returns (attributes, records, combinations) with
    records per box: the counts, volumes and the attribute map {status: {material:
    [attributes]}}."""
    points = np.asarray(points, dtype=float)
    tetrahedra = np.asarray(tetrahedra, dtype=int)
    materials = np.asarray(materials, dtype=int)
    corners = points[tetrahedra]
    volumes = np.abs(np.einsum("ij,ij->i", corners[:, 1] - corners[:, 0],
                               np.cross(corners[:, 2] - corners[:, 0], corners[:, 3] - corners[:, 0]))) / 6.0
    status = np.zeros((len(tetrahedra), len(boxes)), dtype=int)
    for column, box in enumerate(boxes):
        per_tet = box_membership(points, box)[tetrahedra]
        inside, some = per_tet.all(axis=1), per_tet.any(axis=1)
        status[inside, column] = 1
        status[some & ~inside, column] = 2
    attributes = materials.copy()
    combinations = []
    tagged = np.flatnonzero(status.any(axis=1))
    keys = sorted({(tuple(int(s) for s in status[t]), int(materials[t])) for t in tagged})
    for key in keys:
        statuses, material = key
        attribute = ATTRIBUTE_BASE + len(combinations) + 1
        selection = np.all(status == np.asarray(statuses), axis=1) & (materials == material)
        attributes[selection] = attribute
        combinations.append({"Attribute": attribute, "Material": material, "BoxStatus": list(statuses),
                             "Tetrahedra": int(selection.sum()), "Volume": float(volumes[selection].sum())})
    records = []
    for column, box in enumerate(boxes):
        mapping = {STATUS_INSIDE: {}, STATUS_STRADDLE: {}}
        counts, volume = {STATUS_INSIDE: 0, STATUS_STRADDLE: 0}, {STATUS_INSIDE: 0.0, STATUS_STRADDLE: 0.0}
        for item in combinations:
            code = item["BoxStatus"][column]
            if code == 0:
                continue
            name = STATUS_INSIDE if code == 1 else STATUS_STRADDLE
            mapping[name].setdefault(str(item["Material"]), []).append(item["Attribute"])
            counts[name] += item["Tetrahedra"]
            volume[name] += item["Volume"]
        records.append({**box, "Index": column + 1, "Attributes": mapping, "Tetrahedra": counts, "Volume": volume})
    return attributes, records, combinations


def command_relabel(args):
    import meshio

    mesh = meshio.read(args.mesh, file_format="gmsh")
    boxes = placed_boxes(args.library, args.patches, args.catalog, args.model)
    physical = mesh.cell_data["gmsh:physical"]
    tet_blocks = [k for k, block in enumerate(mesh.cells) if block.type == "tetra"]
    if len(tet_blocks) != 1:
        raise DeviceBoxError(f"{args.mesh}: expected one tetrahedron block, found {len(tet_blocks)}")
    k = tet_blocks[0]
    attributes, records, combinations = relabel_tetrahedra(mesh.points, mesh.cells[k].data, physical[k], boxes)
    physical[k] = attributes.astype(physical[k].dtype)
    names = {int(value[0]): name for name, value in mesh.field_data.items() if int(value[1]) == 3}
    for item in combinations:
        label = "_".join(f"box{n + 1}{'in' if code == 1 else 'straddle'}" for n, code in enumerate(item["BoxStatus"]) if code)
        mesh.field_data[f"{label}_{names.get(item['Material'], 'material' + str(item['Material']))}"] = \
            np.array([item["Attribute"], 3])
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    mesh_path = out / (Path(args.mesh).stem + "-boxes.msh")
    meshio.write(mesh_path, mesh, file_format="gmsh22", binary=False)
    materials = {str(int(value[0])): name for name, value in mesh.field_data.items() if int(value[1]) == 3}
    record = {"Version": 1, "SourceMesh": str(Path(args.mesh).resolve()), "Mesh": str(mesh_path), "Library": str(args.library),
              "Patches": str(args.patches), "Catalog": str(args.catalog), "Boxes": records, "Combinations": combinations,
              "VolumeNames": materials,
              "Rule": "tetrahedra with every vertex inside a placed SupportPoints box -> Inside, with some -> Straddle; one "
                      "attribute (9000 + a running index) per distinct (per-box status, original material) combination so "
                      "that overlapping boxes are each read exactly; the nodes and every other element are unchanged "
                      "(decision 299 (1a): the device box domain energy bracket [E_in, E_in + E_straddle])"}
    (out / "device-boxes.json").write_text(json.dumps(record, indent=2) + "\n")
    for item in records:
        print(f"box {item['Index']} {item['Model']}: inside {item['Tetrahedra'][STATUS_INSIDE]} tets "
              f"({item['Volume'][STATUS_INSIDE]:.1f}), straddle {item['Tetrahedra'][STATUS_STRADDLE]} "
              f"({item['Volume'][STATUS_STRADDLE]:.1f}) -> {item['Attributes']}")
    print(f"-> {mesh_path}")
    return 0


def energy_index(box_index, status):
    return ENERGY_INDEX_BASE * box_index + ENERGY_INDEX_STATUS[status]


def box_config(reference, boxes_record, mesh, order, output):
    """The device config on the relabelled mesh at `order`: every material and domain energy
    total carrying a sub-tagged attribute's original material also lists the sub-tag (the
    totals are unchanged), plus one energy entry per box and status."""
    config = json.loads(json.dumps(reference))
    config["Model"]["Mesh"] = str(mesh)
    config["Problem"]["Output"] = str(output)
    config["Solver"]["Order"] = int(order)
    extra = {}
    for item in boxes_record["Combinations"]:
        extra.setdefault(int(item["Material"]), []).append(int(item["Attribute"]))
    for entry in config["Domains"]["Materials"]:
        added = [a for material in entry["Attributes"] for a in extra.get(int(material), [])]
        entry["Attributes"] = list(entry["Attributes"]) + added
    postprocessing = config["Domains"].setdefault("Postprocessing", {})
    energies = postprocessing.setdefault("Energy", [])
    for entry in energies:
        added = [a for material in entry["Attributes"] for a in extra.get(int(material), [])]
        entry["Attributes"] = list(entry["Attributes"]) + added
    for box in boxes_record["Boxes"]:
        for status, mapping in box["Attributes"].items():
            attributes = sorted(a for values in mapping.values() for a in values)
            if attributes:
                energies.append({"Index": energy_index(box["Index"], status), "Attributes": attributes})
    return config


def command_config(args):
    reference = json.loads(Path(args.config).read_text())
    boxes_record = json.loads(Path(args.boxes).read_text())
    out = Path(args.output)
    written = {}
    for order in args.orders:
        directory = out / f"p{order}"
        directory.mkdir(parents=True, exist_ok=True)
        config = box_config(reference, boxes_record, args.mesh or boxes_record["Mesh"], order, directory / "postpro")
        (directory / "config.json").write_text(json.dumps(config, indent=2) + "\n")
        written[f"p{order}"] = str(directory / "config.json")
    (out / "device-box-configs.json").write_text(json.dumps({"Version": 1, "Boxes": str(args.boxes), "Configs": written,
                                                             "Output": {f"p{order}": str(out / f"p{order}" / "postpro") for order in args.orders}},
                                                            indent=2) + "\n")
    print(f"{len(written)} configs -> {out}")
    return 0


def domain_energies(postpro, source=1):
    """{energy index: E_elec[index] (J)} of one source from domain-E.csv (+ the total under 0)."""
    for row in _rows(Path(postpro) / "domain-E.csv"):
        if int(float(row["i"])) == source:
            result = {0: float(row["E_elec (J)"])}
            for key, value in row.items():
                if key.startswith("E_elec[") and key.endswith("]"):
                    result[int(key[len("E_elec["):-1])] = float(value)
            return result
    raise DeviceBoxError(f"{postpro}/domain-E.csv has no source {source}")


def box_energy_record(boxes_record, runs, twin=None):
    """Per box and order: Inside, Straddle, the bracket [In, In + Straddle] and its central
    value; the p-step between consecutive orders; when `twin` gives the thin twin's domain
    energy under the device trace per order ({"p4": E, "p5": E}), the twin's step and the
    ratios twin / device (central and bracket)."""
    orders = sorted(runs, key=lambda name: int(name[1:]))
    record = {"Version": 1, "Orders": orders, "Boxes": [],
              "Rule": "device = the registration device's thin solve at each order on ONE c0 mesh with the placed box "
                      "sub-tagged in the volume: the box domain energy bracket [Inside, Inside + Straddle] (central = Inside + "
                      "Straddle / 2) and its step between consecutive orders; twin = the thin coupon's domain energy under "
                      "the device trace (E_elec of the dense thin runs x EnergyScale). INFORMATION ONLY (decision 299 (1a)): the "
                      "device c0 mesh resolves the sheet edges at micrometres, the twin at its 2-nm cutoff; the domain energy "
                      "of a sheet-edge field converges like h^1, so the coarse device reading is not a pass / fail quantity"}
    for box in boxes_record["Boxes"]:
        entry = {"Model": box["Model"], "Index": box["Index"], "Device": {}, "Twin": None}
        for name in orders:
            energies = domain_energies(runs[name])
            inside = energies.get(energy_index(box["Index"], STATUS_INSIDE), 0.0)
            straddle = energies.get(energy_index(box["Index"], STATUS_STRADDLE), 0.0)
            entry["Device"][name] = {"Inside": inside, "Straddle": straddle, "Central": inside + 0.5 * straddle,
                                     "Bracket": [inside, inside + straddle], "Total": energies[0]}
        steps = {}
        for low, high in zip(orders, orders[1:]):
            a, b = entry["Device"][low], entry["Device"][high]
            steps[f"{low}->{high}"] = {"Central": b["Central"] - a["Central"],
                                       "RelativeCentral": (b["Central"] - a["Central"]) / b["Central"] if b["Central"] else None,
                                       "Inside": b["Inside"] - a["Inside"]}
        entry["DeviceStep"] = steps
        if twin and box["Model"] in twin:
            values = twin[box["Model"]]
            entry["Twin"] = {"Energy": values,
                             "Step": {f"{low}->{high}": values[high] - values[low] for low, high in zip(orders, orders[1:])
                                      if low in values and high in values},
                             "TwinOverDevice": {name: {"Central": values[name] / entry["Device"][name]["Central"]
                                                       if entry["Device"][name]["Central"] else None,
                                                       "Bracket": sorted([values[name] / b if b else float("inf")
                                                                          for b in entry["Device"][name]["Bracket"]])}
                                                for name in orders if name in values}}
        record["Boxes"].append(entry)
    return record


def command_read(args):
    boxes_record = json.loads(Path(args.boxes).read_text())
    runs = {}
    for item in args.run:
        if "=" not in item:
            raise DeviceBoxError("--run needs pORDER=POSTPRO")
        name, path = item.split("=", 1)
        runs[name] = path
    twin = json.loads(Path(args.twin).read_text()) if args.twin else None
    record = box_energy_record(boxes_record, runs, twin)
    Path(args.output).write_text(json.dumps(record, indent=2) + "\n")
    for box in record["Boxes"]:
        print(box["Model"], {name: f"{v['Central']:.4e} [{v['Bracket'][0]:.4e}, {v['Bracket'][1]:.4e}]"
                             for name, v in box["Device"].items()}, box["DeviceStep"])
    return 0


def add_arguments(parser):
    commands = parser.add_subparsers(dest="device_box_command", required=True)
    relabel = commands.add_parser("relabel", help="sub-tag the device mesh's tetrahedra by the placed boxes")
    relabel.add_argument("--mesh", required=True, help="the registration device's thin mesh (Gmsh 2.2)")
    relabel.add_argument("--library", required=True, help="the process library placing the models (SupportPoints)")
    relabel.add_argument("--patches", required=True, help="surface-response-patches.csv of the device run on that library")
    relabel.add_argument("--catalog", required=True, help="palace.json of that run (SurfaceResponse.ModelCatalog)")
    relabel.add_argument("--model", action="append", required=True, help="a spatial model name (repeatable)")
    relabel.add_argument("--output", required=True)
    relabel.set_defaults(func=command_relabel)
    config = commands.add_parser("config", help="the device config at the requested orders on the relabelled mesh")
    config.add_argument("--config", required=True, help="the device run config (the c0 thin solve with the library)")
    config.add_argument("--boxes", required=True, help="device-boxes.json written by relabel")
    config.add_argument("--orders", type=int, nargs="+", default=[4, 5])
    config.add_argument("--mesh", help="the relabelled mesh's path as the solver sees it (default: the relabel record's)")
    config.add_argument("--output", required=True)
    config.set_defaults(func=command_config)
    read = commands.add_parser("read", help="the device box domain energies of the runs")
    read.add_argument("--boxes", required=True)
    read.add_argument("--run", action="append", required=True, metavar="pORDER=POSTPRO")
    read.add_argument("--twin", help="JSON {model: {pORDER: twin domain energy under the device trace}}")
    read.add_argument("--output", required=True)
    read.set_defaults(func=command_read)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(parser)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
