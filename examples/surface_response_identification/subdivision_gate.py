#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Collinear-subdivision invariance gate (supervisor decisions 212 / 213; VALIDATION-PLAN
(h)-8). Inserting collinear vertices on a chord is a geometric no-op, so the identification
must read a mesh and the same mesh with every edge subdivided identically — whatever the arc /
corner rule. Every synthetic layout mesh of the existing gates (the stress / decision-82 /
fillet suites, the arc-cluster set, the stack suite) and the joint-only layouts of
``synthetic_layouts.subdivision_suite`` are identified as meshed and as

  * ``refine-0.5``: every mesh edge split at its midpoint (refine_msh2, 1 -> 8 tetrahedra);
  * ``refine-0.37``: every mesh edge split at 0.37 of its length (a second, uneven spacing);
  * ``rotate-37``: the mesh rotated by 37 deg about the plan-view normal with a seeded
    renumbering (permute_msh2): the canonical, coordinate-sorted numbering gives every closed
    loop another start vertex (the closed-loop scan's seed), without any subdivision;
  * ``rotate-37-refine-0.5``: both.

The gate passes for a layout when every variant has the same GeometryDigest, the same
multiset of features (Type, Signature, Chirality, Length to 1e-5), the same exclusion classes
(class, count, length) and the same vertex types (type, turn) as the base mesh.

    python3 -m surface_response_identification.subdivision_gate --output DIR \\
        --meshes DIR [DIR ...] [--subdivision-suite-meshes DIR] [--palace P] [--np 1] \\
        [--fractions 0.5 0.37] [--rotate-degrees 37] [--only substring] [--jobs 2]
"""

import argparse
import collections
import concurrent.futures
import json
import os
import sys
import time

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "surface_response_identification"

from . import synthetic_layouts as S  # noqa: E402
from .msh2 import read_msh2  # noqa: E402
from .permute_msh2 import permute_file  # noqa: E402
from .preflight_config import preflight_config  # noqa: E402
from .refine_msh2 import refine_file  # noqa: E402

ROTATION_SEED = 11


def layout_index():
    """Layout name -> layout for every synthetic suite (ports and frame normals of the
    preflight configuration)."""
    index = {}
    for suite in (S.stress_suite, S.decision_82_suite, S.fillet_suite, S.arc_cluster_suite, S.stack_suite, S.subdivision_suite):
        for lay in suite():
            index[lay["Name"]] = lay
    return index


def config_arguments(mesh_path, lay):
    mesh = read_msh2(mesh_path)
    attributes = sorted({int(t) for t in mesh.physical_tags(2)})
    ports = [p["Attribute"] for p in (lay or {}).get("Ports", [])]
    return dict(
        ground=[S.GROUND] if S.GROUND in attributes else [],
        terminals=[[S.SECOND_CONDUCTOR]] if S.SECOND_CONDUCTOR in attributes else [],
        sa=[S.SUBSTRATE_AIR] if S.SUBSTRATE_AIR in attributes else [],
        library=S.DEFAULT_LIBRARY,
        ports=ports,
        frame_normal=(lay or {}).get("FrameNormal"),
    )


def content(manifest_path):
    """Coordinate-free content of a manifest: digest, features, exclusions, vertices."""
    with open(manifest_path) as source:
        ident = json.load(source)["Identification"]
    features = sorted((f["Type"], json.dumps(f["Signature"], sort_keys=True), f["Chirality"], round(f["Length"], 5)) for f in ident["Features"])
    exclusions = sorted((e["Class"], e["Count"], round(e["Length"], 5)) for e in ident["Exclusions"])
    vertices = sorted((v["Type"], round(v.get("TurnDegrees", 0.0), 5)) for v in ident["Vertices"])
    arcs = sorted((a.get("Kind"), round(a.get("RadiusOverR", 0.0), 6), a.get("Joints")) for a in ident.get("Arcs", []))
    return {"Digest": ident["GeometryDigest"], "Features": features, "Exclusions": exclusions, "Vertices": vertices, "Arcs": arcs}


def identify(palace, mesh_path, arguments, directory, ranks):
    os.makedirs(directory, exist_ok=True)
    manifest = os.path.join(directory, "postpro", "surface-response-requirements.json")
    if not os.path.exists(manifest):
        config = preflight_config(mesh_path, **arguments, output_directory=os.path.join(directory, "postpro"))
        config_path = os.path.join(directory, "config.json")
        with open(config_path, "w") as target:
            json.dump(config, target, indent=1)
        S.run_preflight(palace, config_path, ranks, os.path.join(directory, "palace.log"))
    return manifest if os.path.exists(manifest) else None


def variants(fractions, rotate_degrees):
    out = [("base", 0.0, None)]
    out += [(f"refine-{f:g}", 0.0, f) for f in fractions]
    if rotate_degrees:
        out.append((f"rotate-{rotate_degrees:g}", rotate_degrees, None))
        if fractions:
            out.append((f"rotate-{rotate_degrees:g}-refine-{fractions[0]:g}", rotate_degrees, fractions[0]))
    return out


def run_layout(name, mesh_path, lay, args):
    started = time.time()
    directory = os.path.join(args.output, name)
    os.makedirs(directory, exist_ok=True)
    arguments = config_arguments(mesh_path, lay)
    row = {"Layout": name, "Mesh": mesh_path, "Variants": {}, "Pass": True}
    base = None
    for label, rotation, fraction in variants(args.fractions, args.rotate_degrees):
        variant_mesh = mesh_path
        variant_directory = os.path.join(directory, label)
        os.makedirs(variant_directory, exist_ok=True)
        if rotation:
            rotated = os.path.join(variant_directory, "rotated.msh2")
            if not os.path.exists(rotated):
                permute_file(mesh_path, rotated, seed=ROTATION_SEED, rotate_degrees=rotation)
            variant_mesh = rotated
        if fraction is not None:
            refined = os.path.join(variant_directory, "refined.msh2")
            if not os.path.exists(refined):
                refine_file(variant_mesh, refined, levels=1, fraction=fraction)
            variant_mesh = refined
        manifest = identify(args.palace, variant_mesh, arguments, variant_directory, args.np)
        if manifest is None:
            row["Variants"][label] = {"Error": "no manifest"}
            row["Pass"] = False
            continue
        this = content(manifest)
        entry = {"Digest": this["Digest"], "FeatureCount": len(this["Features"]), "Classes": dict(collections.Counter(t for t, _, _, _ in this["Features"])), "Arcs": len(this["Arcs"])}
        if label == "base":
            base = this
        else:
            same = {k: this[k] == base[k] for k in ("Digest", "Features", "Exclusions", "Vertices", "Arcs")} if base else {}
            entry["Same"] = same
            entry["Pass"] = bool(base) and all(same.values())
            if not entry["Pass"] and base:
                only_base = sorted(set(base["Features"]) - set(this["Features"]))
                only_this = sorted(set(this["Features"]) - set(base["Features"]))
                entry["OnlyBase"] = [f"{t} chir {c} L {l} {s[:90]}" for t, s, c, l in only_base][:6]
                entry["OnlyVariant"] = [f"{t} chir {c} L {l} {s[:90]}" for t, s, c, l in only_this][:6]
                entry["ArcsBase"], entry["ArcsVariant"] = base["Arcs"], this["Arcs"]
            row["Pass"] = row["Pass"] and entry["Pass"]
        row["Variants"][label] = entry
        for name_ in ("rotated.msh2", "refined.msh2"):
            path = os.path.join(variant_directory, name_)
            if os.path.exists(path) and not args.keep_meshes:
                os.remove(path)
    row["Seconds"] = round(time.time() - started, 1)
    return row


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", required=True)
    parser.add_argument("--meshes", nargs="*", default=[], help="directories of synthetic gate meshes (*.msh2)")
    parser.add_argument("--subdivision-suite-meshes", default=None, help="directory of the joint-only subdivision_suite meshes (generated with Gmsh when missing)")
    parser.add_argument("--palace", default=S.DEFAULT_PALACE)
    parser.add_argument("--np", type=int, default=1)
    parser.add_argument("--fractions", type=float, nargs="*", default=[0.5, 0.37])
    parser.add_argument("--rotate-degrees", type=float, default=37.0)
    parser.add_argument("--only", default=None)
    parser.add_argument("--jobs", type=int, default=1, help="layouts identified concurrently (each at --np ranks)")
    parser.add_argument("--keep-meshes", action="store_true")
    parser.add_argument("--julia", default="julia")
    args = parser.parse_args(argv)
    os.makedirs(args.output, exist_ok=True)
    index = layout_index()
    meshes = []
    if args.subdivision_suite_meshes:
        suite = S.subdivision_suite()
        missing = [lay for lay in suite if not os.path.exists(os.path.join(args.subdivision_suite_meshes, lay["Name"] + ".msh2"))]
        if missing:
            rc, seconds, log = S.generate_meshes(missing, args.subdivision_suite_meshes, args.julia)
            print(json.dumps({"GeneratedSubdivisionSuite": [lay["Name"] for lay in missing], "ReturnCode": rc, "Seconds": round(seconds, 1), "Log": log}), flush=True)
        meshes += [(lay["Name"], os.path.join(args.subdivision_suite_meshes, lay["Name"] + ".msh2")) for lay in suite]
    for directory in args.meshes:
        for file in sorted(os.listdir(directory)):
            if file.endswith(".msh2"):
                meshes.append((file[:-5], os.path.join(directory, file)))
    meshes = [(n, p) for n, p in meshes if os.path.exists(p) and (args.only is None or args.only in n)]
    rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = [pool.submit(run_layout, name, path, index.get(name), args) for name, path in meshes]
        for future in futures:
            row = future.result()
            rows.append(row)
            print(json.dumps({k: row[k] for k in ("Layout", "Pass", "Seconds")} | {"Variants": {k: v.get("Pass", v.get("Digest")) for k, v in row["Variants"].items()}}), flush=True)
            if not row["Pass"]:
                for label, entry in row["Variants"].items():
                    if label != "base" and not entry.get("Pass"):
                        print("   FAIL", label, json.dumps({k: entry.get(k) for k in ("Same", "OnlyBase", "OnlyVariant", "ArcsBase", "ArcsVariant", "Error")}, default=str)[:600], flush=True)
    passed = sum(1 for r in rows if r["Pass"])
    labels = [v[0] for v in variants(args.fractions, args.rotate_degrees)][1:]
    print(f"SUBDIVISION GATE {passed} / {len(rows)} layouts PASS (variants {', '.join(labels)}: identical GeometryDigest, features, exclusions, vertices, arcs)")
    with open(os.path.join(args.output, "results.json"), "w") as target:
        json.dump({"Variants": labels, "Rows": rows, "Pass": passed, "Total": len(rows)}, target, indent=1, default=str)
    return 0 if passed == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
