#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Write the Palace electrostatic configurations of a validation window from its polygon set
and a mesh manifest (VALIDATION-PLAN.md section (c); stage-1 preparation, decision 207).

Two kinds of run share one polygon set (SCHEMA.md): the conductor labels and the `Terminals`
order of the set fix the Ground / Terminal boundaries and their order on both sides.

  fabricated  the reference on the fabricated mesh of `mesh_polygon_window.jl` (manifest
              `OUT.json`, key `attributes`): the single-transmon reference conventions
              (`single-transmon-finitemetal-anisotropic-20260824/finitemetal-electrostatic-
              anisotropic-r10nm-t5um-p4.json`): Ground = ground_air + ground_substrate,
              Terminal i = <label>_air + <label>_substrate, the 2-nm SA / MS / MA layers as
              surface integrals (SA on substrate_air; MS on the metal-substrate shells; MA on
              the metal-air shells; a separate SA index on substrate_backside ONLY when the
              mesh carries attribute 9, i.e. `Vacuum` > 0 — the truncated window boxes have
              none, review d198 m5), fixed mesh (Refinement.MaxIts 0), BoomerAMG / CG 1e-10,
              MaxIts 5000, Order 4 by default (the pilot ladder: --order 5 for r5 / p5).
  thin        the thin window mesh of `mesh_thin_window.jl` (manifest key `Attributes`) with
              the response-correction library: the accepted transmon verification's settings
              (`combined-verification-20260930/configs/T-comb.json` = protocol path T: p5,
              12 AMR solves, SurfaceMortar 2, FixedTrace, SolveTol 1e-6; `P4-comb.json` = the
              p4 path; the `*-preflight.json` variant = geometry only: p2, MaxIts 0,
              Collocated). Ground = the ground sheets of every plane + the ground bump shells
              (`bump_surface`), Terminal i = the label's sheets + its bump shells (`bump_<label>`;
              a `bump_<label>` naming no terminal is refused), SA = the gap sheets, one MS + one MA target per
              metal plane with the plane's facing normal as EdgeFrameNormal (the flip-chip L2
              faces down; lane W's window preflight convention). A plane without a gap sheet
              (one metal sheet over the whole window, S5 / S6's L2) has no SA on that plane; its
              metal sheets keep MS / MA as PLAIN interfaces (no AutomaticEdges, as the
              fabricated writer's) outside TargetInterfaces: an edgeless sheet owns no metal
              perimeter, which a response target requires, and has nothing to correct (raw =
              ft = sc; decisions 329 / 341: the reference's MS / MA integrate every metal shell,
              so a window MS / MA total = the corrected targets + the raw non-target sheets).
              The library path is a parameter (--library; the geometry-only seed by default for
              the preflight).

    python3 write_window_es_configs.py fabricated --polygon-set S1p.json --manifest S1p_r10nm_t5um.json \\
        --mesh /path/on/the/cluster/S1p_r10nm_t5um.msh2 --postpro /path/postpro --output ref-r10-p4.json [--order 4]
        [--estimator-cheap]   # the fixed-mesh reference never refines: loose estimator tolerance / few iterations
    python3 write_window_es_configs.py thin --polygon-set S1p.json --manifest sct002-S1p-thin.json \\
        --mesh /path/sct002-S1p-thin.msh2 --library /path/process-library.json --path T|P4|preflight \\
        --postpro /path/postpro --output thin-T.json
"""

import argparse
import json
import os

# The fabrication record of the transmon reference (`process-library.json` `Fabrication`):
# 2-nm layers, eps SA 4.0 / MS 11.45 / MA 10.0, substrate 11.45; the loss tangents are the
# recorded configs' (electrostatic participations do not depend on them).
LAYER_THICKNESS_UM = 0.002
LAYERS = {
    "SA": {"Permittivity": 4.0, "LossTan": 0.002},
    "MS": {"Permittivity": 11.45, "LossTan": 0.0003},
    "MA": {"Permittivity": 10.0, "LossTan": 0.03},
}
SUBSTRATE_PERMITTIVITY = 11.45
GROUND = "ground"
# Decision 88(1): the identification's matching radius (the library's MatchingRadius is
# authoritative when the library carries one).
DEFAULT_RADIUS_UM = 1.9
DEFAULT_SEED_LIBRARY = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "surface_response_identification",
    "seeds",
    "preflight_seed_r1p9.json",
)

# Protocol paths of the thin + library run (combined-verification-20260930/configs).
THIN_PATHS = {
    "T": {"Order": 5, "MaxIts": 11, "MaxSize": 400000000, "TraceCoupling": "SurfaceMortar"},
    "P4": {"Order": 4, "MaxIts": 15, "MaxSize": 100000000, "TraceCoupling": "SurfaceMortar"},
    "preflight": {"Order": 2, "MaxIts": 0, "MaxSize": None, "TraceCoupling": "Collocated"},
}


def load_json(path):
    with open(path) as source:
        return json.load(source)


def terminal_labels(polygon_set):
    """The non-ground conductor labels in the attribute order of the set (`Terminals`, else
    sorted), checked against the polygons as the mesher does."""
    labels = set()
    for plane in polygon_set["Planes"]:
        for polygon in plane["Polygons"]:
            labels.add(polygon["Conductor"])
    for bump in polygon_set.get("Bumps", []):
        labels.add(bump.get("Conductor", GROUND))
    labels.discard(GROUND)
    terminals = list(polygon_set.get("Terminals", sorted(labels)))
    if sorted(terminals) != sorted(labels) or len(set(terminals)) != len(terminals):
        raise ValueError(f"Terminals {terminals} must list every non-ground label {sorted(labels)} exactly once")
    return terminals


def dielectric_entry(index, attributes, kind, automatic_edges=None, frame_normal=None):
    entry = {
        "Index": index,
        "Attributes": sorted(attributes),
        "Type": kind,
        "Thickness": LAYER_THICKNESS_UM,
        "Permittivity": LAYERS[kind]["Permittivity"],
        "LossTan": LAYERS[kind]["LossTan"],
    }
    if automatic_edges is not None:
        # The identification's edge patches (the library radius as the last EdgeDistances
        # value), as in every window preflight and the transmon verification configs.
        entry.update({"AutomaticEdges": True, "LocalizeEdgeEnergy": False, "SaveLocalEdgeEnergy": False, "EdgeDistances": [automatic_edges]})
    if frame_normal is not None:
        entry["EdgeFrameNormal"] = [float(x) for x in frame_normal]
    return entry


def backside_against_vacuum(polygon_set):
    """Whether the mesher puts a substrate backside against vacuum (`substrate_backside` faces):
    with two planes when `Vacuum.Below` or `Vacuum.Above` > 0; with ONE plane only on the side
    opposite its facing (an `up` plane's `Vacuum.Above` is the vacuum over the metal, not a
    backside; the OSC windows: `Above` 1000, `Below` 0 -> no backside), as `z_stack` in
    PolygonWindowMesh.jl decides."""
    vacuum = polygon_set.get("Vacuum", {"Below": 0.0, "Above": 0.0})
    below = vacuum.get("Below", 0.0) > 0.0
    above = vacuum.get("Above", 0.0) > 0.0
    planes = polygon_set["Planes"]
    if len(planes) == 1:
        return below if planes[0]["Facing"] == "up" else above
    return below or above


def fabricated_config(polygon_set, manifest, mesh, postpro, order=4, linear_max_its=5000, estimator_cheap=False, verbose=2):
    """The fabricated reference configuration (transmon conventions) from the mesher's
    manifest: `attributes` (name -> attribute / dimension) and `surface_attribute_counts`
    (the attributes that carry faces; `substrate_backside` is listed in the table but has no
    faces when `Vacuum` is 0)."""
    table = {name: int(entry["attribute"]) for name, entry in manifest["attributes"].items()}
    counts = {int(k): int(v) for k, v in manifest.get("surface_attribute_counts", {}).items()}

    def attribute(name):
        if name not in table:
            raise ValueError(f"manifest attribute table has no {name!r}")
        if counts and counts.get(table[name], 0) <= 0:
            raise ValueError(f"attribute {table[name]} ({name}) carries no faces in the mesh")
        return table[name]

    terminals = terminal_labels(polygon_set)
    ground = [attribute("ground_air"), attribute("ground_substrate")]
    terminal_groups = [[attribute(f"{label}_air"), attribute(f"{label}_substrate")] for label in terminals]
    metal_air = [attribute("ground_air")] + [group[0] for group in terminal_groups]
    metal_substrate = [attribute("ground_substrate")] + [group[1] for group in terminal_groups]
    dielectrics = [
        dielectric_entry(1, [attribute("substrate_air")], "SA"),
        dielectric_entry(2, metal_substrate, "MS"),
        dielectric_entry(3, metal_air, "MA"),
    ]
    backside = table.get("substrate_backside")
    if backside is not None and counts.get(backside, 0) > 0:
        # The transmon box carries vacuum beyond the backside: its SA is reported separately
        # (index 4) so the fabrication-response comparison does not include it.
        if not backside_against_vacuum(polygon_set):
            raise ValueError("the mesh carries substrate_backside faces but the polygon set has Vacuum 0")
        dielectrics.append(dielectric_entry(4, [backside], "SA"))
    elif backside_against_vacuum(polygon_set):
        raise ValueError("the polygon set has Vacuum > 0 but the mesh carries no substrate_backside faces")
    linear = {"Type": "BoomerAMG", "KSPType": "CG", "Tol": 1.0e-10, "MaxIts": linear_max_its}
    if estimator_cheap:
        # A fixed-mesh reference (MaxIts 0) uses the error indicator for nothing, yet the
        # recorded references spent 10-21 % of their wall time in an estimator flux solve
        # that never converged (10,000 iterations at 1e-6; 13880958: 881 of 4,193 s). The
        # estimator cannot be switched off; this makes it cheap (plan section (c), m15).
        linear.update({"EstimatorTol": 0.1, "EstimatorMaxIts": 100})
    return {
        "Problem": {"Type": "Electrostatic", "Verbose": verbose, "Output": postpro},
        "Model": {"Mesh": mesh, "L0": 1.0e-6, "CrackInternalBoundaryElements": False, "Refinement": {"MaxIts": 0}},
        "Domains": {
            "Materials": [
                {"Attributes": [table["substrate"]], "Permittivity": SUBSTRATE_PERMITTIVITY},
                {"Attributes": [table["vacuum"]], "Permittivity": 1.0},
            ],
            "Postprocessing": {"Energy": [{"Index": 1, "Attributes": [table["substrate"]]}, {"Index": 2, "Attributes": [table["vacuum"]]}]},
        },
        "Boundaries": {
            "Ground": {"Attributes": sorted(ground)},
            "Terminal": [{"Index": i + 1, "Attributes": sorted(group)} for i, group in enumerate(terminal_groups)],
            "Postprocessing": {"Dielectric": dielectrics},
        },
        "Solver": {
            "Order": order,
            "Device": "CPU",
            "Electrostatic": {"Save": 0},
            "Linear": linear,
        },
    }


def thin_config(polygon_set, manifest, mesh, library, path, postpro, radius=DEFAULT_RADIUS_UM, verbose=1):
    """The thin + library configuration from the thin window manifest (`Attributes`: name ->
    attribute; sheets `ground_<plane>`, `<label>_<plane>`, `gap_<plane>`, bump shells
    `bump_surface` (ground) / `bump_<label>` (terminal), volumes `substrate_<plane lower-case>`,
    `vacuum`; the table lists exactly the groups the mesh carries)."""
    if path not in THIN_PATHS:
        raise ValueError(f"path {path!r} must be one of {sorted(THIN_PATHS)}")
    settings = THIN_PATHS[path]
    table = {name: int(value) for name, value in manifest["Attributes"].items()}
    terminals = terminal_labels(polygon_set)
    plane_names = [plane["Name"] for plane in polygon_set["Planes"]]
    facing = {plane["Name"]: plane["Facing"] for plane in polygon_set["Planes"]}
    substrates = sorted(table[f"substrate_{name.lower()}"] for name in plane_names)
    ground = sorted(table[f"ground_{name}"] for name in plane_names if f"ground_{name}" in table)
    if not ground:
        raise ValueError("no ground sheet on any plane")
    if "bump_surface" in table:
        if "surface" in terminals:
            raise ValueError("terminal label 'surface' is ambiguous with the ground bump group bump_surface")
        ground.append(table["bump_surface"])
    bump_groups = {name[len("bump_"):]: value for name, value in table.items() if name.startswith("bump_") and name != "bump_surface"}
    unknown_bumps = sorted(set(bump_groups) - set(terminals))
    if unknown_bumps:
        raise ValueError(f"bump groups {['bump_' + label for label in unknown_bumps]} name no terminal of {terminals}")
    terminal_groups = []
    for label in terminals:
        group = sorted(table[f"{label}_{name}"] for name in plane_names if f"{label}_{name}" in table)
        if not group:
            raise ValueError(f"terminal {label!r} has no sheet on any plane of the thin mesh")
        if label in bump_groups:
            group = sorted(group + [bump_groups[label]])
        terminal_groups.append(group)
    gaps = sorted(table[f"gap_{name}"] for name in plane_names if f"gap_{name}" in table)
    if not gaps:
        raise ValueError("no gap sheet on any plane (no metal edges)")
    # SA = the gap sheets of the planes that have one; a plane that is one metal sheet over the
    # whole window (S5 / S6's L2) contributes no SA. Its metal sheets still carry MS / MA, as
    # PLAIN interfaces outside TargetInterfaces: a response target needs AutomaticEdges and at
    # least one physical metal-perimeter segment (GetInterfaceMetalEdgeSegmentIndices refuses an
    # edgeless sheet), and an edgeless sheet has nothing to correct, its raw energy is the
    # quantity (decisions 329 / 341; the earlier rule dropped all three interfaces).
    dielectrics = [dielectric_entry(1, gaps, "SA", automatic_edges=radius)]
    target_interfaces = [1]
    for name in plane_names:
        sheets = [table[key] for key in [f"ground_{name}"] + [f"{label}_{name}" for label in terminals] if key in table]
        if not sheets:
            raise ValueError(f"plane {name!r} has no metal sheet in the thin mesh")
        edged = f"gap_{name}" in table
        normal = [0.0, 0.0, 1.0] if facing[name] == "up" else [0.0, 0.0, -1.0]
        for kind in ("MS", "MA"):
            index = len(dielectrics) + 1
            if edged:
                dielectrics.append(dielectric_entry(index, sheets, kind, automatic_edges=radius, frame_normal=normal))
                target_interfaces.append(index)
            else:
                dielectrics.append(dielectric_entry(index, sheets, kind))
    refinement = {"MaxIts": settings["MaxIts"], "UniformLevels": 0, "SerialUniformLevels": 0}
    if settings["MaxIts"] > 0:
        refinement.update({
            "Tol": 0.01, "UpdateFraction": 0.7, "MaximumImbalance": 1.15, "MaxNCLevels": 8, "Nonconformal": True,
            "SaveAdaptIterations": True, "SaveAdaptMesh": False, "MaxSize": settings["MaxSize"],
        })
    return {
        "Problem": {"Type": "Electrostatic", "Verbose": verbose, "Output": postpro, "OutputFormats": {"Paraview": False, "GridFunction": False}},
        "Model": {"Mesh": mesh, "L0": 1.0e-6, "Refinement": refinement},
        "Domains": {
            "Materials": [
                {"Attributes": substrates, "Permittivity": SUBSTRATE_PERMITTIVITY},
                {"Attributes": [table["vacuum"]], "Permittivity": 1.0},
            ],
            "Postprocessing": {"Energy": [{"Index": 1, "Attributes": substrates}, {"Index": 2, "Attributes": [table["vacuum"]]}]},
        },
        "Boundaries": {
            "Ground": {"Attributes": sorted(ground)},
            "Terminal": [{"Index": i + 1, "Attributes": group} for i, group in enumerate(terminal_groups)],
            "Postprocessing": {"Dielectric": dielectrics},
        },
        "Solver": {
            "Order": settings["Order"],
            "Electrostatic": {
                "Save": 0 if path == "preflight" else 1,
                "ResponseCorrection": {
                    "Library": library,
                    "TargetInterfaces": target_interfaces,
                    "UnmatchedPolicy": "Warn",
                    "CorrectionMode": "Both",
                    "TranslationalDomainCorrection": "FixedTrace",
                    "TraceCoupling": settings["TraceCoupling"],
                    "MortarOversampling": 2,
                    "SolveTol": 1.0e-6,
                    "PatchConstruction": "Features",
                },
            },
            "Linear": {"Type": "BoomerAMG", "KSPType": "CG", "Tol": 1.0e-10, "MaxIts": 1000, "EstimatorTol": 0.1, "EstimatorMG": True},
        },
    }


def library_radius(library, default=DEFAULT_RADIUS_UM):
    if library and os.path.exists(library):
        return float(load_json(library).get("MatchingRadius", default))
    return default


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("kind", choices=["fabricated", "thin"])
    parser.add_argument("--polygon-set", required=True, help="the window's polygon set (SCHEMA.md)")
    parser.add_argument("--manifest", required=True, help="the mesh manifest (fabricated OUT.json or the thin *-thin.json)")
    parser.add_argument("--mesh", required=True, help="Model.Mesh as the run will see it (the cluster path)")
    parser.add_argument("--postpro", required=True, help="Problem.Output")
    parser.add_argument("--output", required=True, help="configuration file to write")
    parser.add_argument("--order", type=int, default=4, help="fabricated: Solver.Order (4 = the pilot; 5 = the r5 / p5 ladder step)")
    parser.add_argument("--linear-max-its", type=int, default=5000, help="fabricated: Linear.MaxIts (recorded 5000)")
    parser.add_argument("--estimator-cheap", action="store_true", help="fabricated: EstimatorTol 0.1 / EstimatorMaxIts 100 (plan m15)")
    parser.add_argument("--library", help="thin: the response-correction library (default: the geometry-only seed)")
    parser.add_argument("--path", default="T", choices=sorted(THIN_PATHS), help="thin: protocol path T (p5, 12 solves), P4 or preflight")
    parser.add_argument("--radius", type=float, help="thin: EdgeDistances when the library has no MatchingRadius (default 1.9)")
    args = parser.parse_args(argv)
    polygon_set = load_json(args.polygon_set)
    manifest = load_json(args.manifest)
    if args.kind == "fabricated":
        config = fabricated_config(polygon_set, manifest, args.mesh, args.postpro, order=args.order, linear_max_its=args.linear_max_its, estimator_cheap=args.estimator_cheap)
    else:
        library = args.library or DEFAULT_SEED_LIBRARY
        radius = args.radius if args.radius is not None else library_radius(library)
        config = thin_config(polygon_set, manifest, args.mesh, library, args.path, args.postpro, radius=radius)
    output = os.path.abspath(args.output)
    os.makedirs(os.path.dirname(output), exist_ok=True)
    with open(output, "w") as target:
        json.dump(config, target, indent=2)
        target.write("\n")
    print(output)


if __name__ == "__main__":
    main()
