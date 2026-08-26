#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Discover a library-independent surface-response requirement closure."""

import argparse
import copy
import hashlib
import json
import shlex
import subprocess
from pathlib import Path

import prepare_surface_response_coupons as planner


def load_json(path):
    return json.loads(path.read_text())


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def response_section(config):
    solver = config.setdefault("Solver", {})
    if config.get("Problem", {}).get("Type") == "Electrostatic":
        electrostatic = solver.setdefault("Electrostatic", {})
        response = electrostatic.get("ResponseCorrection")
        if response is None:
            raise ValueError(
                "Electrostatic closure discovery requires "
                "Solver.Electrostatic.ResponseCorrection"
            )
        return response
    response = solver.get("SurfaceResponseCorrection")
    if response is None:
        raise ValueError(
            "Maxwell closure discovery requires Solver.SurfaceResponseCorrection"
        )
    return response


def resolve_library(config_path, response):
    path = Path(response["Library"]).expanduser()
    return path.resolve() if path.is_absolute() else (config_path.parent / path).resolve()


def unique_interfaces(requirement):
    result = []
    seen = set()
    for entry in requirement.get("Interfaces", []):
        key = (int(entry.get("Slot", 0)), entry["Type"])
        if key in seen:
            continue
        seen.add(key)
        result.append(
            {"Slot": key[0], "Type": key[1], "Coupon": len(result) + 1}
        )
    result.sort(key=lambda item: (item["Slot"], item["Type"]))
    for index, entry in enumerate(result, start=1):
        entry["Coupon"] = index
    return result


def conductor_references(edges, separation=1.0):
    count = max((int(edge.get("Conductor", 1)) for edge in edges), default=1)
    references = []
    for conductor in range(1, count + 1):
        points = [
            edge.get("Point")
            for edge in edges
            if int(edge.get("Conductor", 1)) == conductor and "Point" in edge
        ]
        if points:
            reference = [float(value) for value in points[0]]
        else:
            reference = [float(conductor - 1) * separation, 0.0, 0.0]
        # Preflight requires distinct references but does not read response matrices.
        if reference in references:
            reference[0] += conductor * max(separation, 1.0)
        references.append(reference)
    return references


def placeholder_model(requirement, matching_radius):
    signature = planner.coupon_signature(requirement)
    digest = hashlib.sha256(
        json.dumps(signature, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()[:16]
    topology = requirement["Topology"]
    geometry = requirement.get("Geometry", {})
    boundary = requirement.get("BoundaryCondition", {"Type": "PEC"})
    model = {
        "Name": f"__preflight_placeholder_{digest}",
        "Topology": topology,
        "FabricatedMatrix": "__preflight_dummy_fabricated.csv",
        "ThinMatrix": "__preflight_dummy_thin.csv",
        "FabricatedSurfaceMatrix": "__preflight_dummy_fabricated_surface.csv",
        "ThinSurfaceMatrix": "__preflight_dummy_thin_surface.csv",
        "BasisPoints": "__preflight_dummy_basis.csv",
        "Interfaces": unique_interfaces(requirement),
        "BoundaryLawQualification": {
            "Version": 1,
            "Status": "Unqualified",
            "Calibration": "GeometryDiscoveryOnly",
            "FrequencyUniversal": False,
        },
    }

    if topology == "IsolatedEdge":
        model["Reference"] = [0.0, 0.0, 0.0]
        model["BoundaryCondition"] = boundary
    elif topology in (
        "SameConductorGap",
        "SameConductorStrip",
        "DifferentConductorGap",
    ):
        separation = float(geometry["Separation"])
        model["Separation"] = separation
        model["SeparationTolerance"] = max(1.0e-9 * matching_radius, 1.0e-12)
        model["BoundaryCondition"] = boundary
        if topology == "DifferentConductorGap":
            model["ConductorReferences"] = [
                [-0.5 * separation, 0.0, 0.0],
                [0.5 * separation, 0.0, 0.0],
            ]
        else:
            model["Reference"] = [0.0, 0.0, 0.0]
    elif topology == "ParallelEdgeCluster":
        edges = copy.deepcopy(geometry["Edges"])
        for edge in edges:
            if isinstance(edge.get("Offset"), list):
                edge["Offset"] = float(edge["Offset"][0])
            if isinstance(edge.get("GapDirection"), list):
                edge["GapDirection"] = int(edge["GapDirection"][0])
        model["Edges"] = edges
        model["EdgeOffsetTolerance"] = max(1.0e-9 * matching_radius, 1.0e-12)
        model["BoundaryCondition"] = boundary
        references = conductor_references(edges, matching_radius)
        if len(references) == 1:
            model["Reference"] = references[0]
        else:
            model["ConductorReferences"] = references
    elif topology in ("ConvexCorner", "ConcaveCorner"):
        model["Angle"] = float(geometry["AngleDegrees"])
        model["AngleTolerance"] = 1.0e-6
        model["CornerRadius"] = float(geometry.get("CornerRadius", 0.0))
        model["CornerRadiusTolerance"] = max(
            1.0e-9 * matching_radius, 1.0e-12
        )
        model["Reference"] = [0.0, 0.0, 0.0]
        model["BoundaryCondition"] = boundary
    elif topology == "SpatialEdgeCluster":
        edges = copy.deepcopy(geometry["Edges"])
        model["Edges"] = edges
        model["EdgePositionTolerance"] = max(
            1.0e-4 * matching_radius, 1.0e-12
        )
        model["EdgeAngleTolerance"] = 1.0e-3
        references = conductor_references(edges, matching_radius)
        if len(references) == 1:
            model["Reference"] = references[0]
        else:
            model["ConductorReferences"] = references
        # Multi-conductor signatures need the exact mask to disambiguate ownership.
        # For a single conductor, edge geometry is sufficient for virtual closure and
        # avoids round-trip sensitivity in canonicalized sampled face loops. The full
        # PlanViewFacets remain in the final requirement fingerprint.
        if len(references) > 1:
            for key in ("PlanViewBoundary", "MaskRegularization"):
                if key in geometry:
                    model[key] = copy.deepcopy(geometry[key])
    elif topology in ("Endpoint", "Junction"):
        model["Reference"] = [0.0, 0.0, 0.0]
        model["BoundaryCondition"] = boundary
        if topology == "Junction":
            model["ArmAngles"] = geometry["ArmAnglesDegrees"]
            model["ArmAngleTolerance"] = 1.0e-6
        for key in ("PlanViewBoundary", "MaskRegularization"):
            if key in geometry:
                model[key] = copy.deepcopy(geometry[key])
    else:
        raise ValueError(f"Unsupported placeholder topology {topology}")
    return digest, model


def run_preflight(palace, config_path, log_path):
    command = [str(palace), "--surface-response-preflight", str(config_path)]
    print("+ " + shlex.join(command), flush=True)
    with log_path.open("w") as stream:
        subprocess.run(command, check=True, stdout=stream, stderr=subprocess.STDOUT)


def restore_source_status(manifest, source_library, placeholder_requirements):
    result = copy.deepcopy(manifest)
    result["Library"]["Name"] = source_library.get(
        "Name", Path(result["Library"]["Path"]).stem
    )
    result["Library"]["Path"] = str(source_library["__SourcePath"])
    counts = {"Exact": 0, "Interpolated": 0, "Missing": 0}
    lengths = {"Exact": 0.0, "Interpolated": 0.0, "Missing": 0.0}
    for requirement in result["Requirements"]:
        selected = requirement.get("SelectedModels", [])
        placeholders = [
            model.get("Name")
            for model in selected
            if model.get("Name") in placeholder_requirements
        ]
        if placeholders:
            # Use the requirement which created the selected virtual model. The production
            # matcher may describe the same match in model-local coordinates on later
            # passes; feeding that transformed description back to the coupon generator
            # would create a different signature than the virtual model that established
            # closure.
            source = placeholder_requirements[placeholders[0]]
            for key in ("Topology", "Geometry", "BoundaryCondition", "Interfaces"):
                if key in source:
                    requirement[key] = copy.deepcopy(source[key])
                else:
                    requirement.pop(key, None)
            requirement["Status"] = "Missing"
            requirement["Reason"] = (
                "Missing from source library after exhaustive geometry discovery"
            )
            requirement.pop("SelectedModels", None)
            requirement.pop("NormalizedLibraryDistance", None)
        status = requirement["Status"]
        counts[status] += int(requirement["Count"])
        lengths[status] += float(requirement.get("TotalEdgeLength", 0.0))
    result["Summary"] = {"Counts": counts, "TotalEdgeLengths": lengths}
    result["Complete"] = counts["Missing"] == 0
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--palace", type=Path, required=True)
    parser.add_argument("--max-passes", type=int, default=8)
    args = parser.parse_args()
    if args.max_passes <= 0:
        parser.error("--max-passes must be positive")

    config_path = args.config.expanduser().resolve()
    config = load_json(config_path)
    response = response_section(config)
    source_path = resolve_library(config_path, response)
    source_library = load_json(source_path)
    source_library["__SourcePath"] = source_path
    virtual_library = copy.deepcopy(source_library)
    virtual_library.pop("__SourcePath", None)
    # Version 2 is sufficient for virtual multi-conductor references. Preserve version 3
    # when fabrication metadata is present, but do not force older libraries to invent it.
    virtual_library["Version"] = max(2, int(virtual_library.get("Version", 0)))
    virtual_library["Name"] = f"{virtual_library.get('Name', source_path.stem)}-closure"

    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    known = set()
    placeholder_names = set()
    placeholder_requirements = {}
    history = []
    final_manifest = None
    for pass_index in range(1, args.max_passes + 1):
        pass_root = output / f"pass-{pass_index:02d}"
        pass_root.mkdir(parents=True, exist_ok=True)
        library_path = pass_root / "process-library.json"
        write_json(library_path, virtual_library)
        pass_config = copy.deepcopy(config)
        pass_config["Problem"]["Output"] = str(pass_root / "postpro")
        pass_response = response_section(pass_config)
        pass_response["Library"] = str(library_path)
        pass_response["UnmatchedPolicy"] = "Warn"
        pass_config_path = pass_root / "config.json"
        write_json(pass_config_path, pass_config)
        run_preflight(args.palace, pass_config_path, pass_root / "preflight.log")
        manifest = load_json(pass_root / "postpro" / "surface-response-requirements.json")
        final_manifest = manifest

        added = []
        for requirement in manifest["Requirements"]:
            if requirement["Status"] != "Missing":
                continue
            digest, model = placeholder_model(
                requirement, float(manifest["Library"]["MatchingRadius"])
            )
            if digest in known:
                continue
            known.add(digest)
            placeholder_names.add(model["Name"])
            placeholder_requirements[model["Name"]] = copy.deepcopy(requirement)
            virtual_library["Models"].append(model)
            added.append({"Id": digest, "Topology": requirement["Topology"]})
        history.append(
            {
                "Pass": pass_index,
                "Summary": manifest["Summary"],
                "AddedPlaceholders": added,
            }
        )
        if not added:
            if manifest["Summary"]["Counts"]["Missing"]:
                raise RuntimeError(
                    "Geometry closure stalled with missing requirements; see "
                    f"{pass_root / 'postpro/surface-response-requirements.json'}"
                )
            break
    else:
        raise RuntimeError(f"Geometry closure did not converge in {args.max_passes} passes")

    source_library["__SourcePath"] = source_path
    final_manifest = restore_source_status(
        final_manifest, source_library, placeholder_requirements
    )
    manifest_path = output / "surface-response-requirements.json"
    write_json(manifest_path, final_manifest)
    write_json(
        output / "closure-history.json",
        {
            "Version": 1,
            "SourceConfig": str(config_path),
            "SourceLibrary": str(source_path),
            "Passes": history,
            "PlaceholderCount": len(placeholder_names),
            "CompleteAgainstSourceLibrary": final_manifest["Complete"],
        },
    )
    print(manifest_path)
    print(json.dumps(final_manifest["Summary"], indent=2))


if __name__ == "__main__":
    main()
