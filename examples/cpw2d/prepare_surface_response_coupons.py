#!/usr/bin/env python3

"""Convert a Palace surface-response preflight manifest into a coupon work plan."""

import argparse
import json
import math
from pathlib import Path


PAIRED_TOPOLOGIES = {
    "SameConductorGap",
    "DifferentConductorGap",
    "SameConductorStrip",
}
CORNER_TOPOLOGIES = {"ConvexCorner", "ConcaveCorner"}
SPATIAL_TOPOLOGIES = {
    "ParallelEdgeCluster",
    "SpatialEdgeCluster",
    "Endpoint",
    "Junction",
}


def load_json(path):
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError as error:
        raise ValueError(f"Invalid JSON in {path}: {error}") from error


def slug(value):
    return (
        f"{value:.12g}".replace("-", "m").replace(".", "p")
        if isinstance(value, float)
        else str(value).lower()
    )


def coupon_id(requirement):
    topology = requirement["Topology"]
    geometry = requirement.get("Geometry", {})
    parts = [topology.lower()]
    for key in ("Separation", "AngleDegrees", "CornerRadius", "ArmCount", "EdgeCount"):
        if key in geometry:
            parts.append(f"{key.lower()}-{slug(geometry[key])}")
    return "_".join(parts)


def preparation(requirement):
    topology = requirement["Topology"]
    geometry = requirement.get("Geometry", {})
    if topology == "IsolatedEdge" or topology in PAIRED_TOPOLOGIES:
        return {
            "Method": "StraightEdgeBuilder",
            "Driver": "examples/cpw2d/build_surface_process_library.py",
        }
    if topology in CORNER_TOPOLOGIES:
        angle = float(geometry.get("AngleDegrees", math.nan))
        if math.isclose(angle, 90.0, abs_tol=1.0e-8):
            return {
                "Method": "CornerPrototype",
                "MeshGenerator":
                    "examples/cpw3d_surface/corner_coupon/mesh_corner_coupon.jl",
                "ResponseGenerator":
                    "examples/cpw3d_surface/corner_coupon/generate_corner_response.py",
            }
        return {
            "Method": "Manual",
            "Reason": "The current corner coupon generator supports 90-degree corners",
        }
    if topology in SPATIAL_TOPOLOGIES:
        return {
            "Method": "Manual",
            "Reason": (
                "A canonical 3D mesh generator is not yet available for this spatial "
                "edge signature"
            ),
        }
    return {"Method": "Manual", "Reason": "Unknown coupon topology"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--include-matched",
        action="store_true",
        help="Include exact and interpolated requirements as validation coupons",
    )
    args = parser.parse_args()

    manifest_path = args.manifest.expanduser().resolve()
    manifest = load_json(manifest_path)
    if manifest.get("Version") != 1 or not isinstance(
        manifest.get("Requirements"), list
    ):
        raise ValueError(f"{manifest_path} is not a version-1 requirements manifest")

    library_path = Path(manifest["Library"]["Path"]).expanduser()
    if not library_path.is_absolute():
        library_path = (manifest_path.parent / library_path).resolve()
    library = load_json(library_path)

    coupons = []
    for requirement in manifest["Requirements"]:
        if not args.include_matched and requirement.get("Status") != "Missing":
            continue
        coupon = {
            "Id": coupon_id(requirement),
            "Topology": requirement["Topology"],
            "Geometry": requirement.get("Geometry", {}),
            "Interfaces": requirement.get("Interfaces", []),
            "BoundaryCondition": requirement.get("BoundaryCondition", {}),
            "DeviceOccurrences": requirement.get("Count", 0),
            "DeviceEdgeLength": requirement.get("TotalEdgeLength", 0.0),
            "CoverageStatus": requirement.get("Status"),
            "Preparation": preparation(requirement),
        }
        if "SelectedModels" in requirement:
            coupon["SelectedModels"] = requirement["SelectedModels"]
        if "Reason" in requirement:
            coupon["CoverageReason"] = requirement["Reason"]
        coupons.append(coupon)

    coupons.sort(key=lambda coupon: coupon["Id"])
    duplicate_ids = {
        coupon["Id"] for coupon in coupons if sum(
            other["Id"] == coupon["Id"] for other in coupons
        ) > 1
    }
    for index, coupon in enumerate(coupons):
        if coupon["Id"] in duplicate_ids:
            coupon["Id"] += f"_{index + 1:03d}"

    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    plan = {
        "Version": 1,
        "SourceManifest": str(manifest_path),
        "ProcessLibrary": str(library_path),
        "Fabrication": library.get("Fabrication"),
        "MatchingRadius": manifest["Library"]["MatchingRadius"],
        "Coupons": coupons,
        "Summary": {
            "CouponCount": len(coupons),
            "AutomaticallyPreparatable": sum(
                coupon["Preparation"]["Method"] != "Manual" for coupon in coupons
            ),
            "Manual": sum(
                coupon["Preparation"]["Method"] == "Manual" for coupon in coupons
            ),
        },
    }
    destination = output / "coupon-plan.json"
    destination.write_text(json.dumps(plan, indent=2) + "\n")
    print(destination)
    if plan["Summary"]["Manual"]:
        print(
            f"Warning: {plan['Summary']['Manual']} spatial coupon requirement(s) "
            "still need a dedicated 3D mesh generator"
        )


if __name__ == "__main__":
    main()
