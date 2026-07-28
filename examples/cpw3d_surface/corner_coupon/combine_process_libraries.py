#!/usr/bin/env python3

"""Combine response models into one portable fabrication-process library."""

import argparse
import json
import math
import re
import shutil
from pathlib import Path


PATH_FIELDS = (
    "FabricatedMatrix",
    "ThinMatrix",
    "FabricatedSurfaceMatrix",
    "ThinSurfaceMatrix",
    "BasisPoints",
)

PATH_NAMES = {
    "FabricatedMatrix": "fabricated-domain-response-matrix.csv",
    "ThinMatrix": "thin-domain-response-matrix.csv",
    "FabricatedSurfaceMatrix": "fabricated-surface-response-matrix.csv",
    "ThinSurfaceMatrix": "thin-surface-response-matrix.csv",
    "BasisPoints": "basis-points.csv",
}


def model_directory(index, name):
    slug = re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-")
    return Path("models") / f"{index:03d}-{slug}"


def load_library(path):
    path = path.expanduser().resolve()
    data = json.loads(path.read_text())
    if data.get("Version") not in (1, 2, 3):
        raise ValueError(f"{path} is not a supported process library")
    if data["Version"] >= 3:
        layers = data.get("Fabrication", {}).get("InterfaceLayers")
        if not isinstance(layers, dict):
            raise ValueError(
                f"{path} is version 3 but has no Fabrication.InterfaceLayers metadata"
            )
    models = data.get("Models")
    if not isinstance(models, list):
        raise ValueError(f"{path} has no Models array")
    radius = float(data["MatchingRadius"])
    if radius <= 0.0:
        raise ValueError(f"{path} has a nonpositive MatchingRadius")
    return path, data


def merge_metadata(first, second, path="Fabrication"):
    if isinstance(first, dict) and isinstance(second, dict):
        result = {}
        for key in first.keys() | second.keys():
            if key not in first:
                result[key] = second[key]
            elif key not in second:
                result[key] = first[key]
            else:
                result[key] = merge_metadata(
                    first[key], second[key], f"{path}.{key}"
                )
        return result
    if (
        isinstance(first, (int, float))
        and not isinstance(first, bool)
        and isinstance(second, (int, float))
        and not isinstance(second, bool)
    ):
        if math.isclose(float(first), float(second), rel_tol=1.0e-12, abs_tol=0.0):
            return first
    elif first == second:
        return first
    raise ValueError(f"Input libraries contain different {path} metadata")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--name", default="combined-fabrication-process")
    parser.add_argument("libraries", type=Path, nargs="+")
    args = parser.parse_args()

    destination = args.output.expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    loaded = [load_library(path) for path in args.libraries]
    matching_radius = float(loaded[0][1]["MatchingRadius"])
    tolerance = 1.0e-12 * matching_radius
    for path, library in loaded[1:]:
        if abs(float(library["MatchingRadius"]) - matching_radius) > tolerance:
            raise ValueError(
                f"{path} uses MatchingRadius={library['MatchingRadius']}; "
                f"expected {matching_radius}"
            )

    version = max(library["Version"] for _, library in loaded)
    fabrication = [library.get("Fabrication") for _, library in loaded]
    known_fabrication = [entry for entry in fabrication if entry is not None]
    merged_fabrication = None
    for entry in known_fabrication:
        merged_fabrication = (
            entry
            if merged_fabrication is None
            else merge_metadata(merged_fabrication, entry)
        )
    if known_fabrication and len(known_fabrication) != len(fabrication):
        print(
            "Warning: combining libraries with and without Fabrication metadata; "
            "the output is downgraded to version 2 until its process is recorded"
        )
        version = min(version, 2)

    result = {
        "Version": version,
        "Name": args.name,
        "MatchingRadius": matching_radius,
        "Models": [],
    }
    if known_fabrication and len(known_fabrication) == len(fabrication):
        result["Fabrication"] = merged_fabrication
    names = set()
    index = 0
    for source_path, library in loaded:
        source_root = source_path.parent
        default_depth = library.get("CouponDepth")
        for source_model in library["Models"]:
            index += 1
            model = dict(source_model)
            name = model.get("Name", "")
            if not name or name in names:
                raise ValueError(
                    f"Response model name {name!r} is empty or repeated"
                )
            names.add(name)
            if default_depth is not None and "CouponDepth" not in model:
                model["CouponDepth"] = default_depth

            relative_directory = model_directory(index, name)
            model_destination = destination / relative_directory
            model_destination.mkdir(parents=True, exist_ok=True)
            for field in PATH_FIELDS:
                if field not in model:
                    continue
                source = Path(model[field])
                if not source.is_absolute():
                    source = source_root / source
                source = source.resolve()
                if not source.is_file():
                    raise FileNotFoundError(
                        f"{name} field {field} does not exist: {source}"
                    )
                target = model_destination / PATH_NAMES[field]
                shutil.copy2(source, target)
                model[field] = str(relative_directory / target.name)
            result["Models"].append(model)

    if not result["Models"]:
        raise ValueError("Combined library contains no response models")

    output = destination / "process-library.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"Wrote {len(result['Models'])} models at R={matching_radius:g} "
        f"to {output}"
    )


if __name__ == "__main__":
    main()
