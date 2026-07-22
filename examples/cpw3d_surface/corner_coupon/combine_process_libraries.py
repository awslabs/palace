#!/usr/bin/env python3

"""Combine response models into one portable fabrication-process library."""

import argparse
import json
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
    if data.get("Version") != 1:
        raise ValueError(f"{path} is not a version-1 process library")
    if not data.get("Models"):
        raise ValueError(f"{path} contains no response models")
    radius = float(data["MatchingRadius"])
    if radius <= 0.0:
        raise ValueError(f"{path} has a nonpositive MatchingRadius")
    return path, data


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

    result = {
        "Version": 1,
        "Name": args.name,
        "MatchingRadius": matching_radius,
        "Models": [],
    }
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

    output = destination / "process-library.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"Wrote {len(result['Models'])} models at R={matching_radius:g} "
        f"to {output}"
    )


if __name__ == "__main__":
    main()
