#!/usr/bin/env python3

"""Prepare the coarse transmon for geometry-only surface-response preflight."""

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "transmon_surface_coarse.json",
    )
    parser.add_argument(
        "--library",
        type=Path,
        default=ROOT / "transmon_surface_process_seed.json",
    )
    parser.add_argument("--postpro", type=Path)
    args = parser.parse_args()

    source = args.config.expanduser().resolve()
    library = args.library.expanduser().resolve()
    output = args.output.expanduser().resolve()
    postpro = (
        args.postpro.expanduser().resolve()
        if args.postpro
        else output.parent / "transmon-surface-preflight"
    )

    config = json.loads(source.read_text())
    mesh = Path(config["Model"]["Mesh"])
    if not mesh.is_absolute():
        mesh = source.parent / mesh
    config["Model"]["Mesh"] = str(mesh.resolve())
    config["Problem"]["Output"] = str(postpro)
    config["Solver"]["SurfaceResponseCorrection"] = {
        "Library": str(library),
        "TargetInterfaces": [1, 2, 3],
        "UnmatchedPolicy": "Warn",
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(config, indent=2) + "\n")
    print(output)
    print(f"palace --surface-response-preflight {output}")


if __name__ == "__main__":
    main()
