#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Prepare fabrication process-ablation variants of one spatial coupon signature.

Given a generated coupon directory (containing coupon.json), emit ready-to-run
command scripts for isolated process effects at identical geometry/discretization:

  metal-only:    finite metal thickness, zero overetch
  overetch-only: nominal-zero metal thickness (10% of nominal), full overetch
  combined:      the production process

Comparing the three fabricated/thin responses identifies which fabrication
operation dominates a large local F/T gain.
"""

import argparse
import json
import shlex
from pathlib import Path


VARIANTS = {
    "metal-only": {"metal_scale": 1.0, "overetch_scale": 0.0},
    "overetch-only": {"metal_scale": 0.1, "overetch_scale": 1.0},
    "combined": {"metal_scale": 1.0, "overetch_scale": 1.0},
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("coupon", type=Path, help="Generated coupon directory")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True, help="Palace repository root")
    parser.add_argument("--python", default="python3")
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--julia-project", type=Path)
    parser.add_argument("--palace", type=Path, required=True)
    parser.add_argument("--ranks", type=int, default=64)
    parser.add_argument("--radius", type=float, required=True)
    parser.add_argument("--metal-thickness", type=float, default=0.1)
    parser.add_argument("--overetch-depth", type=float, default=0.05)
    parser.add_argument("--sidewall-angle", type=float, default=90.0)
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--ring-size", type=int, default=8)
    parser.add_argument("--lc-fine", type=float, default=0.025)
    parser.add_argument("--lc-far", type=float, default=0.4)
    parser.add_argument("--mesh-order", type=int, default=2)
    parser.add_argument("--substrate-permittivity", type=float, default=11.45)
    parser.add_argument("--sa-thickness", type=float, default=0.002)
    parser.add_argument("--sa-permittivity", type=float, default=4.0)
    parser.add_argument("--ms-thickness", type=float, default=0.002)
    parser.add_argument("--ms-permittivity", type=float, default=11.45)
    parser.add_argument("--ma-thickness", type=float, default=0.002)
    parser.add_argument("--ma-permittivity", type=float, default=10.0)
    args = parser.parse_args()

    coupon = args.coupon.expanduser().resolve()
    coupon_json = coupon / "coupon.json"
    if not coupon_json.is_file():
        raise FileNotFoundError(coupon_json)
    coupon_id = json.loads(coupon_json.read_text())["Id"]
    repo = args.repo.expanduser().resolve()
    generator = repo / "examples/cpw3d_surface/spatial_coupon/generate_spatial_response.py"
    mesher = repo / "examples/cpw3d_surface/spatial_coupon/mesh_spatial_coupon.jl"
    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)

    script_lines = ["#!/bin/bash", "set -euo pipefail", ""]
    manifest = {"Version": 1, "Coupon": str(coupon), "Id": coupon_id, "Variants": {}}
    for name, scales in VARIANTS.items():
        metal = args.metal_thickness * scales["metal_scale"]
        overetch = args.overetch_depth * scales["overetch_scale"]
        root = output / name
        root.mkdir(parents=True, exist_ok=True)
        generate = [
            args.python,
            str(generator),
            str(coupon_json),
            "--output", str(root),
            "--radius", args.radius,
            "--metal-thickness", metal,
            "--overetch-depth", overetch,
            "--sidewall-angle", args.sidewall_angle,
            "--top-rounding", 0.0,
            "--trench-rounding", 0.0,
            "--ring-size", args.ring_size,
            "--order", args.order,
            "--model-name", f"{coupon_id[:40]}-{name}",
            "--substrate-permittivity", args.substrate_permittivity,
            "--sa-thickness", args.sa_thickness,
            "--sa-permittivity", args.sa_permittivity,
            "--ms-thickness", args.ms_thickness,
            "--ms-permittivity", args.ms_permittivity,
            "--ma-thickness", args.ma_thickness,
            "--ma-permittivity", args.ma_permittivity,
        ]
        signature = [str(value) for value in generate + ["--signature-only"]]
        julia_project = (
            ["--project=" + str(args.julia_project.expanduser().resolve())]
            if args.julia_project
            else []
        )
        meshes = []
        for kind in ("thin", "fabricated"):
            mesh = [
                args.julia,
                *julia_project,
                str(mesher),
                str(root / "mesh-signature.csv"),
                kind,
                str(root / f"spatial_{kind}.msh"),
                "--mask", str(root / "plan-view-mask.csv"),
                "--boundary", str(root / "plan-view-boundary.csv"),
                "--radius", args.radius,
                "--metal-thickness", metal,
                "--overetch", overetch,
                "--sidewall-angle", args.sidewall_angle,
                "--top-radius", 0.0,
                "--bottom-radius", 0.0,
                "--lc-fine", args.lc_fine,
                "--lc-far", args.lc_far,
                "--mesh-order", args.mesh_order,
            ]
            meshes.append([str(value) for value in mesh])
        full = [
            str(value)
            for value in generate
            + [
                "--thin-mesh", str(root / "spatial_thin.msh"),
                "--fabricated-mesh", str(root / "spatial_fabricated.msh"),
            ]
        ]
        solves = [
            [
                str(args.palace.expanduser().resolve()),
                "-np", str(args.ranks),
                str(root / f"spatial_{kind}.json"),
            ]
            for kind in ("thin", "fabricated")
        ]
        script_lines.append(f"# ---- {name}: metal={metal} overetch={overetch}")
        for command in [signature, *meshes, full, *solves]:
            script_lines.append(shlex.join(command))
        script_lines.append("")
        manifest["Variants"][name] = {
            "MetalThickness": metal,
            "OveretchDepth": overetch,
            "Root": str(root),
        }

    script = output / "run-process-ablation.sh"
    script.write_text("\n".join(script_lines) + "\n")
    script.chmod(0o755)
    (output / "process-ablation.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(script)
    print(output / "process-ablation.json")


if __name__ == "__main__":
    main()
