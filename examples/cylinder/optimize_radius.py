#!/usr/bin/env python3
"""
Gradient-based optimization of cylinder radius to match a target eigenfrequency.

Uses Palace's analytic mesh sensitivity (dual-number QFunction differentiation)
to compute df/dR, then performs gradient descent on the radius parameter.

The mesh is morphed by scaling (x,y) coordinates proportionally to R_new/R_old.
This is the simplest possible shape optimization loop — one design parameter,
analytical mesh deformation, no remeshing.

Usage:
    python3 optimize_radius.py --target 3.0 --palace /path/to/palace-arm64.bin
"""

import argparse
import json
import os
import re
import shutil
import struct
import subprocess
import sys


def read_gmsh_binary_nodes(msh_path):
    """Read node coordinates from a binary Gmsh v2 .msh file."""
    with open(msh_path, "rb") as f:
        content = f.read()
    marker = b"$Nodes\n"
    idx = content.find(marker) + len(marker)
    nl = content.find(b"\n", idx)
    n_nodes = int(content[idx:nl].decode().strip().split()[0])
    idx = nl + 1
    nodes = {}
    for _ in range(n_nodes):
        tag = struct.unpack_from("<i", content, idx)[0]
        x, y, z = struct.unpack_from("<3d", content, idx + 4)
        nodes[tag] = (x, y, z)
        idx += 4 + 24
    return nodes, content, n_nodes


def write_scaled_mesh(original_msh, output_msh, scale_xy):
    """Write a mesh with (x,y) coordinates scaled by scale_xy."""
    with open(original_msh, "rb") as f:
        content = bytearray(f.read())
    marker = b"$Nodes\n"
    idx = content.find(marker) + len(marker)
    nl = content.find(b"\n", idx)
    n_nodes = int(content[idx:nl].decode().strip().split()[0])
    idx = nl + 1
    for _ in range(n_nodes):
        # Skip tag (4 bytes)
        idx += 4
        x, y, z = struct.unpack_from("<3d", content, idx)
        struct.pack_into("<3d", content, idx, x * scale_xy, y * scale_xy, z)
        idx += 24
    with open(output_msh, "wb") as f:
        f.write(content)


def run_palace(palace_bin, config_path, work_dir):
    """Run Palace and extract eigenfrequency and df/dR from output."""
    result = subprocess.run(
        [palace_bin, config_path],
        cwd=work_dir,
        capture_output=True,
        text=True,
        timeout=300,
    )
    output = result.stdout + result.stderr

    # Extract first eigenfrequency
    freq_match = re.search(r"Re\{f\} \(GHz\).*?\n\s+1,\s+([\+\-\d.eE]+)", output)
    if not freq_match:
        print("ERROR: Could not find eigenfrequency in Palace output")
        print(output[-500:])
        return None, None
    f_ghz = float(freq_match.group(1))

    # Extract df/dR from chain contraction diagnostic
    chain_match = re.search(r"Chain contraction.*?df/dR = ([\+\-\d.eE]+) GHz/cm", output)
    if not chain_match:
        print("WARNING: Could not find chain contraction in output")
        return f_ghz, None
    df_dR = float(chain_match.group(1))

    return f_ghz, df_dR


def main():
    parser = argparse.ArgumentParser(description="Optimize cylinder radius for target frequency")
    parser.add_argument("--target", type=float, default=3.0, help="Target frequency in GHz")
    parser.add_argument("--palace", type=str,
                        default="../../build/palace-build/palace-arm64.bin",
                        help="Path to Palace binary")
    parser.add_argument("--config", type=str, default="cavity_pec_meshsens.json",
                        help="Palace config file")
    parser.add_argument("--mesh", type=str, default="mesh/cylinder_hex.msh",
                        help="Original mesh file")
    parser.add_argument("--max-iter", type=int, default=10, help="Maximum iterations")
    parser.add_argument("--tol", type=float, default=1e-4, help="Convergence tolerance (GHz)")
    args = parser.parse_args()

    work_dir = os.path.dirname(os.path.abspath(args.config)) or "."
    mesh_original = os.path.join(work_dir, args.mesh)
    mesh_working = os.path.join(work_dir, "mesh/cylinder_hex_opt.msh")

    # Modify config to use the working mesh
    with open(os.path.join(work_dir, args.config)) as f:
        config = json.load(f)
    config["Model"]["Mesh"] = "mesh/cylinder_hex_opt.msh"
    config_opt = os.path.join(work_dir, "cavity_opt.json")
    with open(config_opt, "w") as f:
        json.dump(config, f, indent=2)

    # Read original mesh to get R
    nodes, _, _ = read_gmsh_binary_nodes(mesh_original)
    import math
    R_original = max(math.sqrt(v[0]**2 + v[1]**2) for v in nodes.values())
    R_current = R_original
    f_target = args.target

    print(f"=== Gradient-Based Cylinder Radius Optimization ===")
    print(f"Target frequency: {f_target:.6f} GHz")
    print(f"Initial radius:   {R_original:.6f} cm")
    print(f"Tolerance:        {args.tol:.1e} GHz")
    print()

    for iteration in range(args.max_iter):
        # Scale mesh
        scale = R_current / R_original
        write_scaled_mesh(mesh_original, mesh_working, scale)

        # Run Palace
        f_ghz, df_dR = run_palace(args.palace, "cavity_opt.json", work_dir)
        if f_ghz is None:
            print("Palace run failed!")
            return 1

        residual = f_ghz - f_target
        print(f"Iter {iteration:2d}: R = {R_current:.6f} cm, f = {f_ghz:.6f} GHz, "
              f"residual = {residual:+.6e} GHz, df/dR = {df_dR:.6e} GHz/cm")

        if abs(residual) < args.tol:
            print(f"\nConverged! Final R = {R_current:.6f} cm, f = {f_ghz:.6f} GHz")
            return 0

        if df_dR is None or abs(df_dR) < 1e-15:
            print("No gradient available!")
            return 1

        # Newton step: R_new = R_old - (f - f_target) / (df/dR)
        delta_R = -residual / df_dR
        # Limit step size to 10% of current R
        max_step = 0.1 * R_current
        if abs(delta_R) > max_step:
            delta_R = max_step * (1 if delta_R > 0 else -1)
        R_current += delta_R

    print(f"\nDid not converge in {args.max_iter} iterations")
    return 1


if __name__ == "__main__":
    sys.exit(main())
