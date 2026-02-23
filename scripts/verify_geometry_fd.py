#!/usr/bin/env python3
"""Geometry-level FD verification for cylinder cavity eigenfrequency sensitivity.

Stretches the cylinder mesh in L (z-axis) and R (radial) directions by a small delta,
re-solves the eigenvalue problem, and computes df/dL and df/dR via forward differences.

For a uniform-fill cavity: (df/dL * L)/f + (df/dR * R)/f = -1.0 (energy partition).

Usage:
    cd examples/cylinder
    python3 ../../scripts/verify_geometry_fd.py

Requires: Palace binary at ../../build/bin/palace
"""

import json
import math
import os
import shutil
import struct
import subprocess
import sys

PALACE = os.path.join(os.path.dirname(__file__), "..", "build", "bin", "palace")
DELTA = 1e-4  # relative perturbation


def read_gmsh_binary_nodes(filename):
    """Read node coordinates from binary Gmsh v2.2 file."""
    with open(filename, "rb") as f:
        content = f.read()
    idx = content.find(b"$Nodes")
    idx = content.index(b"\n", idx) + 1
    end = content.index(b"\n", idx)
    nv = int(content[idx:end])
    idx = end + 1
    verts = []
    for _ in range(nv):
        struct.unpack_from("<i", content, idx)  # skip tag
        x, y, z = struct.unpack_from("<3d", content, idx + 4)
        verts.append((x, y, z))
        idx += 28
    return verts, content


def write_gmsh_scaled(content, output_file, scale_fn):
    """Write Gmsh binary with scaled node coordinates."""
    buf = bytearray(content)
    idx = buf.find(b"$Nodes")
    idx = buf.index(b"\n", idx) + 1
    end = buf.index(b"\n", idx)
    nv = int(buf[idx:end])
    idx = end + 1
    for _ in range(nv):
        x, y, z = struct.unpack_from("<3d", buf, idx + 4)
        x2, y2, z2 = scale_fn(x, y, z)
        struct.pack_into("<3d", buf, idx + 4, x2, y2, z2)
        idx += 28
    with open(output_file, "wb") as f:
        f.write(buf)


def run_palace(config_file):
    """Run Palace and return mode 1 eigenfrequency in GHz."""
    result = subprocess.run(
        [PALACE, config_file], capture_output=True, text=True, timeout=120
    )
    output_dir = None
    with open(config_file) as f:
        cfg = json.load(f)
    output_dir = cfg["Problem"]["Output"]
    with open(os.path.join(output_dir, "eig.csv")) as f:
        for line in f:
            if line.strip().startswith("1.00e+00"):
                return float(line.strip().split(",")[1].strip())
    raise RuntimeError("Could not parse eigenfrequency from output")


def main():
    mesh_file = "mesh/cylinder_hex.msh"
    base_config = "cavity_pec_sens.json"

    if not os.path.exists(mesh_file):
        print(f"Error: {mesh_file} not found. Run from examples/cylinder/", file=sys.stderr)
        sys.exit(1)

    # Read baseline config to get order and material properties
    with open(base_config) as f:
        cfg = json.load(f)

    # Run baseline
    print("Running baseline solve...")
    f_base = run_palace(base_config)
    print(f"  f_base = {f_base:.10f} GHz")

    # Read mesh geometry
    verts, content = read_gmsh_binary_nodes(mesh_file)
    zs = [v[2] for v in verts]
    rs = [math.sqrt(v[0] ** 2 + v[1] ** 2) for v in verts]
    L = max(zs) - min(zs)
    R = max(rs)
    print(f"  Cylinder: L = {L:.4f}, R = {R:.4f} (mesh units)")

    os.makedirs("mesh_temp", exist_ok=True)

    # Stretch L (z-axis)
    write_gmsh_scaled(
        content,
        "mesh_temp/stretched_z.msh",
        lambda x, y, z: (x, y, z * (1 + DELTA)),
    )
    cfg_L = dict(cfg)
    cfg_L["Problem"] = dict(cfg["Problem"], Output="postpro/fd_L")
    cfg_L["Model"] = dict(cfg["Model"], Mesh="mesh_temp/stretched_z.msh")
    # Remove sensitivity options for FD runs
    cfg_L["Solver"] = dict(cfg["Solver"])
    cfg_L["Solver"]["Eigenmode"] = {
        k: v
        for k, v in cfg["Solver"]["Eigenmode"].items()
        if k not in ("Sensitivity", "MeshSensitivity")
    }
    with open("fd_L.json", "w") as f:
        json.dump(cfg_L, f, indent=2)
    print("Running L-perturbed solve...")
    f_L = run_palace("fd_L.json")

    # Stretch R (radial)
    write_gmsh_scaled(
        content,
        "mesh_temp/stretched_r.msh",
        lambda x, y, z: (x * (1 + DELTA), y * (1 + DELTA), z),
    )
    cfg_R = dict(cfg)
    cfg_R["Problem"] = dict(cfg["Problem"], Output="postpro/fd_R")
    cfg_R["Model"] = dict(cfg["Model"], Mesh="mesh_temp/stretched_r.msh")
    cfg_R["Solver"] = dict(cfg["Solver"])
    cfg_R["Solver"]["Eigenmode"] = {
        k: v
        for k, v in cfg["Solver"]["Eigenmode"].items()
        if k not in ("Sensitivity", "MeshSensitivity")
    }
    with open("fd_R.json", "w") as f:
        json.dump(cfg_R, f, indent=2)
    print("Running R-perturbed solve...")
    f_R = run_palace("fd_R.json")

    # Compute FD derivatives
    dfdL = (f_L - f_base) / (L * DELTA)
    dfdR = (f_R - f_base) / (R * DELTA)
    ratio_L = dfdL * L / f_base
    ratio_R = dfdR * R / f_base

    print(f"\n{'=' * 60}")
    print(f"GEOMETRY FD VERIFICATION (delta = {DELTA})")
    print(f"{'=' * 60}")
    print(f"f_base    = {f_base:.10f} GHz")
    print(f"f(L+dL)   = {f_L:.10f} GHz")
    print(f"f(R+dR)   = {f_R:.10f} GHz")
    print(f"\ndf/dL     = {dfdL:+.6f} GHz/mesh_unit")
    print(f"df/dR     = {dfdR:+.6f} GHz/mesh_unit")
    print(f"\n(df/dL * L) / f = {ratio_L:+.6f}  (axial energy fraction)")
    print(f"(df/dR * R) / f = {ratio_R:+.6f}  (radial energy fraction)")
    print(f"Sum             = {ratio_L + ratio_R:+.6f}  (expect -1.000)")

    # Cleanup
    shutil.rmtree("mesh_temp", ignore_errors=True)
    for f in ["fd_L.json", "fd_R.json"]:
        if os.path.exists(f):
            os.remove(f)
    for d in ["postpro/fd_L", "postpro/fd_R"]:
        shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    main()
