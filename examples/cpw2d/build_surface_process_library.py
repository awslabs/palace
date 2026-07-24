#!/usr/bin/env python3

"""Build and package local surface responses for one fabrication process."""

import argparse
import json
import math
import shlex
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
EDGE_MESH_GENERATOR = ROOT / "mesh" / "mesh_edge_coupon.jl"
EDGE_RESPONSE_GENERATOR = ROOT / "generate_edge_response.py"
PAIR_MESH_GENERATOR = ROOT / "mesh" / "mesh_edge_pair_coupon.jl"
PAIR_RESPONSE_GENERATOR = ROOT / "generate_edge_pair_response.py"
LIBRARY_COMBINER = (
    ROOT.parent
    / "cpw3d_surface"
    / "corner_coupon"
    / "combine_process_libraries.py"
)
MATRIX_FILES = (
    "domain-response-matrix.csv",
    "surface-response-matrix.csv",
)
TOPOLOGIES = (
    ("same-conductor-gap", "same", None, "same_conductor_gaps"),
    (
        "different-conductor-gap",
        "different",
        "--different-conductors",
        "different_conductor_gaps",
    ),
    ("same-conductor-strip", "strip", "--strip", "same_conductor_strips"),
)


def run(command):
    print("+ " + shlex.join(str(value) for value in command), flush=True)
    subprocess.run([str(value) for value in command], check=True)


def width_slug(width):
    return f"{width:.12g}".replace(".", "p")


def normalize_widths(widths, matching_radius, name):
    if widths is None:
        return []
    result = sorted(set(widths))
    if any(width <= 0.0 for width in result):
        raise ValueError(f"{name} separations must be positive")
    maximum = 2.0 * matching_radius
    if any(width > maximum for width in result):
        raise ValueError(
            f"{name} separations cannot exceed 2R = {maximum:g} um; "
            "wider neighborhoods use the isolated-edge model"
        )
    return result


def response_complete(directory, solves):
    return all(
        (directory / "postpro" / solve / filename).is_file()
        for solve in solves
        for filename in MATRIX_FILES
    )


def load_json(path):
    try:
        return json.loads(path.read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return None


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2) + "\n")


def julia_command(args):
    command = [args.julia]
    if args.julia_project:
        command.append(f"--project={args.julia_project.expanduser().resolve()}")
    return command


def process_mesh_options(args):
    return [
        "--radius",
        args.matching_radius,
        "--metal-thickness",
        args.metal_thickness,
        "--overetch",
        args.overetch,
        "--sidewall-angle",
        args.sidewall_angle,
        "--top-radius",
        args.top_radius,
        "--bottom-radius",
        args.bottom_radius,
        "--lc-fine",
        args.lc_fine,
        "--lc-far",
        args.lc_far,
        "--mesh-order",
        args.mesh_order,
    ]


def material_options(args):
    return [
        "--substrate-permittivity",
        args.substrate_permittivity,
        "--sa-thickness",
        args.sa_thickness,
        "--sa-permittivity",
        args.sa_permittivity,
        "--ms-thickness",
        args.ms_thickness,
        "--ms-permittivity",
        args.ms_permittivity,
        "--ma-thickness",
        args.ma_thickness,
        "--ma-permittivity",
        args.ma_permittivity,
    ]


def pair_mesh_command(args, kind, width, mode, output):
    return julia_command(args) + [
        PAIR_MESH_GENERATOR,
        kind,
        width,
        output,
        mode,
        *process_mesh_options(args),
    ]


def pair_response_command(args, width, topology_flag, directory):
    command = [
        sys.executable,
        PAIR_RESPONSE_GENERATOR,
        "--separation",
        width,
        "--radius",
        args.matching_radius,
        "--metal-thickness",
        args.metal_thickness,
        "--basis-size",
        args.basis_size,
        "--samples",
        args.samples,
        "--order",
        args.order,
        "--coupon-depth",
        args.coupon_depth,
        "--separation-tolerance",
        args.separation_tolerance,
        "--library-name",
        args.name,
        "--output",
        directory,
        *material_options(args),
    ]
    if topology_flag:
        command.append(topology_flag)
    return command


def edge_mesh_command(args, kind, output):
    return julia_command(args) + [
        EDGE_MESH_GENERATOR,
        kind,
        output,
        *process_mesh_options(args),
    ]


def edge_response_command(args, directory):
    return [
        sys.executable,
        EDGE_RESPONSE_GENERATOR,
        "--radius",
        args.matching_radius,
        "--metal-thickness",
        args.metal_thickness,
        "--basis-size",
        args.basis_size,
        "--samples",
        args.samples,
        "--order",
        args.order,
        "--coupon-depth",
        args.coupon_depth,
        "--library-name",
        args.name,
        "--output",
        directory,
        *material_options(args),
    ]


def palace_command(args, config):
    command = [args.palace]
    command.extend(["--serial"] if args.ranks == 1 else ["-np", args.ranks])
    command.append(config)
    return command


def coupon_spec(args, topology, width):
    return {
        "Version": 1,
        "Topology": topology,
        "Separation": width,
        "MatchingRadius": args.matching_radius,
        "Fabrication": {
            "MetalThickness": args.metal_thickness,
            "OveretchDepth": args.overetch,
            "SidewallAngle": args.sidewall_angle,
            "TopRadius": args.top_radius,
            "BottomRadius": args.bottom_radius,
            "SubstratePermittivity": args.substrate_permittivity,
            "InterfaceLayers": {
                "SA": {
                    "Thickness": args.sa_thickness,
                    "Permittivity": args.sa_permittivity,
                },
                "MS": {
                    "Thickness": args.ms_thickness,
                    "Permittivity": args.ms_permittivity,
                },
                "MA": {
                    "Thickness": args.ma_thickness,
                    "Permittivity": args.ma_permittivity,
                },
            },
        },
        "Mesh": {
            "FineSize": args.lc_fine,
            "FarSize": args.lc_far,
            "Order": args.mesh_order,
        },
        "Response": {
            "BasisSize": args.basis_size,
            "Samples": args.samples,
            "Order": args.order,
            "CouponDepth": args.coupon_depth,
            "SeparationTolerance": args.separation_tolerance,
        },
    }


def isolated_spec(args):
    spec = coupon_spec(args, "isolated-edge", 0.0)
    del spec["Separation"]
    del spec["Response"]["SeparationTolerance"]
    return spec


def build_isolated(args, work):
    directory = work / "isolated-edge"
    directory.mkdir(parents=True, exist_ok=True)
    thin_mesh = directory / "edge_thin.msh"
    fabricated_mesh = directory / "edge_fabricated.msh"
    spec_path = directory / "coupon-spec.json"
    spec = isolated_spec(args)
    reusable = load_json(spec_path) == spec

    if not reusable or not thin_mesh.is_file():
        run(edge_mesh_command(args, "thin", thin_mesh))
    if not reusable or not fabricated_mesh.is_file():
        run(edge_mesh_command(args, "fabricated", fabricated_mesh))

    run(edge_response_command(args, directory))
    write_json(spec_path, spec)
    configs = [directory / "edge_thin.json", directory / "edge_fabricated.json"]
    if args.prepare_only:
        return directory / "process-library.json", configs, False

    solves = ("edge_thin", "edge_fabricated")
    if reusable and response_complete(directory, solves):
        print("Reusing completed isolated-edge response")
    else:
        for config in configs:
            run(palace_command(args, config))
        if not response_complete(directory, solves):
            raise RuntimeError("Palace did not write all isolated-edge response matrices")
    return directory / "process-library.json", configs, True


def build_coupon(args, topology, mode, topology_flag, width, work):
    directory = work / topology / f"w{width_slug(width)}um"
    directory.mkdir(parents=True, exist_ok=True)
    thin_mesh = directory / "edge_pair_thin.msh"
    fabricated_mesh = directory / "edge_pair_fabricated.msh"
    spec_path = directory / "coupon-spec.json"
    spec = coupon_spec(args, topology, width)
    reusable = load_json(spec_path) == spec

    if not reusable or not thin_mesh.is_file():
        run(pair_mesh_command(args, "thin", width, mode, thin_mesh))
    if not reusable or not fabricated_mesh.is_file():
        run(pair_mesh_command(args, "fabricated", width, mode, fabricated_mesh))

    run(pair_response_command(args, width, topology_flag, directory))
    write_json(spec_path, spec)
    configs = [
        directory / "edge_pair_thin.json",
        directory / "edge_pair_fabricated.json",
    ]
    if args.prepare_only:
        return directory / "process-library.json", configs, False

    solves = ("edge_pair_thin", "edge_pair_fabricated")
    if reusable and response_complete(directory, solves):
        print(f"Reusing completed {topology} response at w={width:g} um")
    else:
        for config in configs:
            run(palace_command(args, config))
        if not response_complete(directory, solves):
            raise RuntimeError(
                f"Palace did not write all response matrices for {topology} "
                f"at w={width:g} um"
            )
    return directory / "process-library.json", configs, True


def add_process_metadata(path, args, width_sets):
    library = json.loads(path.read_text())
    library["Version"] = max(3, int(library.get("Version", 0)))
    library["Fabrication"] = {
        "LengthUnit": "um",
        "MetalThickness": args.metal_thickness,
        "OveretchDepth": args.overetch,
        "SidewallAngleDegrees": args.sidewall_angle,
        "TopRoundingRadius": args.top_radius,
        "BottomRoundingRadius": args.bottom_radius,
        "SubstratePermittivity": args.substrate_permittivity,
        "InterfaceLayers": {
            "SA": {
                "Thickness": args.sa_thickness,
                "Permittivity": args.sa_permittivity,
            },
            "MS": {
                "Thickness": args.ms_thickness,
                "Permittivity": args.ms_permittivity,
            },
            "MA": {
                "Thickness": args.ma_thickness,
                "Permittivity": args.ma_permittivity,
            },
        },
    }
    library["PairedEdgeSampling"] = width_sets
    write_json(path, library)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate, solve, and package isolated- and paired-edge responses. "
            "All geometric dimensions are in microns."
        )
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--palace", type=Path)
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--julia-project", type=Path)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--matching-radius", type=float, default=2.0)
    parser.add_argument("--metal-thickness", type=float, default=0.1)
    parser.add_argument("--overetch", type=float, default=0.05)
    parser.add_argument("--sidewall-angle", type=float, default=80.0)
    parser.add_argument("--top-radius", type=float, default=0.01)
    parser.add_argument("--bottom-radius", type=float, default=0.01)
    parser.add_argument("--lc-fine", type=float, default=0.002)
    parser.add_argument("--lc-far", type=float, default=0.05)
    parser.add_argument("--mesh-order", type=int, default=2)
    parser.add_argument("--basis-size", type=int, default=96)
    parser.add_argument("--samples", type=int, default=1200)
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--coupon-depth", type=float, default=1055.0)
    parser.add_argument("--substrate-permittivity", type=float, default=11.47)
    parser.add_argument("--sa-thickness", type=float, default=0.002)
    parser.add_argument("--sa-permittivity", type=float, default=4.0)
    parser.add_argument("--ms-thickness", type=float, default=0.002)
    parser.add_argument("--ms-permittivity", type=float, default=11.47)
    parser.add_argument("--ma-thickness", type=float, default=0.002)
    parser.add_argument("--ma-permittivity", type=float, default=10.0)
    parser.add_argument("--separation-tolerance", type=float, default=1.0e-3)
    parser.add_argument(
        "--separations",
        type=float,
        nargs="+",
        help="Default separation sweep for all three paired topologies",
    )
    parser.add_argument("--same-conductor-gaps", type=float, nargs="*")
    parser.add_argument("--different-conductor-gaps", type=float, nargs="*")
    parser.add_argument("--same-conductor-strips", type=float, nargs="*")
    parser.add_argument(
        "--base-library",
        type=Path,
        action="append",
        default=[],
        help="Existing isolated-edge or corner library to include",
    )
    parser.add_argument(
        "--skip-isolated",
        action="store_true",
        help="Do not generate an isolated-edge response",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Generate meshes and configs without running Palace or packaging matrices",
    )
    args = parser.parse_args()

    if args.matching_radius <= 0.0:
        parser.error("--matching-radius must be positive")
    if args.ranks <= 0:
        parser.error("--ranks must be positive")
    if args.basis_size <= 0 or args.basis_size % 2:
        parser.error("--basis-size must be a positive even integer")
    if args.samples < args.basis_size:
        parser.error("--samples must be at least --basis-size")
    if args.mesh_order <= 0 or args.order <= 0:
        parser.error("FEM orders must be positive")
    material_values = (
        args.substrate_permittivity,
        args.sa_thickness,
        args.sa_permittivity,
        args.ms_thickness,
        args.ms_permittivity,
        args.ma_thickness,
        args.ma_permittivity,
    )
    if any(not math.isfinite(value) or value <= 0.0 for value in material_values):
        parser.error("substrate and interface-layer properties must be finite and positive")
    if not args.prepare_only and args.palace is None:
        parser.error("--palace is required unless --prepare-only is used")

    for _, _, _, attribute in TOPOLOGIES:
        specific = getattr(args, attribute)
        widths = specific if specific is not None else args.separations
        try:
            setattr(
                args,
                attribute,
                normalize_widths(widths, args.matching_radius, attribute),
            )
        except ValueError as error:
            parser.error(str(error))
    if (
        args.skip_isolated
        and not args.base_library
        and not any(getattr(args, attribute) for _, _, _, attribute in TOPOLOGIES)
    ):
        parser.error("no response model or base library was requested")
    return args


def library_has_topology(path, topology):
    path = path.expanduser().resolve()
    library = load_json(path)
    if library is None:
        raise ValueError(f"Unable to read base process library: {path}")
    return any(model.get("Topology") == topology for model in library.get("Models", []))


def main():
    args = parse_args()
    output = args.output.expanduser().resolve()
    work = output / "coupons"
    destination = output / "library"
    work.mkdir(parents=True, exist_ok=True)
    libraries = [path.expanduser().resolve() for path in args.base_library]
    manifest = {
        "Version": 1,
        "Name": args.name,
        "MatchingRadius": args.matching_radius,
        "PreparedOnly": args.prepare_only,
        "Coupons": [],
    }
    width_sets = {}

    base_has_isolated = any(
        library_has_topology(path, "IsolatedEdge") for path in libraries
    )
    if not args.skip_isolated and not base_has_isolated:
        library, configs, complete = build_isolated(args, work)
        libraries.append(library)
        manifest["Coupons"].append(
            {
                "Topology": "isolated-edge",
                "Library": str(library),
                "Configs": [str(path) for path in configs],
                "Complete": complete,
            }
        )
        write_json(output / "build-manifest.json", manifest)
    elif base_has_isolated:
        print("Reusing IsolatedEdge from a base process library")

    for topology, mode, topology_flag, attribute in TOPOLOGIES:
        widths = getattr(args, attribute)
        width_sets[topology] = widths
        for width in widths:
            library, configs, complete = build_coupon(
                args, topology, mode, topology_flag, width, work
            )
            libraries.append(library)
            manifest["Coupons"].append(
                {
                    "Topology": topology,
                    "Separation": width,
                    "Library": str(library),
                    "Configs": [str(path) for path in configs],
                    "Complete": complete,
                }
            )
            write_json(output / "build-manifest.json", manifest)

    if args.prepare_only:
        print(
            f"Prepared {len(manifest['Coupons'])} response coupons under {work}"
        )
        print("Run this command again with --palace and without --prepare-only to solve")
        return

    run(
        [
            sys.executable,
            LIBRARY_COMBINER,
            "--output",
            destination,
            "--name",
            args.name,
            *libraries,
        ]
    )
    library_path = destination / "process-library.json"
    add_process_metadata(library_path, args, width_sets)
    manifest["Library"] = str(library_path)
    write_json(output / "build-manifest.json", manifest)
    print(f"Complete fabrication-process library: {library_path}")
    if not args.base_library:
        print(
            "Warning: no --base-library was supplied; this library contains straight "
            "edge models only and cannot replace detected corner neighborhoods."
        )


if __name__ == "__main__":
    main()
