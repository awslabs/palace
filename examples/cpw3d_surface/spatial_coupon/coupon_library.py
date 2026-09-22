#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""coupon-library: the consolidated command line of the coupon mesh library
(supervisor decision 48): `build` (the mesh path) and `qualify` (the physics path).

  coupon_library.py build [--register CASE_ID=SOURCE_DIR ... --footprint {bound,producer-default}
                           --inventory-status STATUS [--mesh-recipe PATH] [--provenance TEXT]] [--no-thin]
                          [--device PALACE_CONFIG --palace PATH [--device-output DIR] [--ring-size N]
                           [--cap-triangulation METHOD]]
                          [--case ID ...] [--jobs N] [--build-limit N] [--root DIR] [--output PATH] [--manifest PATH]
  coupon_library.py qualify --build-record library-build.json --reference <campaign dir or none>
                            --remote HOST:ROOT [--orders p5] --controls p3,p5 --max-jobs N
                            [--job-policy speed|frugal|fixed [--fixed-jobs N]]
                            [--frozen-binary-sha256 HEX] [--reducer-block-size N] [--case ID ...] [--root DIR] [--dry-run]

`build --device` maps a device layout to coupon source directories first
(device_coupons.py: discovery closure by Palace geometry preflights -> the planner's
spatial coupons -> generate_spatial_response.py --basis-only -> content-hashed source
directories with provenance; the other families are recorded out of scope) and
registers every one of them (footprint producer-default, InventoryStatus DeviceDerived,
the manifest's shared mesh recipe; the device basis's box caps are Delaunay-triangulated
by default - decision 57 -, `--cap-triangulation ear-clipping` selects the gallery
producer's caps); `--build-limit N` builds the N smallest by the
pre-build estimate and records the rest registered-unbuilt with their estimates.
`build --register` registers the given source directories as manifest cases (register_case.py:
source SHA256s, the automated two-pass contract derivation, idempotent by content; the
footprint declaration is mandatory and applies to every directory registered by the
call).  Both then build the selected cases (default: every case of the manifest, the newly
registered ones included) as a job pool (run_gmsh_only_matrix.py: headroom gate ->
run_gmsh_only_case.py DAG -> audits -> verify_canonical_case_entries.py, fail closed per
case) and writes library-build.json.  A registration that fails closed stops the build
before any mesh is made (exit 1) and leaves its record under the case's work directory.

`qualify` runs the physics of every passed coupon of a library-build.json, against its
graded_v2 reference or on its own (--reference none) (qualify/qualify_library.py: the
run config derived from the case's own sources at the recipe's PhysicsRun Order / Tol
on the hash-verified identity mesh, estimate gate, plan with pinned digests, submission
under the user job cap - a coupon's main-order sources split into N worker jobs plus one
reducer job on the archive union under the recorded job policy (speed / frugal / fixed;
decision 61b) -, read-only monitoring, fetch / digest verification / matrix
validation / recorded archive deletion, the frozen machine-readable class gates) and
writes library-qualification.json, qualification-gates.json and process-library.json;
--dry-run writes the plans / configs / estimates / gates without contacting anything.
"""
import argparse
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import device_coupons  # noqa: E402
import register_case  # noqa: E402
import run_gmsh_only_matrix  # noqa: E402
sys.path.insert(0, str(HERE / "qualify"))
import qualify_library  # noqa: E402


def parse_registration(value):
    if "=" not in value:
        raise argparse.ArgumentTypeError("--register expects CASE_ID=SOURCE_DIR")
    case_id, directory = value.split("=", 1)
    if not case_id or not directory:
        raise argparse.ArgumentTypeError("--register expects CASE_ID=SOURCE_DIR")
    return case_id, Path(directory)


def build_parser():
    parser = argparse.ArgumentParser(prog="coupon-library", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    build = commands.add_parser("build", help="register source directories and build the library cases")
    run_gmsh_only_matrix.add_arguments(build)
    build.add_argument("--register", action="append", default=[], type=parse_registration, metavar="CASE_ID=SOURCE_DIR",
                       help="register this source directory before building (repeatable)")
    build.add_argument("--footprint", choices=register_case.FOOTPRINT_DECLARATIONS,
                       help="etch footprint declaration of every registered directory (mandatory with --register)")
    build.add_argument("--inventory-status", choices=register_case.INVENTORY_STATUSES,
                       help="inventory status of every registered directory (mandatory with --register)")
    build.add_argument("--mesh-recipe", help="repository path of the frozen mesh recipe for registered "
                                             "directories without their own mesh-recipe.json")
    build.add_argument("--provenance", help="text appended to the Provenance of every registered case")
    build.add_argument("--no-thin", action="store_true",
                       help="register the fabricated cases only; by default every registered directory (--register and "
                            "--device) also registers and builds its thin counterpart <case>-thin (decision 66: the "
                            "device correction needs both responses of every model)")
    build.add_argument("--work", type=Path, help="parent of the registration work directories")
    build.add_argument("--register-jobs", type=int, default=device_coupons.DEFAULT_REGISTER_JOBS,
                       help="device coupons whose labels-only probe / contract derivation run at once (default "
                            f"{device_coupons.DEFAULT_REGISTER_JOBS}; the manifest append stays serial; decision 62(2))")
    build.add_argument("--device", type=Path, help="a device's Palace config (its ResponseCorrection Library is the "
                                                  "process seed): discovery -> source directories -> registration")
    build.add_argument("--palace", type=Path, help="Palace executable for the discovery preflights (with --device)")
    build.add_argument("--device-output", type=Path, help="output of the device adapter (default ROOT/device)")
    build.add_argument("--ring-size", type=int, default=device_coupons.DEFAULT_RING_SIZE,
                       help="trace-basis ring size of the device coupons (the planner's default)")
    build.add_argument("--cap-triangulation", choices=device_coupons.CAP_TRIANGULATIONS,
                       default=device_coupons.DEFAULT_CAP_TRIANGULATION,
                       help="matching-box cap triangulation of the device basis: delaunay (default; no needle ears, "
                            "decisions 54b / 57) or ear-clipping (the gallery producer's)")
    qualify = commands.add_parser("qualify", help="physics qualification of the built coupons against their references")
    qualify_library.add_arguments(qualify)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "qualify":
        return qualify_library.run_from_args(args)
    if args.command != "build":
        parser.error(f"unknown command {args.command}")
    if args.register and (args.footprint is None or args.inventory_status is None):
        parser.error("--register requires --footprint and --inventory-status (fail closed without a footprint "
                     "declaration)")
    if args.device is not None and args.palace is None:
        parser.error("--device requires --palace (the discovery preflights)")
    registered = []
    extra = None
    if args.device is not None:
        if args.root is None:
            parser.error("--device requires --root (the device adapter writes under ROOT/device)")
        device_output = args.device_output or (args.root / "device")
        try:
            device_record = device_coupons.prepare_device_sources(
                args.device, palace=args.palace, output=device_output, manifest_path=args.manifest, ring_size=args.ring_size,
                cap_triangulation=args.cap_triangulation, python=args.python)
            device_coupons.register_device_sources(device_record, manifest_path=args.manifest, mesh_recipe=args.mesh_recipe,
                                                   work=(args.work or device_output / "register"), python=args.python,
                                                   julia=args.julia, jobs=args.register_jobs, thin=not args.no_thin)
        except device_coupons.DeviceAdapterError as error:
            print(f"DEVICE_ADAPTER_FAILED: {error}", file=sys.stderr)
            return 1
        for coupon in device_record["Coupons"]:
            status = (coupon["Registration"] or {}).get("Status")
            if status in (register_case.STATUS_REGISTERED, register_case.STATUS_REUSED):
                registered.append(coupon["Case"])
            thin_status = (coupon.get("ThinRegistration") or {}).get("Status")
            if thin_status in (register_case.STATUS_REGISTERED, register_case.STATUS_REUSED):
                registered.append(coupon["ThinCase"])
        extra = {"Device": {key: device_record[key] for key in ("Device", "ProcessLibrary", "Discovery", "Plan", "TraceBasis",
                                                                "OutOfScope", "MeshRecipe", "Output")},
                 "DeviceCoupons": device_record["Coupons"]}
        if args.case is None:
            args.case = []
    for case_id, directory in args.register:
        pairs = [(case_id, "fabricated", None)]
        if not args.no_thin:
            pairs.append((register_case.thin_case_id(case_id), "thin", case_id))
        for registered_id, kind, fabricated_case in pairs:
            try:
                record = register_case.register(
                    registered_id, directory, footprint=args.footprint, inventory_status=args.inventory_status,
                    manifest_path=args.manifest, mesh_recipe=args.mesh_recipe, provenance=args.provenance,
                    work=(args.work / registered_id) if args.work is not None else None, python=args.python, julia=args.julia,
                    kind=kind, fabricated_case=fabricated_case)
            except register_case.RegistrationError as error:
                print(f"REGISTRATION_FAILED {registered_id}: {error}", file=sys.stderr)
                return 1
            print(f"{record['Status'].upper()} {registered_id}: {record['Message']}", flush=True)
            registered.append(registered_id)
    if args.case is not None:
        args.case = list(dict.fromkeys(args.case + registered))
    if args.device is not None and not args.case:
        print("DEVICE_ADAPTER_FAILED: no device coupon registered", file=sys.stderr)
        return 1
    return run_gmsh_only_matrix.run_build(args, extra=extra)


if __name__ == "__main__":
    sys.exit(main())
