#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""coupon-library: the consolidated command line of the coupon mesh library
(supervisor decision 48): `build` (the mesh path) and `qualify` (the physics path).

  coupon_library.py build [--register CASE_ID=SOURCE_DIR ... --footprint {bound,producer-default}
                           --inventory-status STATUS [--mesh-recipe PATH] [--provenance TEXT]]
                          [--case ID ...] [--jobs N] [--root DIR] [--output PATH] [--manifest PATH]
  coupon_library.py qualify --build-record library-build.json --reference <campaign dir or none>
                            --remote HOST:ROOT --orders p4 --controls p3,p5 --max-jobs N
                            --frozen-binary-sha256 HEX [--case ID ...] [--root DIR] [--dry-run]

`build` registers the given source directories as manifest cases (register_case.py:
source SHA256s, the automated two-pass contract derivation, idempotent by content; the
footprint declaration is mandatory and applies to every directory registered by the
call), then builds the selected cases (default: every case of the manifest, the newly
registered ones included) as a job pool (run_gmsh_only_matrix.py: headroom gate ->
run_gmsh_only_case.py DAG -> audits -> verify_canonical_case_entries.py, fail closed per
case) and writes library-build.json.  A registration that fails closed stops the build
before any mesh is made (exit 1) and leaves its record under the case's work directory.

`qualify` runs the physics of every passed coupon of a library-build.json against its
graded_v2 reference (qualify/qualify_library.py: configs from the reference config on
the hash-verified identity mesh, estimate gate, plan with pinned digests, submission
under the user job cap, read-only monitoring, fetch / digest verification / matrix
validation / recorded archive deletion, the frozen machine-readable class gates) and
writes library-qualification.json, qualification-gates.json and process-library.json;
--dry-run writes the plans / configs / estimates / gates without contacting anything.
"""
import argparse
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
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
    build.add_argument("--work", type=Path, help="parent of the registration work directories")
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
    registered = []
    for case_id, directory in args.register:
        try:
            record = register_case.register(
                case_id, directory, footprint=args.footprint, inventory_status=args.inventory_status,
                manifest_path=args.manifest, mesh_recipe=args.mesh_recipe, provenance=args.provenance,
                work=(args.work / case_id) if args.work is not None else None, python=args.python, julia=args.julia)
        except register_case.RegistrationError as error:
            print(f"REGISTRATION_FAILED {case_id}: {error}", file=sys.stderr)
            return 1
        print(f"{record['Status'].upper()} {case_id}: {record['Message']}", flush=True)
        registered.append(case_id)
    if args.case is not None:
        args.case = list(dict.fromkeys(args.case + registered))
    return run_gmsh_only_matrix.run_build(args)


if __name__ == "__main__":
    sys.exit(main())
