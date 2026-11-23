#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Prepare the conductor-tagged thin transmon for consistent electrostatic correction.

The default configuration solves the input mesh once (AMR ``MaxIts`` 0). ``--amr-max-its N``
enables N cycles of field-driven AMR with the Refinement block of the recorded prism-library
device runs (transmon_final_library_device_20260830, 2026-09-06..08: Tol 0.01, UpdateFraction
0.7, MaximumImbalance 1.15, MaxNCLevels 8, Nonconformal, SaveAdaptIterations); those runs also
used ``--trace-coupling SurfaceMortar --mortar-oversampling 2``. The geometry-driven
``EdgeRefinement`` tube is an explicit opt-in (``--edge-refinement RADIUS ELEMENTS_PER_RADIUS``):
tubes of Radius 0.2 and ElementsPerRadius 3 (the blocks transmon_surface_amr.json carried
until P4 S4-a) grow the device mesh by ~3.5x per pass over 35 mm of metal edge (138k -> 18.3M
elements after four of about nine passes, measured 2026-09-23), so they are never applied
implicitly.

FORMS (the v1 cut; developer_form.py): by default the response-correction block is the V1 FORM,
post-processing only (no self-consistent solve: ``--correction-mode SelfConsistent`` / ``Both``,
SolveTol and TranslationalDomainCorrection are parked; no mortar: ``--trace-coupling
SurfaceMortar`` and ``--mortar-oversampling`` are parked); ``--developer`` writes the RECORDED
DEVELOPER FORM (the recorded settings: ``--correction-mode`` Both by default, SolveTol 1e-6,
TranslationalDomainCorrection FixedTrace, MortarOversampling 2, with "Developer": true). A parked flag without
``--developer`` is refused. Values from the v1 form are not comparable with records made before
this cut until the benchmark's collocated twins are read.
"""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
IDENTIFICATION = ROOT.parent / "surface_response_identification"
if str(IDENTIFICATION) not in sys.path:
    sys.path.insert(0, str(IDENTIFICATION))
import developer_form  # noqa: E402

# The recorded settings of the response-correction block (the recorded prism-library device runs):
# the developer form writes them, the v1 form less the parked options.
RECORDED_CORRECTION_MODE = "Both"
RECORDED_TRACE_COUPLING = "Collocated"
RECORDED_MORTAR_OVERSAMPLING = 2

# Refinement block of the recorded prism-library device runs (config-p5-exact-amr*.json).
PRISM_ERA_REFINEMENT = {
    "Tol": 0.01,
    "UpdateFraction": 0.7,
    "MaximumImbalance": 1.15,
    "MaxNCLevels": 8,
    "Nonconformal": True,
    "SaveAdaptIterations": True,
    "SaveAdaptMesh": False,
}


def refinement_block(args):
    """Model.Refinement for the produced config: a single solve unless --amr-max-its > 0."""
    if args.amr_max_its == 0:
        return {"MaxIts": 0, "UniformLevels": 0, "SerialUniformLevels": 0}
    block = dict(PRISM_ERA_REFINEMENT)
    block["MaxIts"] = args.amr_max_its
    block["Tol"] = args.amr_tol
    block["SaveAdaptMesh"] = args.save_adapt_mesh
    if args.amr_max_size > 0:
        block["MaxSize"] = args.amr_max_size
        # The predictive MaxSize stop, stated: off by default in the integrated Palace (USER 1636 (4)).
        block["MaxSizePredicted"] = True
    block["UniformLevels"] = 0
    block["SerialUniformLevels"] = 0
    return block


def edge_distances(matching_radius, edge_refinement):
    """EdgeDistances ending at the library matching radius; the tube radius must be one of them."""
    if edge_refinement is None or edge_refinement[0] == matching_radius:
        return [matching_radius]
    return [edge_refinement[0], matching_radius]


def parse_edge_refinement(parser, values):
    if values is None:
        return None
    radius, elements_per_radius = values
    if radius <= 0.0 or radius > 2.0:
        parser.error("edge-refinement radius must be positive and at most the matching radius 2.0")
    if elements_per_radius <= 0 or elements_per_radius != int(elements_per_radius):
        parser.error("edge-refinement elements per radius must be a positive integer")
    return radius, int(elements_per_radius)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--postpro", type=Path, required=True)
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--substrate-permittivity", type=float, default=11.45)
    parser.add_argument(
        "--source", type=Path, default=ROOT / "transmon_surface_coarse.json"
    )
    parser.add_argument(
        "--correction-mode",
        choices=("PostprocessOnly", "SelfConsistent", "Both"),
        help=f"the response-correction mode (the v1 form: PostprocessOnly; the recorded form: "
             f"{RECORDED_CORRECTION_MODE}; SelfConsistent and Both need --developer)",
    )
    parser.add_argument(
        "--trace-coupling",
        choices=("Collocated", "SurfaceMortar"),
        help=f"default {RECORDED_TRACE_COUPLING}; SurfaceMortar needs --developer",
    )
    parser.add_argument("--mortar-oversampling", type=int, choices=range(1, 9),
                        help=f"the mortar quadrature (default {RECORDED_MORTAR_OVERSAMPLING}); needs --developer")
    parser.add_argument(
        "--amr-max-its", type=int, default=0,
        help="AMR cycles with the recorded prism-era Refinement block (default 0: one solve)",
    )
    parser.add_argument("--amr-tol", type=float, default=PRISM_ERA_REFINEMENT["Tol"])
    parser.add_argument(
        "--amr-max-size", type=int, default=0,
        help="stop refining once the solved global unknowns exceed this count (0: no limit)",
    )
    parser.add_argument(
        "--save-adapt-mesh", action="store_true",
        help="write every adapted mesh (<mesh stem>.mesh in the postprocessing directory; .meshgz on a zlib build)",
    )
    parser.add_argument(
        "--edge-refinement", nargs=2, type=float, metavar=("RADIUS_UM", "ELEMENTS_PER_RADIUS"),
        help="opt-in geometry-driven tube refinement of every target interface (see the module docstring)",
    )
    developer_form.add_developer_argument(parser)
    args = parser.parse_args()
    # A parked flag needs --developer (the flags alone are judged, in a block with a Library).
    flags = {"Library": "", "CorrectionMode": args.correction_mode, "TraceCoupling": args.trace_coupling,
             "MortarOversampling": args.mortar_oversampling}
    parked = developer_form.parked_block_settings({k: v for k, v in flags.items() if v is not None})
    if parked and not args.developer:
        parser.error(f"{parked} select a parked developer option: add --developer")
    if args.order <= 0 or args.substrate_permittivity <= 0.0:
        parser.error("order and substrate permittivity must be positive")
    if args.amr_max_its < 0 or args.amr_tol <= 0.0 or args.amr_max_size < 0:
        parser.error("AMR iterations and size must be non-negative and the tolerance positive")
    if args.amr_max_its == 0 and (args.save_adapt_mesh or args.amr_max_size > 0):
        parser.error("--save-adapt-mesh and --amr-max-size require --amr-max-its")
    edge_refinement = parse_edge_refinement(parser, args.edge_refinement)

    config = json.loads(args.source.expanduser().resolve().read_text())
    config["Problem"] = {
        "Type": "Electrostatic",
        "Verbose": 1,
        "Output": str(args.postpro.expanduser().resolve()),
        "OutputFormats": {"Paraview": False, "GridFunction": False},
    }
    config["Model"]["Mesh"] = str(args.mesh.expanduser().resolve())
    config["Model"]["Refinement"] = refinement_block(args)
    config["Domains"]["Materials"] = [
        {"Attributes": [1], "Permittivity": args.substrate_permittivity},
        {"Attributes": [2], "Permittivity": 1.0},
    ]
    postprocessing = config["Boundaries"]["Postprocessing"]
    for dielectric in postprocessing["Dielectric"]:
        dielectric["EdgeDistances"] = edge_distances(2.0, edge_refinement)
        dielectric.pop("EdgeRefinement", None)
        if edge_refinement is not None:
            dielectric["EdgeRefinement"] = {
                "Radius": edge_refinement[0],
                "ElementsPerRadius": edge_refinement[1],
                "OuterRadiusFactor": 2.0,
                "CoreIndicatorWeight": 0.0,
            }
        dielectric["LocalizeEdgeEnergy"] = False
        dielectric["SaveLocalEdgeEnergy"] = False
        if dielectric["Type"] in ("MS", "MA"):
            dielectric["Attributes"] = [5, 6, 7, 9]
        if dielectric["Type"] == "MS":
            dielectric["Permittivity"] = args.substrate_permittivity
    config["Boundaries"] = {
        # Port patches 6 and 7 replace part of the metal sheet and must remain in the
        # grounded conductor union. Attribute 9 is the separately tagged transmon island.
        "Ground": {"Attributes": [5, 6, 7]},
        "Terminal": [{"Index": 1, "Attributes": [9]}],
        "Postprocessing": postprocessing,
    }
    config["Solver"] = {
        "Order": args.order,
        "Electrostatic": {
            "Save": 1,
            "ResponseCorrection": (developer_form.developer_block if args.developer else developer_form.v1_block)({
                "Library": str(args.library.expanduser().resolve()),
                "TargetInterfaces": [1, 2, 3],
                "UnmatchedPolicy": "Error",
                "CorrectionMode": args.correction_mode or RECORDED_CORRECTION_MODE,
                "TranslationalDomainCorrection": "FixedTrace",
                "TraceCoupling": args.trace_coupling or RECORDED_TRACE_COUPLING,
                "MortarOversampling": args.mortar_oversampling or RECORDED_MORTAR_OVERSAMPLING,
                "SolveTol": 1.0e-6,
            }),
        },
        "Linear": {
            "Type": "BoomerAMG",
            "KSPType": "CG",
            "Tol": 1.0e-10,
            "MaxIts": 1000,
            "EstimatorTol": 1.0e-1,
            "EstimatorMG": True,
        },
    }
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(config, indent=2) + "\n")
    print(f"prepare_surface_response_electrostatic.py: the {developer_form.form_name(args.developer)}", file=sys.stderr)
    print(output)


if __name__ == "__main__":
    main()
