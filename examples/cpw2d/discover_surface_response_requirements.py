#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Discover a library-independent surface-response requirement closure."""

import argparse
import copy
import hashlib
import json
import math
import shlex
import subprocess
import sys
from pathlib import Path

import numpy as np

import prepare_surface_response_coupons as planner

IDENTIFICATION_TOOLS = Path(__file__).resolve().parent.parent / "surface_response_identification"
if str(IDENTIFICATION_TOOLS) not in sys.path:
    sys.path.insert(0, str(IDENTIFICATION_TOOLS))
import signature_library  # noqa: E402

SPATIAL_COUPON_TOOLS = Path(__file__).resolve().parent.parent / "cpw3d_surface" / "spatial_coupon"
if str(SPATIAL_COUPON_TOOLS) not in sys.path:
    sys.path.insert(0, str(SPATIAL_COUPON_TOOLS))
import cluster_signature_geometry  # noqa: E402
import device_coupons  # noqa: E402

# The binary's own declaration that a requirement has NO placement in this version
# (surfaceresponseoperator.cpp, the mirror-formed status loop: a mirror-formed CurvedEdge
# configuration is "the 2D curvature family's consumer: no mirror placement in this version",
# decision 557). Such a requirement stays Missing whatever model the library carries, so a
# placeholder can never close it: the closure records it OutOfScope
# UnplaceableInThisVersion (device_coupons.UNPLACEABLE_METHOD) by name instead of stalling
# (decisions 664 / 669 (3) (c): the R4 39e4 attempt-2 stall on 4 mirror-formed CurvedEdge keys).
UNPLACEABLE_DECLARATION = "no mirror placement in this version"

# The two discovery-side RECORD classes of the readiness hardening-2 lane (USER decision
# 1000 (b), decisions 1012 / 1013 F-2 / F-4a, 1018): a Missing requirement that can get no
# placeholder is recorded by name and the closure continues (RAW-BY-FALLBACK: no coupon is
# built for it, the library is built for the rest) instead of the discovery dying.
#  - UNBOXABLE_METHOD: the key is an UnboxableFeature signature ("Unboxable": true, decision
#    282: no box satisfies the face rules), which no builder makes and the library reader
#    refuses to serve (surfaceresponseoperator.cpp ReadProcessLibrary); the census F-2.
#  - UNBUILDABLE_METHOD: the generator's own geometry refuses the key by name
#    (cluster_signature_geometry.SignatureGeometryError: an arc joint in neither built range
#    of the coupon mesher, decision 673 JointTurnGap; the census F-4a / F-4b).
UNBOXABLE_METHOD = "Unboxable"
UNBUILDABLE_METHOD = "UnbuildableJoint"
UNBOXABLE_RULE = ("decision 282 / 1018: an UnboxableFeature signature (\"Unboxable\": true) is a Missing "
                  "placeholder no builder makes and no library may serve; the discovery gives it no "
                  "placeholder, lists it RAW-BY-FALLBACK and continues")
UNBUILDABLE_RULE = ("decision 673 / 1018: the coupon generator refuses the key by name (SignatureGeometryError); "
                    "the discovery records the refusal, gives it no placeholder, lists it RAW-BY-FALLBACK and "
                    "continues (the mesher's joint classes: OI-MESHER-JOINT-CLASS-GAP / OI-MESHER-BUILT-JOINT-TILT)")


def unboxable_requirement(requirement):
    """True for a version-2 requirement keyed by an UnboxableFeature signature."""
    signature = requirement.get("Signature")
    return isinstance(signature, dict) and bool(signature.get("Unboxable", False))


def requirement_record(requirement, **extra):
    record = {
        "Hash": requirement["Hash"],
        "Topology": requirement["Topology"],
        "Count": requirement.get("Count"),
        "Instances": requirement.get("Instances"),
        "TotalEdgeLength": requirement.get("TotalEdgeLength"),
    }
    record.update(extra)
    return record


def unplaceable_declaration(requirement):
    """The binary's Note declaring the requirement unplaceable in this version (None when the
    requirement's Notes carry no such declaration)."""
    for note in requirement.get("Notes") or []:
        if isinstance(note, str) and UNPLACEABLE_DECLARATION in note:
            return note
    return None


def load_json(path):
    return json.loads(path.read_text())


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def response_section(config):
    solver = config.setdefault("Solver", {})
    if config.get("Problem", {}).get("Type") == "Electrostatic":
        electrostatic = solver.setdefault("Electrostatic", {})
        response = electrostatic.get("ResponseCorrection")
        if response is None:
            raise ValueError(
                "Electrostatic closure discovery requires "
                "Solver.Electrostatic.ResponseCorrection"
            )
        return response
    response = solver.get("SurfaceResponseCorrection")
    if response is None:
        raise ValueError(
            "Maxwell closure discovery requires Solver.SurfaceResponseCorrection"
        )
    return response


def resolve_library(config_path, response):
    path = Path(response["Library"]).expanduser()
    return path.resolve() if path.is_absolute() else (config_path.parent / path).resolve()


def unique_interfaces(requirement):
    result = []
    seen = set()
    for entry in requirement.get("Interfaces", []):
        key = (int(entry.get("Slot", 0)), entry["Type"])
        if key in seen:
            continue
        seen.add(key)
        result.append(
            {"Slot": key[0], "Type": key[1], "Coupon": len(result) + 1}
        )
    result.sort(key=lambda item: (item["Slot"], item["Type"]))
    for index, entry in enumerate(result, start=1):
        entry["Coupon"] = index
    return result


def conductor_references(edges, separation=1.0):
    count = max((int(edge.get("Conductor", 1)) for edge in edges), default=1)
    references = []
    for conductor in range(1, count + 1):
        points = [
            edge.get("Point")
            for edge in edges
            if int(edge.get("Conductor", 1)) == conductor and "Point" in edge
        ]
        if points:
            reference = [float(value) for value in points[0]]
        else:
            reference = [float(conductor - 1) * separation, 0.0, 0.0]
        # Preflight requires distinct references but does not read response matrices.
        if reference in references:
            reference[0] += conductor * max(separation, 1.0)
        references.append(reference)
    return references


def spatial_support_points(geometry, matching_radius, fabrication):
    """Reproduce the production spatial-coupon matching box without meshing it."""
    edges = geometry.get("Edges", [])
    if not edges:
        raise ValueError("SpatialEdgeCluster has no complete edge geometry")
    normal = np.asarray(edges[0]["ProcessNormal"], dtype=float)
    normal /= np.linalg.norm(normal)
    gap = np.asarray(edges[0]["GapDirection"], dtype=float)
    gap /= np.linalg.norm(gap)
    axis_y = np.cross(normal, gap)
    axis_y /= np.linalg.norm(axis_y)
    frame = np.vstack((gap, axis_y, normal))
    radius = float(matching_radius)
    points = []
    local_edges = []
    for entry in edges:
        point = frame @ np.asarray(entry["Point"], dtype=float)
        local_gap = frame @ np.asarray(entry["GapDirection"], dtype=float)
        local_normal = frame @ np.asarray(entry["ProcessNormal"], dtype=float)
        tangent = np.cross(local_gap, local_normal)
        begin, end = (float(value) for value in entry["Interval"])
        tolerance = 1.0e-10 * radius
        if begin <= -radius + tolerance:
            begin -= 2.0 * radius
        if end >= radius - tolerance:
            end += 2.0 * radius
        for coordinate in (begin, end):
            boundary = point + coordinate * tangent
            for side in (-1.0, 1.0):
                points.append(boundary + side * radius * local_gap)
        local_edges.append(point)
    points = np.asarray(points)
    lower = np.min(points, axis=0) - radius
    upper = np.max(points, axis=0) + radius
    metal_thickness = float(fabrication.get("MetalThickness", 0.0))
    overetch_depth = float(fabrication.get("OveretchDepth", 0.0))
    lower[2] = min(
        lower[2], min(point[2] for point in local_edges) - radius - overetch_depth
    )
    upper[2] = max(
        upper[2], max(point[2] for point in local_edges) + radius + metal_thickness
    )
    local_support = np.asarray(
        [
            (x, y, z)
            for x in (lower[0], upper[0])
            for y in (lower[1], upper[1])
            for z in (lower[2], upper[2])
        ]
    )
    support = local_support @ frame
    tolerance = max(
        1.0e-10 * radius, 64.0 * np.finfo(float).eps
    )
    exponent = math.floor(math.log10(tolerance))
    step = 10.0**exponent
    decimals = max(0, -exponent)
    support = np.asarray(
        [[round(float(value), decimals) for value in point] for point in support]
    )
    support[np.abs(support) < 0.5 * step] = 0.0
    return support.tolist()


def mirror_formed_placeholder_edges(requirement, matching_radius):
    """The Signature-keyed Edges and the ENTRY `MirrorFormed` record of the discovery placeholder
    of a SpatialEdgeCluster requirement carrying a MirrorFormedContract (decision 599): the
    operator places a mirror-formed configuration on its real half only through a model whose
    Edges carry the contract's per-Edge Weight (surfaceresponseoperator.cpp fails closed on a
    signature-only model, decision 557 (4)), so the placeholder is built as the generator's
    stamp builds the real coupon's entry (device_coupons.stamp_signature_model): the contract
    validated field by field (fail closed by name), the Edges from the signature's chorded
    portions (cluster_signature_geometry.model_edges), Weight 1.0 real / 0.0 image and the
    `MirrorFormed` record (device_coupons.mirror_formed_entry). Returns (edges, record)."""
    contract = device_coupons.validate_mirror_formed_contract(requirement)
    edges = cluster_signature_geometry.model_edges(requirement, matching_radius)
    portions = cluster_signature_geometry.portions_from_signature(
        requirement["Signature"], matching_radius
    )
    if len(portions) != len(edges):
        raise ValueError(
            f"requirement {requirement['Hash'][:12]}: {len(edges)} model edges for "
            f"{len(portions)} chorded portions"
        )
    record, weights = device_coupons.mirror_formed_entry(
        contract, [portion["Portion"] for portion in portions], matching_radius
    )
    for edge, weight in zip(edges, weights):
        edge["Weight"] = weight
    return edges, record


def placeholder_model(requirement, matching_radius, fabrication=None):
    if "Signature" in requirement:
        # A version-2 record (SURFACE-RESPONSE-IDENTIFICATION.md (d)): the model is keyed by
        # the record's canonical Signature (its representative over the tolerance group), which
        # the matcher accepts for every topology; no version-1 site geometry exists or is needed.
        digest = requirement["Hash"][:16]
        model = signature_library.signature_model(
            requirement["Topology"], requirement["Signature"], float(matching_radius),
            f"__preflight_placeholder_{digest}", "__preflight_dummy")
        model["BoundaryLawQualification"] = {
            "Version": 1,
            "Status": "Unqualified",
            "Calibration": "GeometryDiscoveryOnly",
            "FrequencyUniversal": False,
        }
        if (
            requirement["Topology"] == "SpatialEdgeCluster"
            and "MirrorFormedContract" in requirement
        ):
            # A mirror-formed cluster requirement (impl-B5 CONTRACT.md; routed on the
            # CONTRACT, never on the MirrorFormed flag, decision 562): the operator applies a
            # model to the configuration's real half only when the model carries Signature-
            # keyed Edges with the contract's weights, so the placeholder carries them (and
            # the interface mappings its Edges' InterfaceSlots bind) like the real coupon's
            # entry; a signature-only placeholder would leave the key Missing on every pass.
            edges, record = mirror_formed_placeholder_edges(requirement, float(matching_radius))
            model["Edges"] = edges
            model["MirrorFormed"] = record
            # The reader's rule (ReadProcessLibrary: one conductor reference per canonical
            # conductor label of the model's Edges and PlanViewBoundary components; decision
            # 1018 Q2, the census F-3): the placeholder carries no PlanViewBoundary, so its
            # references follow its Edges, not the signature's Context (a context-only
            # conductor is the generator's business: it writes the PlanViewBoundary the
            # reader counts, decision 649 (1)).
            edge_conductors = max(int(edge["Conductor"]) for edge in edges)
            model.pop("ConductorReferences", None)
            if edge_conductors > 1:
                model["ConductorReferences"] = [
                    [0.0, 0.0, -float(matching_radius) * (k + 1)] for k in range(edge_conductors)
                ]
            model["FabricatedSurfaceMatrix"] = "__preflight_dummy_fabricated_surface.csv"
            model["ThinSurfaceMatrix"] = "__preflight_dummy_thin_surface.csv"
            model["Interfaces"] = unique_interfaces(requirement)
        elif requirement["Topology"] != "SpatialEdgeCluster":
            # The library reader binds a spatial model's interface mappings (and with them
            # its surface matrices) to the InterfaceSlot of its stored Edges; a signature-only
            # placeholder has none (the matching key is the Signature, whose portions carry
            # the interface types), so it declares neither.
            model["FabricatedSurfaceMatrix"] = "__preflight_dummy_fabricated_surface.csv"
            model["ThinSurfaceMatrix"] = "__preflight_dummy_thin_surface.csv"
            model["Interfaces"] = unique_interfaces(requirement)
        return digest, model
    signature = planner.coupon_signature(requirement)
    digest = hashlib.sha256(
        json.dumps(signature, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()[:16]
    topology = requirement["Topology"]
    geometry = requirement.get("Geometry", {})
    boundary = requirement.get("BoundaryCondition", {"Type": "PEC"})
    model = {
        "Name": f"__preflight_placeholder_{digest}",
        "Topology": topology,
        "FabricatedMatrix": "__preflight_dummy_fabricated.csv",
        "ThinMatrix": "__preflight_dummy_thin.csv",
        "FabricatedSurfaceMatrix": "__preflight_dummy_fabricated_surface.csv",
        "ThinSurfaceMatrix": "__preflight_dummy_thin_surface.csv",
        "BasisPoints": "__preflight_dummy_basis.csv",
        "Interfaces": unique_interfaces(requirement),
        "BoundaryLawQualification": {
            "Version": 1,
            "Status": "Unqualified",
            "Calibration": "GeometryDiscoveryOnly",
            "FrequencyUniversal": False,
        },
    }

    if topology == "IsolatedEdge":
        model["Reference"] = [0.0, 0.0, 0.0]
        model["BoundaryCondition"] = boundary
    elif topology in (
        "SameConductorGap",
        "SameConductorStrip",
        "DifferentConductorGap",
    ):
        separation = float(geometry["Separation"])
        model["Separation"] = separation
        model["SeparationTolerance"] = max(1.0e-9 * matching_radius, 1.0e-12)
        model["BoundaryCondition"] = boundary
        if topology == "DifferentConductorGap":
            model["ConductorReferences"] = [
                [-0.5 * separation, 0.0, 0.0],
                [0.5 * separation, 0.0, 0.0],
            ]
        else:
            model["Reference"] = [0.0, 0.0, 0.0]
    elif topology == "ParallelEdgeCluster":
        edges = copy.deepcopy(geometry["Edges"])
        for edge in edges:
            if isinstance(edge.get("Offset"), list):
                edge["Offset"] = float(edge["Offset"][0])
            if isinstance(edge.get("GapDirection"), list):
                edge["GapDirection"] = int(edge["GapDirection"][0])
        model["Edges"] = edges
        model["EdgeOffsetTolerance"] = max(1.0e-9 * matching_radius, 1.0e-12)
        model["BoundaryCondition"] = boundary
        references = conductor_references(edges, matching_radius)
        if len(references) == 1:
            model["Reference"] = references[0]
        else:
            model["ConductorReferences"] = references
    elif topology in ("ConvexCorner", "ConcaveCorner"):
        model["Angle"] = float(geometry["AngleDegrees"])
        model["AngleTolerance"] = 1.0e-6
        model["CornerRadius"] = float(geometry.get("CornerRadius", 0.0))
        model["CornerRadiusTolerance"] = max(
            1.0e-9 * matching_radius, 1.0e-12
        )
        model["Reference"] = [0.0, 0.0, 0.0]
        model["BoundaryCondition"] = boundary
    elif topology == "SpatialEdgeCluster":
        edges = copy.deepcopy(geometry["Edges"])
        model["Edges"] = edges
        model["EdgePositionTolerance"] = max(
            1.0e-4 * matching_radius, 1.0e-12
        )
        model["EdgeAngleTolerance"] = 1.0e-3
        references = conductor_references(edges, matching_radius)
        if len(references) == 1:
            model["Reference"] = references[0]
        else:
            model["ConductorReferences"] = references
        # Virtual closure must use the same exact mask and matching support as the
        # production model. Omitting either can make ownership change after the real
        # coupon is inserted and expose a new overlap only after expensive generation.
        for key in ("PlanViewBoundary", "MaskRegularization"):
            if key in geometry:
                model[key] = copy.deepcopy(geometry[key])
        model["SupportPoints"] = spatial_support_points(
            geometry, matching_radius, fabrication or {}
        )
    elif topology in ("Endpoint", "Junction"):
        model["Reference"] = [0.0, 0.0, 0.0]
        model["BoundaryCondition"] = boundary
        if topology == "Junction":
            model["ArmAngles"] = geometry["ArmAnglesDegrees"]
            model["ArmAngleTolerance"] = 1.0e-6
        for key in ("PlanViewBoundary", "MaskRegularization"):
            if key in geometry:
                model[key] = copy.deepcopy(geometry[key])
    else:
        raise ValueError(f"Unsupported placeholder topology {topology}")
    return digest, model


def run_preflight(palace, config_path, log_path):
    command = [str(palace), "--surface-response-preflight", str(config_path)]
    print("+ " + shlex.join(command), flush=True)
    with log_path.open("w") as stream:
        subprocess.run(command, check=True, stdout=stream, stderr=subprocess.STDOUT)


def restore_source_status(manifest, source_library, placeholder_requirements):
    result = copy.deepcopy(manifest)
    result["Library"]["Name"] = source_library.get(
        "Name", Path(result["Library"]["Path"]).stem
    )
    result["Library"]["Path"] = str(source_library["__SourcePath"])
    counts = {"Exact": 0, "Interpolated": 0, "Missing": 0}
    lengths = {"Exact": 0.0, "Interpolated": 0.0, "Missing": 0.0}
    for requirement in result["Requirements"]:
        selected = requirement.get("SelectedModels", [])
        placeholders = [
            model.get("Name")
            for model in selected
            if model.get("Name") in placeholder_requirements
        ]
        if placeholders:
            # Use the requirement which created the selected virtual model. The production
            # matcher may describe the same match in model-local coordinates on later
            # passes; feeding that transformed description back to the coupon generator
            # would create a different signature than the virtual model that established
            # closure.
            source = placeholder_requirements[placeholders[0]]
            for key in ("Topology", "Geometry", "BoundaryCondition", "Interfaces"):
                if key in source:
                    requirement[key] = copy.deepcopy(source[key])
                else:
                    requirement.pop(key, None)
            requirement["Status"] = "Missing"
            requirement["Reason"] = (
                "Missing from source library after exhaustive geometry discovery"
            )
            requirement.pop("SelectedModels", None)
            requirement.pop("NormalizedLibraryDistance", None)
        status = requirement["Status"]
        counts[status] += int(requirement["Count"])
        lengths[status] += float(requirement.get("TotalEdgeLength", 0.0))
    result["Summary"] = {"Counts": counts, "TotalEdgeLengths": lengths}
    if "GeometryFallbacks" in manifest.get("Summary", {}):
        # The geometry fallbacks of the final preflight (USER decision 1000 (b)) survive the
        # status restoration: the designer's manifest names its RAW-BY-FALLBACK features.
        result["Summary"]["GeometryFallbacks"] = copy.deepcopy(
            manifest["Summary"]["GeometryFallbacks"]
        )
    result["Complete"] = counts["Missing"] == 0
    return result


def omitted_prefix(requirement, omit_requirements):
    """The `--omit-requirement` prefix a version-2 record's Hash starts with (None: the
    requirement is not omitted; version-1 records carry no Hash and are never omitted)."""
    digest = requirement.get("Hash", "")
    for prefix in omit_requirements:
        if digest.startswith(prefix):
            return prefix
    return None


def discover(config_path, output, palace, max_passes=8, omit_requirements=()):
    """The closure loop: geometry preflights with the source library plus one signature
    placeholder per Missing requirement until nothing new is Missing; writes
    output/surface-response-requirements.json (the final manifest with the placeholders
    restored to Missing) and output/closure-history.json, returns the final manifest.

    `omit_requirements`: Hash prefixes of Missing requirements that get NO placeholder (an
    omitted requirement stays Missing in every pass and in the final manifest, is listed
    under OmittedRequirements of the manifest and the history, and does not count as a
    stalled closure); a prefix that matches no requirement fails closed.

    A requirement the binary DECLARES UNPLACEABLE (its Notes carry UNPLACEABLE_DECLARATION)
    and that stays Missing with its placeholder in the virtual library is recorded under
    UnplaceableRequirements of the manifest and the history (OutOfScope
    device_coupons.UNPLACEABLE_METHOD, the declaration quoted) and does not count as a
    stalled closure either; it stays Missing in the final manifest with that Reason. A key
    whose Notes carry the declaration but that turns Exact once its placeholder exists (a
    real instance of the same key beside the unplaceable one) closes as any other."""
    if max_passes <= 0:
        raise ValueError("max_passes must be positive")
    omit_requirements = [str(prefix) for prefix in omit_requirements]
    if any(not prefix for prefix in omit_requirements):
        raise ValueError("an omitted requirement needs a non-empty Hash prefix")
    config_path = Path(config_path).expanduser().resolve()
    config = load_json(config_path)
    response = response_section(config)
    source_path = resolve_library(config_path, response)
    source_library = load_json(source_path)
    source_library["__SourcePath"] = source_path
    virtual_library = copy.deepcopy(source_library)
    virtual_library.pop("__SourcePath", None)
    # Version 2 is sufficient for virtual multi-conductor references. Preserve version 3
    # when fabrication metadata is present, but do not force older libraries to invent it.
    virtual_library["Version"] = max(2, int(virtual_library.get("Version", 0)))
    virtual_library["Name"] = f"{virtual_library.get('Name', source_path.stem)}-closure"
    virtual_library["ExhaustiveSpatialClosure"] = True

    output = Path(output).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    known = set()
    placeholder_names = set()
    placeholder_requirements = {}
    omitted = {}
    unplaceable = {}
    unboxable = {}
    unbuildable = {}
    history = []
    final_manifest = None
    for pass_index in range(1, max_passes + 1):
        pass_root = output / f"pass-{pass_index:02d}"
        pass_root.mkdir(parents=True, exist_ok=True)
        library_path = pass_root / "process-library.json"
        write_json(library_path, virtual_library)
        pass_config = copy.deepcopy(config)
        pass_config["Problem"]["Output"] = str(pass_root / "postpro")
        pass_response = response_section(pass_config)
        pass_response["Library"] = str(library_path)
        pass_response["UnmatchedPolicy"] = "Warn"
        # Geometry fallbacks (USER decision 1000 (b), supervisor 1001 Q2): the discovery
        # preflights run with the EXPLICIT RecordRaw leaf, so a design with a geometry
        # fallback (a feature the placement refuses, an alias key mismatch, ...) still
        # closes: the fallback is RECORDED (Summary.GeometryFallbacks, listed RAW-BY-FALLBACK
        # below) and the library is built for the rest; a bare preflight (FailClosed default)
        # would exit non-zero instead. One bad feature never kills a whole library build.
        pass_response["GeometryFailurePolicy"] = "RecordRaw"
        pass_config_path = pass_root / "config.json"
        write_json(pass_config_path, pass_config)
        run_preflight(palace, pass_config_path, pass_root / "preflight.log")
        manifest = load_json(pass_root / "postpro" / "surface-response-requirements.json")

        # Discovery deliberately uses wider geometry components to expose every support
        # overlap before coupon generation. Validate the same virtual library with strict
        # production component construction as well; the returned library must satisfy
        # both views, not merely the wider discovery partition.
        production_config = copy.deepcopy(pass_config)
        production_config["Problem"]["Output"] = str(pass_root / "production-postpro")
        production_response = response_section(production_config)
        production_response["UnmatchedPolicy"] = "Error"
        production_response["GeometryFailurePolicy"] = "RecordRaw"
        production_config_path = pass_root / "production-config.json"
        write_json(production_config_path, production_config)
        run_preflight(
            palace,
            production_config_path,
            pass_root / "production-preflight.log",
        )
        production_manifest = load_json(
            pass_root / "production-postpro" / "surface-response-requirements.json"
        )
        final_manifest = production_manifest

        added = []
        stalled = 0
        # The placeholders both preflights of this pass were run with (a placeholder added
        # from the discovery view is not in the production view's library yet).
        placed_before_pass = set(known)
        for current_manifest in (manifest, production_manifest):
            for requirement in current_manifest["Requirements"]:
                if requirement["Status"] != "Missing":
                    continue
                prefix = omitted_prefix(requirement, omit_requirements)
                if prefix is not None:
                    omitted.setdefault(
                        requirement["Hash"],
                        {
                            "Hash": requirement["Hash"],
                            "Prefix": prefix,
                            "Topology": requirement["Topology"],
                            "Count": requirement.get("Count"),
                            "Instances": requirement.get("Instances"),
                            "TotalEdgeLength": requirement.get("TotalEdgeLength"),
                        },
                    )
                    continue
                if unboxable_requirement(requirement):
                    # F-2 (decision 1018 Q3): no placeholder exists for an Unboxable key (the
                    # reader refuses one); recorded, never a stall.
                    unboxable.setdefault(
                        requirement["Hash"],
                        requirement_record(
                            requirement,
                            EdgeCount=requirement.get("Signature", {}).get("EdgeCount"),
                            MirrorFormed=requirement.get("MirrorFormed"),
                            OutOfScope=UNBOXABLE_METHOD,
                            Rule=UNBOXABLE_RULE,
                            Pass=pass_index,
                        ),
                    )
                    continue
                if requirement["Hash"] in unbuildable:
                    continue
                try:
                    digest, model = placeholder_model(
                        requirement,
                        float(current_manifest["Library"]["MatchingRadius"]),
                        virtual_library.get("Fabrication", {}),
                    )
                except cluster_signature_geometry.SignatureGeometryError as error:
                    # F-4a / F-4b (decision 1018 Q3): the generator refuses the key by name;
                    # recorded with the message, no placeholder, never a stall.
                    unbuildable[requirement["Hash"]] = requirement_record(
                        requirement,
                        MirrorFormed=requirement.get("MirrorFormed"),
                        OutOfScope=UNBUILDABLE_METHOD,
                        Error=str(error),
                        Rule=UNBUILDABLE_RULE,
                        Pass=pass_index,
                    )
                    continue
                if digest in known:
                    if digest not in placed_before_pass:
                        continue
                    # Still Missing with its placeholder in the virtual library: the binary
                    # either declares the requirement unplaceable (recorded, not a stall) or
                    # the placeholder failed to close it (a stall).
                    declaration = unplaceable_declaration(requirement)
                    if declaration is None:
                        stalled += 1
                    else:
                        unplaceable.setdefault(
                            requirement["Hash"],
                            {
                                "Hash": requirement["Hash"],
                                "Topology": requirement["Topology"],
                                "Count": requirement.get("Count"),
                                "Instances": requirement.get("Instances"),
                                "TotalEdgeLength": requirement.get("TotalEdgeLength"),
                                "MirrorFormed": requirement.get("MirrorFormed"),
                                "OutOfScope": device_coupons.UNPLACEABLE_METHOD,
                                "Declaration": declaration,
                                "Placeholder": model["Name"],
                                "Pass": pass_index,
                            },
                        )
                    continue
                known.add(digest)
                placeholder_names.add(model["Name"])
                placeholder_requirements[model["Name"]] = copy.deepcopy(requirement)
                virtual_library["Models"].append(model)
                added.append(
                    {"Id": digest, "Topology": requirement["Topology"]}
                )
        history.append(
            {
                "Pass": pass_index,
                "Summary": manifest["Summary"],
                "DiscoverySummary": manifest["Summary"],
                "ProductionSummary": production_manifest["Summary"],
                "AddedPlaceholders": added,
                "GeometryFallbacks": geometry_fallbacks_record(production_manifest),
            }
        )
        if not added:
            # Missing requirements other than the omitted ones (which never get a
            # placeholder) and the unplaceable ones (which the binary declares so with
            # their placeholder present) stall the closure.
            if stalled:
                raise RuntimeError(
                    "Geometry closure stalled with missing requirements; see "
                    f"{pass_root / 'postpro/surface-response-requirements.json'} and "
                    f"{pass_root / 'production-postpro/surface-response-requirements.json'}"
                )
            break
    else:
        raise RuntimeError(f"Geometry closure did not converge in {max_passes} passes")
    unmatched = sorted(
        set(omit_requirements) - {record["Prefix"] for record in omitted.values()}
    )
    if unmatched:
        raise RuntimeError(
            "omitted requirement prefixes matching no Missing requirement of the device: "
            + ", ".join(unmatched)
        )

    source_library["__SourcePath"] = source_path
    final_manifest = restore_source_status(
        final_manifest, source_library, placeholder_requirements
    )
    omitted_records = [omitted[key] for key in sorted(omitted)]
    final_manifest["OmittedRequirements"] = omitted_records
    unplaceable_records = [unplaceable[key] for key in sorted(unplaceable)]
    for requirement in final_manifest["Requirements"]:
        record = unplaceable.get(requirement.get("Hash"))
        if record is not None and requirement["Status"] == "Missing":
            requirement["OutOfScope"] = device_coupons.UNPLACEABLE_METHOD
            requirement["Reason"] = (
                f"{device_coupons.UNPLACEABLE_METHOD}: declared by the binary, no placeholder "
                f"can close it ({record['Declaration']})"
            )
    final_manifest["UnplaceableRequirements"] = unplaceable_records
    unboxable_records = [unboxable[key] for key in sorted(unboxable)]
    unbuildable_records = [unbuildable[key] for key in sorted(unbuildable)]
    for requirement in final_manifest["Requirements"]:
        if requirement["Status"] != "Missing":
            continue
        record = unboxable.get(requirement.get("Hash"))
        if record is not None:
            requirement["OutOfScope"] = UNBOXABLE_METHOD
            requirement["Reason"] = f"{UNBOXABLE_METHOD}: {UNBOXABLE_RULE}"
        record = unbuildable.get(requirement.get("Hash"))
        if record is not None:
            requirement["OutOfScope"] = UNBUILDABLE_METHOD
            requirement["Reason"] = f"{UNBUILDABLE_METHOD}: {record['Error']}"
    if unboxable_records:
        final_manifest["UnboxableRequirements"] = unboxable_records
        final_manifest["Summary"]["UnboxableRequirements"] = {
            "Count": len(unboxable_records),
            "TotalEdgeLength": sum(float(r.get("TotalEdgeLength") or 0.0) for r in unboxable_records),
        }
    if unbuildable_records:
        final_manifest["UnbuildableRequirements"] = unbuildable_records
        final_manifest["Summary"]["UnbuildableRequirements"] = {
            "Count": len(unbuildable_records),
            "TotalEdgeLength": sum(float(r.get("TotalEdgeLength") or 0.0) for r in unbuildable_records),
        }
    manifest_path = output / "surface-response-requirements.json"
    write_json(manifest_path, final_manifest)
    geometry_fallbacks = geometry_fallbacks_record(final_manifest)
    closure_history = {
        "Version": 1,
        "SourceConfig": str(config_path),
        "SourceLibrary": str(source_path),
        "Passes": history,
        "PlaceholderCount": len(placeholder_names),
        "OmittedRequirements": omitted_records,
        "UnplaceableRequirements": unplaceable_records,
        "CompleteAgainstSourceLibrary": final_manifest["Complete"],
    }
    if unboxable_records:
        closure_history["UnboxableRequirements"] = unboxable_records
    if unbuildable_records:
        closure_history["UnbuildableRequirements"] = unbuildable_records
    if unboxable_records or unbuildable_records:
        print(describe_discovery_fallbacks(unboxable_records, unbuildable_records), flush=True)
    if geometry_fallbacks is not None:
        closure_history["GeometryFallbacks"] = geometry_fallbacks
        print(describe_geometry_fallbacks(geometry_fallbacks), flush=True)
    write_json(output / "closure-history.json", closure_history)
    return final_manifest


def describe_discovery_fallbacks(unboxable_records, unbuildable_records):
    """The RAW-BY-FALLBACK listing of the requirements the DISCOVERY itself could give no
    placeholder (decision 1018: Unboxable keys, generator-refused keys): printed with the
    closure beside the binary's GeometryFallbacks listing (never silent, USER 1000 (b))."""
    records = [(UNBOXABLE_METHOD, r) for r in unboxable_records] + [
        (UNBUILDABLE_METHOD, r) for r in unbuildable_records
    ]
    lines = [
        f"DISCOVERY FALLBACKS: {len(records)} Missing requirement(s) got NO placeholder and stay "
        "RAW-BY-FALLBACK (no coupon is built for them; the library is built for the rest)"
    ]
    for method, record in records:
        lines.append(
            "  RAW-BY-FALLBACK {method} {hash} ({topology}; instances {count}; length {length:.6e})".format(
                method=method,
                hash=record["Hash"],
                topology=record.get("Topology"),
                count=record.get("Count"),
                length=float(record.get("TotalEdgeLength") or 0.0),
            )
            + (f": {record['Error']}" if "Error" in record else "")
        )
    return "\n".join(lines)


def geometry_fallbacks_record(manifest):
    """The manifest's GeometryFallbacks record (Identification.Diagnostics, with the Summary
    leaf's RawShareByLength) or None when no converted site fired (the normal case)."""
    record = manifest.get("Identification", {}).get("Diagnostics", {}).get("GeometryFallbacks")
    if record is None:
        return None
    return {
        "Count": record.get("Count"),
        "Classes": record.get("Classes", {}),
        "RawShareByLength": record.get("RawShareByLength", {}),
        "Entries": record.get("Entries", []),
    }


def describe_geometry_fallbacks(record):
    """The RAW-BY-FALLBACK listing printed with the closure (USER decision 1000 (b): never
    silent): every feature the preflight recorded as a geometry fallback, with its class."""
    lines = [
        f"GEOMETRY FALLBACKS (RecordRaw): {record['Count']} recorded; the listed features are "
        "RAW-BY-FALLBACK (no coupon is built for them; the library is built for the rest)"
    ]
    for entry in record["Entries"]:
        lines.append(
            "  {cls} ({what}; features {features}; models {models}; length {length:.6e}): "
            "{reason}".format(
                cls=entry.get("Class"),
                what=entry.get("What"),
                features=entry.get("Features", []),
                models=entry.get("Models", []),
                length=float(entry.get("Length", 0.0)),
                reason=entry.get("Reason", ""),
            )
        )
    for type_name, share in sorted(record.get("RawShareByLength", {}).items()):
        lines.append(
            f"  raw share by length {type_name}: {100.0 * float(share.get('Share', 0.0)):.3f} %"
        )
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--palace", type=Path, required=True)
    parser.add_argument("--max-passes", type=int, default=8)
    parser.add_argument(
        "--omit-requirement",
        action="append",
        default=[],
        metavar="HASH_PREFIX",
        help="give NO placeholder to the Missing requirement(s) whose Hash starts with this "
        "prefix (repeatable): they stay Missing in the final manifest and are listed under "
        "its OmittedRequirements; a prefix matching no requirement fails closed",
    )
    args = parser.parse_args()
    if args.max_passes <= 0:
        parser.error("--max-passes must be positive")
    final_manifest = discover(
        args.config,
        args.output,
        args.palace,
        max_passes=args.max_passes,
        omit_requirements=args.omit_requirement,
    )
    print(args.output.expanduser().resolve() / "surface-response-requirements.json")
    print(json.dumps(final_manifest["Summary"], indent=2))


if __name__ == "__main__":
    main()
