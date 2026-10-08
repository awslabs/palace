#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""`coupon-library build --device`: from a device's Palace config to registered
coupon-library cases (supervisor decisions 48 step 1 and 52, the discovery ->
source-directory adapter).

1. Discovery: examples/cpw2d/discover_surface_response_requirements.py on the device
   config (its ResponseCorrection / SurfaceResponseCorrection Library is the process
   seed: version-3 Fabrication metadata) with the Palace executable - geometry
   preflights only, never a mesh or a solve - gives the requirement closure
   (surface-response-requirements.json: every Missing requirement with its topology,
   geometry, interfaces).
2. Plan: prepare_surface_response_coupons.plan_from_manifest routes every requirement
   to its builder; the SpatialEdgeCluster coupons (Preparation.Method SpatialCoupon)
   are this library's; every other family (corners, straight edges) is recorded
   out of scope with its method, never dropped silently.
3. Source directory per spatial coupon: the planner's canonical plan-view boundary
   and mask regularization (the same rule its execute path applies), then
   generate_spatial_response.py --basis-only with the seed's fabrication - the
   signature files, the trace basis (basis-contract.json, trace-vertices.csv,
   trace-triangles.csv, basis-points.csv, zero-trace.csv, conductor-N.csv) and the
   process-library.json of the model; process.toml from the fabrication; a
   provenance.json naming the device config, the discovery closure, the requirement
   and every tool.  The device basis triangulates the matching-box caps with Delaunay
   flips by default (--cap-triangulation delaunay, supervisor decision 57: the
   ear-clipped caps produced the needle triangles that drove 20-40% of the mesh cost;
   ear-clipping, the gallery producer's, stays an explicit option and the gallery
   references keep it).  The directory is named by the content hash of its bound source
   files (`spatial-<edge count>-edge-<hash12>`): the same device geometry always maps
   to the same case, and register_case.py reuses a case whose digests it already holds.
4. Registration through register_case.register (footprint declared producer-default:
   the device path binds no retained-etch.csv; InventoryStatus DeviceDerived; the mesh
   recipe every trace-basis case of the manifest binds, or --mesh-recipe).

Near-key reuse (USER decision 428 = Option A; decision 431; DESIGN v2; nearkey_reuse.py): with
`--nearkey-reuse default | fallback` (default `off`) every spatial coupon's generated basis is
offered to the near-key rule BEFORE registration: a qualifying library donor (the four calibrated
structure keys of RuleVersion v1, every item of DESIGN 1.2, the transplant gates, the predicted
bound inside the policy) gives a REUSED model under output/reused/<case>-reused/ (status
ReusedResponse; never registered / built here; `nearkey_reuse.py assemble` adds it to a
library) and the coupon is recorded reused; a refusal is recorded and the coupon follows the
normal path. `default` needs the rule file's DefaultActivation record (validation pair 5, shipped
beside the rule file since decision 459; fail closed when unreadable);
`fallback` needs --nearkey-fallback-approval and a --nearkey-fallback-stop-record
HASH_PREFIX=PATH (the exact key's registration / build STOP record) per requirement it may
reuse; every other requirement follows the normal path.

usage: device_coupons.py DEVICE_CONFIG --palace PATH --output DIR [--manifest PATH]
       [--mesh-recipe REPOSITORY_PATH] [--ring-size N] [--cap-triangulation METHOD] [--register]
"""
import argparse
import concurrent.futures
import copy
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
CPW2D = HERE.parents[1] / "cpw2d"
for path in (str(HERE), str(CPW2D)):
    if path not in sys.path:
        sys.path.insert(0, path)
import cluster_signature_geometry  # noqa: E402
import deterministic_math  # noqa: E402  (surface_response_identification, on sys.path through cluster_signature_geometry)
import general_mesh_manifest  # noqa: E402
import nearkey_predictor  # noqa: E402
import nearkey_reuse  # noqa: E402
import prepare_surface_response_coupons as planner  # noqa: E402
import register_case  # noqa: E402
from refreeze_manifest_tools import PRODUCTION_MANIFEST  # noqa: E402

DISCOVERY = CPW2D / "discover_surface_response_requirements.py"
GENERATOR = HERE / "generate_spatial_response.py"
# Round-3 class (10) (decisions 492 / 493 / 510; DESIGN-part-G G.10.4 rule 3): provenance.json (outside
# the content hash) names the float-serialisation rule of the generated sources, so a case id's
# generation rule is self-describing without changing the id of a straight coupon.
FLOAT_SERIALISATION_MODULE = Path(deterministic_math.__file__).resolve()
SPATIAL_METHOD = "SpatialCoupon"
DEFAULT_RING_SIZE = 16   # prepare_surface_response_coupons --spatial-ring-size default
DEFAULT_CAP_TRIANGULATION = "delaunay"   # generate_spatial_response --cap-triangulation (decisions 54b / 57)
# generate_spatial_response --cap-interior-spacing (x R): interior cap hats within R of the claimed portions at the
# ring spacing (decision 112(b): the ring-only JJ coupons interpolated the cap potential across the whole cap).
DEFAULT_CAP_INTERIOR_SPACING = 1.0
CAP_TRIANGULATIONS = ("ear-clipping", "delaunay")
INVENTORY_STATUS = "DeviceDerived"
# The bound source roles whose digests name a device coupon's directory (the manifest's
# content identity of a case: register_case.source_digests without the derived contract).
CONTENT_ROLES = {**register_case.REQUIRED_SOURCE_FILES,
                 **{role: name for role, name in register_case.OPTIONAL_SOURCE_FILES.items() if role != "Provenance"}}
DEVICE_RECORD = "device-coupons.json"


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


class DeviceAdapterError(ValueError):
    """A fail-closed stop of the device -> coupon mapping (the missing piece is named)."""


def run_discovery(device_config, output, palace, python=sys.executable, omit_requirements=()):
    """The closure manifest of the device (discover_surface_response_requirements.py);
    `omit_requirements` = its --omit-requirement Hash prefixes (no placeholder, stays Missing)."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    command = [python, str(DISCOVERY), str(device_config), "--output", str(output), "--palace", str(palace)]
    for prefix in omit_requirements:
        command += ["--omit-requirement", str(prefix)]
    with open(output / "discovery.log", "w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    manifest_path = output / "surface-response-requirements.json"
    if result.returncode != 0 or not manifest_path.is_file():
        raise DeviceAdapterError(f"discovery failed (rc {result.returncode}); see {output / 'discovery.log'}")
    return manifest_path


def shared_mesh_recipe(manifest):
    """The repository path of the mesh recipe every trace-basis case of the manifest
    binds (fail closed when they differ or none exists)."""
    paths = {case["Source"]["Files"]["MeshRecipe"].get("RepositoryPath") for case in manifest["Cases"]
             if "BasisContract" in case["Source"]["Files"]}
    paths.discard(None)
    if len(paths) != 1:
        raise DeviceAdapterError(f"the manifest's trace-basis cases bind {len(paths)} mesh recipe paths {sorted(paths)}: "
                                 f"pass --mesh-recipe")
    return paths.pop()


def coupon_geometry(coupon, radius, parameters=None):
    """The coupon written for the generator. A version-2 record (Geometry.Signature, no
    Edges: the v2 cluster contract) is built in the canonical frame from its signature alone
    (cluster_signature_geometry.cluster_coupon: edge rows from the portions, the plan-view
    mask from the arrangement of the extended chains); a version-1 record keeps the
    requirement with the planner's canonical plan-view boundary and mask regularization
    (prepare_surface_response_coupons.spatial_spec)."""
    generated = copy.deepcopy(coupon)
    geometry = generated["Geometry"]
    if coupon["Topology"] == "SpatialEdgeCluster" and "Signature" in geometry and not geometry.get("Edges"):
        if parameters is None:
            raise DeviceAdapterError(f"coupon {coupon['Id']}: the signature geometry needs the process parameters")
        try:
            generated, _ = cluster_signature_geometry.cluster_coupon(
                generated, radius, parameters["metal_thickness"], parameters["overetch"])
        except cluster_signature_geometry.SignatureGeometryError as error:
            raise DeviceAdapterError(f"coupon {coupon['Id']}: no coupon geometry from its signature: {error}") from error
        return generated
    facets = geometry.get("PlanViewFacets", [])
    if not facets:
        raise DeviceAdapterError(f"coupon {coupon['Id']} carries no PlanViewFacets: the Gmsh-only recipe needs the "
                                 f"plan-view mask and boundary (plan-view-mask.csv, plan-view-boundary.csv)")
    process_axis = 1 if coupon["Topology"] == "SpatialEdgeCluster" else 2
    canonical = planner.canonical_plan_view_boundary(facets, radius, process_axis)
    if "PlanViewBoundary" in geometry:
        if planner.unclassified_plan_view_boundary(geometry["PlanViewBoundary"]) != canonical:
            raise DeviceAdapterError(f"coupon {coupon['Id']}: Palace PlanViewBoundary does not match its exported facets")
    else:
        geometry["PlanViewBoundary"] = canonical
    geometry["MaskRegularization"] = {"Version": 1, "PhysicalBoundary": "TaperAndRound", "ContinuationBoundary": "Vertical"}
    return {"Topology": coupon["Topology"], "Geometry": geometry, "Interfaces": coupon["Interfaces"],
            "BoundaryCondition": coupon["BoundaryCondition"]}


def stamp_signature_model(library_path, coupon, radius, generated=None, mirror_formed=None):
    """A version-2 coupon's model carries the record's canonical Signature (the matcher's
    exact key for a SpatialEdgeCluster) and, as Edges, the claimed portions exactly (what
    ModelClusterSignature canonicalises to the feature's frame; the generator wrote the
    lengthened rows of the mesh geometry there). Version-1 coupons are left as written.
    `mirror_formed` (the validated MirrorFormedContract of a mirror-formed cluster, decision
    557) stamps the model's `MirrorFormed` record and every Edge's Weight (1.0 real / 0.0 image)."""
    geometry = coupon["Geometry"]
    if coupon["Topology"] != "SpatialEdgeCluster" or "Signature" not in geometry or geometry.get("Edges"):
        if mirror_formed is not None:
            raise DeviceAdapterError(f"coupon {coupon.get('Id')}: a mirror-formed contract needs a version-2 signature coupon")
        return False
    library_path = Path(library_path)
    library = json.loads(library_path.read_text())
    if len(library.get("Models", [])) != 1:
        raise DeviceAdapterError(f"{library_path}: the generator wrote {len(library.get('Models', []))} models, not one")
    model = library["Models"][0]
    model["Signature"] = geometry["Signature"]
    model["Edges"] = cluster_signature_geometry.model_edges(coupon, radius)
    if mirror_formed is not None:
        portions = cluster_signature_geometry.portions_from_signature(geometry["Signature"], radius)
        if len(portions) != len(model["Edges"]):
            raise DeviceAdapterError(f"{library_path}: {len(model['Edges'])} model edges for {len(portions)} chorded portions")
        record, weights = mirror_formed_entry(mirror_formed, [portion["Portion"] for portion in portions], radius)
        for edge, weight in zip(model["Edges"], weights):
            edge["Weight"] = weight
        model["MirrorFormed"] = record
    built = (generated or {}).get("Geometry", {})
    if built.get("SupportBox") is not None:
        # A device-plan coupon (decision 282): the generator's box (the signature's Box,
        # mesh units, canonical frame), the context rows (the mesher's / config's further
        # conductors) and the foreign edges the run config excludes from the within-R
        # accounting (EdgeExcludeSegments, rule B5; case_inputs.derive).
        model["SupportBox"] = list(built["SupportBox"])
        model["ForeignEdges"] = [list(segment) for segment in built.get("ForeignEdges", [])]
        model["ContextEdges"] = [dict(row) for row in built["Edges"] if row.get("Context")]
        model["ContextEdgeCount"] = len(model["ContextEdges"])
    library_path.write_text(json.dumps(library, indent=2) + "\n")
    return True


def generate_spatial_response_default_span_cap():
    import generate_spatial_response
    return generate_spatial_response.DEFAULT_SUPPORT_SPAN_CAP_OVER_R


def generate_sources(coupon, work, *, radius, parameters, ring_size, cap_triangulation=DEFAULT_CAP_TRIANGULATION,
                     cap_interior_spacing=DEFAULT_CAP_INTERIOR_SPACING, python=sys.executable, support_span_cap=None,
                     mirror_formed=None):
    """generate_spatial_response.py --basis-only into `work`; returns the generator command.
    `support_span_cap` (x R) raises the generator's matching-support span bound for this
    coupon alone (--support-span-cap; None = the generator default). `mirror_formed` (the
    validated MirrorFormedContract) stamps the model's real / image split after the generation;
    the generator itself reads the full signature unchanged."""
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    coupon_path = work / "coupon.json"
    generated = coupon_geometry(coupon, radius, parameters)
    coupon_path.write_text(json.dumps(generated, indent=2) + "\n")
    command = [python, str(GENERATOR), str(coupon_path), "--output", str(work), "--radius", str(radius),
               "--metal-thickness", str(parameters["metal_thickness"]), "--overetch-depth", str(parameters["overetch"]),
               "--sidewall-angle", str(parameters["sidewall_angle"]), "--top-rounding", str(parameters["top_radius"]),
               "--trench-rounding", str(parameters["bottom_radius"]), "--ring-size", str(ring_size),
               "--cap-triangulation", cap_triangulation, "--cap-interior-spacing", str(cap_interior_spacing),
               "--order", "1", "--model-name", coupon["Id"], "--basis-only"]
    if support_span_cap is not None:
        command += ["--support-span-cap", str(support_span_cap)]
    command += [str(item) for item in planner.material_options(parameters)]
    with open(work / "generate.log", "w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        failure = work / "generation-failure.json"
        reason = json.loads(failure.read_text()).get("Reason") if failure.is_file() else f"rc {result.returncode}"
        raise DeviceAdapterError(f"coupon {coupon['Id']}: the spatial generator stopped: {reason} (see {work / 'generate.log'})")
    stamp_signature_model(work / "process-library.json", coupon, radius, generated, mirror_formed=mirror_formed)
    return command


def write_process_toml(path, parameters, radius):
    process = {"Units": "um", "Radius": float(radius), "MetalThickness": parameters["metal_thickness"],
               "Overetch": parameters["overetch"], "SidewallAngle": parameters["sidewall_angle"],
               "TopRounding": parameters["top_radius"], "TrenchRounding": parameters["bottom_radius"]}
    Path(path).write_text("\n".join(f"{key} = {json.dumps(value)}" for key, value in process.items()) + "\n")
    return process


def content_hash(directory):
    """SHA-256 over (role, digest) of every bound source file present (sorted by role)."""
    directory = Path(directory)
    digests = {role: sha256(directory / name) for role, name in sorted(CONTENT_ROLES.items()) if (directory / name).is_file()}
    payload = json.dumps(digests, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest(), digests


OMITTED_METHOD = "OmittedRequirement"


# The mirror-formed cluster REQUIREMENT CONTRACT (decision 557; round-3 DESIGN 4.4 (MF), part G G.C.1 /
# X11; the spec: coupon-accuracy-assessment-20260913/mesher-design-round3-20261007/impl-B5/CONTRACT.md).
# The identification EMITS a mirror-formed SpatialEdgeCluster as an ordinary Missing requirement
# whose Signature is the FULL symmetric signature of the extended chain (real + image portions),
# flagged `MirrorFormed: true` and carrying `MirrorFormedContract` {Version, Planes, Frame,
# RealPortions, RealLengthOverR, ImageLengthOverR, RealFeatures, ExtendedFeature, Rule}.  This
# adapter CONSUMES it, routing on the CONTRACT (decision 562 MAJOR-1): the flag `MirrorFormed:
# true` alone is NOT "a mirror-formed key" - the C++ sets it on every feature whose Mirror record
# is not Continued (formed OR touched by an Unmerged configuration) and ORs it over the group, so
# a REAL spatial key touched by a plane carries it; such a flag-only requirement keeps today's
# path (built as the ordinary real coupon it is) and is recorded as information.  A requirement
# CARRYING `MirrorFormedContract`: by default (`refuse`) it is recorded out of scope (today's
# bytes); with `--mirror-formed admit` the contract is validated (fail closed by field), the
# coupon is generated from the full signature unchanged (an image portion is an ordinary claim
# of the generator) and the library ENTRY (the model of process-library.json) is stamped with
# the real / image split: every model Edge gets Weight 1 (real) / 0 (image) so the placement
# applies the model on the real half only, never an image claim.
MIRROR_FORMED_MODES = ("refuse", "admit")
MIRROR_FORMED_FLAG_ONLY_RULE = ("decision 562 MAJOR-1: MirrorFormed: true without a MirrorFormedContract is a REAL key touched by a "
                                "mirror configuration (the identification's flag means formed OR touched, ORed over the group; "
                                "surfaceresponseoperator.cpp instance.mirror_formed / surfaceresponsemirror.cpp Status Unmerged on "
                                "the touched real features): built as the ordinary real coupon, recorded as information")
MIRROR_FORMED_NOT_ADMITTED_METHOD = "MirrorFormedNotAdmitted"
MIRROR_FORMED_CONTRACT_VERSION = 1
# The contract's RealLengthOverR / ImageLengthOverR must agree with the signature's own portion
# lengths within the identification's arc-fit tolerance (the serialised ends are 1e-6 R quantised;
# an arc's length is read on the serialised circle).
MIRROR_FORMED_LENGTH_TOLERANCE_OVER_R = cluster_signature_geometry.ARC_FIT_TOLERANCE_OVER_R
MIRROR_FORMED_FRAME_TOLERANCE = 1.0e-9
MIRROR_FORMED_ENTRY_RULE = ("decision 557 / round-3 DESIGN 4.4 (MF): the library entry of a mirror-formed cluster built from "
                            "the FULL symmetric signature; RealPortions index Signature.Portions (0-based), EdgePortions names "
                            "the portion of every model Edge (an arc portion's chords share it), Edges[].Weight is 1.0 on a "
                            "real portion and 0.0 on an image portion: the placement applies the model on the real half only; "
                            "RealLengthFraction = RealLengthOverR / (RealLengthOverR + ImageLengthOverR)")


def signature_key_hash(signature):
    """sha256 of the identification's Key of a signature (the nlohmann compact dump with sorted
    keys = json.dumps(signature, separators=(',', ':'), sort_keys=True)): the requirement Hash."""
    return hashlib.sha256(json.dumps(signature, separators=(",", ":"), sort_keys=True).encode()).hexdigest()


def signature_portion_lengths_over_R(signature):
    """The length over R of every Signature.Portions entry: a straight portion its chord, an arc
    portion its length on the serialised circle (|sweep| x radius about the serialised centre
    through the serialised midpoint; 2 pi r for a closed circle)."""
    lengths = []
    for entry in signature["Portions"]:
        p = [float(v) for v in entry["P"]]
        if "Arc" in entry:
            arc = [float(v) for v in entry["Arc"]]
            sweep, _ = cluster_signature_geometry._arc_sweep(p[:2], p[2:], arc[:2], arc[2:])
            lengths.append(abs(sweep) * math.hypot(p[0] - arc[0], p[1] - arc[1]))
        else:
            lengths.append(math.hypot(p[2] - p[0], p[3] - p[1]))
    return lengths


def _contract_error(coupon, text):
    return DeviceAdapterError(f"coupon {coupon.get('Id')}: MirrorFormedContract {text} (impl-B5 CONTRACT.md section 2)")


def validate_mirror_formed_contract(coupon):
    """The MirrorFormedContract of a planned mirror-formed SpatialEdgeCluster coupon, validated
    field by field (CONTRACT.md section 2; every failure a named DeviceAdapterError): returns
    {RealPortions, ImagePortions, Planes, Frame, RealLengthOverR, ImageLengthOverR, RealFeatures,
    ExtendedFeature, Version}."""
    if coupon.get("Topology") != "SpatialEdgeCluster":
        raise _contract_error(coupon, f"is consumed for SpatialEdgeCluster only, not {coupon.get('Topology')!r}")
    if coupon.get("MirrorFormed") is not True:
        raise _contract_error(coupon, "needs MirrorFormed: true on the requirement")
    contract = coupon.get("MirrorFormedContract")
    if not isinstance(contract, dict):
        raise _contract_error(coupon, f"is not an object ({contract!r}): a requirement routed on the contract must carry the "
                                      "identification's record (a flag-only requirement never reaches this validation, decision 562)")
    if contract.get("Version") != MIRROR_FORMED_CONTRACT_VERSION:
        raise _contract_error(coupon, f"Version {contract.get('Version')!r} is not {MIRROR_FORMED_CONTRACT_VERSION}")
    signature = coupon.get("Signature")
    geometry = coupon.get("Geometry") or {}
    if not isinstance(signature, dict) or signature.get("Type") != "SpatialEdgeCluster" or not signature.get("Portions"):
        raise _contract_error(coupon, "needs a SpatialEdgeCluster Signature with Portions on the requirement")
    if geometry.get("Signature") != signature:
        raise _contract_error(coupon, "Geometry.Signature differs from Signature")
    portion_count = len(signature["Portions"])
    if int(geometry.get("EdgeCount", -1)) != portion_count:
        raise _contract_error(coupon, f"Geometry.EdgeCount {geometry.get('EdgeCount')!r} is not the portion count {portion_count}")
    expected_hash = signature_key_hash(signature)
    if coupon.get("Hash") != expected_hash:
        raise _contract_error(coupon, f"Hash {str(coupon.get('Hash'))[:12]} is not sha256(Key) {expected_hash[:12]}")
    real = contract.get("RealPortions")
    if (not isinstance(real, list) or not real or any(isinstance(i, bool) or not isinstance(i, int) for i in real)
            or sorted(set(real)) != real):
        raise _contract_error(coupon, f"RealPortions {real!r} must be a non-empty sorted list of distinct ints")
    if real[0] < 0 or real[-1] >= portion_count:
        raise _contract_error(coupon, f"RealPortions {real!r} index outside the {portion_count} portions")
    if len(real) == portion_count:
        raise _contract_error(coupon, f"RealPortions {real!r} names every portion: not a mirror-formed cluster")
    image = [i for i in range(portion_count) if i not in real]
    planes = contract.get("Planes")
    if (not isinstance(planes, list) or not planes or any(isinstance(k, bool) or not isinstance(k, int) or k < 0 for k in planes)
            or len(set(planes)) != len(planes)):
        raise _contract_error(coupon, f"Planes {planes!r} must be a non-empty list of distinct plane indices")
    frame = contract.get("Frame")
    if not isinstance(frame, dict) or not isinstance(frame.get("Origin"), list) or not isinstance(frame.get("Axes"), list):
        raise _contract_error(coupon, "Frame must carry Origin [3] and Axes [3][3]")
    origin = frame["Origin"]
    axes = frame["Axes"]
    if len(origin) != 3 or len(axes) != 3 or any(not isinstance(axis, list) or len(axis) != 3 for axis in axes):
        raise _contract_error(coupon, "Frame must carry Origin [3] and Axes [3][3]")
    try:
        origin = [float(v) for v in origin]
        axes = [[float(v) for v in axis] for axis in axes]
    except (TypeError, ValueError) as error:
        raise _contract_error(coupon, f"Frame is not numeric: {error}") from error
    if not all(math.isfinite(v) for v in origin + [v for axis in axes for v in axis]):
        raise _contract_error(coupon, "Frame is not finite")
    # Decisions on the frame only (nothing of it is written): scalar 3-vector arithmetic.
    dot3 = lambda a, b: a[0] * b[0] + a[1] * b[1] + a[2] * b[2]  # noqa: E731
    for i, axis in enumerate(axes):
        if abs(math.sqrt(dot3(axis, axis)) - 1.0) > MIRROR_FORMED_FRAME_TOLERANCE:
            raise _contract_error(coupon, f"Frame.Axes[{i}] is not a unit vector")
        for j in range(i):
            if abs(dot3(axes[j], axis)) > MIRROR_FORMED_FRAME_TOLERANCE:
                raise _contract_error(coupon, f"Frame.Axes[{j}] and [{i}] are not orthogonal")
    # CONTRACT.md v3 (decision 584 (2)): the Frame is the identification's canonical frame AS
    # IS - Axes[2] is the process normal for BOTH handedness values - with its handedness
    # explicit: Frame.Chirality (1 / -1; = Features[].Chirality for a chiral key, the recorded
    # frame's own handedness for a mirror-symmetric key of chirality 0 - the O1 wedge, S2p's
    # clusters), Axes[2] = Chirality x (Axes[0] x Axes[1]). A record without Chirality, or whose
    # triple disagrees with it, is refused by name (4 of the 16 production cluster contracts of
    # the rerun-2 windows are left-handed triples: C3 02db9a314b1b, S1b cb82f37fe3ba / 95e080429eb9
    # and O4 03fa3fb9166c).
    chirality = frame.get("Chirality")
    if isinstance(chirality, bool) or chirality not in (1, -1):
        raise _contract_error(coupon, f"Frame.Chirality {chirality!r} must be 1 or -1 (CONTRACT.md v3: Axes[2] = Chirality x "
                                      "(Axes[0] x Axes[1]))")
    x, y = axes[0], axes[1]
    cross = [x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]]
    if any(abs(chirality * cross[d] - axes[2][d]) > MIRROR_FORMED_FRAME_TOLERANCE for d in range(3)):
        raise _contract_error(coupon, f"Frame.Axes disagree with Chirality {chirality} (Axes[2] must be Chirality x (Axes[0] x "
                                      "Axes[1]))")
    lengths = signature_portion_lengths_over_R(signature)
    recorded = {}
    for key, indices in (("RealLengthOverR", real), ("ImageLengthOverR", image)):
        value = contract.get(key)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)) or float(value) <= 0.0:
            raise _contract_error(coupon, f"{key} {value!r} must be a positive finite number")
        own = sum(lengths[i] for i in indices)
        if abs(float(value) - own) > MIRROR_FORMED_LENGTH_TOLERANCE_OVER_R:
            raise _contract_error(coupon, f"{key} {float(value)!r} disagrees with the signature's portions {indices} "
                                          f"({own!r}) by more than {MIRROR_FORMED_LENGTH_TOLERANCE_OVER_R} R")
        recorded[key] = float(value)
    features = contract.get("RealFeatures")
    if not isinstance(features, list) or any(isinstance(f, bool) or not isinstance(f, int) for f in features):
        raise _contract_error(coupon, f"RealFeatures {features!r} must be a list of feature ids")
    extended = contract.get("ExtendedFeature")
    if isinstance(extended, bool) or not isinstance(extended, int):
        raise _contract_error(coupon, f"ExtendedFeature {extended!r} must be a feature id")
    if not isinstance(contract.get("Rule"), str) or not contract["Rule"].strip():
        raise _contract_error(coupon, "Rule text is missing")
    return {"Version": MIRROR_FORMED_CONTRACT_VERSION, "RealPortions": list(real), "ImagePortions": image, "Planes": list(planes),
            "Frame": {"Origin": origin, "Axes": axes, "Chirality": int(chirality)}, "RealLengthOverR": recorded["RealLengthOverR"],
            "ImageLengthOverR": recorded["ImageLengthOverR"], "RealFeatures": list(features), "ExtendedFeature": extended}


def mirror_formed_entry(contract, edge_portions, radius):
    """The ENTRY stamp of a mirror-formed cluster's model (CONTRACT.md section 3): the
    `MirrorFormed` record and the per-Edge Weight list (1.0 real / 0.0 image), from the validated
    contract and the signature Portion index of every model Edge (the chorded portions' order of
    cluster_signature_geometry.portions_from_signature = the order of model_edges)."""
    edge_portions = [int(index) for index in edge_portions]
    real = set(contract["RealPortions"])
    weights = [1.0 if index in real else 0.0 for index in edge_portions]
    if not any(weights) or all(weights):
        raise DeviceAdapterError(f"mirror-formed entry: the model edges carry no real / image split ({edge_portions})")
    total = contract["RealLengthOverR"] + contract["ImageLengthOverR"]
    record = {"Version": MIRROR_FORMED_CONTRACT_VERSION, "RealPortions": list(contract["RealPortions"]),
              "ImagePortions": list(contract["ImagePortions"]), "EdgePortions": edge_portions, "Planes": list(contract["Planes"]),
              "RealLengthOverR": contract["RealLengthOverR"], "ImageLengthOverR": contract["ImageLengthOverR"],
              "RealLengthFraction": contract["RealLengthOverR"] / total, "MatchingRadius": float(radius),
              "RealFeatures": list(contract["RealFeatures"]), "ExtendedFeature": contract["ExtendedFeature"],
              "Rule": MIRROR_FORMED_ENTRY_RULE}
    return record, weights


# Per-requirement options of `build --device` (HASH_PREFIX=VALUE, repeatable): a raised
# matching-support span cap for the generator (--support-span-cap, a single closed feature
# wider than 16R that cannot be split) and a raised element cap for the registered case
# (--element-cap, recorded as GateOverrides.MaximumElements of the fabricated and the thin
# case).  Each option must name exactly one spatial coupon of the discovery (fail closed:
# a prefix matching nothing, an omitted requirement or another family is an error) and
# carries its approval / reason text into the provenance and the manifest.
def parse_requirement_option(text, value_type, name):
    """HASH_PREFIX=VALUE -> (prefix, value); the prefix a non-empty hex token."""
    prefix, separator, value = str(text).partition("=")
    if not separator or not prefix or any(character not in "0123456789abcdef" for character in prefix):
        raise DeviceAdapterError(f"{name} expects HASH_PREFIX=VALUE with a hex prefix, not {text!r}")
    try:
        parsed = value_type(value)
    except ValueError as error:
        raise DeviceAdapterError(f"{name} {text!r}: {error}") from error
    if parsed <= 0:
        raise DeviceAdapterError(f"{name} {text!r}: the value must be positive")
    return prefix, parsed


def requirement_options(options, value_type, name):
    """The parsed HASH_PREFIX=VALUE options as {prefix: value}; a repeated prefix is an error."""
    parsed = {}
    for text in options or ():
        prefix, value = parse_requirement_option(text, value_type, name)
        if prefix in parsed:
            raise DeviceAdapterError(f"{name}: the prefix {prefix} is given twice")
        parsed[prefix] = value
    return parsed


def requirement_option_for(options, coupon_hash):
    """The (prefix, value) of the option naming this requirement Hash, or None."""
    matches = [(prefix, value) for prefix, value in options.items() if str(coupon_hash).startswith(prefix)]
    if len(matches) > 1:
        raise DeviceAdapterError(f"requirement {coupon_hash[:12]} is named by several prefixes {sorted(p for p, _ in matches)}")
    return matches[0] if matches else None


def check_requirement_options_consumed(options, consumed, name):
    unused = sorted(set(options) - set(consumed))
    if unused:
        raise DeviceAdapterError(f"{name}: no spatial coupon of the discovery has a requirement Hash starting with "
                                 f"{unused} (an omitted requirement or another family cannot carry it)")


NEARKEY_MODES = ("off", "default", "fallback")


def prepare_device_sources(device_config, *, palace, output, manifest_path=PRODUCTION_MANIFEST, ring_size=DEFAULT_RING_SIZE,
                           cap_triangulation=DEFAULT_CAP_TRIANGULATION,
                           cap_interior_spacing=DEFAULT_CAP_INTERIOR_SPACING, python=sys.executable, log=print,
                           omit_requirements=(), support_span_caps=(), support_span_cap_reason=None, element_caps=(),
                           element_cap_approval=None, element_cap_reason=None, nearkey_reuse_mode="off",
                           nearkey_fallback_approval=None, nearkey_fallback_stop_records=(), nearkey_rule=None, nearkey_record_roots=(),
                           mirror_formed="refuse"):
    """Steps 1-3: the source directories of every spatial coupon of the device under
    output/sources/<case id>; returns the device record (written to output/device-coupons.json).
    `omit_requirements`: Hash prefixes of Missing requirements the discovery gives no
    placeholder (they stay Missing); an omitted requirement is never built here, whatever its
    family: it is recorded out of scope with Method OmittedRequirement.
    `support_span_caps` / `element_caps`: HASH_PREFIX=VALUE options (parse_requirement_option)
    naming one spatial coupon each - the generator's --support-span-cap (x R) for it, and the
    element cap its registered cases carry as GateOverrides.MaximumElements (register_device_sources)
    - with their reason (and approval) texts, all recorded in the provenance and the record.
    `nearkey_reuse_mode` off | default | fallback (nearkey_reuse.py; the module docstring): a reused
    coupon is recorded under NearKeyReuse and is not registered.
    `mirror_formed` refuse | admit (decisions 557 / 562; the contract block above): a spatial coupon
    whose requirement CARRIES `MirrorFormedContract` is recorded out of scope (refuse, the default)
    or built from its full signature with the validated contract stamped on its entry (admit); a
    requirement flagged `MirrorFormed: true` without a contract is built as the ordinary real coupon
    (recorded under MirrorFormed.FlagOnly)."""
    span_caps = requirement_options(support_span_caps, float, "--support-span-cap")
    caps = requirement_options(element_caps, int, "--element-cap")
    if mirror_formed not in MIRROR_FORMED_MODES:
        raise DeviceAdapterError(f"--mirror-formed must be one of {MIRROR_FORMED_MODES}, not {mirror_formed!r}")
    if span_caps and not (isinstance(support_span_cap_reason, str) and support_span_cap_reason.strip()):
        raise DeviceAdapterError("--support-span-cap needs --support-span-cap-reason (recorded with the coupon)")
    if caps and not all(isinstance(text, str) and text.strip() for text in (element_cap_approval, element_cap_reason)):
        raise DeviceAdapterError("--element-cap needs --element-cap-approval and --element-cap-reason (recorded in the manifest)")
    if nearkey_reuse_mode not in NEARKEY_MODES:
        raise DeviceAdapterError(f"--nearkey-reuse must be one of {NEARKEY_MODES}, not {nearkey_reuse_mode!r}")
    stop_records = {}
    for text in nearkey_fallback_stop_records or ():
        prefix, separator, path = str(text).partition("=")
        if not separator or not prefix or not path or any(character not in "0123456789abcdef" for character in prefix):
            raise DeviceAdapterError(f"--nearkey-fallback-stop-record expects HASH_PREFIX=PATH with a hex prefix, not {text!r}")
        if prefix in stop_records:
            raise DeviceAdapterError(f"--nearkey-fallback-stop-record: the prefix {prefix} is given twice")
        stop_records[prefix] = path
    if nearkey_reuse_mode == "fallback":
        if not (isinstance(nearkey_fallback_approval, str) and nearkey_fallback_approval.strip()):
            raise DeviceAdapterError("--nearkey-reuse fallback needs --nearkey-fallback-approval TEXT (recorded on the reused model)")
        if not stop_records:
            raise DeviceAdapterError("--nearkey-reuse fallback needs a --nearkey-fallback-stop-record HASH_PREFIX=PATH per requirement")
    elif stop_records or nearkey_fallback_approval:
        raise DeviceAdapterError("--nearkey-fallback-approval / --nearkey-fallback-stop-record apply to --nearkey-reuse fallback only")
    record_roots = {}
    for text in nearkey_record_roots or ():
        remote, separator, local = str(text).partition("=")
        if not separator or not remote or not local:
            raise DeviceAdapterError(f"--nearkey-record-root expects REMOTE=LOCAL, not {text!r}")
        record_roots[remote] = local
    rule = None
    if nearkey_reuse_mode != "off":
        if nearkey_reuse_mode == "default" and nearkey_rule is not None:
            raise DeviceAdapterError("--nearkey-reuse default refuses --nearkey-rule: only the shipped, test-pinned rule file activates "
                                     "default reuse (decision 438 (5))")
        try:
            rule = nearkey_predictor.load_rule(nearkey_rule or nearkey_predictor.RULE_FILE)
            if nearkey_reuse_mode == "default" and not nearkey_predictor.default_active(rule)[0]:
                raise DeviceAdapterError("--nearkey-reuse default: the rule carries no DefaultActivation record (validation pair 5): "
                                         "default reuse is not active (DESIGN v2 section 3; decision 428)")
        except nearkey_predictor.NearKeyRuleError as error:
            raise DeviceAdapterError(f"near-key rule: {error}") from error
    consumed_stop = []
    consumed_span, consumed_caps = [], []
    device_config = Path(device_config).resolve()
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(manifest_path).resolve()
    commit = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, text=True).strip()
    log(f"device {device_config}: discovery with {palace}")
    closure = run_discovery(device_config, output / "discovery", palace, python, omit_requirements=omit_requirements)
    closure_manifest = json.loads(closure.read_text())
    omitted = {record["Hash"]: record for record in closure_manifest.get("OmittedRequirements", [])}
    library_path = Path(closure_manifest["Library"]["Path"])
    library = json.loads(library_path.read_text())
    parameters = planner.process_parameters(library)
    # The process seed's MatchingRadius is the process value; the closure manifest's is
    # Palace's re-serialization (L0 round trip: 1.9999999999999998 for 2.0) and must agree.
    radius = float(library["MatchingRadius"])
    reported = float(closure_manifest["Library"]["MatchingRadius"])
    if abs(reported - radius) > 1e-9 * radius:
        raise DeviceAdapterError(f"the closure manifest's MatchingRadius {reported} differs from the process library's {radius}")
    plan = planner.plan_from_manifest(closure, closure_manifest, library_path, library, include_matched=False)
    (output / "coupon-plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    record = {"Version": 1, "Command": "coupon-library build --device", "Commit": commit,
              "Device": {"Config": str(device_config), "SHA256": sha256(device_config)},
              "ProcessLibrary": {"Path": str(library_path), "SHA256": sha256(library_path), "MatchingRadius": radius,
                                 "Fabrication": library.get("Fabrication")},
              "Discovery": {"Manifest": str(closure), "SHA256": sha256(closure), "Summary": closure_manifest["Summary"],
                            "Complete": closure_manifest.get("Complete"),
                            "OmittedRequirements": [omitted[key] for key in sorted(omitted)]},
              "Plan": {"Path": str(output / "coupon-plan.json"), "Summary": plan["Summary"]},
              "TraceBasis": {"RingSize": ring_size, "CapTriangulation": cap_triangulation,
                             "DefaultCapTriangulation": DEFAULT_CAP_TRIANGULATION,
                             "CapInteriorSpacingOverR": cap_interior_spacing,
                             "CapInteriorRule": "generate_spatial_response --cap-interior-spacing: interior cap hats "
                                                "on a grid of this spacing (x R) within R of the claimed portions, "
                                                "appended after the ring knots (InteriorTraceCount); 0 = ring-only "
                                                "(the basis of the lane-J JJ coupons; decision 112(b))",
                             "Rule": "generate_spatial_response.build_matching_surface with the planner's default ring "
                                     "size (prepare_surface_response_coupons --spatial-ring-size), the basis every "
                                     "gallery case was produced with; CapTriangulation delaunay (the device default, "
                                     "decision 57) re-triangulates the two box caps without needle ears, ear-clipping "
                                     "is the gallery producer's (an explicit option; device coupons only, never a "
                                     "gallery reference)"},
              "NearKeyReuse": {"Mode": nearkey_reuse_mode, "RuleVersion": rule["RuleVersion"] if rule else None,
                               "RuleFileSHA256": rule["_sha256"] if rule else None, "Reused": [], "Refused": [],
                               "Rule": "USER decision 428 (Option A) / decision 431 / DESIGN v2: nearkey_reuse.reuse_requirement on the "
                                       "generated "
                                       "basis before registration; a reused coupon is not registered / built; a refusal is recorded"},
              "MirrorFormed": {"Mode": mirror_formed, "ContractVersion": MIRROR_FORMED_CONTRACT_VERSION, "FlagOnly": [],
                               "Rule": "decision 557 / round-3 DESIGN 4.4 (MF) / decision 562: a requirement CARRYING the "
                                       "identification's MirrorFormedContract is built from its FULL symmetric signature and its "
                                       "entry stamped with the real / image split (admit), or recorded out of scope with Method "
                                       f"{MIRROR_FORMED_NOT_ADMITTED_METHOD} (refuse, the default); a contract that fails its "
                                       "validation fails closed in admit mode; a requirement flagged MirrorFormed: true WITHOUT a "
                                       "contract (a real key touched by a mirror configuration) is built as the ordinary real coupon "
                                       "and listed under FlagOnly"},
              "Coupons": [], "OutOfScope": []}
    for coupon in plan["Coupons"]:
        method = coupon["Preparation"]["Method"]
        if coupon.get("Hash") in omitted:
            record["OutOfScope"].append({"Id": coupon["Id"], "Topology": coupon["Topology"], "Method": OMITTED_METHOD,
                                         "Reason": f"omitted by --omit-requirement {omitted[coupon['Hash']]['Prefix']}",
                                         "FamilyMethod": method, "Hash": coupon["Hash"],
                                         "DeviceOccurrences": coupon["DeviceOccurrences"],
                                         "DeviceEdgeLength": coupon["DeviceEdgeLength"],
                                         "Rule": "an omitted requirement got no discovery placeholder and is not built by "
                                                 "this call: it stays Missing against the library (recorded, never silent)"})
            continue
        if method != SPATIAL_METHOD:
            record["OutOfScope"].append({"Id": coupon["Id"], "Topology": coupon["Topology"], "Method": method,
                                         "Reason": coupon["Preparation"].get("Reason"),
                                         "DeviceOccurrences": coupon["DeviceOccurrences"],
                                         "DeviceEdgeLength": coupon["DeviceEdgeLength"],
                                         "Rule": "not a coupon of the Gmsh-only spatial library: built by its own family "
                                                 f"({method})"})
            continue
        contract = None
        flag_only = False
        if "MirrorFormedContract" in coupon:
            # The CONTRACT routes (decision 562 MAJOR-1), never the flag alone.
            if mirror_formed == "refuse":
                record["OutOfScope"].append({"Id": coupon["Id"], "Topology": coupon["Topology"],
                                             "Method": MIRROR_FORMED_NOT_ADMITTED_METHOD,
                                             "Reason": "a mirror-formed cluster requirement (MirrorFormedContract) is built only with "
                                                       "--mirror-formed admit (decision 557)",
                                             "FamilyMethod": method, "Hash": coupon.get("Hash"), "MirrorFormed": coupon.get("MirrorFormed"),
                                             "MirrorFormedContract": coupon["MirrorFormedContract"],
                                             "DeviceOccurrences": coupon["DeviceOccurrences"],
                                             "DeviceEdgeLength": coupon["DeviceEdgeLength"],
                                             "Rule": "a mirror-formed requirement is not built by default: it stays Missing against "
                                                     "the library (recorded, never silent)"})
                continue
            contract = validate_mirror_formed_contract(coupon)
        elif coupon.get("MirrorFormed"):
            flag_only = True
            record["MirrorFormed"]["FlagOnly"].append({"Id": coupon["Id"], "Topology": coupon["Topology"], "Hash": coupon.get("Hash"),
                                                       "Rule": MIRROR_FORMED_FLAG_ONLY_RULE})
        work = output / "work" / coupon["Id"]
        if work.exists():
            shutil.rmtree(work)
        span_cap = requirement_option_for(span_caps, coupon["Hash"])
        element_cap = requirement_option_for(caps, coupon["Hash"])
        if span_cap is not None:
            consumed_span.append(span_cap[0])
        if element_cap is not None:
            consumed_caps.append(element_cap[0])
        # The identification's per-case span-cap allowance (block (b) DESIGN section 3 / A4:
        # the requirement record's SpanCapAllowance, resolved by the quantum near-match on the
        # claims-only signature) is the generator's --support-span-cap unless a CLI option
        # names this requirement explicitly (the CLI option wins, both are recorded).
        allowance = coupon.get("SpanCapAllowance")
        if span_cap is None and allowance is not None:
            span_cap = (f"SpanCapAllowance:{allowance['Label']}", float(allowance["SpanCapOverR"]))
            span_cap_reason = (f"{allowance.get('Reason')} (identification SpanCapAllowance {allowance['Label']}, "
                               f"approval {allowance.get('Approval')}, MatchedQuanta {allowance.get('MatchedQuanta')})")
        else:
            span_cap_reason = support_span_cap_reason
        command = generate_sources(coupon, work, radius=radius, parameters=parameters, ring_size=ring_size,
                                   cap_triangulation=cap_triangulation, cap_interior_spacing=cap_interior_spacing,
                                   python=python, support_span_cap=None if span_cap is None else span_cap[1],
                                   mirror_formed=contract)
        write_process_toml(work / "process.toml", parameters, radius)
        digest, digests = content_hash(work)
        edge_count = int(coupon["Geometry"].get("EdgeCount", len(coupon["Geometry"].get("Edges", []))))
        case_id = f"spatial-{edge_count}-edge-{digest[:12]}"
        directory = output / "sources" / case_id
        near_key = None
        if rule is not None:
            stop_option = requirement_option_for(stop_records, coupon["Hash"]) if nearkey_reuse_mode == "fallback" else None
            if stop_option is not None:
                consumed_stop.append(stop_option[0])
            near_key = nearkey_reuse_for_coupon(coupon, work, case_id, library=library, library_path=library_path, rule=rule,
                                                mode=nearkey_reuse_mode, output=output / "reused",
                                                stop_record_path=None if stop_option is None else stop_option[1],
                                                approval=nearkey_fallback_approval, record_roots=record_roots, log=log)
            (record["NearKeyReuse"]["Reused"] if near_key["Reused"] else record["NearKeyReuse"]["Refused"]).append(
                {key: near_key[key] for key in ("Case", "Requirement", "Reason", "Donor", "ReuseMode", "ModelDirectory")})
        provenance = {
            "Version": 1, "Copyright": "Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.",
            "SPDX-License-Identifier": "Apache-2.0",
            "Origin": "coupon-library build --device (device_coupons.py)", "Commit": commit,
            "Device": record["Device"], "ProcessLibrary": {key: record["ProcessLibrary"][key] for key in ("Path", "SHA256")},
            "Discovery": {key: record["Discovery"][key] for key in ("Manifest", "SHA256")},
            "Requirement": {"Id": coupon["Id"], "Topology": coupon["Topology"], "EdgeCount": edge_count,
                            "Interfaces": coupon["Interfaces"], "BoundaryCondition": coupon["BoundaryCondition"],
                            "DeviceOccurrences": coupon["DeviceOccurrences"], "DeviceEdgeLength": coupon["DeviceEdgeLength"]},
            "Generator": {"Command": command, "RingSize": ring_size, "CapTriangulation": cap_triangulation,
                          "CapInteriorSpacingOverR": cap_interior_spacing,
                          "SupportSpanCapOverR": (None if span_cap is None else
                                                  {"Value": span_cap[1], "Prefix": span_cap[0],
                                                   "Default": generate_spatial_response_default_span_cap(),
                                                   "Reason": span_cap_reason,
                                                   "Allowance": allowance,
                                                   "Rule": "generate_spatial_response --support-span-cap: this coupon's "
                                                           "matching-support plan span may exceed the default bound "
                                                           "(basis-contract.json MatchingSupport records the span); "
                                                           "no other coupon's bound changes"}),
                          "PlanViewBoundary": ("cluster_signature_geometry.cluster_coupon: the coupon in the canonical "
                                               "frame of the version-2 Signature (edge rows from the portions, free ends "
                                               "extended to the box, mask = metal faces of the arrangement of the "
                                               "extended chains), boundary by canonical_plan_view_boundary (process "
                                               "axis 2); the model carries the Signature and the exact portions as Edges"
                                               if "Signature" in coupon["Geometry"] and not coupon["Geometry"].get("Edges")
                                               else "prepare_surface_response_coupons.canonical_plan_view_boundary of the "
                                                    "requirement's PlanViewFacets (process axis 1), MaskRegularization "
                                                    "TaperAndRound / Vertical (the planner's execute rule)")},
            "ContentHash": {"SHA256": digest, "Roles": digests, "Rule": "SHA-256 over the sorted (role, digest) pairs of "
                            "the bound source files; the case id is spatial-<edges>-edge-<first 12 hex>"},
            "FloatSerialisation": {
                "Rule": deterministic_math.FLOAT_SERIALISATION_RULE,
                "Module": {"Path": str(FLOAT_SERIALISATION_MODULE.relative_to(HERE.parents[2])),
                           "SHA256": sha256(FLOAT_SERIALISATION_MODULE)},
                "Statement": "every float written into a content-hashed source is produced by CPython scalar float "
                             "arithmetic (+ - x / sqrt, int, math.hypot; squares by multiplication, no ** / pow / log10) on "
                             "serialised inputs and by deterministic_math's "
                             "correctly rounded sin / cos / tan / atan2 / acos on the arc paths; numpy arrays carry data on "
                             "that path but do no BLAS arithmetic on it; the frame rotation is applied with scalar arithmetic "
                             "(round-3 class (10), decisions 492 / 493 / 510); a straight coupon reaches no transcendental "
                             "function and its id is unchanged"},
            "EtchFootprint": {"Declared": register_case.PRODUCER_DEFAULT_ETCH_FOOTPRINT,
                              "Reason": "the device path binds no retained-etch.csv: the producer-default footprint "
                                        "(3R collars around every conductor loop, cut by the retained mask) is a "
                                        "producer outcome (decision 16), recorded, never inferred"},
            "CopiedSourceFiles": sorted(name for name in CONTENT_ROLES.values() if (work / name).is_file()),
            "ExcludedArtifacts": ["traces/ (per-source basis CSVs: regenerated by qualify from the trace basis; their "
                                  "digests are the basis contract's OutputSourceSHA256)", "coupon.json", "generate.log"]}
        if contract is not None:
            # Decision 557: the mirror-formed admission and the contract it was built under (the
            # model's stamp is inside the content hash; this record is the provenance of it).
            provenance["MirrorFormed"] = {"Mode": mirror_formed, "Contract": coupon["MirrorFormedContract"],
                                          "RealPortions": contract["RealPortions"], "ImagePortions": contract["ImagePortions"],
                                          "Rule": MIRROR_FORMED_ENTRY_RULE}
        elif flag_only:
            # Decision 562: information only - the sources and the id are the ordinary real coupon's.
            provenance["MirrorFormed"] = {"Mode": mirror_formed, "FlagOnly": True, "Contract": None, "Rule": MIRROR_FORMED_FLAG_ONLY_RULE}
        if directory.exists():
            existing, _ = content_hash(directory)
            if existing != digest:
                raise DeviceAdapterError(f"{directory} exists with another content hash {existing[:12]}")
            status = "existing"
        else:
            directory.mkdir(parents=True)
            for name in provenance["CopiedSourceFiles"]:
                shutil.copyfile(work / name, directory / name)
            for name in ("zero-trace.csv", "basis-points.csv"):
                if (work / name).is_file():
                    shutil.copyfile(work / name, directory / name)
            for path in sorted(work.glob("conductor-*.csv")):
                shutil.copyfile(path, directory / path.name)
            (directory / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
            status = "written"
        built_model = json.loads((directory / "process-library.json").read_text())["Models"][0]
        contract = json.loads((directory / "basis-contract.json").read_text())
        record["Coupons"].append({"Case": case_id, "Requirement": coupon["Id"], "EdgeCount": edge_count,
                                  "ContextEdgeCount": int(built_model.get("ContextEdgeCount", 0)),
                                  "SupportBox": built_model.get("SupportBox"),
                                  "Directory": str(directory), "ContentHash": digest, "Status": status,
                                  "Sources": contract["Sources"], "MatchingSupport": contract.get("MatchingSupport"),
                                  "SupportSpanCapOverR": provenance["Generator"]["SupportSpanCapOverR"],
                                  "ElementCapOverride": (None if element_cap is None else
                                                         {"Value": element_cap[1], "Prefix": element_cap[0],
                                                          "Approval": element_cap_approval, "Reason": element_cap_reason}),
                                  "Interfaces": coupon["Interfaces"], "DeviceOccurrences": coupon["DeviceOccurrences"],
                                  "DeviceEdgeLength": coupon["DeviceEdgeLength"], "Registration": None, "NearKeyReuse": near_key,
                                  "MirrorFormed": None if contract is None else built_model.get("MirrorFormed"),
                                  "MirrorFormedFlagOnly": flag_only})
        log(f"{case_id}: source directory {status} ({edge_count} edges, requirement {coupon['Id']})"
            + (f"; near-key REUSED from {near_key['Donor']} ({near_key['ReuseMode']}): not registered"
               if near_key and near_key["Reused"] else ""))
    check_requirement_options_consumed(span_caps, consumed_span, "--support-span-cap")
    check_requirement_options_consumed(caps, consumed_caps, "--element-cap")
    if nearkey_reuse_mode == "fallback":
        check_requirement_options_consumed(stop_records, consumed_stop, "--nearkey-fallback-stop-record")
    record["Manifest"] = str(manifest_path)
    record["Output"] = str(output)
    (output / DEVICE_RECORD).write_text(json.dumps(record, indent=2) + "\n")
    return record


def nearkey_reuse_for_coupon(coupon, work, case_id, *, library, library_path, rule, mode, output, stop_record_path, approval,
                             record_roots=None, log=print):
    """The near-key reuse decision of one planned spatial coupon on its generated basis (`work`): the
    coupon record's NearKeyReuse block {Reused, Case, Requirement, Reason, Donor, ReuseMode, ModelDirectory,
    Record}. In fallback mode a requirement without a STOP record is not offered (recorded)."""
    base = {"Reused": False, "Case": case_id, "Requirement": coupon["Id"], "Reason": None, "Donor": None, "ReuseMode": None,
            "ModelDirectory": None, "Record": None}
    signature = coupon["Geometry"].get("Signature")
    if coupon["Topology"] != "SpatialEdgeCluster" or signature is None:
        return {**base, "Reason": "NotApplicable: no SpatialEdgeCluster Signature (a version-1 record)"}
    if mode == "fallback" and stop_record_path is None:
        return {**base, "Reason": "NoStopRecord: fallback reuse is offered only to a requirement with a --nearkey-fallback-stop-record"}
    entry = json.loads((Path(work) / "process-library.json").read_text())["Models"][0]
    try:
        stop = nearkey_reuse.stop_record_from_path(stop_record_path, entry["Name"], case_id) if stop_record_path else None
        result = nearkey_reuse.reuse_requirement(
            exact_signature=signature, exact_basis_dir=work, exact_model_entry=entry, library=library, library_path=library_path, rule=rule,
            mode=mode, output=output, requirement_key=coupon.get("Hash"), stop_record=stop, approval=approval, record_roots=record_roots,
            exact_interfaces=coupon.get("Interfaces"), exact_boundary_condition=coupon.get("BoundaryCondition"), log=log)
    except (nearkey_reuse.NearKeyReuseError, nearkey_predictor.NearKeyRuleError, nearkey_reuse.detection.NearKeyDetectionError,
            nearkey_reuse.transplant.TransplantError) as error:
        raise DeviceAdapterError(f"coupon {coupon['Id']}: near-key reuse stopped: {error}") from error
    record_path = Path(output) / f"{entry['Name']}-nearkey-reuse.json"
    record_path.write_text(json.dumps(result, indent=1) + "\n")
    if not result["Reused"]:
        return {**base, "Reason": result["Refused"]["Reason"], "Record": str(record_path)}
    return {**base, "Reused": True, "Reason": None, "Donor": result["Model"]["ReusedFrom"]["Donor"],
            "ReuseMode": result["Model"]["ReuseMode"],
            "ModelDirectory": result["ModelDirectory"], "Record": str(record_path)}


DEFAULT_REGISTER_JOBS = 2


REGISTRATION_KEYS = ("Status", "Message", "FixtureVersion", "ContractSHA256", "Scope", "StoppedBy", "Work")
STATUS_NEARKEY_REUSED = "near-key-reused"


def register_device_sources(record, *, manifest_path, mesh_recipe=None, work=None, python=sys.executable, julia=None,
                            log=print, probe=None, jobs=DEFAULT_REGISTER_JOBS, thin=True):
    """Step 4: register every produced source directory (idempotent by content); the
    registration record of each is attached to the device record and written back.
    The manifest-free part (digests, scope, the labels-only probe and the contract
    derivation: register_case.prepare_registration) runs for `jobs` coupons at once;
    the manifest append (register_case.commit_registration) is serial, in the device
    record's coupon order (decision 62(2)).  With `thin` (the default; decision 66) the
    thin counterpart of every registered / reused fabricated coupon (<case>-thin, Kind
    thin, the same source directory) is registered the same way afterwards - the device
    correction needs the fabricated AND the thin response of every model - and recorded
    under ThinCase / ThinRegistration."""
    manifest_path = Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    recipe = mesh_recipe or shared_mesh_recipe(manifest)
    if not isinstance(jobs, int) or isinstance(jobs, bool) or jobs < 1:
        raise DeviceAdapterError(f"the registration pool size must be an integer >= 1, not {jobs!r}")

    def gate_overrides(coupon):
        """The GateOverrides block of a coupon's cases (fabricated and thin alike) from its
        recorded --element-cap, validated against the manifest gate; None without one."""
        override = coupon.get("ElementCapOverride")
        if override is None:
            return None
        try:
            return general_mesh_manifest.element_cap_override(
                override["Value"], manifest["Gates"][general_mesh_manifest.ELEMENT_CAP_GATE],
                approval=override["Approval"], reason=override["Reason"])
        except ValueError as error:
            raise DeviceAdapterError(f"{coupon['Case']}: {error}") from error

    def prepare(coupon, kind="fabricated"):
        case_id = coupon["Case"] if kind == "fabricated" else register_case.thin_case_id(coupon["Case"])
        try:
            return register_case.prepare_registration(
                case_id, coupon["Directory"], footprint=register_case.PRODUCER_DEFAULT_ETCH_FOOTPRINT,
                inventory_status=INVENTORY_STATUS, manifest_path=manifest_path, mesh_recipe=recipe,
                provenance=f"device {record['Device']['Config']} (SHA256 {record['Device']['SHA256'][:12]}...), requirement "
                           f"{coupon['Requirement']}, content hash {coupon['ContentHash'][:12]}",
                work=(Path(work) / case_id) if work is not None else None, python=python, julia=julia,
                kind=kind, fabricated_case=coupon["Case"] if kind == "thin" else None,
                gate_overrides=gate_overrides(coupon),
                **({"probe": probe} if probe is not None else {}))
        except (register_case.RegistrationError, ValueError) as error:
            return error

    def outcome(case_id, prepared):
        try:
            if isinstance(prepared, Exception):
                raise prepared
            return register_case.commit_registration(prepared, manifest_path=manifest_path)
        except (register_case.RegistrationError, ValueError) as error:
            registration = {"Status": register_case.STATUS_FAILED, "Message": str(error)}
            failed = Path(work) / case_id / register_case.REGISTRATION_RECORD if work is not None else None
            if failed is not None and failed.is_file():
                registration = json.loads(failed.read_text())
            return registration

    # a near-key REUSED coupon (decision 428) is not registered / built: its reused model replaces the exact coupon
    reused = [coupon for coupon in record["Coupons"] if (coupon.get("NearKeyReuse") or {}).get("Reused")]
    for coupon in reused:
        coupon["Registration"] = {"Status": STATUS_NEARKEY_REUSED, "Message": f"near-key reused from {coupon['NearKeyReuse']['Donor']} "
                                  f"({coupon['NearKeyReuse']['ReuseMode']}): not registered (DESIGN v2 section 3)"}
        coupon["ThinCase"], coupon["ThinRegistration"] = None, None
        log(f"{coupon['Case']}: {coupon['Registration']['Message']}")
    coupons = [coupon for coupon in record["Coupons"] if coupon not in reused]
    started = time.monotonic()
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(jobs, max(len(coupons), 1))) as pool:
        prepared = list(pool.map(prepare, coupons))
    prepare_seconds = time.monotonic() - started
    for coupon, item in zip(coupons, prepared):
        registration = outcome(coupon["Case"], item)
        coupon["Registration"] = {key: registration.get(key) for key in REGISTRATION_KEYS}
        log(f"{coupon['Case']}: registration {registration.get('Status')} - {registration.get('Message')}")
    thin_seconds = None
    if thin:
        # The thin pair of every fabricated coupon now in the manifest (its prepare needs
        # the fabricated case registered: validate_case_kind binds the pair).
        eligible = [coupon for coupon in coupons if coupon["Registration"]["Status"] in
                    (register_case.STATUS_REGISTERED, register_case.STATUS_REUSED)]
        thin_started = time.monotonic()
        with concurrent.futures.ThreadPoolExecutor(max_workers=min(jobs, max(len(eligible), 1))) as pool:
            thin_prepared = list(pool.map(lambda coupon: prepare(coupon, "thin"), eligible))
        for coupon in coupons:
            coupon["ThinCase"] = register_case.thin_case_id(coupon["Case"]) if coupon in eligible else None
            coupon["ThinRegistration"] = None
        for coupon, item in zip(eligible, thin_prepared):
            registration = outcome(coupon["ThinCase"], item)
            coupon["ThinRegistration"] = {key: registration.get(key) for key in REGISTRATION_KEYS}
            log(f"{coupon['ThinCase']}: registration {registration.get('Status')} - {registration.get('Message')}")
        thin_seconds = time.monotonic() - thin_started
    record["MeshRecipe"] = recipe
    record["RegistrationPool"] = {"Jobs": jobs, "PrepareSeconds": prepare_seconds, "ThinSeconds": thin_seconds,
                                  "TotalSeconds": time.monotonic() - started,
                                  "Rule": "register_case.prepare_registration (digests, scope, labels-only probe, contract "
                                          "derivation) for Jobs coupons at once, then commit_registration serially in "
                                          "coupon order (decision 62(2)); then the thin counterparts likewise (decision 66)"}
    (Path(record["Output"]) / DEVICE_RECORD).write_text(json.dumps(record, indent=2) + "\n")
    return record


def add_requirement_option_arguments(parser):
    """The per-requirement options of `build --device` (device_coupons.py and coupon_library.py build)."""
    parser.add_argument("--support-span-cap", action="append", default=[], metavar="HASH_PREFIX=OVER_R",
                        help="with --device: the spatial coupon whose requirement Hash starts with the prefix is generated "
                             "with this matching-support span cap (x R; generate_spatial_response --support-span-cap, "
                             "default 16) - a single closed feature wider than 16R that cannot be split; repeatable; "
                             "needs --support-span-cap-reason; a prefix matching no spatial coupon fails closed")
    parser.add_argument("--support-span-cap-reason", help="why the span cap is raised (recorded in the provenance)")
    parser.add_argument("--element-cap", action="append", default=[], metavar="HASH_PREFIX=ELEMENTS",
                        help="with --device: the registered cases (fabricated + thin) of the spatial coupon whose requirement "
                             "Hash starts with the prefix carry this element cap as GateOverrides.MaximumElements (above the "
                             "manifest gate, which is unchanged for every other case); repeatable; needs "
                             "--element-cap-approval and --element-cap-reason")
    parser.add_argument("--element-cap-approval", help="who approved the element cap override (recorded in the manifest)")
    parser.add_argument("--element-cap-reason", help="why the element cap is raised (recorded in the manifest)")
    parser.add_argument("--nearkey-reuse", choices=NEARKEY_MODES, default="off",
                        help="with --device: offer every spatial coupon to the near-key reuse rule (USER decision 428, Option A; "
                             "nearkey_reuse.py) before registration: default (needs the rule's DefaultActivation record) or fallback "
                             "(needs --nearkey-fallback-approval and a --nearkey-fallback-stop-record per requirement); off = never")
    parser.add_argument("--nearkey-fallback-approval", help="fallback reuse: the recorded approval text (the supervisor decision)")
    parser.add_argument("--nearkey-fallback-stop-record", action="append", default=[], metavar="HASH_PREFIX=PATH",
                        help="fallback reuse: the registration / build STOP record (JSON) of the requirement whose Hash starts with the "
                             "prefix (repeatable; a prefix matching no spatial coupon fails closed)")
    parser.add_argument("--nearkey-rule", help="the near-key rule file (default: the repository's nearkey-reuse-rule-v1.json; refused "
                                               "in default mode, decision 438 (5))")
    parser.add_argument("--nearkey-record-root", action="append", default=[], metavar="REMOTE=LOCAL",
                        help="map the library's record paths (cluster) to a local mirror: the donors' stored (F) records are REQUIRED "
                             "for the T4 gate in default / fallback (decision 438 (2))")
    parser.add_argument("--mirror-formed", choices=MIRROR_FORMED_MODES, default="refuse",
                        help="with --device: admit = build a mirror-formed cluster requirement (one CARRYING the identification's "
                             "MirrorFormedContract, decisions 557 / 562) from its FULL symmetric signature and stamp its entry with "
                             "the real / image split (Edges[].Weight 1 / 0); refuse (default) = record it out of scope with Method "
                             f"{MIRROR_FORMED_NOT_ADMITTED_METHOD}; a MirrorFormed: true requirement without a contract is built as "
                             "the ordinary real coupon either way")


def requirement_option_kwargs(args):
    return {"support_span_caps": args.support_span_cap, "support_span_cap_reason": args.support_span_cap_reason,
            "element_caps": args.element_cap, "element_cap_approval": args.element_cap_approval,
            "element_cap_reason": args.element_cap_reason, "nearkey_reuse_mode": args.nearkey_reuse,
            "nearkey_fallback_approval": args.nearkey_fallback_approval,
            "nearkey_fallback_stop_records": args.nearkey_fallback_stop_record, "nearkey_rule": args.nearkey_rule,
            "nearkey_record_roots": args.nearkey_record_root, "mirror_formed": args.mirror_formed}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("device_config", type=Path)
    parser.add_argument("--palace", type=Path, required=True, help="Palace executable (geometry preflights only)")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, default=PRODUCTION_MANIFEST)
    parser.add_argument("--mesh-recipe", help="repository path of the frozen mesh recipe (default: the one every "
                                              "trace-basis case of the manifest binds)")
    parser.add_argument("--ring-size", type=int, default=DEFAULT_RING_SIZE)
    parser.add_argument("--cap-triangulation", choices=CAP_TRIANGULATIONS, default=DEFAULT_CAP_TRIANGULATION,
                        help="matching-box cap triangulation of the device basis (generate_spatial_response.py): "
                             "delaunay (default, decision 57) or ear-clipping (the gallery producer's)")
    parser.add_argument("--cap-interior-spacing", type=float, default=DEFAULT_CAP_INTERIOR_SPACING,
                        help="spacing (x R) of the interior cap hats within R of the claimed portions "
                             f"(generate_spatial_response.py --cap-interior-spacing; default {DEFAULT_CAP_INTERIOR_SPACING}; "
                             "0 = the ring-only basis)")
    parser.add_argument("--omit-requirement", action="append", default=[], metavar="HASH_PREFIX",
                        help="Missing requirement(s) whose Hash starts with this prefix get no discovery placeholder and "
                             "are not built (repeatable; recorded out of scope with Method OmittedRequirement)")
    add_requirement_option_arguments(parser)
    parser.add_argument("--register", action="store_true", help="register the produced directories into --manifest")
    parser.add_argument("--work", type=Path, help="parent of the registration work directories")
    parser.add_argument("--julia", default=shutil.which("julia"))
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--register-jobs", type=int, default=DEFAULT_REGISTER_JOBS,
                        help=f"coupons whose labels-only probe / contract derivation run at once (default {DEFAULT_REGISTER_JOBS}; "
                             "the manifest append stays serial)")
    args = parser.parse_args(argv)
    try:
        record = prepare_device_sources(args.device_config, palace=args.palace, output=args.output, manifest_path=args.manifest,
                                        ring_size=args.ring_size, cap_triangulation=args.cap_triangulation,
                                        cap_interior_spacing=args.cap_interior_spacing, python=args.python,
                                        omit_requirements=args.omit_requirement, **requirement_option_kwargs(args))
        if args.register:
            register_device_sources(record, manifest_path=args.manifest, mesh_recipe=args.mesh_recipe, work=args.work,
                                    python=args.python, julia=args.julia, jobs=args.register_jobs)
    except DeviceAdapterError as error:
        print(f"DEVICE_ADAPTER_FAILED: {error}", file=sys.stderr)
        return 1
    print(f"DEVICE {record['Device']['Config']}: {len(record['Coupons'])} spatial coupons, {len(record['OutOfScope'])} out of "
          f"scope; record {args.output / DEVICE_RECORD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
