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

usage: device_coupons.py DEVICE_CONFIG --palace PATH --output DIR [--manifest PATH]
       [--mesh-recipe REPOSITORY_PATH] [--ring-size N] [--cap-triangulation METHOD] [--register]
"""
import argparse
import concurrent.futures
import copy
import hashlib
import json
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
import general_mesh_manifest  # noqa: E402
import prepare_surface_response_coupons as planner  # noqa: E402
import register_case  # noqa: E402
from refreeze_manifest_tools import PRODUCTION_MANIFEST  # noqa: E402

DISCOVERY = CPW2D / "discover_surface_response_requirements.py"
GENERATOR = HERE / "generate_spatial_response.py"
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


def stamp_signature_model(library_path, coupon, radius, generated=None):
    """A version-2 coupon's model carries the record's canonical Signature (the matcher's
    exact key for a SpatialEdgeCluster) and, as Edges, the claimed portions exactly (what
    ModelClusterSignature canonicalises to the feature's frame; the generator wrote the
    lengthened rows of the mesh geometry there). Version-1 coupons are left as written."""
    geometry = coupon["Geometry"]
    if coupon["Topology"] != "SpatialEdgeCluster" or "Signature" not in geometry or geometry.get("Edges"):
        return False
    library_path = Path(library_path)
    library = json.loads(library_path.read_text())
    if len(library.get("Models", [])) != 1:
        raise DeviceAdapterError(f"{library_path}: the generator wrote {len(library.get('Models', []))} models, not one")
    model = library["Models"][0]
    model["Signature"] = geometry["Signature"]
    model["Edges"] = cluster_signature_geometry.model_edges(coupon, radius)
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
                     cap_interior_spacing=DEFAULT_CAP_INTERIOR_SPACING, python=sys.executable, support_span_cap=None):
    """generate_spatial_response.py --basis-only into `work`; returns the generator command.
    `support_span_cap` (x R) raises the generator's matching-support span bound for this
    coupon alone (--support-span-cap; None = the generator default)."""
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
    stamp_signature_model(work / "process-library.json", coupon, radius, generated)
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


def prepare_device_sources(device_config, *, palace, output, manifest_path=PRODUCTION_MANIFEST, ring_size=DEFAULT_RING_SIZE,
                           cap_triangulation=DEFAULT_CAP_TRIANGULATION,
                           cap_interior_spacing=DEFAULT_CAP_INTERIOR_SPACING, python=sys.executable, log=print,
                           omit_requirements=(), support_span_caps=(), support_span_cap_reason=None, element_caps=(),
                           element_cap_approval=None, element_cap_reason=None):
    """Steps 1-3: the source directories of every spatial coupon of the device under
    output/sources/<case id>; returns the device record (written to output/device-coupons.json).
    `omit_requirements`: Hash prefixes of Missing requirements the discovery gives no
    placeholder (they stay Missing); an omitted requirement is never built here, whatever its
    family: it is recorded out of scope with Method OmittedRequirement.
    `support_span_caps` / `element_caps`: HASH_PREFIX=VALUE options (parse_requirement_option)
    naming one spatial coupon each - the generator's --support-span-cap (x R) for it, and the
    element cap its registered cases carry as GateOverrides.MaximumElements (register_device_sources)
    - with their reason (and approval) texts, all recorded in the provenance and the record."""
    span_caps = requirement_options(support_span_caps, float, "--support-span-cap")
    caps = requirement_options(element_caps, int, "--element-cap")
    if span_caps and not (isinstance(support_span_cap_reason, str) and support_span_cap_reason.strip()):
        raise DeviceAdapterError("--support-span-cap needs --support-span-cap-reason (recorded with the coupon)")
    if caps and not all(isinstance(text, str) and text.strip() for text in (element_cap_approval, element_cap_reason)):
        raise DeviceAdapterError("--element-cap needs --element-cap-approval and --element-cap-reason (recorded in the manifest)")
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
        work = output / "work" / coupon["Id"]
        if work.exists():
            shutil.rmtree(work)
        span_cap = requirement_option_for(span_caps, coupon["Hash"])
        element_cap = requirement_option_for(caps, coupon["Hash"])
        if span_cap is not None:
            consumed_span.append(span_cap[0])
        if element_cap is not None:
            consumed_caps.append(element_cap[0])
        command = generate_sources(coupon, work, radius=radius, parameters=parameters, ring_size=ring_size,
                                   cap_triangulation=cap_triangulation, cap_interior_spacing=cap_interior_spacing,
                                   python=python, support_span_cap=None if span_cap is None else span_cap[1])
        write_process_toml(work / "process.toml", parameters, radius)
        digest, digests = content_hash(work)
        edge_count = int(coupon["Geometry"].get("EdgeCount", len(coupon["Geometry"].get("Edges", []))))
        case_id = f"spatial-{edge_count}-edge-{digest[:12]}"
        directory = output / "sources" / case_id
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
                                                   "Reason": support_span_cap_reason,
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
            "EtchFootprint": {"Declared": register_case.PRODUCER_DEFAULT_ETCH_FOOTPRINT,
                              "Reason": "the device path binds no retained-etch.csv: the producer-default footprint "
                                        "(3R collars around every conductor loop, cut by the retained mask) is a "
                                        "producer outcome (decision 16), recorded, never inferred"},
            "CopiedSourceFiles": sorted(name for name in CONTENT_ROLES.values() if (work / name).is_file()),
            "ExcludedArtifacts": ["traces/ (per-source basis CSVs: regenerated by qualify from the trace basis; their "
                                  "digests are the basis contract's OutputSourceSHA256)", "coupon.json", "generate.log"]}
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
                                  "DeviceEdgeLength": coupon["DeviceEdgeLength"], "Registration": None})
        log(f"{case_id}: source directory {status} ({edge_count} edges, requirement {coupon['Id']})")
    check_requirement_options_consumed(span_caps, consumed_span, "--support-span-cap")
    check_requirement_options_consumed(caps, consumed_caps, "--element-cap")
    record["Manifest"] = str(manifest_path)
    record["Output"] = str(output)
    (output / DEVICE_RECORD).write_text(json.dumps(record, indent=2) + "\n")
    return record


DEFAULT_REGISTER_JOBS = 2


REGISTRATION_KEYS = ("Status", "Message", "FixtureVersion", "ContractSHA256", "Scope", "StoppedBy", "Work")


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

    coupons = list(record["Coupons"])
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


def requirement_option_kwargs(args):
    return {"support_span_caps": args.support_span_cap, "support_span_cap_reason": args.support_span_cap_reason,
            "element_caps": args.element_cap, "element_cap_approval": args.element_cap_approval,
            "element_cap_reason": args.element_cap_reason}


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
