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
   and every tool.  The directory is named by the content hash of its bound source
   files (`spatial-<edge count>-edge-<hash12>`): the same device geometry always maps
   to the same case, and register_case.py reuses a case whose digests it already holds.
4. Registration through register_case.register (footprint declared producer-default:
   the device path binds no retained-etch.csv; InventoryStatus DeviceDerived; the mesh
   recipe every trace-basis case of the manifest binds, or --mesh-recipe).

usage: device_coupons.py DEVICE_CONFIG --palace PATH --output DIR [--manifest PATH]
       [--mesh-recipe REPOSITORY_PATH] [--ring-size N] [--register]
"""
import argparse
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
import prepare_surface_response_coupons as planner  # noqa: E402
import register_case  # noqa: E402
from refreeze_manifest_tools import PRODUCTION_MANIFEST  # noqa: E402

DISCOVERY = CPW2D / "discover_surface_response_requirements.py"
GENERATOR = HERE / "generate_spatial_response.py"
SPATIAL_METHOD = "SpatialCoupon"
DEFAULT_RING_SIZE = 16   # prepare_surface_response_coupons --spatial-ring-size default
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


def run_discovery(device_config, output, palace, python=sys.executable):
    """The closure manifest of the device (discover_surface_response_requirements.py)."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    command = [python, str(DISCOVERY), str(device_config), "--output", str(output), "--palace", str(palace)]
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


def coupon_geometry(coupon, radius):
    """The coupon written for the generator: the requirement with the planner's canonical
    plan-view boundary and mask regularization (prepare_surface_response_coupons.spatial_spec)."""
    generated = copy.deepcopy(coupon)
    geometry = generated["Geometry"]
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


def generate_sources(coupon, work, *, radius, parameters, ring_size, python=sys.executable):
    """generate_spatial_response.py --basis-only into `work`; returns the generator command."""
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    coupon_path = work / "coupon.json"
    coupon_path.write_text(json.dumps(coupon_geometry(coupon, radius), indent=2) + "\n")
    command = [python, str(GENERATOR), str(coupon_path), "--output", str(work), "--radius", str(radius),
               "--metal-thickness", str(parameters["metal_thickness"]), "--overetch-depth", str(parameters["overetch"]),
               "--sidewall-angle", str(parameters["sidewall_angle"]), "--top-rounding", str(parameters["top_radius"]),
               "--trench-rounding", str(parameters["bottom_radius"]), "--ring-size", str(ring_size), "--order", "1",
               "--model-name", coupon["Id"], "--basis-only"]
    command += [str(item) for item in planner.material_options(parameters)]
    with open(work / "generate.log", "w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        failure = work / "generation-failure.json"
        reason = json.loads(failure.read_text()).get("Reason") if failure.is_file() else f"rc {result.returncode}"
        raise DeviceAdapterError(f"coupon {coupon['Id']}: the spatial generator stopped: {reason} (see {work / 'generate.log'})")
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


def prepare_device_sources(device_config, *, palace, output, manifest_path=PRODUCTION_MANIFEST, ring_size=DEFAULT_RING_SIZE,
                           python=sys.executable, log=print):
    """Steps 1-3: the source directories of every spatial coupon of the device under
    output/sources/<case id>; returns the device record (written to output/device-coupons.json)."""
    device_config = Path(device_config).resolve()
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(manifest_path).resolve()
    commit = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, text=True).strip()
    log(f"device {device_config}: discovery with {palace}")
    closure = run_discovery(device_config, output / "discovery", palace, python)
    closure_manifest = json.loads(closure.read_text())
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
                            "Complete": closure_manifest.get("Complete")},
              "Plan": {"Path": str(output / "coupon-plan.json"), "Summary": plan["Summary"]},
              "TraceBasis": {"RingSize": ring_size, "Rule": "generate_spatial_response.build_matching_surface with the "
                             "planner's default ring size (prepare_surface_response_coupons --spatial-ring-size); "
                             "the basis every gallery case was produced with"},
              "Coupons": [], "OutOfScope": []}
    for coupon in plan["Coupons"]:
        method = coupon["Preparation"]["Method"]
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
        command = generate_sources(coupon, work, radius=radius, parameters=parameters, ring_size=ring_size, python=python)
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
            "Generator": {"Command": command, "RingSize": ring_size,
                          "PlanViewBoundary": "prepare_surface_response_coupons.canonical_plan_view_boundary of the "
                                              "requirement's PlanViewFacets (process axis 1), MaskRegularization "
                                              "TaperAndRound / Vertical (the planner's execute rule)"},
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
        record["Coupons"].append({"Case": case_id, "Requirement": coupon["Id"], "EdgeCount": edge_count,
                                  "Directory": str(directory), "ContentHash": digest, "Status": status,
                                  "Sources": json.loads((directory / "basis-contract.json").read_text())["Sources"],
                                  "Interfaces": coupon["Interfaces"], "DeviceOccurrences": coupon["DeviceOccurrences"],
                                  "DeviceEdgeLength": coupon["DeviceEdgeLength"], "Registration": None})
        log(f"{case_id}: source directory {status} ({edge_count} edges, requirement {coupon['Id']})")
    record["Manifest"] = str(manifest_path)
    record["Output"] = str(output)
    (output / DEVICE_RECORD).write_text(json.dumps(record, indent=2) + "\n")
    return record


def register_device_sources(record, *, manifest_path, mesh_recipe=None, work=None, python=sys.executable, julia=None,
                            log=print, probe=None):
    """Step 4: register every produced source directory (idempotent by content); the
    registration record of each is attached to the device record and written back."""
    manifest_path = Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    recipe = mesh_recipe or shared_mesh_recipe(manifest)
    for coupon in record["Coupons"]:
        case_id = coupon["Case"]
        try:
            registration = register_case.register(
                case_id, coupon["Directory"], footprint=register_case.PRODUCER_DEFAULT_ETCH_FOOTPRINT,
                inventory_status=INVENTORY_STATUS, manifest_path=manifest_path, mesh_recipe=recipe,
                provenance=f"device {record['Device']['Config']} (SHA256 {record['Device']['SHA256'][:12]}...), requirement "
                           f"{coupon['Requirement']}, content hash {coupon['ContentHash'][:12]}",
                work=(Path(work) / case_id) if work is not None else None, python=python, julia=julia,
                **({"probe": probe} if probe is not None else {}))
        except (register_case.RegistrationError, ValueError) as error:
            registration = {"Status": register_case.STATUS_FAILED, "Message": str(error)}
            failed = Path(work) / case_id / register_case.REGISTRATION_RECORD if work is not None else None
            if failed is not None and failed.is_file():
                registration = json.loads(failed.read_text())
        coupon["Registration"] = {key: registration.get(key) for key in
                                  ("Status", "Message", "FixtureVersion", "ContractSHA256", "Scope", "StoppedBy", "Work")}
        log(f"{case_id}: registration {registration.get('Status')} - {registration.get('Message')}")
    record["MeshRecipe"] = recipe
    (Path(record["Output"]) / DEVICE_RECORD).write_text(json.dumps(record, indent=2) + "\n")
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("device_config", type=Path)
    parser.add_argument("--palace", type=Path, required=True, help="Palace executable (geometry preflights only)")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, default=PRODUCTION_MANIFEST)
    parser.add_argument("--mesh-recipe", help="repository path of the frozen mesh recipe (default: the one every "
                                              "trace-basis case of the manifest binds)")
    parser.add_argument("--ring-size", type=int, default=DEFAULT_RING_SIZE)
    parser.add_argument("--register", action="store_true", help="register the produced directories into --manifest")
    parser.add_argument("--work", type=Path, help="parent of the registration work directories")
    parser.add_argument("--julia", default=shutil.which("julia"))
    parser.add_argument("--python", default=sys.executable)
    args = parser.parse_args(argv)
    try:
        record = prepare_device_sources(args.device_config, palace=args.palace, output=args.output, manifest_path=args.manifest,
                                        ring_size=args.ring_size, python=args.python)
        if args.register:
            register_device_sources(record, manifest_path=args.manifest, mesh_recipe=args.mesh_recipe, work=args.work,
                                    python=args.python, julia=args.julia)
    except DeviceAdapterError as error:
        print(f"DEVICE_ADAPTER_FAILED: {error}", file=sys.stderr)
        return 1
    print(f"DEVICE {record['Device']['Config']}: {len(record['Coupons'])} spatial coupons, {len(record['OutOfScope'])} out of "
          f"scope; record {args.output / DEVICE_RECORD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
