#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Fail-closed manifest and evidence gate for coupon mesh generality."""
import csv
import hashlib
import json
import math
from pathlib import Path

from audit_edge_metric_mesh import analyze
from mesh_array_io import read_mesh
from mesh_stage_contract import STAGE_ORDER, validate_stage_dag
from semantic_mesh_contract import (REQUIRED_ROLES, load_semantic_contract,
                                    validate_feature_topology)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
                                     separators=(",", ":")).encode()).hexdigest()


def _artifact_path(base, item):
    path = Path(item["Path"])
    return path if path.is_absolute() else (base / path).resolve()


def _check_artifact(base, item, description):
    if not isinstance(item, dict) or not item.get("Path") or not item.get("SHA256"):
        raise ValueError(f"{description} binding is incomplete")
    path = _artifact_path(base, item)
    if not path.is_file() or sha256(path) != item["SHA256"]:
        raise ValueError(f"{description} artifact hash mismatch")
    return path


def _validate_mesh(path, contract):
    try:
        analyze(read_mesh(path), contract, require_material_names=True)
    except SystemExit as error:
        raise ValueError("audited mesh is not a readable Gmsh mesh") from error


def _transform_points(points, transform):
    result = []
    for x, y, z in points:
        value = (x, y, z, 1.0)
        result.append([sum(transform[4 * row + column] * value[column]
                           for column in range(4)) for row in range(3)])
    return result


def _same_points(expected, actual, tolerance):
    if not expected or not actual or len(expected) != len(actual):
        return False
    unused = [tuple(float(x) for x in point) for point in actual]
    for point in expected:
        point = tuple(float(x) for x in point)
        match = next((i for i, other in enumerate(unused)
                      if len(point) == len(other) and math.dist(point, other) <= tolerance), None)
        if match is None:
            return False
        unused.pop(match)
    return True


def _finite_number(value, *, nonnegative=False, positive=False):
    good = isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
    return good and (not nonnegative or value >= 0) and (not positive or value > 0)


def _is_rigid_transform(transform, tolerance=1e-12):
    if transform[12:] != [0, 0, 0, 1] and transform[12:] != [0.0, 0.0, 0.0, 1.0]:
        return False
    rotation = [[transform[4 * row + column] for column in range(3)]
                for row in range(3)]
    gram = [[sum(rotation[k][i] * rotation[k][j] for k in range(3))
             for j in range(3)] for i in range(3)]
    determinant = (rotation[0][0] * (rotation[1][1] * rotation[2][2] -
                                     rotation[1][2] * rotation[2][1]) -
                   rotation[0][1] * (rotation[1][0] * rotation[2][2] -
                                     rotation[1][2] * rotation[2][0]) +
                   rotation[0][2] * (rotation[1][0] * rotation[2][1] -
                                     rotation[1][1] * rotation[2][0]))
    return (all(abs(gram[i][j] - float(i == j)) <= tolerance
                for i in range(3) for j in range(3)) and
            abs(determinant - 1.0) <= tolerance)


def validate_manifest(manifest, manifest_path, *, check_available_files=True):
    if manifest.get("Version") != 2 or not isinstance(manifest.get("Cases"), list):
        raise ValueError("Unsupported generality-suite manifest")
    if not manifest["Cases"]:
        raise ValueError("Manifest must declare cases")
    identifiers = [case.get("Id") for case in manifest["Cases"]]
    if any(not value for value in identifiers) or len(set(identifiers)) != len(identifiers):
        raise ValueError("Manifest case identifiers must be nonempty and unique")
    gates = manifest.get("Gates", {})
    required_gates = ("CornerTolerance", "MaximumNormalFactor", "MinimumAchievedAspect",
                      "MaximumCornerAspect", "MinimumNoncornerAspect",
                      "MaximumProtectedMeasureError", "MinimumScaledJacobian",
                      "MaximumJacobianCondition", "MaximumSeconds", "MaximumRSSGiB",
                      "MaximumElements")
    if any(not _finite_number(gates.get(name), nonnegative=True) for name in required_gates):
        raise ValueError("Manifest has missing or invalid mesh gates")
    repository = (manifest_path.parent / manifest["RepositoryRoot"]).resolve()
    tools = manifest.get("Tools")
    if not isinstance(tools, list) or not tools:
        raise ValueError("Manifest must freeze at least one evidence tool")
    tool_hashes = {}
    for tool in tools:
        name, digest = tool.get("Name"), tool.get("SHA256")
        if not name or name in tool_hashes or not digest:
            raise ValueError("Tool names and hashes must be nonempty and unique")
        tool_hashes[name] = digest
        if check_available_files:
            path = repository / tool["Path"]
            if not path.is_file() or sha256(path) != digest:
                raise ValueError(f"frozen tool hash mismatch: {name}")
    stage_tools = manifest.get("StageToolSHA256")
    from mesh_stage_contract import STAGE_TOOLS
    if (not isinstance(stage_tools, dict) or set(stage_tools) != set(STAGE_ORDER) or
            any(not isinstance(stage_tools[stage], dict) or
                set(stage_tools[stage]) != STAGE_TOOLS[stage] or
                any(not isinstance(value, str) or len(value) != 64
                    for value in stage_tools[stage].values())
                for stage in STAGE_ORDER)):
        raise ValueError("Manifest must freeze every stage tool digest")
    matrix = set()
    comparison_kinds = set()
    for case in manifest["Cases"]:
        source = case.get("Source", {})
        files = source.get("Files")
        if not isinstance(files, dict) or any(role not in files for role in REQUIRED_ROLES):
            raise ValueError(f"{case['Id']} lacks a required immutable role")
        if "MeshRecipe" not in files:
            raise ValueError(f"{case['Id']} lacks a frozen mesh recipe")
        variants = case.get("Variants")
        if not isinstance(variants, list) or not variants:
            raise ValueError(f"{case['Id']} has no variants")
        variant_ids = []
        for variant in variants:
            if (not isinstance(variant, dict) or not variant.get("Id") or
                    not isinstance(variant.get("Transform"), list) or
                    len(variant["Transform"]) != 16 or
                    not all(_finite_number(x) for x in variant["Transform"]) or
                    not _is_rigid_transform(variant["Transform"])):
                raise ValueError(f"{case['Id']} has an invalid variant transform")
            variant_ids.append(variant["Id"])
            matrix.add((case["Id"], variant["Id"]))
        if len(set(variant_ids)) != len(variant_ids) or "identity" not in variant_ids:
            raise ValueError(f"{case['Id']} variants must be unique and include identity")
        by_id = {variant["Id"]: variant["Transform"] for variant in variants}
        identity = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
        if by_id["identity"] != identity:
            raise ValueError(f"{case['Id']} identity variant is not the identity transform")
        comparison = case.get("TransformComparison")
        physical_tolerances = (
            "MaximumRelativeVolumeError", "MaximumRelativeSurfaceMeasureError",
            "MaximumProtectedSupportHausdorff", "MaximumProtectedMeasureError",
            "MaximumQualityDistributionRelativeError", "MaximumAnisotropyRelativeError",
            "MaximumComplexityRatio")
        if (not isinstance(comparison, dict) or comparison.get("Reference") not in variant_ids or
                comparison.get("Transformed") not in variant_ids or
                comparison.get("Reference") == comparison.get("Transformed") or
                by_id.get(comparison.get("Reference")) ==
                by_id.get(comparison.get("Transformed")) or
                any(not _finite_number(comparison.get(name), nonnegative=True)
                    for name in physical_tolerances) or
                comparison.get("MaximumComplexityRatio", 0) < 1.0):
            raise ValueError(f"{case['Id']} has no frozen physical covariance comparison")
    comparisons = manifest.get("ScalingComparisons")
    if not isinstance(comparisons, list) or not comparisons:
        raise ValueError("Manifest must declare scaling comparisons")
    for comparison in comparisons:
        kind = comparison.get("Kind")
        refs = (tuple(comparison.get("Reference", [])), tuple(comparison.get("Compared", [])))
        if kind not in ("feature-scaling", "cad-subdivision-sensitivity") or any(ref not in matrix for ref in refs):
            raise ValueError("Invalid scaling comparison")
        if not _finite_number(comparison.get("MaximumNormalizedDOFRatio"), positive=True):
            raise ValueError("Scaling comparison has no positive bound")
        comparison_kinds.add(kind)
    if comparison_kinds != {"feature-scaling", "cad-subdivision-sensitivity"}:
        raise ValueError("Both feature and CAD-subdivision scaling controls are required")
    return repository, tool_hashes, matrix


def _recompute_mesh_measurements(mesh_path, binding, source_paths, bounded):
    """Rerun the frozen producer implementation instead of trusting record JSON."""
    from general_mesh_audit_producer import (bounded_record, complexity_record,
                                             invariants_record, topology_record)
    base = {"Transform": binding["Transform"]}
    # Recover reference and ownership paths from the independently validated
    # embedded stage bindings.
    reference = bounded["seed-generation"]["Artifacts"]["seed-mesh"]["Path"]
    ownership = bounded["final-gmsh-publication"]["Artifacts"][
        "ownership-partition"]["Path"]
    topology = topology_record(dict(base), mesh_path, source_paths["SemanticContract"],
        source_paths["MeshRecipe"], source_paths["Process"], source_paths["Signature"],
        reference, ownership)["Measurements"]
    complexity = complexity_record(dict(base), mesh_path, source_paths["SemanticContract"],
                                   source_paths["MeshRecipe"])["Measurements"]
    invariants = invariants_record(dict(base), mesh_path)["Measurements"]
    return {**topology, **complexity, **invariants}


def _validate_source_transformation(reports, binding, source_paths):
    """Independently validate source-transform inputs, outputs, and metric linkage."""
    from transform_coupon_source_contract import (
        transform_semantic_contract, transformed_supports, validate_rigid_transform)
    stage = reports["source-transformation"]
    role_names = {"source-semantic-contract": "SemanticContract",
                  "source-signature": "Signature", "source-boundary": "Boundary",
                  "source-mask": "Mask"}
    for stage_name, source_role in role_names.items():
        if stage["Inputs"][stage_name]["SHA256"] != binding["InputSHA256"][source_role]:
            raise ValueError("source-transform input differs from immutable source")
    transform_path = Path(stage["Inputs"]["canonical-transform"]["Path"])
    transform = json.loads(transform_path.read_text())
    if isinstance(transform, dict):
        transform = transform.get("Transform")
    if transform != binding["Transform"]:
        raise ValueError("source-transform canonical transform differs from variant")
    matrix = validate_rigid_transform(transform)
    semantic_path = Path(stage["Artifacts"]["transformed-semantic-contract"]["Path"])
    supports_path = Path(stage["Artifacts"]["transformed-supports"]["Path"])
    source_semantic = json.loads(Path(source_paths["SemanticContract"]).read_text())
    expected_semantic = transform_semantic_contract(source_semantic, matrix)
    expected_semantic["SourceSemanticContractSHA256"] = sha256(
        source_paths["SemanticContract"])
    expected_supports = transformed_supports(
        Path(source_paths["Signature"]).parent, matrix,
        signature=source_paths["Signature"], boundary=source_paths["Boundary"],
        mask=source_paths["Mask"])
    if (json.loads(semantic_path.read_text()) != expected_semantic or
            json.loads(supports_path.read_text()) != expected_supports):
        raise ValueError("transformed semantic/support artifact differs from source transform")


def _validate_bound_records(evidence_path, evidence, binding, source_paths):
    records = evidence.get("AuditRecords")
    required = {"bounded-run", "mesh-topology-quality", "mesh-complexity",
                "mesh-invariants", "variant-transform"}
    if (not isinstance(records, list) or len(records) != len(required) or
            {item.get("Kind") for item in records} != required):
        raise ValueError("exactly one record of every required audit kind is required")
    measurements = {}
    digests = {evidence["Mesh"]["SHA256"]}
    expected_dependencies = {
        "bounded-run": {},
        "mesh-topology-quality": {
            role: binding["InputSHA256"][role]
            for role in ("SemanticContract", "MeshRecipe", "Process", "Signature")},
        "mesh-complexity": {
            role: binding["InputSHA256"][role]
            for role in ("MeshRecipe", "SemanticContract")},
        "mesh-invariants": {},
        "variant-transform": {},
    }
    records_by_kind = {}
    for item in records:
        path = _check_artifact(evidence_path.parent, item, "audit record")
        digest = item["SHA256"]
        if digest in digests:
            raise ValueError("audit records must be content-distinct")
        digests.add(digest)
        record = json.loads(path.read_text())
        records_by_kind[item["Kind"]] = record
        expected = {"Version": 1, "Kind": item["Kind"],
                    "CaseId": binding["CaseId"], "Variant": binding["Variant"],
                    "MeshSHA256": evidence["Mesh"]["SHA256"],
                    "InputSHA256": binding["InputSHA256"],
                    "Transform": binding["Transform"],
                    "TransformSHA256": binding["TransformSHA256"]}
        if any(record.get(key) != value for key, value in expected.items()):
            raise ValueError("audit record has stale or mismatched bindings")
        producer = record.get("Producer", {})
        if (producer.get("Name") not in binding["ToolSHA256"] or
                binding["ToolSHA256"][producer["Name"]] != producer.get("SHA256") or
                not record.get("Command") or not isinstance(record.get("Environment"), dict)):
            raise ValueError("audit record producer/command/environment is not frozen")
        for section, value in record.get("Measurements", {}).items():
            if section in measurements:
                raise ValueError("measurement section has multiple producers")
            measurements[section] = value
        dependencies = record.get("Dependencies", {})
        if dependencies != expected_dependencies[item["Kind"]]:
            raise ValueError("audit record dependencies differ from frozen inputs")
        if item["Kind"] == "variant-transform":
            if (record.get("TransformVerified") is not True or
                    evidence.get("IdentityMeshSHA256") != record.get("IdentityMeshSHA256") or
                    evidence.get("IdentitySeedMeshSHA256") !=
                    record.get("IdentitySeedMeshSHA256") or
                    evidence.get("TransformMaximumCoordinateError") !=
                    record.get("TransformMaximumCoordinateError")):
                raise ValueError("variant transform result differs from its bound audit")
        if item["Kind"] == "bounded-run":
            stage_items = record.get("BoundedStageRecords")
            if (not isinstance(stage_items, list) or len(stage_items) != len(STAGE_ORDER) or
                    {stage.get("Stage") for stage in stage_items} != set(STAGE_ORDER)):
                raise ValueError("bounded stage records are incomplete")
            reports, stage_digests = validate_stage_dag(
                {stage["Stage"]: Path(stage["Path"]) for stage in stage_items},
                _artifact_path(evidence_path.parent, evidence["Mesh"]),
                "run_bounded_mesher.py", binding["ToolSHA256"]["run_bounded_mesher.py"],
                binding["StageToolSHA256"])
            if (reports != record.get("BoundedStages") or
                    sorted(stage_digests) != record.get("StageRecordSHA256")):
                raise ValueError("bounded stage DAG differs from its bound producer output")
            if digests & stage_digests:
                raise ValueError("bounded stage and audit artifacts must be content-distinct")
            digests.update(stage_digests)
    bounded = records_by_kind["bounded-run"]["BoundedStages"]
    _validate_source_transformation(bounded, binding, source_paths)
    topology = records_by_kind["mesh-topology-quality"]
    if (topology.get("ReferenceMeshSHA256") !=
            bounded["seed-generation"]["Artifacts"]["seed-mesh"]["SHA256"] or
            topology.get("OwnershipReportSHA256") !=
            bounded["final-gmsh-publication"]["Artifacts"]["ownership-partition"]["SHA256"]):
        raise ValueError("topology audit does not bind staged reference/ownership artifacts")
    mesh_path = _artifact_path(evidence_path.parent, evidence["Mesh"])
    recomputed = _recompute_mesh_measurements(mesh_path, binding, source_paths, bounded)
    recorded_mesh_measurements = {
        key: value for key, value in measurements.items()
        if key not in {"Resources", "PhysicalCovariance"}
    }
    # Variant-transform has no measurement section. Resources are independently
    # derived from the validated stage reports and parsed final tetrahedra below.
    if recomputed != recorded_mesh_measurements:
        raise ValueError("audit measurements differ from an independent producer rerun")
    variant = records_by_kind["variant-transform"]
    identity_path = Path(variant.get("IdentityMeshPath", ""))
    identity_seed_path = Path(variant.get("IdentitySeedMeshPath", ""))
    if (not identity_path.is_file() or sha256(identity_path) !=
            variant.get("IdentityMeshSHA256") or not identity_seed_path.is_file() or
            sha256(identity_seed_path) != variant.get("IdentitySeedMeshSHA256")):
        raise ValueError("variant identity mesh bindings changed")
    from general_mesh_audit_producer import _physical_covariance_report
    matrix = __import__("numpy").asarray(binding["Transform"], dtype=float).reshape(4, 4)
    recomputed_physical = _physical_covariance_report(
        read_mesh(identity_path), read_mesh(mesh_path),
        load_semantic_contract(source_paths["SemanticContract"]), matrix)
    if measurements.get("PhysicalCovariance") != recomputed_physical:
        raise ValueError("physical covariance differs from independent normalization")
    resources = measurements.get("Resources", {})
    expected_resources = {
        "ExitCode": 0,
        "Seconds": sum(report["Seconds"] for report in bounded.values()),
        "PeakRSSGiB": max(report["PeakProcessTreeRSSBytes"]
                          for report in bounded.values()) / 2**30,
    }
    if any(resources.get(key) != value for key, value in expected_resources.items()):
        raise ValueError("resource measurements differ from bounded stage reports")

    required_sections = {"Resources", "ActualVolumeMaterials", "ActualBoundaryAttributes",
                         "ActualAdjacency", "OwnershipClosure", "ActualSemanticCorners",
                         "CornerNeighborhoods", "SubdivisionNeighborhoods", "CutNeighborhoods",
                         "ProtectedSurfaces", "AchievedAnisotropy", "TraceDiagonal",
                         "MeshQuality", "Complexity", "ComparisonInvariants",
                         "PhysicalCovariance"}
    if set(measurements) != required_sections:
        raise ValueError("bound records do not supply the exact measurement schema")
    for section, value in measurements.items():
        if evidence.get(section) != value:
            raise ValueError("normalized measurements differ from bound producer records")
    return digests - {evidence["Mesh"]["SHA256"]}


def audit_manifest_evidence(evidence, gates, contract, binding):
    """Judge actual measurements against a separately frozen semantic contract."""
    failures = []
    required_binding = {
        "CaseId": binding["CaseId"], "Variant": binding["Variant"],
        "TransformSHA256": binding["TransformSHA256"],
        "InputSHA256": binding["InputSHA256"],
        "ProcessSHA256": binding["InputSHA256"]["Process"],
        "SemanticContractSHA256": binding["InputSHA256"]["SemanticContract"],
        "RecipeSHA256": binding["InputSHA256"]["MeshRecipe"],
        "ToolSHA256": binding["ToolSHA256"],
        "StageToolSHA256": binding["StageToolSHA256"],
    }
    if evidence.get("Version") != 3 or any(evidence.get(key) != value
                                            for key, value in required_binding.items()):
        failures.append("provenance-binding")

    expected_materials = sorted(contract["VolumeMaterials"], key=lambda item: item["Attribute"])
    actual_materials = evidence.get("ActualVolumeMaterials")
    expected_labels = sorted(item["Attribute"] for item in contract["BoundaryLabels"])
    actual_labels = evidence.get("ActualBoundaryAttributes")
    if (not actual_materials or not actual_labels or
            sorted(actual_materials, key=lambda item: item.get("Attribute", -1)) != expected_materials or
            sorted(actual_labels) != expected_labels):
        failures.append("exact-labels-materials")
    expected_adjacency = {str(item["Attribute"]): sorted(item["AdjacentMaterials"])
                          for item in contract["BoundaryLabels"]}
    actual_adjacency = evidence.get("ActualAdjacency")
    if (not actual_adjacency or
            {str(key): sorted(value) for key, value in actual_adjacency.items()} !=
            expected_adjacency):
        failures.append("material-adjacency")

    ownership = evidence.get("OwnershipClosure", {})
    if (ownership.get("UnmatchedPolicy") != contract["UnmatchedPolicy"] or
            ownership.get("Unmatched") != 0 or ownership.get("Overlaps") != 0 or
            ownership.get("Exhaustive") is not True):
        failures.append("ownership-exhaustive-closure")
    expected_corners = _transform_points(contract["SemanticCorners"], binding["Transform"])
    if not _same_points(expected_corners, evidence.get("ActualSemanticCorners"),
                        float(gates["CornerTolerance"])):
        failures.append("semantic-corners")
    def neighborhood_failure(name, expected, maximum=None, minimum=None):
        values = evidence.get(name)
        if not isinstance(values, list) or len(values) != len(expected):
            return True
        if not expected:
            return False
        if not _same_points(expected, [item.get("Point") for item in values],
                            float(gates["CornerTolerance"])):
            return True
        aspects = [item.get("MaximumAspect") for item in values]
        return (any(not _finite_number(value, positive=True) for value in aspects) or
                (maximum is not None and any(value > maximum for value in aspects)) or
                (minimum is not None and any(value < minimum for value in aspects)))
    topology = contract["FeatureTopology"]
    subdivisions = _transform_points(topology["CADSubdivisionEndpoints"], binding["Transform"])
    cuts = _transform_points(topology["CutEndpoints"], binding["Transform"])
    if (neighborhood_failure("CornerNeighborhoods", expected_corners,
                             maximum=gates["MaximumCornerAspect"]) or
            neighborhood_failure("SubdivisionNeighborhoods", subdivisions,
                                 minimum=gates["MinimumNoncornerAspect"]) or
            neighborhood_failure("CutNeighborhoods", cuts,
                                 minimum=gates["MinimumNoncornerAspect"])):
        failures.append("semantic-corner-and-endpoint-anisotropy")
    protected = evidence.get("ProtectedSurfaces", {})
    if (not protected.get("Actual") or protected.get("PlaneSupportsMatch") is not True or
            protected.get("TopologyMatches") is not True or
            sorted(protected.get("Actual", [])) != sorted(contract["ProtectedSupports"]) or
            not _finite_number(protected.get("MaximumRelativeMeasureError"),
                               nonnegative=True) or
            protected.get("MaximumRelativeMeasureError", math.inf) >
            gates["MaximumProtectedMeasureError"] or
            not _finite_number(protected.get("MaximumSupportVertexDistance"),
                               nonnegative=True) or
            protected.get("MaximumSupportVertexDistance", math.inf) >
            gates["CornerTolerance"]):
        failures.append("protected-surfaces")

    widths = evidence.get("AchievedAnisotropy", {})
    values = [widths.get(name) for name in ("Transverse1P90", "Transverse2P90",
                                             "NormalTarget", "TangentialP50")]
    if (not isinstance(widths.get("Samples"), int) or widths.get("Samples", 0) <= 0 or
            any(not _finite_number(x, positive=True) for x in values) or
            max(values[:2]) > gates["MaximumNormalFactor"] * values[2] or
            values[3] < gates["MinimumAchievedAspect"] * max(values[:2])):
        failures.append("achieved-anisotropy")
    if evidence.get("TraceDiagonal", {}).get("GlobalDiagonalBands") != 0:
        failures.append("trace-diagonal-overrefinement")

    quality = evidence.get("MeshQuality", {})
    if (not isinstance(quality.get("Samples"), int) or quality.get("Samples", 0) <= 0 or
            quality.get("PositiveOrientation") is not True or
            not _finite_number(quality.get("MinimumScaledJacobian"), nonnegative=True) or
            quality.get("MinimumScaledJacobian", -1) < gates["MinimumScaledJacobian"] or
            not _finite_number(quality.get("MaximumJacobianCondition"), positive=True) or
            quality.get("MaximumJacobianCondition", math.inf) > gates["MaximumJacobianCondition"]):
        failures.append("mesh-quality-jacobian")
    resources = evidence.get("Resources", {})
    resource_names = ("Seconds", "PeakRSSGiB", "Elements")
    if (resources.get("ExitCode") != 0 or
            any(not _finite_number(resources.get(name), nonnegative=True) for name in resource_names) or
            resources.get("Seconds", math.inf) > gates["MaximumSeconds"] or
            resources.get("PeakRSSGiB", math.inf) > gates["MaximumRSSGiB"] or
            resources.get("Elements", math.inf) > gates["MaximumElements"]):
        failures.append("bounded-resources")
    complexity = evidence.get("Complexity", {})
    if (not isinstance(complexity.get("H1DOFs"), int) or complexity.get("H1DOFs", 0) <= 0 or
            not isinstance(complexity.get("FeatureCount"), int) or
            complexity.get("FeatureCount", 0) <= 0 or
            not isinstance(complexity.get("CADSubdivisionCount"), int) or
            complexity.get("CADSubdivisionCount", -1) < 0):
        failures.append("complexity-counts")
    invariants = evidence.get("ComparisonInvariants")
    if (not isinstance(invariants, dict) or not invariants or
            any(not _finite_number(value) for value in invariants.values())):
        failures.append("comparison-invariants")
    physical = evidence.get("PhysicalCovariance", {})
    if (physical.get("ComparisonFrame") != "SourceLocal" or
            physical.get("LabelsMaterialsAdjacencyMatch") is not True):
        failures.append("physical-covariance-contract")
    return failures


def _relative_error(a, b):
    scale = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / scale


def _physical_comparison_failures(reference_evidence, transformed_evidence, comparison):
    """Compare final meshes physically; deterministic topology is diagnostic only."""
    failures = []
    physical = transformed_evidence.get("PhysicalCovariance", {})
    if physical.get("LabelsMaterialsAdjacencyMatch") is not True:
        failures.append("physical labels/material adjacency")
        return failures
    left, right = (physical.get(name, {}) for name in
                   ("ReferenceInvariants", "TransformedInvariants"))
    if left.keys() != right.keys() or not left:
        return ["physical measure keys differ"]
    volume_error = max((_relative_error(left[name], right[name]) for name in left
                        if name.startswith("Volume:")), default=math.inf)
    area_error = max((_relative_error(left[name], right[name]) for name in left
                      if name.startswith("Area:")), default=math.inf)
    if volume_error > comparison["MaximumRelativeVolumeError"]:
        failures.append("material volumes")
    if area_error > comparison["MaximumRelativeSurfaceMeasureError"]:
        failures.append("boundary surface measures")
    protected = physical.get("ProtectedSurfaces", {})
    if (protected.get("PlaneSupportsMatch") is not True or
            protected.get("TopologyMatches") is not True or
            protected.get("MaximumSupportVertexDistance", math.inf) >
            comparison["MaximumProtectedSupportHausdorff"] or
            protected.get("MaximumRelativeMeasureError", math.inf) >
            comparison["MaximumProtectedMeasureError"]):
        failures.append("protected support geometry/topology")
    qualities = [physical.get(name, {}) for name in
                 ("ReferenceQuality", "TransformedQuality")]
    quality_values = []
    for name in ("ScaledJacobianQuantiles", "JacobianConditionQuantiles"):
        if (not all(isinstance(item.get(name), list) for item in qualities) or
                len(qualities[0].get(name, [])) != len(qualities[1].get(name, []))):
            quality_values = [math.inf]
            break
        quality_values.extend(_relative_error(a, b)
                              for a, b in zip(qualities[0][name], qualities[1][name]))
    if (not all(item.get("PositiveOrientation") is True for item in qualities) or
            max(quality_values, default=math.inf) >
            comparison["MaximumQualityDistributionRelativeError"]):
        failures.append("orientation/quality distribution")
    anisotropy_names = ("TangentialP50", "Transverse1P90", "Transverse2P90")
    anisotropy = [item.get("AchievedAnisotropy", {})
                  for item in (reference_evidence, transformed_evidence)]
    anisotropy_error = max((_relative_error(anisotropy[0].get(name, math.inf),
                                            anisotropy[1].get(name, -math.inf))
                            for name in anisotropy_names), default=math.inf)
    if anisotropy_error > comparison["MaximumAnisotropyRelativeError"]:
        failures.append("local-frame anisotropy")
    complexity_values = [
        [item["Complexity"]["H1DOFs"], item["Resources"]["Elements"]]
        for item in (reference_evidence, transformed_evidence)]
    complexity_ratio = max(max(a, b) / max(min(a, b), 1)
                           for a, b in zip(*complexity_values))
    if complexity_ratio > comparison["MaximumComplexityRatio"]:
        failures.append("complexity ratio")
    return failures


def run_manifest(args):
    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text())
    repository, tool_hashes, _ = validate_manifest(manifest, manifest_path)
    cases_by_id = {case["Id"]: case for case in manifest["Cases"]}
    overrides = {}
    for item in args.input:
        if "=" not in item:
            raise ValueError("--input must be CASE=DIRECTORY")
        key, value = item.split("=", 1)
        if key in overrides or key not in cases_by_id:
            raise ValueError("Duplicate or unknown input override: " + key)
        overrides[key] = Path(value).resolve()

    records, sources, preflight_ok = [], {}, True
    for case in manifest["Cases"]:
        record = {"Id": case["Id"], "Passed": False,
                  "Variants": [variant["Id"] for variant in case["Variants"]]}
        try:
            source = case["Source"]
            directory = overrides.get(case["Id"])
            if directory is None and source.get("Directory"):
                candidate = Path(source["Directory"])
                directory = candidate if candidate.is_absolute() else repository / candidate
            if directory is None or not directory.is_dir():
                raise ValueError("required immutable input directory is unavailable")
            hashes, paths = {}, {}
            for role, entry in source["Files"].items():
                expected, name = entry.get("SHA256"), entry.get("Name")
                repository_name = entry.get("RepositoryPath")
                if not expected or (not name and not repository_name) or (name and repository_name):
                    raise ValueError(f"{role} has no unambiguous frozen path and SHA256")
                candidate = Path(repository_name or name)
                if candidate.is_absolute():
                    path = candidate
                elif repository_name:
                    path = repository / candidate
                else:
                    path = directory / candidate
                if not path.is_file() or sha256(path) != expected:
                    raise ValueError(f"immutable {role} hash mismatch")
                hashes[role], paths[role] = expected, path
            contract = load_semantic_contract(paths["SemanticContract"])
            validate_feature_topology(contract, paths["Signature"], paths["Boundary"])
            signature_role = source["SignatureRole"]
            with paths[signature_role].open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            required_columns = set(source["SignatureColumns"])
            if not rows or not required_columns.issubset(rows[0]):
                raise ValueError("empty or malformed edge signature")
            record.update({"InputDirectory": str(directory), "InputSHA256": hashes,
                           "DiscoveredEdgeCount": len(rows),
                           "DiscoveredSlots": sorted({int(row["Slot"]) for row in rows}),
                           "DiscoveredConductors": sorted({int(row["Conductor"]) for row in rows})})
            sources[case["Id"]] = (hashes, paths, contract)
            record["Passed"] = True
        except (KeyError, OSError, TypeError, ValueError) as error:
            record["Error"] = str(error)
            preflight_ok = False
        records.append(record)

    summary = {"Version": 2, "Scope": "Mesh-only geometry-independence gates",
               "Manifest": str(manifest_path), "PreflightPassed": preflight_ok,
               "Cases": records, "Passed": False}
    args.root.mkdir(parents=True, exist_ok=False)
    if not preflight_ok or args.preflight_only:
        summary["Passed"] = preflight_ok and args.preflight_only
        (args.root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        return summary["Passed"]
    if args.audit_root is None:
        raise ValueError("--audit-root is required unless --preflight-only is used")

    evidence_by_key, used_audits, used_meshes = {}, set(), set()
    for case, record in zip(manifest["Cases"], records):
        record["VariantResults"] = []
        hashes, _, contract = sources[case["Id"]]
        for variant in case["Variants"]:
            variant_id = variant["Id"]
            result = {"Variant": variant_id, "Passed": False}
            path = args.audit_root / f"{case['Id']}--{variant_id}.json"
            try:
                evidence = json.loads(path.read_text())
                binding = {"CaseId": case["Id"], "Variant": variant_id,
                           "Transform": variant["Transform"],
                           "TransformSHA256": canonical_sha256(variant["Transform"]),
                           "InputSHA256": hashes, "ToolSHA256": tool_hashes,
                           "StageToolSHA256": manifest["StageToolSHA256"]}
                failures = audit_manifest_evidence(evidence, manifest["Gates"], contract, binding)
                mesh_path = _check_artifact(path.parent, evidence.get("Mesh"), "audited mesh")
                _validate_mesh(mesh_path, contract)
                mesh_digest = evidence["Mesh"]["SHA256"]
                if mesh_digest in used_meshes:
                    raise ValueError("audited meshes must be content-distinct per matrix entry")
                used_meshes.add(mesh_digest)
                record_digests = _validate_bound_records(path, evidence, binding, sources[case["Id"]][1])
                if used_audits & record_digests:
                    raise ValueError("audit records must be content-distinct per matrix entry")
                used_audits.update(record_digests)
                result.update({"AuditEvidence": str(path), "GateFailures": failures,
                               "Passed": not failures})
                evidence_by_key[(case["Id"], variant_id)] = evidence
            except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError) as error:
                result["Error"] = str(error)
            record["VariantResults"].append(result)
        record["Passed"] = all(item["Passed"] for item in record["VariantResults"])

    comparison_failures = []
    for case in manifest["Cases"]:
        comparison = case["TransformComparison"]
        keys = [(case["Id"], comparison[name]) for name in ("Reference", "Transformed")]
        if any(key not in evidence_by_key for key in keys):
            comparison_failures.append(case["Id"] + ": missing transform evidence")
            continue
        reference_evidence, transformed_evidence = (evidence_by_key[key] for key in keys)
        identity_digest = reference_evidence["Mesh"]["SHA256"]
        identity_seed_digest = reference_evidence.get("IdentitySeedMeshSHA256")
        coordinate_error = transformed_evidence.get("TransformMaximumCoordinateError")
        if (reference_evidence.get("IdentityMeshSHA256") != identity_digest or
                transformed_evidence.get("IdentityMeshSHA256") != identity_digest or
                transformed_evidence.get("IdentitySeedMeshSHA256") != identity_seed_digest or
                not _finite_number(coordinate_error, nonnegative=True) or
                coordinate_error > manifest["Gates"]["CornerTolerance"]):
            comparison_failures.append(case["Id"] + ": exact source-seed covariance")
            continue
        comparison_failures.extend(
            case["Id"] + ": " + failure for failure in _physical_comparison_failures(
                reference_evidence, transformed_evidence, comparison))
    scaling_failures = []
    for comparison in manifest["ScalingComparisons"]:
        keys = [tuple(comparison[name]) for name in ("Reference", "Compared")]
        if any(key not in evidence_by_key for key in keys):
            scaling_failures.append(comparison["Id"] + ": missing evidence")
            continue
        complexity = [evidence_by_key[key]["Complexity"] for key in keys]
        normalized = [item["H1DOFs"] / item["FeatureCount"] for item in complexity]
        ratio = max(normalized) / min(normalized)
        if (comparison["Kind"] == "cad-subdivision-sensitivity" and
                (complexity[0]["FeatureCount"] != complexity[1]["FeatureCount"] or
                 complexity[0]["CADSubdivisionCount"] == complexity[1]["CADSubdivisionCount"])):
            scaling_failures.append(comparison["Id"] + ": invalid CAD subdivision control")
        elif ratio > comparison["MaximumNormalizedDOFRatio"]:
            scaling_failures.append(comparison["Id"] + ": H1 DOF scaling")
    summary["TransformComparisonFailures"] = comparison_failures
    summary["ScalingComparisonFailures"] = scaling_failures
    summary["Passed"] = (all(record["Passed"] for record in records) and
                         not comparison_failures and not scaling_failures)
    (args.root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary["Passed"]
