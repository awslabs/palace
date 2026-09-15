#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Bounded meshing regressions independent of a device name or edge count.

Scouts validate geometry/topology/labels, not response accuracy. Each mesher has
its own time/RSS guard. Analytic metal area/volume gates precede 3D meshing.
"""
import argparse
import csv
import json
import hashlib
import math
import os
import shutil
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def area(p):
    return abs(sum(a[0]*b[1]-b[0]*a[1] for a, b in zip(p, p[1:]+p[:1]))) / 2


def perimeter(p):
    return sum(math.dist(a, b) for a, b in zip(p, p[1:]+p[:1]))


def transform(p, angle=0., shift=(0., 0.)):
    c, s = math.cos(angle), math.sin(angle)
    return [(c*x-s*y+shift[0], s*x+c*y+shift[1]) for x, y in p]


def rectangle(x0, y0, x1, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def write_case(root, name, polygons, *, masks=None, multislot=False, thickness=.08,
               etch=.03, radius=.5, rounding=0., sidewall_angle=90., ma_area=None):
    # polygons: conductor, process plane, normal sign, hole flag, CCW points.
    path = root / name
    path.mkdir(parents=True, exist_ok=False)
    edges, loops, facets = [], [], []
    for loop, (conductor, z, sign, hole, points) in enumerate(polygons, 1):
        for vertex, (a, b) in enumerate(zip(points, points[1:]+points[:1]), 1):
            length = math.dist(a, b)
            tx, ty = (b[0]-a[0])/length, (b[1]-a[1])/length
            gx, gy = (ty, -tx) if not hole else (-ty, tx)
            slot = (vertex-1) % 2 if multislot else 0
            edges.append([len(edges)+1, slot, conductor, (a[0]+b[0])/2,
                          (a[1]+b[1])/2, z, gx, gy, 0, tx, ty, 0, sign,
                          -length/2, length/2, 0])
            loops.append([loop, vertex, conductor, z, int(hole), "Physical", *a])
        if not hole and masks is None:
            facets.append((conductor, z, points))
    if masks is not None:
        facets = masks
    data = [
        ("mesh-signature.csv", "Index,Slot,Conductor,Px,Py,Pz,Gx,Gy,Gz,Tx,Ty,Tz,Nz,S0,S1,VertexArm", edges),
        ("plan-view-boundary.csv", "Loop,Vertex,Conductor,Plane,Hole,Class,X,Y", loops),
        ("plan-view-mask.csv", "Facet,Conductor,Plane,X,Y",
         [[i, c, z, *p] for i, (c, z, ps) in enumerate(facets, 1) for p in ps]),
    ]
    for filename, header, rows in data:
        with (path / filename).open("w", newline="") as stream:
            stream.write(header+"\n")
            csv.writer(stream).writerows(rows)
    (path / "process.toml").write_text(
        f'Units = "um"\nRadius = {radius}\nMetalThickness = {thickness}\n'
        f'Overetch = {etch}\nSidewallAngle = {sidewall_angle}\nTopRounding = {rounding}\n'
        'TrenchRounding = 0.0\n')
    expected = {
        "MetalArea": sum((-1 if h else 1)*area(p) for _, _, _, h, p in polygons),
        "MetalPerimeter": sum(perimeter(p) for _, _, _, _, p in polygons),
        "Thickness": thickness, "Rounded": rounding > 0, "ExpectedMAArea": ma_area,
        "Conductors": sorted({c for c, _, _, _, _ in polygons}),
        "Slots": [0, 1] if multislot else [0], "Layers": len({z for _, z, _, _, _ in polygons}),
    }
    (path / "expected.json").write_text(json.dumps(expected, indent=2)+"\n")
    return path


def fixtures(root):
    l_shape = [(0., 0.), (1.2, 0.), (1.2, .3), (.3, .3), (.3, 1.), (0., 1.)]
    cases = [write_case(root, "strip", [(1, 0., 1, False, rectangle(-.6, -.1, .6, .1))]),
             write_case(root, "concave-multislot", [(1, 0., 1, False, l_shape)], multislot=True),
             write_case(root, "rotated-concave", [(1, 0., 1, False, transform(l_shape, .63, (2., -1.)))], multislot=True),
             write_case(root, "two-conductors", [(1, 0., 1, False, rectangle(-.6, -.5, -.1, .5)),
                                                  (2, 0., 1, False, rectangle(.1, -.5, .6, .5))], multislot=True)]
    outer, hole = rectangle(-.8, -.8, .8, .8), rectangle(-.3, -.3, .3, .3)
    masks = [(1, 0., p) for p in [rectangle(-.8, -.8, .8, -.3), rectangle(-.8, .3, .8, .8),
                                  rectangle(-.8, -.3, -.3, .3), rectangle(.3, -.3, .8, .3)]]
    cases.append(write_case(root, "hole", [(1, 0., 1, False, outer), (1, 0., 1, True, hole)], masks=masks))
    cases.append(write_case(root, "opposed-layers", [(1, 0., 1, False, rectangle(-.6, -.4, .6, .4)),
                                                     (1, .6, -1, False, rectangle(-.4, -.6, .4, .6))],
                            multislot=True, thickness=.06, etch=.02))
    cases.append(write_case(root, "t-junction", [(1,0.,1,False,
        [(-.7,0.),(.7,0.),(.7,.3),(.15,.3),(.15,1.),(-.15,1.),(-.15,.3),(-.7,.3)])],multislot=True))
    cases.append(write_case(root, "three-conductors",
        [(i,0.,1,False,rectangle(x-.12,-.6,x+.12,.6))
         for i,x in enumerate((-.5,0.,.5),1)],multislot=True))
    cases.append(write_case(root, "rounded-strip", [(1, 0., 1, False, rectangle(-.6, -.2, .6, .2))],
                            rounding=.005))
    angle=80.; thickness=.08
    inset=thickness/math.tan(math.radians(angle))
    top_width,top_height=1.2-2*inset,.4-2*inset
    ma_area=top_width*top_height+(1.2+top_width+.4+top_height)*thickness/math.sin(math.radians(angle))
    cases.append(write_case(root,"sloped-strip",[(1,0.,1,False,rectangle(-.6,-.2,.6,.2))],
                            sidewall_angle=angle,thickness=thickness,ma_area=ma_area))
    return cases


def measures(path):
    with path.open() as stream:
        return {(int(r["dimension"]), int(r["attribute"])): float(r["measure"])
                for r in csv.DictReader(stream)}


def geometry_gate(expected, values, kind):
    total = lambda family: sum(v for (d, a), v in values.items() if d == 2 and a//1000 == family)
    checks = {}
    if kind == "thin":
        checks["ThinMetalAreaRelativeError"] = abs(total(4)/expected["MetalArea"]-1)
    else:
        checks["MSAreaRelativeError"] = abs(total(5)/expected["MetalArea"]-1)
        ma = expected["ExpectedMAArea"]
        if ma is None:
            ma = expected["MetalArea"] + expected["MetalPerimeter"]*expected["Thickness"]
        if not expected["Rounded"]:
            checks["MAAreaRelativeError"] = abs(total(6)/ma-1)
        elif abs(total(6)/ma-1)<1e-5:
            raise ValueError("Requested rounding did not change the planar metal area")
    if any(v > 1e-7 for v in checks.values()):
        raise ValueError(f"Independent geometric area gate failed: {checks}")
    return checks


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _same_points(expected, actual, tolerance):
    if len(expected) != len(actual):
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


def audit_manifest_evidence(evidence, gates):
    """Apply geometry-independent mesh gates to normalized tool evidence."""
    failures = []
    labels = evidence.get("ExactLabelsMaterials", {})
    for name in ("VolumeAttributes", "BoundaryAttributes"):
        expected_key, actual_key = "Expected" + name, "Actual" + name
        if (expected_key not in labels or actual_key not in labels or
                sorted(labels.get(expected_key, [])) != sorted(labels.get(actual_key, []))):
            if "exact-labels-materials" not in failures:
                failures.append("exact-labels-materials")

    ownership = evidence.get("OwnershipClosure", {})
    if (ownership.get("UnmatchedPolicy") != "Error" or
            ownership.get("Unmatched") != 0 or ownership.get("Overlaps") != 0 or
            ownership.get("Exhaustive") is not True):
        failures.append("ownership-exhaustive-closure")

    corners = evidence.get("SemanticCorners", {})
    if ("Expected" not in corners or "Actual" not in corners or
            not _same_points(corners.get("Expected", []), corners.get("Actual", []),
                             float(gates["CornerTolerance"]))):
        failures.append("semantic-corners")

    protected = evidence.get("ProtectedSurfaces", {})
    if ("Expected" not in protected or "Actual" not in protected or
            "Changed" not in protected or
            sorted(protected.get("Expected", [])) != sorted(protected.get("Actual", [])) or
            protected.get("Changed") != 0):
        failures.append("protected-surfaces")

    widths = evidence.get("AchievedAnisotropy", {})
    transverse = [widths.get("Transverse1P90"), widths.get("Transverse2P90")]
    normal = widths.get("NormalTarget")
    tangent = widths.get("TangentialP50")
    if (not isinstance(widths.get("Samples"), int) or widths.get("Samples", 0) <= 0 or
            any(not isinstance(x, (int, float)) or not math.isfinite(x) for x in
                [*transverse, normal, tangent]) or min(*transverse, normal, tangent) <= 0 or
            max(transverse) > float(gates["MaximumNormalFactor"]) * normal or
            tangent < float(gates["MinimumAchievedAspect"]) * max(transverse)):
        failures.append("achieved-anisotropy")

    covariance = evidence.get("RotationCovariance", {})
    covariance_error = covariance.get("MaximumRelativeInvariantError")
    if (not covariance.get("ComparedVariant") or
            not isinstance(covariance_error, (int, float)) or
            not math.isfinite(covariance_error) or covariance_error < 0 or
            covariance_error > float(gates["RotationTolerance"])):
        failures.append("rotation-covariance")

    diagonal = evidence.get("TraceDiagonal", {})
    if diagonal.get("GlobalDiagonalBands") != 0:
        failures.append("trace-diagonal-overrefinement")

    resources = evidence.get("Resources", {})
    finite_resources = all(isinstance(resources.get(name), (int, float)) and
                           math.isfinite(resources[name])
                           for name in ("Seconds", "PeakRSSGiB", "Elements"))
    if (resources.get("ExitCode") != 0 or not finite_resources or
            resources.get("Seconds", math.inf) > float(gates["MaximumSeconds"]) or
            resources.get("PeakRSSGiB", math.inf) > float(gates["MaximumRSSGiB"]) or
            resources.get("Elements", math.inf) > int(gates["MaximumElements"])):
        failures.append("bounded-resources")
    return failures


def run_manifest(args):
    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("Version") != 1 or not isinstance(manifest.get("Cases"), list):
        raise ValueError("Unsupported generality-suite manifest")
    identifiers = [case.get("Id") for case in manifest["Cases"]]
    if any(not value for value in identifiers) or len(set(identifiers)) != len(identifiers):
        raise ValueError("Manifest case identifiers must be nonempty and unique")
    overrides = {}
    for item in args.input:
        if "=" not in item:
            raise ValueError("--input must be CASE=DIRECTORY")
        key, value = item.split("=", 1)
        if key in overrides:
            raise ValueError("Duplicate input override: " + key)
        overrides[key] = Path(value).resolve()
    unknown = set(overrides) - set(identifiers)
    if unknown:
        raise ValueError("Input overrides name unknown cases: " + ", ".join(sorted(unknown)))

    repository = (manifest_path.parent / manifest["RepositoryRoot"]).resolve()
    records = []
    preflight_ok = True
    for case in manifest["Cases"]:
        record = {"Id": case["Id"], "Passed": False, "Variants": case.get("Variants", [])}
        try:
            source = case["Source"]
            directory = overrides.get(case["Id"])
            if directory is None and source.get("Directory"):
                candidate = Path(source["Directory"])
                directory = candidate if candidate.is_absolute() else repository / candidate
            if directory is None or not directory.is_dir():
                raise ValueError("required immutable input directory is unavailable")
            hashes = {}
            for role, entry in source["Files"].items():
                path = directory / entry["Name"]
                expected = entry.get("SHA256")
                if not expected:
                    raise ValueError(f"{role} has no frozen SHA256")
                if not path.is_file():
                    raise ValueError(f"missing required {role}: {path}")
                actual = sha256(path)
                if actual != expected:
                    raise ValueError(f"immutable {role} hash mismatch")
                hashes[role] = actual
            signature_role = source["SignatureRole"]
            with (directory / source["Files"][signature_role]["Name"]).open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            required_columns = set(source["SignatureColumns"])
            if not rows or not required_columns.issubset(rows[0]):
                raise ValueError("empty or malformed edge signature")
            record.update({"InputDirectory": str(directory), "InputSHA256": hashes,
                           "DiscoveredEdgeCount": len(rows),
                           "DiscoveredSlots": sorted({int(row["Slot"]) for row in rows}),
                           "DiscoveredConductors": sorted({int(row["Conductor"]) for row in rows})})
            record["Passed"] = True
        except (KeyError, OSError, ValueError) as error:
            record["Error"] = str(error)
            preflight_ok = False
        records.append(record)

    summary = {"Version": 1, "Scope": "Mesh-only geometry-independence gates",
               "Manifest": str(manifest_path), "PreflightPassed": preflight_ok,
               "Cases": records, "Passed": False}
    args.root.mkdir(parents=True, exist_ok=False)
    if not preflight_ok:
        (args.root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        return False
    if args.preflight_only:
        summary["Passed"] = True
        (args.root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        return True
    if args.audit_root is None:
        raise ValueError("--audit-root is required unless --preflight-only is used")
    for case, record in zip(manifest["Cases"], records):
        evidence_path = args.audit_root / (case["Id"] + ".json")
        if not evidence_path.is_file():
            record["Passed"] = False
            record["Error"] = "missing mesh audit evidence"
            continue
        evidence = json.loads(evidence_path.read_text())
        failures = audit_manifest_evidence(evidence, manifest["Gates"])
        record["AuditEvidence"] = str(evidence_path)
        record["GateFailures"] = failures
        record["Passed"] = not failures
    summary["Passed"] = all(record["Passed"] for record in records)
    (args.root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary["Passed"]


def run(args):
    if args.manifest:
        return run_manifest(args)
    args.root.mkdir(parents=True, exist_ok=False)
    cases = fixtures(args.root)
    tools=args.root/"tools"
    tools.mkdir()
    sources=["mesh_graded_tet_experiment.jl", "frozen_volume_study.jl", "mesh_spatial_coupon.jl",
             "graded_curve_distance.jl", "graded_size_points.jl", "graded_trace_size.jl",
             "surface_ribbon_constraints.jl", "interface_ownership.jl", "ownership_bernstein.jl",
             "label_interface_patches.jl", "run_bounded_mesher.py",
             "run_general_mesh_suite.py", "audit_mesh_measures.cpp"]
    for name in sources:
        shutil.copy2(HERE/name,tools/name)
    summary = {"Scope": "Geometry/mesh scouts only; no response or library qualification",
               "ToolSHA256": {name: hashlib.sha256((tools/name).read_bytes()).hexdigest() for name in sources},
               "Cases": []}
    for option in ("audit_bin", "measures_bin"):
        binary=getattr(args,option)
        if binary:
            target=tools/option
            shutil.copy2(binary,target)
            setattr(args,option,target)
            summary["ToolSHA256"][option]=hashlib.sha256(target.read_bytes()).hexdigest()
    for case in cases:
        if args.case and case.name not in args.case:
            continue
        expected = json.loads((case / "expected.json").read_text())
        for kind in ("thin", "fabricated"):
            if expected["Rounded"] and kind == "thin":
                continue
            record = {"Case": case.name, "Kind": kind, "Passed": False}
            try:
                env = dict(os.environ, JULIA_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                           TET_GEOMETRY_ORDER="2" if expected["Rounded"] else "1",
                           TET_EDGE_TANGENT_SIZE="0.05", TET_SURFACE_ALGORITHM="5",
                           TET_HXT_QUALITY="0.1", TET_ANISOTROPIC_SURFACE="0")
                common = [args.julia, "--startup-file=no", f"--project={args.julia_project}",
                          str(tools / "mesh_graded_tet_experiment.jl"), str(case), kind]
                geometry = case / f"{kind}-geometry.msh"
                options = ["1.0", "0.02", "0.25", "--process", str(case / "process.toml"),
                           "--element-interface-slots"]
                def bounded(command, log):
                    subprocess.run([sys.executable, str(tools / "run_bounded_mesher.py"),
                                    "--seconds", str(args.seconds), "--memory-gib", "4",
                                    "--log", str(log), "--", *command], env=env, check=True,
                                   timeout=args.seconds+15, stdout=subprocess.DEVNULL)
                bounded(common+[str(geometry), *options, "--geometry-only"], case/f"{kind}-geometry.log")
                reference = Path(str(geometry)+".cad-measures.csv")
                cad_measures=measures(reference)
                record.update(geometry_gate(expected, cad_measures, kind))
                # CAD whole-face labels are not a reference slot partition. Derive
                # the expected label set from the independent signature contract.
                expected_tags=set()
                for dim, attribute in cad_measures:
                    if dim!=2 or attribute==1:
                        continue
                    family=attribute//1000
                    if family in (4,5,6):
                        expected_tags.update(1000*family+100*s+attribute%100 for s in expected["Slots"])
                    elif family==3:
                        expected_tags.update((3100 if attribute>=3100 else 3000)+s for s in expected["Slots"])
                tags_file=case/f"{kind}-expected-interfaces.csv"
                tags_file.write_text("attribute\n"+"".join(f"{a}\n" for a in sorted(expected_tags)))
                mesh = case / f"{kind}.msh"
                bounded(common+[str(mesh), *options, "--reference-measures", str(reference),
                                "--expected-interfaces",str(tags_file)], case/f"{kind}.log")
                cad = measures(reference)
                if args.measures_bin:
                    measured_path=case/f"{kind}-measures.json"
                    subprocess.run([str(args.measures_bin),str(mesh),str(measured_path)],
                                   check=True,timeout=30)
                    integrated=json.loads(measured_path.read_text())
                    for attribute, volume in integrated["MaterialVolumes"].items():
                        if abs(volume/cad[(3,int(attribute))]-1)>1e-6:
                            raise ValueError("Integrated mesh volume differs from CAD")
                    family=lambda a: a if a==1 else (a//1000)*1000+(a%100 if a//1000 in (4,5,6) else (100 if a>=3100 else 0))
                    expected_areas, actual_areas = {}, {}
                    for (dim,attribute),value in cad.items():
                        if dim==2:
                            key=family(attribute)
                            expected_areas[key]=expected_areas.get(key,0.)+value
                    for attribute,value in integrated["BoundaryAreas"].items():
                        key=family(int(attribute))
                        actual_areas[key]=actual_areas.get(key,0.)+value
                    if expected_areas.keys()!=actual_areas.keys():
                        raise ValueError("Mesh and CAD physical families differ")
                    record["MaximumMeshCADAreaRelativeError"]=max(abs(actual_areas[k]/v-1) for k,v in expected_areas.items())
                    record["MeshCADAreaRelativeTolerance"] = 1e-3 if expected["Rounded"] else 1e-6
                    if record["MaximumMeshCADAreaRelativeError"]>record["MeshCADAreaRelativeTolerance"]:
                        raise ValueError("Curved/affine boundary area differs from CAD")
                    record["BoundaryCoverageChecked"]=integrated["BoundaryCoverageChecked"]
                    record["MaximumRelativeQuadratureDifference"]=integrated["MaximumRelativeQuadratureDifference"]
                if args.audit_bin:
                    audit = case / f"{kind}-audit.json"
                    subprocess.run([str(args.audit_bin), str(mesh), str(case/"plan-view-boundary.csv"),
                                    str(int(kind=="fabricated")), str(audit)], check=True, timeout=30,
                                   stdout=subprocess.DEVNULL)
                    result = json.loads(audit.read_text())
                    # This legacy quality audit's distance histograms assume the
                    # original process; only affine Jacobian kappa/counts are used.
                    record["Elements"] = result["Elements"]
                    record["KappaMax"] = result["KappaPercentiles"]["100.000000"]
                with Path(str(mesh)+".interface-partition.csv").open() as stream:
                    rows = list(csv.DictReader(stream))
                record["MaximumSampledAmbiguousAreaFraction"] = max(float(r["ambiguous_fraction"]) for r in rows)
                record["Attributes"] = sorted(int(r["attribute"]) for r in rows)
                record["Passed"] = True
            except (ValueError, KeyError, subprocess.SubprocessError) as error:
                record["Error"] = str(error)
            summary["Cases"].append(record)
            (args.root / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
            print(json.dumps(record), flush=True)
    return all(c["Passed"] for c in summary["Cases"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--julia-project", type=Path)
    parser.add_argument("--audit-bin", type=Path)
    parser.add_argument("--measures-bin", type=Path)
    parser.add_argument("--seconds", type=float, default=45)
    parser.add_argument("--case", action="append")
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--input", action="append", default=[], metavar="CASE=DIRECTORY")
    parser.add_argument("--audit-root", type=Path)
    parser.add_argument("--preflight-only", action="store_true")
    parsed = parser.parse_args()
    if not parsed.manifest and parsed.julia_project is None:
        parser.error("--julia-project is required for the generated fixture suite")
    raise SystemExit(0 if run(parsed) else 1)
