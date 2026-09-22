#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Compare the records of two builds of one coupon case field by field (decision 62 step
3 acceptance: a dedupe / caching change must leave every record bit-identical apart from
timings, resources, tool digests, paths and commands).

  compare_case_records.py BEFORE_CASE_DIR AFTER_CASE_DIR [--output report.json]

Every JSON record the case directory holds (stage reports, audit records, receipts,
censuses, normalized evidence under audits/, the verification report) is compared with
its namesake; the fields listed in IGNORED_KEYS (wall seconds, RSS, SHA-256 of tools and
of records that embed them, absolute paths, commands, environments, timestamps) are
skipped wherever they appear, and every other difference is listed with its JSON path.
The meshes and the ownership CSVs are compared byte for byte.  Exit status 0 when nothing
else differs."""
import argparse
import json
from pathlib import Path
import sys

IGNORED_KEYS = {
    # timings and resources
    "Seconds", "PeakProcessTreeRSSBytes", "PeakRSSGiB", "WallClockSeconds", "WallSeconds", "StartedAt", "FinishedAt",
    "Timestamp", "Resources", "CanonicalBuild", "PlacementPublication", "PeakGiB", "Elapsed", "Started", "Finished",
    "Timings",
    # tool identities and the digests of records that embed them
    "SHA256", "ToolSHA256", "StageToolSHA256", "Producer", "Tool", "ToolCommit", "Commit", "CanonicalBuildSHA256",
    "CanonicalBuildId", "CanonicalBuildIds", "SharedCanonicalBuildId", "CanonicalToolSHA256", "OwnershipAuditorSHA256", "OwnershipRuntimeSHA256", "ReceiptPayloadSHA256", "TransformReceiptSHA256",
    "StageRecordSHA256", "BoundedStageRecords", "BoundedStages", "VariantRecordSHA256", "CanonicalStageReportSHA256",
    "CanonicalArtifactSHA256", "AuditRecords", "Evidence", "Tools", "Runtime", "Mesher", "Publisher", "Stager",
    "CanonicalToolHashes", "Manifest", "ManifestSHA256", "BuildCensusSHA256", "ReferenceMeshSHA256",
    "OwnershipReportSHA256", "OwnershipQuadratureSHA256", "TransformedSemanticSHA256", "TransformedSupportsSHA256",
    "IdentityMeshSHA256", "IdentitySeedMeshSHA256", "MeshSHA256", "OutputMeshSHA256", "ParentLabeledMeshSHA256",
    "CanonicalMeshSHA256", "RadialShells", "Census", "InputSHA256", "SourceInputSHA256", "TransformSHA256",
    # paths, commands, environments
    "Path", "Root", "Directory", "Command", "OwnershipCommand", "Environment", "WorkingDirectory", "Inputs", "Artifacts",
    "IdentityMeshPath", "IdentitySeedMeshPath", "Log", "Stdout", "Stderr", "LaunchStdout", "LaunchStderr", "Remote",
    "Mesh", "Output", "Outputs", "Local", "Record", "Report", "Records", "AuditEvidence", "SemanticContract",
}
# Records whose content is the tooling identity or the resource accounting itself.
IGNORED_FILES = {"canonical-tool-hashes.json", "input-hashes.json"}
BYTE_COMPARED_SUFFIXES = (".msh", ".csv")


def walk(left, right, path, differences):
    if type(left) is not type(right):
        differences.append({"Path": path, "Before": left, "After": right})
        return
    if isinstance(left, dict):
        for key in sorted(set(left) | set(right)):
            if key in IGNORED_KEYS:
                continue
            if key not in left or key not in right:
                differences.append({"Path": f"{path}/{key}", "Before": left.get(key, "<absent>"),
                                    "After": right.get(key, "<absent>")})
                continue
            walk(left[key], right[key], f"{path}/{key}", differences)
    elif isinstance(left, list):
        if len(left) != len(right):
            differences.append({"Path": path, "Before": f"{len(left)} items", "After": f"{len(right)} items"})
            return
        for index, (a, b) in enumerate(zip(left, right)):
            walk(a, b, f"{path}[{index}]", differences)
    elif left != right and not (isinstance(left, float) and left != left and right != right):
        differences.append({"Path": path, "Before": left, "After": right})


def compare(before, after):
    before, after = Path(before), Path(after)
    report = {"Before": str(before), "After": str(after), "Compared": [], "Missing": [], "Differences": [],
              "ByteIdentical": [], "ByteDifferent": []}
    for path in sorted(p for p in before.rglob("*") if p.is_file()):
        relative = path.relative_to(before)
        other = after / relative
        if not other.is_file():
            report["Missing"].append(str(relative))
            continue
        if path.suffix == ".json" and path.name not in IGNORED_FILES:
            try:
                left, right = json.loads(path.read_text()), json.loads(other.read_text())
            except ValueError:
                continue
            differences = []
            walk(left, right, str(relative), differences)
            report["Compared"].append(str(relative))
            report["Differences"].extend(differences)
        elif path.suffix in BYTE_COMPARED_SUFFIXES:
            (report["ByteIdentical"] if path.read_bytes() == other.read_bytes() else report["ByteDifferent"]).append(str(relative))
    report["Identical"] = not report["Differences"] and not report["ByteDifferent"]
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("before", type=Path)
    parser.add_argument("after", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    report = compare(args.before, args.after)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"records compared {len(report['Compared'])}, byte-identical files {len(report['ByteIdentical'])}, "
          f"byte-different {len(report['ByteDifferent'])}, missing {len(report['Missing'])}, "
          f"differences {len(report['Differences'])}: {'IDENTICAL' if report['Identical'] else 'DIFFERENT'}")
    for item in report["Differences"][:50]:
        print(f"  {item['Path']}: {item['Before']!r} -> {item['After']!r}")
    for name in report["ByteDifferent"][:20]:
        print(f"  byte-different: {name}")
    return 0 if report["Identical"] else 1


if __name__ == "__main__":
    sys.exit(main())
