#!/usr/bin/env python3

"""Combine response models into one portable fabrication-process library.

--supersede NAME=<record> (decision 505 Q1) replaces the base (first) library's model NAME by the
model of the same name of one later input, recording the replacement under Supersedes with the
mandatory --supersede-decision TEXT (decision 513 (2)); without the flag a repeated model name is
refused and the output is unchanged. Two models of one canonical Signature (whatever their Names:
the one-reference and the -pieces form of one key) are refused by name before anything is written
(OI-LIBRARY-DUPLICATE-SIGNATURE, decisions 1103 / 1126 (5); signature_library.refuse_duplicate_signatures).
--replace-by-signature NAME=<record> (idfix R0-b, decisions 1103 / 1136 (1) / 1140 (4)): the model NAME of one
later input (its piece form, Name <base>-pieces<P>) REPLACES the base library's model carrying the SAME canonical
Signature, in place (the base slot), under the same mandatory --supersede-decision TEXT; recorded under Supersedes
with Kind "Signature" - a library never holds two models of one signature."""

import argparse
import hashlib
import json
import math
import re
import shutil
import sys
from pathlib import Path

IDENTIFICATION_TOOLS = Path(__file__).resolve().parents[2] / "surface_response_identification"
if str(IDENTIFICATION_TOOLS) not in sys.path:
    sys.path.insert(0, str(IDENTIFICATION_TOOLS))
from signature_library import refuse_duplicate_signatures, signature_hash  # noqa: E402


PATH_FIELDS = (
    "FabricatedMatrix",
    "ThinMatrix",
    "FabricatedSurfaceMatrix",
    "ThinSurfaceMatrix",
    "BasisPoints",
)

PATH_NAMES = {
    "FabricatedMatrix": "fabricated-domain-response-matrix.csv",
    "ThinMatrix": "thin-domain-response-matrix.csv",
    "FabricatedSurfaceMatrix": "fabricated-surface-response-matrix.csv",
    "ThinSurfaceMatrix": "thin-surface-response-matrix.csv",
    "BasisPoints": "basis-points.csv",
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_record(path, library):
    """Provenance of one input library: its path, digest, name and model names."""
    return {"Path": str(path), "SHA256": sha256(path), "Name": library.get("Name"),
            "Version": library["Version"], "Models": [model.get("Name") for model in library["Models"]]}


def model_directory(index, name):
    slug = re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-")
    return Path("models") / f"{index:03d}-{slug}"


def load_library(path):
    path = path.expanduser().resolve()
    data = json.loads(path.read_text())
    if data.get("Version") not in (1, 2, 3):
        raise ValueError(f"{path} is not a supported process library")
    if data["Version"] >= 3:
        layers = data.get("Fabrication", {}).get("InterfaceLayers")
        if not isinstance(layers, dict):
            raise ValueError(
                f"{path} is version 3 but has no Fabrication.InterfaceLayers metadata"
            )
    models = data.get("Models")
    if not isinstance(models, list):
        raise ValueError(f"{path} has no Models array")
    radius = float(data["MatchingRadius"])
    if radius <= 0.0:
        raise ValueError(f"{path} has a nonpositive MatchingRadius")
    return path, data


def model_sha256(model, source_root):
    """Digest of one model as its input library records it: the entry (sorted keys) with every
    matrix / basis / trace-mesh path replaced by that file's sha256, so the digest follows the
    model's numbers, not where they are stored."""
    content = dict(model)

    def digest_of(relative):
        path = Path(relative)
        return {"SHA256": sha256(path if path.is_absolute() else Path(source_root) / path)}
    for field in PATH_FIELDS:
        if isinstance(content.get(field), str):
            content[field] = digest_of(content[field])
    if isinstance(content.get("TraceMesh"), dict):
        content["TraceMesh"] = {key: digest_of(value) for key, value in content["TraceMesh"].items()}
    return hashlib.sha256(json.dumps(content, sort_keys=True).encode()).hexdigest()


def parse_supersede(spec):
    """NAME=<record path> of --supersede; the record must exist and parse as JSON."""
    name, separator, record = spec.partition("=")
    if not separator or not name or not record:
        raise ValueError(f"--supersede {spec!r} is not NAME=<record path>")
    path = Path(record).expanduser()
    if not path.is_file():
        raise ValueError(f"--supersede {name}: the record {path} is missing")
    try:
        json.loads(path.read_text())
    except (OSError, ValueError) as error:
        raise ValueError(f"--supersede {name}: the record {path} is unreadable ({error})")
    return name, {"Path": str(path.resolve()), "SHA256": sha256(path)}


def parse_pending_allow(spec):
    """NAME=DECISION of --pending-allow (R-5, fix round P2; USER decision 1230 Q8 / supervisor 1260): the admission of
    a model whose explicit QualificationStatus is not "Qualified", by Name with the curator's decision text (the loader
    refuses such a model unless the library lists it under AllowPendingModels)."""
    name, separator, decision = spec.partition("=")
    if not separator or not name.strip() or not decision.strip():
        raise ValueError(f"--pending-allow {spec!r} is not NAME=DECISION (a nonempty name and a nonempty decision text)")
    return {"Name": name.strip(), "Decision": decision.strip()}


def allow_pending_models(admissions, models):
    """The AllowPendingModels leaf for the output: every admission must name a DISTINCT model of the output whose
    explicit QualificationStatus is not "Qualified" (a stale or misspelt admission is refused at assembly, never carried
    silently); the status-less models (the 2D planner's record form) need no admission."""
    pending = {model["Name"]: model["QualificationStatus"] for model in models
               if isinstance(model.get("QualificationStatus"), str) and model["QualificationStatus"] != "Qualified"}
    names = [entry["Name"] for entry in admissions]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    if duplicates:
        raise ValueError(f"--pending-allow names a model twice: {', '.join(duplicates)}")
    for entry in admissions:
        if entry["Name"] not in pending:
            raise ValueError(f"--pending-allow {entry['Name']}: not a model of the output with an explicit non-Qualified "
                             f"QualificationStatus (the non-Qualified models: {sorted(pending) or 'none'})")
    return admissions


def warn_unadmitted_pending_models(admissions, models):
    """The non-Qualified models of the output no --pending-allow admits: named on stdout (the loader refuses them at
    load, R-5); the output itself is unchanged (the flag is opt-in)."""
    admitted = {entry["Name"] for entry in admissions}
    unadmitted = [(model["Name"], model["QualificationStatus"]) for model in models
                  if isinstance(model.get("QualificationStatus"), str) and model["QualificationStatus"] != "Qualified"
                  and model["Name"] not in admitted]
    if unadmitted:
        print("Warning: the output carries non-Qualified models no --pending-allow admits; the loader refuses them "
              "(R-5, USER decision 1230 Q8): " + ", ".join(f"{name} ({status})" for name, status in unadmitted))
    return unadmitted


def supersede_decision(text, supersedes):
    """The decision text that rules every --supersede / --replace-by-signature replacement (decision 513 (2)):
    mandatory and non-blank whenever either is given, refused without one."""
    if supersedes and not (isinstance(text, str) and text.strip()):
        raise ValueError("--supersede / --replace-by-signature requires --supersede-decision TEXT (the decision that rules "
                         "the replacement, recorded in Supersedes.Decision)")
    if text is not None and not supersedes:
        raise ValueError("--supersede-decision is given without --supersede / --replace-by-signature")
    return text.strip() if text else None


def apply_supersedes(entries, base_path, supersedes, decision):
    """`entries` = [(source path, library, model)] in input order; for every --supersede NAME the
    BASE library's (the first input's) model NAME is replaced IN PLACE by the model of the same
    name of exactly one later input (which then contributes it nowhere else), so every other
    model keeps its index and directory.  Refuses a name absent from the base or not supplied
    by exactly one later input.  Returns (entries, Supersedes records)."""
    records = []
    for name, record in supersedes:
        base_slots = [i for i, (path, _, model) in enumerate(entries) if path == base_path and model.get("Name") == name]
        if not base_slots:
            raise ValueError(f"--supersede {name}: the base library {base_path} has no model {name!r}")
        new_slots = [i for i, (path, _, model) in enumerate(entries) if path != base_path and model.get("Name") == name]
        if len(new_slots) != 1:
            raise ValueError(f"--supersede {name}: {len(new_slots)} later inputs supply a model {name!r}; exactly one must")
        base_entry, new_entry = entries[base_slots[0]], entries[new_slots[0]]
        records.append({"Name": name, "BaseModelSHA": model_sha256(base_entry[2], base_entry[0].parent),
                        "NewModelSHA": model_sha256(new_entry[2], new_entry[0].parent),
                        "BaseSource": {"Path": str(base_entry[0]), "SHA256": sha256(base_entry[0])},
                        "NewSource": {"Path": str(new_entry[0]), "SHA256": sha256(new_entry[0])},
                        "Record": {"Path": record["Path"], "SHA256": record["SHA256"]}, "Decision": decision,
                        "Rule": SUPERSEDE_RULE})
        entries[base_slots[0]] = new_entry
        del entries[new_slots[0]]
    return entries, records


def apply_signature_replacements(entries, base_path, replacements, decision):
    """`entries` = [(source path, library, model)] in input order; for every --replace-by-signature NAME the model
    NAME of exactly one later input (which must carry a Signature) replaces IN PLACE the base library's (the first
    input's) model carrying the SAME canonical Signature (signature_library.signature_hash: exact equality, no
    tolerance), so every other model keeps its index and directory and the library never holds two models of one
    signature (decision 1103; the one-reference model gives way to its piece form, decisions 1136 (1) / 1140 (4)).
    Refuses a NAME absent from the later inputs or supplied twice, a new model without a Signature, a Signature the
    base holds in no model or in several, and a NAME the base already holds (that is --supersede's case). Returns
    (entries, Supersedes records with Kind "Signature")."""
    records = []
    for name, record in replacements:
        new_slots = [i for i, (path, _, model) in enumerate(entries) if path != base_path and model.get("Name") == name]
        if len(new_slots) != 1:
            raise ValueError(f"--replace-by-signature {name}: {len(new_slots)} later inputs supply a model {name!r}; exactly one must")
        if any(path == base_path and model.get("Name") == name for path, _, model in entries):
            raise ValueError(f"--replace-by-signature {name}: the base library {base_path} already holds a model {name!r} "
                             "(a replacement by name is --supersede)")
        new_entry = entries[new_slots[0]]
        signature = new_entry[2].get("Signature")
        if signature is None:
            raise ValueError(f"--replace-by-signature {name}: the model {name!r} carries no Signature")
        key = signature_hash(signature)
        base_slots = [i for i, (path, _, model) in enumerate(entries)
                      if path == base_path and model.get("Signature") is not None and signature_hash(model["Signature"]) == key]
        if len(base_slots) != 1:
            raise ValueError(f"--replace-by-signature {name}: the base library {base_path} holds {len(base_slots)} model(s) of "
                             f"the Signature {key[:12]} of {name!r}; exactly one must")
        base_entry = entries[base_slots[0]]
        records.append({"Kind": "Signature", "Name": base_entry[2].get("Name"), "NewName": name, "SignatureHash": key,
                        "BaseModelSHA": model_sha256(base_entry[2], base_entry[0].parent),
                        "NewModelSHA": model_sha256(new_entry[2], new_entry[0].parent),
                        "BaseSource": {"Path": str(base_entry[0]), "SHA256": sha256(base_entry[0])},
                        "NewSource": {"Path": str(new_entry[0]), "SHA256": sha256(new_entry[0])},
                        "Record": {"Path": record["Path"], "SHA256": record["SHA256"]}, "Decision": decision,
                        "Rule": REPLACE_BY_SIGNATURE_RULE})
        entries[base_slots[0]] = new_entry
        del entries[new_slots[0]]
    return entries, records


REPLACE_BY_SIGNATURE_RULE = ("--replace-by-signature NAME=<record> (idfix R0-b, decisions 1103 / 1136 (1) / 1140 (4)): the model "
                             "NAME of exactly one later input replaces IN PLACE the base library's model carrying the SAME canonical "
                             "Signature (exact equality of the json sort_keys form, signature_library.signature_hash; its matrices "
                             "copied under the base model's index), recorded under Supersedes with Kind Signature, both Names, both "
                             "model digests, the stored record and the mandatory --supersede-decision text; a library never holds two "
                             "models of one signature (the one-reference model gives way to its piece form)")


SUPERSEDE_RULE = ("--supersede NAME=<record>: the base library's model NAME is replaced in place by the model of the same "
                  "name of exactly one later input (its matrices copied under the base model's index / directory), "
                  "recorded under Supersedes with both model digests (the entry with its files' sha256, model_sha256), "
                  "the stored record and the mandatory --supersede-decision text that rule the replacement; "
                  "without the flag a repeated name is refused and the output is unchanged")


def merge_metadata(first, second, path="Fabrication"):
    if isinstance(first, dict) and isinstance(second, dict):
        result = {}
        for key in first.keys() | second.keys():
            if key not in first:
                result[key] = second[key]
            elif key not in second:
                result[key] = first[key]
            else:
                result[key] = merge_metadata(
                    first[key], second[key], f"{path}.{key}"
                )
        return result
    if (
        isinstance(first, (int, float))
        and not isinstance(first, bool)
        and isinstance(second, (int, float))
        and not isinstance(second, bool)
    ):
        if math.isclose(float(first), float(second), rel_tol=1.0e-12, abs_tol=0.0):
            return first
    elif first == second:
        return first
    raise ValueError(f"Input libraries contain different {path} metadata")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--name", default="combined-fabrication-process")
    parser.add_argument(
        "--corner-interpolation-qualification",
        action="append",
        type=Path,
        default=[],
        help="Passed held-out corner-radius interpolation report to embed",
    )
    parser.add_argument(
        "--allow-empty",
        action="store_true",
        help="Write a metadata-only version-3 library when no model qualified",
    )
    parser.add_argument(
        "--supersede",
        action="append",
        default=[],
        metavar="NAME=RECORD",
        help="replace the base (first) library's model NAME by the model of the same name of one later "
             "input, recording Supersedes {Name, BaseModelSHA, NewModelSHA, Record, Decision}; RECORD is the "
             "stored record (JSON) that rules the replacement; requires --supersede-decision",
    )
    parser.add_argument(
        "--replace-by-signature",
        action="append",
        default=[],
        metavar="NAME=RECORD",
        help="replace the base (first) library's model carrying the SAME canonical Signature as the model NAME of "
             "one later input by that model, in place (idfix R0-b, decisions 1103 / 1136 (1): the one-reference model "
             "gives way to its -pieces<P> form), recording Supersedes {Kind: Signature, Name, NewName, SignatureHash, "
             "BaseModelSHA, NewModelSHA, Record, Decision}; RECORD is the stored record (JSON) that rules the "
             "replacement; requires --supersede-decision",
    )
    parser.add_argument(
        "--supersede-decision",
        metavar="TEXT",
        help="the decision that rules every --supersede / --replace-by-signature replacement (mandatory with either; "
             "recorded as Supersedes.Decision)",
    )
    parser.add_argument(
        "--pending-allow",
        action="append",
        default=[],
        metavar="NAME=DECISION",
        help="admit the output model NAME whose explicit QualificationStatus is not \"Qualified\" (R-5, fix round P2; "
             "USER decision 1230 Q8 / supervisor 1260): recorded as the top-level AllowPendingModels [{Name, Decision}] the "
             "loader reads (it refuses such a model unless listed); every non-Qualified model of the output must be admitted "
             "or the assembly is refused; without the flag the output is unchanged",
    )
    parser.add_argument("libraries", type=Path, nargs="+")
    args = parser.parse_args()
    supersedes = [parse_supersede(spec) for spec in args.supersede]
    admissions = [parse_pending_allow(spec) for spec in args.pending_allow]
    replacements = [parse_supersede(spec) for spec in args.replace_by_signature]
    decision = supersede_decision(args.supersede_decision, supersedes + replacements)

    destination = args.output.expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    loaded = [load_library(path) for path in args.libraries]
    matching_radius = float(loaded[0][1]["MatchingRadius"])
    tolerance = 1.0e-12 * matching_radius
    for path, library in loaded[1:]:
        if abs(float(library["MatchingRadius"]) - matching_radius) > tolerance:
            raise ValueError(
                f"{path} uses MatchingRadius={library['MatchingRadius']}; "
                f"expected {matching_radius}"
            )

    version = max(library["Version"] for _, library in loaded)
    trace_lift_version = max(
        int(library.get("TraceLiftVersion", 0)) for _, library in loaded
    )
    fabrication = [library.get("Fabrication") for _, library in loaded]
    known_fabrication = [entry for entry in fabrication if entry is not None]
    merged_fabrication = None
    for entry in known_fabrication:
        merged_fabrication = (
            entry
            if merged_fabrication is None
            else merge_metadata(merged_fabrication, entry)
        )
    if known_fabrication and len(known_fabrication) != len(fabrication):
        print(
            "Warning: combining libraries with and without Fabrication metadata; "
            "the output is downgraded to version 2 until its process is recorded"
        )
        version = min(version, 2)

    result = {
        "Version": version,
        "Name": args.name,
        "MatchingRadius": matching_radius,
        "ExhaustiveSpatialClosure": True,
        "Models": [],
    }
    if trace_lift_version:
        result["TraceLiftVersion"] = trace_lift_version
    if known_fabrication and len(known_fabrication) == len(fabrication):
        result["Fabrication"] = merged_fabrication
    result["Sources"] = [source_record(path, library) for path, library in loaded]
    names = set()
    index = 0
    interpolation = []
    entries = [(source_path, library, source_model) for source_path, library in loaded for source_model in library["Models"]]
    if supersedes:
        entries, result["Supersedes"] = apply_supersedes(entries, loaded[0][0], supersedes, decision)
    if replacements:
        entries, replaced = apply_signature_replacements(entries, loaded[0][0], replacements, decision)
        result["Supersedes"] = result.get("Supersedes", []) + replaced
    refuse_duplicate_signatures([source_model for _, _, source_model in entries], library=str(destination / "process-library.json"))
    source_digests = {source_path: sha256(source_path) for source_path, _ in loaded}
    for source_path, library, source_model in entries:
        source_root = source_path.parent
        combined_from = {"Path": str(source_path), "SHA256": source_digests[source_path]}
        default_depth = library.get("CouponDepth")
        if (
            len(source_model.get("ConductorReferences", [])) > 1
            and int(library.get("TraceLiftVersion", 0)) < 2
        ):
            raise ValueError(
                f"{source_path} contains multiconductor model "
                f"{source_model.get('Name')!r} without TraceLiftVersion >= 2"
            )
        index += 1
        model = dict(source_model)
        name = model.get("Name", "")
        if not name or name in names:
            raise ValueError(
                f"Response model name {name!r} is empty or repeated"
            )
        names.add(name)
        if default_depth is not None and "CouponDepth" not in model:
            model["CouponDepth"] = default_depth

        relative_directory = model_directory(index, name)
        model_destination = destination / relative_directory
        model_destination.mkdir(parents=True, exist_ok=True)
        model["CombinedFrom"] = combined_from
        for field in PATH_FIELDS:
            if field not in model:
                continue
            if not isinstance(model[field], str):
                # Palace reads every matrix path as a string (ThinMatrix unconditionally):
                # a model without the file (a qualify NotLoadable entry) cannot be combined.
                raise ValueError(
                    f"{name} field {field} is {model[field]!r}, not a file path"
                    + (f" ({model['NotLoadable'].get('Reason')})" if isinstance(model.get("NotLoadable"), dict) else "")
                )
            source = Path(model[field])
            if not source.is_absolute():
                source = source_root / source
            source = source.resolve()
            if not source.is_file():
                raise FileNotFoundError(
                    f"{name} field {field} does not exist: {source}"
                )
            target = model_destination / PATH_NAMES[field]
            shutil.copy2(source, target)
            model[field] = str(relative_directory / target.name)
        if "TraceMesh" in model:
            trace_mesh = dict(model["TraceMesh"])
            for field, filename in (
                ("Vertices", "trace-vertices.csv"),
                ("Triangles", "trace-triangles.csv"),
            ):
                source = Path(trace_mesh[field])
                if not source.is_absolute():
                    source = source_root / source
                source = source.resolve()
                if not source.is_file():
                    raise FileNotFoundError(
                        f"{name} TraceMesh.{field} does not exist: {source}"
                    )
                target = model_destination / filename
                shutil.copy2(source, target)
                trace_mesh[field] = str(relative_directory / target.name)
            model["TraceMesh"] = trace_mesh
        result["Models"].append(model)
    for _, library in loaded:
        interpolation.extend(library.get("CornerRadiusInterpolation", []))

    for report_path in args.corner_interpolation_qualification:
        report_path = report_path.expanduser().resolve()
        report = json.loads(report_path.read_text())
        if (
            report.get("Study") != "CornerRadiusInterpolation"
            or not report.get("Passed", False)
            or not isinstance(report.get("LibraryRecord"), dict)
        ):
            raise ValueError(
                f"{report_path} is not a passed corner-radius interpolation report"
            )
        interpolation.append(report["LibraryRecord"])

    seen_spans = set()
    if interpolation and version < 3:
        raise ValueError(
            "CornerRadiusInterpolation requires only version-3 input libraries"
        )
    for span in interpolation:
        if not isinstance(span, dict):
            raise ValueError("CornerRadiusInterpolation entries must be objects")
        lower = span.get("LowerModel")
        upper = span.get("UpperModel")
        qualification = span.get("Qualification")
        if (
            lower not in names
            or upper not in names
            or lower == upper
            or not isinstance(qualification, dict)
            or qualification.get("Method") != "HeldOutCoupon"
            or not qualification.get("Passed", False)
        ):
            raise ValueError(
                "CornerRadiusInterpolation requires distinct included models and "
                "a passed held-out qualification"
            )
        if (lower, upper) in seen_spans:
            raise ValueError(
                f"Duplicate corner-radius interpolation span {lower} -> {upper}"
            )
        seen_spans.add((lower, upper))
    if interpolation:
        result["CornerRadiusInterpolation"] = interpolation

    if not result["Models"] and not (
        args.allow_empty and version >= 3 and merged_fabrication is not None
    ):
        raise ValueError("Combined library contains no response models")

    if admissions:
        result["AllowPendingModels"] = allow_pending_models(admissions, result["Models"])
    warn_unadmitted_pending_models(admissions, result["Models"])
    output = destination / "process-library.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"Wrote {len(result['Models'])} models at R={matching_radius:g} "
        f"to {output}"
    )


if __name__ == "__main__":
    main()
