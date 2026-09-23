#!/usr/bin/env python3
"""Rewrite a qualify run's process-library.json as the Palace-loadable Version-3 library
(qualify_library.LIBRARY_RULE) without re-running any job.

The previous library is read only; every model is kept (the merge path of
qualify_library.process_library_entries with no case of its own): the header comes from
the models' source process-library.json files (consistent_header: the writer stops when
they disagree), the matrices / basis points / trace mesh are copied under
ROOT/models/<slug>/, a model without thin matrices gets ThinMatrix null and a NotLoadable
reason.  With --run, the analyzed cases of that finished qualify run (its
library-qualification.json; manifest through its build record) enter as the run's own
cases: a thin case pairs with the previous library's model of the same Name (the thin run
followed the fabricated run, decision 66; the driver's --merge-into without re-running or
re-fetching anything), a fabricated case replaces the model of its Name.  Writes
ROOT/process-library.json and ROOT/process-library-preflight.json
(qualify_library.PREFLIGHT_RULE) and prints the loadable / not-loadable counts.
"""

import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import qualify_library  # noqa: E402


def run_cases(run_path):
    """The records of a finished qualify run with the contexts process_library_entries reads
    (the recorded stage layout; the MA tail record re-read from its file) and the run's
    manifest (through its build record, as run_qualify binds it)."""
    run = json.loads(run_path.read_text())
    build = json.loads(Path(run["BuildRecord"]["Path"]).read_text())
    manifest_path = Path(build["Library"]["Manifest"]["Path"])
    manifest = json.loads(manifest_path.read_text())
    records = [record for record in run["Cases"] if record.get("Stages") and record.get("Qualification")]
    contexts = {}
    for record in records:
        tail_path = (record.get("MATail") or {}).get("Path")
        ma_tail = json.loads(Path(tail_path).read_text()) | {"Path": tail_path} if tail_path else None
        contexts[record["Case"]] = {"layout": record["Stages"], "ma_tail": ma_tail}
    return records, contexts, manifest_path, manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--previous", type=Path, required=True, help="a qualify run's process-library.json (read only)")
    parser.add_argument("--root", type=Path, required=True, help="output root (ROOT/process-library.json, ROOT/models/)")
    parser.add_argument("--run", type=Path, default=None,
                        help="a finished qualify run's library-qualification.json (read only): its analyzed cases enter as "
                             "this run's - a thin case supplies the thin matrices of the previous library's model of its Name")
    args = parser.parse_args(argv)
    root = args.root.expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    records, contexts, manifest_path, manifest = [], {}, None, None
    if args.run is not None:
        records, contexts, manifest_path, manifest = run_cases(args.run.expanduser().resolve())
    library = qualify_library.process_library_entries(records, contexts, manifest_path=manifest_path, manifest=manifest,
                                                      root=root, merge_into=args.previous.expanduser().resolve())
    qualify_library.write_json(root / qualify_library.PROCESS_LIBRARY_RECORD, library)
    qualify_library.write_json(root / qualify_library.PROCESS_LIBRARY_PREFLIGHT_RECORD,
                               qualify_library.preflight_process_library(library))
    loadable = library["Loadable"]
    print(f"{root / qualify_library.PROCESS_LIBRARY_RECORD}: Version {library['Version']} MatchingRadius {library['MatchingRadius']} "
          f"models {len(library['Models'])} loadable {loadable['Models']} not-loadable {len(loadable['NotLoadable'])}; "
          f"preflight variant {root / qualify_library.PROCESS_LIBRARY_PREFLIGHT_RECORD}")
    for model in library["Models"]:
        if model["NotLoadable"] is not None:
            print(f"  {model['Name']}: NotLoadable - {model['NotLoadable']['Reason']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
