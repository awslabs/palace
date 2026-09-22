#!/usr/bin/env python3
"""Rewrite a qualify run's process-library.json as the Palace-loadable Version-3 library
(qualify_library.LIBRARY_RULE) without re-running any job.

The previous library is read only; every model is kept (the merge path of
qualify_library.process_library_entries with no case of its own): the header comes from
the models' source process-library.json files (consistent_header: the writer stops when
they disagree), the matrices / basis points / trace mesh are copied under
ROOT/models/<slug>/, a model without thin matrices gets ThinMatrix null and a NotLoadable
reason.  Writes ROOT/process-library.json and ROOT/process-library-preflight.json
(qualify_library.PREFLIGHT_RULE) and prints the loadable / not-loadable counts.
"""

import argparse
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import qualify_library  # noqa: E402


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--previous", type=Path, required=True, help="a qualify run's process-library.json (read only)")
    parser.add_argument("--root", type=Path, required=True, help="output root (ROOT/process-library.json, ROOT/models/)")
    args = parser.parse_args(argv)
    root = args.root.expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    library = qualify_library.process_library_entries([], {}, manifest_path=None, manifest=None, root=root,
                                                      merge_into=args.previous.expanduser().resolve())
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
