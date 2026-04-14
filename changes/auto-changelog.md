### New Features

  - Added per-PR changelog fragment system in `changes/` directory to eliminate merge
    conflicts on `CHANGELOG.md`. Each PR includes a `changes/<branch-suffix>.md` file
    validated by CI, and fragments are collected into `CHANGELOG.md` at release time via
    `scripts/collect-changelog`.
