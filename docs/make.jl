# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Run with `julia --project make.jl` from within the `docs` folder.
# Output HTML is saved to the `build` folder.

using Documenter

makedocs(
    format=Documenter.HTML(
        # Always use clean URLs so that raw HTML works, view local builds using a local
        # HTTP server with `python3 -m http.server`, for example
        prettyurls=true,
        sidebar_sitename=false,
        collapselevel=2,
        assets=["assets/favicon.ico"]
    ),
    sitename="Palace",
    authors="Sebastian Grimberg, sjg@amazon.com",
    pages=[
        "Home" => "index.md",
        "Quick Start" => "quick.md",
        "install.md",
        "run.md",
        "User Guide" => Any[
            "guide/guide.md",
            "guide/problem.md",
            "guide/model.md",
            "guide/boundaries.md",
            "guide/postprocessing.md",
            "guide/parallelism.md"
        ],
        "Configuration File" => Any[
            "config/config.md",
            "config/problem.md",
            "config/model.md",
            "config/domains.md",
            "config/boundaries.md",
            "config/solver.md"
        ],
        "Features" => Any["features/farfield.md",],
        "Examples" => Any[
            "examples/examples.md",
            "examples/spheres.md",
            "examples/rings.md",
            "examples/antenna.md",
            "examples/cylinder.md",
            "examples/coaxial.md",
            "examples/cpw.md"
        ],
        "For Developers" => Any[
            "developer/notes.md",
            "developer/testing.md",
            "developer/tutorial_add_new_unit_test.md",
            "developer/tutorial_gpu_profiling.md"
        ],
        "reference.md"
    ]
)

deploydocs(
    repo="github.com/awslabs/palace.git",
    devbranch="main",
    push_preview=true,
    forcepush=true
)
