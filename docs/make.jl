# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Run with `julia --project make.jl` from within the `docs` folder.
# Output HTML is saved to the `build` folder.

using Documenter

# For tagged releases, rewrite GitHub links in the documentation source files to point to
# the tagged version instead of the main branch. This ensures that users consulting the
# documentation for a specific release are directed to the correct version of example files,
# configuration files, and other source code.
const PALACE_REPO_URL = "https://github.com/awslabs/palace"

function get_git_ref()
    # In GitHub Actions, GITHUB_REF_TYPE is "tag" for tag pushes and GITHUB_REF_NAME is the
    # tag name (e.g., "v0.13.0"). For branch pushes, GITHUB_REF_TYPE is "branch".
    ref_type = get(ENV, "GITHUB_REF_TYPE", "branch")
    ref_name = get(ENV, "GITHUB_REF_NAME", "main")
    if ref_type == "tag"
        return ref_name
    end
    return "main"
end

"""
    rewrite_github_links(src_dir, git_ref)

Replace GitHub links pointing to the `main` branch with links to `git_ref` in all Markdown
files under `src_dir`. Handles both file links (`blob/main`) and directory links
(`tree/main`). Does nothing when `git_ref` is `"main"`.
"""
function rewrite_github_links(src_dir, git_ref)
    git_ref == "main" && return  # Nothing to do for dev builds
    patterns = [
        "$PALACE_REPO_URL/blob/main" => "$PALACE_REPO_URL/blob/$git_ref",
        "$PALACE_REPO_URL/tree/main" => "$PALACE_REPO_URL/tree/$git_ref"
    ]
    for (root, dirs, files) in walkdir(src_dir)
        for file in files
            endswith(file, ".md") || continue
            filepath = joinpath(root, file)
            content = read(filepath, String)
            new_content = content
            for (old, new) in patterns
                new_content = replace(new_content, old => new)
            end
            if new_content != content
                write(filepath, new_content)
                @info "Rewrote GitHub links in $filepath to point to $git_ref"
            end
        end
    end
end

git_ref = get_git_ref()
@info "Building documentation with GitHub links pointing to: $git_ref"
rewrite_github_links(joinpath(@__DIR__, "src"), git_ref)

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
    authors="The Palace Developers and Maintainers, palace-maint@amazon.com",
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
            "examples/transmon.md",
            "examples/cylinder.md",
            "examples/coaxial.md",
            "examples/cpw.md",
            "examples/cpw2d.md"
        ],
        "faq.md",
        "For Developers" => Any[
            "developer/notes.md",
            "developer/testing.md",
            "developer/spack.md",
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
