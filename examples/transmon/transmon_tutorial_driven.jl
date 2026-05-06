# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script runs Palace driven solver simulations on the transmon example, for
the "Driven Solver: Uniform vs Adaptive" tutorial. Use the accompanying
`transmon_tutorial_driven_plots.py` python script to generate plots once the data is
generated.

The script runs:
1. A uniform frequency sweep (AdaptiveTol = 0.0)
2. Adaptive frequency sweeps at AdaptiveTol ∈ {1e-1, 1e-2, 1e-3, 1e-4, 1e-5}

Both use `transmon_tutorial_driven.json` as a template, modifying only the `AdaptiveTol`
and `Output` fields for each run.

Results are written to:
  - postpro/transmon_tutorial_driven_rom/driven_uniform_reference/
  - postpro/transmon_tutorial_driven_rom/driven_adaptive_1e-1/
  - postpro/transmon_tutorial_driven_rom/driven_adaptive_1e-2/
  - postpro/transmon_tutorial_driven_rom/driven_adaptive_1e-3/
  - postpro/transmon_tutorial_driven_rom/driven_adaptive_1e-4/
  - postpro/transmon_tutorial_driven_rom/driven_adaptive_1e-5/

## Prerequisites

Install required Julia packages:
```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

## How to run

From the repository root:
```bash
julia --project=examples -e 'include("examples/transmon/transmon_tutorial_driven.jl"); generate_transmon_driven_data(num_processors=4)'
```

This requires `palace` to be a runnable command. If it is not, pass the path explicitly:
```bash
julia --project=examples -e 'include("examples/transmon/transmon_tutorial_driven.jl"); generate_transmon_driven_data(palace_exec="build/bin/palace", num_processors=4)'
```

=#

using JSON

"""
    generate_transmon_driven_data(;
                                   palace_exec::String="palace",
                                   num_processors::Integer=1
                                  )

Run Palace driven solver simulations for the "Driven Solver: Uniform and Adaptive" tutorial
on the transmon model.

Executes a uniform frequency sweep followed by adaptive sweeps at AdaptiveTol ∈
{1e-1, 1e-2, 1e-3, 1e-4, 1e-5}. Both use `transmon_tutorial_driven.json` as a template.

# Arguments

  - `palace_exec`     - path or name of the Palace executable
  - `num_processors`  - number of MPI ranks to use
"""
function generate_transmon_driven_data(; palace_exec="palace", num_processors::Integer=1)
    palace_exec_is_path = occursin(Base.Filesystem.path_separator, palace_exec)
    if palace_exec_is_path
        palace_exec = isabspath(palace_exec) ? palace_exec : abspath(palace_exec)
    end

    transmon_dir = @__DIR__

    # Read the base config once (strip // comments, unsupported by JSON.jl)
    base_config = open(joinpath(transmon_dir, "transmon_tutorial_driven.json")) do f
        content = replace(read(f, String), r"//[^\n]*" => "")
        return JSON.parse(content)
    end

    # Helper: write a modified config, run Palace with tee, move log to outdir
    function run_palace(config, outdir, tmp_name, log_name)
        tmp_path = joinpath(transmon_dir, tmp_name)
        log_path = joinpath(transmon_dir, log_name)
        open(tmp_path, "w") do f
            return JSON.print(f, config, 4)
        end
        try
            palace_cmd = Cmd(`$palace_exec -np $num_processors $tmp_name`; dir=transmon_dir)
            run(pipeline(palace_cmd, `tee $log_path`))
        finally
            isfile(tmp_path) && rm(tmp_path)
            isfile(log_path) &&
                isdir(outdir) &&
                mv(log_path, joinpath(outdir, "palace.log"); force=true)
        end
    end

    # Uniform sweep (AdaptiveTol = 0.0 disables adaptive sampling)
    println("Running uniform sweep...")
    uniform_config = deepcopy(base_config)
    uniform_config["Problem"]["Output"] = "postpro/transmon_tutorial_driven_rom/driven_uniform_reference"
    uniform_config["Solver"]["Driven"]["AdaptiveTol"] = 0.0
    uniform_outdir = joinpath(
        transmon_dir,
        "postpro",
        "transmon_tutorial_driven_rom",
        "driven_uniform_reference"
    )
    run_palace(
        uniform_config,
        uniform_outdir,
        "transmon_tutorial_driven_uniform_tmp.json",
        "transmon_tutorial_driven_uniform.log"
    )

    # Adaptive sweeps
    adaptive_tols = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]

    for tol in adaptive_tols
        exp_val = round(Int, log10(tol))
        tol_label = "1e$(exp_val)"
        println("Running adaptive sweep with AdaptiveTol = $(tol_label)...")

        config = deepcopy(base_config)
        config["Problem"]["Output"] = "postpro/transmon_tutorial_driven_rom/driven_adaptive_$(tol_label)"
        config["Solver"]["Driven"]["AdaptiveTol"] = tol

        outdir = joinpath(
            transmon_dir,
            "postpro",
            "transmon_tutorial_driven_rom",
            "driven_adaptive_$(tol_label)"
        )
        run_palace(
            config,
            outdir,
            "transmon_tutorial_driven_adaptive_tmp_$(tol_label).json",
            "transmon_tutorial_driven_adaptive_$(tol_label).log"
        )
    end

    return
end
