# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script runs Palace driven solver simulations on the CPW lumped-port example,
comparing a uniform frequency sweep against adaptive sweeps at multiple tolerances.

The script runs:
1. A uniform frequency sweep (cpw_lumped_uniform.json, as-is)
2. Adaptive frequency sweeps at AdaptiveTol ∈ {1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5}

Results are written to:
  - postpro/lumped_uniform/             (uniform sweep, output path from the JSON)
  - postpro/driven_adaptive_1e1/        (adaptive, tol = 10)
  - postpro/driven_adaptive_1e0/        (adaptive, tol = 1)
  - postpro/driven_adaptive_1e-1/       (adaptive, tol = 0.1)
  - postpro/driven_adaptive_1e-2/       (adaptive, tol = 1e-2)
  - postpro/driven_adaptive_1e-3/       (adaptive, tol = 1e-3)
  - postpro/driven_adaptive_1e-4/       (adaptive, tol = 1e-4)
  - postpro/driven_adaptive_1e-5/       (adaptive, tol = 1e-5)

## Prerequisites

Install required Julia packages:
```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

## How to run

From the repository root:
```bash
julia --project=examples -e 'include("examples/cpw/cpw_lumped_driven.jl"); generate_cpw_lumped_driven_data(num_processors=4)'
```

This requires `palace` to be a runnable command. If it is not, pass the path explicitly:
```bash
julia --project=examples -e 'include("examples/cpw/cpw_lumped_driven.jl"); generate_cpw_lumped_driven_data(palace_exec="build/bin/palace", num_processors=4)'
```

## Output

Use the accompanying `cpw_lumped_driven.ipynb` notebook to visualize the results as
convergence curves comparing the adaptive sweeps against the uniform reference.
=#

using JSON

"""
    generate_cpw_lumped_driven_data(;
                                     palace_exec::String="palace",
                                     num_processors::Integer=1
                                    )

Run Palace driven solver simulations for the CPW lumped-port example.

Executes a uniform frequency sweep (`cpw_lumped_uniform.json`) followed by
adaptive sweeps at `AdaptiveTol` ∈ {1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5},
derived from `cpw_lumped_adaptive.json` with the tolerance and output directory
modified for each run.

# Arguments

  - `palace_exec`     - path or name of the Palace executable
  - `num_processors`  - number of MPI ranks to use
"""
function generate_cpw_lumped_driven_data(; palace_exec="palace", num_processors::Integer=1)
    palace_exec_is_path = occursin(Base.Filesystem.path_separator, palace_exec)
    if palace_exec_is_path
        palace_exec = isabspath(palace_exec) ? palace_exec : abspath(palace_exec)
    end

    cpw_dir = @__DIR__

    # 1. Run uniform sweep as-is
    println("Running uniform sweep...")
    run(
        Cmd(
            `$palace_exec -np $num_processors cpw_lumped_uniform_convergence.json`;
            dir=cpw_dir
        )
    )

    # 2. Run adaptive sweeps over a range of tolerances
    adaptive_tols = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]

    # Read the base adaptive config once (strip // comments, unsupported by JSON.jl)
    adaptive_config = open(joinpath(cpw_dir, "cpw_lumped_adaptive_convergence.json")) do f
        content = replace(read(f, String), r"//[^\n]*" => "")
        return JSON.parse(content)
    end

    tmpfile = "cpw_lumped_adaptive_temp.json"
    tmppath = joinpath(cpw_dir, tmpfile)

    for tol in adaptive_tols
        exp_val = round(Int, log10(tol))
        tol_label = "1e$(exp_val)"
        println("Running adaptive sweep with AdaptiveTol = $(tol_label)...")

        config = deepcopy(adaptive_config)
        config["Problem"]["Output"] = "postpro/driven_adaptive_$(tol_label)"
        config["Solver"]["Driven"]["AdaptiveTol"] = tol

        open(tmppath, "w") do f
            return JSON.print(f, config, 4)
        end

        outdir = joinpath(cpw_dir, "postpro", "driven_adaptive_$(tol_label)")
        logpath = joinpath(cpw_dir, "cpw_lumped_adaptive_$(tol_label).log")
        try
            palace_cmd = Cmd(`$palace_exec -np $num_processors $tmpfile`; dir=cpw_dir)
            run(pipeline(palace_cmd, `tee $logpath`))
            isfile(logpath) && mv(logpath, joinpath(outdir, "palace.log"); force=true)
        finally
            isfile(tmppath) && rm(tmppath)
        end
    end

    return
end
