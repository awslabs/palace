# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script runs Palace driven solver simulations on the CPW lumped-port example, for
the "Driven Solver: Uniform vs Adaptive" tutorial. Use the accompanying
`cpw_tutorial_lumped_driven_plots.py` python script to generate plots once the data is
generated.

The script runs:
1. A uniform frequency sweep (cpw_tutorial_lumped_uniform.json)
2. Adaptive frequency sweeps at AdaptiveTol [1e-1, 1e-2, 1e-3, 1e-4, 1e-5] (using
   cpw_tutorial_lumped_adaptive.json as a template)

Results are written to:
  - postpro/tutorial_driven_rom/lumped_uniform/ 
  - postpro/tutorial_driven_rom/driven_adaptive_1e-1/       (adaptive, tol = 0.1)
  - postpro/tutorial_driven_rom/driven_adaptive_1e-2/       (adaptive, tol = 1e-2)
  - postpro/tutorial_driven_rom/driven_adaptive_1e-3/       (adaptive, tol = 1e-3)
  - postpro/tutorial_driven_rom/driven_adaptive_1e-4/       (adaptive, tol = 1e-4)
  - postpro/tutorial_driven_rom/driven_adaptive_1e-5/       (adaptive, tol = 1e-5)

## Prerequisites

Install required Julia packages:
```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

## How to run

From the repository root:
```bash
julia --project=examples -e 'include("examples/cpw/cpw_tutorial_lumped_driven.jl"); generate_cpw_lumped_driven_data(num_processors=4)'
```

This requires `palace` to be a runnable command. If it is not, pass the path explicitly:
```bash
julia --project=examples -e 'include("examples/cpw/cpw_tutorial_lumped_driven.jl"); generate_cpw_lumped_driven_data(palace_exec="build/bin/palace", num_processors=4)'
```

=#

using JSON

"""
    generate_cpw_lumped_driven_data(;
                                     palace_exec::String="palace",
                                     num_processors::Integer=1
                                    )

Run Palace driven solver simulations for the "Driven Solver: Uniform and Adaptive" tutorial.

  - Uniform reference just runs cpw_tutorial_lumped_uniform.json.
  - Adaptive solvers uses cpw_tutorial_lumped_adaptive.json but modifies AdaptiveTol and
    output directory.

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

    # AdaptiveTol values to run; 0.0 triggers Palace's uniform driven solver and is
    # written to the `driven_uniform_reference` output directory used by the plot script.
    adaptive_tols = [0.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5]

    # Read the base adaptive config once (strip // comments, unsupported by JSON.jl)
    adaptive_config = open(joinpath(cpw_dir, "cpw_tutorial_lumped_adaptive.json")) do f
        content = replace(read(f, String), r"//[^\n]*" => "")
        return JSON.parse(content)
    end

    for tol in adaptive_tols
        if tol == 0.0
            tol_label = "uniform_reference"
            outdir_rel = "postpro/tutorial_driven_rom/driven_uniform_reference"
        else
            exp_val = round(Int, log10(tol))
            tol_label = "1e$(exp_val)"
            outdir_rel = "postpro/tutorial_driven_rom/driven_adaptive_$(tol_label)"
        end
        println("Running driven sweep with AdaptiveTol = $(tol)...")

        tmp_file = "cpw_tutorial_lumped_adaptive_tmp_$(tol_label).json"
        tmp_path = joinpath(cpw_dir, tmp_file)

        config = deepcopy(adaptive_config)
        config["Problem"]["Output"] = outdir_rel
        config["Solver"]["Driven"]["AdaptiveTol"] = tol

        open(tmp_path, "w") do f
            return JSON.print(f, config, 4)
        end

        outdir = joinpath(
            cpw_dir,
            "postpro",
            "tutorial_driven_rom",
            "driven_adaptive_$(tol_label)"
        )
        log_path = joinpath(cpw_dir, "cpw_tutorial_lumped_adaptive_$(tol_label).log")
        try
            palace_cmd = Cmd(`$palace_exec -np $num_processors $tmp_file`; dir=cpw_dir)
            run(pipeline(palace_cmd, `tee $log_path`))
        finally
            isfile(tmp_path) && rm(tmp_path)
            isfile(log_path) &&
                isdir(outdir) &&
                mv(log_path, joinpath(outdir, "palace.log"); force=true)
        end
    end

    return
end
