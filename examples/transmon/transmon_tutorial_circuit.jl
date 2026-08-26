# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# Transmon Circuit-Synthesis Tutorial — data generator

Generate all Palace output used by `transmon_tutorial_circuit_plots.jl`:

1. One uniform driven reference sweep.
2. Adaptive sweeps without circuit synthesis at AdaptiveTol ∈ {1e-1, 1e-3, 1e-5}.
3. Adaptive sweeps with circuit synthesis at the same tolerances.

All runs use `transmon_tutorial_driven.json` as a template. Results are written below
`examples/transmon/postpro/transmon_tutorial_circuit/`.

From the repository root:

```bash
julia --project=examples -e 'include("examples/transmon/transmon_tutorial_circuit.jl"); generate_transmon_circuit_data(num_processors=4)'
```

This requires `palace` to be a runnable command. Pass `palace_exec` explicitly when it
is not on `PATH`. The uniform reference is a multi-hour run.
=#

using JSON

const CIRCUIT_TUTORIAL_TOLS = [1e-1, 1e-3, 1e-5]

circuit_tol_label(tol) = "1e$(round(Int, log10(tol)))"

"""
    generate_transmon_circuit_data(;
        palace_exec="palace",
        num_processors=1,
    )

Run the uniform, adaptive, and circuit-synthesis simulations used by the transmon
circuit-synthesis feature guide.

# Arguments

  - `palace_exec`: Path or name of the Palace executable.
  - `num_processors`: Number of MPI ranks passed to Palace's `-np` option.
"""
function generate_transmon_circuit_data(; palace_exec="palace", num_processors::Integer=1)
    num_processors > 0 || throw(ArgumentError("num_processors must be positive"))

    palace_exec_is_path = occursin(Base.Filesystem.path_separator, palace_exec)
    if palace_exec_is_path
        palace_exec = isabspath(palace_exec) ? palace_exec : abspath(palace_exec)
    end

    transmon_dir = @__DIR__
    output_root = joinpath(transmon_dir, "postpro", "transmon_tutorial_circuit")

    # JSON.jl does not accept the // comments allowed by Palace's configuration parser.
    base_config = open(joinpath(transmon_dir, "transmon_tutorial_driven.json")) do f
        return JSON.parse(replace(read(f, String), r"//[^\n]*" => ""))
    end

    function run_palace(config, run_name)
        outdir = joinpath(output_root, run_name)
        config["Problem"]["Output"] =
            joinpath("postpro", "transmon_tutorial_circuit", run_name)

        tmp_name = "transmon_tutorial_circuit_$(run_name)_tmp.json"
        tmp_path = joinpath(transmon_dir, tmp_name)
        log_path = joinpath(transmon_dir, "transmon_tutorial_circuit_$(run_name).log")
        output_log = joinpath(outdir, "palace.log")

        open(tmp_path, "w") do f
            return JSON.print(f, config, 4)
        end

        try
            palace_cmd = Cmd(`$palace_exec -np $num_processors $tmp_name`; dir=transmon_dir)
            run(pipeline(palace_cmd, `tee $log_path`))
        finally
            isfile(tmp_path) && rm(tmp_path)
            isfile(log_path) && isdir(outdir) && mv(log_path, output_log; force=true)
        end
        return outdir
    end

    println("Running uniform reference sweep...")
    uniform = deepcopy(base_config)
    uniform["Solver"]["Driven"]["AdaptiveTol"] = 0.0
    uniform["Solver"]["Driven"]["AdaptiveCircuitSynthesis"] = false
    uniform["Solver"]["Driven"]["Restart"] = 1
    run_palace(uniform, "driven_uniform_reference")

    for tol in CIRCUIT_TUTORIAL_TOLS
        label = circuit_tol_label(tol)

        println(
            "Running adaptive sweep without circuit synthesis at AdaptiveTol = $label..."
        )
        adaptive = deepcopy(base_config)
        adaptive["Solver"]["Driven"]["AdaptiveTol"] = tol
        adaptive["Solver"]["Driven"]["AdaptiveCircuitSynthesis"] = false
        adaptive["Solver"]["Driven"]["Restart"] = 1
        run_palace(adaptive, "driven_adaptive_$label")

        println("Running adaptive sweep with circuit synthesis at AdaptiveTol = $label...")
        circuit = deepcopy(base_config)
        circuit["Solver"]["Driven"]["AdaptiveTol"] = tol
        circuit["Solver"]["Driven"]["AdaptiveCircuitSynthesis"] = true
        circuit["Solver"]["Driven"]["Restart"] = 1
        run_palace(circuit, "driven_circuit_$label")
    end

    println("Circuit-synthesis tutorial data written to $output_root")
    return output_root
end
