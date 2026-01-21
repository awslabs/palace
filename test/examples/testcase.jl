# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Test
using CSV
using DataFrames
using JSON

# custom_tests is a dictionary that maps filenames to functions of the signature
# (new_data, ref_data) -> None, where arbitrary tests can be implemented
function testcase(
    testdir,
    testconfig,
    testpostpro;
    palace=["palace"],
    np=1,
    rtol=1.0e-6,
    atol=1.0e-18,
    excluded_columns=[],
    custom_tests=Dict(),
    paraview_fields=true,
    gridfunction_fields=false,
    skip_rowcount=false,
    generate_data=true,
    device="CPU",
    linear_solver="Default",
    eigen_solver="Default"
)
    if isempty(testdir)
        @info "$testdir/ is empty, skipping tests"
        return
    end

    # Check that palace executable exists
    if isnothing(Base.Sys.which(first(palace)))
        error("Executable `$(first(palace))` not found. ")
    end

    palacedir     = dirname(dirname(@__DIR__))
    refdir        = joinpath(@__DIR__, "ref", testdir)
    refpostprodir = joinpath(refdir, testpostpro)
    exampledir    = joinpath(palacedir, "examples", testdir)
    postprodir    = joinpath(exampledir, "postpro", testpostpro)
    logdir        = joinpath(exampledir, "log")

    cd(exampledir)

    # Adjust configuration depending on arguments
    config_to_use = testconfig
    temp_config = nothing

    config_content = read(testconfig, String)
    # Strip comments starting with //
    lines = split(config_content, '\n')
    filtered_lines = [split(line, "//")[1] for line in lines]
    clean_content = join(filtered_lines, '\n')

    # Remove trailing commas (handles multi-line cases)
    clean_content = replace(clean_content, r",(\s*[\r\n]*\s*[}\]])" => s"\1")

    config_json = JSON.parse(clean_content)
    if !haskey(config_json, "Solver")
        config_json["Solver"] = Dict()
        config_json["Solver"]["Linear"] = Dict()
    end
    haskey(config_json["Solver"], "Linear") || (config_json["Solver"]["Linear"] = Dict())

    config_json["Solver"]["Device"] = device

    config_json["Solver"]["Linear"]["Type"] = linear_solver
    haskey(config_json["Solver"], "Eigenmode") &&
        (config_json["Solver"]["Eigenmode"]["Type"] = eigen_solver)

    temp_config = tempname() * ".json"
    write(temp_config, JSON.json(config_json, 2))
    config_to_use = temp_config

    if generate_data
        # Cleanup
        rm(postprodir; force=true, recursive=true)
        rm(logdir; force=true, recursive=true)
        mkdir(logdir)

        @testset "Simulation" begin
            # Run the example simulation
            logfile = "log.out"
            errfile = "err.out"

            proc = run(
                pipeline(
                    ignorestatus(`$palace -np $np $config_to_use`);
                    stdout=joinpath(logdir, logfile),
                    stderr=joinpath(logdir, errfile)
                )
            )
            if proc.exitcode != 0
                @warn "Simulation exited with a nonzero exit code"
                if isfile(joinpath(exampledir, logdir, logfile))
                    @warn "Contents of stdout:"
                    println(String(read(joinpath(exampledir, logdir, logfile))))
                end
                if isfile(joinpath(exampledir, logdir, errfile))
                    @warn "Contents of stderr:"
                    println(String(read(joinpath(exampledir, logdir, errfile))))
                end
            end
            @test proc.exitcode == 0

            # Check for GPU availability when device is GPU.
            # Palace does not crash and revert to CPU if a GPU is not available,
            # but this could lead to false positives
            if device == "GPU" && isfile(joinpath(exampledir, logdir, logfile))
                log_content = String(read(joinpath(exampledir, logdir, logfile)))
                if occursin(
                    "Palace must be built with either CUDA or HIP support for GPU device usage, reverting to CPU!",
                    log_content
                )
                    @error "GPU was requested but Palace was not built with GPU support"
                    @test false
                end
            end
        end
    end

    # Clean up temporary config file
    if !isnothing(temp_config)
        rm(temp_config; force=true)
    end

    @testset "Results" begin
        # Test that directories were created
        @test isdir(postprodir)
        (~, dirs, files) = first(walkdir(postprodir))
        (~, ~, filesref) = first(walkdir(refpostprodir))
        metafiles = filter(x -> last(splitext(x)) != ".csv", files)
        if (paraview_fields && gridfunction_fields)
            @test length(dirs) == 2 && "paraview" in dirs && "gridfunction" in dirs ||
                  (@show dirs; false)
        elseif (paraview_fields)
            @test length(dirs) == 1 && first(dirs) == "paraview" || (@show dirs; false)
        elseif (gridfunction_fields)
            @test length(dirs) == 1 && first(dirs) == "gridfunction" || (@show dirs; false)
        end
        # When using AMR, `iterationN` directories are created
        @test length(dirs) >= 1 && last(dirs) == "paraview" || (@show dirs; false)
        @test length(metafiles) == 1 && first(metafiles) == "palace.json" ||
              (@show metafiles; false)
        @test length(filter(x -> last(splitext(x)) == ".csv", files)) == length(filesref) ||
              (@show filesref; false)

        # Helper to extract the stdout and stderr files and dump their contents.
        # Useful when debugging a failure
        function logdump(data...)
            @show data
            logfile = "log.out"
            errfile = "err.out"
            if isfile(joinpath(exampledir, logdir, logfile))
                @warn "Contents of stdout:"
                println(String(read(joinpath(exampledir, logdir, logfile))))
            end
            if isfile(joinpath(exampledir, logdir, errfile))
                @warn "Contents of stderr:"
                println(String(read(joinpath(exampledir, logdir, errfile))))
            end
            return false
        end

        # Test the simulation outputs
        for file in filesref
            data    = CSV.File(joinpath(postprodir, file); header=1) |> DataFrame
            dataref = CSV.File(joinpath(refpostprodir, file); header=1) |> DataFrame
            if !skip_rowcount
                @test nrow(data) == nrow(dataref) || logdump(data, dataref)
            end
            data = data[1:min(nrow(data), nrow(dataref)), :]

            # Check the number of columns matches, before removing any excluded columns
            @test ncol(data) == ncol(dataref)
            for col ∈ excluded_columns
                select!(data, Not(Cols(contains(col))))
                select!(dataref, Not(Cols(contains(col))))
            end
            rename!(data, strip.(names(data)))
            rename!(dataref, strip.(names(dataref)))

            @test names(data) == names(dataref) || logdump(names(data), names(dataref))

            if file in keys(custom_tests)
                custom_tests[file](data, dataref)
            else
                test = isapprox.(data, dataref; rtol=rtol, atol=atol)
                for (row, rowdataref, rowdata) in
                    zip(eachrow(test), eachrow(dataref), eachrow(data))
                    for (rowcol, rowcoldataref, rowcoldata) in
                        zip(pairs(row), pairs(rowdataref), pairs(rowdata))
                        if !last(rowcol)
                            ref_val = last(rowcoldataref)
                            new_val = last(rowcoldata)
                            abs_err = abs(new_val - ref_val)
                            rel_err = abs_err / max(abs(ref_val), abs(new_val), eps())
                            @warn string(
                                "Regression test error in file '$(file)' at ",
                                "row $(rownumber(row)), column $(strip(string(first(rowcol)))): ",
                                "$(ref_val) ≉ $(new_val), ",
                                "abs_err = $(abs_err) (atol=$(atol)), ",
                                "rel_err = $(rel_err) (rtol=$(rtol))"
                            )
                        end
                        @test last(rowcol)
                    end
                end
            end
        end
    end
end
