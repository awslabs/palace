# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Test
using CSV
using DataFrames

function testcase(
    testdir,
    testconfig,
    testpostpro;
    palace="palace",
    np=1,
    rtol=1.0e-6,
    atol=1.0e-18,
    excluded_columns=[],
    skip_rowcount=false
)
    if isempty(testdir)
        @info "$testdir/ is empty, skipping tests"
        return
    end
    palacedir     = dirname(dirname(@__DIR__))
    refdir        = joinpath(@__DIR__, "ref", testdir)
    refpostprodir = joinpath(refdir, testpostpro)
    exampledir    = joinpath(palacedir, "examples", testdir)
    postprodir    = joinpath(exampledir, "postpro", testpostpro)
    logdir        = joinpath(exampledir, "log")

    # Cleanup
    rm(postprodir; force=true, recursive=true)
    rm(logdir; force=true, recursive=true)
    mkdir(logdir)
    cd(exampledir)

    @testset "Simulation" begin
        # Run the example simulation
        logfile = "log.out"
        errfile = "err.out"
        proc = run(
            pipeline(
                ignorestatus(`$palace -np $np $testconfig`);
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
    end

    @testset "Results" begin
        # Test that directories were created
        @test isdir(postprodir)
        (~, dirs, files) = first(walkdir(postprodir))
        (~, ~, filesref) = first(walkdir(refpostprodir))
        metafiles = filter(x -> last(splitext(x)) != ".csv", files)
        @test length(dirs) == 1 && first(dirs) == "paraview" || (@show dirs; false)
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
            test = isapprox.(data, dataref; rtol=rtol, atol=atol)
            for (row, rowdataref, rowdata) in
                zip(eachrow(test), eachrow(dataref), eachrow(data))
                for (rowcol, rowcoldataref, rowcoldata) in
                    zip(pairs(row), pairs(rowdataref), pairs(rowdata))
                    if !last(rowcol)
                        @warn string(
                            "Regression test error at ",
                            "row $(rownumber(row)), column $(strip(string(first(rowcol)))): ",
                            "$(last(rowcoldataref)) ≉ $(last(rowcoldata))"
                        )
                    end
                    @test last(rowcol)
                end
            end
        end
    end
end
