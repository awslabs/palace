# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using CSV
using DataFrames
using Test

function testcase(testdir, testconfig, testpostpro; np=1, rtol=1.0e-6, atol=1.0e-18)
    if isempty(testdir)
        @info "$testdir/ is empty, skipping tests"
        return
    end
    palacedir     = dirname(@__DIR__)
    refdir        = joinpath(@__DIR__, "ref", testdir)
    refpostprodir = joinpath(refdir, testpostpro)
    exampledir    = joinpath(palacedir, "examples", testdir)
    postprodir    = joinpath(exampledir, "postpro", testpostpro)
    logdir        = joinpath(exampledir, "log")

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
                ignorestatus(`palace -np $np -wdir $exampledir $testconfig`);
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
        (~, ~, filesref) = first(walkdir(refpostprodir))
        (~, dirs, files) = first(walkdir(postprodir))
        metafiles = filter(x -> last(splitext(x)) != ".csv", files)
        @test length(dirs) == 1 && first(dirs) == "paraview"
        @test length(metafiles) == 1 && first(metafiles) == "palace.json"
        @test length(filter(x -> last(splitext(x)) == ".csv", files)) == length(filesref)

        # Test the simulation outputs
        for file in filesref
            dataref = CSV.File(joinpath(refpostprodir, file); header=1) |> DataFrame
            data    = CSV.File(joinpath(postprodir, file); header=1) |> DataFrame

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
