# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# Electrostatic MMS documentation data generator

Run the existing serial electrostatic MMS tests with structured reporting enabled, validate
the records, and write the CSV consumed by `generate_electrostatic_mms_plots.jl`.

From the repository root:

```bash
julia --project=examples docs/generate_electrostatic_mms_data.jl \
  --test-exe build/palace-build/test/unit/palace-unit-tests
```
=#

using CSV
using DataFrames

const DOCS_DIR = @__DIR__
const DATASETS = ["cartesian-smooth", "curved-polynomial", "affine-polynomial"]
const ELEMENTS = ["hexahedron", "curved-tetrahedron", "tetrahedron"]
const COLUMNS = [
    "dataset",
    "solution",
    "element",
    "order",
    "resolution_kind",
    "resolution",
    "h",
    "primary_dofs",
    "potential_l2_error",
    "electric_l2_error"
]

default_test_exe() =
    joinpath(DOCS_DIR, "..", "build", "palace-build", "test", "unit", "palace-unit-tests")
default_output_path() =
    joinpath(DOCS_DIR, "src", "assets", "verification", "electrostatic_mms.csv")

function collect_data(test_exe)
    command = addenv(
        `$test_exe "[electrostatic][mms][analytic]" --skip-benchmarks`,
        "PALACE_MMS_REPORT" => "1"
    )
    output = read(pipeline(command, stderr=stderr), String)
    records = [
        line[(length("MMS_DATA,") + 1):end] for
        line in eachline(IOBuffer(output)) if startswith(line, "MMS_DATA,")
    ]
    length(records) == 14 ||
        error("Expected 14 electrostatic MMS data rows, found $(length(records))")

    csv = join([join(COLUMNS, ","); records], "\n")
    data = CSV.read(IOBuffer(csv), DataFrame)
    names(data) == COLUMNS || error("Unexpected electrostatic MMS data columns")
    count(==("cartesian-smooth"), data.dataset) == 9 ||
        error("Expected 9 Cartesian convergence rows")
    count(==("curved-polynomial"), data.dataset) == 3 ||
        error("Expected 3 curved convergence rows")
    count(==("affine-polynomial"), data.dataset) == 2 ||
        error("Expected 2 affine exactness rows")
    all(data.primary_dofs .> 0) || error("All H1 DoF counts must be positive")
    all(data.potential_l2_error .> 0) || error("All potential errors must be positive")
    all(data.electric_l2_error .> 0) || error("All electric-field errors must be positive")

    data.dataset_order = [findfirst(==(value), DATASETS) for value in data.dataset]
    data.element_order = [findfirst(==(value), ELEMENTS) for value in data.element]
    sort!(data, [:dataset_order, :order, :resolution, :element_order])
    select!(data, Not([:dataset_order, :element_order]))
    return data
end

function main(; test_exe=default_test_exe(), output_path=default_output_path())
    isfile(test_exe) || error("Unit-test executable not found: $test_exe")
    data = collect_data(test_exe)
    mkpath(dirname(output_path))
    CSV.write(output_path, data)
    println("Wrote electrostatic MMS verification data to $output_path")
    return output_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    let test_exe = default_test_exe(), output_path = default_output_path()
        args = copy(ARGS)
        while !isempty(args)
            arg = popfirst!(args)
            if arg == "--test-exe"
                test_exe = popfirst!(args)
            elseif arg == "--out"
                output_path = popfirst!(args)
            else
                error("Unknown argument: $arg")
            end
        end
        main(; test_exe, output_path)
    end
end
