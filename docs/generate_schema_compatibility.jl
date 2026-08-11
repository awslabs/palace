# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Updates the schema compatibility table in docs/src/developer/notes.md from the
# machine-readable data in scripts/schema/schema-compatibility.json.

using JSON

const COMPATIBILITY_PATH =
    joinpath(@__DIR__, "..", "scripts", "schema", "schema-compatibility.json")
const NOTES_PATH = joinpath(@__DIR__, "src", "developer", "notes.md")
const TABLE_BEGIN = "<!-- BEGIN GENERATED SCHEMA COMPATIBILITY TABLE -->"
const TABLE_END = "<!-- END GENERATED SCHEMA COMPATIBILITY TABLE -->"

function markdown_row(values::Vector{String}, widths::Vector{Int})::String
    return "| " *
           join([rpad(value, width) for (value, width) in zip(values, widths)], " | ") *
           " |"
end

function compatibility_table(entries::Vector)::String
    isempty(entries) && error("Schema compatibility data must contain at least one entry")

    rows = Vector{String}[]
    for (index, entry) in enumerate(entries)
        fields = ["schema_version", "first_palace_release", "notes"]
        for field in fields
            haskey(entry, field) ||
                error("Compatibility entry $index is missing \"$field\"")
            entry[field] isa String ||
                error("Compatibility entry $index field \"$field\" must be a string")
        end
        push!(
            rows,
            [
                "`$(entry["schema_version"])`",
                "`$(entry["first_palace_release"])`",
                entry["notes"]
            ]
        )
    end

    headers = ["Schema version", "First *Palace* release", "Notes"]
    widths = [maximum(length(row[i]) for row in [headers, rows...]) for i = 1:3]
    separators = [":" * repeat("-", width - 2) * ":" for width in widths]
    separators[3] = ":" * repeat("-", widths[3] - 1)

    lines = [markdown_row(headers, widths), markdown_row(separators, widths)]
    append!(lines, markdown_row(row, widths) for row in rows)
    return join(lines, "\n")
end

function generate_schema_compatibility_table(;
    compatibility_path::String=COMPATIBILITY_PATH,
    notes_path::String=NOTES_PATH
)
    data = JSON.parsefile(compatibility_path)
    haskey(data, "schema_versions") ||
        error("Compatibility data is missing the \"schema_versions\" array")
    data["schema_versions"] isa Vector ||
        error("Compatibility data field \"schema_versions\" must be an array")

    source = read(notes_path, String)
    length(findall(TABLE_BEGIN, source)) == 1 ||
        error("Expected exactly one compatibility table start marker in $notes_path")
    length(findall(TABLE_END, source)) == 1 ||
        error("Expected exactly one compatibility table end marker in $notes_path")

    generated = "$TABLE_BEGIN\n$(compatibility_table(data["schema_versions"]))\n$TABLE_END"
    pattern = Regex("(?s)" * TABLE_BEGIN * ".*?" * TABLE_END)
    write(notes_path, replace(source, pattern => generated))
    @info "Generated schema compatibility table in $notes_path"
end

if abspath(PROGRAM_FILE) == @__FILE__
    generate_schema_compatibility_table()
end
