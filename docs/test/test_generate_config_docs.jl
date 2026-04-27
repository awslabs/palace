# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Tests for docs/generate_config_docs.jl
using Test

include(joinpath(@__DIR__, "..", "generate_config_docs.jl"))

@testset "load_schema" begin
    schema = load_schema("problem")
    @test haskey(schema, "properties")
    @test haskey(schema["properties"], "Type")
end

@testset "resolve_ref" begin
    schema = load_schema("boundaries")
    resolved = resolve_ref(schema, "#/\$defs/Attributes")
    @test resolved["type"] == "array"
end

@testset "inline_ref merges keys" begin
    root = load_schema("boundaries")
    raw = root["properties"]["PEC"]["properties"]["Attributes"]
    merged = inline_ref(root, raw)
    @test merged["type"] == "array"
    @test !haskey(merged, "\$ref")
end

@testset "collect_fields" begin
    schema = load_schema("problem")
    fields = collect_fields(schema, schema, String[])
    names = [f.name for f in fields]
    @test "Type" in names
    @test "Verbose" in names
    @test "Output" in names
    @test "OutputFormats" in names
end

@testset "field metadata" begin
    schema = load_schema("problem")
    fields = collect_fields(schema, schema, String[])
    type_field = first(f for f in fields if f.name == "Type")
    @test type_field.required == true
    @test !isempty(type_field.enum_values)
    @test "Eigenmode" in type_field.enum_values
    @test haskey(type_field.enum_descriptions, "Eigenmode")
    verbose_field = first(f for f in fields if f.name == "Verbose")
    @test verbose_field.default == "1"
    @test verbose_field.required == false
end

@testset "sub_schemas populated for nested object" begin
    schema = load_schema("problem")
    fields = collect_fields(schema, schema, String[])
    fmt_field = first(f for f in fields if f.name == "OutputFormats")
    @test length(fmt_field.sub_schemas) == 1
    @test fmt_field.sub_schemas[1].first == "OutputFormats"
    @test haskey(fmt_field.sub_schemas[1].second, "properties")
end

@testset "path propagated correctly" begin
    schema = load_schema("problem")
    fields = collect_fields(schema, schema, ["Problem"])
    verbose = first(f for f in fields if f.name == "Verbose")
    @test verbose.path == ["Problem", "Verbose"]
end

@testset "extract_discriminator_key handles plain const and oneOf+const" begin
    # Plain const (pre-enrichment style)
    plain = OrderedDict{String, Any}(
        "properties" => OrderedDict{String, Any}(
            "Type" => OrderedDict{String, Any}("const" => "Linear")
        )
    )
    @test extract_discriminator_key(plain, "fallback") == "Linear"

    # oneOf+const (post-enrichment style)
    enriched = OrderedDict{String, Any}(
        "properties" => OrderedDict{String, Any}(
            "Type" => OrderedDict{String, Any}(
                "oneOf" => [
                    OrderedDict{String, Any}(
                        "const" => "Linear",
                        "description" => "Linear spacing"
                    ),
                    OrderedDict{String, Any}(
                        "const" => "Log",
                        "description" => "Log spacing"
                    )
                ]
            )
        )
    )
    @test extract_discriminator_key(enriched, "fallback") == "Linear"

    # No Type field → fallback
    no_type = OrderedDict{String, Any}("properties" => OrderedDict{String, Any}())
    @test extract_discriminator_key(no_type, "fallback") == "fallback"
end

@testset "render_fields_block! basic" begin
    schema = load_schema("problem")
    fields = collect_fields(schema, schema, ["Problem"])
    leaf = filter(f -> isempty(f.sub_schemas), fields)
    buf = IOBuffer()
    render_fields_block!(buf, leaf)
    md = String(take!(buf))
    # One @raw html block for all leaf fields
    @test count("```@raw html", md) == 1
    @test contains(md, "palace-config")
    # Verbose field
    @test contains(md, "id=\"config-problem-verbose\"")
    @test contains(md, "href=\"#config-problem-verbose\"")
    @test contains(md, "config-type")
    @test contains(md, "integer")
    @test contains(md, "config-default")
    @test contains(md, "default:")
    @test contains(md, "<code>1</code>")
    # Type field (required, enum)
    @test contains(md, "id=\"config-problem-type\"")
    @test contains(md, "config-required")
    @test contains(md, "config-enum")
    @test contains(md, "\"Eigenmode\"")
    @test contains(md, "\"Driven\"")
end

@testset "render_schema_section! produces valid markdown" begin
    schema = load_schema("problem")
    buf = IOBuffer()
    render_schema_section!(buf, schema, schema, ["Problem"], 1)
    md = String(take!(buf))
    @test contains(md, "@id config-problem")
    @test contains(md, "id=\"config-problem-type\"")
    @test contains(md, "id=\"config-problem-verbose\"")
    @test contains(md, "## ")
    @test contains(md, "@id config-problem-outputformats")
    @test contains(md, "Key path:")
    # All leaf fields in one block, not one block per field
    @test count("```@raw html", md) < count("id=\"config-", md)
end

@testset "breadcrumb uses exact key names" begin
    @test breadcrumb(["Solver", "Eigenmode"]) == "solver[\"Eigenmode\"]"
    @test breadcrumb(["Boundaries", "LumpedPort"]) == "boundaries[\"LumpedPort\"]"
    @test breadcrumb(String[]) == ""
end

@testset "anchor_id lowercases all segments" begin
    @test anchor_id(["Solver", "Eigenmode", "Target"]) == "config-solver-eigenmode-target"
    @test anchor_id(["Boundaries", "LumpedPort"]) == "config-boundaries-lumpedport"
    @test anchor_id(String[]) == "config"  # no trailing dash
end

@testset "generate_file writes output" begin
    tmp = mktempdir()
    generate_file("problem"; output_dir=tmp)
    outfile = joinpath(tmp, "problem.md")
    @test isfile(outfile)
    content = read(outfile, String)
    @test contains(content, "Amazon.com")
    @test contains(content, "id=\"config-problem-verbose\"")
    @test contains(content, "id=\"config-problem-type\"")
    @test contains(content, "@id config-problem")
    @test !contains(content, "id=\"config-\"")
end

@testset "collect_toc_entries" begin
    schema = load_schema("problem")
    entries = collect_toc_entries(schema, schema, ["Problem"])
    @test entries[1].title == "Problem"
    @test entries[1].id == "config-problem"
    @test entries[1].depth == 0
    # OutputFormats is a sub-section
    sub = filter(e -> e.depth == 1, entries)
    @test any(e -> contains(e.id, "outputformats"), sub)
end

@testset "generate_all writes reference.md" begin
    tmp = mktempdir()
    generate_all(output_dir=tmp)
    outfile = joinpath(tmp, "reference.md")
    @test isfile(outfile)
    content = read(outfile, String)
    @test contains(content, "config-reference")
    @test contains(content, "config-toc")
    # All five sections present
    for name in ["problem", "model", "domains", "boundaries", "solver"]
        @test contains(content, "config-$(name)")
    end
    # TOC links to sections
    @test contains(content, "href=\"#config-problem\"")
    @test contains(content, "href=\"#config-solver\"")
    # Field entries present
    @test contains(content, "id=\"config-problem-type\"")
end

@testset "collect_fields preserves JSON order within required/optional groups" begin
    schema = load_schema("solver")
    fields = collect_fields(schema, schema, String[])
    names = [f.name for f in fields]
    order_idx = i -> findfirst(==(i), names)
    @test order_idx("Order") < order_idx("PartialAssemblyOrder")
end
