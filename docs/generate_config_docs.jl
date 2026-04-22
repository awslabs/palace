# Auto-generates docs/src/config/*.md from scripts/schema/config/*.json
# Run via: julia --project docs/generate_config_docs.jl
# Or included in docs/make.jl

using JSON
using OrderedCollections
using CommonMark
using JuliaFormatter

const SCHEMA_DIR = joinpath(@__DIR__, "..", "scripts", "schema", "config")
const ROOT_SCHEMA_PATH = joinpath(@__DIR__, "..", "scripts", "schema", "config-schema.json")
const OUTPUT_DIR = joinpath(@__DIR__, "src", "config")

# Use OrderedDict throughout so JSON key order is preserved end-to-end.
# JSON.parsefile with dicttype=OrderedDict returns OrderedDict{String,Any} at top level
# but nested dicts come back as OrderedDict{Any,Any}. We use AbstractDict as the type
# parameter throughout so all levels are handled uniformly.
const SchemaDict = AbstractDict

# --- Data types ---

struct FieldDoc
    name::String              # exact JSON key, e.g. "Target"
    title::String             # human-readable name from "title" field, falls back to name
    path::Vector{String}      # full key path from root, e.g. ["Solver", "Eigenmode", "Target"]
    type_str::String          # rendered type string
    default::Union{String, Nothing}
    required::Bool
    description::String
    enum_values::Vector{String}
    enum_descriptions::Dict{String, String}
    constraints::String
    advanced::Bool
    deprecated::Bool
    is_array::Bool            # true when this field is an array-of-objects (sub_schemas are item schemas)
    sub_schemas::Vector{Pair{String, AbstractDict}}  # (title_or_name => resolved_schema)
end

# --- Schema loading ---

"""
Load a top-level schema dict by section name (e.g. "problem").
Uses dicttype=OrderedDict to preserve JSON key order throughout.
"""
function load_schema(name::String)::AbstractDict
    path = joinpath(SCHEMA_DIR, "$name.json")
    return JSON.parsefile(path; dicttype=OrderedDict)
end

# --- $ref resolution ---

"""
Resolve a local \$ref string (e.g. "#/\$defs/Attributes") against the root schema.
"""
function resolve_ref(root::AbstractDict, ref::String)::AbstractDict
    startswith(ref, "#/") || error("Only local \$ref supported, got: $ref")
    parts = split(ref[3:end], "/")
    node = root
    for part in parts
        node = node[part]
    end
    return node
end

"""
Return a schema node with any top-level \$ref merged in.
Keys on the inline schema override keys from the \$ref target.
"""
function inline_ref(root::AbstractDict, schema::AbstractDict)::AbstractDict
    haskey(schema, "\$ref") || return schema
    base = resolve_ref(root, schema["\$ref"])
    merged = merge(base, schema)
    delete!(merged, "\$ref")
    return merged
end

# --- Type string ---

"""
Produce a compact human-readable type string from a schema node.
"""
function type_string(root::AbstractDict, schema::AbstractDict)::String
    schema = inline_ref(root, schema)

    if haskey(schema, "oneOf")
        branches = schema["oneOf"]
        # Const-enum: all branches have "const", none have "properties"
        if all(haskey(b, "const") && !haskey(b, "properties") for b in branches)
            return branches[1]["const"] isa Number ? "number" : "string"
        end
        # Validation-only oneOf: branches carry only "required" constraints with no type
        # info (e.g. "must have Attributes OR Elements"). Ignore for type resolution and
        # fall through to the explicit "type" field below.
        validation_only = all(
            !haskey(b, "const") && !haskey(b, "type") && !haskey(b, "properties") for
            b in branches
        )
        if !validation_only
            return join(unique(type_string(root, b) for b in branches), " or ")
        end
    end

    if haskey(schema, "anyOf")
        return join(unique(type_string(root, b) for b in schema["anyOf"]), " or ")
    end

    t = get(schema, "type", "any")
    if t == "array"
        items_schema = get(schema, "items", OrderedDict{String, Any}())
        inner = type_string(root, items_schema)
        min_i = get(schema, "minItems", nothing)
        max_i = get(schema, "maxItems", nothing)
        if min_i !== nothing && min_i == max_i
            return "[$inner × $min_i]"
        end
        return "[$inner, ...]"
    end
    return t
end

# --- Constraint string ---

function constraints_string(schema::AbstractDict)::String
    parts = String[]
    for (key, sym) in [
        ("minimum", "≥"),
        ("exclusiveMinimum", ">"),
        ("maximum", "≤"),
        ("exclusiveMaximum", "<")
    ]
        haskey(schema, key) && push!(parts, "$sym $(schema[key])")
    end
    return join(parts, ", ")
end

# --- Schema node classification ---

"""
Extract the discriminator key for a variant schema in a discriminated union.
Handles both plain const ({const: "Linear"}) and oneOf+const after enrichment
({oneOf: [{const: "Linear", description: "..."}, ...]}).
Falls back to `fallback` if no discriminator is found.
"""
function extract_discriminator_key(variant::AbstractDict, fallback::String)::String
    type_field = get(
        get(variant, "properties", OrderedDict{String, Any}()),
        "Type",
        OrderedDict{String, Any}()
    )
    haskey(type_field, "const") && return string(type_field["const"])
    branches = get(type_field, "oneOf", Any[])
    if !isempty(branches) && haskey(first(branches), "const")
        return string(first(branches)["const"])
    end
    return fallback
end

"""
Return true if this oneOf is a discriminated-union of objects, not a const-enum
and not a validation-only oneOf (branches with only "required" constraints).
"""
function is_object_dispatch(branches::Vector)::Bool
    isempty(branches) && return false
    all(haskey(b, "const") for b in branches) && return false
    all(!haskey(b, "properties") && !haskey(b, "type") for b in branches) && return false
    return any(haskey(b, "properties") for b in branches)
end

# --- Field collection ---

"""
Collect FieldDoc entries for all properties at one schema level.
Preserves JSON key order (OrderedDict). Required fields come before optional ones
within the preserved order.
"""
function collect_fields(
    schema::AbstractDict,
    root::AbstractDict,
    path::Vector{String}
)::Vector{FieldDoc}
    props = get(schema, "properties", OrderedDict{String, Any}())
    required_names = Set(get(schema, "required", String[]))
    results = FieldDoc[]

    for (name, raw_field) in props
        field = inline_ref(root, raw_field)
        fpath = vcat(path, [name])

        ftitle = get(field, "title", name)
        ftype = type_string(root, field)
        fdefault = haskey(field, "default") ? JSON.json(field["default"]) : nothing
        fdesc = get(field, "description", "")
        fconstraints = constraints_string(field)
        fadvanced = get(field, "x-palace-advanced", false)
        fdeprecated = get(field, "x-palace-deprecated", false)
        fenum = String[]
        fenum_desc = Dict{String, String}()

        # Detect oneOf+const enum
        oneOf_branches = get(field, "oneOf", Any[])
        is_const_enum =
            !isempty(oneOf_branches) &&
            all(haskey(b, "const") && !haskey(b, "properties") for b in oneOf_branches)
        if is_const_enum
            for branch in oneOf_branches
                val = string(branch["const"])
                push!(fenum, val)
                if haskey(branch, "description")
                    fenum_desc[val] = string(branch["description"])
                end
            end
        end

        # Collect child sub-schemas that become sub-sections
        sub_schemas = Pair{String, AbstractDict}[]
        fis_array = false

        t = get(field, "type", "")
        if t == "object" && haskey(field, "properties")
            push!(sub_schemas, name => field)
        elseif t == "array" && haskey(field, "items")
            items = inline_ref(root, field["items"])
            if get(items, "type", "") == "object" && haskey(items, "properties")
                fis_array = true
                # Use item title as key; empty string signals "unnamed single-item array"
                item_title = get(items, "title", "")
                push!(sub_schemas, item_title => items)
            elseif haskey(items, "oneOf") && is_object_dispatch(items["oneOf"])
                fis_array = true
                for variant in items["oneOf"]
                    resolved = inline_ref(root, variant)
                    haskey(resolved, "properties") || continue
                    variant_key = extract_discriminator_key(resolved, name)
                    push!(sub_schemas, variant_key => resolved)
                end
            end
        elseif !is_const_enum &&
               haskey(field, "oneOf") &&
               is_object_dispatch(field["oneOf"])
            for variant in field["oneOf"]
                resolved = inline_ref(root, variant)
                haskey(resolved, "properties") || continue
                variant_key = extract_discriminator_key(resolved, name)
                push!(sub_schemas, variant_key => resolved)
            end
        end

        push!(
            results,
            FieldDoc(
                name,
                ftitle,
                fpath,
                ftype,
                fdefault,
                name in required_names,
                fdesc,
                fenum,
                fenum_desc,
                fconstraints,
                fadvanced,
                fdeprecated,
                fis_array,
                sub_schemas
            )
        )
    end

    # Stable sort: required fields before optional, preserving relative JSON order within groups
    required_fields = filter(f -> f.required, results)
    optional_fields = filter(f -> !f.required, results)
    return vcat(required_fields, optional_fields)
end

# --- Anchor IDs ---

"""
Build a Documenter.jl anchor ID from a JSON key path.
All segments are lowercased and joined with dashes, prefixed with "config-".
Empty path returns "config" (no trailing dash).
"""
function anchor_id(path::Vector{String})::String
    isempty(path) && return "config"
    return "config-" * join(lowercase.(path), "-")
end

"""
Render a clickable JSON Pointer breadcrumb as an inline @raw html snippet.

path        — key path segments, e.g. ["Solver","Eigenmode"]
array_index — when >= 0, appended as a final /N segment (RFC 6901 array element notation).
For a named-variant path the last segment is a discriminator key rather than
a real JSON key, so drop_last=true strips it before appending /N.
array_anchor — when non-empty, the /N index segment is rendered as a link to this anchor.
Used for named-variant sections so the index links back to the variant heading.

Examples:
path=["Solver","Eigenmode"],   array_index=-1  →  /Solver/Eigenmode
path=["Boundaries","Impedance"], array_index=0  →  /Boundaries/Impedance/0  (plain text)
path=["Solver","Driven","Samples","Point"], array_index=1, drop_last=true,
array_anchor="config-solver-driven-samples-point"  →  /Solver/Driven/Samples/<a>1</a>

Each segment is a hyperlink to the anchor of its prefix path.
"""
function render_keypath(
    path::Vector{String};
    array_index::Int=-1,
    drop_last::Bool=false,
    array_anchor::String=""
)::String
    real_path = drop_last ? path[1:(end - 1)] : path
    buf = IOBuffer()
    println(buf, "```@raw html")
    print(buf, "<p class=\"config-keypath\"><em>Path:</em> ")
    for (i, seg) in enumerate(real_path)
        id = anchor_id(real_path[1:i])
        print(buf, "<code>/</code><a href=\"#$id\"><code>$(html_escape(seg))</code></a>")
    end
    if array_index >= 0
        if !isempty(array_anchor)
            print(
                buf,
                "<code>/</code><a href=\"#$array_anchor\"><code>$(array_index)</code></a>"
            )
        else
            print(buf, "<code>/$(array_index)</code>")
        end
    end
    println(buf, "</p>")
    println(buf, "```")
    return String(take!(buf))
end

# --- HTML helpers ---

html_escape(s::String) =
    replace(s, "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")

const _CM_PARSER = Parser()

"""
Render a description string as inline HTML.
Supports full CommonMark syntax (`code`, **bold**, [text](url)) plus
Documenter cross-references: [text](@ref anchor-id) → <a href="#anchor-id">text</a>.
Returns inner HTML with surrounding <p> tags stripped.
"""
function render_description(s::String)::String
    isempty(s) && return ""
    md = replace(s, r"\[([^\]]+)\]\(@ref ([^\)]+)\)" => s"[\1](#\2)")
    buf = IOBuffer()
    html(buf, _CM_PARSER(md))
    result = String(take!(buf))
    return replace(result, r"^<p>(.*)</p>\n?$"s => s"\1")
end

"""
Sanitize a description for plain Markdown output (section body text).
Converts Documenter cross-references [text](@ref anchor) → [text](#anchor)
so Documenter does not try to resolve them against the Julia docstring index.
"""
sanitize_md(s::String) = replace(s, r"\[([^\]]+)\]\(@ref ([^\)]+)\)" => s"[\1](#\2)")

# --- Markdown rendering ---

const COPYRIGHT_HEADER = """
```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

"""

"""
Render one <dt>/<dd> pair into buf (no surrounding <dl> tags).
For fields with sub-schemas (objects/arrays of objects), only a link to the
sub-section is rendered — no description, type badge, or default.
"""
function render_field_entry!(buf::IOBuffer, f::FieldDoc)
    id = anchor_id(f.path)

    if isempty(f.sub_schemas)
        # Leaf: id on dt + self-link
        print(buf, "  <dt id=\"$id\">")
        print(buf, "<a href=\"#$id\"><code>\"$(html_escape(f.name))\"</code></a>")
    else
        # Sub-section: no id here — the anchor lives on the heading below
        print(buf, "  <dt>")
        print(buf, "<code>\"$(html_escape(f.name))\"</code>")
    end
    # All fields get the full badge row
    print(buf, " <span class=\"config-type\">$(html_escape(f.type_str))</span>")
    if f.required
        print(buf, " <span class=\"config-required\">required</span>")
    else
        default_val = something(f.default, "—")
        print(
            buf,
            " <span class=\"config-default\">default: <code>$(html_escape(default_val))</code></span>"
        )
    end
    if !isempty(f.constraints)
        print(
            buf,
            " <span class=\"config-constraint\"><code>$(html_escape(f.constraints))</code></span>"
        )
    end
    if f.advanced
        print(buf, " <span class=\"config-advanced\">advanced</span>")
    end
    if f.deprecated
        print(buf, " <span class=\"config-deprecated\">deprecated</span>")
    end
    println(buf, "</dt>")

    println(buf, "  <dd>")
    if isempty(f.sub_schemas)
        # Leaf: description + enum values
        if !isempty(f.description)
            println(buf, "    <p>$(render_description(f.description))</p>")
        end
        if !isempty(f.enum_values)
            println(buf, "    <dl class=\"config-enum\">")
            for v in f.enum_values
                println(buf, "      <dt><code>\"$(html_escape(v))\"</code></dt>")
                desc = get(f.enum_descriptions, v, "")
                println(buf, "      <dd>$(render_description(desc))</dd>")
            end
            println(buf, "    </dl>")
        end
    else
        # Sub-section: link to the heading
        println(buf, "    <p><a href=\"#$id\">See full reference ↓</a></p>")
    end
    return println(buf, "  </dd>")
end

"""
Emit a single @raw html block containing one <dl> for all leaf fields.
"""
function render_fields_block!(buf::IOBuffer, fields::Vector{FieldDoc})
    isempty(fields) && return
    println(buf, "```@raw html")
    println(buf, "<dl class=\"palace-config\">")
    for f in fields
        render_field_entry!(buf, f)
    end
    println(buf, "</dl>")
    println(buf, "```")
    return println(buf)
end

# --- TOC collection ---

"""
A single entry in the table of contents: display title, anchor id, indent depth.
"""
struct TocEntry
    title::String
    id::String
    depth::Int   # 0 = top-level section, 1 = subsection
end

"""
Walk the schema tree and collect TocEntry values for the TOC.
Collects the top-level section (depth=0) and all direct H2 subsections (depth=1).
Discriminated-union parent headings count as depth=1; their variants are skipped
(too granular for a TOC).
"""
function collect_toc_entries(
    schema::AbstractDict,
    root::AbstractDict,
    path::Vector{String}
)::Vector{TocEntry}
    entries = TocEntry[]
    display_title = get(schema, "title", isempty(path) ? "Configuration" : path[end])
    push!(entries, TocEntry(display_title, anchor_id(path), 0))

    fields = collect_fields(schema, root, path)
    for f in filter(f -> !isempty(f.sub_schemas), fields)
        if !f.is_array
            _, sub_schema = f.sub_schemas[1]
            sub_title = get(sub_schema, "title", f.name)
            push!(entries, TocEntry(sub_title, anchor_id(f.path), 1))
        else
            raw_field = get(schema["properties"], f.name, OrderedDict{String, Any}())
            parent_title = get(inline_ref(root, raw_field), "title", f.name)
            push!(entries, TocEntry(parent_title, anchor_id(f.path), 1))
        end
    end
    return entries
end

"""
Render a TOC as a single @raw html block.
Each section is a bold link on its own line; subsections follow on the next line
as inline links separated by " · ".
"""
function render_toc(all_entries::Vector{Vector{TocEntry}})::String
    buf = IOBuffer()
    println(buf, "```@raw html")
    println(buf, "<dl class=\"config-toc\">")
    for section_entries in all_entries
        top = section_entries[1]
        subs = section_entries[2:end]
        println(buf, "  <dt><a href=\"#$(top.id)\">$(html_escape(top.title))</a></dt>")
        if !isempty(subs)
            sub_links = join(
                ["<a href=\"#$(s.id)\">$(html_escape(s.title))</a>" for s in subs],
                " · "
            )
            println(buf, "  <dd>$sub_links</dd>")
        end
    end
    println(buf, "</dl>")
    println(buf, "```")
    println(buf)
    return String(take!(buf))
end

"""
Emit the badge span(s) for a FieldDoc as an inline @raw html snippet.
"""
function render_field_badges(f::FieldDoc)::String
    buf = IOBuffer()
    println(buf, "```@raw html")
    print(buf, "<span class=\"config-section-badges\">")
    print(buf, "<span class=\"config-type\">$(html_escape(f.type_str))</span>")
    if f.required
        print(buf, " <span class=\"config-required\">required</span>")
    else
        default_val = something(f.default, "—")
        print(
            buf,
            " <span class=\"config-default\">default: <code>$(html_escape(default_val))</code></span>"
        )
    end
    if !isempty(f.constraints)
        print(
            buf,
            " <span class=\"config-constraint\"><code>$(html_escape(f.constraints))</code></span>"
        )
    end
    if f.advanced
        print(buf, " <span class=\"config-advanced\">advanced</span>")
    end
    if f.deprecated
        print(buf, " <span class=\"config-deprecated\">deprecated</span>")
    end
    println(buf, "</span>")
    println(buf, "```")
    return String(take!(buf))
end

# --- Recursive section renderer ---

"""
Render a schema object node as a Markdown section, recursing into all sub-schemas.
All leaf fields at each level are emitted in a single @raw html <dl> block.
"""
function render_schema_section!(
    buf::IOBuffer,
    schema::AbstractDict,
    root::AbstractDict,
    path::Vector{String},
    heading_level::Int,
    field_doc::Union{FieldDoc, Nothing}=nothing;
    array_index::Int=-1   # >= 0 for named variant sections: last path segment is a
    # discriminator key (not a real JSON key), replaced by /N index
)
    hdr = "#" ^ heading_level
    display_title = get(schema, "title", isempty(path) ? "Configuration" : path[end])
    id = anchor_id(path)

    println(buf, "$hdr [$display_title](@id $id)\n")
    if !isempty(path)
        drop = array_index >= 0
        # For named-variant sections, link the /N index to this section's own anchor
        anchor = array_index >= 0 ? id : ""
        print(
            buf,
            render_keypath(
                path;
                array_index=array_index,
                drop_last=drop,
                array_anchor=anchor
            )
        )
        println(buf)
    end
    if field_doc !== nothing
        print(buf, render_field_badges(field_doc))
        println(buf)
    end
    # Prefer the field-level description (covers array items where schema is the items object)
    desc =
        (field_doc !== nothing && !isempty(field_doc.description)) ? field_doc.description :
        get(schema, "description", "")
    if !isempty(desc)
        println(buf, "$(sanitize_md(desc))\n")
    end

    fields            = collect_fields(schema, root, path)
    subsection_fields = filter(f -> !isempty(f.sub_schemas), fields)

    # Render all fields (leaves get full badge row; sub-section fields get a link)
    render_fields_block!(buf, fields)

    return render_subsections!(buf, subsection_fields, schema, root, heading_level)
end

"""
Render sub-sections for a list of FieldDocs that have sub_schemas.
Extracted so both render_schema_section! and the unnamed-array-item path share
the same array-aware dispatch logic.
"""
function render_subsections!(
    buf::IOBuffer,
    subsection_fields::Vector{FieldDoc},
    parent_schema::AbstractDict,
    root::AbstractDict,
    heading_level::Int
)
    for f in subsection_fields
        if !f.is_array
            # Plain object field: one section at the field's path
            _, sub_schema = f.sub_schemas[1]
            render_schema_section!(buf, sub_schema, root, f.path, heading_level + 1, f)
        else
            # Array-of-objects field: render the array section first (badges, description),
            # then render each item variant as a sub-section.
            array_id = anchor_id(f.path)
            raw_field = get(parent_schema["properties"], f.name, OrderedDict{String, Any}())
            parent_field = inline_ref(root, raw_field)
            array_title = get(parent_field, "title", f.name)
            sub_hdr = "#" ^ (heading_level + 1)
            println(buf, "$sub_hdr [$array_title](@id $array_id)\n")
            print(buf, render_keypath(f.path))
            println(buf)
            print(buf, render_field_badges(f))
            println(buf)
            if !isempty(f.description)
                println(buf, "$(sanitize_md(f.description))\n")
            end

            item_key, item_schema = f.sub_schemas[1]
            if length(f.sub_schemas) == 1 && isempty(item_key)
                # Special case: single unnamed item type (e.g. boundaries["Impedance"]).
                # The item schema has no title in the JSON, so we don't add a new heading —
                # that would be noise for a single homogeneous array. Instead we emit the
                # item's key path (/path/0) and fields inline at the same heading level,
                # making clear these are per-element fields. Any nested sub-sections in the
                # item are recursed into via render_subsections! so they get the same treatment.
                print(buf, render_keypath(f.path; array_index=0))
                println(buf)
                item_desc = get(item_schema, "description", "")
                if !isempty(item_desc)
                    println(buf, "$(sanitize_md(item_desc))\n")
                end
                item_fields = collect_fields(item_schema, root, f.path)
                render_fields_block!(buf, item_fields)
                render_subsections!(
                    buf,
                    filter(fi -> !isempty(fi.sub_schemas), item_fields),
                    item_schema,
                    root,
                    heading_level + 1
                )
            else
                # Multiple named variants (e.g. solver["Driven"]["Samples"]):
                # each variant gets its own sub-section. The anchor uses the variant key
                # for uniqueness; the key path replaces it with a numeric index /0, /1, ...
                for (idx, (variant_key, variant_schema)) in enumerate(f.sub_schemas)
                    variant_path = vcat(f.path, [variant_key])
                    render_schema_section!(
                        buf,
                        variant_schema,
                        root,
                        variant_path,
                        heading_level + 2,
                        nothing;
                        array_index=idx - 1
                    )
                end
            end
        end
    end
end

# --- File generation ---

const SECTIONS = ["problem", "model", "domains", "boundaries", "solver"]

"""
Generate a single docs/src/config/reference.md containing all config sections.
Each schema becomes an H1 section. A TOC with top-level + subsection links is
prepended. Individual per-section files are no longer written.
"""
function generate_all(; output_dir::String=OUTPUT_DIR)
    schemas = [(name, load_schema(name)) for name in SECTIONS]
    root_schema = JSON.parsefile(ROOT_SCHEMA_PATH; dicttype=OrderedDict)
    root_required = Set(get(root_schema, "required", String[]))

    # Collect TOC entries from every schema before rendering
    all_toc = [
        collect_toc_entries(schema, schema, [get(schema, "title", titlecase(name))]) for
        (name, schema) in schemas
    ]

    buf = IOBuffer()
    print(buf, COPYRIGHT_HEADER)

    # TOC
    println(buf, "# [Configuration File Reference](@id config-reference)\n")
    print(buf, render_toc(all_toc))

    # All sections
    for (name, schema) in schemas
        root_key = get(schema, "title", titlecase(name))
        section_field = FieldDoc(
            root_key,
            root_key,
            [root_key],
            "object",
            nothing,
            root_key in root_required,
            get(schema, "description", ""),
            String[],
            Dict{String, String}(),
            "",
            false,
            false,
            false,
            Pair{String, AbstractDict}[]
        )
        render_schema_section!(buf, schema, schema, [root_key], 2, section_field)
    end

    outpath = joinpath(output_dir, "reference.md")
    write(outpath, String(take!(buf)))
    format(outpath; verbose=false, format_markdown=true)
    @info "Generated $outpath"
end

# Keep generate_file for tests (writes a single section to a temp file at H1)
function generate_file(name::String; output_dir::String=OUTPUT_DIR)
    schema = load_schema(name)
    root_key = get(schema, "title", titlecase(name))
    buf = IOBuffer()
    print(buf, COPYRIGHT_HEADER)
    render_schema_section!(buf, schema, schema, [root_key], 1)
    write(joinpath(output_dir, "$name.md"), String(take!(buf)))
    @info "Generated $(joinpath(output_dir, name * ".md"))"
end

# --- Entry point ---

if abspath(PROGRAM_FILE) == @__FILE__
    generate_all()
    @info "Done. Generated config/reference.md."
end
