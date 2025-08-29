# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

module ValidateConfig

using JSON
using JSONSchema

export validate_config_file

"""
    validate_config_file(; config_file::AbstractString, schema_file::AbstractString)

Validate the provided configuration file against the provided JSON Schema
"""
function validate_config_file(; config_file::AbstractString, schema_file::AbstractString)
    # Resolve relative paths in schema by running parsefile in its directory
    schema = cd(() -> Schema(JSON.parsefile(schema_file)), dirname(schema_file))
    config = JSON.parse(replace_range_expand(preprocess_string(read(config_file, String))))
    return validate(schema, config)
end

"""
    preprocess_string(str::AbstractString)

Strip C and C++ style comments (//, /* */) from a file using regex, as well as whitespace
and erroneous trailing commas

Correctly ignores whitespace and comments within strings

See tinyurl.com/2s3n8dkr
"""
function preprocess_string(str::AbstractString)
    rgx1 = Regex(
        raw"(([\"'])(?:(?=(\\?))\3.)*?\2)" * raw"|(\/\*(.|[\r\n])*?\*\/)" * raw"|(\/\/.*)"
    )
    rgx2 = Regex(raw"(([\"'])(?:(?=(\\?))\3.)*?\2)" * raw"|(\s+)")
    rgx3 = Regex(raw"(([\"'])(?:(?=(\\?))\3.)*?\2)" * raw"|,+(?=\s*?[\}\]])")
    return replace(
        replace(replace(str, rgx1 => s"\g<1>"), rgx2 => s"\g<1>"),
        rgx3 => s"\g<1>"
    )
end

"""
    replace_range_expand(str::AbstractString)

Replace integer ranges with lists
"""
function replace_range_expand(str::AbstractString)
    function helper(range::AbstractString)
        # Expand the range and handle negative integers too
        res = ""
        for r in eachsplit(range, ",")
            i = findfirst('-', r[2:end])
            if i !== nothing
                j = i + 1
                r0 = parse(Int, r[1:(j - 1)])
                r1 = parse(Int, r[(j + 1):end])
                @assert r0 < r1 "Invalid range bounds in range expansion!"
                for q = r0:r1
                    res *= string(q) * ","
                end
            else
                res *= r * ","
            end
        end
        if last(res) == ','
            res = res[1:(end - 1)]
        end
        return "[" * res * "]"
    end
    # Unfortunate double match below but OK
    # (https://discourse.julialang.org/t/using-replace-with-a-function-of-the-match/41264)
    rgx = r"\[(-?[0-9][\-\,0-9]*[0-9])\]"
    return replace(str, rgx => s -> helper(match(rgx, s).captures[1]))
end

end
