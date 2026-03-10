# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Thin CLI argument abstraction with env variable fallbacks. (Note for Julia
# experts, this struct is not concrete, but we do not care about performance
# here.)
struct ArgConfig
    "Name of the command-line argument"
    name::String
    "Corresponding environment variable"
    env_var::String
    "Default value if there is not environment variable"
    default::Any
    "Description for -h"
    description::String
    """Function that parses the string to whatever is needed for the rest of
    the execution of the script"""
    parser::Function

    # Keyword argument based constructor.
    ArgConfig(; name, env_var, default, description, parser=identity) =
        new(name, env_var, default, description, parser)
end

function parse_args(configs::Vector{ArgConfig})
    args = Dict{String, Any}()

    # Show help if requested.
    if "--help" in ARGS || "-h" in ARGS
        println("Usage: julia runtests.jl [options]")
        println("\nOptions:")
        for config in configs
            println("  --$(config.name)  $(config.description)")
            println("    Environment variable: $(config.env_var)")
            println("    Default: $(config.default)")
            println()
        end
        exit(0)
    end

    # Check for unrecognized flags
    valid_flags = Set("--$(config.name)" for config in configs)
    push!(valid_flags, "--help", "-h")

    for arg in ARGS
        if startswith(arg, "--") && arg ∉ valid_flags
            error("Unrecognized flag: $arg\nUse --help to see available options.")
        end
    end

    # Parse each argument.
    for config in configs
        cli_arg = "--$(config.name)"

        # Check CLI args first (highest precedence).
        if cli_arg in ARGS
            idx = findfirst(==(cli_arg), ARGS)
            if idx !== nothing && idx < length(ARGS)
                # For test-cases, collect all following non-flag arguments
                # (needed to support syntax like `--test-cases spheres rings`
                # instead of `--test-cases "spheres rings"`)
                if config.name == "test-cases"
                    values = String[]
                    for i = (idx + 1):length(ARGS)
                        if startswith(ARGS[i], "--")
                            break
                        end
                        push!(values, ARGS[i])
                    end
                    args[config.name] = isempty(values) ? config.default : values
                else
                    args[config.name] = config.parser(ARGS[idx + 1])
                end
                continue
            end
        end

        # Check environment variable (medium precedence).
        if config.env_var in keys(ENV)
            args[config.name] = config.parser(ENV[config.env_var])
            continue
        end

        # Use default (lowest precedence).
        args[config.name] = config.default
    end

    return args
end
