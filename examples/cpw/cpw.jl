# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using CSV, DataFrames, Measures, CairoMakie, ColorSchemes

"""
    cpw_impedance(;w,s,h,ϵᵣ)

Compute the characteristic impedance of a coplanar wave guide

See p. 259 of H. J. Visser, _Antenna Theory and Applications_, Wiley, Hoboken, NJ, 2012

# Arguments

    - w width of trace [μm]
    - s separation of trace [μm]
    - h height of substrate [μm]
    - ϵᵣ relative permittivity in surface normal direction
"""
function cpw_impedance(; w=30, s=18, h=500, ϵᵣ)
    """
        Helper function for K(k)/K(k')
    """
    function KoverK′(k, k′)
        s =
            x ->
                log(2 * (sqrt(1 + x) + (4 * x)^(1 // 4)) / (sqrt(1 + x) - (4 * x)^(1 // 4)))
        if k >= 1.0 / sqrt(2)
            return s(k) / (2 * π)
        else
            return 2 * π / s(k′)
        end
    end

    k = w / (w + 2 * s)
    k₁ = sinh(π * w / (4 * h)) / sinh(π * (w + 2 * s) / (4 * h))

    k′ = sqrt(1 - k^2)
    k₁′ = sqrt(1 - k₁^2)

    koverk′ = KoverK′(k, k′)
    k₁overk₁′ = KoverK′(k₁, k₁′)

    ϵ_eff = 1 + ((ϵᵣ - 1) / 2) * k₁overk₁′ / koverk′

    return Z₀ = 30 * π / (koverk′ * sqrt(ϵ_eff))
end

"""
    extract_data(base_path::String, cols::Vector{Int}=[1,2,4,6,8])

Extract data from CSV files located in specific subfolders of `base_path`.

# Arguments

  - `base_path::String`: Path to the parent directory containing the subfolders
  - `cols::Vector{Int}`: Indices of columns to extract (default: [1,2,4,6,8])

# Returns

A nested dictionary structure where:

  - Outer key: folder name (e.g., "wave_uniform", "lumped_uniform", etc.)
  - Inner keys: column names from CSV headers
  - Values: Vector{Float64} of data from each column

# Example

```julia
data = extract_data("postpro", [1, 2, 3])
freq = data["lumped_uniform"]["f (GHz)"]
s11_mag = data["lumped_uniform"]["|S[1][1]| (dB)"]    # Define folder structure
```
"""
function extract_data(base_path::String, cols::Vector{Int}=[1, 2, 4, 6, 8])
    # Define folder structure
    folders = ["wave_uniform", "lumped_uniform", "wave_adaptive", "lumped_adaptive"]

    # Initialize result dictionary
    result = Dict{String, Dict{String, Vector{Float64}}}()

    for folder in folders
        try
            # Construct full file path
            file_path = joinpath(base_path, folder, "port-S.csv")

            # Verify file exists
            if !isfile(file_path)
                @warn "File not found: $file_path"
                continue
            end

            # Read CSV file
            df = CSV.read(file_path, DataFrame)

            # Get column names and strip leading whitespace
            col_names = [lstrip(String(name)) for name in names(df)]

            # Validate column indices
            valid_cols = filter(c -> c ≤ length(col_names), cols)

            if isempty(valid_cols)
                @warn "No valid columns found for $folder"
                continue
            end

            # Initialize inner dictionary
            result[folder] = Dict{String, Vector{Float64}}()

            # Extract specified columns
            for col_idx in valid_cols
                if col_idx ≤ length(col_names)
                    col_name = col_names[col_idx]
                    result[folder][col_name] = Vector{Float64}(df[!, col_idx])
                end
            end
        catch e
            @error "Error processing $folder" exception=(e, catch_backtrace())
        end
    end

    return result
end

"""
    plot_s_parameters(path::String, prefix::String="")

Create and save individual S-parameter plots from data in specified path.
Plots S11, S21, S31, and S41 using Dark2_4 colorscheme, with circles for uniform
and lines for adaptive solutions.

# Arguments

  - `path::String`: Path containing the data folders and where plots will be saved
  - `prefix::String`: Optional prefix for saved files (default="")

# Saves

Four PNG files in the specified path:

  - cpw_[prefix_]11.png
  - cpw_[prefix_]21.png
  - cpw_[prefix_]31.png
  - cpw_[prefix_]41.png
"""
function plot_s_parameters(path::String, prefix::String="")
    data = extract_data(path)

    s_params = ["S[1][1]", "S[2][1]", "S[3][1]", "S[4][1]"]
    unicode_labels = ["|S₁₁| dB", "|S₂₁| dB", "|S₃₁| dB", "|S₄₁| dB"]
    colors = ColorSchemes.Dark2_4

    for (s_param, label) in zip(s_params, unicode_labels)
        fig = Figure(size=(600, 400))
        ax = Axis(fig[1, 1])

        scatter!(
            ax,
            data["wave_uniform"]["f (GHz)"],
            data["wave_uniform"]["|$(s_param)| (dB)"],
            color=colors[1],
            markersize=10,
            label="Wave Uniform"
        )

        scatter!(
            ax,
            data["lumped_uniform"]["f (GHz)"],
            data["lumped_uniform"]["|$(s_param)| (dB)"],
            color=colors[2],
            markersize=10,
            label="Lumped Uniform"
        )

        lines!(
            ax,
            data["wave_adaptive"]["f (GHz)"],
            data["wave_adaptive"]["|$(s_param)| (dB)"],
            color=colors[1],
            linewidth=2,
            label="Wave Adaptive"
        )

        lines!(
            ax,
            data["lumped_adaptive"]["f (GHz)"],
            data["lumped_adaptive"]["|$(s_param)| (dB)"],
            color=colors[2],
            linewidth=2,
            label="Lumped Adaptive"
        )

        ax.xlabel = "Frequency (GHz)"
        ax.ylabel = label

        ax.xlabelsize = 18
        ax.ylabelsize = 18
        ax.xticklabelsize = 18
        ax.yticklabelsize = 18

        data_min =
            minimum(minimum(data[folder]["|$(s_param)| (dB)"]) for folder in keys(data))
        data_max =
            maximum(maximum(data[folder]["|$(s_param)| (dB)"]) for folder in keys(data))
        y_min = floor(data_min / 10) * 10
        y_max = ceil(data_max / 10) * 10
        y_max == 0 ? y_max = 1 : y_max = y_max
        ax.limits = (nothing, (y_min, y_max))

        axislegend(ax, position=:lb, textsize=18)

        suffix = s_param[[3, 6]]
        filename = isempty(prefix) ? "cpw-$(suffix).png" : "cpw-$(prefix)-$(suffix).png"

        save(joinpath(path, filename), fig)
    end
end
