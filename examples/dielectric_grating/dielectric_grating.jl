# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using CSV, DataFrames, CairoMakie, ColorSchemes

"""
    plot_s_parameters(path::String="postpro")

Plot Floquet port S-parameters from uniform (scatter) and adaptive (lines) sweeps.

Reads `port-floquet-S.csv` from `path/uniform/` and `path/adaptive/`, then saves
reflection and transmission plots to `path/`.

# Example

```
julia> include("dielectric_grating.jl")
julia> plot_s_parameters()
```
"""
function plot_s_parameters(path::String="postpro")
    uniform_file = joinpath(path, "uniform", "port-floquet-S.csv")
    adaptive_file = joinpath(path, "adaptive", "port-floquet-S.csv")

    have_uniform = isfile(uniform_file)
    have_adaptive = isfile(adaptive_file)
    if !have_uniform && !have_adaptive
        @warn "No data found. Run the simulations first."
        return nothing
    end

    data_u = have_uniform ? CSV.read(uniform_file, DataFrame) : nothing
    data_a = have_adaptive ? CSV.read(adaptive_file, DataFrame) : nothing

    freq_u, modes_u = have_uniform ? get_modes(data_u) : (nothing, Dict())
    freq_a, modes_a = have_adaptive ? get_modes(data_a) : (nothing, Dict())

    all_labels = sort(unique([keys(modes_u)..., keys(modes_a)...]))
    p1_labels = filter(l -> startswith(l, "P1"), all_labels)
    p2_labels = filter(l -> startswith(l, "P2"), all_labels)

    colors = ColorSchemes.Dark2_8
    y_min = -30.0

    for (labels, title, suffix) in [
        (p1_labels, "Reflection (Port 1)", "reflection"),
        (p2_labels, "Transmission (Port 2)", "transmission")
    ]
        fig = Figure(size=(600, 400))
        ax = Axis(fig[1, 1])
        ax.xlabel = "Frequency (GHz)"
        ax.ylabel = "|S| (dB)"
        ax.title = title
        ax.xlabelsize = 18
        ax.ylabelsize = 18
        ax.xticklabelsize = 18
        ax.yticklabelsize = 18
        ax.limits = (nothing, (y_min, 5))

        plotted = Vector{Vector{Float64}}()
        ci = 0
        for label in labels
            ref_vals =
                haskey(modes_a, label) ? modes_a[label] : get(modes_u, label, nothing)
            isnothing(ref_vals) && continue

            # Skip degenerate symmetric modes (identical S-parameter values).
            is_dup = any(plotted) do prev
                length(prev) == length(ref_vals) && all(zip(prev, ref_vals)) do (a, b)
                    return (isnan(a) && isnan(b)) || isapprox(a, b; rtol=1e-3)
                end
            end
            is_dup && continue
            push!(plotted, ref_vals)

            ci += 1
            color = colors[mod1(ci, length(colors))]

            if haskey(modes_a, label)
                vals = modes_a[label]
                mask = .!isnan.(vals) .& (vals .> y_min)
                if any(mask)
                    lines!(
                        ax,
                        freq_a[mask],
                        vals[mask];
                        color=color,
                        linewidth=2,
                        label="$label (adaptive)"
                    )
                end
            end
            if haskey(modes_u, label)
                vals = modes_u[label]
                mask = .!isnan.(vals) .& (vals .> y_min)
                if any(mask)
                    scatter!(
                        ax,
                        freq_u[mask],
                        vals[mask];
                        color=color,
                        markersize=10,
                        label="$label (uniform)"
                    )
                end
            end
        end

        axislegend(ax; position=:lb, labelsize=12, nbanks=3)
        filename = "dielectric_grating_$(suffix).png"
        save(joinpath(path, filename), fig)
        @info "Saved $(joinpath(path, filename))"
    end
    return nothing
end

"""
    get_modes(df::DataFrame)

Extract non-null |S| magnitude columns from a Palace port-floquet-S.csv DataFrame.
Returns `(freq, modes)` where `modes` is a `Dict{String, Vector{Float64}}`.
"""
function get_modes(df::DataFrame)
    freq = df[!, 1]
    modes = Dict{String, Vector{Float64}}()
    for name in names(df)
        m = match(r"\|S\[P(\d+)\((-?\d+);(-?\d+)\)(\w+)\]\[\d+\]\| \(dB\)", name)
        isnothing(m) && continue
        vals = df[!, name]
        all(isnan, vals) && continue
        label = "P$(m[1])($(m[2]),$(m[3]))$(m[4])"
        modes[label] = vals
    end
    return freq, modes
end
