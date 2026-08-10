# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# Transmon Circuit-Synthesis Tutorial — Julia/Makie plot generator

Regenerate the two circuit-synthesis comparison figures from output produced by
`transmon_tutorial_circuit.jl`:

  * `circuit_transmon_domain_e_comparison.svg`
  * `circuit_transmon_port_s_comparison.svg`

From the repository root:

```bash
julia --project=examples examples/transmon/transmon_tutorial_circuit_plots.jl
```

Use `--data-root <dir>` to read another output tree and `--out <dir>` to write the
figures somewhere other than `docs/src/assets/examples`.
=#

using CSV
using CairoMakie
using ColorSchemes
using DataFrames

const TRANSMON_DIR = @__DIR__
const REPO_ROOT = abspath(joinpath(TRANSMON_DIR, "..", ".."))
const CIRCUIT_TOLS = [1e-1, 1e-3, 1e-5]
const PORT_PAIRS = [(1, 1), (2, 1), (1, 2), (2, 2)]
const EXCITATIONS = [1, 2]
const FREQ_LIM = (3.3, 6.7)
const FREQ_TICKS = collect(3.5:1.0:6.5)
const EIGENMODE_FREQS_GHZ = [4.099115457610, 5.603265962190]
const TOL_COLORS =
    [ColorSchemes.plasma[t] for t in range(0.1, 0.9; length=length(CIRCUIT_TOLS))]

const _RE_S_DB = r"\|S\[(\d+)\]\[(\d+)\]\| \(dB\)"
const _RE_E_COL = r"E_(\w+)\[(\d+)\] \(J\)"

default_data_root() = joinpath(TRANSMON_DIR, "postpro", "transmon_tutorial_circuit")
default_outdir() = joinpath(REPO_ROOT, "docs", "src", "assets", "examples")
tol_label(tol) = "1e$(round(Int, log10(tol)))"
mean_sq(v) = sum(abs2, v) / length(v)

function load_port_s(postpro_dir::AbstractString)
    path = joinpath(postpro_dir, "port-S.csv")
    isfile(path) || error("Missing Palace output: $path")
    df = CSV.read(path, DataFrame; normalizenames=false, stripwhitespace=true)
    freq = Float64.(df[!, "f (GHz)"])
    s = Dict{Tuple{Int, Int}, Vector{ComplexF64}}()
    for col in names(df)
        m = match(_RE_S_DB, col)
        m === nothing && continue
        i, j = parse(Int, m.captures[1]), parse(Int, m.captures[2])
        # Passive/reactive ports can produce string-valued `-inf` magnitudes. Those
        # columns are not measurable S-parameter pairs and should not prevent loading
        # the finite port pairs used by the tutorial.
        try
            db = Float64.(df[!, "|S[$i][$j]| (dB)"])
            angle = Float64.(df[!, "arg(S[$i][$j]) (deg.)"])
            all(isfinite, db) && all(isfinite, angle) || continue
            s[(i, j)] = @. 10^(db / 20) * cis(deg2rad(angle))
        catch error
            error isa ArgumentError || error isa MethodError || rethrow()
            continue
        end
    end
    return freq, s
end

function load_domain_e(postpro_dir::AbstractString)
    path = joinpath(postpro_dir, "domain-E.csv")
    isfile(path) || error("Missing Palace output: $path")
    df = CSV.read(path, DataFrame; normalizenames=false, stripwhitespace=true)
    freq = Float64.(df[!, "f (GHz)"])
    energy = Dict{Tuple{String, Int}, Vector{Float64}}()
    for col in names(df)
        m = match(_RE_E_COL, col)
        m === nothing && continue
        try
            energy[(m.captures[1], parse(Int, m.captures[2]))] = Float64.(df[!, col])
        catch error
            error isa ArgumentError || error isa MethodError || rethrow()
            continue
        end
    end
    return freq, energy
end

function check_frequency_grid(reference, candidate, label)
    same =
        length(reference) == length(candidate) &&
        all(isapprox.(reference, candidate; rtol=0.0, atol=1e-12))
    same || error("Frequency grid for $label does not match the uniform reference")
    return
end

function validate_required_series(label, freq, s_parameters, energy; require_nonzero=false)
    problems = String[]
    isempty(freq) && push!(problems, "frequency grid is empty")
    all(isfinite, freq) || push!(problems, "frequency grid contains non-finite values")

    for pair in PORT_PAIRS
        pair_name = "S[$(pair[1])][$(pair[2])]"
        if !haskey(s_parameters, pair)
            push!(problems, "missing or malformed $pair_name")
            continue
        end
        values = s_parameters[pair]
        length(values) == length(freq) ||
            push!(problems, "$pair_name has the wrong number of rows")
        all(isfinite, values) || push!(problems, "$pair_name contains non-finite values")
        require_nonzero &&
            iszero(mean_sq(values)) &&
            push!(problems, "$pair_name has zero RMS normalization")
    end

    for excitation in EXCITATIONS
        key = ("elec", excitation)
        quantity_name = "E_elec[$excitation]"
        if !haskey(energy, key)
            push!(problems, "missing or malformed $quantity_name")
            continue
        end
        values = energy[key]
        length(values) == length(freq) ||
            push!(problems, "$quantity_name has the wrong number of rows")
        all(isfinite, values) ||
            push!(problems, "$quantity_name contains non-finite values")
        require_nonzero &&
            iszero(mean_sq(values)) &&
            push!(problems, "$quantity_name has zero RMS normalization")
    end

    isempty(problems) ||
        error("Invalid circuit-tutorial data in $label:\n  - " * join(problems, "\n  - "))
    return
end

function load_tutorial_data(data_root)
    uniform_dir = joinpath(data_root, "driven_uniform_reference")
    freq, uniform_s = load_port_s(uniform_dir)
    freq_e, uniform_e = load_domain_e(uniform_dir)
    check_frequency_grid(freq, freq_e, "uniform domain energy")
    validate_required_series(uniform_dir, freq, uniform_s, uniform_e; require_nonzero=true)

    adaptive_s = Dict{Float64, Dict{Tuple{Int, Int}, Vector{ComplexF64}}}()
    circuit_s = Dict{Float64, Dict{Tuple{Int, Int}, Vector{ComplexF64}}}()
    adaptive_e = Dict{Float64, Dict{Tuple{String, Int}, Vector{Float64}}}()
    circuit_e = Dict{Float64, Dict{Tuple{String, Int}, Vector{Float64}}}()

    for tol in CIRCUIT_TOLS
        label = tol_label(tol)
        adaptive_dir = joinpath(data_root, "driven_adaptive_$label")
        circuit_dir = joinpath(data_root, "driven_circuit_$label")

        freq_a, adaptive_s[tol] = load_port_s(adaptive_dir)
        freq_ae, adaptive_e[tol] = load_domain_e(adaptive_dir)
        freq_c, circuit_s[tol] = load_port_s(circuit_dir)
        freq_ce, circuit_e[tol] = load_domain_e(circuit_dir)
        for (grid, name) in (
            (freq_a, "adaptive S at $label"),
            (freq_ae, "adaptive energy at $label"),
            (freq_c, "circuit S at $label"),
            (freq_ce, "circuit energy at $label")
        )
            check_frequency_grid(freq, grid, name)
        end
        validate_required_series(adaptive_dir, freq_a, adaptive_s[tol], adaptive_e[tol])
        validate_required_series(circuit_dir, freq_c, circuit_s[tol], circuit_e[tol])
    end
    return (; freq, uniform_s, uniform_e, adaptive_s, adaptive_e, circuit_s, circuit_e)
end

function tutorial_plot_theme()
    return Theme(
        fontsize=18,
        fonts=(
            regular="Times",
            italic="Times Italic",
            bold="Times Bold",
            bold_italic="Times Bold Italic"
        ),
        Axis=(
            xgridvisible=false,
            ygridvisible=false,
            rightspinevisible=false,
            topspinevisible=false,
            spinewidth=0.75,
            spinecolor=(:black, 0.55),
            xtickwidth=0.75,
            ytickwidth=0.75,
            xtickcolor=(:black, 0.55),
            ytickcolor=(:black, 0.55),
            xticksize=3,
            yticksize=3,
            xtickalign=0,
            ytickalign=0,
            xminorticksvisible=false,
            yminorticksvisible=false,
            xticklabelsize=19,
            yticklabelsize=19,
            xlabelsize=20,
            xlabelcolor=(:black, 0.75)
        )
    )
end

function comparison_axis(fig, row, col; bottom, logscale=false, yticks=nothing)
    options = (
        xlabel=bottom ? "Frequency (GHz)" : "",
        xticks=FREQ_TICKS,
        xticklabelsvisible=bottom,
        xticksvisible=bottom,
        bottomspinevisible=bottom,
        yticklabelsvisible=col != 3,
        backgroundcolor=:white
    )
    ax = if logscale
        Axis(fig[row, col]; options..., yscale=log10, yticks)
    else
        Axis(fig[row, col]; options...)
    end
    xlims!(ax, FREQ_LIM...)
    vlines!(ax, EIGENMODE_FREQS_GHZ; color=(:black, 0.16), linewidth=0.65, linestyle=:dot)
    return ax
end

function add_headers!(fig, reference_title)
    style = (fontsize=22, font="Times Bold", color=(:black, 0.82), padding=(0, 0, 3, 3))
    Label(fig[1, 1], reference_title; style...)
    Label(fig[1, 2], "Error vs uniform"; style...)
    Label(fig[1, 3], "Error vs adaptive"; style...)
    return
end

function add_footer!(fig, row)
    Label(
        fig[row, 1:3],
        "RMS-normalized errors  ·  values below the dotted floor are clipped";
        fontsize=17,
        color=(:black, 0.52),
        padding=(0, 0, 2, 0)
    )
    return
end

function annotate_eigenfrequencies!(ax, y)
    for (index, frequency) in enumerate(EIGENMODE_FREQS_GHZ)
        text!(
            ax,
            frequency - 0.08,
            y;
            text=index == 1 ? "f₁" : "f₂",
            align=(:right, :top),
            color=(:black, 0.48),
            fontsize=16
        )
    end
    return
end

function add_error_guides!(uniform_axis, adaptive_axis, row, floor)
    labels = ("10⁻¹", "10⁻³", "10⁻⁵")
    for (tol, color, label) in zip(CIRCUIT_TOLS, TOL_COLORS, labels)
        for ax in (uniform_axis, adaptive_axis)
            hlines!(ax, [tol]; color=(color, 0.58), linewidth=0.8, linestyle=:dash)
            if row == 1
                text!(ax, 6.62, tol; text=label, align=(:right, :top), color, fontsize=18)
            end
        end
    end
    for ax in (uniform_axis, adaptive_axis)
        hlines!(ax, [floor]; color=(:black, 0.36), linewidth=0.7, linestyle=:dot)
    end
    return
end

function plot_errors!(axes, freq, circuit, adaptive, reference, scale, color, floor)
    error_vs_uniform = max.(abs.(circuit .- reference) ./ scale, floor)
    error_vs_adaptive = max.(abs.(circuit .- adaptive) ./ scale, floor)
    lines!(axes[1], freq, error_vs_uniform; color, linewidth=1.45)
    lines!(axes[2], freq, error_vs_adaptive; color, linewidth=1.45)
    return
end

function plot_domain_energy(data, outdir)
    fig = Figure(size=(960, 560), figure_padding=(14, 14, 10, 10), backgroundcolor=:white)
    add_headers!(fig, "Electric energy (nJ)")

    reference_axes = Axis[]
    error_axes = Axis[]
    error_ticks =
        ([1e-11, 1e-8, 1e-5, 1e-2], [L"10^{-11}", L"10^{-8}", L"10^{-5}", L"10^{-2}"])
    floor = 1e-12

    for (row_index, excitation) in enumerate(EXCITATIONS)
        row = row_index + 1
        bottom = row_index == length(EXCITATIONS)
        reference_axis = comparison_axis(fig, row, 1; bottom)
        uniform_axis =
            comparison_axis(fig, row, 2; bottom, logscale=true, yticks=error_ticks)
        adaptive_axis =
            comparison_axis(fig, row, 3; bottom, logscale=true, yticks=error_ticks)
        ylims!(reference_axis, -0.05, 1.6)
        ylims!(uniform_axis, 2e-13, 3e-1)
        ylims!(adaptive_axis, 2e-13, 3e-1)
        push!(reference_axes, reference_axis)
        push!(error_axes, uniform_axis, adaptive_axis)

        key = ("elec", excitation)
        reference = data.uniform_e[key]
        scale = sqrt(mean_sq(reference))
        lines!(
            reference_axis,
            data.freq,
            reference .* 1e9;
            color=colorant"#222222",
            linewidth=1.8
        )
        text!(
            reference_axis,
            0.97,
            0.68;
            text="Exc. $excitation",
            space=:relative,
            align=(:right, :top),
            color=(:black, 0.72),
            fontsize=19
        )
        row_index == 1 && annotate_eigenfrequencies!(reference_axis, 1.55)

        for (tol, color) in zip(CIRCUIT_TOLS, TOL_COLORS)
            plot_errors!(
                (uniform_axis, adaptive_axis),
                data.freq,
                data.circuit_e[tol][key],
                data.adaptive_e[tol][key],
                reference,
                scale,
                color,
                floor
            )
        end
        add_error_guides!(uniform_axis, adaptive_axis, row_index, floor)
    end

    linkxaxes!(reference_axes..., error_axes...)
    linkyaxes!(reference_axes...)
    linkyaxes!(error_axes...)
    rowgap!(fig.layout, 10)
    colgap!(fig.layout, 20)
    add_footer!(fig, 4)

    path = joinpath(outdir, "circuit_transmon_domain_e_comparison.svg")
    save(path, fig)
    return path
end

function plot_s_parameters(data, outdir)
    fig = Figure(size=(960, 1000), figure_padding=(14, 14, 10, 10), backgroundcolor=:white)
    add_headers!(fig, "Scattering magnitude (dB)")

    reference_axes = Axis[]
    error_axes = Axis[]
    error_ticks =
        ([1e-10, 1e-7, 1e-4, 1e-1], [L"10^{-10}", L"10^{-7}", L"10^{-4}", L"10^{-1}"])
    floor = 1e-11

    for (row_index, pair) in enumerate(PORT_PAIRS)
        i, j = pair
        row = row_index + 1
        bottom = row_index == length(PORT_PAIRS)
        reference_axis = comparison_axis(fig, row, 1; bottom)
        uniform_axis =
            comparison_axis(fig, row, 2; bottom, logscale=true, yticks=error_ticks)
        adaptive_axis =
            comparison_axis(fig, row, 3; bottom, logscale=true, yticks=error_ticks)
        ylims!(reference_axis, -22, 1)
        ylims!(uniform_axis, 5e-13, 4e-1)
        ylims!(adaptive_axis, 5e-13, 4e-1)
        push!(reference_axes, reference_axis)
        push!(error_axes, uniform_axis, adaptive_axis)

        reference = data.uniform_s[pair]
        scale = sqrt(mean_sq(reference))
        lines!(
            reference_axis,
            data.freq,
            20 .* log10.(abs.(reference));
            color=colorant"#222222",
            linewidth=1.8
        )
        for ax in (reference_axis, uniform_axis, adaptive_axis)
            text!(
                ax,
                0.03,
                0.93;
                text=L"S_{%$i%$j}",
                space=:relative,
                align=(:left, :top),
                color=(:black, 0.68),
                fontsize=19
            )
        end
        row_index == 1 && annotate_eigenfrequencies!(reference_axis, 0.2)

        for (tol, color) in zip(CIRCUIT_TOLS, TOL_COLORS)
            plot_errors!(
                (uniform_axis, adaptive_axis),
                data.freq,
                data.circuit_s[tol][pair],
                data.adaptive_s[tol][pair],
                reference,
                scale,
                color,
                floor
            )
        end
        add_error_guides!(uniform_axis, adaptive_axis, row_index, floor)
    end

    linkxaxes!(reference_axes..., error_axes...)
    linkyaxes!(reference_axes...)
    linkyaxes!(error_axes...)
    rowgap!(fig.layout, 10)
    colgap!(fig.layout, 20)
    add_footer!(fig, 6)

    path = joinpath(outdir, "circuit_transmon_port_s_comparison.svg")
    save(path, fig)
    return path
end

function main(; data_root=default_data_root(), outdir=default_outdir())
    set_theme!(tutorial_plot_theme())
    mkpath(outdir)
    data = load_tutorial_data(data_root)
    energy_path = plot_domain_energy(data, outdir)
    s_path = plot_s_parameters(data, outdir)
    println("Wrote circuit-synthesis tutorial figures:")
    println("  $energy_path")
    println("  $s_path")
    return (energy_path, s_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    let data_root = default_data_root(), outdir = default_outdir()
        args = copy(ARGS)
        while !isempty(args)
            arg = popfirst!(args)
            if arg == "--data-root"
                isempty(args) && error("--data-root requires a path")
                data_root = popfirst!(args)
            elseif arg == "--out"
                isempty(args) && error("--out requires a path")
                outdir = popfirst!(args)
            else
                error("Unknown argument: $arg")
            end
        end
        main(; data_root, outdir)
    end
end
