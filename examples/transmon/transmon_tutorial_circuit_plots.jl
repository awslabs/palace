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
using Printf

const TRANSMON_DIR = @__DIR__
const REPO_ROOT = abspath(joinpath(TRANSMON_DIR, "..", ".."))
const CIRCUIT_TOLS = [1e-1, 1e-3, 1e-5]
const PORT_PAIRS = [(1, 1), (2, 1), (1, 2), (2, 2)]
const EXCITATIONS = [1, 2]
const FREQ_LIM = (3.3, 6.7)
const FREQ_TICKS = collect(3.5:0.5:6.5)
const EIGENMODE_FREQS_GHZ = [4.099115457610, 5.603265962190]
const TOL_COLORS =
    [ColorSchemes.plasma[t] for t in range(0.1, 0.9; length=length(CIRCUIT_TOLS))]

const _RE_S_DB = r"\|S\[(\d+)\]\[(\d+)\]\| \(dB\)"
const _RE_E_COL = r"E_(\w+)\[(\d+)\] \(J\)"

default_data_root() = joinpath(TRANSMON_DIR, "postpro", "transmon_tutorial_circuit")
default_outdir() = joinpath(REPO_ROOT, "docs", "src", "assets", "examples")
tol_label(tol) = "1e$(round(Int, log10(tol)))"
legend_tol(tol) = @sprintf("1e-%02d", -round(Int, log10(tol)))
mean_sq(v) = sum(abs2, v) / length(v)

function circuit_theme()
    return Theme(
        fontsize=14,
        fonts=(
            regular="Times",
            italic="Times Italic",
            bold="Times Bold",
            bold_italic="Times Bold Italic"
        ),
        Axis=(
            xgridvisible=true,
            xgridcolor=(:black, 0.18),
            xgridstyle=:dot,
            ygridvisible=false,
            rightspinevisible=false,
            topspinevisible=false,
            xticksize=4,
            yticksize=4,
            xtickalign=1,
            ytickalign=1,
            xminorticksvisible=true,
            yminorticksvisible=true,
            xminortickalign=1,
            yminortickalign=1,
            xminorticksize=2.5,
            yminorticksize=2.5,
            spinewidth=1.0,
            xtickwidth=1.0,
            ytickwidth=1.0,
            xticklabelsize=12,
            yticklabelsize=12,
            xlabelsize=13,
            titlesize=14
        ),
        Legend=(framevisible=false, labelsize=12, patchsize=(14, 12), rowgap=2)
    )
end

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

function add_frequency_guides!(ax)
    xlims!(ax, FREQ_LIM...)
    return vlines!(
        ax,
        EIGENMODE_FREQS_GHZ;
        color=(:black, 0.25),
        linewidth=0.6,
        linestyle=:dot
    )
end

function comparison_axis(fig, row, col; bottom, logscale=false, yticks=nothing)
    kwargs = (
        xlabel=bottom ? "f (GHz)" : "",
        xticks=FREQ_TICKS,
        xminorticks=IntervalsBetween(5),
        xticklabelsvisible=bottom,
        yticklabelsvisible=col != 3
    )
    ax = if logscale
        Axis(fig[row, col]; kwargs..., yscale=log10, yticks=yticks)
    else
        Axis(fig[row, col]; kwargs...)
    end
    add_frequency_guides!(ax)
    return ax
end

function plot_domain_energy(data, outdir)
    fig = Figure(size=(864, 427), figure_padding=7)
    Label(
        fig[0, 1:3],
        L"\text{Domain Energy } E_{\mathrm{elec}}";
        fontsize=18,
        padding=(0, 0, 3, 2)
    )
    Label(
        fig[1, 1],
        L"\text{Uniform Driven Solver } E_{\mathrm{elec}}\;\text{(nJ)}";
        fontsize=14
    )
    Label(
        fig[1, 2],
        L"|E_{\mathrm{circuit}}-E_{\mathrm{uniform}}|/\Vert E_{\mathrm{uniform}}\Vert_{\mathrm{RMS}}";
        fontsize=13
    )
    Label(
        fig[1, 3],
        L"|E_{\mathrm{circuit}}-E_{\mathrm{adaptive}}|/\Vert E_{\mathrm{uniform}}\Vert_{\mathrm{RMS}}";
        fontsize=13
    )

    uniform_axes = Axis[]
    error_axes = Axis[]
    for (row_index, excitation) in enumerate(EXCITATIONS)
        row = row_index + 1
        bottom = row_index == length(EXCITATIONS)
        ax_uniform = comparison_axis(fig, row, 1; bottom)
        ax_uniform.yticks = 0.0:0.5:1.5
        ylims!(ax_uniform, -0.05, 1.6)
        push!(uniform_axes, ax_uniform)

        energy_ticks = [1e-11, 1e-8, 1e-5, 1e-2]
        ax_uniform_error =
            comparison_axis(fig, row, 2; bottom, logscale=true, yticks=energy_ticks)
        ax_adaptive_error =
            comparison_axis(fig, row, 3; bottom, logscale=true, yticks=energy_ticks)
        ylims!(ax_uniform_error, 2e-13, 3e-1)
        ylims!(ax_adaptive_error, 2e-13, 3e-1)
        push!(error_axes, ax_uniform_error, ax_adaptive_error)

        key = ("elec", excitation)
        reference = data.uniform_e[key]
        scale = sqrt(mean_sq(reference))
        scatterlines!(
            ax_uniform,
            data.freq,
            reference .* 1e9;
            color=:black,
            linewidth=1.2,
            markersize=2.5
        )
        text!(
            ax_uniform,
            0.97,
            0.95;
            text="Excitation $excitation",
            space=:relative,
            align=(:right, :top),
            fontsize=12
        )

        for (k, tol) in enumerate(CIRCUIT_TOLS)
            circuit = data.circuit_e[tol][key]
            adaptive = data.adaptive_e[tol][key]
            error_uniform = abs.(circuit .- reference) ./ scale
            error_adaptive = abs.(circuit .- adaptive) ./ scale
            color = TOL_COLORS[k]
            scatterlines!(
                ax_uniform_error,
                data.freq,
                error_uniform;
                color,
                linewidth=0.7,
                markersize=2.5,
                label=legend_tol(tol)
            )
            scatterlines!(
                ax_adaptive_error,
                data.freq,
                error_adaptive;
                color,
                linewidth=0.7,
                markersize=2.5,
                label=legend_tol(tol)
            )
            hlines!(ax_uniform_error, [tol]; color, linewidth=1.0, linestyle=:dash)
        end
        for ax in (ax_uniform_error, ax_adaptive_error)
            hlines!(ax, [1e-12]; color=:black, linewidth=1.0, linestyle=:dot)
        end
    end
    linkxaxes!(uniform_axes..., error_axes...)
    linkyaxes!(uniform_axes...)
    linkyaxes!(error_axes...)
    axislegend(error_axes[2]; position=:rt, unique=true, padding=(2, 2, 2, 2))
    rowgap!(fig.layout, 0)
    colgap!(fig.layout, 7)

    path = joinpath(outdir, "circuit_transmon_domain_e_comparison.svg")
    save(path, fig)
    return path
end

function plot_s_parameters(data, outdir)
    fig = Figure(size=(864, 783), figure_padding=7)
    Label(fig[0, 1:3], "S-Parameters"; fontsize=18, padding=(0, 0, 3, 2))
    Label(fig[1, 1], L"\text{Uniform Driven Solver } |S_{ij}|\;\text{(dB)}"; fontsize=14)
    Label(
        fig[1, 2],
        L"|S_{\mathrm{circuit}}-S_{\mathrm{uniform}}|/\Vert S_{\mathrm{uniform}}\Vert_{\mathrm{RMS}}";
        fontsize=13
    )
    Label(
        fig[1, 3],
        L"|S_{\mathrm{circuit}}-S_{\mathrm{adaptive}}|/\Vert S_{\mathrm{uniform}}\Vert_{\mathrm{RMS}}";
        fontsize=13
    )

    uniform_axes = Axis[]
    error_axes = Axis[]
    for (row_index, pair) in enumerate(PORT_PAIRS)
        i, j = pair
        row = row_index + 1
        bottom = row_index == length(PORT_PAIRS)
        ax_uniform = comparison_axis(fig, row, 1; bottom)
        ax_uniform.yticks = -20:5:0
        ylims!(ax_uniform, -22, 1)
        push!(uniform_axes, ax_uniform)

        error_ticks = [1e-10, 1e-7, 1e-4, 1e-1]
        ax_uniform_error =
            comparison_axis(fig, row, 2; bottom, logscale=true, yticks=error_ticks)
        ax_adaptive_error =
            comparison_axis(fig, row, 3; bottom, logscale=true, yticks=error_ticks)
        ylims!(ax_uniform_error, 3e-12, 4e-1)
        ylims!(ax_adaptive_error, 3e-12, 4e-1)
        push!(error_axes, ax_uniform_error, ax_adaptive_error)

        reference = data.uniform_s[pair]
        scale = sqrt(mean_sq(reference))
        scatterlines!(
            ax_uniform,
            data.freq,
            20 .* log10.(abs.(reference));
            color=:black,
            linewidth=1.2,
            markersize=2.5
        )
        for ax in (ax_uniform, ax_uniform_error, ax_adaptive_error)
            text!(
                ax,
                0.03,
                0.94;
                text=L"S_{%$i%$j}",
                space=:relative,
                align=(:left, :top),
                fontsize=13
            )
        end

        for (k, tol) in enumerate(CIRCUIT_TOLS)
            circuit = data.circuit_s[tol][pair]
            adaptive = data.adaptive_s[tol][pair]
            error_uniform = abs.(circuit .- reference) ./ scale
            error_adaptive = abs.(circuit .- adaptive) ./ scale
            color = TOL_COLORS[k]
            scatterlines!(
                ax_uniform_error,
                data.freq,
                error_uniform;
                color,
                linewidth=0.7,
                markersize=2.5,
                label=legend_tol(tol)
            )
            scatterlines!(
                ax_adaptive_error,
                data.freq,
                error_adaptive;
                color,
                linewidth=0.7,
                markersize=2.5,
                label=legend_tol(tol)
            )
            hlines!(ax_uniform_error, [tol]; color, linewidth=1.0, linestyle=:dash)
        end
        for ax in (ax_uniform_error, ax_adaptive_error)
            hlines!(ax, [1e-12]; color=:black, linewidth=1.0, linestyle=:dot)
        end
    end
    linkxaxes!(uniform_axes..., error_axes...)
    linkyaxes!(uniform_axes...)
    linkyaxes!(error_axes...)
    axislegend(error_axes[2]; position=:rt, unique=true, padding=(2, 2, 2, 2))
    rowgap!(fig.layout, 0)
    colgap!(fig.layout, 7)

    path = joinpath(outdir, "circuit_transmon_port_s_comparison.svg")
    save(path, fig)
    return path
end

function main(; data_root=default_data_root(), outdir=default_outdir())
    set_theme!(circuit_theme())
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
