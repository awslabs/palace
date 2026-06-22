# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# Transmon Driven Solver — Julia/Makie plot generator

Generates the transmon plots for the adaptive driven solver guide using CairoMakie.
Differs from the CPW variant in:
  * Two-port full S-matrix (PORT_PAIRS = [(1,1),(2,1),(1,2),(2,2)])
  * Two excitations indexed in domain energy (E_elec[1], E_elec[2])
  * 1×2 layout (energy plots) and 2×2 (S-mat plots)
  * Eigenmode reference lines drawn as faint vertical dotted markers
  * Frequency window 3.3–6.7 GHz, ticks at 3.5..6.5 step 0.5

Run from repo root:

```bash
julia --project=examples examples/transmon/transmon_tutorial_driven_plots.jl
```
=#

using CSV
using DataFrames
using CairoMakie
using ColorSchemes
using Printf

# -------------------------- Configuration ------------------------------------

const TRANSMON_DIR = joinpath(@__DIR__)
const REPO_ROOT = abspath(joinpath(@__DIR__, "..", ".."))
default_outdir() = joinpath(REPO_ROOT, "docs/src/assets/examples")

const UNIFORM_DIR =
    joinpath(TRANSMON_DIR, "postpro/transmon_tutorial_driven_rom/driven_uniform_reference")
const ADAPTIVE_TOLS = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
adaptive_dir(tol) = joinpath(
    TRANSMON_DIR,
    "postpro/transmon_tutorial_driven_rom",
    "driven_adaptive_1e$(round(Int, log10(tol)))"
)

const PORT_PAIRS = [(1, 1), (2, 1), (1, 2), (2, 2)]
const FREQ_LIM = (3.3, 6.7)
const FREQ_TICKS = collect(3.5:0.5:6.5)
const FREQ_MINOR_BETWEEN = 5

# Reference eigenmode frequencies (test/examples/ref/transmon/transmon_coarse).
const EIGENMODE_FREQS_GHZ = [4.099115457610e0, 5.603265962190e0]

function plasma_palette(n::Int)
    cmap = ColorSchemes.plasma
    return [cmap[t] for t in range(0.1, 0.9; length=n)]
end
const TOL_PALETTE = plasma_palette(length(ADAPTIVE_TOLS))

# -------------------------- Theme --------------------------------------------

function makie_theme()
    return Theme(
        fontsize=22,
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
            xticksize=7,
            yticksize=7,
            xtickalign=1,
            ytickalign=1,
            xminorticksvisible=true,
            yminorticksvisible=true,
            xminortickalign=1,
            yminortickalign=1,
            xminorticksize=4,
            yminorticksize=4,
            spinewidth=1.4,
            xtickwidth=1.4,
            ytickwidth=1.4,
            xticklabelsize=20,
            yticklabelsize=20,
            xlabelsize=22,
            ylabelsize=22,
            titlesize=24
        ),
        Legend=(framevisible=false, labelsize=18, patchsize=(20, 20), rowgap=4)
    )
end

# -------------------------- Data loading -------------------------------------

# domain-E columns are E_(quantity)[(excitation)] (J)  e.g. "E_elec[1] (J)".
const _RE_S_DB = r"\|S\[(\d+)\]\[(\d+)\]\| \(dB\)"
const _RE_E_COL = r"E_(\w+)\[(\d+)\] \(J\)"

function load_port_s(postpro_dir::AbstractString)
    df = CSV.read(
        joinpath(postpro_dir, "port-S.csv"),
        DataFrame;
        normalizenames=false,
        stripwhitespace=true
    )
    freq = collect(df[!, "f (GHz)"])
    s = Dict{Tuple{Int, Int}, Vector{ComplexF64}}()
    for col in names(df)
        m = match(_RE_S_DB, col)
        m === nothing && continue
        i, j = parse(Int, m.captures[1]), parse(Int, m.captures[2])
        # Some passive/reactive ports have inf magnitudes — try to coerce, skip on failure.
        db_raw = df[!, "|S[$i][$j]| (dB)"]
        ang_raw = df[!, "arg(S[$i][$j]) (deg.)"]
        try
            db = Float64.(db_raw)
            ang = Float64.(ang_raw)
            (any(isnan, db) || any(isnan, ang) || any(isinf, db) || any(isinf, ang)) &&
                continue
            s[(i, j)] = @. 10^(db / 20) * cis(deg2rad(ang))
        catch
            continue   # parsing failure (e.g. 'inf' string) → skip pair
        end
    end
    return freq, s
end

function load_domain_e(postpro_dir::AbstractString)
    df = CSV.read(
        joinpath(postpro_dir, "domain-E.csv"),
        DataFrame;
        normalizenames=false,
        stripwhitespace=true
    )
    freq = collect(df[!, "f (GHz)"])
    e = Dict{Tuple{String, Int}, Vector{Float64}}()
    for col in names(df)
        m = match(_RE_E_COL, col)
        m === nothing && continue
        e[(m.captures[1], parse(Int, m.captures[2]))] = collect(df[!, col])
    end
    return freq, e
end

function parse_log_pivots(log_path::AbstractString)
    isfile(log_path) || return Vector{Vector{Float64}}()
    content = read(log_path, String)
    re = r"Sampled frequencies \(GHz\):\s*(.*?)\n\s*Sample error"s
    out = Vector{Vector{Float64}}()
    for m in eachmatch(re, content)
        body = replace(m.captures[1], '\n' => ' ')
        nums = [parse(Float64, strip(t)) for t in split(body, ',') if !isempty(strip(t))]
        push!(out, nums)
    end
    return out
end

# -------------------------- Helpers ------------------------------------------

fmt_tol(t) = (e=round(Int, log10(t)); @sprintf("1e-%02d", -e))
mean_sq(v) = sum(abs2, v) / length(v)
collect_pivots(d, tol) = haskey(d, tol) ? reduce(vcat, d[tol]; init=Float64[]) : Float64[]

function add_eigenmode_lines!(ax)
    return vlines!(
        ax,
        EIGENMODE_FREQS_GHZ;
        color=(:black, 0.25),
        linewidth=0.6,
        linestyle=:dot
    )
end

const REF_LW = 2.0
const REF_MS = 8
const ADAPT_LW = 0.8
const ADAPT_MS = 8
const PIVOT_MS = 14
const TOL_LINE_LW = 1.8
const REF_FLOOR_LW = 1.8
const SUPTITLE_SIZE = 26
const SUPTITLE_PAD = (0, 0, 10, 10)
const SLABEL_SIZE = 26
const LEGEND_SIZE = 18
const EXC_LABEL_SIZE = 20

# -------------------------- Plot 1: S-mat uniform ----------------------------

function plot_smat_uniform(freq, sdict, outdir)
    fig = Figure(size=(900, 720))
    Label(
        fig[0, 1:2],
        L"\text{S-Parameters Magnitude } |S_{ij}| \text{ (dB)}";
        fontsize=SUPTITLE_SIZE,
        padding=SUPTITLE_PAD
    )
    axes = Matrix{Axis}(undef, 2, 2)
    for (idx, (i, j)) in enumerate(PORT_PAIRS)
        row, col = (idx - 1) ÷ 2 + 1, (idx - 1) % 2 + 1
        ax = Axis(
            fig[row, col];
            xlabel=row == 2 ? "f (GHz)" : "",
            xticks=FREQ_TICKS,
            xminorticks=IntervalsBetween(5),
            xticklabelsvisible=row == 2,
            yticklabelsvisible=col == 1
        )
        xlims!(ax, FREQ_LIM...)
        ylims!(ax, -25, 1)
        axes[row, col] = ax
        if haskey(sdict, (i, j))
            y = 20 .* log10.(abs.(sdict[(i, j)]))
            scatterlines!(
                ax,
                freq,
                y;
                color=:black,
                linewidth=REF_LW,
                markersize=REF_MS,
                marker=:circle,
                label="Uniform"
            )
        end
        add_eigenmode_lines!(ax)
        text!(
            ax,
            0.97,
            0.93;
            text=L"S_{%$i%$j}",
            space=:relative,
            align=(:right, :top),
            fontsize=SLABEL_SIZE
        )
    end
    linkaxes!(axes...)
    axislegend(axes[1, 2]; position=:rb, framevisible=false, labelsize=LEGEND_SIZE)
    save(joinpath(outdir, "driven_ua_transmon_sparam_uniform.svg"), fig)
    return fig
end

# -------------------------- Plot 2: S-mat adaptive RMS -----------------------

function plot_smat_adaptive_rms(freq, s_unif, s_adapt, freq_adapt, sampled_freqs, outdir)
    fig = Figure(size=(900, 720))
    Label(
        fig[0, 1:2],
        "Normalized Absolute Error in S-Parameters";
        fontsize=SUPTITLE_SIZE,
        padding=(0, 0, 0, 10)
    )
    Label(
        fig[1, 1:2],
        L"|S_{\mathrm{adaptive}} - S_{\mathrm{uniform}}| / \Vert S_{\mathrm{uniform}} \Vert_{\mathrm{RMS}}";
        fontsize=SUPTITLE_SIZE,
        padding=(0, 0, 6, 0)
    )
    rowgap!(fig.layout, 1, 0)
    axes = Matrix{Axis}(undef, 2, 2)
    for (idx, (i, j)) in enumerate(PORT_PAIRS)
        row, col = (idx - 1) ÷ 2 + 1, (idx - 1) % 2 + 1
        ax = Axis(
            fig[row + 1, col];
            xlabel=row == 2 ? "f (GHz)" : "",
            xticks=FREQ_TICKS,
            xminorticks=IntervalsBetween(5),
            xticklabelsvisible=row == 2,
            yticklabelsvisible=col == 1,
            yscale=log10
        )
        xlims!(ax, FREQ_LIM...)
        ylims!(ax, 5e-14, 0.9)
        axes[row, col] = ax
        haskey(s_unif, (i, j)) || continue
        ref = s_unif[(i, j)]
        scale = sqrt(mean_sq(abs.(ref)))
        tol_line_xmax = FREQ_LIM[1] + 0.8 * (FREQ_LIM[2] - FREQ_LIM[1])
        for (k, tol) in enumerate(ADAPTIVE_TOLS)
            haskey(s_adapt, tol) && haskey(s_adapt[tol], (i, j)) || continue
            err = abs.(s_adapt[tol][(i, j)] .- ref) ./ scale
            scatterlines!(
                ax,
                freq_adapt[tol],
                err;
                color=TOL_PALETTE[k],
                linewidth=ADAPT_LW,
                markersize=ADAPT_MS,
                marker=:circle,
                label=fmt_tol(tol)
            )
            # Tolerance band drawn from xmin..0.8*range, matching Python.
            lines!(
                ax,
                [FREQ_LIM[1], tol_line_xmax],
                [tol, tol];
                color=TOL_PALETTE[k],
                linewidth=TOL_LINE_LW,
                linestyle=:dash
            )
            pivots = collect_pivots(sampled_freqs, tol)
            if !isempty(pivots)
                ys = fill(1.3e-13 * (1.7^k), length(pivots))
                scatter!(
                    ax,
                    pivots,
                    ys;
                    color=(TOL_PALETTE[k], 0.75),
                    marker=:diamond,
                    markersize=PIVOT_MS
                )
            end
        end
        hlines!(ax, [1e-12]; color=:black, linewidth=REF_FLOOR_LW, linestyle=:dot)
        add_eigenmode_lines!(ax)
        text!(
            ax,
            0.97,
            0.55;
            text=L"S_{%$i%$j}",
            space=:relative,
            align=(:right, :top),
            fontsize=SLABEL_SIZE
        )
    end
    linkaxes!(axes...)
    axislegend(
        axes[2, 2];
        position=:rt,
        framevisible=false,
        labelsize=LEGEND_SIZE,
        unique=true,
        padding=(2, 2, 2, 2),
        rowgap=1
    )
    save(joinpath(outdir, "driven_ua_transmon_sparam_adaptive_rms.svg"), fig)
    return fig
end

# -------------------------- Plot 3: energy uniform ---------------------------

function plot_energy_uniform(freq, e, outdir)
    fig = Figure(size=(900, 600))
    Label(
        fig[0, 1:2],
        L"\text{Domain Energy } E_{\mathrm{elec}} \text{ (nJ)}";
        fontsize=SUPTITLE_SIZE,
        padding=SUPTITLE_PAD
    )
    axes = Vector{Axis}(undef, 2)
    for (idx, exc) in enumerate([1, 2])
        ax = Axis(
            fig[1, idx];
            xlabel="f (GHz)",
            xticks=FREQ_TICKS,
            xminorticks=IntervalsBetween(5),
            yticklabelsvisible=idx == 1
        )
        xlims!(ax, FREQ_LIM...)
        axes[idx] = ax
        key = ("elec", exc)
        if haskey(e, key)
            scatterlines!(
                ax,
                freq,
                e[key] .* 1e9;
                color=:black,
                linewidth=REF_LW,
                markersize=REF_MS,
                marker=:circle,
                label="Uniform"
            )
        end
        add_eigenmode_lines!(ax)
        text!(
            ax,
            0.97,
            0.97;
            text="Excitation $exc",
            space=:relative,
            align=(:right, :top),
            fontsize=EXC_LABEL_SIZE
        )
    end
    linkaxes!(axes...)
    axislegend(axes[2]; position=:lt, framevisible=false, labelsize=LEGEND_SIZE)
    save(joinpath(outdir, "driven_ua_transmon_energy_uniform.svg"), fig)
    return fig
end

# -------------------------- Plot 4: energy adaptive sweep --------------------

function plot_energy_adaptive_sweep(
    freq,
    e_unif,
    e_adapt,
    freq_adapt,
    sampled_freqs,
    outdir
)
    fig = Figure(size=(900, 720))
    Label(
        fig[0, 1:2],
        "Error in Domain Energy";
        fontsize=SUPTITLE_SIZE,
        padding=(0, 0, 0, 10)
    )
    Label(
        fig[1, 1:2],
        L"|E_{\mathrm{elec,adaptive}} - E_{\mathrm{elec,uniform}}| / \Vert E_{\mathrm{elec,uniform}} \Vert_{\mathrm{RMS}}";
        fontsize=SUPTITLE_SIZE,
        padding=(0, 0, 6, 0)
    )
    rowgap!(fig.layout, 1, 0)
    axes = Vector{Axis}(undef, 2)
    for (idx, exc) in enumerate([1, 2])
        ax = Axis(
            fig[2, idx];
            xlabel="f (GHz)",
            xticks=FREQ_TICKS,
            xminorticks=IntervalsBetween(5),
            yticklabelsvisible=idx == 1,
            yscale=log10
        )
        xlims!(ax, FREQ_LIM...)
        ylims!(ax, 5e-15, 10)
        axes[idx] = ax
        ref_key = ("elec", exc)
        haskey(e_unif, ref_key) || continue
        ref = e_unif[ref_key]
        scale = sqrt(mean_sq(ref))
        tol_line_xmax = FREQ_LIM[1] + 0.8 * (FREQ_LIM[2] - FREQ_LIM[1])
        for (k, tol) in enumerate(ADAPTIVE_TOLS)
            haskey(e_adapt, tol) && haskey(e_adapt[tol], ref_key) || continue
            err = abs.(e_adapt[tol][ref_key] .- ref) ./ scale
            scatterlines!(
                ax,
                freq_adapt[tol],
                err;
                color=TOL_PALETTE[k],
                linewidth=ADAPT_LW,
                markersize=ADAPT_MS,
                marker=:circle,
                label=fmt_tol(tol)
            )
            lines!(
                ax,
                [FREQ_LIM[1], tol_line_xmax],
                [tol, tol];
                color=TOL_PALETTE[k],
                linewidth=TOL_LINE_LW,
                linestyle=:dash
            )
            pivots = collect_pivots(sampled_freqs, tol)
            if !isempty(pivots)
                ys = fill(2e-13 * tol^0.25, length(pivots))
                scatter!(
                    ax,
                    pivots,
                    ys;
                    color=(TOL_PALETTE[k], 0.75),
                    marker=:diamond,
                    markersize=PIVOT_MS
                )
            end
        end
        hlines!(ax, [1e-12]; color=:black, linewidth=REF_FLOOR_LW, linestyle=:dot)
        add_eigenmode_lines!(ax)
        text!(
            ax,
            0.97,
            0.97;
            text="Excitation $exc",
            space=:relative,
            align=(:right, :top),
            fontsize=EXC_LABEL_SIZE
        )
    end
    linkaxes!(axes...)
    axislegend(
        axes[1];
        position=:rt,
        framevisible=false,
        labelsize=LEGEND_SIZE,
        unique=true,
        padding=(2, 2, 2, 2),
        rowgap=1
    )
    save(joinpath(outdir, "driven_ua_transmon_energy_adaptive_sweep.svg"), fig)
    return fig
end

# -------------------------- Driver -------------------------------------------

function main(; outdir=default_outdir())
    set_theme!(makie_theme())
    mkpath(outdir)

    freq_u, s_u = load_port_s(UNIFORM_DIR)
    _, e_u = load_domain_e(UNIFORM_DIR)
    println("Uniform: $(length(freq_u)) frequency points")

    freq_a = Dict{Float64, Vector{Float64}}()
    s_a = Dict{Float64, Dict{Tuple{Int, Int}, Vector{ComplexF64}}}()
    e_a = Dict{Float64, Dict{Tuple{String, Int}, Vector{Float64}}}()
    pivots = Dict{Float64, Vector{Vector{Float64}}}()
    for tol in ADAPTIVE_TOLS
        d = adaptive_dir(tol)
        if !isfile(joinpath(d, "port-S.csv"))
            @warn "Missing data for tol=$tol at $d"
            continue
        end
        freq_a[tol], s_a[tol] = load_port_s(d)
        _, e_a[tol] = load_domain_e(d)
        pivots[tol] = parse_log_pivots(joinpath(d, "palace.log"))
        n = sum(length, pivots[tol]; init=0)
        println("Adaptive tol=$tol: $(length(freq_a[tol])) freq pts ($n pivots)")
    end

    plot_smat_uniform(freq_u, s_u, outdir)
    plot_smat_adaptive_rms(freq_u, s_u, s_a, freq_a, pivots, outdir)
    plot_energy_uniform(freq_u, e_u, outdir)
    plot_energy_adaptive_sweep(freq_u, e_u, e_a, freq_a, pivots, outdir)

    println("Wrote 4 transmon SVGs to $outdir")
    return outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
    let out = default_outdir()
        args = copy(ARGS)
        while !isempty(args)
            a = popfirst!(args)
            if a == "--out"
                out = popfirst!(args)
            end
        end
        main(; outdir=out)
    end
end
