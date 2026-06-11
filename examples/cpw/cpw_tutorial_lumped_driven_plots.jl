# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# CPW Lumped-Port Driven Solver — Julia/Makie plot generator

Generates the CPW plots for the adaptive driven solver guide using CairoMakie, with
matching colour palette (matplotlib `plasma`, sampled at 0.1..0.9), serif typography,
multi-line titles, and tolerance/pivot annotations.

Run from the repository root:

```bash
julia --project=examples examples/cpw/cpw_tutorial_lumped_driven_plots.jl
```

Output SVGs go to `docs/src/assets/examples/`; pass `--out <dir>` to override.
=#

using CSV
using DataFrames
using CairoMakie
using ColorSchemes
using Printf

# -------------------------- Configuration -------------------------------------

const CPW_DIR = joinpath(@__DIR__)
const REPO_ROOT = abspath(joinpath(@__DIR__, "..", ".."))
default_outdir() = joinpath(REPO_ROOT, "docs/src/assets/examples")

const UNIFORM_DIR =
    joinpath(CPW_DIR, "postpro/tutorial_driven_rom/driven_uniform_reference")
const ADAPTIVE_TOLS = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
adaptive_dir(tol) = joinpath(
    CPW_DIR,
    "postpro/tutorial_driven_rom",
    "driven_adaptive_1e$(round(Int, log10(tol)))"
)

# (out_idx, src_idx) port pairs in row-major order matching Python's PORT_PAIRS
const PORT_PAIRS = [(1, 1), (2, 1), (3, 1), (4, 1)]

# matplotlib plasma sampled at np.linspace(0.1, 0.9, 5)
function plasma_palette(n::Int)
    cmap = ColorSchemes.plasma
    return [cmap[t] for t in range(0.1, 0.9; length=n)]
end
const TOL_PALETTE = plasma_palette(length(ADAPTIVE_TOLS))

const GOLDEN_RATIO = (1 + sqrt(5)) / 2

# -------------------------- Theme ---------------------------------------------

# Matches matplotlib rcParams: serif (Times), 150 dpi.
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

# -------------------------- Data loading --------------------------------------

const _RE_S_DB = r"\|S\[(\d+)\]\[(\d+)\]\| \(dB\)"
const _RE_E_COL = r"E_(\w+) \(J\)"

"""
Load port-S.csv → (freq_GHz, Dict{(Int,Int), Vector{ComplexF64}}).
"""
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
        if m !== nothing
            i, j = parse(Int, m.captures[1]), parse(Int, m.captures[2])
            db = collect(df[!, "|S[$i][$j]| (dB)"])
            ang = collect(df[!, "arg(S[$i][$j]) (deg.)"])
            s[(i, j)] = @. 10^(db / 20) * cis(deg2rad(ang))
        end
    end
    return freq, s
end

"""
Load domain-E.csv → (freq_GHz, Dict{String, Vector{Float64}}).
"""
function load_domain_e(postpro_dir::AbstractString)
    df = CSV.read(
        joinpath(postpro_dir, "domain-E.csv"),
        DataFrame;
        normalizenames=false,
        stripwhitespace=true
    )
    freq = collect(df[!, "f (GHz)"])
    e = Dict{String, Vector{Float64}}()
    for col in names(df)
        m = match(_RE_E_COL, col)
        if m !== nothing
            e[m.captures[1]] = collect(df[!, col])
        end
    end
    return freq, e
end

"""
Parse the log block ``Sampled frequencies (GHz): … (until next``Sample errors``).
The block can span multiple lines with whitespace continuations. Returns a vector
of arrays, one per excitation occurrence.
"""
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

"""
Parse `Sample errors: …` blocks (multi-line, terminated by `Total offline phase`).
"""
function parse_log_errors(log_path::AbstractString)
    isfile(log_path) || return Vector{Vector{Float64}}()
    content = read(log_path, String)
    re = r"Sample errors:\s*(.*?)\n\s*Total offline phase"s
    out = Vector{Vector{Float64}}()
    for m in eachmatch(re, content)
        body = replace(m.captures[1], '\n' => ' ')
        nums = Float64[]
        for t in split(body, ',')
            s = strip(t)
            isempty(s) && continue
            push!(nums, parse(Float64, s))
        end
        push!(out, nums)
    end
    return out
end

# -------------------------- Helpers -------------------------------------------

"""
Format scientific tolerance like Python's `f\"{tol:.0e}\"` → `1e-01`.
"""
fmt_tol(t) = (e=round(Int, log10(t)); @sprintf("1e-%02d", -e))

# Style constants matching the Python script
const REF_COLOR = :black
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

# -------------------------- Plot 1: S-mat uniform -----------------------------

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
            xticks=collect(range(2, 32; step=2)),
            xminorticks=IntervalsBetween(2),
            xticklabelsvisible=row == 2,
            yticklabelsvisible=col == 1
        )
        ylims!(ax, -105, 5)
        axes[row, col] = ax
        if haskey(sdict, (i, j))
            y = 20 .* log10.(abs.(sdict[(i, j)]))
            scatterlines!(
                ax,
                freq,
                y;
                color=REF_COLOR,
                linewidth=REF_LW,
                markersize=REF_MS,
                marker=:circle,
                label="Uniform Driven Solver"
            )
        end
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
    axislegend(
        axes[1, 2];
        position=:rb,
        framevisible=false,
        labelsize=LEGEND_SIZE,
        padding=(0, 0, 0, 0)
    )
    save(joinpath(outdir, "driven_ua_cpw_domain_sparam_uniform.svg"), fig)
    return fig
end

# -------------------------- Plot 2: S-mat adaptive pointwise ------------------

function plot_smat_adaptive_pointwise(freq, s_unif, s_adapt, freq_adapt, outdir)
    fig = Figure(size=(900, 720))
    Label(
        fig[0, 1:2],
        "Magnitude of Pointwise Relative Error in S-Parameters";
        fontsize=SUPTITLE_SIZE,
        padding=(0, 0, 0, 10)
    )
    Label(
        fig[1, 1:2],
        L"|S_{\mathrm{adaptive}} - S_{\mathrm{uniform}}| / |S_{\mathrm{uniform}}|";
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
            xticks=collect(range(2, 32; step=2)),
            xminorticks=IntervalsBetween(2),
            xticklabelsvisible=row == 2,
            yticklabelsvisible=col == 1,
            yscale=log10
        )
        xlims!(ax, 1, 33)
        ylims!(ax, 1.1e-13, 900)
        axes[row, col] = ax
        haskey(s_unif, (i, j)) || continue
        ref = s_unif[(i, j)]
        for (k, tol) in enumerate(ADAPTIVE_TOLS)
            haskey(s_adapt, tol) && haskey(s_adapt[tol], (i, j)) || continue
            err = abs.(s_adapt[tol][(i, j)] .- ref) ./ abs.(ref)
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
            hlines!(ax, [tol]; color=TOL_PALETTE[k], linewidth=TOL_LINE_LW, linestyle=:dash)
        end
        hlines!(ax, [1e-12]; color=:black, linewidth=REF_FLOOR_LW, linestyle=:dot)
        text!(
            ax,
            0.97,
            0.97;
            text=L"S_{%$i%$j}",
            space=:relative,
            align=(:right, :top),
            fontsize=SLABEL_SIZE
        )
    end
    linkaxes!(axes...)
    axislegend(
        axes[2, 2];
        position=:cb,
        framevisible=false,
        labelsize=LEGEND_SIZE,
        nbanks=3,
        unique=true,
        padding=(2, 2, 2, 2),
        rowgap=2,
        colgap=10
    )
    save(joinpath(outdir, "driven_ua_cpw_domain_sparam_adaptive_pointwise.svg"), fig)
    return fig
end

# -------------------------- Plot 3: S-mat adaptive RMS -----------------------

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
            xticks=collect(range(2, 32; step=2)),
            xminorticks=IntervalsBetween(2),
            xticklabelsvisible=row == 2,
            yticklabelsvisible=col == 1,
            yscale=log10
        )
        xlims!(ax, 1, 33)
        ylims!(ax, 1.1e-16, 50)
        axes[row, col] = ax
        haskey(s_unif, (i, j)) || continue
        ref = s_unif[(i, j)]
        scale = sqrt(mean_sq(abs.(ref)))
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
            hlines!(ax, [tol]; color=TOL_PALETTE[k], linewidth=TOL_LINE_LW, linestyle=:dash)
            pivots = collect_pivots(sampled_freqs, tol)
            if !isempty(pivots)
                ys = fill(2e-14 * tol^0.25, length(pivots))
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
        text!(
            ax,
            0.97,
            0.97;
            text=L"S_{%$i%$j}",
            space=:relative,
            align=(:right, :top),
            fontsize=SLABEL_SIZE
        )
    end
    linkaxes!(axes...)
    axislegend(
        axes[2, 2];
        position=:cb,
        framevisible=false,
        labelsize=LEGEND_SIZE,
        nbanks=3,
        unique=true,
        padding=(0, 0, 0, 0),
        rowgap=2,
        colgap=8
    )
    save(joinpath(outdir, "driven_ua_cpw_domain_sparam_adaptive_rms.svg"), fig)
    return fig
end

mean_sq(v) = sum(abs2, v) / length(v)
collect_pivots(d, tol) = haskey(d, tol) ? reduce(vcat, d[tol]; init=Float64[]) : Float64[]

# -------------------------- Plot 4: energy uniform ----------------------------

function plot_energy_uniform(freq, e, outdir)
    fig = Figure(size=(900, 600))
    ax = Axis(
        fig[1, 1];
        title=L"\text{Domain Energy } E_{\mathrm{elec}} \text{ (pJ)}",
        titlesize=SUPTITLE_SIZE,
        xlabel="f (GHz)",
        xticks=collect(range(2, 32; step=2)),
        xminorticks=IntervalsBetween(2)
    )
    xlims!(ax, 1, 33)
    scatterlines!(
        ax,
        freq,
        e["elec"] .* 1e12;
        color=REF_COLOR,
        linewidth=REF_LW,
        markersize=REF_MS,
        marker=:circle,
        label="Uniform Driven Solver"
    )
    axislegend(ax; position=:rt, framevisible=false, labelsize=LEGEND_SIZE)
    save(joinpath(outdir, "driven_ua_cpw_domain_energy_uniform.svg"), fig)
    return fig
end

# -------------------------- Plots 5/6: energy adaptive single & sweep ---------

function _plot_energy_adaptive(
    freq,
    e_unif,
    e_adapt,
    freq_adapt,
    sampled_freqs,
    outdir,
    fname;
    tols=ADAPTIVE_TOLS,
    ylim=(2e-14, 3),
    legend_position=:rt,
    legend_horizontal=false
)
    fig = Figure(size=(900, 600))
    ax = Axis(
        fig[1, 1];
        title="Error in Domain Energy",
        subtitle=L"|E_{\mathrm{elec,adaptive}} - E_{\mathrm{elec,uniform}}| / |E_{\mathrm{elec,uniform}}|",
        titlesize=SUPTITLE_SIZE,
        subtitlesize=SUPTITLE_SIZE,
        subtitlefont=:regular,
        xlabel="f (GHz)",
        xticks=collect(range(2, 32; step=2)),
        xminorticks=IntervalsBetween(2),
        yscale=log10
    )
    xlims!(ax, 1, 33)
    ylims!(ax, ylim...)
    ref = e_unif["elec"]
    for (k, tol) in enumerate(tols)
        haskey(e_adapt, tol) || continue
        err = abs.(e_adapt[tol]["elec"] .- ref) ./ abs.(ref)
        scatterlines!(
            ax,
            freq_adapt[tol],
            err;
            color=TOL_PALETTE[k],
            linewidth=ADAPT_LW,
            markersize=ADAPT_MS,
            marker=:circle,
            label="tol=$(fmt_tol(tol))"
        )
        hlines!(ax, [tol]; color=TOL_PALETTE[k], linestyle=:dash, linewidth=TOL_LINE_LW)
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
    hlines!(ax, [1e-12]; color=:black, linestyle=:dot, linewidth=REF_FLOOR_LW)
    if legend_horizontal
        axislegend(
            ax;
            position=legend_position,
            framevisible=false,
            labelsize=LEGEND_SIZE,
            orientation=:horizontal,
            nbanks=1,
            colgap=8,
            padding=(2, 2, 2, 2)
        )
    else
        axislegend(
            ax;
            position=legend_position,
            framevisible=false,
            labelsize=LEGEND_SIZE,
            padding=(2, 2, 2, 2)
        )
    end
    save(joinpath(outdir, fname), fig)
    return fig
end

plot_energy_adaptive_single(args...) = _plot_energy_adaptive(
    args...;
    tols=ADAPTIVE_TOLS[1:1],
    ylim=(2e-14, 3),
    legend_position=:rt
)
plot_energy_adaptive_sweep(args...) = _plot_energy_adaptive(
    args...;
    tols=ADAPTIVE_TOLS,
    ylim=(5e-15, 3),
    legend_position=:ct,
    legend_horizontal=true
)

# -------------------------- Plot 7: convergence curve -------------------------

function plot_convergence_curve(sampled_errors, outdir)
    fig = Figure(size=(900, 600))
    ax = Axis(
        fig[1, 1];
        title="Adaptive Solver Convergence: Error Indicator",
        xlabel="Sample number",
        xticks=collect(1:13),
        yscale=log10
    )
    last_errs = nothing
    # Iterate from tightest tol → loosest so that the looser tols draw on top,
    # matching matplotlib's `zorder=np.log(tol)` ordering in the Python script.
    for (k, tol) in reverse(collect(enumerate(ADAPTIVE_TOLS)))
        errs_list = get(sampled_errors, tol, Vector{Vector{Float64}}())
        isempty(errs_list) && continue
        errs = errs_list[1]
        finite = isfinite.(errs)
        any(finite) || continue
        n = collect(1:length(errs))
        scatterlines!(
            ax,
            n[finite],
            errs[finite];
            color=TOL_PALETTE[k],
            linewidth=1,
            markersize=PIVOT_MS,
            marker=:circle,
            label="tol=$(fmt_tol(tol))"
        )
        hlines!(ax, [tol]; color=TOL_PALETTE[k], linestyle=:dash, linewidth=TOL_LINE_LW)
        if tol == ADAPTIVE_TOLS[end]
            last_errs = errs   # tightest-tol curve, used for the leading dashed segment
        end
    end
    if last_errs !== nothing && length(last_errs) > 2
        scatterlines!(
            ax,
            [1, 2, 3],
            [1.0, 1.0, last_errs[3]];
            color=:lightgrey,
            marker=:star5,
            markersize=10,
            linestyle=:dash,
            linewidth=1
        )
    end
    axislegend(ax; position=:rt, nbanks=2, framevisible=false, labelsize=LEGEND_SIZE)
    save(joinpath(outdir, "driven_ua_cpw_adaptive_convergence_curve.svg"), fig)
    return fig
end

# -------------------------- Driver --------------------------------------------

function main(; outdir=default_outdir())
    set_theme!(makie_theme())
    mkpath(outdir)

    freq_u, s_u = load_port_s(UNIFORM_DIR)
    _, e_u = load_domain_e(UNIFORM_DIR)
    println("Uniform: $(length(freq_u)) frequency points")

    freq_a = Dict{Float64, Vector{Float64}}()
    s_a = Dict{Float64, Dict{Tuple{Int, Int}, Vector{ComplexF64}}}()
    e_a = Dict{Float64, Dict{String, Vector{Float64}}}()
    pivots = Dict{Float64, Vector{Vector{Float64}}}()
    errs = Dict{Float64, Vector{Vector{Float64}}}()
    for tol in ADAPTIVE_TOLS
        d = adaptive_dir(tol)
        if !isfile(joinpath(d, "port-S.csv"))
            @warn "Missing data for tol=$tol at $d"
            continue
        end
        freq_a[tol], s_a[tol] = load_port_s(d)
        _, e_a[tol] = load_domain_e(d)
        pivots[tol] = parse_log_pivots(joinpath(d, "palace.log"))
        errs[tol] = parse_log_errors(joinpath(d, "palace.log"))
        n = sum(length, pivots[tol]; init=0)
        println("Adaptive tol=$tol: $(length(freq_a[tol])) freq pts ($n pivots)")
    end

    plot_smat_uniform(freq_u, s_u, outdir)
    plot_smat_adaptive_pointwise(freq_u, s_u, s_a, freq_a, outdir)
    plot_smat_adaptive_rms(freq_u, s_u, s_a, freq_a, pivots, outdir)
    plot_energy_uniform(freq_u, e_u, outdir)
    plot_energy_adaptive_single(
        freq_u,
        e_u,
        e_a,
        freq_a,
        pivots,
        outdir,
        "driven_ua_cpw_domain_energy_adaptive_single.svg"
    )
    plot_energy_adaptive_sweep(
        freq_u,
        e_u,
        e_a,
        freq_a,
        pivots,
        outdir,
        "driven_ua_cpw_domain_energy_adaptive_sweep.svg"
    )
    plot_convergence_curve(errs, outdir)

    println("Wrote 7 SVGs to $outdir")
    return outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
    let outdir = default_outdir()
        args = copy(ARGS)
        while !isempty(args)
            a = popfirst!(args)
            if a == "--out"
                outdir = popfirst!(args)
            end
        end
        main(; outdir=outdir)
    end
end
