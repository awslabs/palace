# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# Electrostatic MMS documentation plot generator

Regenerate the electrostatic method-of-manufactured-solutions figure from the
committed verification data:

```bash
julia --project=examples docs/generate_electrostatic_mms_plots.jl
```

Use `--data <file>` to read another CSV and `--out <dir>` to write the SVG
somewhere other than `docs/src/assets/verification`.
=#

using CSV
using CairoMakie
using DataFrames

const DOCS_DIR = @__DIR__

default_data_path() =
    joinpath(DOCS_DIR, "src", "assets", "verification", "electrostatic_mms.csv")
default_outdir() = joinpath(DOCS_DIR, "src", "assets", "verification")

const ORDER_COLORS = Dict(1 => "#2166ac", 2 => "#d6604d", 3 => "#1b7837")
const CURVED_COLOR = "#6a3d9a"
const HEX_COLOR = "#b2182b"
const TET_COLOR = "#006837"

function mms_theme()
    return Theme(
        fontsize=18,
        fonts=(
            regular="Times",
            italic="Times Italic",
            bold="Times Bold",
            bold_italic="Times Bold Italic"
        ),
        Axis=(
            xgridcolor=(:black, 0.12),
            ygridcolor=(:black, 0.12),
            rightspinevisible=false,
            topspinevisible=false,
            xminorticksvisible=true,
            yminorticksvisible=true,
            xticklabelsize=17,
            yticklabelsize=17,
            xlabelsize=19,
            ylabelsize=19,
            titlesize=21
        ),
        Legend=(framevisible=false, labelsize=16, patchsize=(25, 20), rowgap=3)
    )
end

expected_power(order, field) = field == :potential ? order + 1 : order

function reference_error(dofs, error, order, field)
    power = expected_power(order, field)
    return error[1] .* (dofs ./ dofs[1]) .^ (-power / 3)
end

function make_axis(fig, col, title, ylims)
    ax = Axis(
        fig[1, col];
        title,
        xlabel="Global H1 degrees of freedom",
        ylabel="Absolute L2 error",
        xscale=log10,
        yscale=log10
    )
    xlims!(ax, 90, 1.6e5)
    ylims!(ax, ylims...)
    return ax
end

function plot_series!(ax, dofs, error; color, marker, label=nothing)
    scatterlines!(
        ax,
        dofs,
        error;
        color,
        marker,
        markersize=12,
        strokecolor=:white,
        strokewidth=1,
        linewidth=2.5,
        label
    )
    return
end

function plot_reference!(ax, dofs, error, order, field, color)
    lines!(
        ax,
        dofs,
        reference_error(dofs, error, order, field);
        color,
        linestyle=:dot,
        linewidth=2
    )
    return
end

function plot_affine_inset!(fig, col, affine, field)
    yticks =
        field == :potential ? ([1e-14, 1e-13], ["1e-14", "1e-13"]) :
        ([1e-13, 1e-12], ["1e-13", "1e-12"])
    inset = Axis(
        fig[1, col];
        width=Relative(0.30),
        height=Relative(0.28),
        halign=0.97,
        valign=0.96,
        tellwidth=false,
        tellheight=false,
        title="affine p = 2",
        titlesize=12,
        xgridvisible=false,
        ygridcolor=(:black, 0.10),
        xticks=([1, 2], ["hex", "tet"]),
        xticklabelsize=11,
        yticks,
        yticklabelsize=10,
        xminorticksvisible=false,
        yminorticksvisible=false,
        rightspinevisible=true,
        topspinevisible=true,
        backgroundcolor=(:white, 0.94),
        yscale=log10
    )

    for (x, element, color, marker) in
        ((1, "hexahedron", HEX_COLOR, :rect), (2, "tetrahedron", TET_COLOR, :utriangle))
        row = only(eachrow(filter(row -> row.element == element, affine)))
        error = field == :potential ? row.potential_l2_error : row.electric_l2_error
        scatter!(inset, [x], [error]; color, marker, markersize=11)
    end

    xlims!(inset, 0.5, 2.5)
    if field == :potential
        ylims!(inset, 1e-14, 1e-13)
    else
        ylims!(inset, 1e-13, 1e-12)
    end
    translate!(inset.blockscene, 0, 0, 100)
    return inset
end

function plot_convergence(data, outdir)
    smooth = filter(row -> row.dataset == "cartesian-smooth", data)
    curved = sort(filter(row -> row.dataset == "curved-polynomial", data), :primary_dofs)
    affine = filter(row -> row.dataset == "affine-polynomial", data)
    fig = Figure(size=(1150, 520), figure_padding=(8, 18, 8, 8))
    ax_v = make_axis(fig, 1, "(a) Potential V", (1e-6, 1.0))
    ax_e = make_axis(fig, 2, "(b) Electric field E = -grad V", (1e-4, 10.0))

    for order = 1:3
        rows = sort(filter(row -> row.order == order, smooth), :primary_dofs)
        dofs = Float64.(rows.primary_dofs)
        potential = Float64.(rows.potential_l2_error)
        electric = Float64.(rows.electric_l2_error)
        color = ORDER_COLORS[order]

        plot_series!(
            ax_v,
            dofs,
            potential;
            color,
            marker=:rect,
            label="smooth Cartesian, p = $order"
        )
        plot_reference!(ax_v, dofs, potential, order, :potential, color)
        plot_series!(ax_e, dofs, electric; color, marker=:rect)
        plot_reference!(ax_e, dofs, electric, order, :electric, color)
    end

    dofs = Float64.(curved.primary_dofs)
    potential = Float64.(curved.potential_l2_error)
    electric = Float64.(curved.electric_l2_error)
    plot_series!(
        ax_v,
        dofs,
        potential;
        color=CURVED_COLOR,
        marker=:diamond,
        label="curved cylinder, p = 2"
    )
    plot_reference!(ax_v, dofs, potential, 2, :potential, CURVED_COLOR)
    plot_series!(ax_e, dofs, electric; color=CURVED_COLOR, marker=:diamond)
    plot_reference!(ax_e, dofs, electric, 2, :electric, CURVED_COLOR)

    plot_affine_inset!(fig, 1, affine, :potential)
    plot_affine_inset!(fig, 2, affine, :electric)
    axislegend(ax_v; position=:lb)
    Label(
        fig[2, 1:2],
        "Dotted guides start at each coarsest point; expected slopes are " *
        "-(p+1)/3 for V and -p/3 for E.";
        fontsize=16,
        color=(:black, 0.62),
        padding=(0, 0, 0, 0)
    )
    colgap!(fig.layout, 35)
    rowsize!(fig.layout, 2, 30)

    save(joinpath(outdir, "electrostatic-mms-convergence.svg"), fig)
    return fig
end

function main(; data_path=default_data_path(), outdir=default_outdir())
    set_theme!(mms_theme())
    mkpath(outdir)
    data = CSV.read(data_path, DataFrame)
    plot_convergence(data, outdir)
    println("Wrote electrostatic MMS convergence SVG to $outdir")
    return outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
    let data_path = default_data_path(), outdir = default_outdir()
        args = copy(ARGS)
        while !isempty(args)
            arg = popfirst!(args)
            if arg == "--data"
                data_path = popfirst!(args)
            elseif arg == "--out"
                outdir = popfirst!(args)
            else
                error("Unknown argument: $arg")
            end
        end
        main(; data_path, outdir)
    end
end
