# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script runs the silicon carbide example `sic_reststrahlen.json` and compares the
infrared reflectance of a silicon carbide half-space computed by Palace with the Fresnel
formula for normal incidence,

    R = |(1 - n)/(1 + n)|²,   n = √ε,   ε(f) = ε∞ + Δε f₀² / (f₀² - f² + i f_γ f),

where ε is the Lorentz permittivity of the configuration file: the oscillator fit of
Spitzer, Kleinman and Walsh (1959) to the ordinary ray of 6H-SiC.

## Prerequisites

This script requires Julia packages. Install them with:

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

## How to run

From the repository root, run:
```bash
julia --project=examples -e 'include("examples/dielectric_grating/sic_reststrahlen.jl"); generate_reststrahlen_data()'
```

This requires `palace` to be a runnable command. If it is not, you can pass the path to the
executable, e.g.
```bash
julia --project=examples -e 'include("examples/dielectric_grating/sic_reststrahlen.jl"); generate_reststrahlen_data(palace_exec="build/bin/palace")'
```

The script will:
1. Run the adaptive frequency sweep of `sic_reststrahlen.json`
2. Print the largest difference of the S-parameters from the closed form, the number of
   full solves of the sweep, and the reflectance at a few wavelengths
3. Save the figures of the documentation to `postpro/sic_reststrahlen_reflectance.png`
   (reflectance against wavelength) and `postpro/sic_reststrahlen_permittivity.png` (the
   permittivity recovered from the reflection coefficient against the Lorentz model)

`plot_reststrahlen_data()` repeats steps 2 and 3 for the output of a previous run.
=#

using CSV
using DataFrames
using JSON
using Measures
using Plots
using Printf

# Geometry of mesh/sic_halfspace.msh with Model.L0 = 1e-6 (μm): port 1 lies in vacuum at
# d₁ from the surface of the silicon carbide, and port 2 inside it at depth d₂.
const d₁ = 0.5e-6
const d₂ = 0.25e-6
const c₀ = 299792458.0

"""
    read_config(file)

Read a configuration file from the example directory, stripping the // comments that
Palace accepts but JSON.jl does not.
"""
function read_config(file)
    return open(joinpath(@__DIR__, file)) do f
        return JSON.parse(replace(read(f, String), r"//[^\n]*" => ""))
    end
end

"""
    lorentz_permittivity(config)

Return the relative permittivity ε(f) (f in GHz) of the silicon carbide in `config`, a
single Lorentz term in Palace's positive-time convention, where a passive medium has
Im ε < 0.
"""
function lorentz_permittivity(config)
    model = config["Domains"]["Materials"][2]["Permittivity"]
    term = only(model["Terms"])
    @assert term["Type"] == "Lorentz"
    ε_inf = model["HighFrequency"]
    Δε = term["DeltaPermittivity"]
    f₀ = term["ResonanceFrequency"]
    f_γ = term["DampingFrequency"]
    return f -> ε_inf + Δε * f₀^2 / complex(f₀^2 - f^2, f_γ * f)
end

"""
    fresnel_sparameters(ε, f)

Floquet port S-parameters of the half-space at the frequency `f` in GHz, for the
permittivity `ε` at that frequency. The reflection coefficient of the surface is
r = (1 - n)/(1 + n), and the vacuum between port 1 and the surface adds the phase
exp(-2i k₀ d₁). Port 2 normalizes the transmitted wave to the power it carries, which in a
lossy medium is proportional to Re n.
"""
function fresnel_sparameters(ε, f)
    n = sqrt(ε)
    k₀ = 2π * f * 1e9 / c₀
    r = (1 - n) / (1 + n)
    S₁₁ = r * exp(-2im * k₀ * d₁)
    S₂₁ = 2 / (1 + n) * sqrt(real(n)) * exp(-im * k₀ * d₁) * exp(-im * k₀ * n * d₂)
    return S₁₁, S₂₁
end

"""
    read_sparameters(outdir)

Read the frequencies (GHz) and the complex TE and TM S-parameters of both ports from
`port-floquet-S.csv`.
"""
function read_sparameters(outdir)
    data = CSV.File(joinpath(outdir, "port-floquet-S.csv"), header=1) |> DataFrame
    names = strip.(String.(propertynames(data)))
    function column(mode)
        i = findfirst(==("|S[$mode][1]| (dB)"), names)
        return exp10.(data[!, i] ./ 20) .* cis.(deg2rad.(data[!, i + 1]))
    end
    return data[!, 1],
    column("P1(0;0)TE"),
    column("P2(0;0)TE"),
    max.(abs.(column("P1(0;0)TM")), abs.(column("P2(0;0)TM")))
end

λ_μm(f) = c₀ / (f * 1e9) * 1e6     # Wavelength in μm for f in GHz
wavenumber(f) = f * 1e9 / c₀ / 100  # cm⁻¹ for f in GHz

"""
    generate_reststrahlen_data(; palace_exec="palace", num_processors::Integer=1)

Run the silicon carbide example and compare the result with the Fresnel formula.

# Arguments

  - palace_exec - executable for Palace
  - num_processors - number of processors to use for the simulation
"""
function generate_reststrahlen_data(; palace_exec="palace", num_processors::Integer=1)
    is_path = occursin(Base.Filesystem.path_separator, palace_exec)
    palace_exec = is_path && !isabspath(palace_exec) ? abspath(palace_exec) : palace_exec
    run(Cmd(`$palace_exec -np $num_processors sic_reststrahlen.json`; dir=@__DIR__))
    return plot_reststrahlen_data()
end

"""
    plot_reststrahlen_data()

Compare the output of `sic_reststrahlen.json` with the Fresnel formula, print the
comparison, and save the figures.
"""
function plot_reststrahlen_data()
    config = read_config("sic_reststrahlen.json")
    outdir = joinpath(@__DIR__, config["Problem"]["Output"])
    ε = lorentz_permittivity(config)
    term = config["Domains"]["Materials"][2]["Permittivity"]
    ε_inf = term["HighFrequency"]
    Δε = term["Terms"][1]["DeltaPermittivity"]
    f_TO = term["Terms"][1]["ResonanceFrequency"]
    f_LO = f_TO * sqrt((ε_inf + Δε) / ε_inf)  # Re ε = 0 without damping

    f, S₁₁, S₂₁, TM = read_sparameters(outdir)
    exact = fresnel_sparameters.(ε.(f), f)
    ΔS₁₁ = abs.(S₁₁ .- first.(exact))
    ΔS₂₁ = abs.(S₂₁ .- last.(exact))
    R = abs2.(S₁₁)
    R_exact(f) = abs2((1 - sqrt(ε(f))) / (1 + sqrt(ε(f))))
    # Permittivity recovered from the reflection coefficient, n = (1 - r)/(1 + r)
    r = S₁₁ .* exp.(2im .* (2π .* f .* 1e9 ./ c₀) .* d₁)
    ε_palace = ((1 .- r) ./ (1 .+ r)) .^ 2
    Δε_rel = abs.(ε_palace .- ε.(f)) ./ abs.(ε.(f))
    metadata = JSON.parsefile(joinpath(outdir, "palace.json"))
    n_solves = metadata["LinearSolver"]["TotalSolves"]

    println()
    @printf(
        "%d frequencies from %.2f to %.2f μm, from %d full solves (%.1f s)\n",
        length(f),
        λ_μm(maximum(f)),
        λ_μm(minimum(f)),
        n_solves,
        metadata["ElapsedTime"]["Durations"]["Total"]
    )
    @printf(
        "Largest |ΔS₁₁| = %.1e at %.3f μm, |ΔS₂₁| = %.1e at %.3f μm, |S_TM| = %.1e\n",
        maximum(ΔS₁₁),
        λ_μm(f[argmax(ΔS₁₁)]),
        maximum(ΔS₂₁),
        λ_μm(f[argmax(ΔS₂₁)]),
        maximum(TM)
    )
    @printf(
        "Largest relative difference of ε from S₁₁: %.1e at %.3f μm\n",
        maximum(Δε_rel),
        λ_μm(f[argmax(Δε_rel)])
    )
    println("    λ (μm)   ν̃ (cm⁻¹)    f (GHz)   R, Palace      R, exact       Difference")
    for f_row ∈ [maximum(f), 29917.86, 29083.6, 26178.6, f_TO, minimum(f)]
        i = findfirst(==(f_row), f)
        @printf(
            "%10.3f %10.1f %10.2f   %.10f   %.10f   %+.1e\n",
            λ_μm(f[i]),
            wavenumber(f[i]),
            f[i],
            R[i],
            R_exact(f[i]),
            R[i] - R_exact(f[i])
        )
    end

    # Plot settings
    plotsz = (800, 500)
    fntsz = 12
    fnt = font(fntsz)
    default(
        size=plotsz,
        palette=:Set1_9,
        dpi=300,
        tickfont=fnt,
        guidefont=fnt,
        legendfontsize=fntsz - 2,
        margin=5mm
    )
    red, blue, grey = palette(:Set1_9)[1], palette(:Set1_9)[2], :grey50
    λ_TO, λ_LO = λ_μm(f_TO), λ_μm(f_LO)

    # Wavenumber on the top axis, at the wavelengths of round values
    function wavenumber_axis!(p, xlims, wavenumbers)
        top = twiny(p)
        plot!(
            top,
            xlims=xlims,
            xticks=(1e4 ./ wavenumbers, string.(wavenumbers)),
            xlabel="Wavenumber  (cm⁻¹)",
            grid=false
        )
        return p
    end

    # Reflectance on the axes of the fit in Spitzer et al. (1959), as redrawn in the NSM
    # Archive: wavelength from 1 to 22 μm and R from 0 to 1.
    λ_fine =
        range(1e6 * c₀ / (1e9 * maximum(f)), 1e6 * c₀ / (1e9 * minimum(f)), length=4000)
    f_fine = 1e-9 * c₀ ./ (1e-6 .* λ_fine)
    xlims = (1.0, 22.0)
    p1 = plot(
        xlabel="Wavelength  (μm)",
        ylabel="Reflectance \$R\$",
        xlims=xlims,
        ylims=(0.0, 1.0),
        xticks=5:5:20,
        yticks=0:0.2:1.0,
        legend=:right
    )
    vspan!(p1, [λ_LO, λ_TO], color=grey, alpha=0.15, linewidth=0, label="Reststrahlen band")
    plot!(p1, λ_fine, R_exact.(f_fine), color=red, label="Fresnel")
    scatter!(
        p1,
        λ_μm.(f),
        R,
        color=:black,
        markersize=2,
        markerstrokewidth=0,
        label="Palace"
    )
    annotate!(p1, λ_TO + 0.2, 0.96, text("TO", fntsz - 2, :left))
    annotate!(p1, λ_LO - 0.2, 0.96, text("LO", fntsz - 2, :right))
    wavenumber_axis!(p1, xlims, [500, 700, 1000, 2000, 5000])
    savefig(p1, joinpath(@__DIR__, "postpro", "sic_reststrahlen_reflectance.png"))

    # Permittivity recovered from the reflection coefficient, on an axis linear around zero
    # and logarithmic beyond: y = asinh(ε/2).
    y(x) = asinh(x / 2)
    ticks = [-500, -100, -10, 0, 10, 100, 500]
    xlims = (8.0, 15.0)
    in_window = xlims[1] .<= λ_μm.(f) .<= xlims[2]
    λ_fine = range(xlims..., length=4000)
    f_fine = 1e-9 * c₀ ./ (1e-6 .* λ_fine)
    p2 = plot(
        xlabel="Wavelength  (μm)",
        ylabel="Relative permittivity",
        xlims=xlims,
        ylims=y.((-800, 800)),
        yticks=(y.(ticks), replace.(string.(ticks), "-" => "−")),
        legend=:bottomleft
    )
    vspan!(p2, [λ_LO, λ_TO], color=grey, alpha=0.15, linewidth=0, label="")
    annotate!(p2, λ_TO + 0.08, y(400), text("TO", fntsz - 2, :left))
    annotate!(p2, λ_LO - 0.08, y(400), text("LO", fntsz - 2, :right))
    hline!(p2, [0.0], color=:black, linewidth=0.5, label="")
    plot!(
        p2,
        λ_fine,
        y.(real.(ε.(f_fine))),
        color=red,
        label="Lorentz, Re \$\\varepsilon\$"
    )
    plot!(
        p2,
        λ_fine,
        y.(imag.(ε.(f_fine))),
        color=blue,
        label="Lorentz, Im \$\\varepsilon\$"
    )
    scatter!(
        p2,
        λ_μm.(f[in_window]),
        y.(real.(ε_palace[in_window])),
        color=:black,
        markersize=2,
        markerstrokewidth=0,
        label="Palace, from \$S_{11}\$"
    )
    scatter!(
        p2,
        λ_μm.(f[in_window]),
        y.(imag.(ε_palace[in_window])),
        color=:black,
        markersize=2,
        markerstrokewidth=0,
        label=""
    )
    wavenumber_axis!(p2, xlims, [700, 800, 900, 1000, 1100, 1200])
    savefig(p2, joinpath(@__DIR__, "postpro", "sic_reststrahlen_permittivity.png"))

    return
end
