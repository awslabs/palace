# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script sweeps the resonance frequency of the Lorentz medium filling the
cylindrical cavity of `cavity_lorentz.json` through the TM010 and TE111 modes, and
compares the polariton eigenfrequencies computed by Palace with the roots of the
dispersion relation

    (Ω² + ε∞ s²)(s² + γ s + ω₀²) + Δε ω₀² s² = 0,   s = iω,

where Ω = √ε∞ ω_c is set by the frequency ω_c of a mode of the cavity filled with ε∞.

## Prerequisites

This script requires Julia packages. Install them with:

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

## How to run

From the repository root, run:
```bash
julia --project=examples -e 'include("examples/cylinder/cavity_polariton.jl"); generate_polariton_data()'
```

This requires `palace` to be a runnable command. If it is not, you can pass the path to the
executable, e.g.
```bash
julia --project=examples -e 'include("examples/cylinder/cavity_polariton.jl"); generate_polariton_data(palace_exec="build/bin/palace")'
```

The script will:
1. Solve the cavity filled with ε∞ alone, for the frequencies of the uncoupled modes on the
   example mesh
2. For each resonance frequency, solve for the lowest modes below and above the resonance
3. Print the computed and analytic polariton frequencies at the TM010 resonance
4. Save a plot of the polariton branches against the resonance frequency to
   `postpro/cavity_polariton.png`
=#

using CSV
using DataFrames
using JSON
using LinearAlgebra
using Measures
using Plots
using Printf

# Cavity of examples/cylinder/mesh with Model.L0 = 1e-2 (cm): radius a and height d = 2a.
const a_cyl = 2.74e-2
const d_cyl = 2 * a_cyl
const c₀ = 299792458.0

"""
    cavity_frequency(p, l, ε)

Frequency in GHz of the TE or TM mode of the PEC cylinder filled with permittivity `ε`, with
transverse wavenumber `p / a` and `l` half wavelengths along the axis [Pozar, (6.62)].
"""
cavity_frequency(p, l, ε) = c₀ / (2π * sqrt(ε)) * hypot(p / a_cyl, l * π / d_cyl) / 1e9

# First zero of J₀ (TM010) and of J₁′ (TE111)
const p₀₁ = 2.404825557695773
const p′₁₁ = 1.841183781340659

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
    run_eigenmode(config, palace_exec, num_processors)

Run Palace on `config` from the example directory through a temporary configuration file,
and return the computed complex eigenfrequencies in GHz.
"""
function run_eigenmode(config, palace_exec, num_processors)
    tmp_path = joinpath(@__DIR__, "cavity_polariton_tmp.json")
    open(tmp_path, "w") do f
        return JSON.print(f, config, 4)
    end
    run(Cmd(`$palace_exec -np $num_processors $tmp_path`; dir=@__DIR__))
    rm(tmp_path)
    eig = CSV.File(
        joinpath(@__DIR__, config["Problem"]["Output"], "eig.csv"),
        header=1,
        normalizenames=true
    )
    data = DataFrame(eig) |> Matrix
    return data[:, 2] .+ im .* data[:, 3]
end

"""
    polariton_frequencies(f_c, f_0; ε_inf, Δε, f_γ)

Complex frequencies in GHz of the two branches into which a cavity mode of frequency `f_c`
(for the cavity filled with `ε_inf`) splits in a Lorentz medium with resonance `f_0`,
strength `Δε`, and damping `f_γ`, all frequencies in GHz. Returns `(lower, upper)`.
"""
function polariton_frequencies(f_c, f_0; ε_inf, Δε, f_γ)
    Ω2 = ε_inf * (2π * f_c)^2
    ω0 = 2π * f_0
    γ = 2π * f_γ
    # (Ω² + ε∞ s²)(s² + γ s + ω₀²) + Δε ω₀² s² as a polynomial in s, highest degree first
    c4 = ε_inf
    c3 = ε_inf * γ
    c2 = Ω2 + ε_inf * ω0^2 + Δε * ω0^2
    c1 = Ω2 * γ
    c0 = Ω2 * ω0^2
    # Roots from the eigenvalues of the companion matrix
    C = [
        -c3/c4 -c2/c4 -c1/c4 -c0/c4
        1 0 0 0
        0 1 0 0
        0 0 1 0
    ]
    s = eigvals(complex(C))
    f = [z / (2π * im) for z in s if imag(z) > 0]  # Re(f) > 0
    sort!(f, by=real)
    return f[1], f[2]
end

"""
    generate_polariton_data(;
                            palace_exec::String="palace",
                            num_processors::Integer=1,
                            resonances=range(2.55, 3.30, length=31)
                           )

Generate the data for the cavity polariton example and visualize the results.

# Arguments

  - palace_exec - executable for Palace
  - num_processors - number of processors to use for the simulation
  - resonances - resonance frequencies of the Lorentz medium to sweep, in GHz
"""
function generate_polariton_data(;
    palace_exec="palace",
    num_processors::Integer=1,
    resonances=range(2.55, 3.30, length=31)
)
    is_path = occursin(Base.Filesystem.path_separator, palace_exec)
    palace_exec = is_path && !isabspath(palace_exec) ? abspath(palace_exec) : palace_exec

    base = read_config("cavity_lorentz.json")
    material = base["Domains"]["Materials"][1]
    ε_inf = material["Permittivity"]["HighFrequency"]
    term = material["Permittivity"]["Terms"][1]
    Δε = term["DeltaPermittivity"]
    f_γ = term["DampingFrequency"]
    f_res = term["ResonanceFrequency"]

    # Uncoupled TM010, TE111, and TM011 modes of the cavity filled with ε∞: analytic, and
    # computed on the example mesh
    f_tm010 = cavity_frequency(p₀₁, 0, ε_inf)
    f_te111 = cavity_frequency(p′₁₁, 1, ε_inf)
    f_tm011 = cavity_frequency(p₀₁, 1, ε_inf)
    bare = deepcopy(base)
    bare["Domains"]["Materials"][1]["Permittivity"] = ε_inf
    bare["Problem"]["Output"] = "postpro/cavity_polariton_bare"
    bare["Solver"]["Eigenmode"]["N"] = 4
    bare["Solver"]["Eigenmode"]["Target"] = 2.5
    f_bare = real.(run_eigenmode(bare, palace_exec, num_processors))
    f_tm010_h, f_te111_h, f_tm011_h = f_bare[1], f_bare[2], f_bare[4]

    # Polariton modes for each resonance: one solve with a target below the lower TM010
    # branch, and one with a target above the gap (f₀, f₀√(1+Δε/ε∞)), where the
    # permittivity is negative, below the upper TM010 branch. Every computed eigenvalue is
    # kept, and matched afterwards to the nearest analytic branch.
    computed = Tuple{Float64, ComplexF64}[]
    for f_0 ∈ resonances
        f_L = f_0 * sqrt(1 + Δε / ε_inf)
        tm_lower, tm_upper = polariton_frequencies(f_tm010, f_0; ε_inf, Δε, f_γ)
        for target ∈ (0.97 * real(tm_lower), 0.5 * (f_L + real(tm_upper)))
            config = deepcopy(base)
            config["Domains"]["Materials"][1]["Permittivity"]["Terms"][1]["ResonanceFrequency"] =
                f_0
            config["Problem"]["Output"] = "postpro/cavity_polariton"
            config["Solver"]["Eigenmode"]["N"] = 3
            config["Solver"]["Eigenmode"]["Target"] = target
            for f ∈ run_eigenmode(config, palace_exec, num_processors)
                push!(computed, (f_0, f))
            end
        end
    end

    # Comparison at the resonance of cavity_lorentz.json, whose lowest four modes are the
    # lower TM010 branch, the TE111 pair, and the lower TM011 branch
    config = deepcopy(base)
    config["Problem"]["Output"] = "postpro/cavity_polariton"
    f_palace = run_eigenmode(config, palace_exec, num_processors)
    fmt(f) = @sprintf("%.6f%+.6fi", real(f), imag(f))
    println()
    println(
        "Mode     Uncoupled, exact  Uncoupled, mesh  Palace                Dispersion, " *
        "exact          Dispersion, mesh"
    )
    for (name, f_c, f_c_h, f) ∈ (
        ("TM010", f_tm010, f_tm010_h, f_palace[1]),
        ("TE111", f_te111, f_te111_h, f_palace[2]),
        ("TM011", f_tm011, f_tm011_h, f_palace[4])
    )
        f_exact = polariton_frequencies(f_c, f_res; ε_inf, Δε, f_γ)[1]
        f_mesh = polariton_frequencies(f_c_h, f_res; ε_inf, Δε, f_γ)[1]
        println(
            rpad(name, 9),
            rpad(@sprintf("%.6f", f_c), 18),
            rpad(@sprintf("%.6f", f_c_h), 17),
            rpad(fmt(f), 22),
            rpad(fmt(f_exact), 27),
            fmt(f_mesh)
        )
    end
    f_split = polariton_frequencies(f_tm010_h, f_res; ε_inf, Δε, f_γ)
    println(
        "TM010 splitting: $(round(real(f_split[2] - f_split[1]), digits=6)) GHz, " *
        "f₀ √(Δε/ε∞) = $(round(f_res * sqrt(Δε / ε_inf), digits=6)) GHz"
    )
    # Nearest analytic branch, using the uncoupled frequencies from the same mesh
    names = ("TM010-", "TM010+", "TE111-", "TE111+")
    function nearest_branch(f_0, f)
        tm = polariton_frequencies(f_tm010_h, f_0; ε_inf, Δε, f_γ)
        te = polariton_frequencies(f_te111_h, f_0; ε_inf, Δε, f_γ)
        candidates = (tm[1], tm[2], te[1], te[2])
        k = argmin(abs.(f .- candidates))
        return names[k], abs(f - candidates[k]) / abs(candidates[k])
    end
    println("Resonance     Computed eigenvalue (GHz)   Nearest branch  Relative difference")
    for (f_0, f) ∈ computed
        name, err = nearest_branch(f_0, f)
        println(
            lpad(round(f_0, digits=4), 9),
            lpad(string(round(f, digits=6)), 30),
            lpad(name, 16),
            lpad(@sprintf("%.1e", err), 21)
        )
    end

    # Plot settings
    plotsz = (800, 760)
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

    f_fine = range(first(resonances), last(resonances), length=400)
    branches(f_c) = [polariton_frequencies(f_c, f; ε_inf, Δε, f_γ) for f in f_fine]
    tm = branches(f_tm010_h)
    te = branches(f_te111_h)
    blue, red, grey = palette(:Set1_9)[2], palette(:Set1_9)[1], :grey50

    p1 = plot(
        xlabel="",
        ylabel="Re \$f\$  (GHz)",
        legend=:topleft,
        xlims=extrema(resonances),
        ylims=(2.50, 3.45)
    )
    plot!(
        p1,
        f_fine,
        f_fine,
        color=grey,
        linestyle=:dot,
        label="Gap edges, \$f_0\$ and \$f_L\$"
    )
    plot!(p1, f_fine, f_fine .* sqrt(1 + Δε / ε_inf), color=grey, linestyle=:dot, label="")
    hline!(p1, [f_tm010_h], color=blue, linestyle=:dash, alpha=0.6, label="Uncoupled TM010")
    hline!(p1, [f_te111_h], color=red, linestyle=:dash, alpha=0.6, label="Uncoupled TE111")
    for k ∈ 1:2
        plot!(
            p1,
            f_fine,
            real.(getindex.(tm, k)),
            color=blue,
            label=(k == 1 ? "TM010" : "")
        )
        plot!(p1, f_fine, real.(getindex.(te, k)), color=red, label=(k == 1 ? "TE111" : ""))
    end
    scatter!(
        p1,
        first.(computed),
        real.(last.(computed)),
        color=:black,
        markersize=3,
        label="Palace"
    )

    p2 = plot(
        xlabel="Resonance frequency \$f_0\$  (GHz)",
        ylabel="Im \$f\$  (MHz)",
        legend=:top,
        xlims=extrema(resonances),
        ylims=(0.0, 14.0)
    )
    hline!(
        p2,
        [1e3 * f_γ / 2],
        color=grey,
        linestyle=:dot,
        label="Oscillator, \$f_\\gamma/2\$"
    )
    hline!(
        p2,
        [1e3 * f_γ / 4],
        color=grey,
        linestyle=:dash,
        label="Equal mixing, \$f_\\gamma/4\$"
    )
    for k ∈ 1:2
        plot!(p2, f_fine, 1e3 .* imag.(getindex.(tm, k)), color=blue, label="")
        plot!(p2, f_fine, 1e3 .* imag.(getindex.(te, k)), color=red, label="")
    end
    scatter!(
        p2,
        first.(computed),
        1e3 .* imag.(last.(computed)),
        color=:black,
        markersize=3,
        label="Palace"
    )

    pp = plot(p1, p2, layout=grid(2, 1, heights=[0.65, 0.35]))
    savefig(pp, joinpath(@__DIR__, "postpro", "cavity_polariton.png"))
    display(pp)

    return
end
