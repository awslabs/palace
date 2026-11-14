# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script runs the dispersive microstrip example and compares its scattering
parameters with an independent reference: the same line solved once per frequency with the
substrate replaced by the constant complex permittivity that the Djordjevic-Sarkar model
takes at that frequency,

    ε_r(f) = ε∞ + K log((f_hi + i f) / (f_lo + i f)),
    Permittivity = Re ε_r(f),  LossTan = -Im ε_r(f) / Re ε_r(f).

The reference uses the constant-material path of Palace, on the same mesh and with the same
ports, so the comparison isolates the frequency-dependent material: its assembly in the
volume and in the port mode solve, and the adaptive frequency sweep built on it.

## Prerequisites

This script requires Julia packages. Install them with:

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

## How to run

From the repository root, run:
```bash
julia --project=examples -e 'include("examples/microstrip/microstrip.jl"); generate_microstrip_data()'
```

This requires `palace` to be a runnable command. If it is not, you can pass the path to the
executable, e.g.
```bash
julia --project=examples -e 'include("examples/microstrip/microstrip.jl"); generate_microstrip_data(palace_exec="build/bin/palace")'
```

The script will:
1. Run `microstrip_driven.json` (uniform sweep) and `microstrip_adaptive.json` (adaptive
   sweep), both with the dispersive substrate
2. Run one single-frequency case with the constant permittivity ε_r(f) for each reference
   frequency (69 by default, about 30 s each on one process), and collect them in
   `postpro/microstrip_reference.csv`
3. Print the comparison at the uniform-sweep frequencies and the largest differences
4. Save the plot to `postpro/microstrip.png`

`plot_microstrip()` redraws the plot from the existing output without running Palace.
=#

using CSV
using DataFrames
using JSON
using Measures
using Plots
using Printf

const c₀ = 299792458.0

# Length of the line between the two port planes x = ±strip_length/2 of mesh/mesh.jl, with
# strip_length = 100 and Model.L0 = 1e-4 m
const ℓ_line = 10.0e-3

# Reference frequencies in GHz: the output step of microstrip_adaptive.json, refined around
# the |S₁₁| minimum near 7.7 GHz so that the reference curve resolves it
const reference_grid =
    sort(unique(round.(vcat(2.0:0.25:14.0, 6.5:0.125:9.0, 7.4:0.05:8.0), digits=4)))

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
    substrate_permittivity(config)

Return the function f ↦ ε_r(f), f in GHz, of the Djordjevic-Sarkar substrate (material 2)
of `config`, from its `Permittivity` object.
"""
function substrate_permittivity(config)
    model = config["Domains"]["Materials"][2]["Permittivity"]
    term = only(model["Terms"])
    @assert term["Type"] == "DjordjevicSarkar"
    ε_inf, K = model["HighFrequency"], term["Strength"]
    f_lo, f_hi = term["LowerFrequency"], term["UpperFrequency"]
    # ε∞ + K log((s + ω_hi) / (s + ω_lo)) at s = iω
    return f -> ε_inf + K * log(complex(f_hi, f) / complex(f_lo, f))
end

"""
    read_sparameters(dir)

Read `port-S.csv` from the output directory `dir`, returning the frequencies in GHz and the
complex S₁₁ and S₂₁.
"""
function read_sparameters(dir)
    data = CSV.read(joinpath(dir, "port-S.csv"), DataFrame)
    rename!(data, strip.(names(data)))
    S(m, p) = 10 .^ (m ./ 20) .* cis.(deg2rad.(p))
    return (
        f=data[!, "f (GHz)"],
        S11=S(data[!, "|S[1][1]| (dB)"], data[!, "arg(S[1][1]) (deg.)"]),
        S21=S(data[!, "|S[2][1]| (dB)"], data[!, "arg(S[2][1]) (deg.)"])
    )
end

"""
    read_substrate_participation(dir)

Read the substrate's share of the electric energy, `p_elec[1]`, from `domain-E.csv` in the
output directory `dir`.
"""
function read_substrate_participation(dir)
    data = CSV.read(joinpath(dir, "domain-E.csv"), DataFrame)
    rename!(data, strip.(names(data)))
    return data[!, "p_elec[1]"]
end

"""
    run_palace(config, palace_exec, num_processors)

Run Palace on `config` from the example directory through a temporary configuration file.
"""
function run_palace(config, palace_exec, num_processors)
    tmp_path = joinpath(@__DIR__, "microstrip_tmp.json")
    open(tmp_path, "w") do f
        return JSON.print(f, config, 4)
    end
    run(Cmd(`$palace_exec -np $num_processors $tmp_path`; dir=@__DIR__))
    rm(tmp_path)
    return
end

"""
    line_parameters(f, S21, ε)

Effective permittivity ε_eff = (c φ / ω ℓ)² and loss per radian 2α/β of the line from its
transmission coefficient S₂₁ at frequency `f` in GHz, with α ℓ = -ln |S₂₁| and β ℓ = φ. The
2π branch of the phase φ = -arg S₂₁ is the one nearest to the phase ω ℓ √(Re ε_r) / c of
a line filled with the substrate, which is unambiguous since ε_eff < Re ε_r.
"""
function line_parameters(f, S21, ε)
    k0ℓ = 2π * f * 1e9 * ℓ_line / c₀
    φ = -angle(S21)
    φ += 2π * round((k0ℓ * sqrt(real(ε)) - φ) / 2π)
    return (φ / k0ℓ)^2, -2 * log(abs(S21)) / φ
end

"""
    generate_microstrip_data(;
                             palace_exec::String="palace",
                             num_processors::Integer=1,
                             reference_frequencies=reference_grid
                            )

Generate the data for the dispersive microstrip example and visualize the results.

# Arguments

  - palace_exec - executable for Palace
  - num_processors - number of processors to use for the simulation
  - reference_frequencies - frequencies in GHz of the constant-permittivity reference
"""
function generate_microstrip_data(;
    palace_exec="palace",
    num_processors::Integer=1,
    reference_frequencies=reference_grid
)
    is_path = occursin(Base.Filesystem.path_separator, palace_exec)
    palace_exec = is_path && !isabspath(palace_exec) ? abspath(palace_exec) : palace_exec

    # The dispersive line: uniform and adaptive sweeps of the supplied configurations
    for file ∈ ("microstrip_driven.json", "microstrip_adaptive.json")
        run(Cmd(`$palace_exec -np $num_processors $file`; dir=@__DIR__))
    end

    # The reference: one constant-permittivity solve per frequency, otherwise identical to
    # the uniform sweep
    base = read_config("microstrip_driven.json")
    ε = substrate_permittivity(base)
    S11, S21 = ComplexF64[], ComplexF64[]
    for f ∈ reference_frequencies
        config = deepcopy(base)
        material = config["Domains"]["Materials"][2]
        material["Permittivity"] = real(ε(f))
        material["LossTan"] = -imag(ε(f)) / real(ε(f))
        config["Problem"]["Output"] = "postpro/reference"
        config["Solver"]["Driven"] =
            Dict("MinFreq" => f, "MaxFreq" => f, "FreqStep" => 1.0, "SaveStep" => 0)
        run_palace(config, palace_exec, num_processors)
        reference = read_sparameters(joinpath(@__DIR__, "postpro", "reference"))
        @assert only(reference.f) ≈ f
        push!(S11, only(reference.S11))
        push!(S21, only(reference.S21))
    end
    CSV.write(
        joinpath(@__DIR__, "postpro", "microstrip_reference.csv"),
        DataFrame(
            "f (GHz)" => reference_frequencies,
            "Re{eps_r}" => real.(ε.(reference_frequencies)),
            "tan(delta)" =>
                -imag.(ε.(reference_frequencies)) ./ real.(ε.(reference_frequencies)),
            "Re{S[1][1]}" => real.(S11),
            "Im{S[1][1]}" => imag.(S11),
            "Re{S[2][1]}" => real.(S21),
            "Im{S[2][1]}" => imag.(S21)
        )
    )

    print_microstrip_comparison()
    plot_microstrip()
    return
end

"""
    read_reference(root)

Read the constant-permittivity reference written by `generate_microstrip_data`.
"""
function read_reference(root)
    data = CSV.read(joinpath(root, "postpro", "microstrip_reference.csv"), DataFrame)
    return (
        f=data[!, "f (GHz)"],
        S11=data[!, "Re{S[1][1]}"] .+ im .* data[!, "Im{S[1][1]}"],
        S21=data[!, "Re{S[2][1]}"] .+ im .* data[!, "Im{S[2][1]}"]
    )
end

"""
    difference_to_reference(sweep, reference)

Largest complex difference max(|ΔS₁₁|, |ΔS₂₁|) between a sweep and the reference at each
frequency of the sweep, which must be a reference frequency.
"""
function difference_to_reference(sweep, reference)
    return map(enumerate(sweep.f)) do (k, f)
        j = findfirst(x -> isapprox(x, f; atol=1e-9), reference.f)
        isnothing(j) && error("No reference at $f GHz")
        return max(
            abs(sweep.S11[k] - reference.S11[j]),
            abs(sweep.S21[k] - reference.S21[j])
        )
    end
end

"""
    print_microstrip_comparison(; root=@__DIR__)

Print the line parameters at the uniform-sweep frequencies and the agreement of both sweeps
with the reference.
"""
function print_microstrip_comparison(; root::AbstractString=@__DIR__)
    ε = substrate_permittivity(read_config("microstrip_driven.json"))
    uniform = read_sparameters(joinpath(root, "postpro", "driven"))
    adaptive = read_sparameters(joinpath(root, "postpro", "adaptive"))
    reference = read_reference(root)
    p_e = read_substrate_participation(joinpath(root, "postpro", "driven"))
    Δ_uniform = difference_to_reference(uniform, reference)
    Δ_adaptive = difference_to_reference(adaptive, reference)
    dB(S) = 20 * log10(abs(S))

    println()
    println(
        " f (GHz)  Re ε_r   tan δ     ε_eff    2α/β      p_e tan δ  |S11| (dB)  " *
        "|S21| (dB)  |ΔS|"
    )
    for (k, f) ∈ enumerate(uniform.f)
        ε_eff, loss = line_parameters(f, uniform.S21[k], ε(f))
        tanδ = -imag(ε(f)) / real(ε(f))
        @printf(
            "%7.2f  %7.4f  %.5f  %7.4f  %.5f  %.5f   %8.2f   %9.4f   %.1e\n",
            f,
            real(ε(f)),
            tanδ,
            ε_eff,
            loss,
            p_e[k] * tanδ,
            dB(uniform.S11[k]),
            dB(uniform.S21[k]),
            Δ_uniform[k]
        )
    end
    fine = reference.f[1]:0.01:reference.f[end]
    println()
    @printf("Largest |ΔS|, uniform sweep vs reference:  %.1e\n", maximum(Δ_uniform))
    @printf("Largest |ΔS|, adaptive sweep vs reference: %.1e\n", maximum(Δ_adaptive))
    k = argmax(abs.(reference.S11))
    @printf(
        "Reference |S11|: largest %.2f dB at %.3f GHz, smallest %.2f dB at %.3f GHz\n",
        dB(reference.S11[k]),
        reference.f[k],
        minimum(dB.(reference.S11)),
        reference.f[argmin(abs.(reference.S11))]
    )
    @printf(
        "Largest |S11|^2 + |S21|^2, uniform sweep: %.4f\n",
        maximum(abs2.(uniform.S11) .+ abs2.(uniform.S21))
    )
    solves =
        JSON.parsefile(joinpath(root, "postpro", "adaptive", "palace.json"))["LinearSolver"]["TotalSolves"]
    println("Full-order solves of the adaptive sweep: $solves")
    return
end

"""
    plot_microstrip(; root=@__DIR__)

Plot the scattering parameters and the line parameters derived from S₂₁ of both dispersive
sweeps against the constant-permittivity reference. Run `generate_microstrip_data` first.
The generated plot is `postpro/microstrip.png`.
"""
function plot_microstrip(; root::AbstractString=@__DIR__)
    ε = substrate_permittivity(read_config("microstrip_driven.json"))
    uniform = read_sparameters(joinpath(root, "postpro", "driven"))
    adaptive = read_sparameters(joinpath(root, "postpro", "adaptive"))
    reference = read_reference(root)
    dB(S) = 20 * log10(abs(S))
    line(sweep) = [line_parameters(f, S, ε(f)) for (f, S) in zip(sweep.f, sweep.S21)]

    # Plot settings
    plotsz = (900, 700)
    fntsz = 12
    fnt = font(fntsz)
    default(
        size=plotsz,
        palette=:Set1_9,
        dpi=300,
        tickfont=fnt,
        guidefont=fnt,
        legendfontsize=fntsz - 2,
        margin=3mm
    )
    blue, grey = palette(:Set1_9)[2], :grey45
    xaxis = (xlims=(1.4, 14.6), xticks=2:2:14)

    # Reference as a line, the uniform sweep as rings, and the adaptive sweep as dots
    function compare!(p, x, y_reference, y_uniform, y_adaptive; labels=false)
        plot!(
            p,
            x(reference),
            y_reference;
            color=blue,
            linewidth=1.8,
            label=labels ? "Reference, constant \$\\varepsilon_r(f)\$" : ""
        )
        scatter!(
            p,
            x(uniform),
            y_uniform;
            markershape=:circle,
            markersize=6,
            markercolor=:white,
            markerstrokecolor=:black,
            markerstrokewidth=1.2,
            label=labels ? "Dispersive, uniform sweep" : ""
        )
        return scatter!(
            p,
            x(adaptive),
            y_adaptive;
            color=:black,
            markersize=2.5,
            markerstrokewidth=0,
            label=labels ? "Dispersive, adaptive sweep" : ""
        )
    end
    freq(sweep) = sweep.f

    p1 = plot(; xaxis..., xformatter=_ -> "", ylabel="\$|S_{21}|\$  (dB)", legend=false)
    compare!(p1, freq, dB.(reference.S21), dB.(uniform.S21), dB.(adaptive.S21))

    p2 = plot(; xaxis..., xformatter=_ -> "", ylabel="\$|S_{11}|\$  (dB)", legend=false)
    compare!(p2, freq, dB.(reference.S11), dB.(uniform.S11), dB.(adaptive.S11))

    # Effective permittivity and loss per radian of the line, with the substrate's. The one
    # legend sits in the empty band between the two permittivities.
    f_fine = range(2.0, 14.0, length=200)
    p3 = plot(;
        xaxis...,
        xlabel="Frequency  (GHz)",
        ylabel="Permittivity",
        legend=:left,
        ylims=(3.88, 4.42)
    )
    plot!(p3, f_fine, real.(ε.(f_fine)); color=grey, linestyle=:dash, label="")
    compare!(
        p3,
        freq,
        first.(line(reference)),
        first.(line(uniform)),
        first.(line(adaptive));
        labels=true
    )
    annotate!(
        p3,
        13.8,
        4.29,
        text("\$\\mathrm{Re}\\,\\varepsilon_r\$", fntsz - 1, grey, :right)
    )
    annotate!(
        p3,
        13.8,
        3.99,
        text("\$\\varepsilon_{\\mathrm{eff}}\$", fntsz - 1, blue, :right)
    )

    p4 = plot(;
        xaxis...,
        xlabel="Frequency  (GHz)",
        ylabel="Loss tangent",
        legend=false,
        ylims=(0.01960, 0.02070),
        yticks=0.0196:0.0002:0.0206,
        yformatter=y -> @sprintf("%.4f", y)
    )
    plot!(p4, f_fine, -imag.(ε.(f_fine)) ./ real.(ε.(f_fine)); color=grey, linestyle=:dash)
    compare!(p4, freq, last.(line(reference)), last.(line(uniform)), last.(line(adaptive)))
    annotate!(p4, 13.8, 0.02049, text("\$\\tan\\delta\$", fntsz - 1, grey, :right))
    annotate!(p4, 13.8, 0.02012, text("\$2\\alpha/\\beta\$", fntsz - 1, blue, :right))

    pp = plot(p1, p2, p3, p4; layout=(2, 2))
    savefig(pp, joinpath(root, "postpro", "microstrip.png"))
    display(pp)
    return
end
