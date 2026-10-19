# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script reproduces the picosecond pulse propagation on the Nb microstrip of
Talanov et al. (2022). It runs the 2D boundary mode problem of the line cross-section over
a frequency sweep, for the Mattis-Bardeen complex conductivity of the Nb films meshed as
domains or modeled as surface impedance boundaries, and for the two-fluid model. The
propagation constant kₙ(ω) gives the propagation of a single-flux-quantum pulse.
Optionally, 3D lines of two lengths between wave ports give kₙ from the ratio of their
transmission coefficients.

## Prerequisites

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

## How to run

From the repository root:
```bash
julia --project=examples -e 'include("examples/superconducting_line/superconducting_line.jl"); generate_superconducting_line_data()'
```

This requires `palace` to be a runnable command, or pass the executable with
`palace_exec`. The script writes `postpro/mode_sweep.csv`, `postpro/waveport_kn.csv` (with
`waveport = true`), and the plots `postpro/superconducting_line-*.png`.
=#

using CSV, DataFrames, JSON, CairoMakie

include(joinpath(@__DIR__, "mesh", "mesh.jl"))

const c0 = 299792458.0

# Two-fluid model matching the Mattis-Bardeen σ1 and λ = 90 nm at 10 GHz.
const two_fluid = Dict("Conductivity" => 8.80e6, "LondonDepth" => 0.090)

"""
    run_palace(config, out; palace_exec, np)

Write `config` with output directory `out` and run Palace on it.
"""
function run_palace(config, out; palace_exec, np)
    config["Problem"]["Output"] = out
    mkpath(joinpath(@__DIR__, dirname(out)))
    open(joinpath(@__DIR__, out * ".json"), "w") do io
        return JSON.print(io, config, 2)
    end
    cmd = `$palace_exec -np $np $(out * ".json")`
    return cd(() -> run(pipeline(cmd; stdout=devnull, stderr=devnull)), @__DIR__)
end

"""
    load_config(file)

Parse a configuration file of this example, which may contain comments.
"""
load_config(file) =
    JSON.parse(replace(read(joinpath(@__DIR__, file), String), r"//[^\n]*" => ""))

"""
    generate_superconducting_line_data(; kwargs...)

Run the boundary mode sweep (and the 3D wave port simulation with `waveport = true`), and
plot the attenuation, phase velocity, and pulses.
"""
function generate_superconducting_line_data(;
    palace_exec="palace",
    np=1,
    freqs=exp10.(range(0.0, 3.0, length=31)),  # GHz
    waveport=false
)
    for (kwargs, file) in [
        ((;), "superconducting_line.msh"),
        ((; surface=true), "superconducting_line_surface.msh"),
        (
            (; surface=true, l=10.0, lc_film=0.075, n_layers=4),
            "superconducting_line_l10.msh"
        ),
        (
            (; surface=true, l=20.0, lc_film=0.075, n_layers=8),
            "superconducting_line_l20.msh"
        )
    ]
        isfile(joinpath(@__DIR__, "mesh", file)) ||
            generate_superconducting_line_mesh(; filename=file, kwargs...)
    end

    # Boundary mode sweep: the Mattis-Bardeen volume and surface models, and the two-fluid
    # volume model.
    volume = load_config("superconducting_line_mode.json")
    surface = load_config("superconducting_line_mode_surface.json")
    tf = deepcopy(volume)
    nb = tf["Domains"]["Materials"][1]
    delete!(nb, "ComplexConductivity")
    merge!(nb, two_fluid)
    data = DataFrame(f=Float64[], model=String[], kn_re=Float64[], kn_im=Float64[])
    for (model, config) in [("MB", volume), ("MB surface", surface), ("two-fluid", tf)]
        for (i, f) in enumerate(freqs)
            config["Solver"]["BoundaryMode"]["Freq"] = f
            out = "postpro/mode_sweep/$(replace(model, " " => "_"))_$i"
            run_palace(config, out; palace_exec=palace_exec, np=np)
            kn = CSV.read(joinpath(@__DIR__, out, "mode-kn.csv"), DataFrame)
            push!(data, (f, model, kn[1, 2], kn[1, 3]))
        end
    end
    CSV.write(joinpath(@__DIR__, "postpro", "mode_sweep.csv"), data)

    if waveport
        # kₙ of the 3D line between wave ports from the transmission coefficients of two
        # lengths l, which cancels the port contributions to S21.
        config = load_config("superconducting_line_waveport.json")
        S21 = map([10.0, 20.0]) do l
            config["Model"]["Mesh"] = "mesh/superconducting_line_l$(Int(l)).msh"
            out = "postpro/waveport_l$(Int(l))"
            run_palace(config, out; palace_exec=palace_exec, np=np)
            S = CSV.read(
                joinpath(@__DIR__, out, "port-S.csv"),
                DataFrame;
                stripwhitespace=true
            )
            rename!(strip, S)
            s21 = @. exp10(S[!, "|S[2][1]| (dB)"] / 20) *
               cis(deg2rad(S[!, "arg(S[2][1]) (deg.)"]))
            return S[!, 1], s21
        end
        kn = @. 1im * log(S21[2][2] / S21[1][2]) / 10.0e-6
        CSV.write(
            joinpath(@__DIR__, "postpro", "waveport_kn.csv"),
            DataFrame(f=S21[1][1], kn_re=real.(kn), kn_im=imag.(kn))
        )
    end

    return plot_superconducting_line(data)
end

"""
    pulse(data, model, L; σ, t)

Propagate the Gaussian pulse exp(-t²/2σ²) over a length L [m] of line, in retarded time
`t` [s] (the group delay at the lowest frequency is removed). kₙ(f) is interpolated
linearly in the effective index, with zero attenuation at f = 0.
"""
function pulse(data, model, L; σ, t)
    d = data[data.model .== model, :]
    f0 = d.f * 1e9
    n_eff = complex.(d.kn_re, d.kn_im) ./ (2π * f0 / c0)
    function interp(ff)
        ff <= f0[1] && return real(n_eff[1])
        i = min(searchsortedlast(f0, ff), length(f0) - 1)
        return n_eff[i] + (n_eff[i + 1] - n_eff[i]) * (ff - f0[i]) / (f0[i + 1] - f0[i])
    end
    f = range(0.0, maximum(f0), length=2001)
    n = interp.(f)
    ω = 2π * f
    V = @. σ * sqrt(2π) * exp(-(ω * σ)^2 / 2 - 1im * (n - real(n_eff[1])) * ω / c0 * L)
    dω = ω[2] - ω[1]
    return [real(sum(@. V * cis(ω * tt))) * dω / π for tt in t]
end

function plot_superconducting_line(data)
    models = ["MB", "MB surface", "two-fluid"]
    fig = Figure(size=(900, 380))
    ax1 = Axis(fig[1, 1], xscale=log10, xlabel="Frequency (GHz)", ylabel="α (dB/mm)")
    ax2 = Axis(fig[1, 2], xscale=log10, xlabel="Frequency (GHz)", ylabel="vₚ / c")
    for model in models
        d = data[data.model .== model, :]
        lines!(ax1, d.f, -20 / log(10) * d.kn_im * 1e-3, label=model)
        lines!(ax2, d.f, 2π * d.f * 1e9 / c0 ./ d.kn_re, label=model)
    end
    axislegend(ax1, position=:lt)
    save(joinpath(@__DIR__, "postpro", "superconducting_line-1.png"), fig)

    # Single-flux-quantum pulse: 1 mV, FWHM 1.88 ps.
    σ = 1.88e-12 / (2 * sqrt(2 * log(2)))
    t = range(-15e-12, 15e-12, length=601)
    fig = Figure(size=(900, 420))
    for (k, model) in enumerate(["MB", "two-fluid"])
        ax = Axis(fig[1, k], title=model, xlabel="Time (ps)", ylabel="Amplitude (mV)")
        for (j, L) in enumerate([0, 8, 16, 32, 64, 128])
            v = pulse(data, model, L * 1e-3; σ=σ, t=t)
            lines!(ax, t * 1e12, v .- 0.2 * (j - 1), label="$L mm")
        end
        k == 1 && axislegend(ax, position=:rt)
    end
    save(joinpath(@__DIR__, "postpro", "superconducting_line-2.png"), fig)
    return fig
end
