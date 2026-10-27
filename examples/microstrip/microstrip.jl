# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using CSV, DataFrames, CairoMakie, JSON

"""
    plot_microstrip(; root=@__DIR__)

Plot the substrate permittivity and the uniform/adaptive S-parameters. Run
`microstrip_driven.json` and `microstrip_adaptive.json` before calling this function.
The generated plot is `postpro/microstrip.png`.
"""
function plot_microstrip(; root::AbstractString=@__DIR__)
    function read_response(name)
        data = CSV.read(joinpath(root, "postpro", name, "port-S.csv"), DataFrame)
        rename!(data, strip.(names(data)))
        return data
    end
    uniform = read_response("driven")
    adaptive = read_response("adaptive")
    material =
        JSON.parsefile(joinpath(root, "microstrip_adaptive.json"))["Domains"]["Materials"][2]["Permittivity"]
    term = only(material["Terms"])
    @assert term["Type"] == "DjordjevicSarkar"
    ε(f) =
        material["HighFrequency"] +
        term["Strength"] *
        log(complex(term["UpperFrequency"], f) / complex(term["LowerFrequency"], f))
    f = range(first(adaptive[!, "f (GHz)"]), last(adaptive[!, "f (GHz)"]); length=300)
    values = ε.(f)

    figure = Figure(size=(1000, 700), fontsize=18)
    permittivity = Axis(figure[1, 1], xlabel="Frequency (GHz)", ylabel="Re εᵣ")
    loss =
        Axis(figure[1, 2], xlabel="Frequency (GHz)", ylabel="Loss tangent −Im εᵣ / Re εᵣ")
    lines!(permittivity, f, real.(values); color=:steelblue, linewidth=2)
    lines!(loss, f, -imag.(values) ./ real.(values); color=:darkorange, linewidth=2)
    for (column, title, position) in (
        ("|S[1][1]| (dB)", "Reflection |S₁₁|", (2, 1)),
        ("|S[2][1]| (dB)", "Transmission |S₂₁|", (2, 2))
    )
        axis = Axis(
            figure[position...],
            xlabel="Frequency (GHz)",
            ylabel="Magnitude (dB)",
            title=title
        )
        lines!(
            axis,
            adaptive[!, "f (GHz)"],
            adaptive[!, column];
            color=:steelblue,
            linewidth=2,
            label="Adaptive"
        )
        scatter!(
            axis,
            uniform[!, "f (GHz)"],
            uniform[!, column];
            color=:darkorange,
            markersize=10,
            label="Uniform"
        )
        axislegend(axis; position=:rb)
    end
    output = joinpath(root, "postpro", "microstrip.png")
    save(output, figure)
    return figure
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_microstrip()
end
