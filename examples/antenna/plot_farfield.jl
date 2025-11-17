# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""
Generate polar plots on the E and B planes and 3D relative radiation pattern
from the Palace farfield output data (`farfield-rE.csv`). Only plot the first
frequency found in the file.

The 3D relative radiation pattern plot depicts a 3D surface whose distance from
the center is proportional to the strength of the electric field.

Usage:
    julia plot_farfield.jl [model=<model_type> file=<filename>]

Arguments:
    model=<model_type>  - Dipole model type: "antenna_short_dipole" or "antenna_halfwave_dipole"
    file=<filename>     - Path to farfield file

Note: If any arguments are provided, both model= and file= must be specified.

Defaults:
    model: antenna_halfwave_dipole
    file: postpro/antenna_halfwave_dipole/farfield-rE.csv

Requires CSV, DataFrames, and CairoMakie.
"""

import CSV, CairoMakie
import DataFrames: DataFrame, rename!

# Tip: you can choose a different Makie backend. For example GLMakie gives you
# interactivity. CairoMakie is best at producing static images, especially on
# headless servers.
Makie = CairoMakie.Makie

"""
    compute_field_magnitude(df)

Compute total electric field magnitude from real and imaginary components.

Returns |E| = sqrt(|Ex|² + |Ey|² + |Ez|²) where |Ec|² = Re{Ec}² + Im{Ec}² as a
Vector.
"""
function compute_field_magnitude(df)
    Ex_re = df[:, "r*Re{E_x}"]
    Ex_im = df[:, "r*Im{E_x}"]
    Ey_re = df[:, "r*Re{E_y}"]
    Ey_im = df[:, "r*Im{E_y}"]
    Ez_re = df[:, "r*Re{E_z}"]
    Ez_im = df[:, "r*Im{E_z}"]

    return sqrt.(
        Ex_re .^ 2 + Ex_im .^ 2 + Ey_re .^ 2 + Ey_im .^ 2 + Ez_re .^ 2 + Ez_im .^ 2
    )
end

"""
    extract_eplane(df, tolerance_deg=1)

Extract E-plane (xz-plane) radiation pattern data.

This is done by identifying points with φ = 0° and φ = 180° within the given
tolerance.

We have to be a little careful here because we need to offset the θs by 180°
to make sure we cover the entire circumference.
"""
function extract_eplane(df, tolerance_deg=1)
    # E-plane: phi = 0° (no offset) and phi = 180°
    angles_list = Float64[]
    magnitude_list = Float64[]

    # phi ≈ 0°, no theta offset.
    mask1 = abs.(df[:, "phi"] .- 0) .< tolerance_deg
    if any(mask1)
        data1 = df[mask1, :]
        append!(angles_list, data1[:, "theta"])
        append!(magnitude_list, compute_field_magnitude(data1))
    end

    # phi ≈ 180°, theta offset by 180°.
    mask2 = abs.(df[:, "phi"] .- 180) .< tolerance_deg
    if any(mask2)
        data2 = df[mask2, :]
        offset_angles = (data2[:, "theta"] .+ 180) .% 360
        append!(angles_list, offset_angles)
        append!(magnitude_list, compute_field_magnitude(data2))
    end

    # Finally, sort by angle so that we have a nice line.
    sort_idx = sortperm(angles_list)
    return angles_list[sort_idx], magnitude_list[sort_idx]
end

"""
    extract_hplane(df, tolerance_deg=1)

Extract H-plane (xy-plane) radiation pattern data.

This is done by identifying points with θ = 90° within the given tolerance.
"""
function extract_hplane(df, tolerance_deg=1)
    # H-plane: theta = 90°.
    mask = abs.(df[:, "theta"] .- 90) .< tolerance_deg

    hplane_data = df[mask, :]
    return hplane_data[:, "phi"], compute_field_magnitude(hplane_data)
end

"""
    compute_db(magnitude)

Convert field magnitude to normalized dB scale and return values normalized to 0
dB maximum.
"""
function compute_db(magnitude)
    # Floor values 10 orders of magnitude below maximum to avoid log(0).
    magnitude = max.(magnitude, maximum(magnitude) * 1e-10)
    db_values = 20 * log10.(magnitude)
    return db_values .- maximum(db_values)
end

"""
    generate_theoretical_dipole(model_type="antenna_halfwave_dipole")

Return angles and patterns for theoretical dipole radiation pattern.

Arguments:
    model_type - "antenna_short_dipole" for electrical current dipole or "antenna_halfwave_dipole" for half-wave dipole

For antenna_short_dipole (z-oriented electrical current dipole):
Based on geosci.xyz far-field analytical solution for z-oriented dipole.
- E-plane (xz): |E_θ|² ∝ |sin(θ)|² (nulls along dipole axis, maxima perpendicular)
- H-plane (xy): |H_φ|² ∝ constant (omnidirectional)

For antenna_halfwave_dipole:
- E-plane: [cos(π/2 * cos(θ)) / sin(θ)]²
- H-plane: omnidirectional (constant)
"""
function generate_theoretical_dipole(model_type="antenna_halfwave_dipole")
    angles = 0:360

    if model_type == "antenna_short_dipole"
        # E-plane: sin(θ)²
        eplane = zeros(length(angles))
        for (i, θ_deg) in enumerate(angles)
            θ_rad = deg2rad(θ_deg)
            sin_θ = sin(θ_rad)
            if abs(sin_θ) > 1e-6
                eplane[i] = abs(sin_θ)
            end
        end
    elseif model_type == "antenna_halfwave_dipole"
        # E-plane: [cos(π/2 * cos(θ)) / sin(θ)]²
        eplane = zeros(length(angles))
        for (i, θ_deg) in enumerate(angles)
            θ_rad = deg2rad(θ_deg)
            sin_θ = sin(θ_rad)
            if abs(sin_θ) > 1e-6
                eplane[i] = abs(cos(π/2 * cos(θ_rad)) / sin_θ)
            end
        end
    else
        error("Unknown model_type: $model_type. Use 'antenna_short_dipole' or 'antenna_halfwave_dipole'")
    end

    # H-plane: omnidirectional
    hplane = ones(length(angles))

    return angles, eplane, angles, hplane
end

"""
    polar_plots(data, label, model_type="halfwave_dipole", filename="farfield_polar.png")

Plot the polar radiation patterns and the expected dipole pattern.
"""
function polar_plots(data, label, model_type="halfwave_dipole", filename="farfield_polar.png")
    e_angles, e_mag = extract_eplane(data)
    h_angles, h_mag = extract_hplane(data)

    theo_angles, theo_eplane, _, theo_hplane = generate_theoretical_dipole(model_type)
    theo_e_db = compute_db(theo_eplane)
    theo_h_db = compute_db(theo_hplane)

    fig = Makie.Figure(size=(800, 400))

    if length(e_mag) > 1
        e_db = compute_db(e_mag)

        # E-plane plot
        ax1 = Makie.PolarAxis(
            fig[1, 1],
            title="E-plane ($label)",
            rticks=-25:5:2,
            radius_at_origin=-25,
            rlimits=(-25, 2),
            theta_0=-π / 2,
            direction=-1,
            rgridcolor=:lightgray,
            thetagridcolor=:lightgray
        )
        Makie.lines!(
            ax1,
            Makie.deg2rad.(theo_angles),
            theo_e_db,
            linewidth=1,
            linestyle=:dash,
            color=:black
        )
        Makie.lines!(ax1, Makie.deg2rad.(e_angles), e_db, linewidth=2, color=:blue)
    else
        println("No point found on the E-plane")
    end

    if length(h_mag) > 1
        h_db = compute_db(h_mag)

        # H-plane plot
        ax2 = Makie.PolarAxis(
            fig[1, 2],
            title="H-plane ($label)",
            rticks=-25:5:2,
            radius_at_origin=-25,
            rlimits=(-25, 2),
            theta_0=-π / 2,
            direction=-1,
            rgridcolor=:lightgray,
            thetagridcolor=:lightgray
        )
        Makie.lines!(
            ax2,
            Makie.deg2rad.(theo_angles),
            theo_h_db,
            linewidth=1,
            linestyle=:dash,
            color=:black
        )
        Makie.lines!(ax2, Makie.deg2rad.(h_angles), h_db, linewidth=2, color=:blue)
    else
        println("No point found on the H-plane")
    end

    Makie.save(filename, fig)
    return println("Saved $filename")
end

"""
    three_d_plot(data, label, filename="farfield_3d.png")

Plot a 3D representation of the strength of the electric field.

The plot represents the normalized magnitude of electric field with a mesh where
the radial distance is proportional to the magnitude itself.
"""
function three_d_plot(data, label, filename="farfield_3d.png")
    E_mag = compute_field_magnitude(data)
    E_mag ./= maximum(E_mag)

    theta_rad = deg2rad.(data[:, "theta"])
    phi_rad = deg2rad.(data[:, "phi"])

    x = E_mag .* sin.(theta_rad) .* cos.(phi_rad)
    y = E_mag .* sin.(theta_rad) .* sin.(phi_rad)
    z = E_mag .* cos.(theta_rad)

    fig3d = Makie.Figure(size=(600, 450))
    ax = Makie.Axis3(fig3d[1, 1], title="Relative E-field magnitude ($label)", aspect=:data)
    Makie.mesh!(ax, x, y, z, color=E_mag)
    Makie.save(filename, fig3d)
    return println("Saved $filename")
end

function main()
    # Set defaults
    model_type = "antenna_halfwave_dipole"
    filename = "postpro/antenna_halfwave_dipole/farfield-rE.csv"

    # Check if any arguments provided
    if length(ARGS) == 0
        # Use defaults - no arguments provided
    elseif length(ARGS) == 2
        # Both arguments must be provided
        has_model = false
        has_file = false

        # Parse named arguments
        for arg in ARGS
            if startswith(arg, "model=")
                model_type = split(arg, "=", limit=2)[2]
                has_model = true
            elseif startswith(arg, "file=")
                filename = split(arg, "=", limit=2)[2]
                has_file = true
            else
                println("Error: Invalid argument '$arg'")
                println("Usage: julia --project plot_farfield.jl [model=<model_type> file=<filename>]")
                println("Valid model types: 'antenna_short_dipole', 'antenna_halfwave_dipole'")
                println("Examples:")
                println("  julia plot_farfield.jl")
                println("  julia plot_farfield.jl model=antenna_short_dipole file=data.csv")
                return
            end
        end

        # Ensure both arguments were provided
        if !has_model || !has_file
            println("Error: If arguments are provided, both model= and file= must be specified")
            println("Usage: julia --project plot_farfield.jl [model=<model_type> file=<filename>]")
            println("Examples:")
            println("  julia plot_farfield.jl")
            println("  julia plot_farfield.jl model=antenna_short_dipole file=data.csv")
            return
        end
    else
        println("Error: If arguments are provided, both model= and file= must be specified")
        println("Usage: julia --project plot_farfield.jl [model=<model_type> file=<filename>]")
        println("Examples:")
        println("  julia plot_farfield.jl")
        println("  julia plot_farfield.jl model=antenna_short_dipole file=data.csv")
        return
    end

    if !(model_type in ["antenna_short_dipole", "antenna_halfwave_dipole"])
        println("Error: Invalid model_type '$model_type'")
        println("Valid options: 'antenna_short_dipole', 'antenna_halfwave_dipole'")
        println("Usage: julia --project plot_farfield.jl [model=<model_type> file=<filename>]")
        return
    end

    if !isfile(filename)
        println("Error: File '$filename' not found")
        println("Usage: julia --project plot_farfield.jl [model=<model_type> file=<filename>]")
        return
    end

    println("Reading farfield data from: $filename")
    println("Using dipole model: $model_type")

    # Read the entire data file into a DataFrame
    df = CSV.read(filename, DataFrame)
    # Remove spaces and units from column names.
    rename!(df, [name => strip(replace(name, r"\s*\([^)]*\)" => "")) for name in names(df)])

    # Process the first frequency/mode we find.
    if "m" in names(df)         # Eigenmode
        mode = df[:, "m"]
        println("Processing mode: $mode")
        label = "m = $mode"
        data = filter(row -> row["m"] == mode, df)
    else
        freq = first(df[:, "f"])
        println("Processing frequency: $(freq) GHz")
        label = "freq = $freq GHz"
        data = filter(row -> row["f"] == freq, df)
    end

    polar_filename = "farfield_polar_$(model_type).png"
    three_d_filename = "farfield_3d_$(model_type).png"

    polar_plots(data, label, model_type, polar_filename)
    three_d_plot(data, label, three_d_filename)
    return nothing
end

if !isinteractive()
    main()
end
