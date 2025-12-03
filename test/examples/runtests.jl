# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Hwloc
using Test

include("testcase.jl")

if "PALACE_TEST" in keys(ENV)
    palace = String.(split(ENV["PALACE_TEST"], ' '))
else
    palace = ["palace"]
end

if isnothing(Base.Sys.which(first(palace)))
    error(
        "Executable `$palace` not found. " *
        "You can customize the path by setting the PALACE_TEST environment variable"
    )
end

if "NUM_PROC_TEST" in keys(ENV)
    numprocs = parse(Int, ENV["NUM_PROC_TEST"])
else
    numprocs = num_physical_cores()
end

if "OMP_NUM_THREADS" in keys(ENV)
    numthreads = parse(Int, ENV["OMP_NUM_THREADS"])
else
    numthreads = 1
end

if "TEST_CASES" in keys(ENV)
    cases = String.(split(ENV["TEST_CASES"], ' '))
else
    cases = [
        "spheres",
        "rings",
        "antenna",
        "cylinder/cavity_pec",
        "cylinder/cavity_impedance",
        "cylinder/waveguide",
        "cylinder/floquet",
        "cylinder/driven_wave",
        "coaxial/open",
        "coaxial/matched",
        "cpw/lumped_uniform",
        "cpw/wave_uniform",
        "cpw/lumped_adaptive",
        "cpw/wave_adaptive",
        "cpw/lumped_eigen",
        "cpw/wave_eigen",
        "adapter/hybrid"
    ]
end

@info "Starting regression tests using `$palace` on $numprocs core$(numprocs > 1 ? "s" : "") $(numthreads > 1 ? "with $(numthreads) threads" : "")"

reltol = 1.0e-4
abstol = 1.0e-16

if "spheres" in cases
    @info "Testing spheres..."
    @time testcase(
        "spheres",
        "spheres.json",
        "";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum"],
        gridfunction_fields=true
    )
end

if "rings" in cases
    @info "Testing rings..."
    @time testcase(
        "rings",
        "rings.json",
        "";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum"]
    )
end

if "cylinder/cavity_pec" in cases
    @info "Testing cylinder/cavity (PEC)..."
    @time testcase(
        "cylinder",
        "cavity_pec.json",
        "cavity_pec";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum", "Mean", "Error (Bkwd.)", "Error (Abs.)"],
        skip_rowcount=true
    )
end

if "cylinder/cavity_impedance" in cases
    @info "Testing cylinder/cavity (impedance)..."
    @time testcase(
        "cylinder",
        "cavity_impedance.json",
        "cavity_impedance";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum", "Mean", "Error (Bkwd.)", "Error (Abs.)"],
        skip_rowcount=true
    )
end

if "cylinder/waveguide" in cases
    @info "Testing cylinder/waveguide (periodic)..."
    @time testcase(
        "cylinder",
        "waveguide.json",
        "waveguide";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum", "Mean", "Error (Bkwd.)", "Error (Abs.)"],
        skip_rowcount=true
    )
end

if "cylinder/floquet" in cases

    # The phase of the eigenmodes is not stable across architectures, so we
    # compare the magnitudes.
    function test_probe_magnitude(new_data, ref_data)
        num_probes = Int((size(new_data)[2] - 1) / 6) # 6 components per probe, 1 column with the mode number
        for i = 0:(num_probes - 1)
            # Compute magnitudes.
            Ex_new = new_data[:, 2 + i * 6] + 1im * new_data[:, 3 + i * 6]
            Ey_new = new_data[:, 4 + i * 6] + 1im * new_data[:, 5 + i * 6]
            Ez_new = new_data[:, 6 + i * 6] + 1im * new_data[:, 7 + i * 6]
            E_mag_new = sqrt.(abs.(Ex_new) .^ 2 + abs.(Ey_new) .^ 2 + abs.(Ez_new) .^ 2)

            Ex_ref = ref_data[:, 2 + i * 6] + 1im * ref_data[:, 3 + i * 6]
            Ey_ref = ref_data[:, 4 + i * 6] + 1im * ref_data[:, 5 + i * 6]
            Ez_ref = ref_data[:, 6 + i * 6] + 1im * ref_data[:, 7 + i * 6]
            E_mag_ref = sqrt.(abs.(Ex_ref) .^ 2 + abs.(Ey_ref) .^ 2 + abs.(Ez_ref) .^ 2)

            # Test magnitudes with relative tolerance.
            @test E_mag_new ≈ E_mag_ref rtol=reltol
        end
        return true
    end

    @info "Testing cylinder/floquet (periodic)..."
    @time testcase(
        "cylinder",
        "floquet.json",
        "floquet";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum", "Mean", "Error (Bkwd.)", "Error (Abs.)"],
        custom_tests=Dict(
            "probe-E.csv" => test_probe_magnitude,
            "probe-B.csv" => test_probe_magnitude
        ),
        skip_rowcount=true
    )
end

reltol = 1.0e-3
abstol = 1.0e-16

if "cylinder/driven_wave" in cases
    @info "Testing cylinder/driven_wave..."
    @time testcase(
        "cylinder",
        "driven_wave.json",
        "driven_wave";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum"]
    )
end

# Coarser test tolerances for driven simulations with ports
reltol = 2.0e-2
abstol = 1.0e-10

if "antenna" in cases
    @info "Testing antenna..."
    @time testcase(
        "antenna",
        "antenna.json",
        "";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=50abstol
    )
end

abstol = 2.0e-12

if "coaxial/open" in cases
    @info "Testing coaxial (open)..."
    @time testcase(
        "coaxial",
        "coaxial_open.json",
        "open";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum"]
    )
end

if "coaxial/matched" in cases
    @info "Testing coaxial (matched)..."
    @time testcase(
        "coaxial",
        "coaxial_matched.json",
        "matched";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum"]
    )
end

if "cpw/lumped_uniform" in cases
    @info "Testing CPW (lumped ports)..."
    @time testcase(
        "cpw",
        "cpw_lumped_uniform.json",
        "lumped_uniform";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum"]
    )
end

if "cpw/wave_uniform" in cases
    @info "Testing CPW (wave ports)..."
    @time testcase(
        "cpw",
        "cpw_wave_uniform.json",
        "wave_uniform";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum"]
    )
end

# Don't check accuracy for adaptive frequency sweep simulations

if "cpw/lumped_adaptive" in cases
    @info "Testing CPW (lumped ports, adaptive)..."
    @time testcase(
        "cpw",
        "cpw_lumped_adaptive.json",
        "lumped_adaptive";
        palace=palace,
        np=numprocs,
        rtol=Inf,
        atol=Inf
    )
end

if "cpw/wave_adaptive" in cases
    @info "Testing CPW (wave ports, adaptive)..."
    @time testcase(
        "cpw",
        "cpw_wave_adaptive.json",
        "wave_adaptive";
        palace=palace,
        np=numprocs,
        rtol=Inf,
        atol=Inf
    )
end

if "cpw/lumped_eigen" in cases

    # The phase of the eigenmodes is not stable across architectures, so we
    # compare the magnitudes.
    function test_farfield(new_data, ref_data)
        # Compute E field magnitudes.
        Ex_new = new_data[:, 4] + 1im * new_data[:, 5]
        Ey_new = new_data[:, 6] + 1im * new_data[:, 7]
        Ez_new = new_data[:, 8] + 1im * new_data[:, 9]
        E_mag_new = sqrt.(abs.(Ex_new) .^ 2 + abs.(Ey_new) .^ 2 + abs.(Ez_new) .^ 2)

        Ex_ref = ref_data[:, 4] + 1im * ref_data[:, 5]
        Ey_ref = ref_data[:, 6] + 1im * ref_data[:, 7]
        Ez_ref = ref_data[:, 8] + 1im * ref_data[:, 9]
        E_mag_ref = sqrt.(abs.(Ex_ref) .^ 2 + abs.(Ey_ref) .^ 2 + abs.(Ez_ref) .^ 2)

        # Test magnitudes with relative tolerance.
        @test E_mag_new ≈ E_mag_ref rtol=reltol
        return true
    end

    @info "Testing CPW (lumped ports, eigenmode)..."
    @time testcase(
        "cpw",
        "cpw_lumped_eigen.json",
        "lumped_eigen";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=[
            "Maximum",
            "Minimum",
            "Mean",
            "Error (Bkwd.)",
            "Error (Abs.)",
            "Re{V[1]} (V)",
            "Im{V[1]} (V)",
            "Re{V[2]} (V)",
            "Im{V[2]} (V)",
            "Re{V[3]} (V)",
            "Im{V[3]} (V)",
            "Re{V[4]} (V)",
            "Im{V[4]} (V)",
            "Re{I[1]} (A)",
            "Im{I[1]} (A)",
            "Re{I[2]} (A)",
            "Im{I[2]} (A)",
            "Re{I[3]} (A)",
            "Im{I[3]} (A)",
            "Re{I[4]} (A)",
            "Im{I[4]} (A)",
            "Q_ext[1]",
            "κ_ext[1] (GHz)",
            "Q_ext[2]",
            "κ_ext[2] (GHz)",
            "Q_ext[3]",
            "κ_ext[3] (GHz)",
            "Q_ext[4]",
            "κ_ext[4] (GHz)"
        ],
        custom_tests=Dict("farfield-rE.csv" => test_farfield),
        skip_rowcount=true
    )
end

if "cpw/wave_eigen" in cases
    @info "Testing CPW (wave ports, eigenmode)..."
    @time testcase(
        "cpw",
        "cpw_wave_eigen.json",
        "wave_eigen";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum", "Mean", "Error (Bkwd.)", "Error (Abs.)"],
        skip_rowcount=true
    )
end

if "adapter/hybrid" in cases
    @info "Testing adapter (hybrid)"
    @time testcase(
        "adapter",
        "hybrid.json",
        "hybrid";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum", "Mean", "Error (Bkwd.)", "Error (Abs.)"],
        skip_rowcount=true
    )
end

if "adapter/slp" in cases
    @info "Testing adapter (slp)"
    @time testcase(
        "adapter",
        "slp.json",
        "slp";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum", "Mean", "Error (Bkwd.)", "Error (Abs.)"],
        skip_rowcount=true
    )
end
