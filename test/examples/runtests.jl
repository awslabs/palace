# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Hwloc
using Test

include("argconfig.jl")
include("testcase.jl")

# Helper function to test farfield data by comparing E field magnitudes
# (phase is not stable across architectures). Column layout differs between driven
# (idx, excitation, theta, phi, rE0_re, rE0_im, ...) and eigenmode
# (idx, f_re, f_im, excitation, theta, phi, rE0_re, rE0_im, ...), so locate the
# real-part of rE0 by counting from the right.
function test_farfield(new_data, ref_data)
    n_components = 3  # Ex, Ey, Ez each contributing (re, im)
    n_complex_cols = 2 * n_components
    rE0_re_col = size(new_data, 2) - n_complex_cols + 1

    Ex_new = new_data[:, rE0_re_col + 0] + 1im * new_data[:, rE0_re_col + 1]
    Ey_new = new_data[:, rE0_re_col + 2] + 1im * new_data[:, rE0_re_col + 3]
    Ez_new = new_data[:, rE0_re_col + 4] + 1im * new_data[:, rE0_re_col + 5]
    E_mag_new = sqrt.(abs.(Ex_new) .^ 2 + abs.(Ey_new) .^ 2 + abs.(Ez_new) .^ 2)

    Ex_ref = ref_data[:, rE0_re_col + 0] + 1im * ref_data[:, rE0_re_col + 1]
    Ey_ref = ref_data[:, rE0_re_col + 2] + 1im * ref_data[:, rE0_re_col + 3]
    Ez_ref = ref_data[:, rE0_re_col + 4] + 1im * ref_data[:, rE0_re_col + 5]
    E_mag_ref = sqrt.(abs.(Ex_ref) .^ 2 + abs.(Ey_ref) .^ 2 + abs.(Ez_ref) .^ 2)

    # Test magnitudes with relative tolerance.
    @test E_mag_new ≈ E_mag_ref rtol=reltol
    return true
end

# Configure arguments.
arg_configs = [
    ArgConfig(
        name="palace-test",
        env_var="PALACE_TEST",
        default=["palace"],
        description="Palace executable path",
        parser=s -> String.(split(s, ' '))
    ),
    ArgConfig(
        name="num-proc-test",
        env_var="NUM_PROC_TEST",
        default=num_physical_cores(),
        description="Number of processes for testing",
        parser=s -> parse(Int, s)
    ),
    ArgConfig(
        name="palace-solver",
        env_var="PALACE_SOLVER",
        default="Default",
        description="Solver to test"
    ),
    ArgConfig(
        name="palace-eigensolver",
        env_var="PALACE_EIGENSOLVER",
        default="Default",
        description="Eigensolver to test"
    ),
    ArgConfig(
        name="omp-num-threads",
        env_var="OMP_NUM_THREADS",
        default=1,
        description="Number of OpenMP threads",
        parser=s -> parse(Int, s)
    ),
    ArgConfig(
        name="palace-device",
        env_var="PALACE_DEVICE",
        default="CPU",
        description="Device to use for testing (CPU or GPU)"
    ),
    ArgConfig(
        name="test-cases",
        env_var="TEST_CASES",
        default=[
            "spheres",
            "rings",
            # "transmon/transmon_coarse", # not included by default because of its cost
            # "transmon/transmon_amr", # not included by default because of its cost
            "antenna/antenna_halfwave_dipole",
            "antenna/antenna_short_dipole",
            "antenna/antenna_halfwave_dipole_surfacecurrent",
            "cylinder/cavity_pec",
            "cylinder/cavity_impedance",
            "cylinder/waveguide",
            "cylinder/floquet",
            "cylinder/driven_wave",
            "dielectric_grating/uniform",
            "coaxial/open",
            "coaxial/matched",
            "coaxial/lumped_wave",
            "cpw/lumped_uniform",
            "cpw/wave_uniform",
            "cpw/lumped_adaptive",
            "cpw/wave_adaptive",
            "cpw/lumped_eigen",
            "cpw/wave_eigen",
            "adapter/hybrid",
            "cavity2d/eigenmode",
            "cavity2d/driven",
            "cavity2d/electrostatic",
            "cavity2d/magnetostatic",
            "cavity2d/transient",
            "cpw2d/thin",
            "cpw2d/thick_impedance",
            "cpw/wave_2dmode"
        ],
        description="Test cases to run",
        parser=s -> String.(split(s, ' '))
    )
]

# Parse arguments
args = parse_args(arg_configs)

palace = args["palace-test"]
numprocs = args["num-proc-test"]
solver = args["palace-solver"]
eigensolver = args["palace-eigensolver"]
numthreads = args["omp-num-threads"]
device = args["palace-device"]
cases = args["test-cases"]

@info "Starting regression tests using `$palace` on $numprocs core$(numprocs > 1 ? "s" : "") $(numthreads > 1 ? "with $(numthreads) threads" : "")"
@info "Running test cases: $(join(cases, ", "))"

start_time = time()

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
        gridfunction_fields=true,
        device=device,
        linear_solver=solver,
        eigen_solver=eigensolver
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
        excluded_columns=["Maximum", "Minimum"],
        device=device,
        linear_solver="Default",
        eigen_solver=eigensolver
    )
end

reltol = 1.0e-2

if "transmon/transmon_coarse" in cases
    @info "Testing single transmon coarse..."
    @time testcase(
        "transmon",
        "transmon_coarse.json",
        "transmon_coarse";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        abs_columns=["κ_ext"],
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
            "Re{I[1]} (A)",
            "Im{I[1]} (A)",
            "Re{I[2]} (A)",
            "Im{I[2]} (A)",
            "Re{I[3]} (A)",
            "Im{I[3]} (A)"
        ],
        skip_rowcount=true,
        gridfunction_fields=true
    )
end

if "transmon/transmon_amr" in cases
    @info "Testing single transmon amr..."
    @time testcase(
        "transmon",
        "transmon_amr.json",
        "transmon_amr";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        abs_columns=["κ_ext"],
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
            "Re{I[1]} (A)",
            "Im{I[1]} (A)",
            "Re{I[2]} (A)",
            "Im{I[2]} (A)",
            "Re{I[3]} (A)",
            "Im{I[3]} (A)"
        ],
        skip_rowcount=true,
        gridfunction_fields=true
    )
end

reltol = 1.0e-4

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
        skip_rowcount=true,
        device=device,
        linear_solver="Default",
        eigen_solver=eigensolver
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
        skip_rowcount=true,
        device=device,
        linear_solver=solver,
        eigen_solver=eigensolver
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
        skip_rowcount=true,
        device=device,
        linear_solver=solver,
        eigen_solver=eigensolver
    )
end

if "cylinder/floquet" in cases
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
        skip_rowcount=true,
        device=device,
        linear_solver=solver,
        eigen_solver=eigensolver
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
        excluded_columns=["Maximum", "Minimum", "Mean"],
        device=device,
        linear_solver=solver,
        eigen_solver=eigensolver
    )
end

# Floquet port S-parameter test: compare only magnitude columns (|S| in dB),
# skipping NaN entries (evanescent modes) and negligible signals (< -200 dB).
function test_floquet_sparams(new_data, ref_data)
    for col_name in names(new_data)
        # Only compare magnitude columns, skip phase columns.
        occursin("|S[", col_name) && occursin("(dB)", col_name) || continue
        for (v_new, v_ref) in zip(new_data[!, col_name], ref_data[!, col_name])
            if isnan(v_new) && isnan(v_ref)
                @test true
            elseif v_ref < -200  # negligible signal, skip
                @test true
            else
                @test v_new ≈ v_ref rtol=reltol atol=abstol
            end
        end
    end
    return true
end

if "dielectric_grating/uniform" in cases
    @info "Testing dielectric_grating/uniform..."
    @time testcase(
        "dielectric_grating",
        "dielectric_grating_uniform.json",
        "uniform";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum"],
        custom_tests=Dict("port-floquet-S.csv" => test_floquet_sparams),
        paraview_fields=false,
        skip_rowcount=true,
        device=device,
        linear_solver=solver,
        eigen_solver=eigensolver
    )
end

# Coarser test tolerances for driven simulations with ports
reltol = 2.0e-2
abstol = 1.0e-10

if "antenna/antenna_halfwave_dipole" in cases
    @info "Testing antenna/antenna_halfwave_dipole..."
    @time testcase(
        "antenna",
        "antenna_halfwave_dipole.json",
        "antenna_halfwave_dipole";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=50abstol,
        device=device,
        linear_solver=solver
    )
end

if "antenna/antenna_short_dipole" in cases
    @info "Testing antenna/antenna_short_dipole..."
    @time testcase(
        "antenna",
        "antenna_short_dipole.json",
        "antenna_short_dipole";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=50abstol,
        custom_tests=Dict("farfield-rE.csv" => test_farfield),
        device=device,
        linear_solver=solver
    )
end

if "antenna/antenna_halfwave_dipole_surfacecurrent" in cases
    @info "Testing antenna/antenna_halfwave_dipole_surfacecurrent..."
    @time testcase(
        "antenna",
        "antenna_halfwave_dipole_surfacecurrent.json",
        "antenna_halfwave_dipole_surfacecurrent";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=50abstol,
        device=device,
        linear_solver=solver,
        paraview_fields=false
    )
end

abstol = 1.0e-11

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
        excluded_columns=["Maximum", "Minimum"],
        device=device,
        linear_solver=solver,
        eigen_solver=eigensolver
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
        excluded_columns=["Maximum", "Minimum"],
        device=device,
        linear_solver=solver,
        eigen_solver=eigensolver
    )
end

if "coaxial/lumped_wave" in cases
    @info "Testing coaxial (lumped + wave port mix)..."
    @time testcase(
        "coaxial",
        "coaxial_lumped_wave.json",
        "lumped_wave";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum"],
        device=device,
        linear_solver=solver,
        eigen_solver=eigensolver
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
        excluded_columns=["Maximum", "Minimum"],
        custom_tests=Dict("farfield-rE.csv" => test_farfield),
        device=device,
        linear_solver=solver,
        eigen_solver=eigensolver
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
        excluded_columns=["Maximum", "Minimum"],
        device=device,
        linear_solver=solver,
        eigen_solver=eigensolver
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
        atol=Inf,
        linear_solver=solver,
        eigen_solver=eigensolver
    )
end

# We skip when on GPUs because of this:
# https://github.com/awslabs/palace/issues/375
if "cpw/wave_adaptive" in cases && device != "GPU"
    @info "Testing CPW (wave ports, adaptive)..."
    @time testcase(
        "cpw",
        "cpw_wave_adaptive.json",
        "wave_adaptive";
        palace=palace,
        np=numprocs,
        rtol=Inf,
        atol=Inf,
        linear_solver=solver,
        eigen_solver=eigensolver
    )
end

if "cpw/lumped_eigen" in cases
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
        skip_rowcount=true,
        device=device,
        linear_solver=solver,
        eigen_solver=eigensolver
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
        skip_rowcount=true,
        device=device,
        linear_solver=solver,
        eigen_solver=eigensolver
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
        skip_rowcount=true,
        device=device,
        linear_solver=solver,
        eigen_solver=eigensolver
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
        skip_rowcount=true,
        device=device,
        linear_solver=solver,
        eigen_solver=eigensolver
    )
end

reltol = 1.0e-4
abstol = 1.0e-16

if "cavity2d/eigenmode" in cases
    @info "Testing cavity2d/eigenmode (2D eigenmode)..."
    @time testcase(
        "cavity2d",
        "cavity2d.json",
        "eigenmode";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum", "Mean", "Error (Bkwd.)", "Error (Abs.)"],
        skip_rowcount=true,
        device=device,
        linear_solver="Default",
        eigen_solver=eigensolver
    )
end

# Coarser test tolerances for 2D driven simulations: the coarse meshes and low-order ports
# make results more sensitive to MPI partitioning and platform differences.
reltol = 2.0e-2
abstol = 1.0e-8

if "cavity2d/driven" in cases
    @info "Testing cavity2d/driven (2D driven)..."
    @time testcase(
        "cavity2d",
        "cavity2d_driven.json",
        "driven";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum"],
        device=device,
        linear_solver="Default",
        eigen_solver=eigensolver
    )
end

reltol = 1.0e-4
abstol = 1.0e-10

if "cavity2d/electrostatic" in cases
    @info "Testing cavity2d/electrostatic (2D electrostatic)..."
    @time testcase(
        "cavity2d",
        "cavity2d_electrostatic.json",
        "electrostatic";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum"],
        device=device,
        linear_solver=solver,
        eigen_solver=eigensolver
    )
end

if "cavity2d/magnetostatic" in cases
    @info "Testing cavity2d/magnetostatic (2D magnetostatic)..."
    @time testcase(
        "cavity2d",
        "cavity2d_magnetostatic.json",
        "magnetostatic";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum"],
        device=device,
        linear_solver="Default",
        eigen_solver=eigensolver
    )
end

if "cavity2d/transient" in cases
    @info "Testing cavity2d/transient (2D transient)..."
    @time testcase(
        "cavity2d",
        "cavity2d_transient.json",
        "transient";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum"],
        device=device,
        linear_solver="Default",
        eigen_solver=eigensolver
    )
end

# 2D boundary mode: impedance depends on voltage path integration, tolerant for
# cross-platform reproducibility.
reltol = 1.0e-2

if "cpw2d/thin" in cases
    @info "Testing cpw2d/thin (2D mode analysis, thin PEC)..."
    @time testcase(
        "cpw2d",
        "cpw2d_thin.json",
        "thin";
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
            "Im{kn} (1/m)",
            "Im{n_eff}"
        ],
        custom_tests=Dict(
            "mode-V.csv" =>
                (data, dataref) -> begin
                    # Compare voltage magnitudes (phase is arbitrary for eigenmodes).
                    re_cols = filter(n -> startswith(n, "Re{V["), names(data))
                    for rc in re_cols
                        idx = match(r"Re\{V\[(\d+)\]\}", rc)[1]
                        ic = "Im{V[$idx]} (V)"
                        mag = sqrt.(data[!, rc] .^ 2 .+ data[!, ic] .^ 2)
                        mag_ref = sqrt.(dataref[!, rc] .^ 2 .+ dataref[!, ic] .^ 2)
                        for (i, (v, vr)) in enumerate(zip(mag, mag_ref))
                            @test isapprox(v, vr; rtol=reltol, atol=abstol) ||
                                  (@warn "|V[$idx]| row $i: $vr ≉ $v"; false)
                        end
                    end
                end
        ),
        skip_rowcount=true,
        device=device,
        linear_solver="Default",
        eigen_solver=eigensolver
    )
end

if "cpw2d/thick_impedance" in cases
    @info "Testing cpw2d/thick_impedance (2D mode analysis, thick impedance)..."
    @time testcase(
        "cpw2d",
        "cpw2d_thick_impedance.json",
        "thick_impedance";
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
            "Im{kn} (1/m)",
            "Im{n_eff}"
        ],
        skip_rowcount=true,
        device=device,
        linear_solver="Default",
        eigen_solver=eigensolver
    )
end

reltol = 1.0e-4

if "cpw/wave_2dmode" in cases
    @info "Testing cpw/wave_2dmode (2D mode analysis from 3D mesh)..."
    @time testcase(
        "cpw",
        "cpw_wave_2dmode.json",
        "wave_2dmode";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum", "Mean", "Error (Bkwd.)", "Error (Abs.)"],
        skip_rowcount=true,
        device=device,
        linear_solver="Default",
        eigen_solver=eigensolver
    )
end

total_time = time() - start_time
@info "Total runtime: $(round(total_time, digits=2)) seconds"
