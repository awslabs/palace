# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Hwloc
using Test

include("argconfig.jl")
include("testcase.jl")

# Helper function to test farfield data by comparing E field magnitudes
# (phase is not stable across architectures)
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
            # "transmon/coarse", # not included by default because of its cost
            # "transmon/amr", # not included by default because of its cost
            "antenna/antenna_halfwave_dipole",
            "antenna/antenna_short_dipole",
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

if "transmon/coarse" in cases
    @info "Testing single transmon coarse..."
    @time testcase(
        "transmon",
        "transmon_coarse.json",
        "transmon_coarse";
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

if "transmon/amr" in cases
    @info "Testing single transmon arm..."
    @time testcase(
        "transmon",
        "transmon_amr.json",
        "transmon_amr";
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
        excluded_columns=["Maximum", "Minimum"],
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

total_time = time() - start_time
@info "Total runtime: $(round(total_time, digits=2)) seconds"
