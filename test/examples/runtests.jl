# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Hwloc

include("testcase.jl")

if "PALACE_TEST" in keys(ENV)
    palace = String.(split(ENV["PALACE_TEST"], ' '))
else
    palace = "palace"
end

if "NUM_PROC_TEST" in keys(ENV)
    numprocs = parse(Int, ENV["NUM_PROC_TEST"])
else
    numprocs = num_physical_cores()
end

if "TEST_CASES" in keys(ENV)
    cases = String.(split(ENV["TEST_CASES"], ' '))
else
    cases = [
        # "spheres",
        # "rings",
        # "cavity/pec",
        # "cavity/impedance",
        # "coaxial/open",
        # "coaxial/matched",
        # "cpw/lumped_uniform",
        # "cpw/wave_uniform",
        # "cpw/lumped_adaptive",
        # "cpw/wave_adaptive",
        "cavity/amr_conformal"
        # "cavity/amr_nonconformal"
    ]
end

@info "Starting regression tests using `$palace` on $numprocs core$(numprocs > 1 ? "s" : "")"

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
        excluded_columns=["Maximum", "Minimum"]
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

if "cavity/pec" in cases
    @info "Testing cavity (PEC)..."
    @time testcase(
        "cavity",
        "cavity_pec.json",
        "pec";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum", "Mean", "Error (Bkwd.)", "Error (Abs.)"],
        skip_rowcount=true
    )
end

if "cavity/impedance" in cases
    @info "Testing cavity (impedance)..."
    @time testcase(
        "cavity",
        "cavity_impedance.json",
        "impedance";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum", "Mean", "Error (Bkwd.)", "Error (Abs.)"],
        skip_rowcount=true
    )
end

if "cavity/amr_conformal" in cases
    @info "Testing cavity (amr conformal)..."
    @time testcase(
        "cavity",
        "cavity_amr_conformal.json",
        "amr/conformal";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum", "Mean", "Error (Bkwd.)", "Error (Abs.)"],
        skip_rowcount=true
    )
end

if "cavity/amr_nonconformal" in cases
    @info "Testing cavity (amr nonconformal)..."
    @time testcase(
        "cavity",
        "cavity_amr_nonconformal.json",
        "amr/nonconformal";
        palace=palace,
        np=numprocs,
        rtol=reltol,
        atol=abstol,
        excluded_columns=["Maximum", "Minimum", "Mean", "Error (Bkwd.)", "Error (Abs.)"],
        skip_rowcount=true
    )
end

# Coarser test tolerances for driven simulations with ports
reltol = 2.0e-2
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
    @info "Testing CPW (lumped ports)"
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
    @info "Testing CPW (wave ports)"
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
    @info "Testing CPW (lumped ports, adaptive)"
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
    @info "Testing CPW (wave ports, adaptive)"
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
