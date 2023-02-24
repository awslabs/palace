# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Hwloc
include("testcase.jl")

if "NUM_PROC_TEST" in keys(ENV)
    numprocs = ENV["NUM_PROC_TEST"]
else
    numprocs = num_physical_cores()
end

reltol = 1.0e-4
abstol = 1.0e-16

@info "Starting regression tests on $numprocs cores"

@info "Testing spheres..."
@time testcase("spheres", "spheres.json", ""; np=numprocs, rtol=reltol, atol=abstol)

@info "Testing rings..."
@time testcase("rings", "rings.json", ""; np=numprocs, rtol=reltol, atol=abstol)

@info "Testing cavity (PEC)..."
@time testcase("cavity", "cavity_pec.json", "pec"; np=numprocs, rtol=reltol, atol=abstol)

@info "Testing cavity (impedance)..."
@time testcase(
    "cavity",
    "cavity_impedance.json",
    "impedance";
    np=numprocs,
    rtol=reltol,
    atol=abstol
)

# Coarser test tolerances for driven simulations with ports
reltol = 2.0e-2
abstol = 2.0e-12

@info "Testing coaxial (open)..."
@time testcase(
    "coaxial",
    "coaxial_open.json",
    "open";
    np=numprocs,
    rtol=reltol,
    atol=abstol
)

@info "Testing coaxial (matched)..."
@time testcase(
    "coaxial",
    "coaxial_matched.json",
    "matched";
    np=numprocs,
    rtol=reltol,
    atol=abstol
)

@info "Testing CPW (lumped ports)"
@time testcase(
    "cpw",
    "cpw_lumped_uniform.json",
    "lumped_uniform";
    np=numprocs,
    rtol=reltol,
    atol=abstol
)

@info "Testing CPW (wave ports)"
@time testcase(
    "cpw",
    "cpw_wave_uniform.json",
    "wave_uniform";
    np=numprocs,
    rtol=reltol,
    atol=abstol
)

# Don't check accuracy for adaptive frequency sweep simulations

@info "Testing CPW (lumped ports, adaptive)"
@time testcase(
    "cpw",
    "cpw_lumped_adaptive.json",
    "lumped_adaptive";
    np=numprocs,
    rtol=Inf,
    atol=Inf
)

@info "Testing CPW (wave ports, adaptive)"
@time testcase(
    "cpw",
    "cpw_wave_adaptive.json",
    "wave_adaptive";
    np=numprocs,
    rtol=Inf,
    atol=Inf
)
