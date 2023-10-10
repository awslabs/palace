# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Test
using Hwloc
using CSV
using DataFrames

include("testcase.jl")

if "PALACE_TEST" in keys(ENV)
    palace = split(ENV["PALACE_TEST"], ' ')
else
    palace = "palace"
end

if "NUM_PROC_TEST" in keys(ENV)
    numprocs = ENV["NUM_PROC_TEST"]
else
    numprocs = num_physical_cores()
end

reltol = 1.0e-4
abstol = 1.0e-16

@info "Starting regression tests using `$palace` on $numprocs cores"

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

@info "Testing cavity (PEC)..."
@time testcase(
    "cavity",
    "cavity_pec.json",
    "pec";
    palace=palace,
    np=numprocs,
    rtol=reltol,
    atol=abstol,
    excluded_columns=["Maximum", "Minimum"],
    skip_rowcount=true
)

@info "Testing cavity (impedance)..."
@time testcase(
    "cavity",
    "cavity_impedance.json",
    "impedance";
    palace=palace,
    np=numprocs,
    rtol=reltol,
    atol=abstol,
    excluded_columns=["Maximum", "Minimum"],
    skip_rowcount=true
)

# Coarser test tolerances for driven simulations with ports
reltol = 2.0e-2
abstol = 2.0e-12

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

# Don't check accuracy for adaptive frequency sweep simulations

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
