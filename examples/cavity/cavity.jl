# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using CSV
using DataFrames
using Dates
using ForwardDiff
using JSON
using Roots
using SpecialFunctions

include(joinpath(@__DIR__, "mesh", "mesh.jl"))

"""
    solve_cavity_resonator(
        params;
        order,
        refinement,
        mesh_type::Integer      = 0,
        radius::Real            = 2.74,
        aspect_ratio::Real      = 1.0,
        num_processors::Integer = 1,
        cleanup_files::Bool     = true
    )

Solve the cavity mode problem, with an automatically generated Gmsh mesh

See also [`generate_cylindrical_cavity_mesh`](@ref)

# Arguments

  - params - dictionary storing the parsed json parameter file
  - order - the polynomial order used in the solution representation
  - geo_order - the polynomial order used in the mesh representation
  - refinement - the level of mesh refinement
  - mesh_type - 0 = tetrahedral mesh, 1 = prism mesh, 2 = hexahedral mesh
  - radius - the radius of the cavity resonator
  - aspect_ratio - the ratio of the DIAMETER of the cavity to the height
  - num_processors - number of processors to use for the simulation
  - cleanup_files - delete temporary mesh and configuration files after simulation
"""
function solve_cavity_resonator(
    params::Dict;
    order::Integer,
    geo_order::Integer,
    refinement::Integer,
    mesh_type::Integer      = 0,
    radius::Real            = 2.74,
    aspect_ratio::Real      = 1.0,
    num_processors::Integer = 1,
    cleanup_files::Bool     = true
)
    @assert refinement >= 0
    @assert order > 0
    @assert geo_order > 0
    @assert mesh_type ∈ [0, 1, 2]

    # Generate a mesh
    cavity_dir = @__DIR__
    file_root = string("cavity_p", order, "_h", refinement)
    mesh_filename = string(file_root, ".msh")
    generate_cylindrical_cavity_mesh(
        filename=mesh_filename,
        refinement=refinement,
        order=geo_order,
        mesh_type=mesh_type,
        radius=radius,
        aspect_ratio=aspect_ratio,
        verbose=0
    )

    # Generate solver parameter file
    params["Solver"]["Order"] = order
    params["Model"]["Mesh"] = joinpath(cavity_dir, "mesh", mesh_filename)
    json_filename = string(file_root, ".json")
    open(joinpath(cavity_dir, json_filename), "w") do f
        return JSON.print(f, params)
    end

    # Call the solver, storing the terminal output
    call_command = Cmd(`palace -np $num_processors $json_filename`, dir=cavity_dir)
    log_file = read(call_command, String)
    # println(log_file)

    # Search through for the DOF count
    start_ind = findfirst("ND", log_file)[end]
    start_ind += findfirst(":", log_file[start_ind:end])[end]
    end_ind = start_ind + findfirst(",", log_file[start_ind:end])[1]
    dof = parse(Int, filter(isdigit, log_file[start_ind:end_ind]))

    # Extract the top two frequency modes
    eig_df = CSV.read(joinpath(cavity_dir, "postpro", "convergence", "eig.csv"), DataFrame)
    eig = Matrix(eig_df[:, 2:end])[:, 1]

    # Clean up the parameter and mesh file
    if cleanup_files
        rm(joinpath(cavity_dir, "mesh", mesh_filename))
        rm(joinpath(cavity_dir, json_filename))
    end

    return dof, eig
end

∂besselj = (ν, x) -> ForwardDiff.derivative(y -> besselj(ν, y), x)

"""
    besselj_roots(ν, n::Integer)::Float64

Compute the n-th root of the Bessel functions of first kind, J_ν

Note: Compare against https://mathworld.wolfram.com/BesselFunctionZeros.html
with besselj_roots.((0:5)', 1:5)
"""
function besselj_roots(ν, n::Integer)::Float64
    upper_bound = 10
    roots = find_zeros(x -> besselj(ν, x), 0, upper_bound)
    while length(roots) < n + 1
        upper_bound *= 2
        roots = find_zeros(x -> besselj(ν, x), 0, upper_bound)
    end
    if roots[1] < 1e-10
        # If the first root is marked as 0, ignore
        popfirst!(roots)
    end
    return roots[n]
end

"""
    ∂besselj_roots(ν, n::Integer)::Float64

Compute the n-th root of the first derivative of the Bessel functions of first kind, J'_ν

Note: Compare against https://mathworld.wolfram.com/BesselFunctionZeros.html with
∂besselj_roots.((0:5)', 1:5)
"""
function ∂besselj_roots(ν, n::Integer)::Float64
    upper_bound = 10
    roots = find_zeros(x -> ∂besselj(ν, x), 0, upper_bound)
    while length(roots) < n + 1
        upper_bound *= 2
        roots = find_zeros(x -> ∂besselj(ν, x), 0, upper_bound)
    end
    if roots[1] < 1e-10
        # If the first root is marked as 0, ignore
        popfirst!(roots)
    end
    return roots[n]
end

"""
    frequency_transverse(n, m, l; ϵᵣ, μᵣ, a, d)

Compute the resonant frequency of the transverse electric and magnetic mode indexed by n, m,
l, in GHz

# Arguments

  - n - mode number in the circumferential direction
  - m - mode number in the radial direction
  - l - mode number in the z direction
  - ϵᵣ - relative electric permittivity
  - μᵣ - relative magnetic permeability
  - a - radius of cavity in centimeters
  - d - height of cavity in centimeters
"""
function frequency_transverse(n, m, l; ϵᵣ, μᵣ, a_cm, d_cm)
    ϵ₀ = 8.8541878176e-12
    μ₀ = 4e-7 * π

    a = a_cm * 0.01
    d = d_cm * 0.01

    c = 1.0 / sqrt(ϵ₀ * μ₀)
    p_nm = besselj_roots(n, m)
    ∂p_nm = ∂besselj_roots(n, m)

    C = (c / (2 * π * sqrt(ϵᵣ * μᵣ)))

    f_M = C * sqrt((p_nm / a)^2 + (l * π / d)^2) / 1e9
    f_E = C * sqrt((∂p_nm / a)^2 + (l * π / d)^2) / 1e9

    return f_E, f_M
end

"""
    generate_cavity_convergence_data(
        p_min::Integer          = 1,
        p_max::Integer          = 3,
        ref_min::Integer        = 0,
        ref_max::Integer        = 3,
        mesh_type::Integer      = 0,
        num_processors::Integer = 1
    )

Generate the data for the cavity convergence study

# Arguments

  - p_min - minimum polynomial order
  - p_max - maximum polynomial order
  - ref_min - minimum number of levels of uniform mesh refinement
  - ref_max - maximum number of levels of uniform mesh refinement
  - mesh_type - 0 = tetrahedral mesh, 1 = prism mesh, 2 = hexahedral mesh
  - num_processors - number of processors to use for the simulation
"""
function generate_cavity_convergence_data(;
    p_min::Integer          = 1,
    p_max::Integer          = 3,
    ref_min::Integer        = 0,
    ref_max::Integer        = 3,
    mesh_type::Integer      = 0,
    num_processors::Integer = 1
)
    # Load the default JSON script (the file contains comments and we need to sanitize them)
    cavity_dir = @__DIR__
    params = open(joinpath(cavity_dir, "cavity_pec.json"), "r") do f
        return JSON.parse(join(getindex.(split.(eachline(f), "//"), 1), "\n"))
    end

    # Update the dictionary
    params["Problem"]["Verbose"] = 2
    params["Problem"]["Output"] = joinpath(cavity_dir, "postpro", "convergence")
    params["Model"]["Refinement"]["UniformLevels"] = 0 # Don't perform any mesh refinement
    params["Solver"]["Eigenmode"]["Save"] = 0 # Don't write any fields to file
    params["Solver"]["Eigenmode"]["N"] = 4 # Look only for the top 4 modes
    params["Solver"]["Eigenmode"]["Tol"] = 1.0e-12
    params["Solver"]["Eigenmode"]["Target"] = 2.0
    params["Solver"]["Eigenmode"]["StartVectorConstant"] = true
    params["Solver"]["Linear"]["Tol"] = 1.0e-14

    # Compute the exact solution for reference
    radius = 2.74
    aspect_ratio = 1 / sqrt(2)
    ~, f_TM_010_true = frequency_transverse(
        0,
        1,
        0;
        ϵᵣ   = 2.08,
        μᵣ   = 1.0,
        a_cm = radius,
        d_cm = aspect_ratio * 2 * radius
    )
    f_TE_111_true, ~ = frequency_transverse(
        1,
        1,
        1;
        ϵᵣ   = 2.08,
        μᵣ   = 1.0,
        a_cm = radius,
        d_cm = aspect_ratio * 2 * radius
    )

    dof = Vector{Vector{Int}}()
    f_TM_010 = Vector{Vector{Float64}}()
    f_TE_111 = Vector{Vector{Float64}}()

    # Generate the data
    for p = p_min:p_max
        push!(dof, eltype(dof)())
        push!(f_TM_010, eltype(f_TM_010)())
        push!(f_TE_111, eltype(f_TE_111)())
        for ref = ref_min:ref_max
            print("p = ", p, ", ref = ", ref, ": ")
            results = solve_cavity_resonator(
                params,
                order=p,
                geo_order=p,
                refinement=ref,
                mesh_type=mesh_type,
                radius=radius,
                aspect_ratio=aspect_ratio,
                num_processors=min(num_processors, 4 * 2^(2 * ref))
            )
            println("Success! $(results[1]) dofs, finished at $(now())")
            push!(dof[end], results[1])

            # The first dominant frequency should be the magnetic, and the second the
            # electric, but just in case we search for the closest
            push!(f_TM_010[end], results[2][argmin(abs.(results[2] .- f_TM_010_true))])
            push!(f_TE_111[end], results[2][argmin(abs.(results[2] .- f_TE_111_true))])
        end
    end

    f_TM_010_rel_error = map(x -> abs.(x .- f_TM_010_true) ./ f_TM_010_true, f_TM_010)
    f_TE_111_rel_error = map(x -> abs.(x .- f_TE_111_true) ./ f_TE_111_true, f_TE_111)

    return dof, f_TM_010_rel_error, f_TE_111_rel_error
end
