# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""
    PalaceSensitivity

Read a Palace sensitivity export and contract its gradients with design velocities, giving
derivatives with respect to design parameters, `dJ/dp = Σₖ gₖ ⋅ V(xₖ)`. The README specifies
the export format.
"""
module PalaceSensitivity

using LinearAlgebra
import JSON

export read_sensitivity, write_sensitivity, objective, derivative, jacobian, check_velocity

const FORMAT = "PalaceSensitivity"
const FORMAT_VERSION = 1
# Entries of `sensitivity.json` and of its objectives that have their own fields.
const EXPORT_KEYS = (
    "Format",
    "Version",
    "ProblemType",
    "LengthUnit",
    "NumPoints",
    "Points",
    "ByteOrder",
    "Objectives",
    "Boundary"
)
const OBJECTIVE_KEYS = ("Index", "Type", "Value", "Unit", "Gradient")

"""
    struct Objective

An objective of a [`SensitivityExport`](@ref): its `index` and `type`, the rest of its
configuration in `spec` (for example `"Terminals"`), its `value` in `unit`, and its
`gradient` with respect to the exported points, a `3 × N` matrix in `unit` per mesh unit.
"""
struct Objective
    index::Int
    type::String
    spec::Dict{String, Any}
    value::Float64
    unit::String
    gradient::Matrix{Float64}
end

"""
    struct SensitivityExport

A sensitivity export, from [`read_sensitivity`](@ref): the exported `points` (a `3 × N`
matrix in mesh units), the mesh unit `length_unit` in meters, the `objectives`, the
`boundary` attributes of each point (or `nothing`), and the other entries of
`sensitivity.json` in `metadata`.
"""
struct SensitivityExport
    problem_type::String
    length_unit::Float64
    points::Matrix{Float64}
    objectives::Vector{Objective}
    boundary::Union{Nothing, Vector{Vector{Int32}}}
    metadata::Dict{String, Any}
end

function Base.show(io::IO, s::SensitivityExport)
    return print(
        io,
        "SensitivityExport(",
        s.problem_type,
        ", ",
        size(s.points, 2),
        " points, objectives ",
        [o.index for o in s.objectives],
        ")"
    )
end

"""
    objective(s::SensitivityExport, index)

The [`Objective`](@ref) of `s` with the given index.
"""
function objective(s::SensitivityExport, index::Integer)
    i = findfirst(o -> o.index == index, s.objectives)
    isnothing(i) && throw(KeyError(index))
    return s.objectives[i]
end

"""
    read_sensitivity(dir)

Read the sensitivity export in `dir` (a run's `<Problem.Output>/sensitivity`) as a
[`SensitivityExport`](@ref).
"""
function read_sensitivity(dir::AbstractString)
    path = joinpath(dir, "sensitivity.json")
    meta = JSON.parsefile(path)
    get(meta, "Format", nothing) == FORMAT || error("$path is not a $FORMAT file")
    version = get(meta, "Version", nothing)
    version == FORMAT_VERSION ||
        error("Unsupported $FORMAT version $version in $path (supported: $FORMAT_VERSION)")
    byte_order = get(meta, "ByteOrder", "LittleEndian")
    byte_order == "LittleEndian" || error("Unsupported byte order $byte_order in $path")
    n = Int(meta["NumPoints"])
    points = _read_vectors(dir, get(meta, "Points", nothing), n)
    objectives = Objective[]
    for o in meta["Objectives"]
        spec = Dict{String, Any}(String(k) => v for (k, v) in o if !(k in OBJECTIVE_KEYS))
        gradient = _read_vectors(dir, get(o, "Gradient", nothing), n)
        push!(
            objectives,
            Objective(
                Int(o["Index"]),
                String(o["Type"]),
                spec,
                Float64(o["Value"]),
                String(get(o, "Unit", "")),
                gradient
            )
        )
    end
    boundary =
        haskey(meta, "Boundary") ? _read_boundary(joinpath(dir, meta["Boundary"]), n) :
        nothing
    metadata = Dict{String, Any}(String(k) => v for (k, v) in meta if !(k in EXPORT_KEYS))
    return SensitivityExport(
        String(meta["ProblemType"]),
        Float64(meta["LengthUnit"]),
        points,
        objectives,
        boundary,
        metadata
    )
end

# N × 3 row-major Float64 is 3 × N column-major. Exports without node gradients have no
# files and no points.
function _read_vectors(dir, file, n)
    if isnothing(file)
        n == 0 || error("Missing points or gradient file in the export in $dir")
        return zeros(3, 0)
    end
    path = joinpath(dir, file)
    filesize(path) == 24n ||
        error("$path has $(filesize(path)) bytes, expected $(24n) for $n points")
    v = Vector{Float64}(undef, 3n)
    read!(path, v)
    return reshape(ltoh.(v), 3, n)
end

function _write_vectors(path, X)
    size(X, 1) == 3 || throw(DimensionMismatch("expected a 3 × N matrix"))
    open(path, "w") do io
        return write(io, htol.(vec(Matrix{Float64}(X))))
    end
end

function _read_boundary(path, n)
    invalid() = error("$path is not a boundary file for $n points")
    filesize(path) >= 8(n + 1) || invalid()
    open(path) do io
        offsets = ltoh.(read!(io, Vector{Int64}(undef, n + 1)))
        (
            offsets[1] == 0 &&
            issorted(offsets) &&
            filesize(path) == 8(n + 1) + 4offsets[end]
        ) || invalid()
        attributes = ltoh.(read!(io, Vector{Int32}(undef, offsets[end])))
        return [attributes[(offsets[k] + 1):offsets[k + 1]] for k = 1:n]
    end
end

"""
    write_sensitivity(dir; problem_type, length_unit, points, objectives, boundary=nothing,
                      metadata=Dict())

Write a sensitivity export to `dir`, for tests and other tools.

`points` is a `3 × N` matrix in mesh units, and `length_unit` the mesh unit in meters. Each
objective is an [`Objective`](@ref) or a `NamedTuple` with its fields (`unit` and `spec` are
optional). `boundary` lists the boundary attributes of each point, and `metadata` adds
entries to `sensitivity.json`. Without points (a `3 × 0` matrix), no gradients are written,
as by *Palace* without `"Shape"`.
"""
function write_sensitivity(
    dir::AbstractString;
    problem_type,
    length_unit,
    points,
    objectives,
    boundary=nothing,
    metadata=Dict{String, Any}()
)
    allunique(o.index for o in objectives) ||
        throw(ArgumentError("objective indices must be unique"))
    any(k -> String(k) in EXPORT_KEYS, keys(metadata)) &&
        throw(ArgumentError("metadata must not contain the entries $EXPORT_KEYS"))
    mkpath(dir)
    n = size(points, 2)
    n > 0 && _write_vectors(joinpath(dir, "points.bin"), points)
    entries = Any[]
    for o in objectives
        size(o.gradient) == (3, n) ||
            throw(DimensionMismatch("the gradient of objective $(o.index) is not 3 × $n"))
        spec = hasproperty(o, :spec) ? o.spec : Dict{String, Any}()
        any(k -> String(k) in OBJECTIVE_KEYS, keys(spec)) &&
            throw(ArgumentError("spec must not contain the entries $OBJECTIVE_KEYS"))
        entry = Dict{String, Any}(
            "Index" => o.index,
            "Type" => o.type,
            "Value" => o.value,
            "Unit" => hasproperty(o, :unit) ? o.unit : ""
        )
        if n > 0
            entry["Gradient"] = "gradient-$(o.index).bin"
            _write_vectors(joinpath(dir, entry["Gradient"]), o.gradient)
        end
        push!(entries, merge(entry, spec))
    end
    meta = Dict{String, Any}(
        "Format" => FORMAT,
        "Version" => FORMAT_VERSION,
        "ProblemType" => problem_type,
        "LengthUnit" => length_unit,
        "NumPoints" => n,
        "ByteOrder" => "LittleEndian",
        "Objectives" => entries
    )
    n > 0 && (meta["Points"] = "points.bin")
    if !isnothing(boundary)
        length(boundary) == n || throw(
            DimensionMismatch("boundary has $(length(boundary)) entries, expected $n")
        )
        offsets = cumsum(vcat(0, length.(boundary)))
        open(joinpath(dir, "boundary.bin"), "w") do io
            write(io, htol.(Int64.(offsets)))
            return write(io, htol.(Int32.(reduce(vcat, boundary; init=Int32[]))))
        end
        meta["Boundary"] = "boundary.bin"
    end
    merge!(meta, metadata)
    open(joinpath(dir, "sensitivity.json"), "w") do io
        return JSON.print(io, meta, 2)
    end
    return dir
end

"""
    jacobian(s::SensitivityExport, velocities; objectives=all, velocity_length_unit=1e-6)

The derivatives `Σₖ gₖ ⋅ V(xₖ)` of the objectives of `s` (rows, in the order of
`objectives`) with respect to the parameters of `velocities` (columns, in iteration order).
A velocity maps a `3 × N` matrix of points to a `3 × N` matrix of velocities, both with
lengths in `velocity_length_unit` meters (the default, μm, is that of DeviceLayout's design
velocities).
"""
function jacobian(
    s::SensitivityExport,
    velocities;
    objectives=[o.index for o in s.objectives],
    velocity_length_unit=1e-6
)
    scale = s.length_unit / velocity_length_unit # mesh unit in velocity units
    X = s.points .* scale
    V = [v(X) ./ scale for v in values(velocities)]
    for (j, Vj) in enumerate(V)
        size(Vj) == size(X) || throw(
            DimensionMismatch(
                "velocity $j returned size $(size(Vj)) for $(size(X, 2)) points"
            )
        )
    end
    rows = collect(objectives)
    J = zeros(length(rows), length(V))
    for (i, index) in enumerate(rows)
        g = objective(s, index).gradient
        for (j, Vj) in enumerate(V)
            J[i, j] = dot(g, Vj)
        end
    end
    return J
end

"""
    derivative(s::SensitivityExport, index, velocity; velocity_length_unit=1e-6)

The derivative of objective `index` of `s` with respect to the parameter of `velocity`. See
[`jacobian`](@ref).
"""
function derivative(
    s::SensitivityExport,
    index::Integer,
    velocity;
    velocity_length_unit=1e-6
)
    return only(jacobian(s, (velocity,); objectives=(index,), velocity_length_unit))
end

"""
    check_velocity(s::SensitivityExport, velocity; fixed=(), rigid=(), rtol=1e-8,
                   velocity_length_unit=1e-6)

Check that `velocity` is admissible for the export `s`, from the boundary attributes of its
points: it must vanish on the nodes of each attribute in `fixed`, such as boundaries that do
not move, and move the nodes of each attribute in `rigid` rigidly, such as ports and lumped
elements, whose impedances depend on their size. Violations larger than `rtol` times the
largest velocity throw an error naming the attribute.
"""
function check_velocity(
    s::SensitivityExport,
    velocity;
    fixed=(),
    rigid=(),
    rtol=1e-8,
    velocity_length_unit=1e-6
)
    isnothing(s.boundary) &&
        error("The export has no boundary attributes: set \"BoundaryTags\" in \"Shape\"")
    scale = s.length_unit / velocity_length_unit
    X = s.points .* scale
    V = velocity(X)
    size(V) == size(X) || throw(
        DimensionMismatch("velocity returned size $(size(V)) for $(size(X, 2)) points")
    )
    tol = rtol * maximum(norm, eachcol(V); init=0.0)
    nodes(attr) = [k for (k, b) in enumerate(s.boundary) if attr in b]
    for attr in fixed
        violation = maximum(k -> norm(V[:, k]), nodes(attr); init=0.0)
        violation <= tol ||
            error("The velocity does not vanish on boundary attribute $attr")
    end
    for attr in rigid
        # Residual of the least-squares rigid motion t + ω × x of the nodes.
        k = nodes(attr)
        Q = X[:, k] .- sum(X[:, k]; dims=2) / max(length(k), 1)
        A = zeros(3length(k), 6)
        for (j, q) in enumerate(eachcol(Q))
            A[(3j - 2):(3j), 1:3] = I(3)
            A[(3j - 2):(3j), 4:6] = [0 q[3] -q[2]; -q[3] 0 q[1]; q[2] -q[1] 0]
        end
        b = vec(V[:, k])
        r = isempty(k) ? b : b - A * (A \ b)
        violation = maximum(norm, eachcol(reshape(r, 3, :)); init=0.0)
        violation <= tol ||
            error("The velocity does not move boundary attribute $attr rigidly")
    end
    return nothing
end

end # module
