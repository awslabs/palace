# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""
    PalaceSensitivity

Read the sensitivity export of a Palace run, and contract it with design velocities.

A run with a `"Sensitivity"` configuration writes, for each objective, its value and its
gradient with respect to the positions of the mesh nodes in an export region.
[`read_sensitivity`](@ref) reads the export, and [`jacobian`](@ref) and
[`derivative`](@ref) contract it with design velocities `V`, which describe how the geometry
moves with a design parameter `p`, giving `dJ/dp = Σₖ gₖ ⋅ V(xₖ)`. Design velocities are any
callable mapping a `3 × N` matrix of points to a `3 × N` matrix of velocities, for example
those of `DeviceLayout.SolidModels.design_velocities`. [`write_sensitivity`](@ref) writes the
same format, for tests and for other tools.

# Format (version 1)

The export is written to `<Problem.Output>/sensitivity/`:

  - `sensitivity.json`: metadata. `"Format": "PalaceSensitivity"`, `"Version": 1`,
    `"ProblemType"`, `"LengthUnit"` (the mesh unit, `Model.L0`, in m), `"NumPoints"`,
    `"Points"` (the file of node positions), `"ByteOrder": "LittleEndian"`, and
    `"Objectives"`: for each, `"Index"`, `"Type"`, its configuration (for example
    `"Terminals"`), `"Value"`, `"Unit"` and `"Gradient"` (the file of its gradient).
    Optional entries include `"Boundary"`, `"Region"`, `"PalaceVersion"`,
    `"MeshFingerprint"`, and scalar sensitivities in `"Materials"` and `"LumpedElements"`.
  - `points.bin`: N × 3 Float64, row-major: the positions of the exported nodes, in mesh
    units.
  - `gradient-<index>.bin`: N × 3 Float64 for each objective: its gradient with respect to
    the node positions, in the objective's unit per mesh unit.
  - `boundary.bin` (optional): the boundary attributes of each point, as N + 1 Int64 offsets
    (0-based, compressed sparse rows) followed by the Int32 attributes; -1 marks the interface
    of a substructuring region.

Each node of the export region appears once. The set of points does not depend on Palace's
partition, but their order does: a gradient is a covector at its points, and is only ever
contracted with velocities evaluated at the same points, never interpolated.
"""
module PalaceSensitivity

using LinearAlgebra
import JSON

export read_sensitivity, write_sensitivity, objective, derivative, jacobian

const FORMAT = "PalaceSensitivity"
const FORMAT_VERSION = 1

"""
    struct Objective
        index::Int
        type::String
        spec::Dict{String, Any}
        value::Float64
        unit::String
        gradient::Matrix{Float64}
    end

One objective of a [`SensitivityExport`](@ref): its index and type as configured, the rest
of its configuration in `spec` (for example `"Terminals"`), its value and the value's unit,
and its gradient with respect to the exported node positions as a `3 × N` matrix, in `unit`
per mesh unit.
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
        problem_type::String
        length_unit::Float64
        points::Matrix{Float64}
        objectives::Vector{Objective}
        boundary::Union{Nothing, Vector{Vector{Int32}}}
        metadata::Dict{String, Any}
    end

The sensitivity export of a Palace run, from [`read_sensitivity`](@ref).

`points` are the exported node positions as a `3 × N` matrix in mesh units, `length_unit`
is the mesh unit in meters (`Model.L0`), and `boundary`, if exported, lists the boundary
attributes of each point. `metadata` holds the rest of `sensitivity.json`, including
material and lumped-element sensitivities.
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

Read the sensitivity export in `dir`, a run's `<Problem.Output>/sensitivity` directory, as a
[`SensitivityExport`](@ref).
"""
function read_sensitivity(dir::AbstractString)
    path = joinpath(dir, "sensitivity.json")
    meta = JSON.parsefile(path)
    get(meta, "Format", nothing) == FORMAT || error("$path is not a $FORMAT file")
    version = meta["Version"]
    version == FORMAT_VERSION ||
        error("Unsupported $FORMAT version $version in $path (supported: $FORMAT_VERSION)")
    get(meta, "ByteOrder", "LittleEndian") == "LittleEndian" ||
        error("Unsupported byte order $(meta["ByteOrder"]) in $path")
    n = Int(meta["NumPoints"])
    points = _read_vectors(joinpath(dir, meta["Points"]), n)
    objectives = Objective[]
    for o in meta["Objectives"]
        spec = Dict{String, Any}(
            String(k) => v for
            (k, v) in o if !(k in ("Index", "Type", "Value", "Unit", "Gradient"))
        )
        push!(
            objectives,
            Objective(
                Int(o["Index"]),
                String(o["Type"]),
                spec,
                Float64(o["Value"]),
                String(get(o, "Unit", "")),
                _read_vectors(joinpath(dir, o["Gradient"]), n)
            )
        )
    end
    boundary =
        haskey(meta, "Boundary") ? _read_boundary(joinpath(dir, meta["Boundary"]), n) :
        nothing
    rest = Dict{String, Any}(
        String(k) => v for
        (k, v) in meta if !(k in ("Points", "Objectives", "Boundary", "NumPoints"))
    )
    return SensitivityExport(
        String(meta["ProblemType"]),
        Float64(meta["LengthUnit"]),
        points,
        objectives,
        boundary,
        rest
    )
end

# N × 3 row-major Float64 is 3 × N column-major.
function _read_vectors(path, n)
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
    open(path) do io
        offsets = ltoh.(read!(io, Vector{Int64}(undef, n + 1)))
        attributes = ltoh.(read!(io, Vector{Int32}(undef, offsets[end])))
        eof(io) || error("$path is longer than its offsets indicate")
        return [attributes[(offsets[k] + 1):offsets[k + 1]] for k = 1:n]
    end
end

"""
    write_sensitivity(dir; problem_type, length_unit, points, objectives, boundary=nothing,
                      metadata=Dict())

Write a sensitivity export to `dir`, in the format that [`read_sensitivity`](@ref) reads.

`points` is a `3 × N` matrix in mesh units and `length_unit` the mesh unit in meters.
`objectives` is a collection of [`Objective`](@ref)s, or of `NamedTuple`s with fields
`index`, `type`, `value` and `gradient` (a `3 × N` matrix), and optionally `unit` and `spec`.
`boundary`, if given, lists the boundary attributes of each point. `metadata` entries are
added to `sensitivity.json`.
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
    mkpath(dir)
    n = size(points, 2)
    _write_vectors(joinpath(dir, "points.bin"), points)
    entries = Any[]
    for o in objectives
        size(o.gradient) == (3, n) ||
            throw(DimensionMismatch("the gradient of objective $(o.index) is not 3 × $n"))
        file = "gradient-$(o.index).bin"
        _write_vectors(joinpath(dir, file), o.gradient)
        entry = Dict{String, Any}(
            "Index" => o.index,
            "Type" => o.type,
            "Value" => o.value,
            "Unit" => hasproperty(o, :unit) ? o.unit : "",
            "Gradient" => file
        )
        hasproperty(o, :spec) && merge!(entry, o.spec)
        push!(entries, entry)
    end
    meta = Dict{String, Any}(
        "Format" => FORMAT,
        "Version" => FORMAT_VERSION,
        "ProblemType" => problem_type,
        "LengthUnit" => length_unit,
        "NumPoints" => n,
        "Points" => "points.bin",
        "ByteOrder" => "LittleEndian",
        "Objectives" => entries
    )
    if !isnothing(boundary)
        length(boundary) == n || throw(
            DimensionMismatch("boundary has $(length(boundary)) entries, expected $n")
        )
        offsets = Int64[0]
        for b in boundary
            push!(offsets, offsets[end] + length(b))
        end
        open(joinpath(dir, "boundary.bin"), "w") do io
            write(io, htol.(offsets))
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

The matrix of derivatives of the objectives of `s` (rows, in the order of `objectives`, by
default all in order) with respect to the parameters whose design velocities are
`velocities` (columns, in iteration order; for a `NamedTuple`, in the order of its keys).

Each velocity maps a `3 × N` matrix of points to a `3 × N` matrix of velocities. Its points
and the lengths of its velocities are in `velocity_length_unit` (in meters; the default, μm,
is the unit of `DeviceLayout.SolidModels.design_velocities`), and are converted from and to
the export's mesh unit. The derivatives are in each objective's unit per unit of the
parameter.
"""
function jacobian(
    s::SensitivityExport,
    velocities;
    objectives=[o.index for o in s.objectives],
    velocity_length_unit=1e-6
)
    # Points in the velocities' unit, and velocity lengths back in mesh units.
    to_velocity_unit = s.length_unit / velocity_length_unit
    X = s.points .* to_velocity_unit
    V = [v(X) ./ to_velocity_unit for v in values(velocities)]
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

The derivative `Σₖ gₖ ⋅ V(xₖ)` of objective `index` of `s` with respect to the parameter whose
design velocity is `velocity`. See [`jacobian`](@ref).
"""
function derivative(
    s::SensitivityExport,
    index::Integer,
    velocity;
    velocity_length_unit=1e-6
)
    return only(jacobian(s, (velocity,); objectives=(index,), velocity_length_unit))
end

end # module
