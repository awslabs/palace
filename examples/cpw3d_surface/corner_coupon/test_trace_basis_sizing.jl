# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
include(joinpath(@__DIR__, "mesh_corner_coupon.jl"))

# The trace-basis knot sizing law (supervisor decisions 328 (2), 783 / 787): around every
# free knot the mesh size is min(the knot's along-ring gap, the across-ring spacing of its
# ring to the neighbouring rings), both read from the trace layout itself. A hat on the
# extra ring levels above the metal has a support of gap x spacing; when gap > spacing the
# along-ring gap alone no longer puts two elements across the hat (the concave 75-degree
# node of decision 783: hat 105 at 79 active nodes vs 80 on the h-2 calibration mesh), so
# the size must read the spacing too. The test builds the sizing fields on a synthetic
# layout (no mesh) and reads every knot's size back from the Threshold / Distance fields.

const R = 2.0
const TOLERANCE = 1.0e-9

# A square ring of `count` vertices at half width `half_width` and height `z`, in perimeter
# order; `conductor_of(index)` marks the PEC knots (not free).
function square_ring_vertices(half_width, z, count; conductor_of = index -> 0)
    rows = String[]
    perimeter = 8.0half_width
    for index = 0:(count - 1)
        s = perimeter * index / count
        # Perimeter coordinate s counterclockwise from (-half_width, 0), as the mesher's
        # square_perimeter_coordinate.
        if s <= half_width
            x, y = -half_width, -s
        elseif s <= 3half_width
            x, y = -half_width + (s - half_width), -half_width
        elseif s <= 5half_width
            x, y = half_width, -half_width + (s - 3half_width)
        elseif s <= 7half_width
            x, y = half_width - (s - 5half_width), half_width
        else
            x, y = -half_width, half_width - (s - 7half_width)
        end
        push!(rows, "$x,$y,$z,1,$(conductor_of(index))")
    end
    return rows
end

function write_trace_vertices(directory, rows)
    open(joinpath(directory, "trace-vertices.csv"), "w") do io
        println(io, "vertex,x,y,z,basis,conductor")
        for (index, row) in enumerate(rows)
            println(io, "$index,$row")
        end
    end
end

# Knot sizes read back from the fields: point coordinates => SizeMin of the Threshold field
# whose Distance sources include the point.
function knot_sizes(fields)
    sizes = Dict{NTuple{3, Float64}, Float64}()
    for threshold in fields
        distance = round(Int, gmsh.model.mesh.field.getNumber(threshold, "InField"))
        size = gmsh.model.mesh.field.getNumber(threshold, "SizeMin")
        for tag in gmsh.model.mesh.field.getNumbers(distance, "PointsList")
            x, y, z = gmsh.model.getValue(0, round(Int, tag), Float64[])
            sizes[(x, y, z)] = size
        end
    end
    return sizes
end

function sizes_on_ring(sizes, half_width, z)
    return [
        size for (point, size) in sizes if
        abs(max(abs(point[1]), abs(point[2])) - half_width) <= TOLERANCE &&
        abs(point[3] - z) <= TOLERANCE
    ]
end

function build_fields(rows; lc_far = 0.3R)
    directory = mktempdir()
    write_trace_vertices(directory, rows)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("trace-basis-sizing")
    try
        fields, smallest =
            trace_basis_size_fields(gmsh.model.occ, directory, R, lc_far, 3)
        return knot_sizes(fields), smallest
    finally
        gmsh.finalize()
    end
end

@testset "trace basis sizing law (decision 787 R-A)" begin
@testset "the knot size reads the across-ring spacing" begin
    # The layout: the box |x|, |y|, |z| <= R. Side rings (half width R) at z = -R (the
    # bottom edge ring, 8 knots: gap R), z = 0 and z = 0.05 R (32 knots each: gap 0.25 R,
    # four PEC knots on z = 0), z = +R (the top edge ring, 32 knots: gap 0.25 R). Cap rings
    # (half width R / 3) at z = +R with a cap centre and at z = -R without one (4 knots
    # each, at the side midpoints: the along-ring gap is the chord sqrt(2) R / 3 across the
    # corner).
    rows = vcat(
        square_ring_vertices(R, -R, 8),
        square_ring_vertices(R, 0.0, 32; conductor_of = index -> index % 8 == 0 ? 1 : 0),
        square_ring_vertices(R, 0.05R, 32),
        square_ring_vertices(R, R, 32),
        square_ring_vertices(R / 3, R, 4),
        square_ring_vertices(R / 3, -R, 4),
        ["0.0,0.0,$R,0,0"]  # the top cap centre: an interpolated vertex, basis 0
    )
    sizes, smallest = build_fields(rows)
    @test length(sizes) == 8 + 28 + 32 + 32 + 4 + 4

    # The hat whose support is gap x spacing with gap > spacing: the two side rings 0.05 R
    # apart with knots 0.25 R apart along the ring get a size <= the spacing (the law of
    # record gave the along-ring gap 0.25 R: this is the assertion that FAILS on the old
    # law).
    for z in (0.0, 0.05R)
        ring = sizes_on_ring(sizes, R, z)
        @test length(ring) == (z == 0.0 ? 28 : 32)
        @test all(size -> size <= 0.05R + TOLERANCE, ring)
        @test all(size -> isapprox(size, 0.05R; atol = TOLERANCE), ring)
    end
    @test isapprox(smallest, 0.05R; atol = TOLERANCE)

    # A knot whose along-ring gap is already the smaller keeps its gap: the top edge ring
    # (gap 0.25 R; spacings 0.95 R to the side ring, 2 R / 3 to the cap ring).
    @test all(size -> isapprox(size, 0.25R; atol = TOLERANCE), sizes_on_ring(sizes, R, R))
    # The bottom cap ring has no cap centre: its only neighbour is the edge ring 2 R / 3
    # away, more than its gap sqrt(2) R / 3, so its size is unchanged.
    @test all(
        size -> isapprox(size, sqrt(2.0) * R / 3; atol = TOLERANCE),
        sizes_on_ring(sizes, R / 3, -R)
    )

    # The edge ring belongs to both groups and takes the smaller spacing: the bottom edge
    # ring (gap R) has the side ring z = 0 at R and the cap ring at 2 R / 3 -> 2 R / 3.
    @test all(size -> isapprox(size, 2R / 3; atol = TOLERANCE), sizes_on_ring(sizes, R, -R))

    # The cap centre counts as the innermost cap level: the top cap ring (gap sqrt(2) R / 3)
    # is R / 3 from the centre -> R / 3.
    @test all(
        size -> isapprox(size, R / 3; atol = TOLERANCE),
        sizes_on_ring(sizes, R / 3, R)
    )

    # free_knot_gaps still reads the along-ring gaps alone.
    vertices = read_trace_vertices(
        let directory = mktempdir()
            write_trace_vertices(directory, rows)
            directory
        end
    )
    gaps = free_knot_gaps(vertices, trace_rings(vertices, R), R)
    @test length(gaps) == 108
    @test isapprox(minimum(values(gaps)), 0.25R; atol = TOLERANCE)
    @test isapprox(maximum(values(gaps)), R; atol = TOLERANCE)
end

@testset "a ring on neither box surface fails closed" begin
    rows = vcat(
        square_ring_vertices(R, -R, 8),
        square_ring_vertices(R, 0.0, 32),
        square_ring_vertices(R, R, 32),
        square_ring_vertices(0.5R, 0.0, 8)  # neither half width R nor height +-R
    )
    @test_throws ErrorException build_fields(rows)
end
end
