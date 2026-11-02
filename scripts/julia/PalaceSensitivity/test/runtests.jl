# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using PalaceSensitivity, Test, LinearAlgebra, Random

@testset "PalaceSensitivity" begin
    rng = Random.Xoshiro(1)

    @testset "Export round trip" begin
        dir = mktempdir()
        n = 57
        X = randn(rng, 3, n)
        G1, G2 = randn(rng, 3, n), randn(rng, 3, n)
        boundary = [Int32.(rand(rng, [Int[], [4], [-1, 4], [7]])) for _ = 1:n]
        write_sensitivity(
            dir;
            problem_type="Electrostatic",
            length_unit=1e-6,
            points=X,
            objectives=[
                (
                    index=1,
                    type="Capacitance",
                    value=1.5e-13,
                    unit="F",
                    gradient=G1,
                    spec=Dict("Terminals" => [25, 25])
                ),
                (
                    index=2,
                    type="Capacitance",
                    value=-2.5e-15,
                    unit="F",
                    gradient=G2,
                    spec=Dict("Terminals" => [25, 26])
                )
            ],
            boundary=boundary,
            metadata=Dict("PalaceVersion" => "test")
        )
        s = read_sensitivity(dir)
        @test s.problem_type == "Electrostatic"
        @test s.length_unit == 1e-6
        @test s.points == X
        @test objective(s, 1).gradient == G1
        @test objective(s, 2).gradient == G2
        @test objective(s, 2).value == -2.5e-15
        @test objective(s, 2).unit == "F"
        @test objective(s, 1).spec["Terminals"] == [25, 25]
        @test s.boundary == boundary
        @test s.metadata == Dict("PalaceVersion" => "test")
        @test_throws KeyError objective(s, 3)
        @test sprint(show, s) ==
              "SensitivityExport(Electrostatic, 57 points, objectives [1, 2])"

        # Invalid exports are rejected.
        write(joinpath(dir, "boundary.bin"), fill(Int64(5), n + 1))
        @test_throws r"not a boundary file" read_sensitivity(dir)
        write(joinpath(dir, "points.bin"), zeros(10))
        @test_throws r"expected 1368" read_sensitivity(dir)
        other = mktempdir()
        write(joinpath(other, "sensitivity.json"), "{\"Format\": \"Other\"}")
        @test_throws r"not a PalaceSensitivity file" read_sensitivity(other)
        write(joinpath(other, "sensitivity.json"), "{\"Format\": \"PalaceSensitivity\"}")
        @test_throws r"version nothing" read_sensitivity(other)

        # So are objectives and metadata that would overwrite each other or the format's
        # entries.
        o = (index=1, type="Capacitance", value=0.0, gradient=G1)
        w(; kwargs...) = write_sensitivity(
            mktempdir();
            problem_type="Electrostatic",
            length_unit=1e-6,
            points=X,
            kwargs...
        )
        @test_throws ArgumentError w(objectives=[o, o])
        @test_throws ArgumentError w(objectives=[merge(o, (spec=Dict("Value" => 1),))])
        @test_throws ArgumentError w(objectives=[o], metadata=Dict("Version" => 2))
    end

    @testset "Exports without node gradients" begin
        dir = mktempdir()
        write_sensitivity(
            dir;
            problem_type="Electrostatic",
            length_unit=1e-6,
            points=zeros(3, 0),
            objectives=[(index=1, type="Capacitance", value=1.0, gradient=zeros(3, 0))]
        )
        @test sort(readdir(dir)) == ["sensitivity.json"]
        s = read_sensitivity(dir)
        @test size(s.points) == (3, 0)
        @test objective(s, 1).value == 1.0
        @test size(objective(s, 1).gradient) == (3, 0)
    end

    @testset "Admissible velocities" begin
        # A unit square of points: attribute 1 on x = 0, 2 on x = 1, and 3 on y = 0.
        X = reduce(hcat, [[x, y, 0.0] for x = 0:0.25:1 for y = 0:0.25:1])
        boundary = [Int32[] for _ in axes(X, 2)]
        for (k, x) in enumerate(eachcol(X))
            x[1] == 0 && push!(boundary[k], 1)
            x[1] == 1 && push!(boundary[k], 2)
            x[2] == 0 && push!(boundary[k], 3)
        end
        dir = mktempdir()
        write_sensitivity(
            dir;
            problem_type="Electrostatic",
            length_unit=1e-6,
            points=X,
            objectives=[(index=1, type="Capacitance", value=0.0, gradient=zero(X))],
            boundary=boundary
        )
        s = read_sensitivity(dir)
        bulge(Y) = vcat((Y[1:1, :] .* (1 .- Y[1:1, :])), zeros(2, size(Y, 2)))
        rotation(Y) = vcat(-Y[2:2, :], Y[1:1, :], zeros(1, size(Y, 2)))
        stretch(Y) = vcat(zeros(1, size(Y, 2)), Y[2:2, :], zeros(1, size(Y, 2)))
        @test isnothing(check_velocity(s, bulge; fixed=(1, 2)))
        @test isnothing(check_velocity(s, rotation; rigid=(1, 2, 3)))
        @test_throws r"vanish on boundary attribute 3" check_velocity(
            s,
            bulge;
            fixed=(1, 3)
        )
        @test_throws r"attribute 2 rigidly" check_velocity(s, stretch; rigid=(3, 2))
        @test isnothing(check_velocity(s, stretch; rigid=(3,)))
    end

    @testset "Contraction with design velocities" begin
        # J(X) = Σ wₖ ⋅ xₖ has gradient wₖ and derivative Σ wₖ ⋅ V(xₖ) along any velocity V.
        n = 200
        X = randn(rng, 3, n)
        W = randn(rng, 3, n)
        A, b = randn(rng, 3, 3), randn(rng, 3)
        V(Y) = A * Y .+ b
        U(Y) = 2 .* Y
        dir = mktempdir()
        write_sensitivity(
            dir;
            problem_type="Electrostatic",
            length_unit=1e-6,
            points=X,
            objectives=[(index=7, type="Capacitance", value=0.0, gradient=W)]
        )
        s = read_sensitivity(dir)
        @test derivative(s, 7, V) ≈ dot(W, A * X .+ b)
        @test jacobian(s, (a=V, b=U)) ≈ [dot(W, A * X .+ b) dot(W, 2X)]
        @test jacobian(s, [V, U]) ≈ jacobian(s, (a=V, b=U))
        # The same export in mm, with velocities in μm.
        dirmm = mktempdir()
        write_sensitivity(
            dirmm;
            problem_type="Electrostatic",
            length_unit=1e-3,
            points=X ./ 1000,
            objectives=[(index=7, type="Capacitance", value=0.0, gradient=W .* 1000)]
        )
        smm = read_sensitivity(dirmm)
        @test derivative(smm, 7, V) ≈ derivative(s, 7, V)
        @test_throws DimensionMismatch jacobian(s, (Y -> zeros(2, 3),))
    end
end
