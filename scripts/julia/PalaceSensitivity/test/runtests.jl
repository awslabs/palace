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
        @test s.metadata["PalaceVersion"] == "test"
        @test_throws KeyError objective(s, 3)
        @test sprint(show, s) ==
              "SensitivityExport(Electrostatic, 57 points, objectives [1, 2])"
        # Truncated binaries and other formats or versions are rejected.
        write(joinpath(dir, "points.bin"), zeros(10))
        @test_throws ErrorException read_sensitivity(dir)
        other = mktempdir()
        write(joinpath(other, "sensitivity.json"), "{\"Format\": \"Other\"}")
        @test_throws ErrorException read_sensitivity(other)
        write(
            joinpath(other, "sensitivity.json"),
            "{\"Format\": \"PalaceSensitivity\", \"Version\": 2}"
        )
        @test_throws ErrorException read_sensitivity(other)
    end

    @testset "Contraction with design velocities" begin
        # A functional that is linear in the node positions, J(X) = Σ wₖ ⋅ xₖ, has gradient
        # wₖ, and its derivative along an affine velocity V(x) = A x + b is Σ wₖ ⋅ (A xₖ + b).
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
        # The same export in mm, contracted with velocities in μm: the points are converted
        # to μm, and the velocity lengths back to mm.
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
