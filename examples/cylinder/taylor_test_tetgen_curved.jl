# End-to-end curved-mesh gradient: EGADS → TetGen → P2 curving → Palace.
#
# This is the full alternative meshing pathway exercised on the cylinder:
#   1. EGADS cylinder primitive (R = 2.74 cm, H = 2.74 cm).
#   2. TetGen graded linear mesh (EGADSTetGen.tetrahedralize_size_field).
#   3. Vertex → CAD classification (EGADSTetGen.classify_tetgen_vertices).
#   4. Continuation curving to a valid order-2 mesh (EGADSTetGen.curve_tetgen_mesh).
#   5. Palace eigensolve + analytic mesh sensitivity (full-node export).
#   6. Chain contraction df/dR = df/dx · dx/dR and Taylor verification.
#
# Expected: TM010 with εr = 2.08 → f ≈ 2.905 GHz, df/dR·R/f ≈ -1, Taylor
# remainder rates ≈ 2. Run from examples/cylinder with `palace` on PATH and
# EGADS_LIBRARY_PATH (or Preferences) set; needs EGADS.jl + TetGen.jl.
#
# Usage: julia --project=<env with EGADS+TetGen> taylor_test_tetgen_curved.jl

using LinearAlgebra, Printf
using EGADS, TetGen

include(joinpath(@__DIR__, "..", "shared", "palace_io.jl"))
include(joinpath(@__DIR__, "..", "shared", "egads_tetgen.jl"))
using .EGADSTetGen

const R = 2.74
const H = 2.74
const MESH = joinpath(@__DIR__, "mesh", "cyl_tetgen_curved.msh2")
const OUTDIR = "postpro/tetgen_curved"

# ─── 1-4: EGADS → TetGen → classify → curve ──────────────────────────────────

ctx = EGADS.EGADSContext()
EGADS.set_out_level!(ctx, 0)
cyl = EGADS.make_cylinder(ctx, [0.0, 0.0, 0.0], [0.0, 0.0, H], R)

println("TetGen mesh...")
out = EGADSTetGen.tetrahedralize_size_field([cyl], "pAYQ";
    meshsize_callback = (x, y, z) -> 0.8, method = :callback,
    max_iterations = 3, angle = 30.0)
println("  ", size(out.pointlist, 2), " vertices, ",
        size(out.tetrahedronlist, 2), " tets")

println("Classify + curve...")
cls = EGADSTetGen.classify_tetgen_vertices(out, cyl)
curved = EGADSTetGen.curve_tetgen_mesh(out, cls, cyl; verbose = true)
@assert curved.t_final == 1.0 "curving did not reach t=1"
@assert curved.min_detj_ratio > 0.0 "curved mesh invalid"

nfaces = length(EGADS.get_body_topos(cyl, EGADS.OClass.FACE))
EGADSTetGen.write_msh2_p2(MESH, curved; tet_attr = 1, tri_attr = m -> m,
    physical_names = vcat([(3, 1, "volume")],
                          [(2, m, "face_$m") for m in 1:nfaces]))
println("  wrote ", MESH)

# ─── 5: Palace with mesh sensitivity ─────────────────────────────────────────

cfg = joinpath(@__DIR__, "tetgen_curved.json")
write(cfg, """{
  "Problem": {"Type": "Eigenmode", "Verbose": 2, "Output": "$OUTDIR"},
  "Model": {"Mesh": "$(relpath(MESH, @__DIR__))", "L0": 1.0e-2},
  "Domains": {"Materials": [{"Attributes": [1], "Permeability": 1.0, "Permittivity": 2.08, "LossTan": 0.0}]},
  "Boundaries": {"PEC": {"Attributes": $(collect(1:nfaces))}},
  "Solver": {
    "Order": 2, "Device": "CPU",
    "Eigenmode": {"N": 2, "Tol": 1.0e-8, "Target": 2.0, "Sensitivity": true, "MeshSensitivity": true},
    "Linear": {"Type": "Default", "KSPType": "GMRES", "Tol": 1.0e-8, "MaxIts": 100}
  }
}""")
println("Palace eigensolve + mesh sensitivity...")
run(pipeline(`palace $cfg`; stdout = devnull, stderr = devnull))

f0 = read_eigenfrequency(joinpath(OUTDIR, "eig.csv"); mode = 1)
sol = read_solb(joinpath(OUTDIR, "mesh_sensitivity.solb"))
sidecar = joinpath(OUTDIR, "mesh_sensitivity_nodes.txt")
@assert isfile(sidecar) "full-node sensitivity sidecar missing"
lines = readlines(sidecar)
nn = parse(Int, split(lines[1])[1])
coords = reduce(vcat, [parse.(Float64, split(lines[1 + i]))' for i in 1:nn])

# ─── 6: chain contraction + Taylor ───────────────────────────────────────────

dxdR = zeros(nn, 3)
for i in 1:nn
    dxdR[i, 1] = coords[i, 1] / R
    dxdR[i, 2] = coords[i, 2] / R
end
dfdR = dot(vec(sol.data), vec(dxdR))
@printf("\n  f₀ = %.9f GHz (analytic TM010: %.6f)\n", f0,
        2.99792458e10 * 2.404825557695773 / (2π * R * sqrt(2.08)) / 1e9)
@printf("  df/dR = %.9f GHz/cm   (df/dR)·R/f = %.6f  (expect -1)\n\n", dfdR, dfdR * R / f0)

# Taylor: pure xy-scaling of every node (matches dx/dR exactly).
epsilons = [1e-1, 1e-2, 1e-3]
rems = Float64[]
rates = Float64[]
base = read(MESH, String)
for (k, ε) in enumerate(epsilons)
    s = 1 + ε / R
    # ASCII msh2: rescale node x, y in the $Nodes block.
    io = IOBuffer()
    in_nodes = false
    for ln in split(base, '\n')
        if ln == "\$Nodes"
            in_nodes = true
            println(io, ln)
        elseif ln == "\$EndNodes"
            in_nodes = false
            println(io, ln)
        elseif in_nodes && count(' ', ln) == 3
            t = split(ln)
            println(io, t[1], " ", parse(Float64, t[2]) * s, " ",
                    parse(Float64, t[3]) * s, " ", t[4])
        else
            println(io, ln)
        end
    end
    pm = joinpath(@__DIR__, "mesh", "cyl_tetgen_curved_eps$k.msh2")
    write(pm, String(take!(io)))
    pc = joinpath(@__DIR__, "tetgen_curved_eps$k.json")
    write(pc, replace(read(cfg, String),
        relpath(MESH, @__DIR__) => relpath(pm, @__DIR__),
        OUTDIR => "postpro/tetgen_curved_$k",
        "\"MeshSensitivity\": true" => "\"MeshSensitivity\": false"))
    run(pipeline(`palace $pc`; stdout = devnull, stderr = devnull))
    fε = read_eigenfrequency(joinpath("postpro/tetgen_curved_$k", "eig.csv"); mode = 1)
    rem = abs(fε - f0 - ε * dfdR)
    push!(rems, rem)
    rate = k >= 2 ? log(rems[k - 1] / rems[k]) / log(epsilons[k - 1] / ε) : NaN
    k >= 2 && push!(rates, rate)
    @printf("  ε=%.0e: f=%.9f  remainder=%.3e  rate=%s\n",
            ε, fε, rem, isnan(rate) ? "—" : @sprintf("%.2f", rate))
    rm(pm; force = true); rm(pc; force = true)
    rm("postpro/tetgen_curved_$k"; recursive = true, force = true)
end

close(ctx)
ok = count(>(1.8), rates) >= 1 && abs(dfdR * R / f0 + 1) < 0.05
println(ok ? "\nStatus: PASS ✓" : "\nStatus: FAIL ✗")
exit(ok ? 0 : 1)
