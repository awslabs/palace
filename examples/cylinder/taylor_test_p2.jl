# Taylor test for CURVED (P2) mesh sensitivity — cylinder TM010 radius derivative.
#
# Verifies the full-node mesh sensitivity export on an order-2 (curved) tet mesh:
#   1. Generates a P2 gmsh tet mesh of the dielectric cylinder (R=2.74, H=2.74 cm).
#   2. Runs Palace with Eigenmode MeshSensitivity → mesh_sensitivity.solb over ALL
#      nodes + mesh_sensitivity_nodes.txt sidecar (mid-edge node coordinates).
#   3. Contracts df/dx with the radial velocity dx/dR = (x/R, y/R, 0) at every node.
#   4. Taylor test: scales all node coordinates (vertices + mid-edge) by (1+ε/R)
#      in xy — the exact perturbation corresponding to dx/dR — and checks
#      |f(R+ε) - f(R) - ε·df/dR| = O(ε²).
#
# Expected: df/dR·R/f ≈ -1 (TM010: f ∝ 1/R) and convergence rates ≈ 2.0.
# Before the full-node export fix, vertex-only truncation gave df/dR·R/f ≈ -0.0007.
#
# The in-log "Dilation identity" line is the stronger, parameter-free check:
# sum(df/dx · x) = -f exactly for any mode on any mesh (Maxwell scale invariance).
#
# Usage: cd examples/cylinder && julia taylor_test_p2.jl
#        (requires `palace` on PATH and the Gmsh julia package for mesh generation)

using LinearAlgebra, Printf

include(joinpath(@__DIR__, "..", "shared", "palace_io.jl"))

const R = 2.74
const MESH = joinpath(@__DIR__, "mesh", "cylinder_taylor_p2.msh")
const CONFIG = joinpath(@__DIR__, "taylor_p2.json")
const OUTDIR = "postpro/taylor_p2_baseline"

# ─── Binary Gmsh v2 node scaling (same as taylor_test.jl) ───────────────────

function find_nodes_offset(buf::Vector{UInt8})
    target = UInt8['$', 'N', 'o', 'd', 'e', 's', '\n']
    pos = 0
    for i in 1:length(buf) - length(target) + 1
        if @view(buf[i:i+length(target)-1]) == target
            pos = i + length(target)
            break
        end
    end
    pos == 0 && error("No \$Nodes section found")
    nl = findnext(==(UInt8('\n')), buf, pos)
    n_nodes = parse(Int, String(buf[pos:nl-1]))
    return nl + 1, n_nodes
end

read_f64(buf, pos) = reinterpret(Float64, buf[pos:pos+7])[1]
write_f64!(buf, pos, val) = (buf[pos:pos+7] .= reinterpret(UInt8, Float64[val]))

function scale_mesh_xy!(buf, n_nodes, pos, xy)
    for _ in 1:n_nodes
        pos += 4  # skip tag
        write_f64!(buf, pos, read_f64(buf, pos) * xy)
        write_f64!(buf, pos + 8, read_f64(buf, pos + 8) * xy)
        pos += 24
    end
end

# ─── Setup ───────────────────────────────────────────────────────────────────

if !isfile(MESH)
    include(joinpath(@__DIR__, "mesh", "mesh.jl"))
    # height = aspect_ratio * 2 * radius; 0.5 → H = R so TM010 (2.91 GHz) is
    # well separated from TE111 (4.40 GHz) and mode 1 is unambiguous.
    generate_cylindrical_mesh(filename=basename(MESH), mesh_type=0, order=2,
                              refinement=1, aspect_ratio=0.5, verbose=2)
end

write(CONFIG, """{
  "Problem": {"Type": "Eigenmode", "Verbose": 2, "Output": "$OUTDIR"},
  "Model": {"Mesh": "$(relpath(MESH, @__DIR__))", "L0": 1.0e-2},
  "Domains": {"Materials": [{"Attributes": [1], "Permeability": 1.0, "Permittivity": 2.08, "LossTan": 0.0}]},
  "Boundaries": {"PEC": {"Attributes": [2, 3, 4]}},
  "Solver": {
    "Order": 2, "Device": "CPU",
    "Eigenmode": {"N": 2, "Tol": 1.0e-8, "Target": 2.0, "Sensitivity": true, "MeshSensitivity": true},
    "Linear": {"Type": "Default", "KSPType": "GMRES", "Tol": 1.0e-8, "MaxIts": 100}
  }
}""")

println("═══ P2 Curved-Mesh Taylor Test: Cylinder Radius ═══\n")
println("Baseline Palace run...")
run(pipeline(`palace $CONFIG`; stdout=devnull, stderr=devnull))

# ─── Contraction over ALL nodes (vertices + mid-edge) ───────────────────────

f0 = read_eigenfrequency(joinpath(OUTDIR, "eig.csv"); mode=1)
sol = read_solb(joinpath(OUTDIR, "mesh_sensitivity.solb"))
sidecar = joinpath(OUTDIR, "mesh_sensitivity_nodes.txt")
@assert isfile(sidecar) "mesh_sensitivity_nodes.txt missing — Palace build lacks the full-node export"
lines = readlines(sidecar)
nn = parse(Int, split(lines[1])[1])
coords = reduce(vcat, [parse.(Float64, split(lines[1 + i]))' for i in 1:nn])
@assert length(vec(sol.data)) == 3nn "solb/sidecar node count mismatch"

dxdR = zeros(nn, 3)
for i in 1:nn
    dxdR[i, 1] = coords[i, 1] / R
    dxdR[i, 2] = coords[i, 2] / R
end
dfdR = dot(vec(sol.data), vec(dxdR))
@printf("  f₀ = %.9f GHz   nodes = %d\n", f0, nn)
@printf("  df/dR = %.9f GHz/cm   (df/dR)·R/f = %.6f  (expect -1)\n\n", dfdR, dfdR * R / f0)

# ─── Taylor convergence ──────────────────────────────────────────────────────

buf0 = read(MESH)
epsilons = [1e-1, 1e-2, 1e-3]
rems = Float64[]
rates = Float64[]
for (k, ε) in enumerate(epsilons)
    buf = copy(buf0)
    pos, n_nodes = find_nodes_offset(buf)
    scale_mesh_xy!(buf, n_nodes, pos, 1 + ε / R)
    pm = joinpath(@__DIR__, "mesh", "cylinder_taylor_p2_eps$k.msh")
    write(pm, buf)
    cfg = replace(read(CONFIG, String),
        relpath(MESH, @__DIR__) => relpath(pm, @__DIR__),
        OUTDIR => "postpro/taylor_p2_$k",
        "\"MeshSensitivity\": true" => "\"MeshSensitivity\": false")
    pc = joinpath(@__DIR__, "taylor_p2_$k.json")
    write(pc, cfg)
    run(pipeline(`palace $pc`; stdout=devnull, stderr=devnull))
    fε = read_eigenfrequency(joinpath("postpro/taylor_p2_$k", "eig.csv"); mode=1)
    rem = abs(fε - f0 - ε * dfdR)
    push!(rems, rem)
    rate = k >= 2 ? log(rems[k-1] / rems[k]) / log(epsilons[k-1] / ε) : NaN
    k >= 2 && push!(rates, rate)
    @printf("  ε=%.0e: f=%.9f  remainder=%.3e  rate=%s\n",
            ε, fε, rem, isnan(rate) ? "—" : @sprintf("%.2f", rate))
    rm(pm; force=true); rm(pc; force=true)
    rm("postpro/taylor_p2_$k"; recursive=true, force=true)
end

ok = count(>(1.8), rates) >= 1 && abs(dfdR * R / f0 + 1) < 0.05
println(ok ? "\nStatus: PASS ✓" : "\nStatus: FAIL ✗")
exit(ok ? 0 : 1)
