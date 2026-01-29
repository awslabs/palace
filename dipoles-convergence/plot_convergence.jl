using Plots
using Printf

# Data structure to hold convergence results
struct ConvergenceData
    p::Int                    # polynomial order
    n::Int                    # mesh refinement level
    nd_unknowns::Int         # Number of ND unknowns
    e_rel_error::Float64     # E-field relative error (selective domain)
    b_rel_error::Float64     # B-field relative error (selective domain)
end

function parse_log_file(filename::String)
    """Parse a convergence log file and extract relevant data."""

    lines = readlines(filename)

    # Extract p and n from filename (e.g., "conv-f200-p2-n7.log" -> p=2, n=7)
    m_p = match(r"p(\d+)", filename)
    m_n = match(r"n(\d+)", filename)
    p = parse(Int, m_p.captures[1])
    n = parse(Int, m_n.captures[1])

    nd_unknowns = 0
    e_rel_error = 0.0
    b_rel_error = 0.0

    in_selective_e = false
    in_selective_b = false

    for (i, line) in enumerate(lines)
        # Extract ND unknowns
        if contains(line, "ND (p =")
            m = match(r"ND \(p = \d+\): (\d+)", line)
            if m !== nothing
                nd_unknowns = parse(Int, m.captures[1])
            end
        end

        # Find selective E-field section
        if contains(line, "--- E-FIELD: Selective Elements")
            in_selective_e = true
            in_selective_b = false
        elseif contains(line, "--- B-FIELD: Selective Elements")
            in_selective_e = false
            in_selective_b = true
        elseif contains(line, "---") && (in_selective_e || in_selective_b)
            in_selective_e = false
            in_selective_b = false
        end

        # Extract relative errors from selective domain
        if in_selective_e && contains(line, "Total:") && contains(line, "||E||")
            m = match(r"Total:\s+.*=\s+([\d.e+-]+)", line)
            if m !== nothing
                e_rel_error = parse(Float64, m.captures[1])
            end
        end

        if in_selective_b && contains(line, "Total:") && contains(line, "||B||")
            m = match(r"Total:\s+.*=\s+([\d.e+-]+)", line)
            if m !== nothing
                b_rel_error = parse(Float64, m.captures[1])
            end
        end
    end

    return ConvergenceData(p, n, nd_unknowns, e_rel_error, b_rel_error)
end

# Parse all log files
log_files = [
    # p1 files
    "conv-f200-p1-n7.log",
    "conv-f200-p1-n10.log",
    "conv-f200-p1-n13.log",
    # p2 files
    "conv-f200-p2-n5.log",
    "conv-f200-p2-n7.log",
    "conv-f200-p2-n10.log",
    "conv-f200-p2-n13.log",
    # p3 files
    "conv-f200-p3-n5.log",
    "conv-f200-p3-n7.log",
    "conv-f200-p3-n10.log"
]

data = ConvergenceData[]
for file in log_files
    if isfile(file)
        push!(data, parse_log_file(file))
        println("Parsed $file: p=$(data[end].p), n=$(data[end].n), ND=$(data[end].nd_unknowns), E_err=$(data[end].e_rel_error), B_err=$(data[end].b_rel_error)")
    else
        @warn "File not found: $file"
    end
end

# Group data by polynomial order
data_by_p = Dict{Int, Vector{ConvergenceData}}()
for d in data
    if !haskey(data_by_p, d.p)
        data_by_p[d.p] = ConvergenceData[]
    end
    push!(data_by_p[d.p], d)
end

# Sort each group by number of unknowns
for p in keys(data_by_p)
    sort!(data_by_p[p], by = d -> d.nd_unknowns)
end

# Define colors and markers for each polynomial order
colors = Dict(1 => :blue, 2 => :red, 3 => :green)
markers = Dict(1 => :circle, 2 => :square, 3 => :diamond)

# Create E-field plot
p_e = plot(
    xscale=:log10,
    yscale=:log10,
    xlabel="DoFs",
    ylabel="Relative L₂ Error",
    title="Electric Field (f=200 MHz)",
    legend=:topright,
    grid=true,
    minorgrid=true,
    size=(800, 600),
    dpi=300)

# Plot E-field data for each polynomial order
for p in sort(collect(keys(data_by_p)))
    p_data = data_by_p[p]
    unknowns = [d.nd_unknowns for d in p_data]
    e_errors = [d.e_rel_error for d in p_data]

    plot!(p_e, unknowns, e_errors,
          label="P=$p",
          marker=markers[p],
          markersize=6,
          linewidth=2,
          color=colors[p])

    # Add reference line with slope -p for E-field
    max_e_idx = argmax(e_errors)
    x0_e = unknowns[max_e_idx]
    y0_e = e_errors[max_e_idx]
    x_ref_e = range(x0_e, maximum(unknowns), length=100)
    slope = -p
    y_ref_e = y0_e .* (x_ref_e ./ x0_e).^slope
    plot!(p_e, x_ref_e, y_ref_e,
          label="slope $slope ",
          linestyle=:dash,
          linewidth=2,
          color=colors[p])
end

# Save E-field plot
savefig(p_e, "convergence_E.png")
println("\nE-field plot saved as convergence_E.png")

# Create B-field plot
p_b = plot(
    xscale=:log10,
    yscale=:log10,
    xlabel="DoFs",
    ylabel="Relative L₂ Error",
    title="Magnetic Field (f=200 MHz)",
    legend=:topright,
    grid=true,
    minorgrid=true,
    size=(800, 600),
    dpi=300)

# Plot B-field data for each polynomial order
for p in sort(collect(keys(data_by_p)))
    p_data = data_by_p[p]
    unknowns = [d.nd_unknowns for d in p_data]
    b_errors = [d.b_rel_error for d in p_data]

    plot!(p_b, unknowns, b_errors,
          label="P=$p",
          marker=markers[p],
          markersize=6,
          linewidth=2,
          color=colors[p])

    # Add reference line with slope -p for B-field
    max_b_idx = argmax(b_errors)
    x0_b = unknowns[max_b_idx]
    y0_b = b_errors[max_b_idx]
    x_ref_b = range(x0_b, maximum(unknowns), length=100)
    slope = -p
    y_ref_b = y0_b .* (x_ref_b ./ x0_b).^slope
    plot!(p_b, x_ref_b, y_ref_b,
          label="slope $slope ",
          linestyle=:dash,
          linewidth=2,
          color=colors[p])
end

# Save B-field plot
savefig(p_b, "convergence_B.png")
println("B-field plot saved as convergence_B.png")

# Print summary table
println("\n" * "="^70)
println("Convergence Study Summary")
println("="^70)
println(@sprintf("%-4s %-6s %-12s %-15s %-15s", "p", "n", "ND Unknowns", "E-field Error", "B-field Error"))
println("-"^70)
for p in sort(collect(keys(data_by_p)))
    for d in data_by_p[p]
        println(@sprintf("%-4d %-6d %-12d %-15.6e %-15.6e", d.p, d.n, d.nd_unknowns, d.e_rel_error, d.b_rel_error))
    end
end
println("="^70)
