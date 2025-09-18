#!/usr/bin/env julia
############################################################
# TestUncertifiedIntegrate.jl
#
# Simple test harness for UncertifiedIntegrate.jl
# - Generates random initial points
# - Integrates until 3 crossings or max T
# - Projects onto y-z plane
# - Calculates distance to nearest Lorenz fixed point
# - Reports min distance from outputs to nearest unstable fixed point
# - Annotates plot with runtime, seed, num_points, T_max
############################################################

using Random
using Dates
using Printf
using LinearAlgebra
using Plots

include("UncertifiedIntegrate.jl")
using .UncertifiedIntegrate

# Defaults
num_points = 50
T_max = 50.0

# -----------------------------
# Parse command-line arguments
# -----------------------------
for arg in ARGS
    if occursin("--num_points=", arg)
        global num_points = parse(Int, split(arg, "=")[2])
    elseif occursin("--T_max=", arg)
        global T_max = parse(Float64, split(arg, "=")[2])
    end
end

println("Generating $num_points points, max integration T = $T_max")

# -----------------------------
# RNG for reproducibility
# -----------------------------
seed = Random.random_seed()
rng = MersenneTwister(seed)

# -----------------------------
# Lorenz unstable fixed points (wing centers)
# -----------------------------
fixed_points = [
    (sqrt(72.0)/2, sqrt(72.0)/2, 27.0),
    (-sqrt(72.0)/2, -sqrt(72.0)/2, 27.0)
]

# -----------------------------
# Distance to nearest fixed point
# -----------------------------
distance_to_nearest_fp(u::NTuple{3,Float64}) = minimum([norm([u[1]-fx[1], u[2]-fx[2], u[3]-fx[3]]) for fx in fixed_points])

# -----------------------------
# Main program
# -----------------------------
println("=== TestUncertifiedIntegrate.jl ===")
t_start = now()
println("Start time: $t_start")

points_final = NTuple{3,Float64}[]
distances = Float64[]

for i in 1:num_points
    u0 = (rand(rng)*40.0 - 20.0, rand(rng)*40.0 - 20.0, rand(rng)*50.0)
    try
        uT = uncertified_integrate(u0, T_max; rng=rng)  # pass RNG
        push!(points_final, uT)
        push!(distances, distance_to_nearest_fp(uT))
    catch e
        @warn "Integration failed for point $u0: $e"
    end
end

min_d = isempty(distances) ? NaN : minimum(distances)

t_end = now()
# Compute runtime in milliseconds
runtime_ms = Dates.value(t_end - t_start)  # returns Int, milliseconds
h  = runtime_ms ÷ 3_600_000
r1 = runtime_ms % 3_600_000
m  = r1 ÷ 60_000
r2 = r1 % 60_000
s  = r2 ÷ 1000
ms = r2 % 1000
runtime_str = @sprintf("%d:%02d:%02d.%03d", h, m, s, ms)
println("Simulation runtime: $runtime_str")

# -----------------------------
# Plot projection onto y-z plane
# -----------------------------
ys = [p[2] for p in points_final]
zs = [p[3] for p in points_final]

scatter(ys, zs, markersize=1,
        xlabel="y", ylabel="z",
        legend=false)

# Annotate plot
annot_text = """
num_points = $num_points
T_max = $T_max
seed = $seed
runtime = $runtime_str
min distance to FP = $(round(min_d,digits=4))
"""
# Add simulation info as a text box in the upper left, slightly below top
x_pos = minimum(ys) + 0.05*(maximum(ys)-minimum(ys))
y_pos = maximum(zs) - 0.05*(maximum(zs)-minimum(zs))

# Add white rectangle background
annotate!(x_pos, y_pos, text("█████████████████████████████", :left, 12, :white))
# Add your text on top
annotate!(
    x_pos, y_pos,
    text(
        "Min distance: $(round(min_d,digits=4)).   Num points: $num_points.\nT_max: $T_max.   Runtime: $runtime_str",
        :left, 10, :brown
    )
)

fname = "UncertifiedIntegrate_yz_projection.png"
savefig(fname)
println("Plot saved → $fname")

# Auto-open plot (cross-platform)
if Sys.isapple()
    run(`open $fname`)
elseif Sys.iswindows()
    run(`start $fname`)
else
    run(`xdg-open $fname`)
end
