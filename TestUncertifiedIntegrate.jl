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
using OrdinaryDiffEq

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
distance_to_nearest_fp(u::NTuple{3,Float64}) =
    minimum([norm([u[1]-fx[1], u[2]-fx[2], u[3]-fx[3]]) for fx in fixed_points])

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
runtime_ms = Dates.value(t_end - t_start)
h  = runtime_ms ÷ 3_600_000
r1 = runtime_ms % 3_600_000
m  = r1 ÷ 60_000
r2 = r1 % 60_000
s  = r2 ÷ 1000
ms = r2 % 1000
runtime_str = @sprintf("%d:%02d:%02d.%03d", h, m, s, ms)
println("Simulation runtime: $runtime_str")

# -----------------------------
# Special trajectory from (0.01,0.01,0)
# -----------------------------
u0_special = [0.01, 0.01, 0.0]
prob_special = ODEProblem(UncertifiedIntegrate.lorenz_ode!, u0_special, (0.0, 17.5))
sol_special = solve(prob_special, Tsit5(); reltol=1e-6, abstol=1e-6)
xs_special = [u[1] for u in sol_special.u]
ys_special = [u[2] for u in sol_special.u]
xplusy_special = [xs_special[i] + ys_special[i] for i in eachindex(xs_special)]
zs_special = [u[3] for u in sol_special.u]

# -----------------------------
# Projection of uncertified points
# -----------------------------
xs = [p[1] for p in points_final]
ys = [p[2] for p in points_final]
xplusy = [xs[i] + ys[i] for i in eachindex(xs)]
zs = [p[3] for p in points_final]

# -----------------------------
# Plot
# -----------------------------
plot(xplusy_special, zs_special, color=:red, lw=1.5, alpha=0.6, label="", xlabel="x + y", ylabel="z", legend=false)
scatter!(xplusy, zs, markersize=1)


# -----------------------------
# Annotate plot
# -----------------------------
annot_text = """
num_points = $num_points
T_max = $T_max
seed = $seed
runtime = $runtime_str
min distance to FP = $(round(min_d,digits=4))
"""
x_pos = minimum(xplusy) + 0.05*(maximum(xplusy) - minimum(xplusy))
y_pos = maximum(zs) - 0.05*(maximum(zs)-minimum(zs))

annotate!(x_pos, y_pos, text("█████████████████████████████", :left, 12, :white))
annotate!(
    x_pos, y_pos,
    text(
        "Min distance: $(round(min_d,digits=4)).   Num points: $num_points.\nT_max: $T_max.   Runtime: $runtime_str",
        :left, 10, :brown
    )
)

# -----------------------------
# Save and auto-open
# -----------------------------
timestamp = Dates.format(now(), "HHMMSS")
fname = "$(num_points)_UncertifiedIntegrate_yz_projection_$(timestamp).png"
savefig(fname)
println("Plot saved → $fname")

if Sys.isapple()
    run(`open $fname`)
elseif Sys.iswindows()
    run(`start $fname`)
else
    run(`xdg-open $fname`)
end
