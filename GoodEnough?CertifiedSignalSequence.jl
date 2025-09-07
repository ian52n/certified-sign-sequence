############################################################
# CertifiedSignalSequence.jl
# Generate rigorously certified symbol sequences
# for the Lorenz system using validated Taylor models
# NOTE: We overapproximate TM reach-sets as Hyperrectangles
#       (no direct poking into Taylor models).
############################################################

using ReachabilityAnalysis
using LazySets
using IntervalArithmetic
using Dates

# --------------------------------
# Lorenz vector field (taylorized)
# --------------------------------
@taylorize function lorenz!(dx, x, p, t)
    σ, ρ, β = 10.0, 28.0, 8/3
    dx[1] = σ * (x[2] - x[1])
    dx[2] = x[1] * (ρ - x[3]) - x[2]
    dx[3] = x[1] * x[2] - β * x[3]
end

# --------------------------------
# Helper: certified box bounds from a reach set that wraps a Hyperrectangle
# --------------------------------
function bounds_from_reachset(rs)
    X = set(rs)                           # should be a Hyperrectangle after overapproximate
    # low/high may be StaticArrays etc.; collect into plain Vectors and convert to Float64
    l = Float64.(collect(LazySets.low(X)))
    h = Float64.(collect(LazySets.high(X)))
    return [(l[1], h[1]),
            (l[2], h[2]),
            (l[3], h[3])]
end

# --------------------------------
# Initial condition set (small box)
# --------------------------------
X0 = Hyperrectangle(
    low  = [1.0,     1.0,     20.0],
    high = [1.0001,  1.0001,  20.0001],
)

# --------------------------------
# Problem setup
# --------------------------------
prob = @ivp(x' = lorenz!(x), dim: 3, x(0) ∈ X0)

println("Starting validated integration...")
println("Initial set: $X0")

# --------------------------------
# Algorithm parameters
# --------------------------------
alg = TMJets21a(abstol=1e-10, orderT=6, orderQ=2, maxsteps=10000)

T_max = 5.0
println("Integration time: $T_max")
println("Starting integration at $(now())...")

try
    sol_tm = solve(prob, T=T_max, alg=alg)
    println("Integration completed successfully!")
    println("Number of reach sets: $(length(sol_tm))")
    println("Time span covered: $(tspan(sol_tm))")

    if length(sol_tm) == 0
        println("❌ Empty solution: no reach sets returned.")
    else
        # Convert the Taylor-model flowpipe into Hyperrectangles for easy, certified bounds.
        sol = overapproximate(sol_tm, Hyperrectangle)

        # --------------------------------
        # Quick trajectory check (sampled)
        # --------------------------------
        println("\nQuick trajectory check:")
        sample_indices = [1, max(1, div(length(sol), 4)), max(1, div(length(sol), 2)),
                          max(1, 3*div(length(sol), 4)), length(sol)]
        sample_indices = unique(filter(i -> 1 <= i <= length(sol), sample_indices))

        z_values = Tuple{Float64,Float64}[]
        for i in sample_indices
            try
                bnds = bounds_from_reachset(sol[i])
                z_low, z_high = bnds[3]
                push!(z_values, (z_low, z_high))
                t_approx = (i-1) * T_max / length(sol)
                println("  t≈$(round(t_approx, digits=2)): z ∈ [$(round(z_low, digits=3)), $(round(z_high, digits=3))]")
            catch e
                println("  Could not extract bounds from reach set $i: $e")
            end
        end

        section = 27.0
        early_stop_reason = nothing

        if isempty(z_values)
            early_stop_reason = "No valid z-intervals were collected from sampled reach sets."
        else
            z_max_seen = maximum(last, z_values)
            z_min_seen = minimum(first, z_values)

            println("\nOverall z range sampled: [$(round(z_min_seen, digits=3)), $(round(z_max_seen, digits=3))]")

            if z_max_seen < section
                early_stop_reason = "z never reaches $section in samples — need longer integration time (try T_max=$(T_max*3))."
            elseif z_min_seen > section
                early_stop_reason = "z always above $section in samples — trajectory already passed the section (adjust X0 or analyze earlier)."
            end
        end

        symbols = Char[]
        max_symbols = 10
        crossings_checked = 0

        if early_stop_reason === nothing
            println("✓ z range includes $section — looking for certified crossings...")

            println("\nLooking for crossings of z = $section (upward only, dz/dt > 0)...")
            for (segment_idx, rs) in enumerate(sol)
                try
                    crossings_checked += 1
                    if segment_idx % 500 == 0
                        println("... checking segment $segment_idx/$(length(sol)), symbols found: $(length(symbols))")
                    end

                    bnds = bounds_from_reachset(rs)
                    x_low, x_high = bnds[1]
                    y_low, y_high = bnds[2]
                    z_low, z_high = bnds[3]

                    # Check if z-interval straddles the section
                    if z_low <= section <= z_high
                        println("Potential crossing at segment $segment_idx: z ∈ [$(round(z_low, digits=3)), $(round(z_high, digits=3))]")

                        # Certified derivative bounds via interval arithmetic: dz/dt = x*y - (8/3) * z
                        xI = IntervalArithmetic.Interval(x_low, x_high)
                        yI = IntervalArithmetic.Interval(y_low, y_high)
                        zI = IntervalArithmetic.Interval(z_low, z_high)
                        dzI = xI * yI - (8/3) * zI
                        dzdt_min = IntervalArithmetic.inf(dzI)
                        dzdt_max = IntervalArithmetic.sup(dzI)

                        println("  x ∈ [$(round(x_low, digits=3)), $(round(x_high, digits=3))]")
                        println("  y ∈ [$(round(y_low, digits=3)), $(round(y_high, digits=3))]")
                        println("  dz/dt ∈ [$(round(dzdt_min, digits=3)), $(round(dzdt_max, digits=3))]")

                        # Only consider certified upward crossings (dz/dt strictly > 0)
                        if dzdt_min > 0
                            println("  ✓ Confirmed upward crossing!")
                            if y_low > 0
                                push!(symbols, 'R')
                                println("  ✓ Symbol: R (y > 0)")
                            elseif y_high < 0
                                push!(symbols, 'L')
                                println("  ✓ Symbol: L (y < 0)")
                            else
                                push!(symbols, '?')
                                println("  ? Symbol: ? (y sign uncertain)")
                            end
                            println("🎯 Certified symbol '$(symbols[end])' at crossing #$(length(symbols))")

                            if length(symbols) >= max_symbols
                                println("Reached $max_symbols certified symbols, stopping.")
                                break
                            end
                        else
                            println("  ✗ Not a certified upward crossing (dz/dt not guaranteed > 0).")
                        end
                    end
                catch e
                    println("Warning: Could not process segment $segment_idx: $e")
                    continue
                end

                if length(symbols) >= max_symbols
                    break
                end
            end
        else
            println("\n❌ $early_stop_reason")
        end

        println("\nCrossings checked: $crossings_checked")
        println("\n" * "="^50)
        println("RESULTS")
        println("="^50)
        println("Certified sequence: $(String(symbols))")
        println("Number of symbols: $(length(symbols))")

        if isempty(symbols)
            println("\nNo crossings found (or search not performed). Possible actions:")
            println("- Increase T_max (e.g., $(T_max*2) or $(T_max*3))")
            println("- Tighten integration (smaller abstol, larger orderT/orderQ)")
            println("- Use a smaller initial set X0 to reduce over-approximation")
        else
            println("\n🎉 SUCCESS! Sequence analysis:")
            println("L count: $(count(==('L'), symbols))")
            println("R count: $(count(==('R'), symbols))")
            println("? count: $(count(==('?'), symbols))")

            if length(symbols) < max_symbols
                println("\nTo get more symbols, try:")
                println("- Increase T_max from $T_max to $(T_max * 2)")
                println("- Or increase max_symbols from $max_symbols")
            end
        end
    end

catch e
    println("ERROR during integration: $e")
    println("\nTroubleshooting suggestions:")
    println("1. Try smaller initial set (shrink X0)")
    println("2. Try different algorithm parameters (increase orderT/orderQ, tighten abstol)")
    println("3. Ensure ReachabilityAnalysis/LazySets/IntervalArithmetic versions are compatible")
    rethrow(e)
end
