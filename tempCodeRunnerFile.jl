############################################################
# CertifiedSignalSequence.jl
# Generate rigorously certified signal sequences
# for the Lorenz system using validated Taylor models
# - Uses MathematicalSystems.InitialValueProblem + BlackBoxContinuousSystem
# - Prints progress every 1000 segments
# - Stops after 10 certified symbols
############################################################

using ReachabilityAnalysis
using TaylorModels
using LazySets
using IntervalArithmetic
using IntervalRootFinding
using MathematicalSystems  # for InitialValueProblem, BlackBoxContinuousSystem

# -------------------------------
# Lorenz vector field definition
# -------------------------------
@taylorize function lorenz!(dx, x, p, t)
    σ, ρ, β = 10.0, 28.0, 8/3
    dx[1] = σ * (x[2] - x[1])         # dx/dt
    dx[2] = x[1] * (ρ - x[3]) - x[2]  # dy/dt
    dx[3] = x[1] * x[2] - β * x[3]    # dz/dt
end

# -------------------------------
# Initial condition set (small box)
# -------------------------------
X0 = Hyperrectangle(
    low  = [1.0,     1.0,     20.0],
    high = [1.0001,  1.0001,  20.0001],
)

# -------------------------------
# Problem setup (no @ivp macro)
# -------------------------------
# Build a system object understood by ReachabilityAnalysis:
sys  = BlackBoxContinuousSystem(lorenz!, 3)
prob = InitialValueProblem(sys, X0)

println("Starting validated integration...")

# -------------------------------
# Validated integration
# -------------------------------
# You can tweak orders/tolerances if steps fail:
alg = TMJets21a(orderT = 18, orderQ = 1, abstol = 1e-12, adaptive = true)
sol = solve(prob, T = 3.0, alg = alg)

# -------------------------------
# Certified section crossing logic
# -------------------------------
section      = 27.0
symbols      = Char[]
max_symbols  = 10

for (segment_idx, seg) in enumerate(sol)
    if segment_idx % 1000 == 0
        println("... reached segment $segment_idx, symbols so far = $(length(symbols))")
    end

    # Overapproximate this segment as a Hyperrectangle, then extract z-interval
    Z        = overapproximate(seg, Hyperrectangle)
    lo, hi   = low(set(Z)), high(set(Z))
    zI       = IntervalArithmetic.Interval(lo[3], hi[3])

    # Does this segment straddle z = section?
    if zI.lo <= section <= zI.hi
        # Taylor model for z(t) on this segment
        tmz = set(seg)[3]

        # Validated roots of z(t) - section = 0 on the segment's local time domain
        # Note: We pass a function t -> tmz(t) and the domain of the Taylor model
        rts = roots(t -> tmz(t) - section, domain(tmz))

        for rt in rts
            t_mid = mid(dom(rt))  # midpoint of the certified root interval

            # Evaluate y(t) and dz/dt at that certified time interval midpoint
            yI  = evaluate(set(seg)[2], t_mid)
            dzI = evaluate(set(seg)[1], t_mid) *
                  evaluate(set(seg)[2], t_mid) -
                  (8/3) * evaluate(set(seg)[3], t_mid)

            # Require upward crossing: dz/dt > 0 in interval sense
            if dzI.lo > 0
                if yI.lo > 0
                    push!(symbols, 'R')
                elseif yI.hi < 0
                    push!(symbols, 'L')
                else
                    push!(symbols, '?')  # cannot certify sign of y at the crossing
                end
                println("Certified symbol $(symbols[end]) at crossing #$(length(symbols)) (segment $segment_idx)")

                if length(symbols) >= max_symbols
                    println("Reached $max_symbols certified symbols, stopping.")
                    break
                end
            end
        end
    end

    if length(symbols) >= max_symbols
        break
    end
end

println("Finished integration.")
println("Certified sequence: ", String(symbols))