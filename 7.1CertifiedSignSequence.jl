############################################################
# 7_1CertifiedSignSequence.jl — MULTI-SEQUENCE (v7.1)
#
# What this version does:
# • Starts from a random point (or user-provided range), performs
#   a fast, *uncertified* burn-in with OrdinaryDiffEq to get near
#   the Lorenz attractor.
# • Runs CHUNKED, *rigorous* certified integration to produce a
#   sign sequence (L/R/?) stopping at the FIRST uncertifiable
#   crossing candidate.
# • Takes a point estimate at the failing slab’s center, performs
#   a random *uncertified* hop, and then generates the next
#   certified sign sequence; repeats until the requested number of
#   sequences is produced.
# • Writes a .txt report in the format the user requested.
#
# Reproducibility change (this edit):
# • All randomness now comes from a *dedicated* RNG instance
#   (MersenneTwister seeded with `--seed`). No use of the global
#   RNG. This makes initial point *and* hop times reproducible.
#
# Key Fixes carried over:
# • Typed local for `last_below_box::Union{Nothing,Hyperrectangle} = nothing`
#   (avoids Union-constructor error).
# • Robust CLI parsing for --init-range and numeric knobs.
############################################################

using ReachabilityAnalysis
using LazySets
using IntervalArithmetic
using OrdinaryDiffEq
using Dates
using Random
using Printf

const IA = IntervalArithmetic
const VERSION_STR = "v7.1"

# ============================
# Lorenz dynamics
# ============================
@taylorize function lorenz!(dx, x, p, t)
    σ, ρ, β = 10.0, 28.0, 8/3
    dx[1] = σ * (x[2] - x[1])
    dx[2] = x[1] * (ρ - x[3]) - x[2]
    dx[3] = x[1] * x[2] - β * x[3]
end

# OrdinaryDiffEq version for uncertified steps
function lorenz_ode!(du, u, p, t)
    σ, ρ, β = 10.0, 28.0, 8/3
    du[1] = σ * (u[2] - u[1])
    du[2] = u[1] * (ρ - u[3]) - u[2]
    du[3] = u[1] * u[2] - β * u[3]
end

# ============================
# Configuration / CLI
# ============================
Base.@kwdef mutable struct Config
    seed::Int
    burn_in::Float64                    = 5.0
    num_sequences::Int                 = 2
    output_path::String                = "CertifiedSignSequence_report.txt"
    # Certified scan knobs
    section::Float64                   = 27.0
    global_t_max::Float64              = 25.0
    chunk_len::Float64                 = 1.0
    target_time_width::Float64         = 0.05
    # uncertified random hop time range
    hop_min::Float64                   = 1.0
    hop_max::Float64                   = 3.0
    # initial random sampling range (if not provided via --init-range)
    init_range::NTuple{6,Float64}      = (-15.0, 15.0, -20.0, 20.0, 10.0, 40.0)
    eps_box::Float64                   = 1e-4  # half-width for initial certified box
    # Global flowpipe algorithm (coarser for speed)
    global_alg                         = TMJets21a(abstol=1e-8, orderT=4, orderQ=1, maxsteps=25000)
    # Local bracketing algorithm (still rigorous)
    local_alg                          = TMJets21a(abstol=1e-11, orderT=6, orderQ=2, maxsteps=120000)
end

function parse_init_range(s::String)::NTuple{6,Float64}
    parts = split(s, ",")
    length(parts) == 6 || error("--init-range must have 6 comma-separated numbers")
    vals = parse.(Float64, strip.(parts))
    return (vals[1], vals[2], vals[3], vals[4], vals[5], vals[6])
end

function parse_kv!(cfg::Config, arg::AbstractString)
    occursin("=", arg) || return cfg
    k, v = split(arg, "=", limit=2)
    k = strip(k); v = strip(v)
    if     k == "--seed"                cfg.seed = parse(Int, v)
    elseif k == "--burn-in"             cfg.burn_in = parse(Float64, v)
    elseif k == "--num-sequences"       cfg.num_sequences = parse(Int, v)
    elseif k == "--output"              cfg.output_path = v
    elseif k == "--target-time-width"   cfg.target_time_width = parse(Float64, v)
    elseif k == "--global-t-max"        cfg.global_t_max = parse(Float64, v)
    elseif k == "--chunk-len"           cfg.chunk_len = parse(Float64, v)
    elseif k == "--hop-min"             cfg.hop_min = parse(Float64, v)
    elseif k == "--hop-max"             cfg.hop_max = parse(Float64, v)
    elseif k == "--eps-box"             cfg.eps_box = parse(Float64, v)
    elseif k == "--init-range"
        cfg.init_range = parse_init_range(v)
    end
    return cfg
end

function get_config()::Config
    seed_default = Int(mod(Dates.now().instant.periods.value, 10^9))
    cfg = Config(seed=seed_default)
    for a in ARGS
        parse_kv!(cfg, a)
    end
    cfg.hop_max ≥ cfg.hop_min || error("--hop-max must be ≥ --hop-min")
    return cfg
end

# ============================
# Helpers for sets/intervals
# ============================
box_of(S) = S isa Hyperrectangle ? S : overapproximate(S, Hyperrectangle)

function bounds_from_reachset(rs)
    X = set(rs) |> box_of
    l = Float64.(collect(LazySets.low(X)))
    h = Float64.(collect(LazySets.high(X)))
    return [(l[1], h[1]), (l[2], h[2]), (l[3], h[3])]
end

function xyz_intervals_from_box(B::Hyperrectangle)
    l = LazySets.low(B); h = LazySets.high(B)
    xI = IA.Interval(Float64(l[1]), Float64(h[1]))
    yI = IA.Interval(Float64(l[2]), Float64(h[2]))
    zI = IA.Interval(Float64(l[3]), Float64(h[3]))
    return xI, yI, zI
end

dz_interval(xI, yI, zI) = xI*yI - (8/3)*zI
hull3(a::IA.Interval, b::IA.Interval, c::IA.Interval) = IA.hull(IA.hull(a, b), c)

# Normalize any tspan-like object to (tlo, thi)
function tspan_tuple(obj)
    ts = tspan(obj)
    if ts isa Tuple
        return (float(ts[1]), float(ts[2]))
    elseif ts isa IA.Interval
        return (float(IA.inf(ts)), float(IA.sup(ts)))
    else
        try
            return (float(first(ts)), float(last(ts)))
        catch
            error("Unsupported tspan type: $(typeof(ts))")
        end
    end
end

# Midpoint of a hyperrectangle
function center_point(B::Hyperrectangle)::NTuple{3,Float64}
    l = LazySets.low(B); h = LazySets.high(B)
    c = ((l[1]+h[1])/2, (l[2]+h[2])/2, (l[3]+h[3])/2)
    return (Float64(c[1]), Float64(c[2]), Float64(c[3]))
end

# Build small box around a point
function box_from_point(u::NTuple{3,Float64}, eps::Float64)::Hyperrectangle
    (x,y,z) = u
    return Hyperrectangle(low=[x-eps, y-eps, z-eps], high=[x+eps, y+eps, z+eps])
end

# ============================
# Local certified bracketing
# ============================
function certify_upward_crossing!(
    lorenz!; B_below::Hyperrectangle, section::Float64,
    t_start_abs::Float64,
    T_guess::Float64 = 0.05,
    T_max_local::Float64 = 0.30,
    target_time_width::Float64,
    max_expand::Int = 8,
    max_bisect::Int = 22,
    alg_refined = TMJets21a(abstol=1e-11, orderT=6, orderQ=2, maxsteps=120000)
)::Tuple{Bool,Float64,Float64,IA.Interval{Float64},IA.Interval{Float64},Float64,Float64}

    _, _, zI0 = xyz_intervals_from_box(B_below)
    if IA.sup(zI0) >= section
        return (false, 0.0, 0.0, IA.Interval(0,0), IA.Interval(0,0), 0.0, 0.0)
    end

    # 1) Expand-forward until strictly above section
    T = T_guess
    T_total = 0.0
    B_lo = B_below
    B_hi = B_below
    found = false

    for _ in 1:max_expand
        prob_loc = @ivp(x' = lorenz!(x), dim: 3, x(0) ∈ B_lo)
        sol_loc_tm = solve(prob_loc, T=T, alg=alg_refined)
        sol_loc = overapproximate(sol_loc_tm, Hyperrectangle)
        B_end = set(sol_loc[end]) |> box_of

        _, _, zI_end = xyz_intervals_from_box(B_end)
        if IA.inf(zI_end) > section
            B_hi = B_end
            found = true
            break
        else
            B_lo = B_end
            T_total += T
            if T_total >= T_max_local
                break
            end
            T = min(T*1.7, T_max_local - T_total)
            if T <= 0
                break
            end
        end
    end
    if !found
        return (false, 0.0, 0.0, IA.Interval(0,0), IA.Interval(0,0), 0.0, 0.0)
    end

    # 2) Bisection with direction check on hull{lo,mid,hi}
    τ_lo = T_total
    τ_hi = T_total + T
    last_dz_min = -Inf
    last_dz_max =  Inf

    mid_cache = Dict{Float64,Hyperrectangle}()

    for _ in 1:max_bisect
        τ_mid = (τ_lo + τ_hi) / 2

        B_mid = get!(mid_cache, τ_mid) do
            prob_mid = @ivp(x' = lorenz!(x), dim: 3, x(0) ∈ B_below)
            sol_mid_tm = solve(prob_mid, T=τ_mid, alg=alg_refined)
            sol_mid = overapproximate(sol_mid_tm, Hyperrectangle)
            set(sol_mid[end]) |> box_of
        end

        xI_lo, yI_lo, zI_lo = xyz_intervals_from_box(B_lo)
        xI_mid, yI_mid, zI_mid = xyz_intervals_from_box(B_mid)
        xI_hi, yI_hi, zI_hi = xyz_intervals_from_box(B_hi)

        xI_all = hull3(xI_lo, xI_mid, xI_hi)
        yI_all = hull3(yI_lo, yI_mid, yI_hi)
        zI_all = hull3(zI_lo, zI_mid, zI_hi)
        dzI_all = dz_interval(xI_all, yI_all, zI_all)
        last_dz_min, last_dz_max = IA.inf(dzI_all), IA.sup(dzI_all)

        if IA.inf(zI_mid) > section
            B_hi = B_mid
            τ_hi = τ_mid
        else
            B_lo = B_mid
            τ_lo = τ_mid
        end

        if (τ_hi - τ_lo) ≤ target_time_width && last_dz_min > 0
            return (true, t_start_abs + τ_lo, t_start_abs + τ_hi, yI_lo, yI_hi, last_dz_min, last_dz_max)
        end
    end

    # Final acceptance if positive and small enough
    xI_lo, yI_lo, zI_lo = xyz_intervals_from_box(B_lo)
    xI_hi, yI_hi, zI_hi = xyz_intervals_from_box(B_hi)
    dzI_all = dz_interval(hull3(xI_lo,xI_hi,xI_lo), hull3(yI_lo,yI_hi,yI_lo), hull3(zI_lo,zI_hi,zI_lo))
    last_dz_min, last_dz_max = IA.inf(dzI_all), IA.sup(dzI_all)

    if last_dz_min > 0 && (τ_hi - τ_lo) ≤ target_time_width
        return (true, t_start_abs + τ_lo, t_start_abs + τ_hi, yI_lo, yI_hi, last_dz_min, last_dz_max)
    else
        return (false, 0.0, 0.0, IA.Interval(0,0), IA.Interval(0,0), last_dz_min, last_dz_max)
    end
end

# ============================
# CHUNKED certified sign sequence
# ============================
Base.@kwdef mutable struct SequenceResult
    symbols::String
    brackets::Vector{Tuple{Float64,Float64}}
    stopped_uncertified::Bool
    stop_tspan::Tuple{Float64,Float64}
    stop_box::Union{Nothing,Hyperrectangle}
    total_time_scanned::Float64
end

function run_certified_sequence(X_start::Hyperrectangle, cfg::Config)::SequenceResult
    println("\n===== Generating certified sign sequence =====")
    symbols = Char[]
    hit_brackets = Vector{Tuple{Float64,Float64}}()
    awaiting_new_cross = true
    last_below_box::Union{Nothing,Hyperrectangle} = nothing  # TYPED LOCAL (fix)
    last_below_tend_abs = 0.0

    stopped_uncertified = false
    stop_tlo_abs = 0.0
    stop_thi_abs = 0.0
    stop_box::Union{Nothing,Hyperrectangle} = nothing

    t_offset = 0.0
    chunk_id = 0
    X_current = X_start

    while t_offset < cfg.global_t_max && !stopped_uncertified
        chunk_id += 1
        T_chunk = min(cfg.chunk_len, cfg.global_t_max - t_offset)

        prob_chunk = @ivp(x' = lorenz!(x), dim: 3, x(0) ∈ X_current)
        sol_tm_chunk = solve(prob_chunk, T=T_chunk, alg=cfg.global_alg)
        sol_chunk = overapproximate(sol_tm_chunk, Hyperrectangle)

        # Actual end time (could be < T_chunk if step cap)
        t0_rel, t1_rel = tspan_tuple(sol_chunk)
        t0_abs = t_offset + t0_rel
        t1_abs = t_offset + t1_rel

        println("— Integrated chunk $(chunk_id): t ∈ [$(round(t0_abs,digits=3)), $(round(t1_abs,digits=3))], sets: $(length(sol_chunk))")

        for (i, rs) in enumerate(sol_chunk)
            bnds = bounds_from_reachset(rs)
            z_low, z_high = bnds[3]
            tlo_rel, thi_rel = tspan_tuple(rs)
            tlo_abs = t_offset + tlo_rel
            thi_abs = t_offset + thi_rel

            # Track when fully below (to arm and provide B_below)
            if z_high < cfg.section
                last_below_box = set(rs) |> box_of
                last_below_tend_abs = thi_abs
                if !awaiting_new_cross
                    awaiting_new_cross = true
                end
                continue
            end

            # Consider straddling slabs only when armed
            if awaiting_new_cross && (z_low <= cfg.section <= z_high)
                if last_below_box === nothing
                    continue
                end

                ok, t_lo_abs, t_hi_abs, yI_lo, yI_hi, dz_min, dz_max = certify_upward_crossing!(
                    lorenz!;
                    B_below=last_below_box, section=cfg.section,
                    t_start_abs=last_below_tend_abs,
                    target_time_width=cfg.target_time_width,
                    alg_refined=cfg.local_alg
                )

                if ok
                    yI_hull = IA.hull(yI_lo, yI_hi)
                    if IA.inf(yI_hull) > 0
                        push!(symbols, 'R')
                        println("   ✓ Certified ABS bracket [$(round(t_lo_abs,digits=6)), $(round(t_hi_abs,digits=6))], width=$(round(t_hi_abs-t_lo_abs,digits=6)) ⇒ R")
                    elseif IA.sup(yI_hull) < 0
                        push!(symbols, 'L')
                        println("   ✓ Certified ABS bracket [$(round(t_lo_abs,digits=6)), $(round(t_hi_abs,digits=6))], width=$(round(t_hi_abs-t_lo_abs,digits=6)) ⇒ L")
                    else
                        push!(symbols, '?')
                        println("   ✓ Certified ABS bracket [$(round(t_lo_abs,digits=6)), $(round(t_hi_abs,digits=6))], width=$(round(t_hi_abs-t_lo_abs,digits=6)) ⇒ ?")
                    end
                    push!(hit_brackets, (t_lo_abs, t_hi_abs))
                    awaiting_new_cross = false
                else
                    println("   ⛔ Could NOT certify this straddling slab. Stopping sequence.")
                    stopped_uncertified = true
                    stop_tlo_abs = tlo_abs
                    stop_thi_abs = thi_abs
                    stop_box = set(rs) |> box_of
                    break
                end
            end
        end

        # Early stop?
        if stopped_uncertified
            break
        end

        # Next chunk starts from terminal set of this chunk
        X_current = set(sol_chunk[end]) |> box_of
        t_offset = t1_abs
    end

    return SequenceResult(
        symbols = String(symbols),
        brackets = hit_brackets,
        stopped_uncertified = stopped_uncertified,
        stop_tspan = (stop_tlo_abs, stop_thi_abs),
        stop_box = stop_box,
        total_time_scanned = t_offset
    )
end

# ============================
# Uncertified steps (burn-in / hop)
# ============================
function uncertified_integrate(u0::NTuple{3,Float64}, T::Float64; reltol=1e-6, abstol=1e-6)
    u = collect(u0)
    prob = ODEProblem(lorenz_ode!, u, (0.0, T))
    sol = solve(prob, Tsit5(); reltol=reltol, abstol=abstol)
    uT = sol.u[end]
    return (Float64(uT[1]), Float64(uT[2]), Float64(uT[3]))
end

# ============================
# Reporting (requested format)
# ============================
# Format: HH:MM:SS.mmm
function fmt_runtime(ms::Int)
    h  = ms ÷ 3_600_000
    r1 = ms % 3_600_000
    m  = r1 ÷ 60_000
    r2 = r1 % 60_000
    s  = r2 ÷ 1_000
    mm = r2 % 1_000
    return @sprintf("%d:%02d:%02d.%03d", h, m, s, mm)
end

function write_report(cfg::Config;
    seed::Int,
    sequences::Vector{String},
    init_pt::NTuple{3,Float64},
    start_points::Vector{NTuple{3,Float64}},
    runtime_ms::Int
)
    open(cfg.output_path, "w") do io
        # Exact header requested by the user:
        println(io, "7_1CertifiedSignSequence.jl")
        println(io, "Program Run Time: ", fmt_runtime(runtime_ms))
        println(io, "Random seed: ", seed)
        println(io, "")
        println(io, "Parameters / Knobs")
        println(io, "  section              = $(cfg.section)")
        println(io, "  global_t_max         = $(cfg.global_t_max)")
        println(io, "  chunk_len            = $(cfg.chunk_len)")
        println(io, "  target_time_width    = $(cfg.target_time_width)")
        println(io, "  burn_in              = $(cfg.burn_in)")
        println(io, "  hop_min..hop_max     = $(cfg.hop_min)..$(cfg.hop_max)")
        println(io, "  eps_box              = $(cfg.eps_box)")
        println(io, "  init_range           = $(cfg.init_range)")
        println(io, "")
        println(io, "Initial random point   = $(init_pt)")

        # Post burn-in / hop points per sequence (requested wording)
        for (k, pt) in enumerate(start_points)
            if k == 1
                # Match sample: "Post burn-in point 1    = (...)"
                println(io, "Post burn-in point $(k)    = $(pt)")
            else
                # Match sample: "Post-burn in point 2 = ...."
                println(io, "Post-burn in point $(k) = $(pt)")
            end
        end

        println(io, "")
        println(io, "Sign sequences separated by newlines:")
        for s in sequences
            println(io, s)
        end
    end
    println("Wrote report to: $(cfg.output_path)")
end

# ============================
# Main
# ============================
function main()
    cfg = get_config()
    t0_ns = time_ns()
    println("Starting run at $(now()). Version $(VERSION_STR)")
    println("Requested sequences: $(cfg.num_sequences)  |  burn-in: $(cfg.burn_in) s  |  output: $(cfg.output_path)")

    # Dedicated RNG for *all* random draws → reproducible results
    rng = MersenneTwister(cfg.seed)

    # Pick initial point with dedicated RNG
    (xl, xh, yl, yh, zl, zh) = cfg.init_range
    u0 = (rand(rng)*(xh-xl)+xl, rand(rng)*(yh-yl)+yl, rand(rng)*(zh-zl)+zl)
    println("Random initial point: $(u0)")

    # Burn-in (uncertified) to get near attractor (deterministic given u0, T, solver/tols)
    u_burn = uncertified_integrate(u0, cfg.burn_in)
    println("Post burn-in point: $(u_burn)")

    sequences = String[]
    start_points = NTuple{3,Float64}[]
    push!(start_points, u_burn)

    X_start = box_from_point(u_burn, cfg.eps_box)

    for seq_idx in 1:cfg.num_sequences
        println("\n===== Generating certified sign sequence $(seq_idx) =====")
        res = run_certified_sequence(X_start, cfg)
        push!(sequences, res.symbols)

        if seq_idx == cfg.num_sequences
            break
        end

        # Determine next starting point:
        # If we stopped due to uncertifiable slab, take its center; otherwise use current X_start center.
        u_est = if res.stopped_uncertified && !(res.stop_box === nothing)
            center_point(res.stop_box)
        else
            center_point(X_start)
        end
        # Random uncertified hop to decorrelate — draw from dedicated RNG
        hop_T = rand(rng) * (cfg.hop_max - cfg.hop_min) + cfg.hop_min
        u_next = uncertified_integrate(u_est, hop_T)
        println("Post-hop point: $(u_next)")

        push!(start_points, u_next)
        X_start = box_from_point(u_next, cfg.eps_box)
    end

    runtime_ms = Int(round((time_ns() - t0_ns) / 1e6))
    write_report(cfg;
        seed = cfg.seed,
        sequences = sequences,
        init_pt = u0,
        start_points = start_points,
        runtime_ms = runtime_ms
    )
end

# Entry
main()
