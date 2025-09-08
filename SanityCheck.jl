# SanityCheck.jl
#
# Runs 7_4CertifiedSignSequence.jl several times with different
# approximation knobs (chunk size, bracket width, root locator),
# keeping seed/burn-in fixed. Compares the first 8 symbols of the
# sign sequence across variants.

using Dates
using Printf

const SCRIPT = joinpath(@__DIR__, "CertifiedSignSequence.jl")
const N_PREFIX = 8
const FIXED_SEED = 930415324  # reproducible initial point + hop

# ---------- helpers ----------

function parse_last_sequences(path::AbstractString)::Vector{String}
    lines = readlines(path)
    marker = "Sign sequences separated by newlines:"
    idx = findlast(x -> occursin(marker, x), lines)
    idx === nothing && error("Marker \"$marker\" not found in $(path).")
    seqs = String[]
    for j in (idx + 1):length(lines)
        l = strip(lines[j])
        if isempty(l) || startswith(l, "7_") || occursin("This file contains multiple runs", l)
            break
        end
        push!(seqs, l)
    end
    return seqs
end

function prefix_n(s::AbstractString, n::Int)
    io = IOBuffer()
    i = 0
    for ch in s
        i += 1
        write(io, ch)
        if i == n
            break
        end
    end
    return String(take!(io))
end

function run_one(label::String, extra_flags::Vector{String})
    outpath = joinpath(@__DIR__, "sanity_out_$(label).txt")

    # NOTE: We intentionally DO NOT set --init-range because 7_4 expects a String
    # and receives a SubString, which throws a MethodError in its parser.
    base_flags = String[
        "--seed=$(FIXED_SEED)",
        "--burn-in=5.0",
        "--num-sequences=1",
        "--eps-box=0.0001",
        "--output=$(outpath)",
    ]

    flags = vcat(base_flags, extra_flags)

    println("[$(Dates.format(now(), "HH:MM:SS"))] Running variant: $label")
    julia = Base.julia_cmd()
    cmd = `$(julia) $SCRIPT $flags`
    try
        run(cmd)  # throws on failure
    catch e
        @printf "[%s] Variant %-12s FAILED: %s\n" Dates.format(now(), "HH:MM:SS") label sprint(showerror, e)
        @printf "  Cmd: %s\n" cmd
        rethrow()
    end

    seqs = parse_last_sequences(outpath)
    isempty(seqs) && error("No sequences found in $(outpath).")
    return (label, seqs[1])
end

# ---------- main ----------

function main()
    println("Sanity check for 7_4CertifiedSignSequence.jl")
    println("Comparing first $(N_PREFIX) symbols across approximation variants.\n")

    variants = [
        ("default",      ["--chunk-len=1.0", "--target-time-width=0.05", "--root-locator=bisect_squeeze", "--global-t-max=25.0"]),
        ("coarse_tw",    ["--chunk-len=1.0", "--target-time-width=0.10", "--root-locator=bisect_squeeze", "--global-t-max=25.0"]),
        ("bigger_chunk", ["--chunk-len=2.0", "--target-time-width=0.10", "--root-locator=bisect_squeeze", "--global-t-max=25.0"]),
        ("bisect_only",  ["--chunk-len=1.0", "--target-time-width=0.05", "--root-locator=bisect",         "--global-t-max=25.0"]),
        ("strict_local", ["--chunk-len=1.0", "--target-time-width=0.02", "--root-locator=bisect_squeeze", "--global-t-max=25.0"]),
    ]

    results = Dict{String,String}()
    for (label, knobs) in variants
        try
            lab, seq = run_one(label, knobs)
            results[lab] = seq
        catch
            # already reported
        end
    end

    if isempty(results)
        println("\nNo successful runs to compare.")
        return
    end

    min_len = minimum(length.(values(results)))
    use_len = min(N_PREFIX, min_len)

    println("\nFirst $(use_len) symbols per variant:")
    println("────────────────────────────────────")
    for (label, seq) in sort(collect(results); by=first)
        @printf("%-12s : %s\n", label, prefix_n(seq, use_len))
    end

    prefixes = unique(prefix_n(seq, use_len) for seq in values(results))
    if length(prefixes) == 1
        println("\n✅ Sanity check PASSED: all variants share the same first $(use_len) symbols.")
    else
        println("\n⚠️  Sanity check WARNING: prefixes differ across variants.")
    end
end

main()
