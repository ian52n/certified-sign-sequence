using TaylorIntegration, IntervalArithmetic
using IntervalArithmetic: ..

# --- Helper struct to hold our results ---
mutable struct CertificationResult
    sequence::Vector{Char}
    is_successful::Bool
end

# 1. Lorenz system definition (unchanged)
function lorenz!(dx, x, p, t)
    σ, ρ, β = p
    dx[1] = σ * (x[2] - x[1])
    dx[2] = x[1] * (ρ - x[3]) - x[2]
    dx[3] = x[1] * x[2] - β * x[3]
    return nothing
end

# 2. Event condition function (unchanged)
function crossing_plane(t, x, params)
    return x[3] - 27.0 # z - 27
end

# --- Main script logic ---
function generate_certified_sequence()
    # Holds the final sequence and success status
    result = CertificationResult(Char[], true)
    
    # 3. Define the CERTIFYING EVENT ACTION
    # This function is called when the interval box `x` crosses z=27.
    function certify_and_record_crossing!(integrator)
        # Stop if we already failed or got our 5 letters
        if !result.is_successful || length(result.sequence) >= 5
            terminate!(integrator)
            return
        end

        x_interval, y_interval, z_interval = integrator.u
        σ, ρ, β = integrator.p
        dz_dt_interval = x_interval * y_interval - β * z_interval

        # CHECK 1: Is the crossing direction guaranteed upwards?
        if inf(dz_dt_interval) <= 0
            println("! CERTIFICATION FAILED: Crossing direction is ambiguous.")
            result.is_successful = false
            terminate!(integrator)
            return
        end

        # CHECK 2: Is the lobe (L/R) guaranteed?
        local symbol::Char
        if inf(y_interval) > 0
            symbol = 'L'
        elseif sup(y_interval) < 0
            symbol = 'R'
        else
            println("! CERTIFICATION FAILED: Lobe is ambiguous.")
            result.is_successful = false
            terminate!(integrator)
            return
        end

        push!(result.sequence, symbol)
        println("✓ Event at t ≈ $(inf(integrator.t)): Certified symbol -> $(symbol)")

        if length(result.sequence) >= 5
            println("--- Successfully certified 5 symbols. Stopping. ---")
            terminate!(integrator)
        end
    end

    # 4. Set up and run the integration using `taylorinteg`
    println("--- Starting Lorenz system integration for a CERTIFIED sequence ---")
    p = (10.0, 28.0, 8/3)
    uncertainty = 1e-12
    x₀ = [
    (-8.0 - uncertainty) .. (-8.0 + uncertainty),
    (-9.0 - uncertainty) .. (-9.0 + uncertainty),
    (25.0 - uncertainty) .. (25.0 + uncertainty)
]
    println("Initial uncertainty box width: $(2*uncertainty)")

    t₀ = 0.0
    T = 200.0
    order = 28
    abstol = 1e-25

    # Call the core `taylorinteg` function, bypassing the broken extension
    taylorinteg(lorenz!, x₀, t₀, T, order, abstol, p,
        event_function=crossing_plane,
        event_action=certify_and_record_crossing!,
        event_order=10) # Using a high order for accurate event location

    # 5. Display the final result
    println("\n--- Integration complete ---")
    if result.is_successful && length(result.sequence) >= 5
        println("✅ Success! Certified 5-letter sequence: ", join(result.sequence))
    else
        println("❌ Failure. Could not certify a 5-letter sequence with the initial uncertainty.")
    end
end

# Run the whole process
generate_certified_sequence()