using TaylorIntegration, IntervalArithmetic

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

# --- Main script logic ---
function generate_certified_sequence()
    # Holds the final sequence and success status
    result = CertificationResult(Char[], true)
    
    # --- Integration parameters ---
    p = (10.0, 28.0, 8/3)
    uncertainty = 1e-12
    x₀ = [
        interval(-8.0 - uncertainty, -8.0 + uncertainty),
        interval(-9.0 - uncertainty, -9.0 + uncertainty),
        interval(25.0 - uncertainty, 25.0 + uncertainty)
    ]

    # --- Manual loop parameters ---
    t = 0.0
    step_size = 0.01
    max_steps = 50000
    
    println("--- Starting Lorenz system integration with manual event detection ---")
    println("Initial uncertainty box width: $(2*uncertainty)")

    for i in 1:max_steps
        # --- MANUAL EVENT DETECTION LOOP ---

        # 1. Store the state *before* the integration step
        x_before = x₀
        
        # 2. Integrate forward for one small step
        # We call taylorinteg without the unsupported event keywords.
        tv, xv = taylorinteg(lorenz!, x_before, t, t + step_size, 28, 1e-25, p)
        
        # 3. Get the state *after* the integration step
        x_after = xv[end]
        t = tv[end]

        # 4. Check if the entire interval box has crossed the z=27 plane upwards
        z_before = x_before[3]
        z_after = x_after[3]
        
        # This is a guaranteed up-crossing if the highest point before the step
        # is below 27 and the lowest point after the step is above 27.
        if sup(z_before) < 27 && inf(z_after) > 27
            println("Plane crossing detected around t ≈ $t")
            
            # --- Certification Checks ---
            x_interval, y_interval, z_interval = x_after
            dz_dt_interval = x_interval * y_interval - p[3] * z_interval

            if inf(dz_dt_interval) <= 0
                println("! CERTIFICATION FAILED: Crossing direction is ambiguous.")
                result.is_successful = false
                break # Exit the loop
            end

            local symbol::Char
            if inf(y_interval) > 0
                symbol = 'L'
            elseif sup(y_interval) < 0
                symbol = 'R'
            else
                println("! CERTIFICATION FAILED: Lobe is ambiguous.")
                result.is_successful = false
                break # Exit the loop
            end

            push!(result.sequence, symbol)
            println("✓ Certified symbol -> $(symbol)")
            
            # Check if we are done
            if length(result.sequence) >= 5
                println("--- Successfully certified 5 symbols. ---")
                break # Exit the loop
            end
        end
        
        # 5. Prepare for the next iteration
        x₀ = x_after
    end

    # --- Final Result Display ---
    println("\n--- Integration complete ---")
    if result.is_successful && length(result.sequence) >= 5
        println("✅ Success! Certified 5-letter sequence: ", join(result.sequence))
    else
        println("❌ Failure. Could not certify a 5-letter sequence.")
        println("   Generated sequence: ", join(result.sequence))
    end
end

# Run the whole process
generate_certified_sequence()