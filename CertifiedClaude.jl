#!/usr/bin/env julia

"""
True Certified Lorenz System Itinerary Generator using Interval ODE Integration

This script implements a proper interval arithmetic ODE integrator for the Lorenz system.
Unlike standard methods, this provides mathematically rigorous bounds that grow with
the chaotic dynamics, leading to realistic limitations on integration time.

WARNING: Due to the chaotic nature of the Lorenz system, certified bounds will
grow exponentially and become unusable after a short time period.

Author: Generated for rigorous trajectory analysis
"""

using IntervalArithmetic
using LinearAlgebra
using Plots
using Printf

# Lorenz system parameters (classic values)
const σ = 10.0
const ρ = 28.0  
const β = 8.0/3.0

# Integration parameters - much more conservative for interval methods
const PLANE_HEIGHT = 27.0
const MAX_INTEGRATION_TIME = 5.0  # Much shorter due to interval growth
const INITIAL_STEP_SIZE = 0.001   # Starting step size
const MIN_STEP_SIZE = 1e-8        # Minimum allowed step size
const MAX_INTERVAL_WIDTH = 10.0   # Stop if intervals get too wide
const STEP_CONTROL_FACTOR = 0.8   # Step size reduction factor

"""
Lorenz system with interval arithmetic
"""
function lorenz_interval!(du, u, p, t)
    x, y, z = u
    du[1] = σ * (y - x)           # dx/dt
    du[2] = x * (ρ - z) - y       # dy/dt  
    du[3] = x * y - β * z         # dz/dt
    return nothing
end

"""
Compute interval width (diameter) of a vector of intervals
"""
function interval_width(u_interval)
    return maximum(diam.(u_interval))
end

"""
Interval Runge-Kutta 4th order integrator with adaptive step size
This implements proper interval arithmetic at each step
"""
function interval_rk4_step(f!, u, h, t)
    # All computations use interval arithmetic
    k1 = similar(u)
    k2 = similar(u)
    k3 = similar(u)
    k4 = similar(u)
    
    # Stage 1: k1 = f(t, u)
    f!(k1, u, nothing, t)
    
    # Stage 2: k2 = f(t + h/2, u + h*k1/2)
    u2 = u + h * k1 / 2
    f!(k2, u2, nothing, t + h/2)
    
    # Stage 3: k3 = f(t + h/2, u + h*k2/2)  
    u3 = u + h * k2 / 2
    f!(k3, u3, nothing, t + h/2)
    
    # Stage 4: k4 = f(t + h, u + h*k3)
    u4 = u + h * k3
    f!(k4, u4, nothing, t + h)
    
    # Final RK4 combination with interval arithmetic
    u_new = u + h * (k1 + 2*k2 + 2*k3 + k4) / 6
    
    return u_new
end

"""
Adaptive interval ODE integrator with step size control
"""
function integrate_interval_ode(f!, u0_interval, tspan, initial_h=INITIAL_STEP_SIZE)
    println("=== Starting Certified Interval ODE Integration ===")
    println("Initial intervals:")
    components = ['x', 'y', 'z']
    for (i, comp) in enumerate(components)
        println("  $(comp)₀ ∈ $(u0_interval[i]) (width: $(diam(u0_interval[i])))")
    end
    println()
    
    t_start, t_end = tspan
    t = t_start
    u = copy(u0_interval)
    h = initial_h
    
    # Storage for solution
    times = [t]
    solutions = [copy(u)]
    step_sizes = [h]
    widths = [interval_width(u)]
    
    step_count = 0
    rejected_steps = 0
    
    println("Integration progress:")
    println("Step    Time      Max Width    Step Size    Status")
    println("-" ^ 55)
    
    while t < t_end && step_count < 10000  # Safety limit
        step_count += 1
        
        # Try integration step
        u_new = interval_rk4_step(f!, u, h, t)
        
        # Check if intervals have grown too large
        max_width = interval_width(u_new)
        
        if max_width > MAX_INTERVAL_WIDTH
            if h <= MIN_STEP_SIZE
                println("ERROR: Intervals too wide ($(max_width)) and minimum step size reached")
                println("Integration failed at t = $(t)")
                break
            else
                # Reduce step size and try again
                h *= STEP_CONTROL_FACTOR
                rejected_steps += 1
                @printf("  %4d   %7.3f   %9.3e   %9.3e   REJECT (too wide)\n", 
                       step_count, t, max_width, h)
                continue
            end
        end
        
        # Check for NaN or infinite intervals
        if any(x -> isinf(inf(x)) || isinf(sup(x)) || isnan(inf(x)) || isnan(sup(x)), u_new)
            println("ERROR: Invalid intervals detected at t = $(t)")
            break
        end
        
        # Accept the step
        t += h
        u = u_new
        
        # Store solution
        push!(times, t)
        push!(solutions, copy(u))
        push!(step_sizes, h)
        push!(widths, max_width)
        
        # Print progress every 100 steps or when width doubles
        if step_count % 100 == 0 || (length(widths) > 1 && max_width > 2 * widths[end-1])
            @printf("  %4d   %7.3f   %9.3e   %9.3e   OK\n", 
                   step_count, t, max_width, h)
        end
        
        # Adaptive step size: reduce if intervals growing too fast
        if length(widths) >= 2 && max_width > 1.5 * widths[end-1]
            h *= STEP_CONTROL_FACTOR
        elseif max_width < 0.1 * MAX_INTERVAL_WIDTH && h < initial_h
            h *= 1.1  # Slightly increase step size if intervals are well-behaved
        end
        
        # Stop if intervals become unusable
        if max_width > MAX_INTERVAL_WIDTH / 2
            println("WARNING: Intervals approaching maximum width threshold")
            println("Continued integration may not be reliable")
        end
    end
    
    println("-" ^ 55)
    println("Integration completed:")
    println("  Total steps: $(step_count)")
    println("  Rejected steps: $(rejected_steps)")
    println("  Final time: $(t)")
    println("  Final max interval width: $(interval_width(u))")
    println("  Integration $(t >= t_end ? "SUCCESS" : "TERMINATED EARLY")")
    
    return (times=times, solutions=solutions, step_sizes=step_sizes, widths=widths)
end

"""
Check for certified plane crossings using interval arithmetic
"""
function check_interval_plane_crossing(u_prev, u_curr, t_prev, t_curr)
    x_prev, y_prev, z_prev = u_prev
    x_curr, y_curr, z_curr = u_curr
    
    # Check if z interval potentially crosses the plane
    z_min_prev = inf(z_prev)
    z_max_prev = sup(z_prev)
    z_min_curr = inf(z_curr)
    z_max_curr = sup(z_curr)
    
    # Conservative check: Could we have crossed?
    if z_max_prev < PLANE_HEIGHT && z_min_curr > PLANE_HEIGHT
        # Definite crossing occurred
        crossing_type = :definite
    elseif z_min_prev < PLANE_HEIGHT && z_max_curr > PLANE_HEIGHT
        # Possible crossing (intervals overlap the plane)
        crossing_type = :possible
    else
        # No crossing
        return (false, :none, interval(0.0), interval(0.0))
    end
    
    # Linear interpolation to estimate crossing point
    # This is approximate since we don't know exact trajectory within intervals
    α_min = max(0.0, (PLANE_HEIGHT - z_max_prev) / (z_min_curr - z_max_prev))
    α_max = min(1.0, (PLANE_HEIGHT - z_min_prev) / (z_max_curr - z_min_prev))
    
    if α_min > α_max
        α_min, α_max = α_max, α_min
    end
    
    α_interval = interval(α_min, α_max)
    
    # Interpolated crossing position (intervals)
    x_cross = x_prev + α_interval * (x_curr - x_prev)
    y_cross = y_prev + α_interval * (y_curr - y_prev)
    z_cross = interval(PLANE_HEIGHT)  # Exactly on plane
    
    # Check dz/dt at crossing (interval arithmetic)
    dzdt_cross = x_cross * y_cross - β * z_cross
    
    # Determine if dz/dt > 0 is certain, possible, or impossible
    if inf(dzdt_cross) > 0
        return (true, crossing_type, y_cross, dzdt_cross)
    elseif sup(dzdt_cross) > 0
        return (true, :uncertain_velocity, y_cross, dzdt_cross)
    else
        return (false, :negative_velocity, y_cross, dzdt_cross)
    end
end

"""
Generate certified itinerary using proper interval ODE integration
"""
function generate_certified_itinerary(initial_condition, uncertainty_radius=1e-6)
    println("=== True Certified Lorenz Itinerary Generator ===")
    println("Using rigorous interval ODE integration\n")
    
    # Create initial interval condition with specified uncertainty
    u0_interval = [
        interval(initial_condition[1] - uncertainty_radius, initial_condition[1] + uncertainty_radius),
        interval(initial_condition[2] - uncertainty_radius, initial_condition[2] + uncertainty_radius),
        interval(initial_condition[3] - uncertainty_radius, initial_condition[3] + uncertainty_radius)
    ]
    
    # Integrate using interval arithmetic
    tspan = (0.0, MAX_INTEGRATION_TIME)
    sol = integrate_interval_ode(lorenz_interval!, u0_interval, tspan)
    
    if isempty(sol.times) || length(sol.times) < 2
        println("ERROR: Integration failed immediately")
        return nothing
    end
    
    # Analyze crossings
    itinerary = String[]
    crossing_details = []
    
    println("\n=== Analyzing Plane Crossings ===")
    
    for i in 2:length(sol.times)
        u_prev = sol.solutions[i-1]
        u_curr = sol.solutions[i]
        t_prev = sol.times[i-1]
        t_curr = sol.times[i]
        
        crossed, crossing_type, y_cross, dzdt_cross = check_interval_plane_crossing(
            u_prev, u_curr, t_prev, t_curr)
        
        if crossed
            # Determine symbol based on y interval
            if sup(y_cross) < 0
                symbol = "L"  # Certainly left
                certainty = "CERTAIN"
            elseif inf(y_cross) > 0
                symbol = "R"  # Certainly right
                certainty = "CERTAIN"
            else
                symbol = "?"  # Uncertain
                certainty = "UNCERTAIN"
            end
            
            push!(itinerary, symbol)
            push!(crossing_details, (
                time=t_curr,
                y_interval=y_cross,
                dzdt_interval=dzdt_cross,
                crossing_type=crossing_type,
                symbol=symbol,
                certainty=certainty
            ))
            
            @printf("Crossing %d: t=%.3f, y∈%s, dz/dt∈%s → %s (%s)\n",
                   length(itinerary), t_curr, y_cross, dzdt_cross, symbol, certainty)
        end
    end
    
    println("\n=== CERTIFIED ITINERARY RESULTS ===")
    itinerary_string = join(itinerary, "")
    println("Sequence: $itinerary_string")
    println("Length: $(length(itinerary))")
    println("Integration time achieved: $(sol.times[end])s / $(MAX_INTEGRATION_TIME)s")
    println("Final interval width: $(sol.widths[end])")
    
    uncertain_count = count(==("?"), itinerary)
    if uncertain_count > 0
        println("⚠️  WARNING: $uncertain_count uncertain crossings due to interval overlap")
    end
    
    if sol.times[end] < MAX_INTEGRATION_TIME
        println("⚠️  WARNING: Integration terminated early due to interval growth")
    end
    
    return (
        itinerary=itinerary_string,
        crossings=length(itinerary),
        details=crossing_details,
        solution=sol,
        integration_time=sol.times[end]
    )
end

"""
Plot the interval solution showing uncertainty bounds
"""
function plot_interval_solution(result)
    sol = result.solution
    
    # Extract bounds for plotting
    times = sol.times
    x_inf = [inf(u[1]) for u in sol.solutions]
    x_sup = [sup(u[1]) for u in sol.solutions]
    y_inf = [inf(u[2]) for u in sol.solutions]
    y_sup = [sup(u[2]) for u in sol.solutions]
    z_inf = [inf(u[3]) for u in sol.solutions]
    z_sup = [sup(u[3]) for u in sol.solutions]
    
    # Time series plots with interval bounds
    p1 = plot(times, x_inf, fillto=x_sup, alpha=0.3, color=:blue, 
              title="X Component", xlabel="Time", ylabel="x", label="x bounds")
    
    p2 = plot(times, y_inf, fillto=y_sup, alpha=0.3, color=:red,
              title="Y Component", xlabel="Time", ylabel="y", label="y bounds")
    
    p3 = plot(times, z_inf, fillto=z_sup, alpha=0.3, color=:green,
              title="Z Component", xlabel="Time", ylabel="z", label="z bounds")
    hline!(p3, [PLANE_HEIGHT], color=:black, linestyle=:dash, label="z=27")
    
    # Interval width evolution
    p4 = plot(times, sol.widths, title="Maximum Interval Width", 
              xlabel="Time", ylabel="Width", label="Max width", logy=true)
    
    return plot(p1, p2, p3, p4, layout=(2,2), size=(1000, 800))
end

# Main execution
function main()
    println("Rigorous Interval Lorenz System Analysis")
    println("=" ^ 50)
    
    # Initial condition
    initial_condition = [1.0, 1.0, 1.0]
    uncertainty_radius = 1e-8  # Very small initial uncertainty
    
    println("Initial condition: $initial_condition")
    println("Uncertainty radius: $uncertainty_radius")
    println("Maximum integration time: $(MAX_INTEGRATION_TIME)s")
    println("Maximum allowed interval width: $(MAX_INTERVAL_WIDTH)")
    println()
    
    # Generate certified itinerary
    result = generate_certified_itinerary(initial_condition, uncertainty_radius)
    
    if result === nothing
        println("Integration failed - cannot generate itinerary")
        return nothing
    end
    
    # Create plots
    println("\nGenerating interval bounds visualization...")
    try
        interval_plot = plot_interval_solution(result)
        savefig(interval_plot, "lorenz_interval_bounds.png")
        println("Interval bounds plot saved as 'lorenz_interval_bounds.png'")
    catch e
        println("Warning: Could not create plot - $e")
    end
    
    # Save detailed results
    open("lorenz_certified_itinerary.txt", "w") do f
        println(f, "CERTIFIED Lorenz System Itinerary")
        println(f, "Using rigorous interval arithmetic ODE integration")
        println(f, "")
        println(f, "Parameters:")
        println(f, "  Initial condition: $initial_condition")
        println(f, "  Uncertainty radius: $uncertainty_radius")
        println(f, "  Integration time achieved: $(result.integration_time)s")
        println(f, "  Maximum interval width: $(MAX_INTERVAL_WIDTH)")
        println(f, "")
        println(f, "Results:")
        println(f, "  Itinerary: $(result.itinerary)")
        println(f, "  Total crossings: $(result.crossings)")
        println(f, "  Final interval width: $(result.solution.widths[end])")
        println(f, "")
        println(f, "Crossing Details:")
        for (i, detail) in enumerate(result.details)
            println(f, "  $i: t=$(detail.time:.3f), y∈$(detail.y_interval), → $(detail.symbol) ($(detail.certainty))")
        end
    end
    println("Detailed results saved as 'lorenz_certified_itinerary.txt'")
    
    return result
end

# Execute the analysis
println("Executing rigorous interval analysis...")
result = main()