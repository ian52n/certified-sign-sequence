module UncertifiedIntegrate
export uncertified_integrate

using OrdinaryDiffEq
using Random

# Lorenz dynamics
function lorenz_ode!(du, u, p, t)
    σ, ρ, β = 10.0, 28.0, 8/3
    du[1] = σ * (u[2] - u[1])
    du[2] = u[1] * (ρ - u[3]) - u[2]
    du[3] = u[1] * u[2] - β * u[3]
end

"""
    uncertified_integrate(u0::NTuple{3,Float64}, T::Float64; rng, reltol=1e-6, abstol=1e-6)

Integrate Lorenz system starting from `u0` until the trajectory has intersected the
z=27 plane with dz/dt>0 at least 5 times (tracking the sign of y). If `T` is reached first,
a new random initial point is drawn from `rng` and retried. Once 5 crossings are reached,
integrates an additional 10 time units and returns the final state.
"""
function uncertified_integrate(u0::NTuple{3,Float64}, T::Float64; rng::AbstractRNG, reltol=1e-6, abstol=1e-6)
    while true
        u = collect(u0)
        crossings = Ref(0)
        last_y_sign = 0.0

        # Condition for z-plane crossing with dz/dt>0
        function condition(u, t, integrator)
            x, y, z = u
            dz = x*y - (8/3)*z
            return z - 27.0
        end

        function affect!(integrator)
            x, y, z = integrator.u
            dz = x*y - (8/3)*z
            if dz > 0
                current_y_sign = sign(y)
                if last_y_sign != 0.0 && current_y_sign != last_y_sign
                    crossings[] += 1
                end
                last_y_sign = current_y_sign
            end
            if crossings[] >= 5
                terminate!(integrator)
            end
        end

        cb = ContinuousCallback(condition, affect!; rootfind=true)

        prob = ODEProblem(lorenz_ode!, u, (0.0, T))
        sol = solve(prob, Tsit5(); callback=cb, reltol=reltol, abstol=abstol, save_everystep=false)

        if crossings[] < 5
            # T reached before 5 crossings → generate new random point
            u0 = (rand(rng)*40-20, rand(rng)*40-20, rand(rng)*50 + 0.1)
            println("T reached before 5 crossings. Retrying with new point: $(u0)")
            continue
        end

        # Integrate an additional 10 time units from the last state
        u_last = sol.u[end]
        prob_extra = ODEProblem(lorenz_ode!, u_last, (0.0, 10.0))
        sol_extra = solve(prob_extra, Tsit5(); reltol=reltol, abstol=abstol, save_everystep=false)
        u_final = sol_extra.u[end]
        return (Float64(u_final[1]), Float64(u_final[2]), Float64(u_final[3]))
    end
end

end # module
