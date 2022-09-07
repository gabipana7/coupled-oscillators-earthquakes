# Tutorial 1 - scalar equations

using DifferentialEquations
f(u,p,t) = 1.01*u
u0 = 1/2
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

using Plots
plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false
plot!(sol.t, t->0.5*exp(1.01t),lw=3,ls=:dash,label="True Solution!")


sol
[t+u for (u,t) in tuples(sol)]


# Example 2 
function lorenz!(du,u,p,t)
     du[1] = 10.0*(u[2]-u[1])
     du[2] = u[1]*(28.0-u[3]) - u[2]
     du[3] = u[1]*u[2] - (8/3)*u[3]
    end

u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz!,u0,tspan)
sol = solve(prob)

# 3D plot in phase space
plot(sol,vars=(1,2,3))

# Timeseries of just the second component
plot(sol,vars=(0,2))


# Parametrized Functions

function parameterized_lorenz!(du,u,p,t)
     du[1] = p[1]*(u[2]-u[1])
     du[2] = u[1]*(p[2]-u[3]) - u[2]
     du[3] = u[1]*u[2] - p[3]*u[3]
end

# Or better:
function parameterized_lorenz!(du,u,p,t)
     x,y,z = u
     σ,ρ,β = p
     du[1] = dx = σ*(y-x)
     du[2] = dy = x*(ρ-z) - y
     du[3] = dz = x*y - β*z
end

u0 = [1.0,0.0,0.0]
tspan = (0.0,1.0)
p = [10.0,28.0,8/3]  # parameters
prob = ODEProblem(parameterized_lorenz!,u0,tspan,p)


# Example 3 : Nonhomogenous ODE
# pendulum
l = 1.0                             # length [m]
m = 1.0                             # mass[m]
g = 9.81                            # gravitational acceleration [m/s²]

function pendulum!(du,u,p,t)
     du[1] = u[2]                    # θ'(t) = ω(t)
     du[2] = -3g/(2l)*sin(u[1]) + 3/(m*l^2)*p(t) # ω'(t) = -3g/(2l) sin θ(t) + 3/(ml^2)M(t)
end

θ₀ = 0.01                           # initial angular deflection [rad]
ω₀ = 0.0                            # initial angular velocity [rad/s]
u₀ = [θ₀, ω₀]                       # initial state vector
tspan = (0.0,10.0)                  # time interval  

M = t->0.1sin(t)                    # external torque [Nm]

prob = ODEProblem(pendulum!,u₀,tspan,M)
sol = solve(prob)

plot(sol,linewidth=2,xaxis="t",label=["θ [rad]" "ω [rad/s]"],layout=(2,1))


# Example 4 : Other types for systems of eqs - matrix

A  = [1. 0  0 -5
      4 -2  4 -3
     -4  0  0  1
      5 -2  2  3]
u0 = rand(4,2)
tspan = (0.0,1.0)
f(u,p,t) = A*u
prob = ODEProblem(f,u0,tspan)

sol = solve(prob)
plot(sol)