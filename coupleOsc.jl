using DifferentialEquations
using Plots


function parametrized_strogatz(du,u,p,t)
    #Coupling constant
    K=p.K
    ##order parameter
    du[1]=(du[2]+du[3]+du[4])/3.0   
    #oscillators
    du[2] = (1.0-abs2(u[2])+p.w2*im)u[2]+K*(u[1]-u[2])
    du[3] = (1.0-abs2(u[3])+p.w3*im)u[3]+K*(u[1]-u[3])
    du[4] = (1.0-abs2(u[4])+p.w4*im)u[4]+K*(u[1]-u[4])
end

u0 = [1.0+0.0im,1.0+0.0im,1.0+0.0im] ##starting points
pushfirst!(u0,sum(u0)) ## starting point of the order parameter(sum of all oscillators)
tspan = (0.0,100.0)
# K=coupling constant 
#wi=natural frequency of oscillator i 
p = (K=0.9,w2=7,w3=10,w4=13)
prob = ODEProblem(parametrized_strogatz,u0,tspan,p)
sol = solve(prob)
p1=plot(sol,vars=(0,[1,2,3,4]),labels=["Order parameter" "oscillator 1" "oscillator 2" "oscillator 3"])


par = (K = 1., N = 10, w = rand(10)) # sent to the following function

function parametrized_strogatz2(du,u,p,t)
    #Coupling constant
    K = p.K
    N = p.N
    w = p.w

    ##centroid
    du[1]=0
    for i in 1:N
        du[1]+=du[i+1]
    end
    du[1]=du[1]/N
    #oscillators
    for i in 1:N-1
        du[i+1]=(1.0-abs2(u[i+1])+w[i+1]*im)*u[i+1]+K*(u[1]-u[i+1])
    end
end

u0 = [1.0+0.0im,1.0+0.0im,1.0+0.0im,1.0+0.0im,1.0+0.0im,1.0+0.0im,1.0+0.0im,1.0+0.0im,1.0+0.0im,1.0+0.0im] ##starting points
pushfirst!(u0,sum(u0)) ## starting point of the order parameter(sum of all oscillators)
tspan = (0.0,100.0)
# K=coupling constant 
#wi=natural frequency of oscillator i 
p = (K = 1., N = 10, w = rand(10))
prob = ODEProblem(parametrized_strogatz2,u0,tspan,p)
sol = solve(prob)


p1=plot(sol,vars=(0,[1,2,3,4,5,6,7,8,9,10]),labels=["Order parameter" "oscillator 1" "oscillator 2" "oscillator 3"])
