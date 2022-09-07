using Plots


N=30
alpha=10
phi=1.5

Nsteps = 2000
a = 0 
b = 10
h = (b-a)/Nsteps

# Initial positions of blocms ( bordered by zeros )
X = zeros(N+2,N+2,N+2)
for i in 2:N+1
    for j in 2:N+1
        for k in 2:N+1
            X[i,j,k] = rand()
        end
    end 
end 

# Heatmap
heat = zeros(N+2,N+2,N+2)

n1 = Nsteps + 1 

solutions = Array{Matrix{Float64},3}(undef,N+2,N+2,N+2)
# Matrix that remembers solutions : NxN array with each element a n1,2 matrix with time on first column and position on second 
#solutions = Matrix{Matrix{Float64}}(undef, N+2,N+2)
for i in 2:N+1
    for j in 2:N+1
        for k in 2:N+1
            solutions[i,j,k] = zeros(n1,3)
        end
    end 
end 

# Create initial conditions ( they can change with each iteration of time m )
for i in 2:N+1
    for j in 2:N+1
        for k in 2:N+1
            #if (X[i,j]+alpha*(4*X[i,j]-X[i,j+1]-X[i,j-1]-X[i+1,j]-X[i-1,j])) >= 1.0
            # Time
            solutions[i,j,k][1,1] = 0
            # Position
            solutions[i,j,k][1,2] = X[i,j,k]
            #Velocity
            solutions[i,j,k][1,3] = 0
        end 
    end
end


# Global iterator for time 
for m in 2:n1


    # Actual euler method to solve differential equation
    for i in 2:N+1
        for j in 2:N+1
            for k in 2:N+1
                # Each element of the NxN matrix is a system of 2 equations than needs solving, represented by a n1x3 matrix
                # First element of the n1x3 matrix is time : t1 = t0 + dt
                solutions[i,j,k][m,1] = a + (m-1) * h
                # Second element of the n1x3 matrix is position : x1 = x0 + dt*u0
                solutions[i,j,k][m,2] = solutions[i,j,k][m-1,2] + h * solutions[i,j,k][m-1,3]
                # Last element of the n1x3 matrix is velocity : u1 = u0 + dt*f(t0,x0)
                solutions[i,j,k][m,3] = solutions[i,j,k][m-1,3] + h * (1/phi - X[i,j,k] -
                        alpha*(6*X[i,j,k]-X[i,j+1,k]-X[i,j-1,k]-X[i+1,j,k]-X[i-1,j,k]-X[i,j,k-1]-X[i,j,k+1]))
            end
        end 
    end 


    # At the end, update each value of X that "moved" with it's new value, since X is used in f
    for i in 2:N+1
        for j in 2:N+1
            for k in 2:N+1
                X[i,j,k] = solutions[i,j,k][m,2]
            end
        end 
    end 


    # For those that do not respect failure condition after new X is computed, change velocity to 0 
    for i in 2:N+1
        for j in 2:N+1
            for k in 2:N+1
                if (X[i,j,k]+alpha*(6*X[i,j,k]-X[i,j+1,k]
                            -X[i,j-1,k]-X[i+1,j,k]-X[i-1,j,k]-X[i,j,k-1]-X[i,j,k+1])) < 1.0
                    # Time
                    #solutions[i,j][m-1,1] = 
                    # Position
                    #solutions[i,j][m-1,2] = X[i,j]
                    # Velocity
                    solutions[i,j,k][m,3] = 0
                else 
                    heat[i,j,k]+=1
                end 
            end
        end
    end

end 


#=
m=0
for i in 2:N+1
    for j in 2:N+1
        plotx = plot(solutions[i,j][:,1],solutions[i,j][:,2],title="Position",
                        label="id=$i,$j",xlabel="t",ylabel="x")
        #savenamex = string(m)*"pos"*"plot"*".png"
        #savefig(plotx,savenamex)
        

        ploty = plot(solutions[i,j][:,1],solutions[i,j][:,3],title="Velocity",
                        label="id=$i,$j",xlabel="t",ylabel="y")
        savename= string(m)*"plot"*".png"

        plotconcat = plot(plotx,ploty,layout = grid(2, 1, heights=[0.5,0.5]))

        savefig(plotconcat,savename)

        m+=1
    end 
end
=#

#print(solutions[2,2])
#print(X[2,3])
#print(heat)

heatmap(heat,title=["alpha=$alpha, phi=$phi"],xlabel="i",ylabel="j",zlabel="k")

htmp = heatmap(heat,title=["alpha=$alpha, phi=$phi"],xlabel="i",ylabel="j",zlabel="k")
savefig(htmp,"heatmap3D$alpha.png")

histvalues = Vector{Float64}()
for i in 2:N+1
    for j in 2:N+1
        for k in 2:N+1
            append!(histvalues,heat[i,j,k])
        end
    end 
end 


hst = histogram(histvalues,bins = 50,title=["alpha=$alpha, phi=$phi"])
savefig(hst,"histogram3D$alpha.png")

#scatter(heat)


