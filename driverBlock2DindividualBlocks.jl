using Plots

for alpha in (10,20,30,50)

    M=10
    N=100
    phi=1.5
    
    Nsteps = 10000
    a = 0 
    b = 10
    h = (b-a)/Nsteps
    
    # Initial positions of blocks ( bordered by zeros )
    X = zeros(M+2,N+2)
    
    for i in 2:M+1
        for j in 2:N+1
            X[i,j] = rand()
        end 
    end 
    
    # Heatmap
    heat = zeros(M+2,N+2)
    
    n1 = Nsteps + 1 
    
    
    # Matrix that remembers solutions : NxN array with each element a n1,2 matrix with time on first column and position on second 
    solutions = Matrix{Matrix{Float64}}(undef,M+2,N+2)
    for i in 2:M+1
        for j in 2:N+1
            solutions[i,j] = zeros(n1,4)
        end 
    end 
    
    
    # Create initial conditions ( they can change with each iteration of time k )
    for i in 2:M+1
        for j in 2:N+1
            #if (X[i,j]+alpha*(4*X[i,j]-X[i,j+1]-X[i,j-1]-X[i+1,j]-X[i-1,j])) >= 1.0
            # Time
            solutions[i,j][1,1] = 0
            # Position
            solutions[i,j][1,2] = X[i,j]
            #Velocity
            solutions[i,j][1,3] = 0
            # Activity
            solutions[i,j][1,4] = 0
        end
    end


    # Global iterator for time 
    for k in 2:n1


        # Actual euler method to solve differential equation
        for i in 2:M+1
            for j in 2:N+1
                # Each element of the NxN matrix is a system of 2 equations than needs solving, represented by a n1x3 matrix
                # First element of the n1x3 matrix is time : t1 = t0 + dt
                solutions[i,j][k,1] = a + (k-1) * h
                # Second element of the n1x3 matrix is position : x1 = x0 + dt*u0
                solutions[i,j][k,2] = solutions[i,j][k-1,2] + h * solutions[i,j][k-1,3]
                # Last element of the n1x3 matrix is velocity : u1 = u0 + dt*f(t0,x0)
                solutions[i,j][k,3] = solutions[i,j][k-1,3] + h * (1/phi - X[i,j] - alpha*(4*X[i,j]-X[i,j+1]-X[i,j-1]-X[i+1,j]-X[i-1,j]))
            end 
        end 


        # At the end, update each value of X that "moved" with it's new value, since X is used in f
        for i in 2:M+1
            for j in 2:N+1
                X[i,j] = solutions[i,j][k,2]
            end 
        end 

        #=
        # Then find the position of the minimum
        for i in 2:N+1
            for j in 2:N+1
                Xmin[i-1,j-1] = X[i,j]
            end
        end
        M = argmin(Xmin)
        m = M[1]
        n = M[2]
        =#
        

        # For those that do not respect failure condition after new X is computed, change velocity to 0 
        for i in 2:M+1
            for j in 2:N+1
                if (X[i,j]+alpha*(4*X[i,j]-X[i,j+1]-X[i,j-1]-X[i+1,j]-X[i-1,j])) < 1.0
                    # Time
                    #solutions[i,j][k-1,1] = 
                    # Position
                    #solutions[i,j][k-1,2] = X[i,j]
                    # Velocity
                    solutions[i,j][k,3] = 0

                    # Stop activity
                    solutions[i,j][k,4] = 0
                else 
                    heat[i,j]+=1
                    solutions[i,j][k,4] = 1
                end 
            end
        end

    end 


    # Save position, velocity for each block
    k=0
    for i in 2:M+1
        for j in 2:N+1
            plotx = plot(solutions[i,j][:,1],solutions[i,j][:,2],title="Position",
                            label="id=$i,$j",xlabel="t",ylabel="x")
            #savenamex = string(k)*"pos"*"plot"*".png"
            #savefig(plotx,savenamex)
            

            ploty = plot(solutions[i,j][:,1],solutions[i,j][:,3],title="Velocity",
                            label="id=$i,$j",xlabel="t",ylabel="y")
            savename= "./alpha"*"$alpha"*"_10x100/"*string(k)*"plot"*".png"

            plotconcat = plot(plotx,ploty,layout = grid(2, 1, heights=[0.5,0.5]))

            savefig(plotconcat,savename)

            k+=1
        end 
    end


end
