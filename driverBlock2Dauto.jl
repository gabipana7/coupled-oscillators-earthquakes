using Plots

for phi in (3.0,0.1)
    for alpha in (1,10,30,50)

        
        N=50
        Nsteps = 10000
        a = 0 
        b = 10
        h = (b-a)/Nsteps

        # Initial positions of blocks ( bordered by zeros )
        X = zeros(N+2,N+2)

        for i in 2:N+1
            for j in 2:N+1
                X[i,j] = rand()
            end 
        end 

        # Heatmap
        heat = zeros(N+2,N+2)


        n1 = Nsteps + 1 

        #=
        # Activity
        activity = Matrix{Vector{Float64}}(undef,N+2,N+2)

        for i in 2:N+1
            for j in 2:N+1
                activity[i,j] = Array{Float64,1}(undef,n1)
            end 
        end
        =#


        # Matrix that remembers solutions : NxN array with each element a n1,2 matrix with time on first column and position on second 
        solutions = Matrix{Matrix{Float64}}(undef,N+2,N+2)

        for i in 2:N+1
            for j in 2:N+1
                solutions[i,j] = zeros(n1,4)
            end 
        end 


        # Create initial conditions ( they can change with each iteration of time k )
        for i in 2:N+1
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


        Xintermediary=zeros(N,N)
        # Find the first minimum
        for i in 2:N+1
            for j in 2:N+1
                Xintermediary[i-1,j-1] = X[i,j]
            end
        end

        minima = zeros(n1,2)
        minima[1,1] = argmin(Xintermediary)[1]
        minima[1,2] = argmin(Xintermediary)[2]

        maxima = zeros(n1,2)
        maxima[1,1] = argmax(Xintermediary)[1]
        maxima[1,2] = argmax(Xintermediary)[2]

        differenceMIN = Vector{Float64}()
        differenceMAX = Vector{Float64}()

        #-------------------------------------- EULER ODE ---------------------------------
        # Global iterator for time 
        for k in 2:n1

            # Actual euler method to solve differential equation
            for i in 2:N+1
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
            # ------------------------------------------------------------------------

            # At the end, update each value of X that "moved" with it's new value, since X is used in f
            for i in 2:N+1
                for j in 2:N+1
                    X[i,j] = solutions[i,j][k,2]
                end 
            end 
            # -----------------------------------------------------------------------

            # For those that do not respect failure condition after new X is computed, change velocity to 0 
            for i in 2:N+1
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
            # ----------------------------------------------------------------------------


            # FIND MINIMA AND MAXIMA POSITION --------------------------------------------
            # Intermediary matrix save for argmin and argmax to work
            Xintermediary=zeros(N,N)
            for i in 2:N+1
                for j in 2:N+1
                    Xintermediary[i-1,j-1] = X[i,j]
                end
            end

            # MINIMA
            minima[k,1] = argmin(Xintermediary)[1]
            minima[k,2] = argmin(Xintermediary)[2]

            dMIN = sqrt((minima[k,1]-minima[k-1,1])^2+(minima[k,2]-minima[k-1,2])^2) 
            append!(differenceMIN,dMIN)

            # MAXIMA
            maxima[k,1] = argmax(Xintermediary)[1]
            maxima[k,2] = argmax(Xintermediary)[2]

            dMAX = sqrt((maxima[k,1]-maxima[k-1,1])^2+(maxima[k,2]-maxima[k-1,2])^2) 
            append!(differenceMAX,dMAX)
            #-----------------------------------------------------------------------------

        end 


        # ------------------------------------STATISTICS -------------------------------------
        savepath = "./phi"*"$phi"*"/"

        # ACTIVE STATIC SITES AT EACH STEP
        activity = zeros(n1,2)
        for k in 2:n1
            for i in 2:N+1
                for j in 2:N+1
                    # ACTIVE
                    activity[k,1] += solutions[i,j][k,4]
                end 
            end
            # INACTIVE
            activity[k,2] = N*N - activity[k,1]
        end
        activity
        counter = convert(Int64,(n1+1)/2)
        activityDifference = Vector{Float64}()
        for k in 2:counter-1
            if activity[2*k,1]-activity[2*k-1,1] != 0
                append!(activityDifference,activity[2*k,1]-activity[2*k-1,1])

            end
        end
        histogram(activityDifference)
        actdif = histogram(activityDifference,xlabel="activity difference",ylabel="count",
                        title=["alpha=$alpha, phi=$phi"])
        savefig(actdif,savepath*"alpha$alpha"*"_activityDifference.png")


        #=
        activityDifference = [activity[k,1]-activity[k,2] for k in 2:n1]
        histogram(activityDifference)

        activityDifference2 = Vector{Float64}()
        for item in activityDifference
            if item != -2500 
                append!(activityDifference2,item)
            end 
        end      
        histogram(activityDifference2,title=["alpha=$alpha, phi=$phi"])   
        actDif = histogram(activityDifference2,title=["alpha=$alpha, phi=$phi"],xlabel="Activity Change",ylabel="Count")
        savefig(actDif,savepath*"alpha$alpha"*")actDif.png")

        io = open(savepath*"alpha$alpha"*"_activityDifference.txt", "w") do io
            for x in activityDifference2
            println(io, x)
            end
        end
        =#


        #=
        # MINIMA MAXIMA HISTOGRAMS
        differenceMIN
        diffMIN = Vector{Float64}()
        for item in differenceMIN
            if item != 0 
                append!(diffMIN,item)
            end 
        end 
        histogram(diffMIN)
        hstmin = histogram(diffMIN,title=["alpha=$alpha, phi=$phi"],xlabel="distance between minima")
        savename = "alpha$alpha"*"_minimaDistance.png"
        savefig(hstmin,savepath*savename)
        =#

        differenceMAX
        diffMAX = Vector{Float64}()
        for item in differenceMAX
            if item != 0 
                append!(diffMAX,item)
            end 
        end 
        histogram(diffMAX)
        hstmax = histogram(diffMAX,title=["alpha=$alpha, phi=$phi"],xlabel="distance between maxima")
        savename = "alpha$alpha"*"_maximaDistance.png"
        savefig(hstmax,savepath*savename)

        # ----------------------------------------------------------------------------------

        #=
        # POSITION - VELOCITY PLOTS PER BLOCK
        k=0
        for i in 2:N+1
            for j in 2:N+1
                plotx = plot(solutions[i,j][:,1],solutions[i,j][:,2],title="Position",
                                label="id=$i,$j",xlabel="t",ylabel="x")
                #savenamex = string(k)*"pos"*"plot"*".png"
                #savefig(plotx,savenamex)
                

                ploty = plot(solutions[i,j][:,1],solutions[i,j][:,3],title="Velocity",
                                label="id=$i,$j",xlabel="t",ylabel="y")
                savename="/home/gabi/A1work/earthquakeNetworks/earthquakeModels/driverBlockJulia/alpha"*"$alpha"*"_100x100/"*string(k)*"plot"*".png"

                plotconcat = plot(plotx,ploty,layout = grid(2, 1, heights=[0.5,0.5]))

                savefig(plotconcat,savename)

                k+=1
            end 
        end
        # ----------------------------------------------------------------------------------
        =#

        # HEATMAP OF MOVING BLOCKS 
        heatmap(heat,title=["alpha=$alpha, phi=$phi"],xlabel="i",ylabel="j")

        htmp = heatmap(heat,title=["alpha=$alpha, phi=$phi"],xlabel="i",ylabel="j")
        savefig(htmp,savepath*"alpha$alpha"*"_heatmap.png")
        # -----------------------------------------------------------------------------------'=


        # TOTAL MOVING TIMES PER ALL BLOCK - CENTRAL THEOREM -> GAUSSIAN
        histvalues = Vector{Float64}()
        for i in 2:N+1
            for j in 2:N+1
                if heat[i,j] != 0
                    append!(histvalues,heat[i,j])
                end
            end 
        end 

        hst = histogram(histvalues,title=["alpha=$alpha, phi=$phi"],xlabel="total activity all blocks",ylabel="count")
        savefig(hst,savepath*"alpha$alpha"*"_activityALL.png")
        # ------------------------------------------------------------------



        # HISTOGRAM OF TOTAL ACTIVITY EXPORTED FOR PYTHON
        #=
        histvalues = Vector{Float64}()
        for i in 2:N+1
            for j in 2:N+1
                if heat[i,j] != 0
                    append!(histvalues,heat[i,j])
                end
            end 
        end 

        io = open("data2D$alpha.txt", "w") do io
            for x in histvalues
            println(io, x)
            end
        end
        =#




        
        # CHOOSE THE HIGHEST ACTIVITY BLOCK 
        m = argmax(heat)[1]
        n = argmax(heat)[2]

        
        # VELOCITY AND POSITION PLOTS FOR HIGHEST ACTIVITY BLOCK
        plotx = plot(solutions[m,n][:,1],solutions[m,n][:,2],title="Position",
                                label="id=$m,$n",xlabel="t",ylabel="x")


        ploty = plot(solutions[m,n][:,1],solutions[m,n][:,3],title="Velocity",
                        label="id=$m,$n",xlabel="t",ylabel="u")


        plotconcat = plot(plotx,ploty,layout = grid(2, 1, heights=[0.5,0.5]),
                            title=["alpha=$alpha, phi=$phi"])

        savefig(plotconcat,savepath*"alpha$alpha"*"_posVel"*".png")
        
        
        # ACTIVITY DISTRIBUTION PER BLOCK 
        #print(solutions[m,n])
        plot(solutions[m,n][:,1],solutions[m,n][:,4])

        changeInActivity = Vector{Float64}()

        for k in 2:n1
            if solutions[m,n][k,4] != solutions[m,n][k-1,4]
                append!(changeInActivity,solutions[m,n][k,1])
            end
        end

        #print(changeInActivity)
        activityTimes = Vector{Float64}()

        counter=0
        for item in changeInActivity
            counter+=1
        end
        counter = counter/2
        counter = convert(Int64,counter)

        for k in 1:counter
            append!(activityTimes,changeInActivity[2*k]-changeInActivity[2*k-1])
        end

        histogram(activityTimes,title=["alpha=$alpha, phi=$phi"],label="id=$m,$n")
        hs = histogram(activityTimes,title=["alpha=$alpha, phi=$phi"],label="id=$m,$n",xlabel="activity times for a block",ylabel="count")
        savefig(hs,savepath*"alpha$alpha"*"_activityONEBLOCK.png")

        size(activityTimes)

        io = open(savepath*"alpha$alpha"*"_activityONEBLOCK.txt", "w") do io
            for x in activityTimes
            println(io, x)
            end
        end


    end
end
