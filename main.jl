# Nonlinear LMPC with CS as terminals set and P as terminal cost
# Author: Ugo Rosolia, Date 7/20/2017

using Ipopt
using JuMP
using JLD
using PyPlot
using JLD
using PyPlot

# Include the following files ======> IMPORTANT: NEED TO CHANGE ROAD PROFILE IN WARM START!!!!
include("classes.jl")             
include("LMPC_models.jl")           # Need to modify this to change road profile
include("SolveLMPCProblem.jl")      
include("ComputeFeasibleTraj.jl")   
include("ComputeCost.jl")           
include("RoadProfile.jl")         

SystemParams = TypeSystemParams()
LMPCparams   = TypeLMPCparams()
LMPCSol      = TypeLMPCSol()

# Initialize System Parameters
SystemParams.g   = 9.81
SystemParams.xF  = [0 3.14 0 2.0]
SystemParams.dt  = 0.05
SystemParams.l   = 1.0
SystemParams.m   = 0.5
SystemParams.b   = 0.01


alpha         = 0.001               # Scaling:      Good Settings = 1
beta          = 1

LMPCparams.N  = 3
LMPCparams.Qt = 10.0                   # Time:         Good Settings = 10.0
LMPCparams.Q1 = beta * alpha * 1       # Torque:       Good settings = 1
LMPCparams.Q2 =        alpha * 0.1     # Acceleration: Good Settings = 0.1
# Initial Conditions;
x0 = [0.0,3.14,0.0,0.0]

# Compute First Feasible Iteration    
x_feasible, u_feasible = Feasible_Traj(SystemParams, x0)

figure()
hold()
plot(x_feasible[2,:], x_feasible[1,:], "--g*")

# Initialize SS and Q function for first feasible iteration
Buffer = 400

SS   = zeros(4, Buffer, 20)
Qfun = zeros(1, Buffer, 20)
time = zeros(20)
time = round(Int64, time)

IndexTime = find(x -> x <=0.0001, abs(x_feasible[4,:]-SystemParams.xF[4]'*ones(1,size(x_feasible)[2])))

time[1] = IndexTime[1] + 1

x_LMPC = zeros(4,Buffer)
u_LMPC = zeros(2,Buffer-1)

x_LMPC[:,1:time[1]]   = x_feasible[:, 1:time[1]] 
u_LMPC[:,1:time[1]-1] = u_feasible[:, 1:time[1]-1]  

it = 1
SS[:, 1:time[it], it]   = x_LMPC[:,1:time[it]]
Qfun[:, 1:time[it], it] = ComputeCost(x_LMPC[:,1:time[it]], u_LMPC[:,1:time[it]], LMPCparams, SystemParams)

# Now start with the Second iteration (The first is for the feasible trajectory)
it = 2
Difference = 1
xWarm = zeros(4, LMPCparams.N+1)
uWarm = zeros(2, LMPCparams.N)

while (abs(Difference) > (1e-1))&&(it<10)
    
    # Vectorize the SS and the Q function
    SSdim = sum(time) # sum(Time) = Total number of time steps for all iterations
    ConvSS   = zeros(4, SSdim)
    ConvQfun = zeros(SSdim)
    lambWarm = zeros(SSdim)
    
    xWarm[:,1:LMPCparams.N+1] = x_LMPC[:,1:LMPCparams.N+1]
    uWarm[:,1:LMPCparams.N]   = u_LMPC[:,1:LMPCparams.N]

    Counter  = 1
    for ii = 1:it-1
        for kk = 1:time[ii]
            ConvSS[:,Counter]  = SS[:, kk, ii]
            ConvQfun[Counter]  = Qfun[1, kk, ii]

            Counter = Counter + 1
        end
    end

    # Here start the iteration
    x_LMPC      = zeros(4, Buffer)
    x_LMPC[:,1] = x0


    u_LMPC      = zeros(2, Buffer)
    cost_LMPC   = ones(1,500)

    # Define the model at the j-th iteration (Need to define it at each iterations as SS and Q function change)
    mdl    = LMPC_Model(LMPCparams,SystemParams, SSdim)
    
    # Enter the time loop for the LMPC at the j-th iteration
    t = 1
    while ((cost_LMPC[t] > (1e-2))&&(t<Buffer-1))

        if t == 1
            #solveLMPCProblem(mdl,LMPCSol, x_LMPC[:,t], ConvSS, ConvQfun) 
            solveLMPCProblem(mdl,LMPCSol, x_LMPC[:,t], ConvSS, ConvQfun, xWarm, uWarm, lambWarm) 
        else
            #solveLMPCProblem(mdl,LMPCSol, x_LMPC[:,t], ConvSS, ConvQfun) 
            solveLMPCProblem(mdl,LMPCSol, x_LMPC[:,t], ConvSS, ConvQfun, xWarm, uWarm, lambWarm) 
        end

        x_LMPC[:,t+1]  = LMPCSol.x[:,2]
        u_LMPC[:,t]    = LMPCSol.u[:,1]

	# Compute the Wart Start
	N = LMPCparams.N 
	xWarm[:, 1:N]   = LMPCSol.x[:,2:N+1]
	uWarm[:, 1:N-1] = LMPCSol.u[:,2:N]

	Counter = 1
	for ii = 1:it-1
		for kk = 1:time[ii]
			if kk < time[ii]
				xWarm[:, N+1] = SS[:, kk+1, ii] * LMPCSol.lamb[Counter]
				if kk == 1
					lambWarm[Counter] = 0
				else
					lambWarm[Counter] = LMPCSol.lamb[Counter-1]
				end
			else
				lambWarm[Counter] = LMPCSol.lamb[Counter-1] + LMPCSol.lamb[Counter]
				xWarm[:, N+1] = SS[:, kk, ii] * LMPCSol.lamb[Counter]
			end
			Counter = Counter + 1
		end
	end
	uWarm[1,N] = 0
	uWarm[2,N] = 0
	#	uWarm[1, N] = (xWarm[1,N+1] - xWarm[1,N])/SystemParams.m*SystemParams.dt + SystemParams.rho*xWarm[1,N]^2 + SystemParams.m*SystemParams.g*sin(theta) 
#
	# Extract Cost
	cost_LMPC[t+1] = LMPCSol.cost 
	println("LMPC cost at step ",t, " of iteration ", it," is ", cost_LMPC[t+1])
        

        t=t+1    
    end

    # Now post process the data after the LMPC has converged
    time[it] = t 

    # Add data to SS and Q function
    SS[:, 1:time[it], it]   = x_LMPC[:,1:time[it]]
    Qfun[:, 1:time[it], it] = ComputeCost(x_LMPC[:,1:time[it]], u_LMPC[:,1:time[it]], LMPCparams, SystemParams)

    Difference = Qfun[1,1,it-1]-Qfun[1,1,it]

    it = it + 1

end

it = it - 1
for i = 1:it
    println(i,"-th itearion cost; ", Qfun[1,1,i])
end
