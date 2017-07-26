function Feasible_Traj(SystemParams::TypeSystemParams, x0::Array{Float64,1})

    
    g   = SystemParams.g
    dt  = SystemParams.dt
    xF  = SystemParams.xF
    l   = SystemParams.l  
    m   = SystemParams.m
    b   = SystemParams.b

    Points = 100
    
    x_feasible = zeros(4, Points+1)
    u_feasible = zeros(2, Points)


    x_feasible[:,1] = x0

    for i = 1:Points
	# Logic to compute input
	if i==1
		u_feasible[2, i] =  10.0
	elseif ((x_feasible[4,i] + dt*x_feasible[3,i]-xF[4])^2<0.0001)&&((x_feasible[3,i])^2>0.001)
		u_feasible[2, i] = -10.0
	else
		u_feasible[2, i] = -0.0
	end
	u_feasible[1, i] =  + m*g*l*( sin( x_feasible[2,i] ) + cos( x_feasible[2,i])*u_feasible[2,i]  ) + b*x_feasible[1,i]

	# Propagate input forward
	x_feasible[1,i+1] = x_feasible[1,i] + dt/(m*l^2)*( u_feasible[1,i]  - m*g*l*( sin( x_feasible[2,i] ) + cos( x_feasible[2,i])*u_feasible[2,i]  ) - b*x_feasible[1,i])
	x_feasible[2,i+1] = x_feasible[2,i] + dt*( x_feasible[1,i] )
	x_feasible[3,i+1] = x_feasible[3,i] + dt*( u_feasible[2,i] )
	x_feasible[4,i+1] = x_feasible[4,i] + dt*( x_feasible[3,i] )


	#println("x_feasible " , x_feasible[4,i])
	#println("dt*x_feasible ", dt*x_feasible[3,i])
	#println("x_next ", x_feasible[4,i+1])
    end
    return x_feasible, u_feasible

end
