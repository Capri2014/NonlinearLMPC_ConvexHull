type LMPC_Model
    mdl::JuMP.Model

    x0::Array{JuMP.NonlinearParameter,1}
    SS::Array{JuMP.NonlinearParameter,2}
    xWarm::Array{JuMP.NonlinearParameter,2}
    uWarm::Array{JuMP.NonlinearParameter,2}
    lambWarm::Array{JuMP.NonlinearParameter,1}
    Qfun::Array{JuMP.NonlinearParameter,1}

    x_Ol::Array{JuMP.Variable,2}
    u_Ol::Array{JuMP.Variable,2}
    lamb::Array{JuMP.Variable,1}

    state_cost::JuMP.NonlinearExpression
    termi_cost::JuMP.NonlinearExpression

    function LMPC_Model(LMPCparams::TypeLMPCparams,SystemParams::TypeSystemParams, SSdim::Int64)
        println("Starting creation")
        model = new()
	
        n          = 4
        d          = 2
	
        N  = LMPCparams.N
	Q1 = LMPCparams.Q1
	Q2 = LMPCparams.Q2
	Qt = LMPCparams.Qt


	dt = SystemParams.dt
        g  = SystemParams.g
	l  = SystemParams.l
	b  = SystemParams.b
        xF = SystemParams.xF
	m  = SystemParams.m

        # Create Model
        mdl = Model(solver = IpoptSolver(print_level=0))

        # Create variables (these are going to be optimized)
        @variable( mdl, x_Ol[1:n,1:(N+1)]) 
        @variable( mdl, u_Ol[1:d,1:N])
        @variable( mdl, lamb[1:SSdim])
	
        @NLparameter(mdl, x0[1:n] == 0)
        @NLparameter(mdl, SS[1:n,1:SSdim] == 0)
        @NLparameter(mdl, Qfun[1:SSdim] == 0)
        @NLparameter(mdl, lambWarm[1:SSdim] == 0)
        @NLparameter(mdl, xWarm[1:n, 1:N+1] == 0)
        @NLparameter(mdl, uWarm[1:d, 1:N] == 0)

	for i = 1:N+1
		setvalue(x_Ol[1,i], getvalue( xWarm[1,i]) )
		setvalue(x_Ol[2,i], getvalue( xWarm[2,i]) )
		setvalue(x_Ol[3,i], getvalue( xWarm[3,i]) )
		setvalue(x_Ol[4,i], getvalue( xWarm[4,i]) )
	end

	for i = 1:N
		setvalue(u_Ol[1,i], getvalue( uWarm[1,i]) )
		setvalue(u_Ol[2,i], getvalue( uWarm[2,i]) )
	end

	for i = 1:SSdim
		setvalue(lamb[i], getvalue( lambWarm[i]) )
	end

        # System Initial Condition
        @NLconstraint(mdl, [i=1:n], x_Ol[i,1] == x0[i])         # initial condition
        println("Initializing model...")

        # System dynamics
        for i=1:N	
	    @NLconstraint(mdl, x_Ol[1,i+1] == x_Ol[1, i] + dt/(m*l^2) * ( u_Ol[1, i] - m*g*l*( sin( x_Ol[2,i] )+cos( x_Ol[2,i] )*u_Ol[2, i] ) - b*x_Ol[1,i]  ))
           
	    @NLconstraint(mdl, x_Ol[2,i+1] == x_Ol[2, i] + dt * ( x_Ol[1, i]) )
	    @NLconstraint(mdl, x_Ol[3,i+1] == x_Ol[3, i] + dt * ( u_Ol[2, i]) )
	    @NLconstraint(mdl, x_Ol[4,i+1] == x_Ol[4, i] + dt * ( x_Ol[3, i]) ) 
        end
       
	# Constraints
        # for i=1:N+1
        #     setupperbound(x_Ol[2,i],  0.2)
        # end

        for i=1:N
            setlowerbound(u_Ol[2,i], -15)
            setupperbound(u_Ol[2,i],  15)
        end

	# Constraints on lambda
        for j=1:SSdim
            setlowerbound(lamb[j,1], 0)
        end
        @NLconstraint(mdl,sum{lamb[j,1], j = 1:SSdim} == 1)

        # Convex Sefe Set Constraint
        @NLconstraint(mdl, x_Ol[1,N+1] == sum{SS[1,k] * lamb[k], k = 1:SSdim})        
        @NLconstraint(mdl, x_Ol[2,N+1] == sum{SS[2,k] * lamb[k], k = 1:SSdim})        
        @NLconstraint(mdl, x_Ol[3,N+1] == sum{SS[3,k] * lamb[k], k = 1:SSdim})        
        @NLconstraint(mdl, x_Ol[4,N+1] == sum{SS[4,k] * lamb[k], k = 1:SSdim})        
        
	# Cost definitions
        # State cost
	@NLexpression(mdl, state_cost, sum{ Q1*(u_Ol[1, j]^2), j=1:N} + sum{ Q2*(u_Ol[2, j]^2), j=1:N} 
 + sum{Qt*10*( x_Ol[1,j]^2 + (x_Ol[2,j]-xF[2])^2 + x_Ol[3,j]^2 + (x_Ol[4,j]-xF[4])^2  )/( (10*( x_Ol[1,j]^2 + (x_Ol[2,j]-xF[2])^2 + x_Ol[3,j]^2 + (x_Ol[4,j]-xF[4])^2) )^2 + 1)^0.5 , j=1:N} ) 

        # Terminal  cost
        @NLexpression(mdl, termi_cost, sum{ Qfun[j] * lamb[j] ,j=1:SSdim})


        # Objective function
        @NLobjective(mdl, Min, state_cost + termi_cost)

        # First solve
        #sol_stat=solve(mdl)
        #println("Finished solve 1: $sol_stat")
        #sol_stat=solve(mdl)
        #println("Finished solve 2: $sol_stat")
        
        model.mdl   = mdl
        model.x0    = x0
        model.Qfun  = Qfun
	model.lamb  = lamb
	model.SS    = SS
        model.x_Ol  = x_Ol
	model.xWarm = xWarm
	model.uWarm = uWarm
        model.u_Ol  = u_Ol
	model.lambWarm   = lambWarm
        model.state_cost = state_cost
        model.termi_cost = termi_cost
        
        return model
    end
end
