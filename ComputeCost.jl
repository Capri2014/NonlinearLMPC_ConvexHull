function ComputeCost(x::Array{Float64,2}, u::Array{Float64,2}, LMPCparams::TypeLMPCparams, SystemParams::TypeSystemParams)

    Points = size(x)[2] 
    Qfun = zeros(Points)
    xF  = SystemParams.xF
    dt  = SystemParams.dt
    g   = SystemParams.g
    m   = SystemParams.m

    Qt  = LMPCparams.Qt
    Q1  = LMPCparams.Q1
    Q2  = LMPCparams.Q2
    
    for i = 1:Points
        index = Points - i + 1
        if i == 1
            w     = x[1, index]
            theta = x[2, index]
            v     = x[3, index]
            s     = x[4, index]


	    Norm2 = w^2 + ( theta - xF[2] )^2 + v^2 + ( s - xF[4] )^2 
	    Qfun[index] = Qt*10*Norm2 / ((10*Norm2)^2 + 1)^0.5  
	    
        else
            w     = x[1, index]
            theta = x[2, index]
            v     = x[3, index]
            s     = x[4, index]


	    Norm2 = w^2 + ( theta - xF[2] )^2 + v^2 + ( s - xF[4] )^2 

	    Qfun[index] = Qfun[index+1] + Qt*10*Norm2/( (10*Norm2)^2  + 1)^0.5 + Q1*u[1,index]^2 + Q2*u[2,index]^2 

        end
    end

    return Qfun
    
end
