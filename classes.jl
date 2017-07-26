type TypeSystemParams          # parameters for MPC solver
    g::Float64
    l::Float64
    dt::Float64
    m::Float64
    b::Float64
    xF::Array{Float64,2}
    TypeSystemParams( g = 9, l = 4,  dt= 0, m=0,b=0, xF = [0 0 0 0] ) = new(g, l, dt, m,b, xF)
end

type TypeLMPCparams          # parameters for MPC solver
    Qt::Float64
    Q1::Float64
    Q2::Float64
    R::Array{Float64,2}
    N::Int64
    TypeLMPCparams( Qt=1, Q1=1, Q2=1, R=[1 1;1 1], N = 0 ) = new(Qt, Q1, Q2, R, N)
end

type TypeLMPCSol          # parameters for MPC solver
    x::Array{Float64,2}
    u::Array{Float64,2}
    lamb::Array{Float64,1}
    cost::Float64
    TypeLMPCSol( x=[1 1;1 1], u=[1 1;1 1], lamb=[1], cost = 0) = new(x, u, lamb, cost)
end
