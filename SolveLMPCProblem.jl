function solveLMPCProblem(mdl::LMPC_Model,LMPCSol::TypeLMPCSol,xCurr::Array{Float64,1},ConvSS::Array{Float64,2},ConvQfun::Array{Float64,1},xWarm::Array{Float64,2},uWarm::Array{Float64,2},lambWarm::Array{Float64,1})

    # Load Parameters
    sol_status::Symbol

    # Update current initial condition, curvature and previous input
    setvalue(mdl.x0,xCurr)
    setvalue(mdl.SS,ConvSS)
    setvalue(mdl.Qfun,ConvQfun)
    setvalue(mdl.xWarm,xWarm)
    setvalue(mdl.uWarm,uWarm)
    setvalue(mdl.lambWarm,lambWarm)
    
    # Solve Problem and return solution
    sol_status  = solve(mdl.mdl)

    LMPCSol.x    = getvalue(mdl.x_Ol)
    LMPCSol.u    = getvalue(mdl.u_Ol)
    LMPCSol.lamb = getvalue(mdl.lamb)

    println("State Cost: ", getvalue(mdl.state_cost))
    println("Terminal Cost: ", getvalue(mdl.termi_cost))
    LMPCSol.cost = getvalue(mdl.state_cost) + getvalue(mdl.termi_cost)

    println("Solved, status = $sol_status")
    nothing
end
