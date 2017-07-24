function RoadProfile(s::Float64, SystemParams::TypeSystemParams)

    xF  = SystemParams.xF
    
    theta = sin(sin( (s - xF[2]) / xF[2] * 4 * 3.14 ))

    return theta 
    
end
