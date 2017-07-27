using PyPlot


dt = SystemParams.dt

figure()
subplot(411)
hold(1)
i = 1
TimePlot  = collect(0:1:time[i]-1)*dt
println(TimePlot)
plot(TimePlot',  SS[1, 1:time[i], i], "-ro" )
 
i = it
TimePlot  = collect(0:1:time[i]-1)*dt
plot(TimePlot',  SS[1, 1:time[i], i], "-go")
grid(1)
subplot(412)
hold(1)
i = 1
TimePlot  = collect(0:1:time[i]-1)*dt
plot(TimePlot',  SS[2, 1:time[i], i], "-ro" )
 
i = it
TimePlot  = collect(0:1:time[i]-1)*dt
plot(TimePlot',  SS[2, 1:time[i], i], "-go")
grid(1)
subplot(413)
hold(1)
i = 1
TimePlot  = collect(0:1:time[i]-1)*dt
plot(TimePlot',  SS[3, 1:time[i], i], "-ro" )
 
i = it
TimePlot  = collect(0:1:time[i]-1)*dt
plot(TimePlot',  SS[3, 1:time[i], i], "-go")
grid(1)
subplot(414)
hold(1)
i = 1
TimePlot  = collect(0:1:time[i]-1)*dt
plot(TimePlot',  SS[4, 1:time[i], i], "-ro" )
 
i = it
TimePlot  = collect(0:1:time[i]-1)*dt
plot(TimePlot',  SS[4, 1:time[i], i], "-go")
grid(1)
 
title("LMPC Steady State")
axis("equal")

