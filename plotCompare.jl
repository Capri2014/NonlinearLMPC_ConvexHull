using PyPlot
using PyCall
using JLD
using Polyhedra
close("all")
@load "alpha1_beta100"

#vrep = SimpleVRepresentation([1 1; 2 2])

figure()
hold(true)
for i = 1:it
	plot(SS[1, 1:time[i], i], SS[2, 1:time[i], i], "-ob")	
end
plot(SS[1, 1, it], SS[2, 1, it], "-ob", label="High Input Cost")	

@load "alpha0.001_beta1"
for i = 1:it
	plot(SS[1, 1:time[i], i], SS[2, 1:time[i], i], "-or")	
end
plot(SS[1, 1, it], SS[2, 1, it], "-or", label="High Time Cost")	
legend()
hold(false)

@load "alpha1_beta100"

figure()
hold(true)
subplot(211)
TimePlot = collect(0:time[it]-2)*dt
plot(TimePlot, u_LMPC[1,1:time[it]-1]', "-ob", label="High Input Cost")
legend()
xlim([0,5])
ylabel("Torque [Nm]")
xlabel("Time [s]")
subplot(212)
plot(TimePlot, u_LMPC[2,1:time[it]-1]', "-bo", label="High Input Cost")
legend()
xlim([0,5])
ylabel("Acceleration [m/s^2]")
xlabel("Time")

@load "alpha0.001_beta1"
subplot(211)
TimePlot = collect(0:time[it]-2)*dt
plot(TimePlot, u_LMPC[1,1:time[it]-1]', "-or", label="High Time Cost")
legend()
subplot(212)
plot(TimePlot, u_LMPC[2,1:time[it]-1]', "-ro", label="High Time Cost")
legend()

hold(false)
