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
ConvHull_x = [-1.898, 0.111, 1.2604, 1.87978, 1.889, 3.2519, 1.5645, 0.3754, 0.1430, -0.1018, -3.0349, -3.0448, -1.898]
ConvHull_y = [3.2419, 3.2771,3.1966, 3.14787, 3.145, 3.0335, 2.8989, 2.809 , 2.8   , 2.8059 ,  3.1358, 3.14   , 3.2419]
plot(ConvHull_x, ConvHull_y, "-ob", label="Minimize Torque")	

@load "alpha0.001_beta1"
for i = 1:it
	plot(SS[1, 1:time[i], i], SS[2, 1:time[i], i], "-or")	
end
ConvHull_x = [-0.3157, -0.005371, 0.919, 0.9586, 3.1837, 1.3937, -1.6494, -3.97361, -4.07926, -0.3157]
ConvHull_y = [3.46142, 3.4619   , 3.416, 3.4135, 3.1027, 3.0581, 3.1398 , 3.33827, 3.35172  , 3.46142]
plot(ConvHull_x, ConvHull_y, "-or", label="Minimum Time Problem")	
legend()
xlabel("Angular Velocity [rad/s]")
ylabel("Angle [rad]")
hold(false)

@load "alpha1_beta100"

figure()
hold(true)
subplot(211)
TimePlot = collect(0:time[it]-2)*dt
plot(TimePlot, u_LMPC[1,1:time[it]-1]', "-ob", label="Minimize Torque")
legend()
xlim([0,5])
ylabel("Torque [Nm]")
xlabel("Time [s]")
subplot(212)
plot(TimePlot, u_LMPC[2,1:time[it]-1]', "-bo", label="Minimize Torque")
legend()
xlim([0,5])
ylabel("Acceleration [m/s^2]")
xlabel("Time [s]")

@load "alpha0.001_beta1"
subplot(211)
TimePlot = collect(0:time[it]-2)*dt
plot(TimePlot, u_LMPC[1,1:time[it]-1]', "-or", label="Minimum Time Problem")
legend()
subplot(212)
plot(TimePlot, u_LMPC[2,1:time[it]-1]', "-ro", label="Minimum Time Problem")
legend()

hold(false)
