using PyPlot
using PyCall
using JLD
@pyimport matplotlib as mpl
@pyimport matplotlib.patches as patches
@pyimport matplotlib.animation as anim # First set up the figure, the axis, and the plot element we want to animate


close("all")

dt = SystemParams.dt

figure()
subplot(411)
hold(1)
i = 1
TimePlot  = collect(0:1:time[i]-1)*dt
println(TimePlot)
plot(TimePlot,  SS[1, 1:time[i], i]', "-ro" )
 
i = it
TimePlot  = collect(0:1:time[i]-1)*dt
plot(TimePlot,  SS[1, 1:time[i], i]', "-go")
grid(1)
subplot(412)
hold(1)
i = 1
TimePlot  = collect(0:1:time[i]-1)*dt
plot(TimePlot,  SS[2, 1:time[i], i]', "-ro" )
 
i = it
TimePlot  = collect(0:1:time[i]-1)*dt
plot(TimePlot,  SS[2, 1:time[i], i]', "-go")
grid(1)
subplot(413)
hold(1)
i = 1
TimePlot  = collect(0:1:time[i]-1)*dt
plot(TimePlot,  SS[3, 1:time[i], i]', "-ro" )
 
i = it
TimePlot  = collect(0:1:time[i]-1)*dt
plot(TimePlot,  SS[3, 1:time[i], i]', "-go")
grid(1)
subplot(414)
hold(1)
i = 1
TimePlot  = collect(0:1:time[i]-1)*dt
plot(TimePlot,  SS[4, 1:time[i], i]', "-ro" )
 
i = it
TimePlot  = collect(0:1:time[i]-1)*dt
plot(TimePlot,  SS[4, 1:time[i], i]', "-go")
grid(1)
 
title("LMPC Steady State")
axis("equal")


figure()
subplot(211)
TimePlot = collect(0:time[i]-2)*dt
plot(TimePlot, u_LMPC[1,1:time[i]-1]', "-ob")
subplot(212)
plot(TimePlot, u_LMPC[2,1:time[i]-1]', "-bo")

figure()
i = 1
TimePlot  = collect(0:1:time[i]-1)*dt
plot(SS[4, 1:time[i], i]',  SS[2, 1:time[i], i]', "-ro")
plot(SS[4, 1:time[i], 1]',  SS[2, 1:time[i], 1]', "-ro", label="1st Iteration")
i = 2
TimePlot  = collect(0:1:time[i]-1)*dt
plot(SS[4, 1:time[i], i]',  SS[2, 1:time[i], i]', "-ko", label="2nd Iteration")
i = 3
TimePlot  = collect(0:1:time[i]-1)*dt
plot(SS[4, 1:time[i], i]',  SS[2, 1:time[i], i]', "-yo", label="3rd Iteration")
i = 4
TimePlot  = collect(0:1:time[i]-1)*dt
plot(SS[4, 1:time[i], i]',  SS[2, 1:time[i], i]', "-mo", label="4th Iteration")
i = 6
TimePlot  = collect(0:1:time[i]-1)*dt
plot(SS[4, 1:time[i], i]',  SS[2, 1:time[i], i]', "-bo", label="6th Iteration")
i = 8
TimePlot  = collect(0:1:time[i]-1)*dt
plot(SS[4, 1:time[i], i]',  SS[2, 1:time[i], i]', "-go", label="8th Iteration")
i = 9
TimePlot  = collect(0:1:time[i]-1)*dt
plot(SS[4, 1:time[i], i]',  SS[2, 1:time[i], i]', "-co", label="9th Iteration")
legend(loc="lower right")
xlabel("x-axis [m]")
ylabel("Angle [rad]")


# Now draaw the car
width  = 0.2
height = 0.1
l      = SystemParams.l

x1 = zeros(time[it]); y1 = zeros(time[it])
x2 = zeros(time[it]); y2 = zeros(time[it])
x3 = zeros(time[it]); y3 = zeros(time[it])
x4 = zeros(time[it]); y4 = zeros(time[it])
xPos = zeros(time[it], 5); yPos = zeros(time[it], 5)

headx = zeros(time[it]); heady = zeros(time[it])
tailx = zeros(time[it]); taily = zeros(time[it])
xVec  = zeros(time[it], 2); yVec = zeros(time[it], 2)

for t = 1:time[it]
	frame = t
	x1[t] = SS[4, frame, it]-width/2; y1[t] = height /2
	x2[t] = SS[4, frame, it]+width/2; y2[t] = height /2
	x3[t] = SS[4, frame, it]+width/2; y3[t] = -height /2
	x4[t] = SS[4, frame, it]-width/2; y4[t] = -height /2
	xPos[t,:] = [x1[t] x2[t] x3[t] x4[t] x1[t]]; yPos[t,:]= [y1[t] y2[t] y3[t] y4[t] y1[t]]
	headx[t] = SS[4, frame, it]  + sin(SS[2, frame, it])*l; heady[t] = height/2 - l*cos(SS[2, frame, it])
	tailx[t] = SS[4, frame, it]                           ; taily[t] = height/2
	xVec[t,:]  = [headx[t] tailx[t]]; yVec[t,:] = [heady[t] taily[t]]
end


fig = figure()
ax = plt
frame = 1#time[it] 
ax[:plot](xPos[frame, :]', yPos[frame, :]', color = "gray")
ax[:plot](xVec[frame, :]', yVec[frame, :]', color ="gray", "-o", linewidth = "2")

frame = 5
ax[:plot](xPos[frame, :]', yPos[frame, :]', "r")
ax[:plot](xVec[frame, :]', yVec[frame, :]', "-ob", linewidth = "4")
xlim([-0.2, 2.2])
ylim([-0.1, 1.5])
xlabel("x-axis [m]")

#Use External Viewer for Animation
pygui(true)

#Construct Figure and Plot Data
fig = figure()
xlim([-0.2, 2.2])
ylim([-0.1, 1.2])
global line1 = ax[:plot]([],[],"r-")[1]
global line2 = ax[:plot]([],[],"-ob", linewidth = "4")[1]

function init()
    global line1
    global line2
    line1[:set_data]([],[])
    line2[:set_data]([],[])
    return (line1,line2,None)
end

function animate(i)
    k = i + 1
    global line1
    global line2
    line1[:set_data](xPos[k,:]', yPos[k,:]')
    line2[:set_data](xVec[k,:]', yVec[k,:]')
    return (line1,line2,None)
end

#Call the animator.
myanim = anim.FuncAnimation(fig, animate, init_func=init, frames=time[it], interval=20)
