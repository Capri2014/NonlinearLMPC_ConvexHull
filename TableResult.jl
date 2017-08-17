using PyPlot
using PyCall
using JLD
using Polyhedra
close("all")

Table = zeros(6,4)

@load "alpha1_beta100"
Table[1,1] = Qfun[1, 1, 1]
Table[1,2] = Qfun[1, 1, 9]
Table[1,3] = 1
Table[1,4] = 100

@load "alpha1_beta10"
Table[2,1] = Qfun[1, 1, 1]
Table[2,2] = Qfun[1, 1, 9]
Table[2,3] = 1
Table[2,4] = 10

@load "alpha1_beta1"
Table[3,1] = Qfun[1, 1, 1]
Table[3,2] = Qfun[1, 1, 9]
Table[3,3] = 1
Table[3,4] = 1


@load "alpha0.1_beta1"
Table[4,1] = Qfun[1, 1, 1]
Table[4,2] = Qfun[1, 1, 9]
Table[4,3] = 0.1
Table[4,4] = 1


@load "alpha0.01_beta1"
Table[5,1] = Qfun[1, 1, 1]
Table[5,2] = Qfun[1, 1, 9]
Table[5,3] = 0.01
Table[5,4] = 1

@load "alpha0.001_beta1"
Table[6,1] = Qfun[1, 1, 1]
Table[6,2] = Qfun[1, 1, 9]
Table[6,3] = 0.001
Table[6,4] = 1

improvement = zeros(6)
for i = 1:6
	improvement[i] = (Table[i,1]-Table[i,2]) / Table[i,1]
end
