using PyPlot
using PyCall
using JLD
using Polyhedra
close("all")
@load "alpha1_beta100"

IterationCost = zeros(9)
for i = 1:9
	IterationCost[i] = Qfun[1,1,i]
end


