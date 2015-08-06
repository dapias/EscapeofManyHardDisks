push!(LOAD_PATH,"./")


include("./input_parameters.jl")
using Simulation
using Visual

#In this file the main functions of the project are called, i.e. *simulation* from main.jl and *visualizate* from visualization.jl.
#Additionally the parameters from input_parameters.jl are called.

radii = ones(N)*r
masses = ones(N)*r
h = 0.1*(L - 2*r)
#h = 2.
simulacion = simulation(tinicial, tmax, N, L, etotal, masses, radii,r, h);
@time visualize(simulacion, N, L, r, h, tmax)
