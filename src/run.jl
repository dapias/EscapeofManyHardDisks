push!(LOAD_PATH,"./")


include("./input_parameters.jl")
using Simulation
using Visual

#In this file the main functions of the project are called, i.e. *simulation* from main.jl and *visualizate* from visualization.jl.
#Additionally the parameters from input_parameters.jl are called.

simulacion = simulation(tinicial, tmax, N, Lx1, Lx2, Ly1, Ly2, etotal, masses, radii);
@time visualize(simulacion, N, Lx1, Lx2, Ly1, Ly2)
