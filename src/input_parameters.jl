"""In this file the user can adjust the values of the parameters to be passed to visualization.jl"""

export N, Lx1, Lx2, Ly1, Ly2, vmin, vmax, tinicial, tmax

N = 5 #Number of particles.
Lx1 = 0 #Position at the x-axis of the "left" vertical wall.
Lx2 = 6.27 #Position at the x-axis of the "right" vertical wall.
Ly1 = 0.  #Position at the y-axis of the "down" vertical wall.
Ly2 = 6.27 #Position at the y-axis of the "up" vertical wall.
etotal = 5
r = 1/2.
mass = 1.
tinicial = 0 #Time when simulation starts.
tmax = 1000 #Time when simulation ends


