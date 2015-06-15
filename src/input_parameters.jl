"""In this file the user can adjust the values of the parameters to be passed to visualization.jl"""

export N, Lx1, Lx2, Ly1, Ly2, vmin, vmax, tinicial, tmax

N = 10 #Number of particles.
Lx1 = 0 #Position at the x-axis of the "left" vertical wall.
Lx2 = 10 #Position at the x-axis of the "right" vertical wall.
Ly1 = 0.  #Position at the y-axis of the "down" vertical wall.
Ly2 = 10 #Position at the y-axis of the "up" vertical wall.
vmin = -5.0  #Minimum velocity for the particles.
vmax = 5.0 #Maximum velocity for the particles.
tinicial = 0 #Time when simulation starts.
tmax = 100 #Time when simulation ends

