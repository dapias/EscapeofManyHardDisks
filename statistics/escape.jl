push!(LOAD_PATH,"../src/")
push!(LOAD_PATH,"../myDataStructures/")

include("../src/HardDiskBilliardSimulation.jl")
include("../myDataStructures/MyCollections.jl")

importall Simulation
using Model
using MyCollections
using Docile

function escape_time_and_energy(parameters)

  t_initial  = parameters[t_initial]
  t_max  = parameters[t_max]
  radius =  parameters[radius]
  mass  = parameters[mass]
  e_total = parameters[e_total]
  L = parameters[L]
  nofrealizations = parameters[nofrealizations]
  N = parameters[N]

  radii = ones(N)*radius
  masses = ones(N)*mass
  h = 0.1*(L - 2*r)

  disks, paredes, posiciones, velocidades, masas, pq, t, tiempo = startsimulation(tinitial, tmax, N, L, e_total, masses, radii, r, h)
  label = 0

  t_escape = Array(Float64,(N, nofrealizations))
  e_escape = Array(Float64,(N, nofrealizations))

  for i in 1:N

    while(!isempty(pq))
      label += 1
      evento, event_time = dequeue!(pq)
      validcollision = validatecollision(evento)

      if validcollision
        updatelabels(evento,label)
        moveparticles(disks,event_time-t)
        t = event_time
        push!(tiempo,t)

        rcondition = ((evento.referencedisk.r[1] > (L - (r + h/2.))) && (evento.referencedisk.r[2] > (L - (r + h/2.)))) #Si al momento de chocar el disco está en esta región escapa
        vcondition = velocitycondition(evento.referencedisk, L, h)

        if (rcondition && vcondition)
          tescape = t + 2*evento.referencedisk.radius/norm(evento.referencedisk.v)
          evento.referencedisk.v = [0.,0.]
          evento.referencedisk.insidetable = false
          println(tescape)
        else
          collision(evento.referencedisk,evento.diskorwall)
        end

        #       updateanimationlists(disks, posiciones,velocidades,N)
        futurecollisions!(evento.referencedisk, evento.diskorwall, disks, paredes, t, tmax, pq,label)
      end
    end
  end

  push!(tiempo, tmax)
  posiciones, velocidades, tiempo, disks, masas

end
