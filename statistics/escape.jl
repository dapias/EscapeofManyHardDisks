push!(LOAD_PATH,"../src/")
push!(LOAD_PATH,"../myDataStructures/")

using Simulation
using Model
using MyCollections
using Docile

function escape_time_and_energy(parameters)

  t_initial  = parameters["t_initial"]
  t_max  = parameters["t_max"]
  radius =  parameters["radius"]
  mass  = parameters["mass"]
  e_total = parameters["e_total"]
  L = parameters["L"]
  nofrealizations = parameters["nofrealizations"]
  N = parameters["N"]

  radii = ones(N)*radius
  masses = ones(N)*mass
  h = 0.1*(L - 2*radius)

  t_escape = Array(Float64,(N, nofrealizations))
  e_escape = Array(Float64,(N, nofrealizations))


#   t_escape = zeros(N, nofrealizations)
#   e_escape = zeros(N, nofrealizations)


  for i in 1:nofrealizations
    disks, paredes, posiciones, velocidades, masas, pq, t, tiempo = startsimulation(t_initial, t_max, N, L, e_total, masses, radii, radius, h)
    label = 0
    escaped_disks = 0

    while(!isempty(pq))
      label += 1
      evento, event_time = dequeue!(pq)
      validcollision = validatecollision(evento)

      if validcollision
        updatelabels(evento,label)
        moveparticles(disks,event_time-t)
        t = event_time
        push!(tiempo,t)

        rcondition = ((evento.referencedisk.r[1] > (L - (radius + h/2.))) && (evento.referencedisk.r[2] > (L - (radius + h/2.)))) #Si al momento de chocar el disco está en esta región escapa
        vcondition = velocitycondition(evento.referencedisk, L, h)

        if (rcondition && vcondition)
          escaped_disks += 1
          escapetime = t + 2*evento.referencedisk.radius/norm(evento.referencedisk.v)
          evento.referencedisk.insidetable = false
          t_escape[escaped_disks, i] = escapetime
          escapeenergy = energy(evento.referencedisk)
          e_escape[escaped_disks, i] = escapeenergy
          if escaped_disks == N
            break
          end
        else
          collision(evento.referencedisk,evento.diskorwall)
        end

        #       updateanimationlists(disks, posiciones,velocidades,N)
        futurecollisions!(evento.referencedisk, evento.diskorwall, disks, paredes, t, t_max, pq,label)
      end
    end
  end

  t_escape, e_escape

end

function energy(d::Disk)
  d.mass*norm(d.v)^2./2.
end

# import YAML
# parameters = YAML.load(open("parameters.yaml"))
# escape_time_and_energy(parameters)
