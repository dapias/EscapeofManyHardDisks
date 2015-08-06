include("./Model.jl")
include("../myDataStructures/MyCollections.jl")

module Simulation

VERSION < v"0.4-" && using Docile

importall Model
importall MyCollections
export simulation, energy


@doc doc"""Contains the main loop of the project. The PriorityQueue is filled at each step with Events associated
to the collider Disk(s); and the element with the highest physical priority (lowest time) is removed
from the Queue and ignored if it is physically meaningless. The loop goes until the last Event is removed
from the Data Structure, which is delimited by the maximum time(tmax)."""->
function simulation(tinitial::Float64, tmax::Float64, N::Int, L::Float64, etotal::Float64, masses::Array{Float64,1}, radii::Array{Float64,1}, r::Float64, h::Float64)
  disks, paredes, posiciones, velocidades, masas, pq, t, tiempo = startsimulation(tinitial, tmax, N, L, etotal, masses, radii, r, h)
  label = 0
  while(!isempty(pq))
    label += 1
    evento, event_time = dequeue!(pq)
    validcollision = validatecollision(evento)
    if validcollision
      updatelabels(evento,label)
      moveparticles(disks,event_time-t)
      t = event_time
      #      t = evento.tiempo
      push!(tiempo,t)

      rcondition = ((evento.referencedisk.r[1] > (L - (r + h/2.))) && (evento.referencedisk.r[2] > (L - (r + h/2.)))) #Si al momento de chocar el disco está en esta región escapa
      vcondition = velocitycondition(evento.referencedisk, L, h)
      if (rcondition && vcondition)
        evento.referencedisk.r += evento.referencedisk.v*0.2
        tescape = t + 2*evento.referencedisk.radius/norm(evento.referencedisk.v)
        evento.referencedisk.v = [0.,0.]
        evento.referencedisk.insidetable = false
        println(tescape)
      else
        collision(evento.referencedisk,evento.diskorwall)
      end
      updateanimationlists(disks, posiciones,velocidades,N)
      futurecollisions!(evento.referencedisk, evento.diskorwall, disks, paredes, t, tmax, pq,label)
    end
  end
  push!(tiempo, tmax)
  posiciones, velocidades, tiempo, disks, masas
end

function velocitycondition(d::Disk, L::Float64, h::Float64)
  vx = d.v[1]
  vy = d.v[2]
  m = vy/vx
  theta = atan(m)
  r = d.radius
  Lp = L - (r + h/2.)
  condition1 = (vx > 0. && vy > 0.)
  a = Lp + sin(theta)*r - (L - d.r[2] - cos(theta)*r)/m
  b = L - sin(theta)*r  + (d.r[2] - cos(theta)*r - Lp)/m
  condition2 = a < d.r[1] < b
  condition1 && condition2
end


function startsimulation(tinitial, tmax, N, L, etotal, masses, radii,r, h)
  disks = createdisks(N,L,etotal,masses,radii)
  paredes = createwalls(L,r, h)
  posiciones = [disk.r for disk in disks]
  velocidades = [disk.v for disk in disks]
  #masas = [disk.mass for disk in disks]
  pq = PriorityQueue()
  enqueue!(pq,Event(Disk([0.,0.],[0.,0.],1.0),Disk([0.,0.],[0.,0.],1.0), 0),0.)
  pq = initialcollisions!(disks,paredes,tinitial,tmax, pq)
  evento, t = dequeue!(pq)
  #t = evento.tiempo
  tiempo = [t]
  return disks, paredes, posiciones, velocidades, masses, pq, t, tiempo
end

@doc """Calculates the next collision with a wall of a Disk and put it in the Priority Queue with a given label"""->
function collisions_with_wall!(disk::Disk,paredes::Array{Wall,1},tinitial::Float64, tmax::Float64, pq, label::Int)
  tiempo = zeros(4)
  indice = 1
  for pared in paredes
    dt = dtcollision(disk, pared)
    tiempo[indice] = dt
    indice += 1
  end
  dt,k = findmin(tiempo)
  if tinitial + dt < tmax
    enqueue!(pq,Event(disk, paredes[k], label),tinitial+dt)
  end
end

@doc """Calculates the initial feasible Events and push them into the PriorityQueue with label
equal to 0"""->
function initialcollisions!(disks::Array, paredes::Array, tinitial::Number, tmax::Number, pq)

  #Puts etiqueta by default equal to zero

  for i in 1:length(disks)
    collisions_with_wall!(disks[i],paredes,tinitial, tmax, pq, 0)

    for j in i+1:length(disks) #Numero de pares sin repetición N(N-1)/2
      dt = dtcollision(disks[i], disks[j])
      if tinitial + dt < tmax
        enqueue!(pq,Event(disks[i], disks[j],0),tinitial+dt)
      end
    end

  end
  pq
end


@doc """Updates the PriorityQueue pushing into it all the feasible Events that can occur after the collision
of a Disk with a Wall"""->
function futurecollisions!(disk::Disk, wall::Wall, disks, paredes, tinitial, tmax, pq, etiqueta )

  #The wall parameter is not used but it is passed to take advantage of the multiple dispatch for futurecollisions!.
  if disk.insidetable
    collisions_with_wall!(disk,paredes,tinitial, tmax, pq, etiqueta)
    #   tiempo = Float64[]
    for p in disks
      if (disk != p && p.insidetable)
        dt = dtcollision(disk, p)
        if tinitial + dt < tmax
          enqueue!(pq,Event(disk, p, etiqueta),tinitial+dt)
        end
      end
    end
  end
  pq
end

@doc doc"""Updates the PriorityQueue pushing into it all the possible Events that can occur after the collision
of two Disks."""->
function futurecollisions!(disk1::Disk, disk2::Disk, disks, paredes, tinitial, tmax, pq, etiqueta)

  if disk1.insidetable
    collisions_with_wall!(disk1,paredes,tinitial, tmax, pq, etiqueta)
    for p in disks
      if (disk1 != p && p.insidetable)
        dt = dtcollision(disk1, p)
        if tinitial + dt < tmax
          enqueue!(pq,Event(disk1, p, etiqueta),tinitial+dt)
        end
      end
    end
  end

  if disk2.insidetable
    collisions_with_wall!(disk2,paredes,tinitial, tmax, pq, etiqueta)
    for p in disks
      if (((disk1 != p) && (disk2 != p)) && p.insidetable)
        dt = dtcollision(disk2, p)
        if tinitial + dt < tmax
          enqueue!(pq,Event(disk2, p, etiqueta),tinitial+dt)
        end
      end
    end
  end
  pq
end


@doc """Returns true if the event was predicted after the lastcollision label of the Disk(s)"""->
function validatecollision(event::Event)
  validcollision = false
  function validate(d::Disk)
    if (event.predictedcollision >= event.diskorwall.lastcollision)
      validcollision = true
    end
  end
  function validate(w::Wall)
    validcollision  = true
  end

  if event.predictedcollision >= event.referencedisk.lastcollision
    validate(event.diskorwall)
  end
  validcollision
end

@doc """Update the lastcollision label of the Disk(s) with the label of the loop"""->
function updatelabels(evento::Event,label)
  function update(d1::Disk,d2::Disk)
    evento.diskorwall.lastcollision = label
    evento.referencedisk.lastcollision = label
  end

  function update(d1::Disk,w::Wall)
    evento.referencedisk.lastcollision = label
  end

  update(evento.referencedisk,evento.diskorwall)
end



function moveparticles(disks, delta_t)
  for disk in disks
    move(disk,delta_t)
  end
end


function updateanimationlists(disks, posiciones,velocidades,N)
  for i in 1:N
    push!(posiciones, disks[i].r)
    push!(velocidades, disks[i].v)
  end
end


@doc """Calculates the total energy (kinetic) of the system."""->
function energy(masas,velocidades)
  e = 0.
  for i in 1:length(masas)
    e += masas[i]*norm(velocidades[i])^2/2.
  end
  e
end


end

