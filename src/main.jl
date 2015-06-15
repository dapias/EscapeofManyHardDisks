include("./objects.jl")
include("./timesandrules.jl")
include("./creatingobjects.jl")

module HardDisksSimulation

VERSION < v"0.4-" && using Docile

importall Objects
importall Rules
importall Init
import Base.isless
export simulation, energy

#This allow to use the PriorityQueue providing a criterion to select the priority of an Event.
isless(e1::Event, e2::Event) = e1.tiempo < e2.tiempo

function collisions_with_wall!(particula::Disk,paredes::Array{Wall,1},tinicial, tmax, pq,etiqueta)
    tiempo = zeros(4)
    indice = 1
    for pared in paredes
        dt = dtcollision(particula, pared)
        tiempo[indice] = dt
        indice += 1
    end
    dt,k = findmin(tiempo)
    if tinicial + dt < tmax
        Collections.enqueue!(pq,Event(tinicial+dt, particula, paredes[k], etiqueta),tinicial+dt)
    end
end

@doc doc"""Calculates the initial feasible Events and push them into the PriorityQueue with label
equal to 0"""->
function initialcollisions!(particulas::Array, paredes::Array, tinicial::Number, tmax::Number, pq)
    #Puts etiqueta by default equal to zero
    for i in 1:length(particulas)

        collisions_with_wall!(particulas[i],paredes,tinicial, tmax, pq, 0)

        for j in i+1:length(particulas) #Numero de pares sin repetición N(N-1)/2
            dt = dtcollision(particulas[i], particulas[j])
            if tinicial + dt < tmax
                Collections.enqueue!(pq,Event(tinicial+dt, particulas[i], particulas[j],0),tinicial+dt)
            end
        end
    end
    pq
end





@doc doc"""Updates the PriorityQueue pushing into it all the feasible Events that can occur after the collision
of a Disk with a Wall"""->
function futurecollisions!(particula::Disk, wall::Wall, particulas, paredes, tinicial, tmax, pq, etiqueta )
    #The wall parameter is not used but it is passed to take advantage of the multiple dispatch for futurecollisions!.

    collisions_with_wall!(particula,paredes,tinicial, tmax, pq, etiqueta)

    tiempo = Float64[]
    for p in particulas
        if particula != p
            dt = dtcollision(particula, p)
            if tinicial + dt < tmax
                Collections.enqueue!(pq,Event(tinicial+dt, particula, p, etiqueta),tinicial+dt)
            end
        end
    end
    pq
end

@doc doc"""Updates the PriorityQueue pushing into it all the possible Events that can occur after the collision
of two Disks."""->
function futurecollisions!(particula1::Disk, particula2::Disk, particulas, paredes, tinicial, tmax, pq, etiqueta)

    collisions_with_wall!(particula1,paredes,tinicial, tmax, pq, etiqueta)
    collisions_with_wall!(particula2,paredes,tinicial, tmax, pq, etiqueta)

    tiempo = Float64[]
    for p in particulas
        if (particula1 != p) & (particula2 != p)
            dt = dtcollision(particula1, p)
            if tinicial + dt < tmax
                Collections.enqueue!(pq,Event(tinicial+dt, particula1, p, etiqueta),tinicial+dt)
            end
        end
    end

    tiempo = Float64[]
    for p in particulas
        if (particula1 != p) & (particula2 != p)
            dt = dtcollision(particula2, p)
            if tinicial + dt < tmax
                Collections.enqueue!(pq,Event(tinicial+dt, particula2, p, etiqueta),tinicial+dt)
            end
        end
    end
    pq
end

@doc doc"""Calculates the total energy (kinetic) of the system."""->
function energy(masas,velocidades)
    e = 0.
    for i in 1:length(masas)
        e += masas[i]*norm(velocidades[i])^2/2.
    end
    e
end

function startingsimulation(tinicial, tmax, N, Lx1, Lx2, Ly1, Ly2, vmin, vmax)
    particulas = createdisks(N,Lx1,Lx2,Ly1,Ly2,vmin,vmax)
    paredes = createwalls(Lx1,Lx2,Ly1,Ly2)
    posiciones = [particula.r for particula in particulas]
    velocidades = [particula.v for particula in particulas]
    masas = [particula.mass for particula in particulas]
    pq = Collections.PriorityQueue()
    Collections.enqueue!(pq,Event(0.0, Disk([0.,0.],[0.,0.],1.0),Disk([0.,0.],[0.,0.],1.0), 0),0.)
    pq = initialcollisions!(particulas,paredes,tinicial,tmax, pq)
    evento = Collections.dequeue!(pq)
    t = evento.tiempo
    tiempo = [evento.tiempo]
    return particulas, paredes, posiciones, velocidades, masas, pq, t, tiempo
end

@doc """Returns true if the event was predicted after the last collision label of the Disk(s)"""->
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



function moveparticles(particulas, delta_t)
    for particula in particulas
        move(particula,delta_t)
    end
end


function updateanimationlists(particulas, posiciones,velocidades,N)
    for i in 1:N
        push!(posiciones, particulas[i].r)
        push!(velocidades, particulas[i].v)
    end
end


@doc doc"""Contains the main loop of the project. The PriorityQueue is filled at each step with Events associated
to the collider Disk(s); and the element with the highest physical priority (lowest time) is removed
from the Queue and ignored if it is physically meaningless. The loop goes until the last Event is removed
from the Data Structure, which is delimited by the maximum time(tmax)."""->
function simulation(tinicial, tmax, N, Lx1, Lx2, Ly1, Ly2, vmin, vmax)
    particulas, paredes, posiciones, velocidades, masas, pq, t, tiempo = startingsimulation(tinicial, tmax, N, Lx1, Lx2, Ly1, Ly2, vmin, vmax)
    label = 0
    while(!isempty(pq))
        label += 1
        evento = Collections.dequeue!(pq)
        validcollision = validatecollision(evento)
        if validcollision == true
            updatelabels(evento,label)
            moveparticles(particulas,evento.tiempo-t)
            t = evento.tiempo
            push!(tiempo,t)
            collision(evento.referencedisk,evento.diskorwall)
            updateanimationlists(particulas, posiciones,velocidades,N)
            futurecollisions!(evento.referencedisk, evento.diskorwall, particulas, paredes, t, tmax, pq,label)
        end
    end
    push!(tiempo, tmax)
    posiciones, velocidades, tiempo, particulas, masas
end


end


# @doc doc"""Contains the main loop of the project. The PriorityQueue is filled at each step with Events associated
# to the collider Disk(s); and at the same time the element with the highest physical priority (lowest time) is removed
# from the Queue and ignored if it is physically meaningless. The loop goes until the last Event is removed
# from the Data Structure, which is delimited by the maximum time(tmax)."""->
# function simulation(tinicial, tmax, N, Lx1, Lx2, Ly1, Ly2, vmin, vmax)
#     particulas, paredes, posiciones, velocidades, masas, pq, t, tiempo = startingsimulation(tinicial, tmax, N, Lx1, Lx2, Ly1, Ly2, vmin, vmax)
#     label = 0

#     while(!isempty(pq))
#         label += 1
#         evento = Collections.dequeue!(pq)
#         if (evento.predictedcollision >= evento.referencedisk.lastcollision)
#             if typeof(evento.diskorwall) == Disk
#                 if (evento.predictedcollision >= evento.diskorwall.lastcollision)
#                     evento.diskorwall.lastcollision = label
#                     evento.referencedisk.lastcollision = label
#                     for particula in particulas
#                         move(particula,evento.tiempo - t)
#                     end
#                     t = evento.tiempo
#                     push!(tiempo,t)
#                     collision(evento.referencedisk,evento.diskorwall)
#                     for i in 1:N
#                         push!(posiciones, particulas[i].r)
#                         push!(velocidades, particulas[i].v)
#                     end
#                     futurecollisions!(evento.referencedisk, evento.diskorwall, particulas, paredes, t, tmax, pq,label)
#                 end
#             else
#                 evento.referencedisk.lastcollision = label
#                 for particula in particulas
#                     move(particula,evento.tiempo - t)
#                 end
#                 t = evento.tiempo
#                 push!(tiempo,t)
#                 collision(evento.referencedisk,evento.diskorwall)
#                 for i in 1:N
#                     push!(posiciones, particulas[i].r)
#                     push!(velocidades, particulas[i].v)
#                 end
#                 futurecollisions!(evento.referencedisk, particulas, paredes, t, tmax, pq, label)
#             end
#         end
#     end
#     push!(tiempo, tmax)
#     posiciones, velocidades, tiempo, particulas, masas
# end


# end
