module Model

VERSION < v"0.4-" && using Docile

export Object, Wall, Disk, VerticalWall, HorizontalWall, Event
export createdisks, createwalls
export move, dtcollision, collision

abstract Object
abstract Wall <: Object

@doc """Type with attributes position(r), velocity, radius, mass and two labels (lastcollision and insidetable). Both position and velocity
are vectors. The label lastcollision corresponds to the cycle within the main loop in which the Disk suffered
its last collision (see simulation in *Simulation.jl*).  """->
type Disk <:Object
  r::Array{Float64,1}
  v::Array{Float64,1}
  radius::Float64
  mass::Float64
  lastcollision ::Int
  insidetable::Bool
end

#Disk(r,v,radius) = Disk(r,v,radius,1.0,0) #Fixed mass of 1.0 and label equal to 0 (by default)
#Disk(r,v,radius, mass) = Disk(r,v,radius, mass ,0) #Label equal to 0 (by default)
Disk(r,v,radius) = Disk(r,v,radius,1.0,0,true)
Disk(r,v,radius, mass) = Disk(r,v,radius, mass ,0, true) #Label equal to 0 (by default)


@doc """Type with attributes x and y. x corresponds to its horizontal position in a Cartesian Plane
(just a number) and y represents its initial and final height in the Plane (Array of length equal to 2)."""  ->
type VerticalWall <:Wall
  x :: Float64
  y :: Array{Float64,1}
end

@doc """Type with attributes x and y. x corresponds to its horizontal extension in a Cartesian plane
(initial and final position -Array of length equal to 2- and y corresponds to its vertical position
(a number).""" ->
type HorizontalWall <:Wall
  x :: Array{Float64,1}
  y :: Float64
end


@doc """Type with attributes  referencedisk (collider1), diskorwall(collider2) and predictedcollision(label). This label makes reference to the cycle
within the main loop in which the event was predicted (see simulation in *Simulation.jl*)."""->
type Event
  referencedisk::Disk
  diskorwall ::Object
  predictedcollision :: Int
end

@doc """Creates a random Array(Vector) of *dimension* dim with limits: liminf, limsup""" ->
function randuniform(a, b, c::Int=1)
  a + rand(c)*(b - a)
end

@doc doc"""Evaluates if two Disks are overlapped""" ->
function overlap(p1::Disk, p2::Disk)
  deltar = p1.r - p2.r
  r = norm(deltar)
  return r < (p1.radius + p2.radius)
end

@doc """Creates a Disk enclosed in the square box with size L (coordinates [0.,L]*[0.,L]); and with velocity equal to [0.,0.]"""->
function createdisk(L::Float64, mass::Float64, radius::Float64)
  deltax = L
  deltay = L
  cotainfx = radius
  cotasupx = L - radius
  cotainfy = radius
  cotasupy = L - radius
  x = randuniform(cotainfx, cotasupx)
  y = randuniform(cotainfy, cotasupy)
  v = [0.,0.]
  p = Disk([x,y],v,radius, mass)
  p
end

@doc """Sample in a microcanonical way the energy hypersurface of a system of *N = length(masses)* hard disks. The arguments needed are the
total energy and an array of masses. The procedure is taken from  Dumont, R., J. Chem. Phys. 95 (12), 1991. pages: 9172-"""->
function microcanonicalsampling(etotal::Float64, masses::Array{Float64,1})
  #Esta función que calcula la energía del sistema dado un arreglo de momentos y uno de masas (el primero con 2 veces la dimensión del segundo)
  function energy(moments::Array{Float64,1}, masses::Array{Float64,1})
    psquare = moments .^ 2./2.
    e = 0.
    for i in 1:length(psquare)
      j = ceil(i/2)  #El 2 corresponde a dividir px y py.
      e += psquare[i]/masses[j]
    end
    e
  end

  #Uso un vector con distribución normal que muestrearía la superficie de una hiperesfera.
  q = randn(2*length(masses))
  #Calculo el factor de correción por el hecho de que las masas no son iguales. Entonces la superficie es como una hiperelipse.
  N = (etotal/energy(q, masses))^0.5
  #Calculo los momentos con la distribución correcta
  moments = N*q
end


@doc """Creates N Disks enclosed in the square box of size L (coordinates [0,L]*[0,L]; and with a random
velocity generated with the microcanonical distribution."""->
function createdisks(N::Int, L::Float64, etotal::Float64, masses::Array{Float64,1}, radii::Array{Float64,1})
  p = createdisk(L, masses[1], radii[1])
  #particulas = Array(typeof(p), N)
  particulas = [p]
  for i in 2:N
    superposition = true
    while(superposition)
      p = createdisk(L, masses[i], radii[i])
      arreglo = [false]
      for particula in particulas
        test = overlap(particula, p)
        push!(arreglo,test)
      end
      k = findin(arreglo,true)
      if k == []
        superposition = false
      end
    end
    push!(particulas,p)
  end

  moments = microcanonicalsampling(etotal, masses)

  for i in 1:length(particulas)
    particulas[i].v = moments[2i-1:2i]/particulas[i].mass
  end
  particulas
end

@doc """Creates the square box in which Disks will be enclosed. The box has a hole
in the right side corner of horizontal and vertical length of r+h/2 . """->
function createwalls(L::Float64, r::Float64, h::Float64)
  arreglo = Array(Wall,4)
  Lp = L -(r + h/2.)
  arreglo[1] = VerticalWall(0.,[0.,L])
  arreglo[2] = VerticalWall(L,[0.,Lp])
  arreglo[3] = HorizontalWall([0.,L],0.)
  arreglo[4] = HorizontalWall([0.,Lp],L)
  arreglo
end

@doc """Update the position of the Disk through moving it as a free particle with velocity Disk.v
during a time interval dt"""->
move(p::Disk, dt::Real) = p.r += p.v * dt

@doc doc"""Calculates the time of collision between the Disk and a VerticalWall"""->
#Hacer test cuando la pared no esté acotada por números positivos.
function dtcollision(p::Disk, V::VerticalWall)
    dt = Inf
    if p.v[1] > 0
        if p.r[1] < V.x
            dt = (V.x - (p.r[1] + p.radius))/p.v[1]
        end
    elseif p.v[1] < 0
        if p.r[1] > V.x
            dt = ((p.r[1] - p.radius) - V.x)/-p.v[1]
        end
    end
    dt
end


@doc """Calculates the time of collision between the Disk and a HorizontallWall"""->
#Hacer test cuando la pared no esté acotada por números positivos.
function dtcollision(p::Disk, H::HorizontalWall)
    dt = Inf
    if p.v[2] > 0
        if p.r[2] < H.y
            dt = (H.y - (p.r[2] + p.radius))/p.v[2]
        end
    elseif p.v[2] < 0
        if p.r[2] > H.y
            dt = ((p.r[2] - p.radius) - H.y)/-p.v[2]
        end
    end
    dt
end

@doc """Calculates the time of collision between two Disks."""->
function dtcollision(p1::Disk,p2::Disk)
    deltar = p1.r - p2.r
    deltav = p1.v - p2.v
    rdotv = dot(deltar, deltav)
    rcuadrado = dot(deltar,deltar)
    vcuadrado = dot(deltav, deltav)
    if rdotv >= 0
        return Inf
    end
    d = (rdotv)^2 -(vcuadrado)*(rcuadrado - (p1.radius + p2.radius)^2)
    if d < 0
        return Inf
    end
    #dt = min((-rdotv+ sqrt(d))/vcuadrado, (-rdotv - sqrt(d))/vcuadrado)
    dt = (rcuadrado - (p1.radius + p2.radius)^2)/(-rdotv + sqrt(d))
    return dt
end

@doc """Update the velocity vector of a disk (Disk.v) after it collides with a VerticalWall."""->
function collision(p1::Disk, V::VerticalWall )
        p1.v = [-p1.v[1], p1.v[2]]
end

@doc doc"""Update the velocity vector of a disk (Disk.v) after it collides with a HorizontallWall."""->
function collision(p1::Disk, H::HorizontalWall )
    p1.v = [p1.v[1],-p1.v[2]]
end

@doc """Update the velocity vector of two Disks after they collides."""->
function collision(p1::Disk, p2::Disk)
    deltar = p1.r - p2.r
    deltav = p1.v - p2.v
    h = dot(deltar,deltav)
    sigma = p1.radius+p2.radius
    J = 2*p1.mass*p2.mass*h/(sigma*(p1.mass + p2.mass))
    p1.v -= J*deltar/(sigma*p1.mass)
    p2.v += J*deltar/(sigma*p2.mass)
end

end
