module Init

VERSION < v"0.4-" && using Docile

importall Objects


export createdisks, createwalls

@doc """Creates an random Array(Vector) of *dimension* dim with limits: liminf, limsup""" ->
function randuniform(a, b, c=1)
  a + rand(c)*(b - a)
end

@doc doc"""Evaluates if two Disks are overlapped""" ->
function solape(p1::Disk, p2::Disk)
  deltar = p1.r - p2.r
  r = norm(deltar)
  return r < (p1.radius + p2.radius)
end

@doc doc"""Creates a Disks enclosed in the box with boundaries at Lx1, Lx2, Ly1, Ly2; and with velocity equal to [0.,0.]"""->
function createdisk(Lx1, Lx2, Ly1, Ly2, mass, radius)
  deltax = Lx2 - Lx1
  deltay = Ly2 - Ly1
  cotainfx = Lx1 + radius
  cotasupx = Lx2 - radius
  cotainfy = Ly1 + radius
  cotasupy = Ly2 - radius
  x = randuniform(cotainfx, cotasupx)
  y = randuniform(cotainfy, cotasupy)
  v = [0.,0.]
  p = Disk([x,y],v,radius, mass)
  p
end

@doc """Sample in a microcanonical way the energy hypersurface of a system of *N = length(masses)* hard disks. The arguments needed are the
total energy and an array of masses. The procedure is taken from  Dumont, R., J. Chem. Phys. 95 (12), 1991. pages: 9172-"""->
function microcanonicalsampling(etotal, masses)
  #Creo esta función que calcula la energía del sistema dado un arreglo de momentos y uno de masas (el primero con 2 veces la dimensión del segundo)
  function energy(moments, masses)
    psquare = moments .^ 2./2.
    e = 0
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


@doc """Creates N Disks enclosed in the box with boundaries at Lx1, Lx2, Ly1, Ly2; and with a random
velocity generated with the microcanonical distribution."""->
function createdisks(N, Lx1, Lx2, Ly1, Ly2, etotal, masses, radii)
  p = createdisk(Lx1, Lx2, Ly1, Ly2, masses[1], radii[1])
  #particulas = Array(typeof(p), N)
  particulas = [p]
  for i in 2:N
    overlap = true
    while(overlap)
      p = createdisk(Lx1, Lx2, Ly1, Ly2, masses[i], radii[i])
      arreglo = [false]
      for particula in particulas
        test = solape(particula, p)
        push!(arreglo,test)
      end
      k = findin(arreglo,true)
      if k == []
        overlap = false
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

@doc doc"""Creates the box in which Disks will be enclosed. Its horizontal boundaries are at Lx1 and Lx2
(Lx1 < Lx2). Its vertical boundaries are at Ly1 and Ly2 (Ly1 < Ly2)."""->
function createwalls(Lx1,Lx2,Ly1,Ly2)
  arreglo = Array(Wall,4)
  arreglo[1] = VerticalWall(Lx1,[Ly1,Ly2])
  arreglo[2] = VerticalWall(Lx2,[Ly1,Ly2])
  arreglo[3] = HorizontalWall([Lx1,Lx2],Ly1)
  arreglo[4] = HorizontalWall([Lx1,Lx2],Ly2)
  arreglo
end

end
