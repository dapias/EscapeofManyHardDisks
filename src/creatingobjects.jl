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

@doc doc"""Creates a Disks enclosed in the box with boundaries at Lx1, Lx2, Ly1, Ly2; and with a random velocity
between vmim and vmax"""->
function createdisk(Lx1, Lx2, Ly1, Ly2, N, mass, radius)
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

function microcanonicalsampling(etotal, masses)
  function energy(moments, masses)
    moments = moments .^ 2./2.
    e = 0
    for i in 1:length(moments)
      j = ceil(i/2)  #El 2 corresponde a dividir px y py.
      e += moments[i]/masses[j]
    end
    e
  end

  #Muestreo microcanónico (tomado de Dumont, R., J. Chem. Phys. 95 (12) 9172-)
  #Uso un vector con distribución normal.
  q = randn(2*length(masses))
  #Calculo el factor de cambio de coordenadas
  N = (etotal/energy(q, masses))^0.5
  #Calculo los momentos con la distribución correcta
  moments = N*q

end


@doc """Creates N Disks enclosed in the box with boundaries at Lx1, Lx2, Ly1, Ly2; and with a random
velocity between vmin and vmax"""->
function createdisks(N, Lx1, Lx2, Ly1, Ly2, etotal, masses, radii)
  p = createdisk(Lx1, Lx2, Ly1, Ly2, N, masses[1], radii[1])
  #particulas = Array(typeof(p), N)
  particulas = [p]
  for i in 2:N
    overlap = true
    while(overlap)
      p = createdisk(Lx1, Lx2, Ly1, Ly2, N, masses[i], radii[i])
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
