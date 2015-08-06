include("../src/Model.jl")
include("../src/Simulation.jl")
using FactCheck
#push!(LOAD_PATH, "C:\\Users\\marisol\\Documentos\\GitHub\\Hard-Disk-Gas\\src")
#println(LOAD_PATH)

using Model
using Simulation

facts("Disk tests") do
  D = Disk([1.0, 4.0], [2., 1.5],1.0,1.0, 0, true)
  @fact D.r => [1.0,4.0]
  @fact D.v => [2.0,1.5]
  @fact D.radius => 1.0
  @fact D.mass => 1.0
  @fact D.lastcollision => 0
  @fact D.insidetable => true
end

facts("VerticalWall tests") do
  VW = VerticalWall(2.,[3.,6.])
  @fact VW.x => 2.0
  @fact VW.y => [3.0,6.0]
end

#Checaré la condición de escape. De acuerdo al archivo calculationangles.mw el ángulo de escape está entre 0.601 y 0.970 radianes,
facts("Valid escape test") do
  L = 6.27
  radius = 1/2.
  h = 0.1*(L-2*radius)
  mass = 1.0
  # Estos valores de vx y vy definen un ángulo de 0.524 radianes. Es decir, la partícula no escapa.
  vx = 0.5
  vy = sqrt(0.75)
  D = Disk([L-radius, L-radius], [vx, vy],radius, mass, 0, true)
  @fact velocitycondition(D, L, h) => false

  # Estos valores de vx y vy definen un ángulo de 0.785 radianes. Es decir, la partícula escapa.
  vx = 1.
  vy = 1.
  K = Disk([L-radius, L-radius], [vx, vy],radius, mass, 0, true)
  @fact velocitycondition(K, L, h) => true
end




