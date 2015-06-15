include("../src/objects.jl")
using FactCheck
#push!(LOAD_PATH, "C:\\Users\\marisol\\Documentos\\GitHub\\Hard-Disk-Gas\\src")
#println(LOAD_PATH)

using Objects

facts("Disk tests") do
  D = Disk([1.0, 4.0], [2., 1.5],1.0,1.0, 0)
    @fact D.r => [1.0,4.0]
    @fact D.v => [2.0,1.5]
    @fact D.radius => 1.0
    @fact D.mass => 1.0
    @fact D.lastcollision => 0
end

facts("VerticalWall tests") do
    VW = VerticalWall(2.,[3.,6.])
    @fact VW.x => 2.0
    @fact VW.y => [3.0,6.0]
end





