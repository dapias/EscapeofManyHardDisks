using Dates
import YAML

include("makehdf5.jl")

### Aquí se definen las variables que se van a utilizar para guardar los archivos.
parameters = YAML.load(open("parameters.yaml"))
t_max = parameters["t_max"]
nofrealizations = parameters["nofrealizations"]

s = string(now())
hour = s[end-7:end]
foldername = "Data/$(today())/$hour/"
filename = "Ensemble_size:$nofrealizations""t_max:""$t_max"
filepath = "$foldername""$filename"

try
  mkdir("./Data/$(today())/")
end
try
  mkdir("./Data/$(today())/$hour")
end
##############

##Aquí se ejecutan las funciones que crean los archivos.
createhdf5(filepath, parameters)
data(filepath, parameters)
println("Look for the file at the \"$foldername\" folder with the name \"$filename\"")
#println(Pkg.status())
################