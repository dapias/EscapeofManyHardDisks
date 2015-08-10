include("escape.jl")

using HDF5
using Docile
using YAML

function createhdf5(filepath, parameters)
  h5open("$filepath.hdf5", "w") do file
    attrs(file)["Nofrealizations"] = nofrealizations
    for (key,value) in parameters
      attrs(file)[string(key)] = value
    end
  end
end

function data(filepath, parameters)
  file = h5open("$filepath.hdf5", "r+")

  tdata,edata = escape_time_and_energy(parameters)

  file["tdata"] = tdata
  file["edata"] = edata

  close(file)

end