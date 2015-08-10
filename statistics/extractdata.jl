using PyPlot
using HDF5
using StatsBase

function getN(filepath)
  file = h5open("$filepath.hdf5","r")
  N = read(attrs(file)["N"])
  close(file)
  N
end


function gettdata(filepath)
  file = h5open("$filepath.hdf5","r")
  tdata = read(file["tdata/"])
  close(file)
  tdata
end

function getedata(filepath)
  file = h5open("$filepath.hdf5","r")
  edata = read(file["edata/"])
  close(file)
  edata
end


#Debo tener cuidado con los puntos en los que el array no ha sido reemplazado y contiene ceros. La estadística puede salir mal por esto.
#Debo pensar cóno descartar esos ceros.
function plotehistograms(filepath, foldername)
  edata = getedata(filepath)
  N = getN(filepath)
  for i in 1:N
    e = edata[i,:][:]
    e = e[e .> 1.e-100]
    try
      plt.clf()
      plt.hist(e, 20)
      plt.gca()[:set_xlabel](L"$Energy$")
      plt.gca()[:set_ylabel](L"$Frequency$")
      savefig("./$(foldername)/energy_plots/ehistogramparticle$i.png")
    end
  end
end

#Falta poner etiquetas en los ejes.
function plottdata(filepath, foldername)
  tdata = gettdata(filepath)
  N = getN(filepath)
  plt.clf()
  plt.gca()[:set_yscale]("log")
  plt.gca()[:set_xscale]("log")
  plt.gca()[:set_xlabel](L"$t$")
  plt.gca()[:set_ylabel](L"$P_n(t)$")
  for i in 1:N
    t = tdata[i,:][:]
    t = t[t .> 1.e-100]   #Quita elementos que no quiero en el arreglo
    try                   #Trata de calcular las funciones y graficar cuando length(t) es diferente de cero
      empirica = ecdf(t)
      max = findmax(t)[1]
      intervalo = [0.:.01:max]
      Ps = 1 - empirica(intervalo)
      plt.plot(intervalo, Ps, ".--", label = "n = $i")
    end
  end
  handles, labels = plt.gca()[:get_legend_handles_labels]()
  plt.gca()[:legend](handles, labels, loc =3)

  savefig("./$(foldername)/time_plots/tcumulative.png")

end

