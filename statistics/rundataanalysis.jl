include("extractdata.jl")

try
  mkdir("./$(foldername)/time_plots")
  mkdir("./$(foldername)/energy_plots")
end

plotehistograms(filepath, foldername)
plottdata(filepath, foldername)

println("Look for the plots into the folders \"./$(foldername)time_plots\" and \"./$(foldername)energy_plots\"")


