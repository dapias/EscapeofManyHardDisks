module Visual

using Simulation
using PyPlot
using PyCall

VERSION < v"0.4-" && using Docile

pygui(true)

@pyimport matplotlib.patches as patch
@pyimport matplotlib.lines as lines
@pyimport matplotlib.animation as animation



export visualize

@doc doc"""This function generates the animation of the data generated by the function *simulacionanimada* in main.jl"""->
function visualize(simulation_results, N, L, r, h, t_max)

  posiciones, velocidades, tiempo, particulas, masas = simulation_results

  pos = [[posiciones[k] for k in j:N:length(posiciones)] for j in 1:N];
  vel = [[velocidades[k] for k in j:N:length(posiciones)] for j in 1:N];

  fig = plt.figure()
  ax = fig[:add_axes]([0.05, 0.05, 0.9, 0.9])
  c = patch.Circle(pos[1][1],particulas[1].radius) #En pos[1][1] el primer 1 se refiere a la particula, en tanto que el
  #segundo se refiere al evento.
  c[:set_color]((rand(),rand(),rand()))
  circulos = [c]
  ax[:add_patch](c)

  for k in 2:N
    c = patch.Circle(pos[k][1],particulas[k].radius)
    c[:set_color]((rand(),rand(),rand()))
    push!(circulos,c)
    ax[:add_patch](c)
  end


  energy_text = plt.text(0.02,0.9,"",transform=ax[:transAxes])
  plt.gca()[:set_aspect]("equal")

  drawwalls(ax, L, r,h)
  ax[:set_xlim](0., 2*L)
  ax[:set_ylim](0., 2*L)

  function animate(i)

    z = [i/10 > t for t in tiempo]
    k = findfirst(z,false) - 1

    if k == 0
      for j in 1:N
        circulos[j][:center] = (pos[j][1][1], pos[j][1][2])
      end
      #circulos[2][:center] = (pos2[1][1],pos2[1][2])

    else
      #if tiempo[k] < i/10 < tiempo[k+1]
      for j in 1:N
        circulos[j][:center] = (pos[j][k][1] + vel[j][k][1]*(i/10-tiempo[k]), pos[j][k][2] + vel[j][k][2]*(i/10-tiempo[k]))

        #circulos[2][:center] = (pos2[k][1] + vel2[k][1]*(i/10-tiempo[k]), pos2[k][2] + vel2[k][2]*(i/10-tiempo[k]))
      end

      energy_text[:set_text]("energy = $(energy(masas, [vel[j][k] for j in 1:N]))")

    end
    return (circulos,)
  end


  anim = animation.FuncAnimation(fig, animate, frames=int(t_max*10), interval=20, blit=false, repeat = false)
end


function drawwalls(ax,L,r,h)

  ##Notation for lines: from x1,y1, to x2, y2, Line2D([x1,x2],[y1,y2])
  line1 = lines.Line2D([0.,L],[0.,0.])
  line2 = lines.Line2D([0.,0.],[0.,L])
  line3 = lines.Line2D([L,L],[0.,L - (r+h/2.)])
  line4 = lines.Line2D([0.,L - (r+h/2.)],[L , L])
  ax[:add_line](line1)
  ax[:add_line](line2)
  ax[:add_line](line3)
  ax[:add_line](line4)

end

end

