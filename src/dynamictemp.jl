module dynamictemp

# module to find the dynamic path of temperature, emissions and meat consumption
# when temeprature is endogenised

cd("C:\\Users\\lucil\\OneDrive\\Documents\\GitHub\\Master_thesis\\src")
include("autarkyeq.jl")
#include("tradeeq.jl")
include("subsistconst.jl")

using Gadfly

# PARAMETERS
# Environment parameters
k1 = 20.0 # kg of GHG emissions per kg of meat produced
k2 = .9 # Natural regeneration of GHG
k3 = 1e-3 # coefficient of how temperature rises with GHG emissions


function dyn(nper::Int64, L, T0, sig)

  Nct = length(L) # number of countries

    #if Nct == 1
     q_m(Temp::Float64) = autarkyeq.autarky_eq(L, Temp, sig)["Optimal meat consumption"]
     q_c(Temp::Float64) = autarkyeq.autarky_eq(L, Temp, sig)["Optimal net cereal consumption"]
     subcons(Temp::Float64) = subsistconst.subsist(L, Temp, sig)
    #else
    # q_m(Temp::Vector)  = tradeeq.trade_eq(L, Temp, sig)["Optimal meat consumption"]
    # subcons(Temp::Vector) = subsistconst.subsist(L, Temp, sig)
    #end

  # Temperature dynamics
  e = zeros(nper)                 # to be replaced
  GHG = zeros(nper)               # to be replaced
  qm = zeros(nper)                # to be replaced
  qc = zeros(nper)                # to be replaced
   #if Nct == 1
    Temp = vcat(T0, ones(nper-1)) # to be replaced
      for t in 2:nper
        e[t] = k1 * q_m(Temp[t-1])      # agri emissions at t due to agri prod at t-1
        GHG[t] = k2 * GHG[t-1] + e[t]   # tot emissions at t: natural + agri emissions
        Temp[t] = T0 + k3 * GHG[t]      # temperature at t: consequence of emissions at t
      end
      qm = collect(q_m(Temp[t]) for t in 1:nper)
      qc = collect(q_c(Temp[t]) for t in 1:nper)
      sc = collect(subcons(Temp[t]) for t in 1:nper)
  # else
  #  Temp = Matrix(Nct, nper)        # to be replaced
  #  Temp[:,1] = T0                  # at t=0, temepratures are the original T0
  #    for t in 2:nper
  #      e[t] = k1 * q_m(Temp[:,t-1])    # agri emissions at t due to agri prod at t-1
  #      GHG[t] = k2 * GHG[t-1] + e[t]   # tot emissions at t: natural + agri emissions
  #      Temp[:,t] = T0 + k3 * GHG[t]    # temperature at t: consequence of emissions at t
  #   end
  #   qm = collect(q_m(Temp[:,t]) for t in 1:nper)
  #   sc = collect(subcons(Temp[:,t]) for t in 1:nper)
  # end

   return Dict("Temperatures" => Temp,
               "Agri emissions" => e,
               "Global GHG" => GHG,
               "q_m" => qm,
               "q_c" => qc,
               "Subsistence constraint" => sc)

end #function dyn

function runall(nper, L, T0, sig)
  res = dyn(nper, L, T0, sig)
  println("-------------------------------------------")
  println("Solutions when Lmax = ", L, "  intitial temperatue = ", T0, "  and sigma = ", sig)
  println("Evolution of temperature: ", res["Temperatures"])
  println("Evolution of agriculture emissions: ", res["Agri emissions"])
  println("Evolution of global GHG emissions: ", res["Global GHG"])
  println("Evolution of meat demand: ", res["q_m"])
  println("Evolution of food security (1 -> satisfied / 0 -> not): ", res["Subsistence constraint"])
  println("-------------------------------------------")

end

function plotall(nper, L, T0, sig)
  println("World where land endowment = $L, initial temperature = $T0, and sigma $sig")
  set_default_plot_size(20cm, 15cm)
  grid = collect(1:nper)
  res = dyn(nper, L, T0, sig)
    plot1 = plot(x=collect(1:nper), y=res["Temperatures"], Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Time"), Guide.ylabel("Temperature"), Guide.title("Evolution of temperature"))
    plot2 = plot(x=collect(1:nper), y=res["Agri emissions"], Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Time"), Guide.ylabel("Agri emission"), Guide.title("Evolution of Agri emissions"))
    plot3 = plot(x=collect(1:nper), y=res["Global GHG"], Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Time"), Guide.ylabel("Tot GHG"), Guide.title("Evolution of tot GHG"))
    plot4 = plot(x=collect(1:nper), y=res["q_m"], Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Time"), Guide.ylabel("Meat prod"), Guide.title("Evolution of Meat prod"))
    plot5 = plot(x=collect(1:nper), y=res["q_c"], Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Time"), Guide.ylabel("Net cereal C°"), Guide.title("Evolution of net cereal consumption"))
    plot6 = plot(x=collect(1:nper), y=res["Subsistence constraint"], Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Time"), Guide.ylabel("Subs constraint"), Guide.title("Evolution of Subsistence constraint"))
  stack = gridstack([plot1 plot2; plot3 plot4; plot5 plot6])
  display(stack)
end

end #module
