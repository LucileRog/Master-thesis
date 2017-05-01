module autarkyplots

# Module to plot results from the autarky equilibrium by making Lmax, t and sigma
# vary
# We use the analytic solutions here because we have them

include("autarkyeq.jl")

using Gadfly

set_default_plot_size(20cm, 15cm)
nevals = 100

## Land endowment, fixing t=15.0 and sigma = 1.0
# Land endowment and cereal consumption
  gridland = linspace(1e-4, 1e6, nevals) # setting a grid for land endowment
  yland = [autarkyeq.autarky_eq(gridland[i], 15.0, .9) for i in 1:nevals] #corresponding y using function autarky_eq

  landplot1 = plot(x = gridland, y = collect(yland[i]["Optimal net cereal consumption"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Land endowment"), Guide.ylabel("Opt q_c"), Guide.title("Fixing t=15.0 and sigma = .9"))
  landplot2 = plot(x = gridland, y = collect(yland[i]["Optimal meat consumption"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Land endowment"), Guide.ylabel("Opt q_m"), Guide.title("Fixing t=15.0 and sigma = .9"))
  landplot3 = plot(x = gridland, y = collect(yland[i]["Optimal Relative Demand"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Land endowment"), Guide.ylabel("Opt RD"), Guide.title("Fixing t=15.0 and sigma = .9"))
  landplot4 = plot(x = gridland, y = collect(yland[i]["Optimal price"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Land endowment"), Guide.ylabel("Opt price"), Guide.title("Fixing t=15.0 and sigma = .9"))

  landplots = gridstack([landplot1 landplot2; landplot3 landplot4])

## Temperatures, fixing Lmax=1000 and sigma = 1.0
  gridtemp = linspace(5.0, 25.0, nevals) # setting a grid for land endowment
  ytemp = [autarkyeq.autarky_eq(1000.0, gridtemp[i], .9) for i in 1:nevals] #corresponding y using function autarky_eq

  tempplot1 = plot(x = gridtemp, y = collect(ytemp[i]["Optimal net cereal consumption"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Temperatures"), Guide.ylabel("Opt q_c"), Guide.title("Fixing Lmax=1000 and sigma = .9"))
  tempplot2 = plot(x = gridtemp, y = collect(ytemp[i]["Optimal meat consumption"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Temperatures"), Guide.ylabel("Opt q_m"), Guide.title("Fixing Lmax=1000 and sigma = .9"))
  tempplot3 = plot(x = gridtemp, y = collect(ytemp[i]["Optimal Relative Demand"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Temperatures"), Guide.ylabel("Opt RD"), Guide.title("Fixing Lmax=1000 and sigma = .9"))
  tempplot4 = plot(x = gridtemp, y = collect(ytemp[i]["Optimal price"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Temperatures"), Guide.ylabel("Opt price"), Guide.title("Fixing Lmax=1000 and sigma = .9"))

  tempplots = gridstack([tempplot1 tempplot2; tempplot3 tempplot4])

# CES, fixing Lmax=1000 and t = 15.0
  gridsigma = linspace(1e-2, 1-1e-2, nevals) # setting a grid for land endowment
  ysigma = [autarkyeq.autarky_eq(1000.0, 15.0, gridsigma[i]) for i in 1:nevals] #corresponding y using function autarky_eq

  sigmaplot1 = plot(x = gridsigma, y = collect(ysigma[i]["Optimal net cereal consumption"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("Opt q_c"), Guide.title("Fixing Lmax=1000 and t = 15.0"))
  sigmaplot2 = plot(x = gridsigma, y = collect(ysigma[i]["Optimal meat consumption"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("Opt q_m"), Guide.title("Fixing Lmax=1000 and t = 15.0"))
  sigmaplot3 = plot(x = gridsigma, y = collect(ysigma[i]["Optimal Relative Demand"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("Opt RD"), Guide.title("Fixing Lmax=1000 and t = 15.0"))
  sigmaplot4 = plot(x = gridsigma, y = collect(ysigma[i]["Optimal price"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("Opt price"), Guide.title("Fixing Lmax=1000 and t = 15.0"))
  sigmaplot5 = plot(x = gridsigma, y = collect(ysigma[i]["Percentage of land to crop"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("% land to crop"), Guide.title("Fixing Lmax=1000 and t = 15.0"))
  sigmaplot6 = plot(x = gridsigma, y = collect(ysigma[i]["Percentage of feed"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("% crop to feed"), Guide.title("Fixing Lmax=1000 and t = 15.0"))

  sigmaplots = gridstack([sigmaplot1 sigmaplot2; sigmaplot3 sigmaplot4; sigmaplot5 sigmaplot6])

function plotall()
  display(landplots)
  display(tempplots)
  display(sigmaplots)
end

end #module
