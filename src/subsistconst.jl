module subsistconst

# module to find check whether the subsistence consumption constraint is
# satisfied or not

pwd()
include("autarkyeq.jl")
using Gadfly

# PARAMETERS
mu_m = 1370.0 # coefficient of nutrition for meat (here beef in cal/kg)
          # https://www.fatsecret.com/calories-nutrition/usda/ground-beef-%2895%25-lean---5%25-fat%29?portionid=63122&portionamount=100.000
mu_c = 1980.0 # coefficient of nutrition for cereals (here wheat in cal/kg)
          # https://www.fatsecret.com/calories-nutrition/usda/wheat-(sprouted)
Q = 1.15e6*mu_c + 112*mu_m # subsistence level, calibrated so that the subsistence constraint is
           # not satisfied when t is around 12.5 for a country where
           # land endowment=1000.0 ha and sig=1.0
density = 1.0 # average density of France in person per hectare

function subsist(Lmax, t, sigma)
  # To call different results
  Nct = length(Lmax)
    #if Nct == 1
      res = autarkyeq.autarky_eq(Lmax, t, sigma)
      pop = Lmax*density #population
    # else
    #   res = tradeeq.trade_eq(Lmax, t, sigma)
    #   Nct = length(Lmax)
    #   Lsum = sum(Lmax[i] for i in 1:Nct)
    #   pop = Lsum*density #population
    # end
  # Optimal consumptions of cereals and meat
  qc = res["Optimal net cereal consumption"]
  qm = res["Optimal meat consumption"]
  # In order to plot this as a graph, I make the function return 0 and 1
  if (mu_m*qm + mu_c*qc)*pop >= Q
    return 1
  else
    return 0
  end
end

# function sigthreshold(L, t)
#   a = 10.0 # coef of conversion of crop into meat
#   b = 1e-4 # coef of conversion of land into meat
#   gamma1 = 3800/15 # coef of conversion of land into crop
#                   # http://data.worldbank.org/indicator/AG.YLD.CREL.KG
#   gamma2 = 600 # slope coef of how crop productivity decreases with distance to
#                # optimal temperature
#   # Optimal temeprature
#   Topt = 15.0
#
#   function d(t)
#     abs(Topt - t)
#   end
#
#   dist = d(t)
#
#   # crop productivity
#   function theta(d::Float64)
#     if gamma1 * Topt - gamma2 *d > 0.0
#       return gamma1 * Topt - gamma2 *d
#     else
#       return 0.0
#     end
#   end
#
#   p = theta(dist)/a + 1/b
#   return 1/log(p)*(log(mu_m-(Q*p)/(theta(dist)*L))-log((Q)/(theta(dist)*L)-mu_c))
#
# end

################################################################################
#################################### PLOTS #####################################
################################################################################

set_default_plot_size(20cm, 10cm)
nevals = 100

## Land endowment, fixing t=15.0 and sigma = 1.0
# Land endowment and cereal consumption
gridland = linspace(1e-4, 1e6, nevals) # setting a grid for land endowment
yland = [subsist(gridland[i], 15.0, 1.0) for i in 1:nevals] #corresponding y using function autarky_eq
landplot = plot(x = gridland, y = yland, Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Land endowment"), Guide.ylabel("Subs const"), Guide.title("Fixing t=15.0 and sigma = 1.0"))

## Temperatures, fixing Lmax=1000 and sigma = 1.0
gridtemp = linspace(5.0, 25.0, 1e2) # setting a grid for land endowment
ytemp = [subsist(1000.0, gridtemp[i], 1.0) for i in 1:nevals] #corresponding y using function autarky_eq
tempplot = plot(x = gridtemp, y = ytemp, Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Temperatures"), Guide.ylabel("Subs const"), Guide.title("Fixing Lmax=1000 and sigma = 1.0"))

## CES, fixing Lmax=1000 and t = 15.0
gridsigma = linspace(0.0, 25.0, 1e2) # setting a grid for land endowment
ysigma = [subsist(1000.0, 15.0, gridsigma[i]) for i in 1:nevals] #corresponding y using function autarky_eq
sigmaplot = plot(x = gridsigma, y = ysigma, Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("Subs const"), Guide.title("Fixing Lmax=1000 and t = 15.0"))

function plotall()
  stack = gridstack([landplot tempplot sigmaplot])
  display(stack)
end

end #module
