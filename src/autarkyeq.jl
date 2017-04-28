
module autarkyeq

# module to find the equilibrium quantities and prices with CES, land endowment
# and temperatures as arguments, in case of autarky

using Gadfly
using Roots
using JuMP
using NLopt
using AmplNLWriter, CoinOptServices

# PARAMETERS

# Technology parameters (do not vary from one country to another)
a = 10.0 # coef of conversion of crop into meat
b = 1e-4 # coef of conversion of land into meat
gamma1 = 3800/15 # coef of conversion of land into crop
                # http://data.worldbank.org/indicator/AG.YLD.CREL.KG
gamma2 = 600 # slope coef of how crop productivity decreases with distance to
             # optimal temperature
# Optimal temeprature
Topt = 15.0

function autarky_eq(Lmax::Float64, t::Float64, sigma::Float64)

  # distance to optimal temperature for crops
  function d(t)
    abs(Topt - t)
  end

  dist = d(t)

  # crop productivity
  function theta(d::Float64)
    if gamma1 * Topt - gamma2 *d > 0.0
      return gamma1 * Topt - gamma2 *d
    else
      return 0.0
    end
  end

  ##############################################################################
  ########################### Numerical solutions ##############################
  ##############################################################################

  ######################################
  ############### DEMAND ###############
  ######################################

  # Demand maximizes utility under the budget constraint
  # Recall that 0 < sigma < 1
    function maxuty(p)
      if gamma1 * Topt - gamma2 * dist > 0.0
        ############## PROBLEM #############
        # ne fonctionne que pour sigma <.3, sinon too high power
        # mais a déjà fonctionné correctement avec sigma = .1 et below
        m = Model(solver=BonminNLSolver())
        @variable(m, q_c >= 0.0)
        @variable(m, q_m >= 0.0)
        @NLobjective(m, Max, (q_c^((sigma-1)/sigma) + q_m^((sigma-1)/sigma))^(sigma/(sigma-1)))
        @NLconstraint(m, q_c + p*q_m == (gamma1 * Topt - gamma2 *dist)*Lmax)
        solve(m)
        return res = (getvalue(q_m), getvalue(q_c))
      else
        return (0.0, 0.0)
      end
    end

  ######################################
  ############### SUPPLY ###############
  ######################################

  # Optimal price, from profit = 0
  popt = fzeros(p -> p - theta(dist)/a - 1/b)[1]

  # maximum production of cereals possible
  maxqc = theta(dist)*Lmax

  # L_c, obtained by meat producers equating inputs
  function num_Lc(q_c::Float64)
    if maxqc == 0.0
      return 0.0
    else
      return fzeros(L_c -> a*(Lmax - L_c) - b*(theta(dist)*L_c - q_c), 0.0, maxqc)[1]
    end
  end

  # Gives Q_m
  function num_Qm(q_c::Float64)
    min( a*(Lmax-num_Lc(q_c)) , b*(theta(dist)*num_Lc(q_c)-q_c))
  end

  # Relative demand with a CES utility function
  num_RD = maxuty(popt)[1]/maxuty(popt)[2]

  function num_RS(q_c)
    return num_Qm(q_c)/q_c
  end

  # Optimal q_c from RD = RS
  if maxqc == 0.0
    num_qc = 0.0
  else
    num_qc = fzeros(q_c -> num_RD - num_RS(q_c), 0.0, maxqc)[1]
  end

  # Optimal q_m
  num_qm = num_Qm(num_qc)


  ##############################################################################
  ########################### Analytical solutions #############################
  ##############################################################################

  # Land allocated to cereals given net consumption of cereals
  function L_c(q_c::Float64)
    (a*Lmax + b*q_c)/(b*theta(dist) + a)
  end

  # Meat production as a function of net consumption of cereals
  function Q_m(q_c::Float64)
    a * (Lmax - L_c(q_c))
  end

  # Analytical solutions from RD = RS
  function q_c(sigma)
    (theta(dist) * Lmax) / (popt^(1-sigma)+1)
  end
  function q_m(sigma)
    (theta(dist) * Lmax) * (popt^(-sigma)) / (popt^(1-sigma)+1)
  end

  totcrops = theta(dist)*L_c(q_c(sigma))
  K= totcrops - q_c(sigma)

return Dict("Optimal Relative Demand" => q_m(sigma)/q_c(sigma),
            "Optimal Relative Demand num" => num_RD,
            "Optimal price" => popt,
            "Optimal net cereal consumption" => q_c(sigma),
            "Optimal net c C° num" => num_qc,
            "Optimal meat consumption" => q_m(sigma),
            "Optimal m C° num" => num_qm,
            "Percentage of feed" => (K/totcrops)*100.0,
            "Percentage of land to crop" => (L_c(q_c(sigma))/Lmax)*100.0,
            "Percentage of land to crop num" => (num_Lc(num_qc)/Lmax)*100.0)

end #function autarky_eq

###################################################################
# Plot how the equilibrium evolves with country-specific parameters
###################################################################

# We use the analytic solutions because we have them

set_default_plot_size(20cm, 15cm)
nevals = 100

## Land endowment, fixing t=15.0 and sigma = 1.0
# Land endowment and cereal consumption
gridland = linspace(1e-4, 1e6, nevals) # setting a grid for land endowment
yland = [autarky_eq(gridland[i], 15.0, .9) for i in 1:nevals] #corresponding y using function autarky_eq

landplot1 = plot(x = gridland, y = collect(yland[i]["Optimal net cereal consumption"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Land endowment"), Guide.ylabel("Opt q_c"), Guide.title("Fixing t=15.0 and sigma = .9"))
landplot2 = plot(x = gridland, y = collect(yland[i]["Optimal meat consumption"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Land endowment"), Guide.ylabel("Opt q_m"), Guide.title("Fixing t=15.0 and sigma = .9"))
landplot3 = plot(x = gridland, y = collect(yland[i]["Optimal Relative Demand"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Land endowment"), Guide.ylabel("Opt RD"), Guide.title("Fixing t=15.0 and sigma = .9"))
landplot4 = plot(x = gridland, y = collect(yland[i]["Optimal price"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Land endowment"), Guide.ylabel("Opt price"), Guide.title("Fixing t=15.0 and sigma = .9"))
#landplot5 = plot(x = gridland, y = collect(yland[i]["Percentage of land to crop"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Land endowment"), Guide.ylabel("Per Land to crop"), Guide.title("Fixing t=15.0 and sigma = 1.0"))
#landplot6 = plot(x = gridland, y = collect(yland[i]["Percentage of feed"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Land endowment"), Guide.ylabel("Per crop to feed"), Guide.title("Fixing t=15.0 and sigma = 1.0"))

landplots = gridstack([landplot1 landplot2; landplot3 landplot4])

## Temperatures, fixing Lmax=1000 and sigma = 1.0
gridtemp = linspace(5.0, 25.0, nevals) # setting a grid for land endowment
ytemp = [autarky_eq(1000.0, gridtemp[i], .9) for i in 1:nevals] #corresponding y using function autarky_eq

tempplot1 = plot(x = gridtemp, y = collect(ytemp[i]["Optimal net cereal consumption"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Temperatures"), Guide.ylabel("Opt q_c"), Guide.title("Fixing Lmax=1000 and sigma = .9"))
tempplot2 = plot(x = gridtemp, y = collect(ytemp[i]["Optimal meat consumption"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Temperatures"), Guide.ylabel("Opt q_m"), Guide.title("Fixing Lmax=1000 and sigma = .9"))
tempplot3 = plot(x = gridtemp, y = collect(ytemp[i]["Optimal Relative Demand"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Temperatures"), Guide.ylabel("Opt RD"), Guide.title("Fixing Lmax=1000 and sigma = .9"))
tempplot4 = plot(x = gridtemp, y = collect(ytemp[i]["Optimal price"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Temperatures"), Guide.ylabel("Opt price"), Guide.title("Fixing Lmax=1000 and sigma = .9"))
#tempplot5 = plot(x = gridtemp, y = collect(yland[i]["Percentage of land to crop"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Temperatures"), Guide.ylabel("Per Land to crop"), Guide.title("Fixing Lmax=1000 and sigma = 1.0"))
#tempplot6 = plot(x = gridtemp, y = collect(yland[i]["Percentage of feed"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Temperatures"), Guide.ylabel("Per crop to feed"), Guide.title("Fixing Lmax=1000 and sigma = 1.0"))

tempplots = gridstack([tempplot1 tempplot2; tempplot3 tempplot4])

# CES, fixing Lmax=1000 and t = 15.0
gridsigma = linspace(.3, 1-1e-4, nevals) # setting a grid for land endowment
ysigma = [autarky_eq(1000.0, 15.0, gridsigma[i]) for i in 1:nevals] #corresponding y using function autarky_eq

sigmaplot1 = plot(x = gridsigma, y = collect(ysigma[i]["Optimal net cereal consumption"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("Opt q_c"), Guide.title("Fixing Lmax=1000 and t = 15.0"))
sigmaplot2 = plot(x = gridsigma, y = collect(ysigma[i]["Optimal meat consumption"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("Opt q_m"), Guide.title("Fixing Lmax=1000 and t = 15.0"))
sigmaplot3 = plot(x = gridsigma, y = collect(ysigma[i]["Optimal Relative Demand"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("Opt RD"), Guide.title("Fixing Lmax=1000 and t = 15.0"))
sigmaplot4 = plot(x = gridsigma, y = collect(ysigma[i]["Optimal price"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("Opt price"), Guide.title("Fixing Lmax=1000 and t = 15.0"))
sigmaplot5 = plot(x = gridsigma, y = collect(ysigma[i]["Percentage of land to crop"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("% land to crop"), Guide.title("Fixing Lmax=1000 and t = 15.0"))
sigmaplot6 = plot(x = gridsigma, y = collect(ysigma[i]["Percentage of feed"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("% crop to feed"), Guide.title("Fixing Lmax=1000 and t = 15.0"))

sigmaplots = gridstack([sigmaplot1 sigmaplot2; sigmaplot3 sigmaplot4; sigmaplot5 sigmaplot6])

function runall(Lmax, t, sigma)
  res = autarky_eq(Lmax, t, sigma)
  println("-------------------------------------------")
  println("Solutions when Lmax = ", Lmax, "  t = ", t, "  and sigma = ", sigma)
  println("Numerical solution for q_c: ", res["Optimal net c C° num"])
  println("Numerical solution for q_m: ", res["Optimal m C° num"])
  println("Percent of cereals used as feed: ", res["Percentage of feed"])
  println("Percentage of land for crop production: ", res["Percentage of land to crop"])
  println("Optimal relative price: ", res["Optimal price"])
  println("-------------------------------------------")
end

function plotall()
display(landplots)
display(tempplots)
display(sigmaplots)
end

end #module
