
module autarkyeq

# module to find the equilibrium quantities and prices with CES, land endowment
# and temperatures as arguments, in case of autarky

using Gadfly
using Roots
using JuMP
using NLopt
using Ipopt
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

  ######################################
  ############### SUPPLY ###############
  ######################################

    # distance to optimal temperature for crops
    function dist(t)
      abs(Topt - t)
    end

    # crop productivity
    function theta(d::Float64)
      if gamma1 * Topt - gamma2 *d > 0.0
        return gamma1 * Topt - gamma2 *d
      else
        return 0.0
      end
    end

    # Optimal price, from profit = 0
    function popt(d::Float64)
      theta(d)/a + 1/b
    end

    ######## Analytical solutions #########

    # Land allocated to cereals given net consumption of cereals
    function L_c(q_c::Float64, Lmax::Float64, d::Float64)
      (a*Lmax + b*q_c)/(b*theta(d) + a)
    end

    # Meat production as a function of net consumption of cereals
    function Q_m(q_c::Float64, Lmax::Float64, d::Float64)
      min( a*(Lmax-L_c(q_c, Lmax, d)) , b*(theta(d)*L_c(q_c, Lmax, d)-q_c))
    end

  ######################################
  ############### DEMAND ###############
  ######################################

  ######## Analytical solution #########

  function q_c(Lmax::Float64, d::Float64, sigma::Float64)
    (theta(d) * Lmax) / (popt(d)^(1-sigma)+1)
  end
  function q_m(Lmax::Float64, d::Float64, sigma::Float64)
    (theta(d) * Lmax) * (popt(d)^(-sigma)) / (popt(d)^(1-sigma)+1)
  end

  ######## Numerical solution ##########

  function maxuty(p::Float64, Lmax::Float64, d::Float64, sigma::Float64)
    if gamma1 * Topt - gamma2 * d > 0.0
      m = Model(solver=IpoptSolver())
      @variable(m, q_c >= 0.0)
      @variable(m, q_m >= 0.0)
      @NLobjective(m, Max, (q_c^((sigma-1)/sigma) + q_m^((sigma-1)/sigma))^(sigma/(sigma-1)))
      @NLconstraint(m, q_c + p*q_m == (gamma1 * Topt - gamma2 *d)*Lmax)
      solve(m)
      return res = (getvalue(q_m), getvalue(q_c))
    else
      return (0.0, 0.0)
    end
  end

  ######################################
  ############ EQUILIBRIUM #############
  ######################################

  function RSup(q_c::Float64, Lmax::Float64, d::Float64)
    Q_m(q_c, Lmax, d)/q_c
  end

  function RDem(q_c::Float64, q_m::Float64)
    q_m/q_c
  end




function autarky_eq(Lmax::Float64, t::Float64, sigma::Float64)

  # Define the distance to the optimal temperature
  d = dist(t)

  ######## Analytical solutions #########

  qc = q_c(Lmax, d, sigma)
  qm = q_m(Lmax, d, sigma)
  RD = RDem(qc, qm)
  RS = RSup(qc, Lmax, d)

  ######## Numerical solutions ##########

  p      = popt(d)
  res    = maxuty(p, Lmax, d, sigma)

  # if maxqc == 0.0
  #   num_qc = 0.0
  # else
  #   num_qc = fzeros(q_c -> num_RD - num_RS(q_c), 0.0, maxqc)[1]
  # end

  num_qm = res[1]
  num_qc = res[2]
  num_RD = RDem(num_qc, num_qm)
  num_RS = RSup(num_qc, Lmax, d)


  # Useful values
  Lc       = L_c(qc, Lmax, d)
  num_Lc   = L_c(num_qc, Lmax, d)
  totcrops = theta(d)*Lc
  K        = totcrops - qc
  maxqm    = Q_m(0.0, Lmax, d)


return Dict("Optimal Relative Demand"        => RD,
            "Optimal Relative Demand num"    => num_RD,
            "Optimal price"                  => p,
            "Optimal net cereal consumption" => qc,
            "Optimal net c C° num"           => num_qc,
            "Optimal meat consumption"       => qm,
            "Optimal m C° num"               => num_qm,
            "Percentage of feed"             => (K/totcrops)*100.0,
            "Percentage of land to crop"     => (Lc/Lmax)*100.0,
            "Percentage of land to crop num" => (num_Lc/Lmax)*100.0,
            "Max qm"                         => maxqm)

end #function autarky_eq

###################################################################
# Plot how the equilibrium evolves with country-specific parameters
###################################################################

# We use the analytic solutions because we have them

# set_default_plot_size(20cm, 15cm)
# nevals = 100
#
# ## Land endowment, fixing t=15.0 and sigma = 1.0
# # Land endowment and cereal consumption
# gridland = linspace(1e-4, 1e6, nevals) # setting a grid for land endowment
# yland = [autarky_eq(gridland[i], 15.0, .9) for i in 1:nevals] #corresponding y using function autarky_eq
#
# landplot1 = plot(x = gridland, y = collect(yland[i]["Optimal net cereal consumption"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Land endowment"), Guide.ylabel("Opt q_c"), Guide.title("Fixing t=15.0 and sigma = .9"))
# landplot2 = plot(x = gridland, y = collect(yland[i]["Optimal meat consumption"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Land endowment"), Guide.ylabel("Opt q_m"), Guide.title("Fixing t=15.0 and sigma = .9"))
# landplot3 = plot(x = gridland, y = collect(yland[i]["Optimal Relative Demand"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Land endowment"), Guide.ylabel("Opt RD"), Guide.title("Fixing t=15.0 and sigma = .9"))
# landplot4 = plot(x = gridland, y = collect(yland[i]["Optimal price"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Land endowment"), Guide.ylabel("Opt price"), Guide.title("Fixing t=15.0 and sigma = .9"))
# #landplot5 = plot(x = gridland, y = collect(yland[i]["Percentage of land to crop"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Land endowment"), Guide.ylabel("Per Land to crop"), Guide.title("Fixing t=15.0 and sigma = 1.0"))
# #landplot6 = plot(x = gridland, y = collect(yland[i]["Percentage of feed"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Land endowment"), Guide.ylabel("Per crop to feed"), Guide.title("Fixing t=15.0 and sigma = 1.0"))
#
# landplots = gridstack([landplot1 landplot2; landplot3 landplot4])
#
# ## Temperatures, fixing Lmax=1000 and sigma = 1.0
# gridtemp = linspace(5.0, 25.0, nevals) # setting a grid for land endowment
# ytemp = [autarky_eq(1000.0, gridtemp[i], .9) for i in 1:nevals] #corresponding y using function autarky_eq
#
# tempplot1 = plot(x = gridtemp, y = collect(ytemp[i]["Optimal net cereal consumption"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Temperatures"), Guide.ylabel("Opt q_c"), Guide.title("Fixing Lmax=1000 and sigma = .9"))
# tempplot2 = plot(x = gridtemp, y = collect(ytemp[i]["Optimal meat consumption"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Temperatures"), Guide.ylabel("Opt q_m"), Guide.title("Fixing Lmax=1000 and sigma = .9"))
# tempplot3 = plot(x = gridtemp, y = collect(ytemp[i]["Optimal Relative Demand"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Temperatures"), Guide.ylabel("Opt RD"), Guide.title("Fixing Lmax=1000 and sigma = .9"))
# tempplot4 = plot(x = gridtemp, y = collect(ytemp[i]["Optimal price"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Temperatures"), Guide.ylabel("Opt price"), Guide.title("Fixing Lmax=1000 and sigma = .9"))
# #tempplot5 = plot(x = gridtemp, y = collect(yland[i]["Percentage of land to crop"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Temperatures"), Guide.ylabel("Per Land to crop"), Guide.title("Fixing Lmax=1000 and sigma = 1.0"))
# #tempplot6 = plot(x = gridtemp, y = collect(yland[i]["Percentage of feed"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("Temperatures"), Guide.ylabel("Per crop to feed"), Guide.title("Fixing Lmax=1000 and sigma = 1.0"))
#
# tempplots = gridstack([tempplot1 tempplot2; tempplot3 tempplot4])
#
# # CES, fixing Lmax=1000 and t = 15.0
# gridsigma = linspace(1e-2, 1-1e-2, nevals) # setting a grid for land endowment
# ysigma = [autarky_eq(1000.0, 15.0, gridsigma[i]) for i in 1:nevals] #corresponding y using function autarky_eq
#
# sigmaplot1 = plot(x = gridsigma, y = collect(ysigma[i]["Optimal net cereal consumption"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("Opt q_c"), Guide.title("Fixing Lmax=1000 and t = 15.0"))
# sigmaplot2 = plot(x = gridsigma, y = collect(ysigma[i]["Optimal meat consumption"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("Opt q_m"), Guide.title("Fixing Lmax=1000 and t = 15.0"))
# sigmaplot3 = plot(x = gridsigma, y = collect(ysigma[i]["Optimal Relative Demand"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("Opt RD"), Guide.title("Fixing Lmax=1000 and t = 15.0"))
# sigmaplot4 = plot(x = gridsigma, y = collect(ysigma[i]["Optimal price"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("Opt price"), Guide.title("Fixing Lmax=1000 and t = 15.0"))
# sigmaplot5 = plot(x = gridsigma, y = collect(ysigma[i]["Percentage of land to crop"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("% land to crop"), Guide.title("Fixing Lmax=1000 and t = 15.0"))
# sigmaplot6 = plot(x = gridsigma, y = collect(ysigma[i]["Percentage of feed"] for i in 1:nevals), Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("% crop to feed"), Guide.title("Fixing Lmax=1000 and t = 15.0"))
#
# sigmaplots = gridstack([sigmaplot1 sigmaplot2; sigmaplot3 sigmaplot4; sigmaplot5 sigmaplot6])

function runall(Lmax::Float64, t::Float64, sigma::Float64)
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
# display(landplots)
# display(tempplots)
# display(sigmaplots)
end

end #module
