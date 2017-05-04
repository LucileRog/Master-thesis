
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
    a      = 1e4      # how many kg of meat produced with one unit of land (ha)
    b      = 1/6.5    # how many kg of meat produced with one unit of cereals (kgs)
                      # http://www.nature.com/nature/journal/v418/n6898/full/nature01014.html
    gamma1 = 3800/15  # how many kgs of cereals with one unit of land (ha)
                      # http://data.worldbank.org/indicator/AG.YLD.CREL.KG
    gamma2 = 600      # slope coef of how crop productivity decreases with
                      # distance to optimal temperature
    # Optimal temeprature
    Topt   = 15.0

  ######################################
  ############### SUPPLY ###############
  ######################################

    # distance to optimal temperature for crops
    function dist(t)
      abs(Topt - t)
    end

    # crop productivity
    function theta(d::Float64)
      max(gamma1 * Topt - gamma2 *d, 0.0)
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


end #module
