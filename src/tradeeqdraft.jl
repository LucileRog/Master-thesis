module tradeeq

# module to find the equilibrium quantities and prices with CES, land endowment
# and temperatures as arguments, in case of autarky

using Gadfly
using Roots
using NLopt
using Optim
using Ipopt
# import ApproXD: getBasis, BSpline
# using Distributions, ApproxFun, FastGaussQuadrature


# PARAMETERS
a = 10.0         # coef of conversion of crop into meat
b = 1e-4         # coef of conversion of land into meat
gamma1 = 3800/15 # coef of conversion of land into crop
                 # http://data.worldbank.org/indicator/AG.YLD.CREL.KG
gamma2 = 600     # slope coef of how crop productivity decreases with distance
                 # to the optimal temperature
Topt = 15.0      # Optimal temeprature


  ######################################
  ############### SUPPLY ###############
  ######################################

  # distance to optimal temperature for crops
    function d(t)
      abs(Topt - t)
    end

  # Crop productivity at the country level
    function theta(d::Vector)
      theta = ones(length(d))
      for i in 1:length(d)
        if gamma1 * Topt - gamma2 *d[i] > 0.0
          theta[i] = gamma1 * Topt - gamma2 *d[i]
        else
          theta[i] = 0.0
        end
      end
      return theta
    end

  # True crop productivity at the world level
    function theta_w(L_c::Float64, ct::Matrix)
      Nct = length(ct[:,1])
      res=0.0 # to be replaced
      if L_c <= ct[1,1]
        res = ct[1,2]
      end
      for i in 2:Nct
        if sum(ct[j,1] for j in 1:i-1) <= L_c <= sum(ct[j,1] for j in 1:i)
          res = sum((ct[j,1]/L_c)*ct[j,2] for j in 1:i-1) + ((L_c-sum(ct[j,1] for j in 1:i-1))/L_c)*ct[i,2]
        end
      end
      return res
    end

  # L_c obtained by meat producers equating inputs
    function L_c(q_c::Float64, ct::Matrix)
      Nct = length(ct[:,1])
      Lsum = sum(ct[i,1] for i in 1:Nct)
      fzeros(L_c -> a*(Lsum - L_c) - b*(theta_w(L_c, ct)*L_c - q_c), 1.0, Lsum)[1]
    end
    # 20. seconds

  # Gives Q_m as a function of q_c
    function Q_m(q_c::Float64, ct::Matrix)
      Nct = length(ct[:,1])
      Lsum = sum(ct[i,1] for i in 1:Nct)
      min( a*(Lsum-L_c(q_c, ct)) , b*(theta_w(L_c(q_c, ct), ct)*L_c(q_c, ct)-q_c) )
    end
    # 0.36 seconds

  # Optimal price, from profit = 0
    function popt(q_c::Float64, ct::Matrix)
      theta_w(L_c(q_c, ct), ct)/a + 1/b
    end
    # about 2.7 seconds


  ######################################
  ############### DEMAND ###############
  ######################################

  # We assume that sigma is the same in each country

    ## Analytical ######################

      function RD(q_c::Float64, ct::Matrix)
        sig = ct[1,3]
        popt(q_c, ct) ^ (-sig)
      end

    ## Numerical ######################

      # Demand maximizes utility under the budget constraint
      # Recall that 0 < sigma < 1

      # q[1] = meat
      # q[2] = cereals

      function objfun(q::Vector, grad::Vector, ct::Matrix)
        sig = ct[1,3]
        if length(grad) > 0
          grad[1] = q[1]^(-1/sig) * (q[1]^((sig-1)/sig) + q[2]^((sig-1)/sig))^(1/(sig-1))
          grad[2] = q[2]^(-1/sig) * (q[1]^((sig-1)/sig) + q[2]^((sig-1)/sig))^(1/(sig-1))
        end
        obj = (q[1]^((sig-1)/sig) + q[2]^((sig-1)/sig))^(sig/(sig-1))
        return grad, obj
      end

      function constr(q::Vector, grad::Vector, ct::Matrix)
        Nct = length(ct[:,1])
        if length(grad) > 0
  				 grad[1] = theta_w(L_c(q[2], ct), ct)/a + 1/b
   				grad[2] = 1
   			end
     		constr = q[1]*(theta_w(L_c(q[2], ct),ct)/a + 1/b) + q[2] - sum(ct[i,2]*ct[i, 1] for i in 1:Nct)
     		return grad, constr
      end

      function maxuty(ct::Matrix)
  		    opt = NLopt.Opt(:LD_SLSQP,2)
  		      # define the type optimization
          NLopt.max_objective!(opt,(q,g)->objfun(q,g,ct)[2])
          # define the bounds, note that consumption can not be negative
  		      lower_bounds!(opt,[0.0 ; 0.0])
  		      xtol_rel!(opt,1e-9)
  		      ftol_rel!(opt,1e-9)
  	      # define the constraint
          NLopt.equality_constraint!(opt, (q,g)-> constr(q,g,ct)[2],1e-9)
          (optf,optq,ret) = NLopt.optimize(opt, rand(2))
  		  return optq
  	   end

    # I use NLopt because the multiple occurence of theta_w() in the constraint
    # would have complicated a lot the expression in @NLconstraint


    ######################################
    ############ EQUILIBRIUM #############
    ######################################

    ## Analytical solutions
     # where the analytical part comes from the maximuzation of utility

      function RDRS(q_c::Float64, ct::Matrix)
        RD(q_c,ct) - Q_m(q_c,ct)/q_c
      end
      # 0.68 seconds

      function hRDRS(q_c::Float64, ct::Matrix)
        (RD(q_c,ct) - Q_m(q_c,ct)/q_c)*1e8
      end
      # 0.71 seconds

      function abshRDRS(q_c::Float64, ct::Matrix)
        abs((RD(q_c, ct) - Q_m(q_c,ct)/q_c)*1e8)
      end
      # plot(x = linspace(0.0, maxqc, 100), y = [abshRDRS(linspace(0.0, maxqc, 100)[i]) for i in 1:100], Geom.point)

      ##### Note on the method for finding the equilibrium
      ## The function RD - RS was extremely flat towards zero which made the function
        # fzeros() very slow (about 780 seconds to find the equilibrium)
      ## A first tentative was to raise the "power" of RD - RS by multiplying the
        # result by 1e8, but it did not increase the speed dramatically
      ## The solution I chose was to take the absolute value of the function, so that
        # the minimum would be where the function hit 0, I use the package Optim,
        # and now the equilibrium is found in about 20 seconds


function trade_eq(Lmax::Vector, t::Vector, sigma::Vector)

  ##############################################################################
  ########################### WORLD CHARACTERISTICS ############################
  ##############################################################################

  # Number of countries
    Nct = length(Lmax)
  # Sum of all countries land endowments
    Lsum = sum(Lmax[i] for i in 1:Nct)
  # Distance from optimal temperature for all countries
    dist = collect(d(t[i]) for i in 1:Nct)

  # Sorting countries by decreasing theta(d)
    countries = hcat(Lmax, theta(dist), sigma)
    ct = sortrows(countries, by=x->x[2], rev=true) # dim = Ncountries x 2

  # Define some useful maximal values
    maxqc = theta_w(Lsum, ct)*Lsum
    htheta = ct[1,2] # the maximal tehta because ct are sorted by decreasing theta
    maxpopt = htheta/a + 1/b

  ## Equilibrium quantities

     # Analytical solution
        if Optim.optimize(q_c -> abshRDRS(q_c, ct), 0.0, maxqc).minimum < 1e-3
          qc = Optim.optimize(q_c -> abshRDRS(q_c, ct), 0.0, maxqc).minimizer
          qm = Q_m(qc, ct)
        else
          qc = 0.0
          qm = 0.0
        end

     # Numerical solution
       num_qm = maxuty(ct)[1]
       num_qc = maxuty(ct)[2]

  # Useful values
    Lc       = L_c(qc, ct)
    num_Lc   = L_c(num_qc, ct)
    totcrops = theta_w(Lc, ct)*Lc
    K        = totcrops - qc


  return Dict("Optimal Relative Demand" => RD(qc, ct),
              "Optimal Relative Demand num" => num_qm/num_qc,
              "Optimal price" => popt(qc, ct),
              "Optimal price num" => popt(num_qc, ct),
              "Optimal net cereal consumption" => qc,
              "Optimal net c C° num" => num_qc,
              "Optimal meat consumption" => qm,
              "Optimal m C° num" => num_qm,
              "Percentage of feed" => (K/totcrops)*100.0,
              "Percentage of land to crop" => (Lc/Lsum)*100.0,
              "Percentage of land to crop num" => (num_Lc/Lsum)*100.0)

end # function trade_eq



function runall(Lmax, t, sigma)
  res = trade_eq(Lmax, t, sigma)
  println("-------------------------------------------")
  println("Solutions when Lmax = ", Lmax, "  t = ", t, "  and sigma = ", sigma)
  println("Numerical solution for q_c: ", res["Optimal net c C° num"])
  println("Numerical solution for q_m: ", res["Optimal m C° num"])
  println("Percent of cereals used as feed: ", res["Percentage of feed"])
  println("Percentage of land for crop production: ", res["Percentage of land to crop"])
  println("Optimal relative price: ", res["Optimal price"])
  println("-------------------------------------------")
end


end # end
