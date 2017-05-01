module tradeeq

# module to find the equilibrium quantities and prices with CES, land endowment
# and temperatures as arguments, in case of autarky

using Gadfly
using Roots
using NLopt
using Optim
using Ipopt
using DataFrames

include("autarkyeq.jl")



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
      if length(fzeros(L_c -> a*(Lsum - L_c) - b*(theta_w(L_c, ct)*L_c - q_c), 1.0, Lsum)) == 1
       return fzeros(L_c -> a*(Lsum - L_c) - b*(theta_w(L_c, ct)*L_c - q_c), 1.0, Lsum)[1]
     else
       0.0
     end
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
        Nct = length(ct[:,1])
          if ct[:,1] == zeros(Nct) || ct[:,2] == zeros(Nct)
            return (0.0, 0.0)
          else
    		    opt = NLopt.Opt(:LD_SLSQP,2) # define the type optimization
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
        (RD(q_c,ct) - Q_m(q_c,ct)/q_c)*1e10
      end
      # 0.71 seconds

      function abshRDRS(q_c::Float64, ct::Matrix)
        abs((RD(q_c, ct) - Q_m(q_c,ct)/q_c)*1e10)
      end
      # plot(x = linspace(0.0, maxqc, 100), y = [abshRDRS(linspace(0.0, maxqc, 100)[i]) for i in 1:100], Geom.point)

      ##### Note on the method for finding the equilibrium
      ## For some values of Lmax, t, and sigma,  RD - RS was extremely flat
        # towards zero which made the function fzeros() very slow (about 780
        # seconds to find the equilibrium)
      ## A first tentative was to raise the "power" of RD - RS by multiplying the
        # result by 1e10, but it did not increase the speed dramatically
      ## The solution I chose was to take the absolute value of the function, so that
        # the minimum would be where the function hit 0, I use the package Optim,
        # and now the equilibrium is found in about 20 seconds


# A function to generate a sorted list of countries
    function sortct(Lmax::Vector, t::Vector, sigma::Vector)
      # Number of countries
        Nct = length(Lmax)
      # Sum of all countries land endowments
        Lsum = sum(Lmax[i] for i in 1:Nct)
      # Distance from optimal temperature for all countries
        dist = collect(d(t[i]) for i in 1:Nct)
      # Sorting countries by decreasing theta(d)
        countries = hcat(Lmax, theta(dist), sigma)
        return sortrows(countries, by=x->x[2], rev=true) # dim = Ncountries x 2
    end

function trade_eq(Lmax::Vector, t::Vector, sigma::Vector)

  ##############################################################################
  ########################### WORLD CHARACTERISTICS ############################
  ##############################################################################
  #(Lmax,t,sigma) = ([1200.0, 800.0], [15.0, 14.0], [.9,.9])

  # Number of countries
    Nct = length(Lmax)
  # Sum of all countries land endowments
    Lsum = sum(Lmax[i] for i in 1:Nct)

  # Sorted countries
    ct = sortct(Lmax, t, sigma)

  # Define some useful maximal values
    maxqc   = theta_w(Lsum, ct)*Lsum
    htheta  = ct[1,2] # the maximal tehta because ct are sorted by decreasing theta
    maxpopt = htheta/a + 1/b

  # Equilibrium quantities

      # Analytical solution
        # For some values of Lmax, t and sigma, RDRS() is very flat towards 0,
        # We consruct a gris of qc
        gridqc = linspace(0.0, maxqc, 100)
        # The criterion is that the absolute value of the last 95% of points on
        # the grid are lower than 1e-2
        crit = [[abs(RDRS(gridqc[i], ct)) for i in 5:100 ] .< [1e-2 for i in 5:100]][1]
        vectrue = [true for i in 5:100]
        # if this is true, we apply the method to find the root explained above
         if crit == vectrue
            if Optim.optimize(q_c -> abshRDRS(q_c, ct), 0.0, maxqc).minimum < 1e-1
              qc = Optim.optimize(q_c -> abshRDRS(q_c, ct), 0.0, maxqc).minimizer
              qm = Q_m(qc, ct)
            else
              qc = 0.0
              qm = 0.0
            end
         else
        # otherwise, we simply run the fzeros() function
          if length(fzeros(q_c -> RDRS(q_c, ct), 0.0, maxqc)) == 1
            qc = fzeros(q_c -> RDRS(q_c, ct), 0.0, maxqc)[1]
            qm = Q_m(qc, ct)
          else
            qc = 0.0
            qm = 0.0
          end
        end

        ## Note on time/efficiency
         # I tested the efficiency of a fasttrade_eq() function that would
         # implement this method vs. a slowtrade_eq() that would only use
         # fzeros() and compared their results:
         # for qc : a difference of 1e-2 (qc around 1e6)
         # for qm : a difference of 1e-6 (qm around 1e3)
         # @time fasttrade_eq = 6 sec / @time slowtrade_eq = 53 sec
         # I conclude that it is better to choose the fasttrade_eq because the
         # efficiency loss is rather small compared to the large gain in time.

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
              "Total land to crops" => Lc,
              "Total crop production" => totcrops,
              "Percentage of feed" => (K/totcrops)*100.0,
              "Percentage of land to crop" => (Lc/Lsum)*100.0,
              "Percentage of land to crop num" => (num_Lc/Lsum)*100.0)

end # function trade_eq


function runall(Lmax::Vector, t::Vector, sigma::Vector)
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


function datafr(Lmax::Vector, t::Vector, sigma::Vector)
  # Lmax = [1200.0, 800.0, 1600.0]
  # t = [15.0, 14.0, 17.5]
  # sigma = [.9,.9,.9]

  ct = sortct(Lmax, t, sigma)      # Matrix of countries, sorted (L, theta, sig)
  Nct = length(ct[:,1])
  countries = hcat(Lmax, t, sigma) # Matrix of countries, not sorted (L, t, sig)
  worldres = trade_eq(Lmax, t, sigma)

  # Temperature among sorted countries
    Temp = zeros(Nct)
    for i in 1:Nct
      for j in 1:Nct
        if ct[i,1] == countries[j,1]
          Temp[i] = countries[j,2]
        end
      end
    end

  # Share of lands to crops per country
    worldLc = worldres["Total land to crops"]
    shareLc = zeros(Nct) # to be replaced
    if worldLc <= ct[1,1]
      shareLc[1] = 100.0
    end
    for i in 2:Nct
      if sum(ct[j,1] for j in 1:i-1) <= worldLc <= sum(ct[j,1] for j in 1:i)
        for k in 1:Nct-1
        shareLc[k] = ct[k,1]/worldLc
        end
        shareLc[i] = 1 - sum(shareLc[j] for j in 1:i-1)
      end
    end
    shareLc

  # Share of cereal production
    worldcprod = worldres["Total crop production"]
    sharecprod = zeros(Nct)
    for i in 1:Nct
      sharecprod[i] = ct[i,2]*(shareLc[i]*worldLc)/worldcprod
    end
    sharecprod

  # Share of meat production
    # Recall that Q_m = q_m
    worldmprod = worldres["Optimal meat consumption"]
    sharemprod = zeros(Nct)
    for i in 1:Nct
      sharemprod[i] = a*(ct[i,1]-shareLc[i]*worldLc)
    end
    sharemprod

  # Imports/exports
    # cereals
      worldccons = worldres["Optimal net cereal consumption"]
      tradec = zeros(Nct)
      autqc  = zeros(Nct)
      for i in 1:Nct
        autqc[i] = autarkyeq.autarky_eq(ct[i,1], Temp[i], ct[i,3])["Optimal net cereal consumption"]
        tradec[i] = (1/Nct)*worldccons - autqc[i]
      end
      tradec

    # meat
      worldmcons = worldres["Optimal meat consumption"]
      tradem = zeros(Nct)
      autqm  = zeros(Nct)
      for i in 1:Nct
        autqm[i] = autarkyeq.autarky_eq(ct[i,1], Temp[i], ct[i,3])["Optimal meat consumption"]
        tradem[i] = (1/Nct)*worldmcons - autqm[i]
      end
      tradem

return DataFrame(Land_endowment =ct[:,1],
                 Temperature=Temp,
                 CES=ct[:,3],
                 Crop_yields=ct[:,2],
                 Share_land_to_crops=shareLc,
                 Share_cereal_prod=sharecprod,
                 Share_meat_prod=sharemprod,
                 Trade_in_cereals=tradec,
                 Trade_in_meat=tradem)

end


end # end module
