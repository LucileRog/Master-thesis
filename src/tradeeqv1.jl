module tradeeq

# module to find the equilibrium quantities and prices with CES, land endowment
# and temperatures as arguments, in case of autarky

using Gadfly
using Roots
using NLopt
using Optim
using Ipopt


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
  ############# COUNTRIES ##############
  ######################################

  # A function to generate a sorted list of countries
      function sortct(Lmax::Vector, t::Vector, sigma::Vector)
        # Number of countries
          Nct = length(Lmax)
        # Sum of all countries land endowments
          Lsum = sum(Lmax[i] for i in 1:Nct)
        # Distance from optimal temperature for all countries
          dist = collect(d(t[i]) for i in 1:Nct)
        # Sorting countries by decreasing theta(d)
          countries = hcat(Lmax, theta(dist), sigma, t, dist)
          ct = sortrows(countries, by=x->x[2], rev=true) # 1st sort by decreasing thetas
      end


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


  # Numerical expression for theta_w

    function obj_thetaw(L_cR::Float64, L_c::Float64, ct::Matrix)
      -( (L_cR/L_c)*ct[1,2] + (1-L_cR/L_c)*ct[2,2] )
    end

    function num_theta_w(L_c::Float64, ct::Matrix)
      Nct = length(ct[:,1])
      if Nct == 2
        res = Optim.optimize(x -> obj_thetaw(x, L_c, ct), 0.0, min(ct[1,1], L_c))
        return ([res.minimizer, L_c-res.minimizer], -res.minimum)
      elseif Nct > 2
        m = Model(solver=IpoptSolver())
        @variable(m, 0.0 <= Lc[i=1:Nct] <= ct[i,1])
        @NLobjective(m, Max, sum((Lc[i]/L_c)*ct[i,2] for i in 1:Nct))
        @NLconstraint(m, sum(Lc[i] for i in 1:Nct) == L_c)
        solve(m)
        return res = (getvalue(Lc), getobjectivevalue(m))
      end
    end

  # True crop productivity at the world level
    # function theta_w(L_c::Float64, ct::Matrix)
    #   Nct = length(ct[:,1]) # countries are sorted by decreasing theta
    #   res=0.0 # to be replaced
    #   alldiff = [ct[i,2] != ct[i+1,2] for i in 1:Nct-1]
    #   # for 2 countries
    #   if Nct == 2
    #     if alldiff == trues(Nct-1) # if all thetas are different
    #       if L_c <= ct[1,1]
    #         res = ct[1,2] # theta_w = highest theta
    #       elseif ct[1,1] <= L_c <= ct[1,1] + ct[2,1]
    #         res = (ct[1,1]/L_c)*ct[1,2] + (1-ct[1,1]/L_c)*ct[2,2] # weighted sum
    #       end
    #     else # if both countries have same theta, theta_w = theta_1 = theta_2
    #       res = ct[1,2] # or ct[2,2]
    #     end # same theta
    #   # If there are more than 2 countries
    #   else
    #     if alldiff == trues(Nct-1) # if all thetas are different
    #       if L_c <= ct[1,1]
    #         res = ct[1,2] # theta_w = highest theta
    #       end
    #       for i in 2:Nct # if not, a weighted sum of thetas where highest theta has highest share
    #         if sum(ct[j,1] for j in 1:i-1) <= L_c <= sum(ct[j,1] for j in 1:i)
    #           res = sum((ct[j,1]/L_c)*ct[j,2] for j in 1:i-1) + ((L_c-sum(ct[j,1] for j in 1:i-1))/L_c)*ct[i,2]
    #         end
    #       end
    #     end # all thetas are different
    #   end
    #   return res
    # end

  # L_c obtained by meat producers equating inputs
    function L_c(q_c::Float64, ct::Matrix)
      Nct = length(ct[:,1])
      Lsum = sum(ct[i,1] for i in 1:Nct)
        if length(fzeros(L_c -> a*(Lsum - L_c) - b*(num_theta_w(L_c, ct)[2]*L_c - q_c), 1.0, Lsum)) == 1
          return fzeros(L_c -> a*(Lsum - L_c) - b*(num_theta_w(L_c, ct)[2]*L_c - q_c), 1.0, Lsum)[1]
        else
         0.0
        end
    end

    function num_shares(q_c::Float64, ct::Matrix)
      Nct = length(ct[:,1])
      Lsum = sum(ct[i,1] for i in 1:Nct)
      Lcw = L_c(q_c, ct)
      Lmw = Lsum - Lcw
      Lc = num_theta_w(L_c(q_c, ct), ct)[1]
      shLc = zeros(Nct)
      shLm = zeros(Nct)
      for i in 1:Nct
        shLc[i] = Lc[i]/Lcw
        shLm[i] = (ct[i,1] - Lc[i])/Lmw
      end
      return (shLc, shLm)
    end

    # function returning land allocation
    # function shares(q_c::Float64, ct::Matrix)
    #   # recall countries are first sorted by decreasing thetas, then by incresing
    #   # land endowments
    #   Nct = length(ct[:,1])
    #   Lsum = sum(ct[i,1] for i in 1:Nct)
    #   worldLc = L_c(q_c, ct)
    #   worldLm = Lsum - L_c(q_c, ct)
    #   sh = zeros(Nct, 2) # to be replaced (sharesLc, sharesQm)
    #   alldiff = [ct[i,2] != ct[i+1,2] for i in 1:Nct-1]
    #   # for 2 countries
    #   if Nct == 2
    #     if alldiff == trues(Nct-1) # if all thetas are different
    #       if worldLc <= ct[1,1]       # tot lands required for crops < Lmax of highest theta
    #         sh[1,1] = 1.0         # only highest theta's land used for crops
    #         sh[1,2] = (ct[1,1]-worldLc)/worldLm
    #       else                        # if highest theta's land not enough
    #         sh[1,1] = ct[1,1]/worldLc      # highest theta produces as much crops as possible
    #         sh[2,1] = 1 - ct[1,1]/worldLc  # lowest theta produces the rest
    #       end
    #     else # if both countries have same theta
    #       if ct[1,1] == ct[2,1] # if countries are symmetric
    #         sh[:,1] = [1/2, 1/2] # equal share for each type of good
    #       else # if countries have different land endowments
    #         if worldLc/2 <= ct[1,1] # if lowest land country has more than half the resuired land for crops
    #           sh[:,1] = [1/2, 1/2]  # land to crops equally allocated across countries
    #           for i in 1:2
    #             sh[i,2] = (ct[i,1]/2)/worldLm # the rest of each country goes to meat
    #           end
    #         else # if half the required lands for crops is larger than the smallest land endowment
    #           sh[1,1] = ct[1,1]/worldLc # smallest land country takes as much as it can
    #           sh[2,1] = 1.0 - sh[1,1]
    #         end
    #       end
    #     end # same theta
    #   end
    #   for i in 1:Nct
    #     sh[i,2] = (ct[i,1]-sh[i,1]*worldLc)/worldLm
    #   return sh
    # end

  # Optimal price, from profit = 0
    function popt(q_c::Float64, ct::Matrix)
      Nct = length(ct[:,1])
      sharem = num_shares(q_c, ct)[2]
      for i in 1:Nct
        return (1/a)*sum(sharem[i]*ct[i,2]) + 1/b
      end
    end
    # about 2.7 seconds


  ######################################
  ############### DEMAND ###############
  ######################################

  # We assume that sigma is the same in each country

      # Demand maximizes utility under the budget constraint
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
  				grad[1] = popt(q[2], ct)
   				grad[2] = 1
   			end
     		constr = q[1]*popt(q[2], ct) + q[2] - sum(ct[i,2]*ct[i, 1] for i in 1:Nct)
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



function trade_eq(Lmax::Vector, t::Vector, sigma::Vector)

  ##############################################################################
  ########################### WORLD CHARACTERISTICS ############################
  ##############################################################################

  # Number of countries
    Nct = length(Lmax)
  # Sum of all countries land endowments
    Lsum = sum(Lmax[i] for i in 1:Nct)

  # Sorted countries
    ct = sortct(Lmax, t, sigma)

  # Define some useful maximal values
    maxqc   = num_theta_w(Lsum, ct)[2]*Lsum
    htheta  = ct[1,2] # the maximal tehta because ct are sorted by decreasing theta
    maxpopt = htheta/a + 1/b

  # Equilibrium quantities
    qm = maxuty(ct)[1]
    qc = maxuty(ct)[2]

  # Useful values
    Lc       = L_c(qc, ct)
    totcrops = num_theta_w(Lc, ct)[2]*Lc
    K        = totcrops - qc


  return Dict("Number of countries"            => Nct,
              "Optimal Relative Demand"        => qm/qc,
              "Optimal price"                  => popt(qc, ct),
              "Optimal net cereal consumption" => qc,
              "Optimal meat consumption"       => qm,
              "Total land to crops"            => Lc,
              "Total crop production"          => totcrops,
              "Percentage of feed"             => (K/totcrops)*100.0,
              "Percentage of land to crop"     => (Lc/Lsum)*100.0)
end # function trade_eq



end # end module
