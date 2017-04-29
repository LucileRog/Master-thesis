# module to find the equilibrium quantities and prices with CES, land endowment
# and temperatures as arguments, in case of autarky

using Gadfly
using Roots
using NLopt
using Ipopt
# import ApproXD: getBasis, BSpline
using Distributions, ApproxFun, FastGaussQuadrature


# PARAMETERS
a = 10.0         # coef of conversion of crop into meat
b = 1e-4         # coef of conversion of land into meat
gamma1 = 3800/15 # coef of conversion of land into crop
                 # http://data.worldbank.org/indicator/AG.YLD.CREL.KG
gamma2 = 600     # slope coef of how crop productivity decreases with distance
                 # to the optimal temperature
Topt = 15.0      # Optimal temeprature

# For tests
Lmax = [1000.0, 700.0, 1400.0]
t = [14.5, 15.2, 16.5]
sigma = [.9,.9,.9]

##############################################################################
########################### WORLD CHARACTERISTICS ############################
##############################################################################

# Sum of all countries land endowments
  Lsum = sum(Lmax[i] for i in 1:length(Lmax))

# distance to optimal temperature for crops
  function d(t)
    abs(Topt - t)
  end

dist = collect(d(t[i]) for i in 1:length(t))

# crop productivity
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

  # For theta_w, we only have an analytical expression, for which we coded
  # the approximation but it is too slow to be used

  # Sorting countries by decreasing theta(d)
  countries = hcat(Lmax, theta(dist))
  ct = sortrows(countries, by=x->x[2], rev=true) # dim = Ncountries x 2
  Nct = length(Lmax) #Ncountries

  # True crop productivity at the world level
  function theta_w(L_c::Float64)
    if L_c <= ct[1,1]
      return ct[1,2]
    end
    for i in 2:Nct
      if sum(ct[j,1] for j in 1:i-1) <= L_c <= sum(ct[j,1] for j in 1:i)
        return sum((ct[j,1]/L_c)*ct[j,2] for j in 1:i-1) + ((L_c-sum(ct[j,1] for j in 1:i-1))/L_c)*ct[i,2]
      end
    end
  end

######################################
############### SUPPLY ###############
######################################

# Define some useful maximal values
maxqc = theta_w(Lsum)*Lsum
htheta = ct[1,2] # the maximal tehta because ct are sorted by decreasing theta
maxpopt = htheta/a + 1/b

# L_c obtained by meat producers equating inputs
  function L_c(q_c::Float64)
  fzeros(L_c -> a*(Lsum - L_c) - b*(theta_w(L_c)*L_c - q_c), 1.0, Lsum)[1]
  end
  # 20. seconds

# Gives Q_m as a function of q_c
  function Q_m(q_c::Float64)
    min( a*(Lsum-L_c(q_c)) , b*(theta_w(L_c(q_c))*L_c(q_c)-q_c) )
  end
  # 0.36 seconds

# Optimal price, from profit = 0
  function popt(q_c::Float64)
    theta_w(L_c(q_c))/a + 1/b
  end
  # about 2.7 seconds


######################################
############### DEMAND ###############
######################################



    ############################################################################

      function maxuty(p::Float64)
        sig = sigma[1]
        if sum((gamma1 * Topt - gamma2 * dist[i]) for i in 1:Nct) == 0.0
          return (0.0, 0.0)
        else
          ############## PROBLEM #############
          # ne fonctionne que pour sigma <.3, sinon too high power
          # mais a déjà fonctionné correctement avec sigma = .1 et below
          m = Model(solver=IpoptSolver())
          @variable(m, q[i=1:2] >= 0.0)
          @NLobjective(m, Max, (q[1]^((sig-1)/sig) + q[2]^((sig-1)/sig))^(sig/(sig-1)))
          @NLconstraint(m, == sum((gamma1 * Topt - gamma2 *dist[i])*Lmax[i] for i in 1:Nct))
          solve(m)
          return res = (getvalue(q))
        end
      end

      function num_q_m(q_c::Float64)
        if isnan(maxuty(popt(q_c))[1]) == false
          return maxuty(popt(q_c))[1]
        else
          return 0.0
        end
      end

      function num_q_c(q_c::Float64)
        if isnan(maxuty(popt(q_c))[2]) == false
          return maxuty(popt(q_c))[2]
        else
          return 0.0
        end
      end

      function num_RD(q_c::Float64)
        if isnan(num_q_c(q_c)) == false
          return num_q_c(q_c)/num_q_m(q_c)
        else
          return 0.0
        end
      end


######################################
############ EQUILIBRIUM #############
######################################


  function num_RDRS(q_c::Float64)
    num_RD(q_c) - Q_m(q_c)/q_c
  end
  gridqc = linspace(0.0, maxqc, 100)
  plot(x = gridqc, y = [num_RDRS(gridqc[i]) for i in 1:100], Geom.point)

  function num_hRDRS(q_c::Float64)
    (num_RD(q_c) - Q_m(q_c)/q_c)*1e8
  end

  function num_abshRDRS(q_c::Float64)
    abs((num_RD(q_c) - Q_m(q_c)/q_c)*1e8)
  end
  plot(x = linspace(0.0, maxqc, 100), y = [num_abshRDRS(linspace(0.0, maxqc, 100)[i]) for i in 1:100], Geom.point)

    if optimize(num_abshRDRS, 0.0, maxqc).minimum < 1e-3
      num_qc = optimize(num_abshRDRS, 0.0, maxqc).minimizer
      num_qm = Q_m(qc)
    else
      num_qc = 0.0
      num_qm = 0.0
    end

    ##############################################################################
    ####################### Approximation of theta(L_c) ##########################
    ##############################################################################

    # Approximation of crop productivity at the world level
    # WARNING : so far only defined for a world with 3 countries
    # function approx_thetaw(L_c)
    #   ub,lb = (Lsum, 1.0)
    #   deg = 3
    #   # Countries
    #   countries = hcat(Lmax, theta(dist))
    #   ct = sortrows(countries, by=x->x[2], rev=true) # dim = Ncountries x 2
    #   N = length(Lmax) #Ncountries
    #   # myknots with knot multiplicity at 0
    #   kink = ones(N-1)
    #     for n in 1:N-1
    #       kink[n] = sum(ct[i,1] for i in 1:n)
    #     end
    #   nknots = 5
    #   linsp = zeros(N,nknots)
    #     linsp[1,:] = linspace(1e-3, ct[1,1], nknots)
    #     for n in 2:N
    #       linsp[n,:] = linspace(sum(ct[i,1] for i in 1:n-1), sum(ct[i,1] for i in 1:n), nknots)
    #     end
    #   myknots = vcat(linsp[1,:], kink[1], linsp[2,:], kink[2], linsp[3,:])
    #   knots = sort(myknots, rev=false)
    #   # get coefficients for each case
    #   params = BSpline(knots,deg)
    #   nevals = 5 * params.numKnots
    #   eval_points = collect(linspace(lb,ub,nevals))
    #   c = getBasis(eval_points,params) \ [theta_w(eval_points[i]) for i in 1:nevals]
    #   return (getBasis([L_c],params)*c)[1]
    # end


Lc       = L_c(qc)
totcrops = theta_w(Lc)*Lc
K        = totcrops - qc

figure_eq = plot(xintercept=[RD(sigma)], yintercept=[popt], Geom.vline, Geom.hline)

return Dict("Optimal Relative Demand" => RD(sigma),
          "Optimal price" => popt,
          "Max theta" => theta_w(Lc),
          "Optimal net cereal consumption" => q_c(popt),
          "Optimal meat consumption" => q_m(popt),
          "Optimal land for crops" => Lc,
          "Optimal meat prod" => Q_m(q_c(sigma)),
          "Percentage of feed" => (K/totcrops)*100.0)
