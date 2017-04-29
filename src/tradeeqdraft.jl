module tradeeqdraft

using Gadfly
using Roots
using NLopt

# PARAMETERS
a = 10.0         # coef of conversion of crop into meat
b = 1e-4         # coef of conversion of land into meat
gamma1 = 3800/15 # coef of conversion of land into crop
                 # http://data.worldbank.org/indicator/AG.YLD.CREL.KG
gamma2 = 600     # slope coef of how crop productivity decreases with distance
                 # to the optimal temperature
Topt = 15.0      # Optimal temeprature

Lmax = [1200.0, 800.0, 1600.0]
t = [15.0, 14.0, 17.5]
sigma = [.9,.9,.9]


function num_tradeeq(Lmax::Vector, t::Vector, sigma::Vector)
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



      function objfun(q::Vector, grad::Vector)
        sig = sigma[1]
        if length(grad) > 0
          grad[1] = q[1]^(-1/sig) * (q[1]^((sig-1)/sig) + q[2]^((sig-1)/sig))^(1/(sig-1))
          grad[2] = q[2]^(-1/sig) * (q[1]^((sig-1)/sig) + q[2]^((sig-1)/sig))^(1/(sig-1))
        end
        obj = (q[1]^((sig-1)/sig) + q[2]^((sig-1)/sig))^(sig/(sig-1))
        return grad, obj
      end

      function constr(q::Vector, grad::Vector)
        if length(grad) > 0
            grad[1] = theta_w(L_c(q[2]))/a + 1/b
            grad[2] = 1
        end
        constr = q[1]*(theta_w(L_c(q[2]))/a + 1/b) + q[2] - sum((gamma1*Topt - gamma2*dist[i])*Lmax[i] for i in 1:Nct)
        return grad, constr
      end

      function maxuty()
        opt = NLopt.Opt(:LD_SLSQP,2)
          # define the type optimization
        NLopt.max_objective!(opt,(q,g)->objfun(q,g)[2])
        # define the bounds, note that consumption can not be negative
          NLopt.lower_bounds!(opt,[0.0 ; 0.0])
          NLopt.xtol_rel!(opt,1e-9)
          NLopt.ftol_rel!(opt,1e-9)
        # define the constraint
        NLopt.equality_constraint!(opt, (q,g)-> constr(q,g)[2],1e-9)
        (optf,optq,ret) = NLopt.optimize(opt, rand(2))
          return optq
      end



num_qm = maxuty()[1]
num_qc = maxuty()[2]

num_Lc   = L_c(num_qc)
totcrops = theta_w(num_Lc)*num_Lc
K        = totcrops - num_qc


return Dict("Optimal Relative Demand num" => num_qm/num_qc,
            "Optimal price num" => popt(num_qc),
            "Optimal net c C° num" => num_qc,
            "Optimal meat consumption num" => num_qm,
            "Optimal m C° num" => num_qm,
            "Percentage of feed" => (K/totcrops)*100.0,
            "Percentage of land to crop num" => (num_Lc/Lsum)*100.0)

end







end
