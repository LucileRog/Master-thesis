module tradeeq

# module to find the equilibrium quantities and prices with CES, land endowment
# and temperatures as arguments, in case of autarky

using Gadfly
using Roots
using ApproxFun
using JuMP
using Optim
using AmplNLWriter, CoinOptServices
import ApproXD: getBasis, BSpline
using Distributions, ApproxFun, FastGaussQuadrature


# PARAMETERS
a = 10.0         # coef of conversion of crop into meat
b = 1e-4         # coef of conversion of land into meat
gamma1 = 3800/15 # coef of conversion of land into crop
                 # http://data.worldbank.org/indicator/AG.YLD.CREL.KG
gamma2 = 600     # slope coef of how crop productivity decreases with distance
                 # to the optimal temperature
Topt = 15.0      # Optimal temeprature


function trade_eq(Lmax::Vector, t::Vector, sigma::Vector)

  # For tests
  Lmax = [1000.0, 700.0, 1400.0]
  t = [14.5, 15.2, 16.5]
  sigma = [.9,.9,.9]


# Sum of all countries land endowments
  Lsum = sum(Lmax[i] for i in 1:length(Lmax))

# in the 2 country case, each vector Lmax, t and sigma is of dimension 2

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

  # For the trade equilibrium, we have an analytical solution only for theta
  # L_c and Q_m are computed numerically in each case

  ##############################################################################
  ########################### Analytical solutions #############################
  ##############################################################################

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
  plot(theta_w, 0.0, Lsum)

  # Define some maximal values
  maxqc = theta_w(Lsum)*Lsum
  htheta = ct[1,2] # the maximal tehta because ct are sorted by decreasing theta
  maxpopt = htheta/a + 1/b

  # L_c obtained by meat producers equating inputs
  function L_c(q_c::Float64)
  fzeros(L_c -> a*(Lsum - L_c) - b*(theta_w(L_c)*L_c - q_c), 1.0, Lsum)[1]
  end
  # about 20 seconds

  # Optimal price, from profit = 0
  function popt(q_c::Float64)
    popt = fzeros(p -> p - theta_w(L_c(q_c))/a - 1/b, 0.0, maxpopt)[1]
  end
  # about 2.7 seconds

  # Relative demand with a CES utility function
  function RD(q_c)
    popt(q_c) ^ (-sigma[1])
  end
  # about 0.4 seconds

  # Gives Q_m
  function Q_m(q_c::Float64)
    min( a*(Lsum-L_c(q_c)) , b*(theta_w(L_c(q_c))*L_c(q_c)-q_c) )
  end
  # 0.36 seconds


  # Optimal q_c from RD = RS
  # By plotting q_c -> RD(q_c) - Q_m(q_c)/q_c, we notice that the function is
  # extremely flat towards 0, which slows a lot the computation
  function RDRS(q_c::Float64)
    RD(q_c) - Q_m(q_c)/q_c
  end
  # 0.68 seconds

  function hRDRS(q_c::Float64)
    (RD(q_c) - Q_m(q_c)/q_c)*1e8
  end
  # 0.71 seconds
  plot(x = linspace(0.0, maxqc, 100), y = [hRDRS(linspace(0.0, maxqc, 100)[i]) for i in 1:100], Geom.point)
  # 70 seconds

  qc = fzeros(x -> hRDRS(x), 0.0, maxqc)
  #

  function abshRDRS(q_c::Float64)
    abs((RD(q_c) - Q_m(q_c)/q_c)*1e8)
  end

  qc = optimize(abshRDRS, 0.0, LBFGS(), Optim.Options(autodiff = true))
  #

  # Approximate the function
  # function approx_RDRS(q_c::Float64)
  #   ub,lb = (maxqc, 0.0)
	#   nknots = 13
  #   deg = 3
	#   bs1 = BSpline(nknots,deg,lb,ub)	# equally spaced knots in [lb,ub]
	#   nevals = 5 * bs1.numKnots # get nBasis < nEvalpoints
	# # scaled knots
  #   G(k,s) = GeneralizedPareto(k,s,0)
  #   pf(k,s) = quantile(GeneralizedPareto(k,s,0),linspace(0.05,cdf(G(0.5,1),5),6))
	#   myknots = vcat(-reverse(pf(0.5,1)),0.0,pf(0.5,1))
  #   bs2 = BSpline(myknots,deg)
	# # get coefficients
	#   eval_points = collect(linspace(lb,ub,nevals))
	# #c1 = getBasis(eval_points,bs1) \ runge(eval_points)
	#   c2 = getBasis(eval_points,bs2) \ [RDRS(eval_points[i]) for i in 1:length(eval_points)]
  # return getBasis([q_c],bs1) * c2
  # end

  if length(fzeros(q_c -> RDRS(q_c), 0.0, maxqc)) == 1
  qc = fzeros(q_c -> RDRS(q_c), 0.0, maxqc)[1]
  end
  if length(fzeros(q_c -> RDRS(q_c), 0.0, maxqc)) == 0
  qc = 0.0
  end

  # Optimal q_m
  qm = num_Qm(qc)

  ##############################################################################
  ########################### Numerical solutions ##############################
  ##############################################################################

  # Approximation of crop productivity at the world level
  # WARNING : so far only defined for a world with 3 countries
  function approx_thetaw(L_c)
    ub,lb = (Lsum, 1.0)
    deg = 3
    # Countries
    countries = hcat(Lmax, theta(dist))
    ct = sortrows(countries, by=x->x[2], rev=true) # dim = Ncountries x 2
    N = length(Lmax) #Ncountries
    # myknots with knot multiplicity at 0
    kink = ones(N-1)
      for n in 1:N-1
        kink[n] = sum(ct[i,1] for i in 1:n)
      end
    nknots = 5
    linsp = zeros(N,nknots)
      linsp[1,:] = linspace(1e-3, ct[1,1], nknots)
      for n in 2:N
        linsp[n,:] = linspace(sum(ct[i,1] for i in 1:n-1), sum(ct[i,1] for i in 1:n), nknots)
      end
    myknots = vcat(linsp[1,:], kink[1], linsp[2,:], kink[2], linsp[3,:])
    knots = sort(myknots, rev=false)
    # get coefficients for each case
    params = BSpline(knots,deg)
    nevals = 5 * params.numKnots
    eval_points = collect(linspace(lb,ub,nevals))
    c = getBasis(eval_points,params) \ [theta_w(eval_points[i]) for i in 1:nevals]
    return (getBasis([L_c],params)*c)[1]
  end

  grid = linspace(1e-3, Lsum, 100)
  plot(x=grid, y=[maxtheta(grid[i]) for i in 1:100], Geom.point,
    Theme(default_color=colorant"darkblue", default_point_size=1.5pt,
    highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("L_c,w"),
    Guide.ylabel("MAximal theta"), Guide.title("Numerical solution for theta(L_c)"))
  # We notice that the function shows kinks each time cereals start to be
  # produced in a new country => ApproXD package and use of Bsplines


  Lc       = L_c(q_c(sigma))
  totcrops = theta_w(Lc)*Lc
  K        = totcrops - q_c(sigma)

  figure_eq = plot(xintercept=[RD(sigma)], yintercept=[popt], Geom.vline, Geom.hline)

return Dict("Optimal Relative Demand" => RD(sigma),
            "Optimal price" => popt,
            "Max theta" => theta_w(Lc),
            "Max theta num" => approx_maxtheta(Lc),
            "Optimal net cereal consumption" => q_c(sigma),
            "Optimal meat consumption" => q_m(sigma),
            "Optimal land for crops" => Lc,
            "Optimal meat prod" => Q_m(q_c(sigma)),
            "Percentage of feed" => (K/totcrops)*100.0)

end # function trade_eq


function runall()
  Lmax = [1000.0, 800.0, 500.0]
  t = [15.0, 14.0, 17.5]
  sigma = ones(length(Lmax))
  res = trade_eq(Lmax, t, sigma)
  println("Solution for q_c: ", res["Optimal net cereal consumption"])
  println("Solution for q_m: ", res["Optimal meat consumption"])
  println("Percent of cereals used as feed: ", res["Percentage of feed"])
  #display(landplots)
  #display(tempplots)
  #display(sigmaplots)
end

end
