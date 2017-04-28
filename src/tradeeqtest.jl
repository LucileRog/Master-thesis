using Gadfly
using Roots
using JuMP
using AmplNLWriter, CoinOptServices
import ApproXD: getBasis, BSpline
using Distributions, ApproxFun, FastGaussQuadrature


a = 10.0         # coef of conversion of crop into meat
b = 1e-4         # coef of conversion of land into meat
gamma1 = 3800/15 # coef of conversion of land into crop
                 # http://data.worldbank.org/indicator/AG.YLD.CREL.KG
gamma2 = 600     # slope coef of how crop productivity decreases with distance
                 # to the optimal temperature
Topt = 15.0      # Optimal temeprature

Lmax = [1300.0, 100.0, 400.0]
t = [14.2, 15.3, 12.7]
sigma = ones(3)
Nct   = length(Lmax)
Lsum  = sum(Lmax[i] for i in 1:Nct)

# in the 2 country case, each vector Lmax, t and sigma is of dimension 2

# distance to optimal temperature for crops
function d(t)
  abs(Topt - t)
end

dist = collect(d(t[i]) for i in 1:Nct)

# crop productivity
function theta(d::Vector)
  theta = ones(Nct)
  for i in 1:Nct
    if gamma1 * Topt - gamma2 *d[i] > 0.0
      theta[i] = gamma1 * Topt - gamma2 *d[i]
    else
      theta[i] = 0.0
    end
  end
  return theta
end

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

# Says infeasible, probably because there are not as many constraints as
# variables, but the result seems intuitive although I can't find a mathematical
# proof for it.
function maxtheta(L_c::Float64)
  m = Model()
  @variable(m, 0.0 <= L[i=1:Nct] <= Lmax[i])
  @objective(m, Max, (1/L_c)*sum(L[i]*theta(dist)[i] for i in 1:Nct))
  @constraint(m, sum(L[i] for i in 1:Nct) == L_c)
  solve(m)
  res = (getvalue(L), getobjectivevalue(m))
  return res[2]
end

grid = linspace(1, Lsum, 100)
plot(x=grid, y=[maxtheta(grid[i]) for i in 1:100], Geom.point,
  Theme(default_color=colorant"darkblue", default_point_size=1.5pt,
  highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("L_c,w"),
  Guide.ylabel("Maximal theta"), Guide.title("Numerical solution for theta(L_c)"))
# We notice that the function shows kinks each time cereals start to be
# produced in a new country => ApproXD package and use of Bsplines

function approx_maxtheta(L_c)
  # compare 2 knot vectors with runge's function
  ub,lb = (Lsum, 1e-3)
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
  #myknots = vcat(collect(vcat(vcat(linsp[n,:], kink[n]) for n in 1:N), linspace(ct[N,1], Lsum, nknots)))
    ############## PROBLEM #############
    # je n'arrive pas à faire ça en loop sur les n
  myknots = vcat(linsp[1,:], kink[1], linsp[2,:], kink[2], linsp[3,:])
  knots = sort(myknots, rev=false)
  params = BSpline(knots,deg)
  nevals = 5 * params.numKnots # get nBasis < nEvalpoints

  # get coefficients for each case
  eval_points = collect(linspace(lb,ub,nevals))
    ############## PROBLEM #############
    # obtient que des NAN
  c = getBasis(eval_points,params) \ [theta_w(eval_points[i]) for i in 1:nevals]
  return (getBasis([L_c],params)*c)[1]
end

grid = linspace(1, Lsum, 100)
plot(x=grid, y=[approx_maxtheta(grid[i]) for i in 1:100], Geom.point,
  Theme(default_color=colorant"darkblue", default_point_size=1.5pt,
  highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("L_c,w"),
  Guide.ylabel("Maximal theta"), Guide.title("Numerical solution for theta(L_c)"))
