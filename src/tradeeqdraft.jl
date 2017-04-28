module tradeeqdraft

using Gadfly
using Roots
using Distributions, ApproxFun, FastGaussQuadrature

# PARAMETERS
a = 10.0         # coef of conversion of crop into meat
b = 1e-4         # coef of conversion of land into meat
gamma1 = 3800/15 # coef of conversion of land into crop
                 # http://data.worldbank.org/indicator/AG.YLD.CREL.KG
gamma2 = 600     # slope coef of how crop productivity decreases with distance
                 # to the optimal temperature
Topt = 15.0      # Optimal temeprature

function theta_w(L_c, Lmax, t, sigma)

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

# ordering countries by theta(d)
countries = hcat(Lmax, theta(dist))
ct = sortrows(countries, by=x->x[2], rev=true) # dim = Ncountries x 2
N = length(Lmax) #Ncountries

  # countries are characterized by a lend endowment and a crop productivity
  res=0.0 # to be replaced
  if L_c <= ct[1,1]
    res = ct[1,2]
  end
  for i in 2:N
    if sum(ct[j,1] for j in 1:i-1) <= L_c <= sum(ct[j,1] for j in 1:i)
      res = sum((ct[j,1]/L_c)*ct[j,2] for j in 1:i-1) + ((L_c-sum(ct[j,1] for j in 1:i-1))/L_c)*ct[i,2]
    end
  end
  return res

end

function RDRS(q_c, Lmax, t, sigma)
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

  # rdering countries by theta(d)
  countries = hcat(Lmax, theta(dist))
  ct = sortrows(countries, by=x->x[2], rev=true) # dim = Ncountries x 2
  N = length(Lmax) #Ncountries

  function theta_w(L_c::Float64)
    # countries are characterized by a lend endowment and a crop productivity
    res=0.0 # to be replaced
    if L_c <= ct[1,1]
      res = ct[1,2]
    end
    for i in 2:N
      if sum(ct[j,1] for j in 1:i-1) <= L_c <= sum(ct[j,1] for j in 1:i)
        res = sum((ct[j,1]/L_c)*ct[j,2] for j in 1:i-1) + ((L_c-sum(ct[j,1] for j in 1:i-1))/L_c)*ct[i,2]
      end
    end
    return res
  end

  # L_c obtained by meat producers equating inputs
  function L_c(q_c::Float64)
  fzeros(L_c -> a*(Lsum - L_c) - b*(theta_w(L_c)*L_c - q_c), 0.0, Lsum)[1]
  end

  maxqc = theta_w(Lsum)*Lsum
  htheta = ct[1,2] #findmax(theta(dist))[1]

  # Optimal price, from profit = 0
  function popt(q_c::Float64)
    popt = fzeros(p -> p - theta_w(L_c(q_c))/a - 1/b, 0.0, htheta/a + 1/b)[1]
  end

  # Relative demand with a CES utility function
  function RD(q_c)
    popt(q_c) ^ (-sigma[1])
  end

  # Gives Q_m
  function Q_m(q_c::Float64)
    a * (Lsum - L_c(q_c))
  end

  # Optimal q_c from RD = RS
  # By plotting q_c -> RD(q_c) - Q_m(q_c)/q_c, we notice that the function is
  # extremely flat towards 0, which slows a lot the computation
    return  RD(q_c) - Q_m(q_c)/q_c

end

function plotall(Lmax, t, sigma)

Lsum = sum(Lmax[i] for i in 1:length(Lmax))
gridLc = linspace(0.0, Lsum, 100)
ythetaw = [theta_w(gridLc[i], Lmax, t, sigma) for i in 1:100]
plottheta = plot(x=gridLc, y=ythetaw, Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("L_c,w"), Guide.ylabel("Theta_w"), Guide.title("Theta_w : world crop productivity"))

gridqc = linspace(0.0, 1e6, 100)
yRDRS = [RDRS(gridqc[i], Lmax, t, sigma) for i in 1:100]
plotRDRS = plot(x=gridqc, y=yRDRS, Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("q_c,w"), Guide.ylabel("RD - RS"), Guide.title("RD = RS, flat curve"))
display(plottheta)
display(plotRDRS)
end


end
