module subsistconst

# module to find check whether the subsistence consumption constraint is
# satisfied or not for autarky and trade cases

include("autarkyeq.jl")
include("tradeeq.jl")

using Gadfly

# PARAMETERS
mu_m = .26         # coefficient of nutrition for meat (here beef in cal/kg)
                   # https://www.google.fr/search?q=proteins+in+beef&oq=proteins+in+beef&aqs=chrome..69i57j0l5.6851j0j7&sourceid=chrome&ie=UTF-8
mu_c = .14         # coefficient of nutrition for cereals (here wheat in cal/kg)
                   # https://www.google.fr/search?q=proteins+in+beef&oq=proteins+in+beef&aqs=chrome..69i57j0l5.6851j0j7&sourceid=chrome&ie=UTF-8#q=proteins+in+wheat
Q = (.8*75.0)*365  # subsistence level of proteins in a year, taking the average
                   # weight at 75kg
                   # http://www.health.harvard.edu/blog/how-much-protein-do-you-need-every-day-201506188096

function subsist(Lmax, t, sigma)

  # To call different results
  Nct = length(Lmax)
    if Nct == 1
      res = autarkyeq.autarky_eq(Lmax, t, sigma)
    elseif Nct > 1
      res = tradeeq.trade_eq(Lmax, t, sigma)
    end

  # Optimal consumptions of cereals and meat, optimal price
  qc    = res["Optimal net cereal consumption"]
  qm    = res["Optimal meat consumption"]
  popt  = res["Optimal price"]
  theta = 0.0 # to be replaced
    if Nct == 1
      theta = autarkyeq.theta(autarkyeq.dist(t))
    elseif Nct > 1
      ct = tradeeq.sortct(Lmax, t, sigma)
      Lc = res["Total land to crops"]
      Lsum = sum(Lmax[i] for i in 1:Nct)
      res = tradeeq.num_theta_w(Lc,ct)
    end
  theta

  # # sigma star, the lower bound for satisfaction of the constraint
  # if Nct == 1
  #   sigstar = (1/log(popt))*(log(mu_m-(Q*popt)/theta*Lmax)-log(Q/theta*Lmax - mu_c))
  # else
  #   sigstar = (1/log(popt))*(log(mu_m-(Q*popt)/theta*Lsum)-log(Q/theta*Lsum - mu_c))
  # end
  # if isnan(sigstar)
  #   return sigstar = NaN
  # else
  #   return sigstar
  # end

  # In order to plot this as a graph, I make the function return 0 and 1
  if (mu_m*qm + mu_c*qc) >= Q
    foodcons= 1
  else
    foodcons= 0
  end

  return foodcons, res["Optimal price"] #, sigstar

end # function subsist


################################################################################
#################################### PLOTS #####################################
################################################################################

function plotallautarky(Lmax::Float64, T::Float64)
  set_default_plot_size(20cm, 10cm)
  nevals = 100
  gridsigma = linspace(1e-2, 1-1e-2, nevals) # setting a grid for land endowment
    ysigma = [subsist(Lmax, T, gridsigma[j])[1] for j in 1:nevals] #corresponding y using function autarky_eq
    sigmaplot = plot(x = gridsigma, y = ysigma, Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("Food sec const"), Guide.title("Fixing Lmax=$Lmax and T = $T"))
  println("Optimal price = ", subsist(Lmax, T, .5)[2])
  # println("Sigma star = ", subsist(Lmax, T, .5)[3])
  display(sigmaplot)
end

function plotalltrade(Lmax::Vector, T::Vector)
  set_default_plot_size(20cm, 10cm)
  nevals = 30
  gridsigma = linspace(1e-1, 1-1e-1, nevals) # setting a grid for land endowment
    ysigma = [subsist(Lmax, T, [gridsigma[j] for i in 1:length(Lmax)])[1] for j in 1:nevals] #corresponding y using function autarky_eq
    sigmaplot = plot(x = gridsigma, y = ysigma, Theme(default_color=colorant"darkblue", default_point_size=1.5pt, highlight_width=0.05pt), Geom.point, Geom.line, Guide.xlabel("CES"), Guide.ylabel("Food sec const"), Guide.title("Fixing Lmax=$Lmax and T = $T"))
  println("Optimal price = ", subsist(Lmax, T, [.5 for i in 1:length(Lmax)])[2])
  # println("Sigma star = ", subsist(Lmax, T, [.5 for i in 1:length(Lmax)])[3])
  display(sigmaplot)
end

end #module
