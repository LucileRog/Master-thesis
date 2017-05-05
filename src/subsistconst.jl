module subsistconst

# module to find check whether the subsistence consumption constraint is
# satisfied or not for autarky and trade cases

include("autarkyeq.jl")
include("tradeeq.jl")

using Gadfly
using ApproxFun

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
  p     = res["Optimal price"]
    if mu_m*qm + mu_c*qc >= Q
      foodcons = 1
    else
      foodcons = 0
    end

  # The threshold sigma star, the interval of p
  if Nct == 1
    thet  = autarkyeq.theta(autarkyeq.dist(t))
      if (Q*p)/(thet*Lmax)-mu_m > 0 && mu_c-Q/(thet*Lmax) > 0 # if sigstar defined
          sigstar = (1/log(p))*(log((Q*p)/(thet*Lmax) - mu_m)-log(mu_c - Q/(thet*Lmax)))
          if p > (thet*Lmax)*(mu_m+mu_c)/Q -1 && p > exp(log((Q*p)/(thet*Lmax) - mu_m)-log(mu_c - Q/(thet*Lmax)))
            pinter = true
          else
            pinter = false
          end
        else
          sigstar = NaN
          pinter = NaN
      end
  else
    sigstar = nothing
    pinter  = nothing
  end

  return foodcons, p, sigstar, pinter

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
  println("Price over both lowerbounds: ", subsist(Lmax, T, .5)[4])
  println("Sigma star = ", round(subsist(Lmax, T, .5)[3], 3))
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
