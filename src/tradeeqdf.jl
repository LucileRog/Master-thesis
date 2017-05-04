module tradeeqdf

# Module to print dataframes summarizing the results for the trade case
# both at the country-level and world level

include("autarkyeq.jl")
include("tradeeq.jl")

using DataFrames

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


# distance to optimal temperature for crops
  function d(t)
    abs(Topt - t)
  end

# Crop productivity at the country level
  function theta(d::Vector)
    max([gamma1 * Topt - gamma2 *d[i] for i in 1:length(d)], 0.0)
  end


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

function datafr(Lmax::Vector, t::Vector, sigma::Vector)
  ct = sortct(Lmax, t, sigma)      # Matrix of countries, sorted (L, theta, sig)
  Nct = length(ct[:,1])
  countries = hcat(Lmax, t, sigma) # Matrix of countries, not sorted (L, t, sig)
  worldres = tradeeq.trade_eq(Lmax, t, sigma)
  qc = worldres["Optimal net cereal consumption"]
  worldLc = worldres["Total land to crops"]

  # Share of lands to crops per country
    shareLc = [round(tradeeq.num_shares(qc, ct)[1][i], 3) for i in 1:Nct]

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
    sharemprod = [round(tradeeq.num_shares(qc, ct)[2][i], 3) for i in 1:Nct] # because technologies are the same

  # Imports/exports
    # cereals
      worldccons = worldres["Optimal net cereal consumption"]
      tradec = zeros(Nct)
      for i in 1:Nct
        tradec[i] = (sharecprod[i] - 1/Nct)*worldccons
      end
      tradec

    # meat
      worldmcons = worldmprod
      tradem = zeros(Nct)
      for i in 1:Nct
        tradem[i] = (sharemprod[i] - 1/Nct)*worldmcons
      end
      tradem

    # Relative price and relative demand
      popt_aut = zeros(Nct) # to be replaced
      RD_aut   = zeros(Nct) # to be replaced
      for i in 1:Nct
        popt_aut[i] = autarkyeq.autarky_eq(ct[i,1], ct[i,4], ct[i,3])["Optimal price"]
        RD_aut[i]   = autarkyeq.autarky_eq(ct[i,1], ct[i,4], ct[i,3])["Optimal Relative Demand"]
      end

return DataFrame(Land_endowment      = ct[:,1],
                 Temperature         = ct[:,4],
                 CES                 = ct[:,3],
                 Crop_yields         = ct[:,2],
                 Share_cereal_prod   = sharecprod,
                 Share_meat_prod     = sharemprod,
                 Trade_in_cereals    = tradec,
                 Trade_in_meat       = tradem,
                 Rel_price_autarky   = popt_aut,
                 Rel_demand_autarky  = RD_aut)

end


function eqdatafr(Lmax::Vector, t::Vector, sigma::Vector)
  res = tradeeq.trade_eq(Lmax, t, sigma)
  return DataFrame(Net_cereal_cons   = res["Optimal net cereal consumption"],
                   Meat_cons         = res["Optimal meat consumption"],
                   Per_land_to_crops = res["Percentage of land to crop"],
                   Opt_price         = res["Optimal price"],
                   Rel_demand        = res["Optimal Relative Demand"])
end


function runall(Lmax::Vector, T::Vector, sigma::Vector)
  # data frames
    df1 = datafr(Lmax, T, sigma)
    df2 = eqdatafr(Lmax, T, sigma)
  # print
    println("World where Lmax = $Lmax, Temperature = $T, and CES=$sigma")
    println(df1)
    println(df2)
end






end #module
