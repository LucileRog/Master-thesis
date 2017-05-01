module tradeeqdf

include("autarkyeq.jl")
include("tradeeq.jl")

using DataFrames

# PARAMETERS

# Technology parameters (do not vary from one country to another)
a = 10.0 # coef of conversion of crop into meat
b = 1e-4 # coef of conversion of land into meat
gamma1 = 3800/15 # coef of conversion of land into crop
                # http://data.worldbank.org/indicator/AG.YLD.CREL.KG
gamma2 = 600 # slope coef of how crop productivity decreases with distance to
             # optimal temperature
# Optimal temeprature
Topt = 15.0


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

# A function to generate a sorted list of countries
    function sortct(Lmax::Vector, t::Vector, sigma::Vector)
      # Number of countries
        Nct = length(Lmax)
      # Sum of all countries land endowments
        Lsum = sum(Lmax[i] for i in 1:Nct)
      # Distance from optimal temperature for all countries
        dist = collect(d(t[i]) for i in 1:Nct)
      # Sorting countries by decreasing theta(d)
        countries = hcat(Lmax, theta(dist), sigma, t)
        return sortrows(countries, by=x->x[2], rev=true) # dim = Ncountries x 2
    end

function datafr(Lmax::Vector, t::Vector, sigma::Vector)

  ct = sortct(Lmax, t, sigma)      # Matrix of countries, sorted (L, theta, sig)
  Nct = length(ct[:,1])
  countries = hcat(Lmax, t, sigma) # Matrix of countries, not sorted (L, t, sig)
  worldres = tradeeq.trade_eq(Lmax, t, sigma)

  # Share of lands to crops per country
    worldLc = worldres["Total land to crops"]
    shareLc = zeros(Nct) # to be replaced
    # variable checking wether thetas are different for all countries or not
    # because ct[:,2] are ordered, we can check 2 by 2
    alldiff = [ct[i,2] != ct[i+1,2] for i in 1:Nct-1]

    # for 2 countries
    if Nct == 2
      if alldiff == trues(Nct-1) # if all thetas are different
        if worldLc <= ct[1,1]
          shareLc[1] = 100.0
        else
          shareLc[1] = ct[1,1]/worldLc
          shareLc[2] = 1 - ct[1,1]/worldLc
        end
      else # if both countries have same theta
        ct = sortrows(ct, by=x->x[1]) # sort again, by decreasing Lmax
        if worldLc/2 <= ct[1,1]
          shareLc[1] = 1/2
          shareLc[2] = 1/2
        else # if total land required for crop production per country is larger
             # than the smallest land endowment
          shareLc[1] = ct[1,1]/worldLc
          shareLc[2] = 1 - ct[1,1]/worldLc
        end
      end # same theta
    # Nct = 2
    # PROBLEM : I do not manage to loop over countries
    else # for Nct countries
      if alldiff == trues(Nct-1) #all thetas are different
        # if all thetas are different, then the country with higher theta will
        # have the largest share, then second the second largest etc.
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
      else # if some thetas are the same
        shareLc = [NaN for i in 1:Nct]
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
      if worldmprod != 0.0
        sharemprod[i] = (a*(ct[i,1]-shareLc[i]*worldLc))/worldmprod
      else
        sharemprod[i]
      end
    end
    sharemprod

  # Imports/exports
    # cereals
      worldccons = worldres["Optimal net cereal consumption"]
      tradec = zeros(Nct)
      autqc  = zeros(Nct)
      for i in 1:Nct
        autqc[i] = autarkyeq.autarky_eq(ct[i,1], ct[i,4], ct[i,3])["Optimal net cereal consumption"]
        tradec[i] = autqc[i] - (1/Nct)*worldccons
      end
      tradec

    # meat
      worldmcons = worldres["Optimal meat consumption"]
      tradem = zeros(Nct)
      autqm  = zeros(Nct)
      for i in 1:Nct
        autqm[i] = autarkyeq.autarky_eq(ct[i,1], ct[i,4], ct[i,3])["Optimal meat consumption"]
        tradem[i] = autqm[i] - (1/Nct)*worldmcons
      end
      tradem

return DataFrame(Land_endowment      = ct[:,1],
                 Temperature         = ct[:,4],
                 CES                 = ct[:,3],
                 Crop_yields         = ct[:,2],
                 Share_land_to_crops = shareLc,
                 Share_cereal_prod   = sharecprod,
                 Share_meat_prod     = sharemprod,
                 Trade_in_cereals    = tradec,
                 Trade_in_meat       = tradem)

end


function eqdatafr(Lmax::Vector, t::Vector, sigma::Vector)
  res = tradeeq.trade_eq(Lmax, t, sigma)
  return DataFrame(Net_cereal_cons   = res["Optimal net cereal consumption"],
                   Meat_cons         = res["Optimal meat consumption"],
                   Tot_land_to_crops = res["Total land to crops"],
                   Per_land_to_crops = res["Percentage of land to crop"])
end


function runall()
  # menu
    Lmax = linspace(500.0, 3000.0, 6)
    t = linspace(10.0, 20.0, 9)
    sigma = [[.1 for i in 1:6], [.5 for i in 1:6], [.9 for i in 1:6]]
  # data frames
    # 2 identical countries
      df2a = datafr([Lmax[3] for i in 1:2], [t[4] for i in 1:2], sigma[3][1:2])
      df2a_eq = eqdatafr([Lmax[3] for i in 1:2], [t[4] for i in 1:2], sigma[2][1:2])
    # 2 countries, different Lmax, same t
      df2b = datafr([Lmax[3],Lmax[4]], [t[4] for i in 1:2], sigma[3][1:2])
      df2b_eq = eqdatafr([Lmax[3],Lmax[4]], [t[4] for i in 1:2], sigma[2][1:2])
    # 2 countries, same Lmax, different t
      df2c = datafr([Lmax[3] for i in 1:2], [t[5], t[6]], sigma[3][1:2])
      df2c_eq = eqdatafr([Lmax[3] for i in 1:2], [t[5], t[6]], sigma[2][1:2])
    # 3 countries, different in Lmax and t
      df3 = datafr([Lmax[3], Lmax[5], Lmax[6]], [t[4], t[7], t[5]], sigma[3][1:3])
      df3_eq = eqdatafr([Lmax[3], Lmax[5], Lmax[6]], [t[4], t[7], t[5]], sigma[3][1:3])

  # print all
  println("2 indentical countries")
  println(df2a)
  println(df2a_eq)
  println("--------------------------------------------------------------")
  println("--------------------------------------------------------------")
  println("2 countries with different land endowments but same temperature")
  println(df2b)
  println(df2b_eq)
  println("--------------------------------------------------------------")
  println("--------------------------------------------------------------")
  println("2 countries with same land endowment but different temperatures")
  println(df2c)
  println(df2c_eq)
  println("--------------------------------------------------------------")
  println("--------------------------------------------------------------")
  println("3 countries with different land endowments and temperatures")
  println(df3)
  println(df3_eq)
end

end #module
