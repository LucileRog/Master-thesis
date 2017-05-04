module Modeltests

using Base.Test

################################################################################
############################## TESTS AUTARKYEQ #################################
################################################################################

@testset "Autarky equilibrium tests" begin

include("../src/autarkyeq.jl")
tol = 1e-8
ntest = 5
grid = hcat(linspace(1.0, 1000.0, ntest), linspace(5.0, 25.0, ntest), linspace(.3, 1-1e-4, ntest))
res = collect(autarkyeq.autarky_eq(grid[i,1], grid[i,2], grid[i,3]) for i in 1:ntest)

  @testset "Basic values" begin
  zeroland = autarkyeq.autarky_eq(0.0, 15.0, .9)
    @test zeroland["Optimal net cereal consumption"] == 0.0
    @test zeroland["Optimal meat consumption"] == 0.0
  extremetemp = autarkyeq.autarky_eq(0.0, 30.0, .9)
    @test zeroland["Optimal net cereal consumption"] == 0.0
    @test zeroland["Optimal meat consumption"] == 0.0
  end

  @testset "RD, Analytical vs. Numerical" begin
    for i in 1:ntest
      if isnan(res[i]["Optimal Relative Demand"]) == true
        @test isnan(res[i]["Optimal Relative Demand num"]) == true
      else
        @test abs(res[i]["Optimal Relative Demand"] - res[i]["Optimal Relative Demand num"]) < tol
      end
    end
  end

  @testset "L_c(q_c), Analytical vs. Numerical" begin
    for i in 1:ntest
      if isnan(res[i]["Percentage of land to crop"]) == true
        @test isnan(res[i]["Percentage of land to crop num"]) == true
      else
        @test abs(res[i]["Percentage of land to crop"] - res[i]["Percentage of land to crop"]) < tol
      end
    end
  end

  @testset "Optimal q_c and q_m, Analytical vs. Numerical" begin
  tol = 1e-1
    for i in 1:ntest
      @test abs(res[i]["Optimal net cereal consumption"] - res[i]["Optimal net c C° num"]) < tol
    end
    for i in 1:ntest
      @test abs(res[i]["Optimal meat consumption"] - res[i]["Optimal m C° num"]) < tol
    end
  end

end # autarky tests set


################################################################################
######################### TESTS TRADEEQ & TRADEEQDF ############################
################################################################################


@testset "Trade equilibrium tests" begin

include("../src/autarkyeq.jl")
include("../src/tradeeq.jl")
include("../src/tradeeqdf.jl")
tol = 1e-8
ntest = 3
worlds = (hcat([1200.0, 650.0], [15.4, 14.5], [.1,.1]), hcat([1200.0, 800.0], [17.4, 14.5], [.9,.9]), hcat([300.0, 1600.0], [13.9, 14.3], [.3,.3]))
res = collect(tradeeq.trade_eq(worlds[i][:,1], worlds[i][:,2], worlds[i][:,3]) for i in 1:ntest)

  @testset "Basic values" begin
  zeroland = tradeeq.trade_eq([0.0,0.0], [13.6,15.4], [.6,.6])
    @test zeroland["Optimal net cereal consumption"] == 0.0
    @test zeroland["Optimal meat consumption"] == 0.0
  extremetemp = tradeeq.trade_eq([1200.0,1960.0], [4.0,30.0], [.6,.6])
    @test zeroland["Optimal net cereal consumption"] == 0.0
    @test zeroland["Optimal meat consumption"] == 0.0
  end

  @testset "Maximum aggregate theta" begin
  tol = 1e-5
    ct1 = tradeeq.sortct([400.0, 400.0], [15.0, 15.0], [.9, .9])
    ct2 = tradeeq.sortct([100.0, 400.0], [15.0, 15.0], [.9, .9])
    ct3 = tradeeq.sortct([400.0, 400.0], [15.0, 14.0], [.9, .9])
    gridLc = [100.0, 400.0, 800.0]
    # ct1 : 2 countries identical
      for i in 1:3
        @test tradeeq.num_theta_w(gridLc[i], ct1)[1] == [gridLc[i]/2, gridLc[i]/2]
      end
    # ct2 : different thetas
      @test tradeeq.num_theta_w(gridLc[1], ct2)[1] == [gridLc[1]/2, gridLc[1]/2]
      for i in 2:3
        @test tradeeq.num_theta_w(gridLc[i], ct2)[1] == [100, gridLc[i]-100]
      end
    # ct3 : different land endowments
      for i in 1:2
        @test [tradeeq.num_theta_w(gridLc[i], ct3)[1] .- [gridLc[i], 0] .< tol][1] == trues(2)
      end
      @test [tradeeq.num_theta_w(gridLc[2], ct3)[1] .- [gridLc[3]/2, gridLc[3]/2] .< tol][1] == trues(2)
  end

  @testset "Imports = Exports" begin
      resdf = collect(tradeeqdf.datafr(worlds[i][:,1], worlds[i][:,2], worlds[i][:,3]) for i in 1:ntest)
        for i in 1:ntest
          @test abs(sum(resdf[i][7][j] for j in 1:length(resdf[i][7]))) < 10.0
          @test abs(sum(resdf[i][8][j] for j in 1:length(resdf[i][8]))) < 1e-7
        end
  end

end # Trade equilibrium tests

################################################################################
######################### TESTS TRADEEQ & TRADEEQDF ############################
################################################################################

@testset "Subsistence consumption condition" begin

include("../src/subsistconst.jl")

  @testset "Basic values" begin
    @test subsistconst.subsist(0.0, 15.0, .9)[1] == 0.0
    @test subsistconst.subsist(1000.0, 5.0, .9)[1] == 0.0
    @test subsistconst.subsist([0.0,0.0], [15.0,15.0], [.9,.9])[1] == 0.0
    @test subsistconst.subsist([1000.0,1000.0], [5.0,5.0], [.9,.9])[1] == 0.0
  end

end # subsistence const


end # Module
