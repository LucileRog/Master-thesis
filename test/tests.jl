module Modeltests

include("C:/Users/lucil/OneDrive/Documents/GitHub/Master_thesis/src/autarkyeq.jl")

#using autarkyeq, tradeeq, dynamictemp
using Base.Test


@testset "Autarky equilibrium tests" begin
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


@testset "Trade equilibrium tests" begin

end

@testset "Dynamic temperature tests" begin

end


end # Module
