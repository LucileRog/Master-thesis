module dynamictemp_test

# module to find the dynamic path of temperature, emissions and meat consumption
# when temeprature is endogenised

include("autarkyeq.jl")
#include("src/tradeeq.jl")
include("subsistconst.jl")

using Gadfly
using Interpolations
using Optim
using Roots


# PARAMETERS

# Technology parameters (do not vary from one country to another)
a      = 1e4      # how many kg of meat produced with one unit of land (ha)
b      = 1/10     # how many kg of meat produced with one unit of cereals (kgs)
gamma1 = 3800/15  # how many kgs of cereals with one unit of land (ha)
                  # http://data.worldbank.org/indicator/AG.YLD.CREL.KG
gamma2 = 600      # slope coef of how crop productivity decreases with distance to
                  # optimal temperature
# Optimal temeprature
Topt   = 15.0

# Temperature dynamics parameters
k1     = 200.0 # kg of GHG emissions per kg of meat produced
k2     = .9   # Natural regeneration of GHG
k3     = 1e-2 # coefficient of how temperature rises with GHG emissions
beta   = 0.95 # intertemporal discount rate
GHG0   = 100

function dynamictemp(Lmax, T0, sigma)

  Lmax = 1000.0
  T0 = 16.0
  sigma = .9

  res    = autarkyeq.autarky_eq(Lmax, T0, sigma)

  ######################################
  ############### DEMAND ###############
  ######################################

    maxqm  = res["Max qm"] # upper bound for qm
    n      = 50  # number of grid points (for n periods)
    N_iter = 3000  # number of iterations

    Tgrid  = linspace(1e-4, maxqm-1e-4, n+1)  # equispaced grid

    # Define theta as a function of T

        function theta(T)
          max(gamma1*Topt - gamma2*abs(Topt-T) , 0.0)
        end

        function theta_prim(T)
          if  Topt * (1 - gamma1/gamma2) < T < Topt
             gamma1
          elseif Topt * (1 + gamma1/gamma2) > T > Topt
            - gamma1
          else
            0.0
          end
        end

        # Define the utility function

    # Analytical solution

    using Roots

    ufun(qc::Float64, qm::Float64)    = (qc^((sigma-1)/sigma) + qm^((sigma-1)/sigma))^(sigma/(sigma-1))
    uprim_c(qc::Float64, qm::Float64) = qc^(-1/sigma) * ufun(qc, qm)^(1/sigma)
    uprim_m(qc::Float64, qm::Float64) = qm^(-1/sigma) * ufun(qc, qm)^(1/sigma)
    # qc as a function of qm from the supply side
    f(qm::Float64, T::Float64)          = theta(T)*Lmax - qm*(theta(T)/a + 1/b)
    fprim(qm::Float64, T::Float64)      = - (theta(T)/a + 1/b)

    n=6
    qm   = [res["Optimal meat consumption"] for i in 1:n]
    GHG  = [GHG0 for i in 1:n]
    Temp = [T0 for t in 1:n]
      for t in 2:n
        Temp[t] = T0 + k3*(k2* GHG[t-1] + k1*qm[t-1])
        obj(x, t) = uprim_c(f(qm[t-1],Temp[t-1]),qm[t-1])/uprim_c(f(x, Temp[t]),x)*(theta(Temp[t-1])/a + 1/b) + uprim_m(f(qm[t-1],Temp[t-1]),qm[t-1])/uprim_c(f(x, Temp[t]),x)+ k1*k3*beta*(Lmax-x/a)*theta_prim(Temp[t])
        if length(fzeros(x -> obj(x,t), 1e-13, maxqm-1e-13)) == 1
          qm[t] = fzeros(x -> obj(x,t), 1e-13, maxqm-1e-13)[1]
        else
          qm[t] = 0.0
        end
      end

    # Numerical solution

# Utility function
n = 1
qmgrid = linspace(1-6, maxqm, 100)
ufun(qc::Float64, qm::Float64)    = (qc^((sigma-1)/sigma) + qm^((sigma-1)/sigma))^(sigma/(sigma-1))
uprim_c(qc::Float64, qm::Float64) = qc^(-1/sigma) * ufun(qc, qm)^(1/sigma)
uprim_m(qc::Float64, qm::Float64) = qm^(-1/sigma) * ufun(qc, qm)^(1/sigma)
# qc as a function of qm from the supply side
f(qm::Float64)           = theta(T)*Lmax - qm*(theta(T)/a + 1/b)
fprim(qm::Float64)       = - (theta(T)/a + 1/b)
# Temperature at t+1 as a function of previous meat consumption
Temp = [T0 for t in 1:n]
GHG  = zeros(n)
qm   = zeros(n)
for t in 2:n
  GHG[t] = k2^t * GHG0 + sum(k2^(j-2)*k1*qm[j-1] for j in 2:t)
end
for t in 2:n-1
  Temp[t+1] = T0 + k3 * (k2* GHG[t-1] + k1*qm[t])
end


    function policy_iter(grid,qm0)
        qm1  = zeros(length(grid))     # next guess
        pol_fun = interpolate((collect(grid),), qm0, Gridded(Linear()) )

        # loop over current states
        # of current capital
        for t in 1:length(grid)
            objective(x) = uprim_c(f(x), x)*fprim(x) + uprim_m(f(x), x) + k1*k3
            qm1[t] = fzero(objective, 1e-10, maxqm-1e-10)
        end
        return qm1
    end

    function PFI()
        qm_init = qmgrid
        for iter in 1:N_iter
            qm_next = policy_iter(qmgrid,qm_init)
            # check convergence
            # if maxabs(qm_init.-qm_next) < tol
            #     perrors = maxabs(qm_next.-c_star(kgrid))
            #     println("PFI:")
            #     println("Found solution after $iter iterations")
            #     println("max policy function error = $perrors")
            #     return c_next
            # elseif iter==N_iter
            #     warn("No solution found after $iter iterations")
            #     return c_next
            # end
            qm_init = qm_next  # update guess
        end
        return qm_init
    end

    v = PFI()
    plot(x = kgrid, y =v, Geom.point, Geom.line )




alpha     = 0.65
beta      = 0.95
grid_max  = 2  # upper bound of capital grid
n         = 150  # number of grid points
N_iter    = 3000  # number of iterations
kgrid     = 1e-6:(grid_max-1e-6)/(n-1):grid_max  # equispaced grid
f(x) = x^alpha  # defines the production function f(k)
tol = 1e-9

ab        = alpha * beta
c1        = (log(1 - ab) + log(ab) * ab / (1 - ab)) / (1 - beta)
c2        = alpha / (1 - ab)
# optimal analytical values
v_star(k) = c1 .+ c2 .* log(k)
k_star(k) = ab * k.^alpha
c_star(k) = (1-ab) * k.^alpha
ufun(x) = log(x)


function policy_iter(grid,c0,u_prime,f_prime)

    c1  = zeros(length(grid))     # next guess
    pol_fun = interpolate((collect(grid),), c0, Gridded(Linear()) )

    # loop over current states
    # of current capital
    for (i,k) in enumerate(grid)
        objective(c) = u_prime(c) - beta * u_prime(pol_fun[f(k)-c]) * f_prime(f(k)-c)
        c1[i] = fzero(objective, 1e-10, f(k)-1e-10)
    end
    return c1
end

uprime(x) = 1.0 ./ x
fprime(x) = alpha * x.^(alpha-1)

function PFI()
    c_init = kgrid
    for iter in 1:N_iter
        c_next = policy_iter(kgrid,c_init,uprime,fprime)
        # check convergence
        if maxabs(c_init.-c_next) < tol
            perrors = maxabs(c_next.-c_star(kgrid))
            println("PFI:")
            println("Found solution after $iter iterations")
            println("max policy function error = $perrors")
            return c_next
        elseif iter==N_iter
            warn("No solution found after $iter iterations")
            return c_next
        end
        c_init = c_next  # update guess
    end
end



      #   function bellman_operator(grid,v0)
      #
      #       v1   = vcat(0.0, zeros(n))         # next guess, t=0, t=1, ... t=n
      #       Temp = collect(T0 for i in 1:n+1)    # Temperature at each period
      #       qm  = vcat(0.0, zeros(n))          # consumption policy function (meat consumption)
      #       GHG  = vcat(GHG0, zeros(n))        # total GHG emissions
      #
      #       Interp = interpolate((collect(grid),), v0, Gridded(Linear()) )
      #
      #       # loop over current states of temperature
      #       for i in 2:n
      #           GHG[i] = k2^i * GHG0 + sum(k2^(j-2)*k1*qm[j-1] for j in 2:i)
      #           objective(qm) = - ( (((theta(Temp[i])*Lmax)-qm*(theta(Temp[i])/a + 1/b))^((sigma-1)/sigma) + qm^((sigma-1)/sigma))^(sigma/(sigma-1)) + beta * Interp[Topt + k3 *(k2*GHG[i-1] + k1*qm)])
      #           # find max of ojbective between [0,k^alpha]
      #           res = Optim.optimize(objective, 1e-6, maxqm)
      #           Temp[i] = Topt + k3 *(k2*GHG[i-1] + k1*res.minimizer) # T_(t+1)
      #           qm[i] = res.minimizer # qm_t
      #           v1[i] = -res.minimum  # v_t
      #       end
      #       return (v1, qm, Temp)   # return both value and policy function
      #   end
      #
      #   function VFI()
      #       v_init     = vcat(0.0, zeros(n)) # initial guess
      #       qm_init    = vcat(0.0, zeros(n)) # initial guess
      #       temp_init  = vcat(0.0, zeros(n)) # initial guess
      #       for iter in 1:N_iter
      #           v_next    = bellman_operator(Tgrid,v_init)  # returns a tuple: (v1,qm, Temp)
      #           v_init    = v_next[1]  # update guess
      #           qm_init  = v_next[2]
      #           temp_init = v_next[3]
      #       end
      #       return (v_init, qm_init, temp_init)
      #   end
      #
      # v    = VFI()[1]
      # qm  = VFI()[2]
      # temp = VFI()[3]
      # plot(x=Tgrid,y=[v[i] for i in 1:n+1], Geom.point, Geom.line)
      # plot(x=Tgrid,y=[qm[i] for i in 1:n+1], Geom.point, Geom.line)
      # plot(x=Tgrid,y=[temp[i] for i in 1:n+1], Geom.point, Geom.line)










end #module
