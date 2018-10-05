# Particle balance in plasma

![Concentration Evolution](https://github.com/tungli/F5170-python/blob/master/6_Balance/reactions.svg)

This chapter is about reactions in plasma.
Plasmas usually contain a large number of species.
These can be electrons, different elements or molecules, as well as ionized and excited species.
It can be very useful be able to predict the concentration of species in time and space.

Here, we will restrict our description to time evolution and assume two different temperatures - one for electrons, one for the heavy particles (gas temperature).
The reaction scheme we will use is one for argon plasma at atmospheric pressure.
In the picture below you will find the reactions and reaction rate temperature dependencies (third column), the data are from [this article](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.85.056404).

![scheme](https://github.com/tungli/F5170-python/blob/master/6_Balance/rate_table.png)

To transform a set of chemical reactions to a set of ordinary differential equations (ODEs) describing the evolution of every species here is a general approach:
 * Take the *r*th equation of *M* reactions in total. It has the form of ![reac](http://mathurl.com/ycnjqt5p.png) where *a* and *b* are the stoichiometric coefficients (always nonnegative) of species *X*. Generally, there are two rate constants - forward and backward rates.
 * The ODE for the *i*th species is:

![masterODE](http://mathurl.com/yd4yxekn.png)

we used `[]` to represent the concentration.

This may look scary, so here is an example with one-way reactions (like the ones for argon):

![ex](http://mathurl.com/ycv8c7w5.png)

**Note**:
A reaction with 3 reactants does not necessarily represent a three-body collision.
If you have two reactions:

![twobody](http://mathurl.com/yaecpw9l.png)

and the second reaction is faster (more probable), then we may write this as:

[!threebody](http://mathurl.com/yadrlr8h.png)


## Exercises
>  **Exercise 1**
>  * What is the unit of the time derivative of the concentration (source term)? What must be the unit of the rate constant for a reaction involving two reactants and what for a reaction involving three reactants?
>  
>  **Exercise 2**
>  If the argon gas was without excited and ionized species, the concentration of argon could be calculated from the equation of state for ideal gas.
>  * Calculate the actual concentration of ground state argon in plasma given the gas temperature, gas pressure, and the concentration of excited and ionized species (Ar<sup>+</sup>, Ar<sup>\*</sup>, Ar<sup>+</sup><sub>2</sub>)
>  
>  **Exercise 3**
>  Look at the [reaction table picture](https://github.com/tungli/F5170-python/blob/master/6_Balance/rate_table.png), reactions 1, 3, 4.
>  * What is the unit of the constant in the exponential? What is its physical meaning?
>  * Which reaction do you think will be the most important ionization channel?
>  
>  **Exercise 4**
>  * Using the [reaction table picture](https://github.com/tungli/F5170-python/blob/master/6_Balance/rate_table.png) express the time derivatives of concentration (source terms) of species: e, Ar<sup>+</sup>, Ar<sup>\*</sup>, Ar<sup>+</sup><sub>2</sub>.
>  
>  **Exercise 5**
>  * Complete the [script](https://github.com/tungli/F5170-julia/blob/master/6_Balance/odesolve.jl). In particular, you need to:
>    1. Write the expressions for the source terms (`SArs`,`SArp`,`SAr2p` and `Se`).
>    2. Create the functions for rate constants `f_k[1-11](Te,Tg)`
>    3. Run the script and decide if the output makes sense.
>  
>  **Exercise 6**
>  Try increasing and decreasing the electron temperature and answer the following questions:
>  * How does the equilibrium electron density change and why?
>  * How does the ignition time change?
>  * What is the dominant ion in the ignition phase and what is the dominant ion when the plasma stabilizes? Does the dominance of the two ions change with electron temperature?
>  
>  **Exercise 7**
>  Run the program for your chosen value of electron temperature and for pressures of 10<sup>3</sup>, 10<sup>4</sup>, 10<sup>5</sup> and 10<sup>6</sup> Pa and answer the following questions:
>  * How does the equilibrium electron density change and why?
>  * How does the ignition time change?
>  
>  **Advanced Exercise**
>  * Modify the program in the previous section so that it solves the system of equations for several values of *T*<sub>e</sub>*/p* and plots the steady-state number densities as functions of *T*<sub>e</sub>*/p*. You can do this by adding a for loop. What makes this exercise difficult is the fact that the ignition time changes quite quickly with electron temperature and pressure. Therefore, the time interval has to be updated in each iteration according to the current value of *T*<sub>e</sub> and *p*, otherwise the solution will take very long.


## Implementation
The details of the script are discussed below.
```julia
using DifferentialEquations
using Plots
pyplot()
using LaTeXStrings

p = 1e5
kB = 1.38e-23
Tg = 400.0
Te = 2.0
```
We import the ODE solver and a plotting package with LaTeX strings.
Then we define some constants.

```julia
tsteps = 1000
tout = 10.0.^range(-11,stop=-6,length=tsteps)
tspan = (minimum(tout),maximum(tout))

nArs_init = 1e12
nArp_init = 1e12
nAr2p_init = 1e12
ne_init = nArp_init + nAr2p_init
initial = [nArs_init, nArp_init, nAr2p_init, ne_init]
```
To solve the system of ODEs, we need to define the time span of the integration and the initial conditions.

```julia
function odefun!(dy,y,params,t)
    nArs = y[1]
    nArp = y[2]
    nAr2p = y[3]
    ne = y[4]
    nAr = p/(kB*Tg) - nArs - 0.5*nAr2p - nArp


    k1 = f_k1(Te,Tg)    
    k2 = f_k2(Te,Tg)
    k3 = f_k3(Te,Tg)
    k4 = f_k4(Te,Tg)
    k5 = f_k5(Te,Tg)
    k6 = f_k6(Te,Tg)
    k7 = f_k7(Te,Tg)
    k8 = f_k8(Te,Tg)
    k9 = f_k9(Te,Tg)
    k10 = f_k10(Te,Tg)
    k11 = f_k11(Te,Tg)

    dy[1] = k1*ne*nAr - k2*ne*nArs - k4*ne*nArs + k6*ne*nAr2p - 2*k10*nArs^2 + k11*nArs*nAr
    dy[2] = 1.0
    dy[3] = 1.0
    dy[4] = 1.0
    return nothing
end
```
This is the first derivatives function.
The functions for the rate constants are missing, you will need to implement them, as well as the right-hand sides for the derivatives.

```julia
prob = ODEProblem(odefun!,initial,tspan,Vector{Float64}())
sol = solve(prob)
```
We use the DifferentialEquations.jl package to solve the problem numerically.

```julia
nAr = p/(kB*Tg) .- [ i[1] - 0.5*i[3] - i[2] for i in sol.u ]
plot(sol,yscale=:log10,xscale=:log10,vars=(0,1),label=L"Ar^*")
plot!(sol,yscale=:log10,xscale=:log10,vars=(0,2),label=L"Ar^+")
plot!(sol,yscale=:log10,xscale=:log10,vars=(0,3),label=L"Ar_2^+")
plot!(sol,yscale=:log10,xscale=:log10,vars=(0,4),label=L"e^-")
plot!(range(tspan[1],stop=tspan[2],length=length(sol)),nAr,label="Ar")
```
Finally, we plot our results using log-log plots.



