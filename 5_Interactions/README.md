# Particle interactions in plasma

This chapter is mostly about numerical integration and interpolation applied to the study of distribution functions.

## Numerical integration
We will be using the trapezoidal rule here.
This is not the best choice but (hopefully) it will be enough.
If you are interested in better algorithms read about the [Newton-Cotes rules](https://en.wikipedia.org/wiki/Newton%E2%80%93Cotes_formulas) or the [Gaussian quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature)

To numerically integrate a function, you could just sum the functional values evaluated in the center of some intervals multiplied by the width of the interval.
This corresponds to dividing the function into rectangles.
This is called the midpoint rule. This is not optimal, because you would need a large number of intervals and, therefore, function evaluations to get a good accuracy.

The trapeziodal rule is a simple improvement.
Instead of dividing into rectangles use trapezoids.

I have not found a Julia package which implements the trapezoidal rule for numeric integration, so I have created my own function:
```julia
function trapz(y::AbstractVector; x=1:length(y))
    #Trapezoidal rule
    dx = diff(x)
    sum([ (y[i-1] + y[i])/2*dx[i-1] for i in 2:length(y)])
end
```

Here is a simple example of usage:
```julia
f(x) = 1/(1 + x^2)

for i in [100,1000,10000,100000]
    x = range(0,stop=1,length=i)
    y = f.(x)
    integral = trapz(y,x=x)
    println(4*integral)
end
println(pi)
```

## Maxwell-Boltzmann distribution
Now, you should apply numerical integration to Maxwell-Boltzmann distribution to calculate mean values and compare these with the theoretical values.

## Exercises
>  **Exercise 1**
>  * Complete this [script](https://github.com/tungli/F5170-julia/blob/master/5_Interactions/maxwell.jl). The script should plot the Maxwell-Boltzmann distribution with repsect the velocity magnitude and a few points of interest, namely:
>    1. The mean speed
>    2. The mean squared speed (norm of the velocity vector)
>    3. The most probable speed
>  * Calculate these values analytically and numerically. Compare them and explain any differences.
>  * Look at the plot. Which of the 3 speeds is the lowest. Does their order change with temperature?
>  * Change the number of points of integration. How does it affect the result of the integration?
>  
>  **Exercise 2**
>  * The mass of nitrogen molecules is 28 a.m.u. and their number density at ambient conditions is approximately *1.7e25* m<sup>−3</sup>. How many nitrogen molecules in your room are faster than 50, 500, 1000, 2500, 5000 and 10000 m/s?
>  
>  **Exercise 3**
>  The figure below shows three Maxwell-Boltzmann distributions.
>  * Assume that the distributions differ only in temperature. Which of the distributions has the highest temperature and which the lowest?
>  * Assume that the distributions differ only in particle mass. Which of the distributions has the highest mass and which the lowest?
>  
>  ![dists](https://github.com/tungli/F5170-python/blob/master/5_Interactions/dists.svg)
>  

```julia
using Plots
pyplot()

v = range(0,stop=4000,length=500)

function Fv(v)
    amu = 1.66053904e-27
    m = 40.0*amu
    kB = 1.3806485e-23
    T = 1000.0
    a = m/kB/T
    sqrt(2/pi)*v^2*a^(3.0/2.0)*exp(-a/2*v^2)
end

v_mean =
v_sq_mean =
v_mp =
v_mean_an =
v_sq_mean_an =
v_mp_an =

plot(v,Fv.(v),label="Maxwell-Boltzmann distribution")
scatter!([v_mean],[Fv(v_mean)])
scatter!([v_sq_mean],[Fv(v_sq_mean)])
```

## Time evolution of a distribution function
In the previous section you analyzed the distribution function of velocity magnitude.
In this section you will be using 3 distribution functions for 3 Cartesian components of the velocity vector and you will look at a simple time evolution of a system described by distrubtion functions.

Our simple model is based on the [Boltzmann transport equation](https://en.wikipedia.org/wiki/Boltzmann_equation)

![boltzmannEq](http://mathurl.com/y9qsxtzt.png)

We will assume homogenity in spatial coordinates, zero external force and for the collision term we will assume the following form ([Krook](https://en.wikipedia.org/wiki/Bhatnagar%E2%80%93Gross%E2%80%93Krook_operator)):

![krook](http://mathurl.com/y77eakb6.png)

where *ν<sub>m</sub>* is the collision frequency.

Our kinetic equation therefore simplifies greatly, in fact, it can be integrated analytically:

![kinetic](http://mathurl.com/ya2qhqxq.png)


## Exercises
>  **Exercise 4**
>  If all goes well, this [script](https://github.com/tungli/F5170-julia/blob/master/5_Interactions/evol.jl) should make a *gif* file with an animation of the distribution time evolution.
>  Answer the following questions:
>  * What is the physical meaning of the initial condition for the distribution function?
>  * What kind of particles could the distribution functions describe?
>  * The collision frequency is *5e8* Hz, which is a resonable value. What is the time necessary for reaching the equilibrium?
>  * Try increasing and decreasing the collision frequency in the script. What happens with the time necessary for reaching equilibrium and why?
  

```julia
using Plots
pyplot()

v = range(-2000,stop=2000,length=200)
t = range(0,stop=1e-8,length=100)
amu = 1.66053904e-27
m = 40.0*amu
kB = 1.3806485e-23
T = 1000.0
coll_freq = 5e8
```
This is just the usual importing and defining constants.

```julia
gaussian(x,x0,sigma) = exp(-(x-x0)^2/(2*sigma^2))/sqrt(2*pi)/sigma

fx0 = gaussian.(v,500.0,10.0)
fy0 = zeros(length(v))
fz0 = zeros(length(v))

v_sq = sqrt(trapz((fx0+fy0+fz0).*v.^2,x=v))

Teq = m*v_sq^2/kB/3

fx_eq = sqrt.(m/(2*pi*kB*Teq))*exp.(-m*v.^2/(2*kB*Teq))
fy_eq = sqrt.(m/(2*pi*kB*Teq))*exp.(-m*v.^2/(2*kB*Teq))
fz_eq = sqrt.(m/(2*pi*kB*Teq))*exp.(-m*v.^2/(2*kB*Teq))

calc_distribution(feq,f0,t,coll_f) = feq + (f0 - feq)*exp(-t*coll_f) 
```
Here we define the initial and equilibrium distribution functions and a function which calculates the distribution function for any time.

```julia
@gif for i in t
    px = plot(v,calc_distribution(fx_eq,fx0,i,coll_freq),ylim=([0,maximum(fx0)]))
    py = plot(v,calc_distribution(fy_eq,fy0,i,coll_freq),ylim=([0,maximum(fx0)]))
    pz = plot(v,calc_distribution(fz_eq,fz0,i,coll_freq),ylim=([0,maximum(fx0)]))
    plot(px,py,pz,layout=(1,3),legend=false)
end
```
These just the things for the animation.
You can add some labels or change the colors, etc.

Here is an example of one frame from the animation:

![Evolution](https://github.com/tungli/F5170-python/blob/master/5_Interactions/evol_ex.svg)

## Rate coefficients
If you have a system described by a distribution function you have access to macroscopic properties of the system.
Here, we will demonstrate this by calculating the rate constants of reactions taking place in plasma.
Rate constants are calculated from collisional cross sections *σ* of a reaction indexed *r* as:

![rateintegral](http://mathurl.com/yde22kcc.png)

Cross sections of most reactions do not have functional representations - they are available only as tabulated data.

You will find the data you need in this directory, here are the [cross sections](https://github.com/tungli/F5170-python/blob/master/5_Interactions/sigmaion.dat) to interpolate.
The first column of the data are speeds in m/s, the second column are the cross sections in m<sup>2</sup>.
The cross section are those of argon ionization by electron impact:

![arelimpioni](http://mathurl.com/ybcgs3t8.png)


>  ## Exercises
>  **Exercise 5**
>  * Run this [scipt](https://github.com/tungli/F5170-julia/blob/master/5_Interactions/interpol_rates.jl). It calculates the rate constant from the cross section data using two different distribution function - a Maxwell-Boltzmann and a nearly monoenergetic distribution (delta function approximated by a Gaussian function). You will notice that the mean velocity is the same for both distributions but the rate constants differ. Provide an explanation.
>  
>  **Exercise 6**
>  * Modify the script so that it plots cross section as function of speed.
>  * Answer the following questions:
>    1. What is the ionization threshold in electronvolts?
>    2. Where does the cross section reach the maximum value?
>  
>  **Exercise 7**
>  * Run the script for several electron temperatures.
>  * What happens with the rate coefficients with increasing electron temperature?
>  * Is there electron temperature for which the rate coefficient for the nearly monoenergetic beam exceeds the Maxwell-Boltzmann coefficient? Provide an explanation why this is/is not possible
>  
>  **Advanced Exercise**
>  * Rewrite the scipt so that it uses electron energy rather than electron speed. Using electron energy is more common in plasma physics.
  

**(Potentially) Important Note:**
I have decided on using the [Dierckx.jl](https://github.com/kbarbary/Dierckx.jl) package in favour of [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl).
At the time I was writing this, the package available _**by using `Pkg.add()`**_ was still for older version of Julia.
Nevertheless, on github the version was corrected, therefore you can just clone [the repository](https://github.com/kbarbary/Dierckx.jl) by using:
```julia
#Pkg.rm("Dierckx")
Pkg.add("https://github.com/kbarbary/Dierckx.jl")
#Pkg.resolve()
```
Hopefully this will work. Contact me if it does not or use a different interpolation method.

```julia
using Dierckx
using DelimitedFiles

data = readdlm("sigmaion.dat",',')
x = data[:,1]
y = data[:,2]

function Fv(v)
    m = 9.109e-31
    kB = 1.3806485e-23
    T = 11604.0*15.0
    a = m/kB/T
    sqrt(2/pi)*v^2*a^(3.0/2.0)*exp(-a/2*v^2)
end

gaussian(x,x0,sigma) = exp(-(x-x0)^2/(2*sigma^2))/sqrt(2*pi)/sigma

function trapz(y::AbstractVector; x=1:length(y))
    #Trapezoidal rule
    dx = diff(x)
    sum([ (y[i-1] + y[i])/2*dx[i-1] for i in 2:length(y)])
end

v = range(0.0,stop=1e7,length=10^6)
cr_sec = Spline1D(x,y,k=3,bc="zero")

@show kr_MB = trapz(v.*Fv.(v).*cr_sec.(v),x=v)
@show vmean = trapz(Fv.(v).*v,x=v)
@show kr_beam = trapz(v.*gaussian.(v,vmean,vmean/100).*cr_sec.(v),x=v)
```

The scipt involves importing some packages, getting the data from a file, interpolating and integrating some functions.

The interpolation takes place in this line:
```julia
cr_sec = Spline1D(x,y,k=3,bc="zero")
```
Notice the arguments - `k=3` means to use the 3<sup>rd</sup> order (cubic) interpolation, `bc="zero"` is to allow extrapolation with all the extrapolated values being zero.
The `cr_sec` object will act like a function - meaning you can write `cr_sec(0.5e7)` to get the value at `0.5e7`.

Here is an example plot with the interpolation:

![Rates Interpolation](https://github.com/tungli/F5170-python/blob/master/5_Interactions/rates_interp.svg)


