# Particle drift

![Drift](https://github.com/tungli/F5170-python/blob/master/3_Motion/drift.svg)

## Mathematical formulation

The exercises want us to solve the equation of motion for a charged particle in a electromagnetic field.
The equation of motion is a set of 3 (x,y,z) ordinary differential equations (ODE) of second order.
To solve a system of ODEs, we usually transform the set of higher order equations to first order equations.  We will create more variables by doing this but that's fine.
We express the first order derivatives as functions of variables and other first order derivatives.
We can then create a function in Julia which for a given set of values of variables evaluates the derivatives.
This is the function that a general ODE solver needs.

In our problem we treat positions and velocities as variables and express their derivatives as functions of these variables:

![Equation of Motion](http://mathurl.com/yaan82k2.png)

where *F* is the Lorentz force.

**Note:** There is another way! Here a general approach to ODE is shown but to integrate second order ODEs it is usually not necessary to transform to first order set. 
Another way is, for example, to use the *leapfrog* algorithm. 
If you are not comfortable with using a black-box solver, this is the way for you. Try googling "Boris leapfrog" and implementing it in your code.

## Exercises

>  **Exercise 1**
>  * Run the [script](https://github.com/tungli/F5170-python/blob/master/3_Motion/motion.py).  
>  * Configure the velocity, position and fields as you want.  
>  * Configure the mass and charge to that of an electron.  
>  * What kind of drift do you observe?  
>  * What is the direction of the drift for an electron a what for a positron?  
>  
>  **Exercise 2**
>  * Now configure the parameters for a proton
>  * How many times do you have to increase/decrease the time scale for the plot of the trajectory to be comparable to that of an electron
>  * Compare the amplitudes of the oscillation and the magnitudes of the drift velocities for proton and electron
>  
>  **Exercise 3**
>  * Study a charged particle in the following field and with the following velocity:
>  
> ![prob3](http://mathurl.com/ycp4a5wj.png)
>  
>  * Do you observe any drift? If yes, for what parameters did you use? What is the direction of the drift velocity for an electron and a positron?
>  * Try to match your observations with theoretical predictions.
>  
>  **Advanced exercise**
>  * Study both an electron and a proton in an electric field which varies harmonically with time and in uniform magnetic field. Try different frequencies.
>  * How do they react to the field?
>  * Compare the effect for various frequencies 


## Implementation
Copy the [motion.jl](https://github.com/tungli/F5170-julia/blob/master/3_Motion/motion.jl) file and run Julia REPL.
Download the packages DifferentialEquations and Plots using the Julia package manager (`using Pkg`, then `Pkg.add("Plots")`,...).
Type `include("PATH TO motion.jl")`.

We will be using the brilliant **DifferentialEquations.jl** package.
Take a look at the [tutorial](http://docs.juliadiffeq.org/latest/tutorials/ode_example.html#Defining-Parameterized-Functions-1).
By following the same procedure we solve our problem. Let us look at it in parts:
```julia
using Plots
pyplot()

using DifferentialEquations
using LinearAlgebra
```
These are the packages we will be using here.
It might take a while for Julia to load these (especially the very first time) - both are huge packages.
Once loaded, the code will run fast.

```julia
function derivs!(dy,y,p,t)
    dy[1:3] = y[4:6]
    dy[4:6] = p[1]/p[2]*(p[3] + cross(u[4:6],p[4]))
end
```
This is the function that updates the first order derivatives. If you open the file `motion.jl` you can notice that in Julia we can make use Unicode characters to make the code look very similar to the mathematical formulation (I do not know how to use Unicode in Markdown here).
Also notice the `!` in the name of our function - this is not mandatory but in Julia the functions that change their inputs have a `!` at the end.

```julia
y0 = [0.0,0,0,0,0,0]
p = [1.0,1.0,[0.0,1,0],[0.0,0,1]]
tspan = (0.0,10.0)
prob = ODEProblem(derivs!,y0,tspan,p)
```
Here, we define the problem as DifferentialEquations.jl needs it.
The first three values of the state vector `y0` are the positions (x,y,z), the other three are the velocities.
In the parameter vector `p` we have the charge, the mass of the particle and then two the two field vectors.
`tspan` is a `Tuple` with the initial time and final time.
Finally, `prob` is the object that encapsulates the entire mathematical problem.

```julia
sol = solve(prob)
```
We solve the problem by calling the `solve()` function on the `ODEProblem` object.
Take a look at the `sol` object.
It first seems that it contains only a very limited number of points from the solutions.
This is fine, since DifferentialEquations.jl by default only solves for the necessary points and interpolates in-between these point to obtain the whole solution.

```julia
plot(sol,vars=(1,2))
```
This line is possible because of **Plots.jl** knows how to handle the `sol` object.
Another way is shown below, where we extract the x,y,z vectors and plot those.
```julia
t = range(tspan[1],length=1000,stop=tspan[2])
a = sol(t)'[:,1:3]
x = a[:,1]
y = a[:,2]
z = a[:,3]

plot(x,y)
```
To configure the plot to your liking take a look at some [examples](http://docs.juliaplots.org/latest/) and try them out.
You can, for example, make a GIF to see the particle trajectory moving.

# Van Allen radiation belt

![Van Allen](https://github.com/tungli/F5170-python/blob/master/3_Motion/van_allen.svg)

Earth's [Van Allen radiation belts](https://en.wikipedia.org/wiki/Van_Allen_radiation_belt) is a good example of charged particle moving in a field and suitable for numerical simulation.

## Exercises
>  **Exercise 4**
>  Let us assume that the magnetic moment of the Earth is accurately described by a single magnetic dipole moment and orient our frame of reference in such a way that this the magnetic moment is **_m_** = *(0,0,M)* in Cartesian coordinates.
>  * Express components of the magnetic field **_B_** in Cartesian coordinates.
>  * What is the value of *M* if the geomagnetic field at the equator is *3.12e-5* T?
>
>  **Exercise 5**
>  * Analyze the motion of the high-energy proton in the geomagnetic field. What are the three components of the motion?
>
>  **Exercise 6**
>  * Try changing the initial position of the proton
>  * What is the maximum initial distance from Earth's center for which the proton still has a stable (bound) trajectory?
>  * What is the minimum initial distance for which the proton does not hit the surface of the Earth?
>
>  **Advanced exercise**
>  * Replace the proton with an electron and try to find suitable initial conditions for a stable (bound) trajectory. Think before implementing it. Will the magnetic field required for an electron be higher or lower? How does the drift differ from that of the proton?

Since the same equations and solvers apply here, let us jump to the implementation right away.

## Implementation
Take a look at the [script](https://github.com/tungli/F5170-julia/blob/master/3_Motion/van_allen.jl).
You can see that the implementation is roughly the same as with the particle drift study.
Some differences here are:
 - We assume electric field is zero
 - Earth's magnetic field can be approximated as if produced by a (huge) magnetic dipole.

You will need to implement the `bfield()` function that calculates the magnetic field from the "Earth dipole" at a certain position and calculate the magnetic dipole moment.
Good luck!


