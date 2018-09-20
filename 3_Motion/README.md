# Particle drift
This sections covers exercises 3.1, 3.2, 3.3 and 3.1.a.

## Matematical formulation

The exercises want us to solve the equation of motion for a charged particle in a electromagnetic field.
The equation of motion is a set of 3 (x,y,z) ordinary differential equations (ODE) of second order.
To solve a system of ODEs, we usually transform the set of higher order equations to first order equations.
We will create more variables by doing this but that's fine.
We express the first order derivatives as functions of variables and other first order derivatives.
We can then create a function in Julia which for a given set of values of variables evaluates the derivatives.
This is the function that a general ODE solver needs.

In our problem we treat positions and velocities as variables and express their derivatives as functions of these variables:

![Equation of Motion](http://mathurl.com/yaan82k2.png)

where *F* is the Lorentz force.

**Note:** There is another way! Here a general approach to ODE is shown but to integrate second order ODEs it is usually not necessary to transform to first order set. 
Another way is, for example, to use the *leapfrog* algorithm. 
If you are not comfortable with using a blackbox solver, this is the way for you. Try googling "Boris leapfrog" and implementing it in your code.

## Implementation
Download the [motion.jl](https://github.com/tungli/F5170-julia/blob/master/3_Motion/motion.jl) file and run Julia REPL.
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
This is fine, since DifferentialEquations.jl by default only solves for the necessary points and interpolates inbetween these point to obtain the whole solution.

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
Since the same equations and solvers apply here, let us jump to the implementation right away.

## Implementation
The script is named `van_allen.py`.
You can see that the implementation is basically the same as with the particle drift study.
Some differences here are:
 - We assume electric field is zero
 - Earth's magnetic field can be approximated as if produced by a (huge) magnetic dipole.


