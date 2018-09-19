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
We will be using the brilliant **DifferentialEquations.jl** package. 
Take a look at the [tutorial](http://docs.juliadiffeq.org/latest/tutorials/ode_example.html#Defining-Parameterized-Functions-1).
By following the same procedure (using the `@ode_def` macro (take a look at [macro-related magic](https://docs.julialang.org/en/v1/manual/metaprogramming/)) we solve our problem. Let us look at it in parts:
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
With this macro we define the ODEs. If you open the file `motion.jl` you can notice that in Julia we can make use Unicode characters to make the code look very similar to the mathematical formulation (I do not know how to use Unicode in Markdown here).

```julia
y0 = [0.0,0,0,0,0,0]
p = [1.0,1.0,[0.0,1,0],[0.0,0,1]]
tspan = range(0,stop=10,length=1000)
prob =  ODEProblem(f,y0,tspan,p)
```


# Van Allen radiation belt
Since the same equations and solvers apply here, let us jump to the implementation right away.

## Implementation
The script is named `van_allen.py`.
You can see that the implementation is basically the same as with the particle drift study.
Some differences here are:
 - We assume electric field is zero
 - Earth's magnetic field can be approximated as if produced by a (huge) magnetic dipole.

```python
from scipy import integrate as itg
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
```
Here some extra things are needed for 3D plots.

```python
def derivs(y,t,E,q,m,magM):
    d = np.zeros(np.shape(y))
    d[0:3] = y[3:6]
    d[3:6] = q/m*(E + np.cross(y[3:6],bfield(y[0:3],magM)))
    return d

def bfield(r,M):
    #This function calculates the magnetic field from 
    #position vector and magnetic dipole moment
    #Good luck!
```
We replace the magnetic field vector with a magnetic field vector generating function (which you are meant to write yourself).

```python
E = np.array([0.0,0.0,0.0])
magM = np.array([0.0,0.0,8.10e22])
q = 1.60217662e-19
m = 1.6726219e-27 
c = 3e8
re = 6.38e6
Ek_ev = 5e7

#Velocity
vr = c/np.sqrt(1.0+m*c**2/Ek_ev/np.abs(q))
vp = 0.0
vt = np.pi/4
v = np.array([vr*np.sin(vt)*np.cos(vp),vr*np.sin(vt)*np.sin(vp),vr*np.cos(vt)])
#Position
rr = 2.5*re
rp = 0.0
rt = np.pi/2
r = np.array([rr*np.sin(rt)*np.cos(rp),rr*np.sin(rt)*np.sin(rp),rr*np.cos(rt)])
#State vector
y0 = np.array([r[0],r[1],r[2],v[0],v[1],v[2]])

ti = 0.0
tf = 25.0
num_points = 10000
t = np.linspace(ti,tf,num_points)
```
Constants and initial values. Spherical coordinates for velocity and position vectors are used, which are then converted to Cartesian coordinates.
This is convenient because of the spherical nature of the Earth.


```python
res = itg.odeint(derivs,y0,t,args=(E,q,m,magM))

x = [ i[0]/re for i in res ]
y = [ i[1]/re for i in res ]
z = [ i[2]/re for i in res ]
```
Integration, followed by extraction of coordinates.

```python
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
xs = np.cos(u)*np.sin(v)
ys = np.sin(u)*np.sin(v)
zs = np.cos(v)
ax.plot_surface(xs, ys, zs, cmap=cm.plasma)
#Note: the colormap "plasma" is the only acceptable colormap when doing plasma physics!
#...although "jet" can be justified in some cases.
```
Plot the ball...

```python
ax.plot(x, y, z, lw=1.0, c="b")
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")
#ax.set_xlim((-3.0,3.0))
#ax.set_ylim((-3.0,3.0))
#ax.set_zlim((-3.0,3.0))

plt.show(block=True)
```
...and the trajectory.









