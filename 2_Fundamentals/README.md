# Julia Basics

Julia is very intuitive and you will find everything you will need and more in the official [documentation](https://docs.julialang.org/en/v1/)
Here are just a few things to get you started.

To download a package use:
```julia
using Pkg
Pkg.add("Plots")
```
to use it:
```julia
using Plots
```
[Plots.jl](https://github.com/JuliaPlots/Plots.jl) offers a powerful plotting interface, although it is a bit slow to start. 
You can select the backend for plots - call `pyplot()` or `gr()` after `using Plots` line.

Using this in Julia REPL (run the `julia` executable) can be helpful:
```julia
?range
```

Julia arrays are indexed from one:
```julia
a = [1,2,3]
println(a[1])

A = [1 2 3; 4 5 6]
println(a[:,1])
```
There is an `end` keyword:
```julia
println(a[2:end])
println(a[1:(end-1)])
```
To use matrix mutliplication you don't need special functions:
```
A = rand(3,3)
a = [1,2,0]
println(A*a)
println(a'*A)
```
Vectors are columns and `'` is the transpose function.

Functions in Julia are usually **not** written for vectors.
It is better to write functions for single values, then use *broadcasting* to apply it to vector.
Broadcasting is done by placing a dot between the function name and the argument parenthesis:
```julia
a = 1
sin(a)
a = [1,2,3]
sin.(a)
```

Now try to do the exercises 2.1 and 2.2 from the MATLAB material.

# Extra problems

## Vector product
Write a function that accepts 5 arguments: *q,m,v,E,B*, and calculates the Lorentz force.
Browse the internet to find the function which calculates the [vector product](https://en.wikipedia.org/wiki/Cross_product) of two 3-element vectors and implement it in your function.
Test your function. 

## Debye length
One of the important parameters of a plasma is the [Debye length](https://en.wikipedia.org/wiki/Debye_length).
In plasma it can be calculated from electron number density and electron temperature:

![Debye](http://mathurl.com/y876kcbb.png)

First, implement a function which calculates the Debye length.
Then, download the `Discharges.txt` file. 
Ideally, you want to have the file in the same directory where you are working.

Working with delimited files in Julia is easy.
Take a look at the `readdlm()` and `parse()` functions:
```julia
?readdlm
?parse
```
Use it to import the data from `Discharges.txt`.
To convert the type `String` into `Float64` use the `parse()` function.

## Larmor radius
What will be the trajectory of a charged partice with an initial velocity in homogeneous magnetic field?
Try to guess what the following parameter, called *Larmor radius*, represents in the trajectory:

![Larmor](http://mathurl.com/ybs37jkj.png)

Download the `result.txt` file which contains the data from such trajectory.
In this file there are 3 columns -- 1. column is time, 2.column is the *x*-position and 3. column the *y*-position.

Using `readdlm()` import the data to Julia.

Now, plot *x* as a function of time and *y* as a function of *x*.

Suggest how to get an estimate of the Larmor radius from the data.
Implement the method.  

Assume the charged particle in the trajectory was an electron with an intial velocity of 1 m/s.
Calculate the magnetic field.


