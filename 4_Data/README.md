# Data processing
In this chapter, you will learn to use _*for*_-loops and _*if*_ statement blocks.

You will understand how to do those just by going over the examples, so I will not go into detail here.

The data files you will be working here (`csk[1-3].dat`) contain data for the cross sections of collision (in m<sup>2</sup>) of electrons with argon atom species for different temperatures (in Kelvins).
In particular, `csk1.dat` contains the data for elastic collision:

![coll](http://mathurl.com/ycnhzk89.png)

`csk2.dat` for argon excitation:

![exci](http://mathurl.com/ybd2s6ql.png)

`csk3.dat` for argon ionization:

![ioni](http://mathurl.com/ydyupuzm.png)

The scripts referenced in the exercises are described in subsequent sections.

## Exercises
>  **Exercise 1**
>  * Modify [this script](https://github.com/tungli/F5170-julia/blob/master/4_Data/simple_plot.jl) so that it plots:
>  
>  ![sinc](http://mathurl.com/y983ysyp.png)
>  
>  **Exercise 2**
>  [This script](https://github.com/tungli/F5170-julia/blob/master/4_Data/data_plot.jl) takes the data in `csk[1-3].dat` and plots them.
>  * Run the script with the data in the same directory
>  * You can see that the magnitudes are too different and the plot is not very practical. Your task is to modify it so that:
>    1. The electron energy is expressed in electronvolts (eV).
>    2. The y-axis scale is logarithmic
>  
>  The logarithmic plot should look similar to this:
>  ![Data](https://github.com/tungli/F5170-python/blob/master/4_Data/data_plot.svg)
>  
>  * Look at your plot -- what are the excitation and ionization thresholds?
>  
>  **Exercise 3**
>  * Modify the script so that it converts the *x*-data from Kelvin to electronvolt and saves each cross-section to tab-delimited file named `csevN.dat` with `N` being the corresponding file number.
>  
>  **Exercise 4**
>  * Find the function in Julia which calculates the inverse of a matrix. 
>  * Define a [singular matrix](http://mathworld.wolfram.com/SingularMatrix.html) and try calculating its inverse. What happens?
>  * Run [this script](https://github.com/tungli/F5170-julia/blob/master/4_Data/inverse_matrix.jl) and verify that it works correctly by testing it on various matrices.
>  * Add another `elif` statement so that it displays a warning when the [matrix rank](http://mathworld.wolfram.com/MatrixRank.html) is greater than 10.
>  
>  **Advanced exercise**
>  In this directory you will also find data files `adv_csk[1-3].dat`. Two of these files include electron temperature in electronvolts (eV) while one of them includes electron energy in Kelvin.
>  * Modify the [script](https://github.com/tungli/F5170-python/blob/master/4_Data/data_plot.py) using *if-else* statements so that it decides which files should be converted and which not.
>  * Such a script could be very useful when processing thousands of similar files. However, what are the limitations of your program?



Let us take a look at the scrips you will need for the exercises.

## Simple plot
```julia
using Plots
using LaTeXStrings
pyplot()

colors =  [:red,:green,:blue,:magenta,:black,:yellow,:cyan]
x = range(0,stop=1,length=1000)

plot()
for i in 1:length(colors)
    plot!(x,x.^(i-1),c=colors[i],label = latexstring("x^{$(i)}"))
end
Plots.gui()
```
In the beginning we import a package for plotting and a package that allows us to use LaTeX inside strings.
We set PyPlot as our backend for plots.
Then, we define a `Vector` of `Symbol`s.
These are arguments used in `plot()` to define the color of the line.
Next, we define the *x* coordinate for the plot.

Then we create a for-loop which iterates over the numbers 1,2,3...
Inside the for-loop we make a simple plot with a legend.
Notice that in the Plots.jl package, `plot()` creates a new plot, while `plot!()` adds to the current plot.

**Note:** Here is one other way we could have designed the for loop -- using the `enumerate()` function.
```julia
plot()
for (i,c) in enumerate(colors):
    plot!(x,x^(i-1),c=c,label = latexstring("x^{$(i)}"))
end
```
Maybe an interesting thing is to notice is `c=c` - while the first `c` is a defined keyword inside the `plot()` function, the second `c` is the variable containing a `Symbol`


After the for-loop we tell Julia it can show us the figure.

## Loading and plotting a matrix

```julia
using Plots
pyplot()

colors = [:red,:green,:blue,:magenta,:black,:yellow,:cyan]

plot()
for i in 1:3
    filename = "csk$(i).dat"
    f = open(filename)
    data = readlines(filename)
    close(f)

    x = Vector{Float64}()
    y = Vector{Float64}()
    for i in data
        a = parse.(Float64,split(i,'\t'))
        push!(x,a[1])
        push!(y,a[2])
    end
    
    #Removal of x,y pair if y = 0
    t = .!(iszero.(y))
    x = x[t]
    y = y[t]

    plot!(x,y,c=colors[i],label=filename,xlim=([0, 1e6]))
end
Plots.gui()
```
The first steps are same as in the "Simple Plot" case.
In the *for*-loop, we define `filename` using a string and the iteration variable `i`.
In Julia, the symbol `$()` places a variable into the string.
Then we continue by opening the file (creating a stream variable `f`) and loading everything in an `Vector` of `String`s.
We close the file stream since we have everything in the variable `data` now.
To separate the first and the second column, we define empty arrays `x` and `y`.
Then we loop over the data lines a push the values into the arrays.
We use `parse` to convert `String`s into `Float`s.

Another (shorter) option to load a data file is to use:
```julia
using DelimitedFiles

data = readdlm("csk1.dat")
x = data[:,1]
y = data[:,2]
```
The `readdlm()` function does the split part as well as the conversion to floats.

Then I added removal of value where `y` is zero so that you can do a log-plot without errors.

## Inverse matrix
```julia
using LinearAlgebra

function myinv(mat::Matrix)
    s = size(mat)
    if s[1] != s[2]
        println("Not a square matrix!")
        return nothing
    elseif det(mat) == 0
        println("Matrix is singular!")
        return nothing
    else
        return inv(mat)
    end
end
```
We created a more sophisticated version of matrix inversion.
It notifies us when the matrix is not square or when its determinant is zero and in these cases, it does not try to call the matrix inversion funcion.
Couple of things to notice:
 * `size()` returns a `Tuple` (an object similar to `Array`) with the array's size.
 * `if`,`elseif` are followed by a condition. A condition has a value of `true` or `false`.
 * `end` is used to end the if-block
 * `return` keyword exits the function with the desired value.

 

