# Data processing

In this chapter, you will learn to use _*for*_-loops and _*if*_ statement blocks.

You will understand how to do those just by going over the examples, so I will not go into detail here.
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

 

