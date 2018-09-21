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
