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
    filter!(x->x!=0.0,y)

    plot!(x,y,c=colors[i],label=filename,xlim=([0, 1e6]))
end
Plots.gui()
