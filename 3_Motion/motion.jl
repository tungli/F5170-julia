using Plots
pyplot()

using DifferentialEquations
using LinearAlgebra

function derivs!(dy,y,p,t)
    dy[1:3] = y[4:6]
    dy[4:6] = p[1]/p[2]*(p[3] + y[4:6] × p[4])
end

y0 = [0.0,0,0,0,0,0]
p = [1.0,1.0,[0.0,1,0],[0.0,0,1]]
tspan = (0.0,10.0)# range(0,stop=10,length=1000)
prob =  ODEProblem(derivs!,y0,tspan,p)

sol = solve(prob)

t = range(tspan[1],length=1000,stop=tspan[2])
a = sol(t)'[:,1:3]

x = a[:,1]
y = a[:,2]
z = a[:,3]

plot(x,y)

