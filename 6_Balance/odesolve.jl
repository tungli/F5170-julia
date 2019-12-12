using DifferentialEquations
using Plots
pyplot()
using LaTeXStrings

p = 1e5
kB = 1.38e-23
Tg = 400.0
Te = 2.0

tsteps = 1000
tout = 10.0.^range(-11,stop=-6,length=tsteps)
tspan = (minimum(tout),maximum(tout))

nArs_init = 1e12
nArp_init = 1e12
nAr2p_init = 1e12
ne_init = nArp_init + nAr2p_init
initial = [nArs_init, nArp_init, nAr2p_init, ne_init]

function odefun!(dy,y,params,t)
    nArs = y[1]
    nArp = y[2]
    nAr2p = y[3]
    ne = y[4]
    nAr = p/(kB*Tg) - nArs - 0.5*nAr2p - nArp


    k1 = f_k1(Te,Tg)    
    k2 = f_k2(Te,Tg)
    k3 = f_k3(Te,Tg)
    k4 = f_k4(Te,Tg)
    k5 = f_k5(Te,Tg)
    k6 = f_k6(Te,Tg)
    k7 = f_k7(Te,Tg)
    k8 = f_k8(Te,Tg)
    k9 = f_k9(Te,Tg)
    k10 = f_k10(Te,Tg)
    k11 = f_k11(Te,Tg)

    dy[1] = k1*ne*nAr - k2*ne*nArs - k4*ne*nArs + k6*ne*nAr2p - 2*k10*nArs^2 - k11*nArs*nAr
    dy[2] = 1.0
    dy[3] = 1.0
    dy[4] = 1.0

    return nothing
end

#You need to define these functions
f_k1(Te, Tg) = 
f_k2(Te, Tg) = 
f_k3(Te, Tg) = 
f_k4(Te, Tg) = 
f_k5(Te, Tg) = 
f_k6(Te, Tg) = 
f_k7(Te, Tg) = 
f_k8(Te, Tg) = 
f_k9(Te, Tg) = 
f_k10(Te, Tg) =
f_k11(Te, Tg) =

prob = ODEProblem(odefun!,initial,tspan,Vector{Float64}())
sol = solve(prob)

nAr = p/(kB*Tg) .- [ i[1] - 0.5*i[3] - i[2] for i in sol.u ]
plot(sol,yscale=:log10,xscale=:log10,vars=(0,1),label=L"Ar^*")
plot!(sol,yscale=:log10,xscale=:log10,vars=(0,2),label=L"Ar^+")
plot!(sol,yscale=:log10,xscale=:log10,vars=(0,3),label=L"Ar_2^+")
plot!(sol,yscale=:log10,xscale=:log10,vars=(0,4),label=L"e^-")
plot!(range(tspan[1],stop=tspan[2],length=length(sol)),nAr,label="Ar")

