using Plots
pyplot()

v = range(-2000,stop=2000,length=200)
t = range(0,stop=1e-8,length=100)
amu = 1.66053904e-27
m = 40.0*amu
kB = 1.3806485e-23
T = 1000.0
coll_freq = 5e8

gaussian(x,x0,sigma) = exp(-(x-x0)^2/(2*sigma^2))/sqrt(2*pi)/sigma

fx0 = gaussian.(v,500.0,10.0)
fy0 = zeros(length(v))
fz0 = zeros(length(v))

function trapz(y::AbstractVector; x=1:length(y))
    #Trapezoidal rule
    dx = diff(x)
    sum([ (y[i-1] + y[i])/2*dx[i-1] for i in 2:length(y)])
end

v_sq = sqrt(trapz((fx0+fy0+fz0).*v.^2,x=v))

Teq = m*v_sq^2/kB/3

fx_eq = sqrt.(m/(2*pi*kB*Teq))*exp.(-m*v.^2/(2*kB*Teq))
fy_eq = sqrt.(m/(2*pi*kB*Teq))*exp.(-m*v.^2/(2*kB*Teq))
fz_eq = sqrt.(m/(2*pi*kB*Teq))*exp.(-m*v.^2/(2*kB*Teq))

calc_distribution(feq,f0,t,coll_f) = feq + (f0 - feq)*exp(-t*coll_f) 

@gif for i in t
    px = plot(v,calc_distribution(fx_eq,fx0,i,coll_freq),ylim=([0,maximum(fx0)]))
    py = plot(v,calc_distribution(fy_eq,fy0,i,coll_freq),ylim=([0,maximum(fx0)]))
    pz = plot(v,calc_distribution(fz_eq,fz0,i,coll_freq),ylim=([0,maximum(fx0)]))
    plot(px,py,pz,layout=(1,3),legend=false)
end
