using Plots
pyplot()

v = range(0,stop=4000,length=500)

function Fv(v)
    amu = 1.66053904e-27
    m = 40.0*amu
    kB = 1.3806485e-23
    T = 1000.0
    a = m/kB/T
    sqrt(2/pi)*v^2*a^(3.0/2.0)*exp(-a/2*v^2)
end

v_mean = 
v_sq_mean =
v_mp = 
v_mean_an = 
v_sq_mean_an = 
v_mp_an =

plot(v,Fv.(v),label="Maxwell-Boltzmann distribution")
scatter!([v_mean],[Fv(v_mean)])
scatter!([v_sq_mean],[Fv(v_sq_mean)])





