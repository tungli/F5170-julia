using Dierckx
using DelimitedFiles

data = readdlm("sigmaion.dat",',')
x = data[:,1]
y = data[:,2]

function Fv(v)
    m = 9.109e-31
    kB = 1.3806485e-23
    T = 11604.0*15.0
    a = m/kB/T
    sqrt(2/pi)*v^2*a^(3.0/2.0)*exp(-a/2*v^2)
end

gaussian(x,x0,sigma) = exp(-(x-x0)^2/(2*sigma^2))/sqrt(2*pi)/sigma

function trapz(y::AbstractVector; x=1:length(y))
    #Trapezoidal rule
    dx = diff(x)
    sum([ (y[i-1] + y[i])/2*dx[i-1] for i in 2:length(y)])
end

v = range(0.0,stop=1e7,length=10^6)
cr_sec = Spline1D(x,y,k=3,bc="zero")

@show kr_MB = trapz(v.*Fv.(v).*cr_sec.(v),x=v)
@show vmean = trapz(Fv.(v).*v,x=v)
@show kr_beam = trapz(v.*gaussian.(v,vmean,vmean/100).*cr_sec.(v),x=v)



