include("testmodule.jl")
using ClassicalOrthogonalPolynomials, SpecialFunctions, LinearAlgebra, Plots
α=0.0;β=0.0;b=-1;p=2.0;N=1000
ϵ=1.0
x=Inclusion(-1..1);
y=testmodule.y2x(p);
S=Jacobi(α,β);

# compute every component
I1=testmodule.OpI22_stable(α,β,b,p,0.5,N)[1:N+1,1:N+1]; # I^0.5
I2=(testmodule.OpI(α,β,b,p)^2)[1:N+1,1:N+1]; # I^2
X=testmodule.OpX(α,β,b,p)[1:N+1,1:N+1]; # x
fv=Complex{Float64}.(vcat(testmodule.OpC(BigFloat(α),β,b,p,5)\[ϵ*im^BigFloat(1.5)/gamma(BigFloat(1.5));0;0;1;0;1], zeros(N-5))) # coefficients of f
gv=vcat(testmodule.OpC(α,β,b,p,1)\[0;1],zeros(N-1)) # coefficients of 1

# solve system
@time temp=[ϵ*im^1.5*I1-X*I2 fv;0.5*S[1,1:N+1]'*I2 2]\[zeros(N+1);1]
vv=temp[1:end-1]
a=temp[end]
uv=I2*vv

# compute solutions on grid
ygrid=0:0.00002:1
xgrid=testmodule.y2x.(ygrid,p)
Qgrid=(1 .+ ygrid).^b .* S[ygrid,1:550]
sol=Qgrid*uv[1:550] .+ a*(xgrid.+1);

# plot solutions
plot(xgrid, [real.(sol),imag.(sol)], xlims=(-1,1), xlabel="\$x\$", ylabel="\$u(x)\$", size=(700,500))

# plot coefficients
plot(1:N+1, abs.(uv), yaxis=:log, xlabel="coefficient index", ylabel="absolute value of coefficients", size=(700,500))

plot(1:N+1, [real.(uv),imag.(uv)], xlabel="coefficient index", ylabel="coefficients", size=(700,500), labels=["real" "imag"])

cond([ϵ*im^1.5*I1-X*I2 fv;0.5*S[1,1:N+1]'*I2 2])