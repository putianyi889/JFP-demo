## (1-I^0.5+I-I^1.5+I^2)u=1
include("testmodule.jl")
using ClassicalOrthogonalPolynomials, ContinuumArrays, LinearAlgebra, SpecialFunctions, Plots
pgfplotsx()
α=0.0;β=0.0;b=0;p=2.0;N=30;
x=Inclusion(-1..1);
y=testmodule.y2x(p);
S=Jacobi(α,β);

# compute every component
I1=testmodule.OpI22_stable(α,β,b,p,0.5,N)[1:N+1,1:N+1]; # I^0.5
I2=testmodule.OpI(α,β,b,p,N+1); # I^1
I3=testmodule.OpI22_stable(α,β,b,p,1.5,N)[1:N+1,1:N+1]; # I^1.5
I4=(testmodule.OpI(α,β,b,p)^2)[1:N+1,1:N+1] # I^2
f=ones(x);
fv=(S \ ((1 .+ x).^(-b) .* f))[1:N+1]; # coefficients of f

# solve system
uv=(I-I1+I2-I3+I4)\fv;

# compute solution on grid
ygrid=-1:0.01:1 
xgrid=testmodule.y2x.(ygrid,p)
Qgrid=(1 .+ ygrid).^b .* S[ygrid,1:N+1]
sol=Qgrid*uv;

# plot solution
plot(xgrid, sol, xlims=(-1,1), xlabel="\$x\$", legend=false, size=(400,300), framestyle=:box, ylabel="\$u(x)\$")

# plot coefficients
plot(1:N+1, abs.(uv), yaxis=:log, xlabel="coefficient index", legend=false, size=(400,300), framestyle=:box, ylabel="absolute value of coefficients")


## u-erfc(sqrt(1+x))I^0.5u=1
include("testmodule.jl")
using ClassicalOrthogonalPolynomials, ContinuumArrays, LinearAlgebra, SpecialFunctions, Plots
pgfplotsx()
α=0.0;β=0.0;b=0;p=2.0;N=30;
x=Inclusion(-1..1);
y=testmodule.y2x(p);
S=Jacobi(α,β);

# compute every component
iop=testmodule.OpI22_stable(α,β,b,p,0.5,N)[1:N+1,1:N+1]; # I^0.5
m=erfc.(sqrt.(1 .+ y))
M=(S\(m.*S))[1:N+1,1:N+1]; # multiplication operator
f=ones(x);
fv=(S \ ((1 .+ x).^(-b) .* f))[1:N+1]; # coefficients of f

# solve system
uv=(I-M*iop)\fv;

# compute solution on grid
ygrid=-1:0.01:1
xgrid=testmodule.y2x.(ygrid,p)
Qgrid=(1 .+ ygrid).^b .* S[ygrid,1:N+1]
sol=Qgrid*uv;

# plot solution
plot(xgrid, sol, xlims=(-1,1), xlabel="\$x\$", legend=false, size=(400,300), framestyle=:box, ylabel="\$u(x)\$")

# plot coefficients
plot(1:N+1, abs.(uv), yaxis=:log, xlabel="coefficient index", legend=false, size=(400,300), framestyle=:box, ylabel="absolute value of coefficients")

# Bagley–Torvik
include("testmodule.jl")
using ClassicalOrthogonalPolynomials, SpecialFunctions, LinearAlgebra, Plots
pgfplotsx()
α=0.0;β=0.0;b=-1;p=2.0;N=30
x=Inclusion(-1..1);
y=testmodule.y2x(p);
S=Jacobi(α,β);

# compute every component
I1=testmodule.OpI22_stable(α,β,b,p,1.5,N)[1:N+1,1:N+1]; # I^1.5
I2=(testmodule.OpI(α,β,b,p)^2)[1:N+1,1:N+1]; # I^2
fvrl=vcat(testmodule.OpC(α,β,1)\[sqrt(2)/gamma(1/2);1], zeros(N-1)) # coefficients of f
fvc=vcat(testmodule.OpC(α,β,1)\[0;1], zeros(N-1))
gv=vcat(testmodule.OpC(α,β,3)\[0;0;1/gamma(3/2)/sqrt(2);0.5], zeros(N-3)) # coefficients of g

# solve system
temprl=[I+I1+I2 gv; 0.5*S[1,1:N+1]'*I2 2] \ [-fvrl; -1] # [v;a]
tempc=[I+I1+I2 gv; 0.5*S[1,1:N+1]'*I2 2] \ [-fvc; -1]
vvrl=temprl[1:end-1]; arl=temprl[end]
vvc=tempc[1:end-1]; ac=tempc[end]
uvrl=I2*vvrl # coefficients of I^2*v
uvc=I2*vvc

# compute solutions on grid
ygrid=-1:0.01:1
xgrid=testmodule.y2x.(ygrid,p)
Qgrid=(1 .+ ygrid).^b .* S[ygrid,1:N+1]
solrl=Qgrid*uvrl .+ arl*(xgrid.+1) .+ 1;
solc=Qgrid*uvc .+ ac*(xgrid.+1) .+ 1;

# plot solutions
plot(xgrid, [solrl solc], xlims=(-1,1), xlabel="\$x\$", size=(400,300), framestyle=:box, ylabel="\$u(x)\$", labels=["RL type" "Caputo type"], legend=:bottomleft, linestyle=[:solid :dash])

# plot coefficients of I^2*v
plot(1:N+1, abs.([uvrl uvc]), yaxis=:log, xlabel="coefficient index", labels=["RL type" "Caputo type"], size=(400,300), framestyle=:box, ylabel="absolute value of coefficients", legend=:bottomleft, linestyle=[:solid :dash])