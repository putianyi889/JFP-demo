include("testmodule.jl")
using Plots, LinearAlgebra, ClassicalOrthogonalPolynomials, SpecialFunctions, LaTeXStrings
pgfplotsx()
α=-0.5;β=-0.5;b=-0.5;p=3;μ=1/3;N=50;
xaxis=Inclusion(-1..1);xgrid=-1:0.01:1;
yaxis=testmodule.y2x(float(p));ygrid=testmodule.x2y.(xgrid,p);
fx=(1 .+ xaxis).^0.5; fy=fx[yaxis];
S=Jacobi(α,β);
Q=(1 .+ ygrid).^b .* S[ygrid,1:N+1];
fv=(S \ ((1 .+ Inclusion(-1..1)).^(-b) .* fy))[1:N+1];
print("computing mittleff:");
@time sol=gamma(1.5).*sqrt.(xgrid.+1).*Float64.(testmodule.mittleff_matlab.(μ,1.5,-(xgrid.+1).^BigFloat(μ)));
print("computing iop:");
@time iop=testmodule.OpI22_stable(α,β,b,p,μ,N);
err=zeros(N);
uv=(I+iop[1:N+1,1:N+1])\fv;
for n in 1:N
    err[n]=maximum(abs.(sol[2:end]-Q[2:end,1:n]*uv[1:n]))
end

# error and coefficients
plot(ylims=(1e-16,1), yaxis=:log, yticks=10. .^(0:-4:-16), xlims=(0,40), legend=:topright, size=(300,250), xlabel="truncation size / coefficient index")
plot!(1:N,err, label="maximum error")
plot!([1],[err[1]], label="on -1:0.01:1", color=:transparent)
plot!(1:N,abs.(uv[1:N]), label="absolute value", linestyle=:dash, color=:red)
plot!([1],[err[1]], label="of coefficient", color=:transparent)

# uniform convergence
plot(xgrid, sol-Q[:,1:N+1]*uv, size=(300,250), xlabel="\$x\$", legend=false, yticks=(1e-16*(-3:3:6),[k==0 ? k : "$(k)e-16" for k=-3:3:6]), xlims=(-1,1), framestyle=:box)

# solution
plot(xgrid,sol, xlabel="\$x\$", size=(300,250), framestyle=:box, xlims=(-1,1), ylims=(0,0.7), legend=false, yticks=0.1*(0:7), ylabel="\$u(x)\$")