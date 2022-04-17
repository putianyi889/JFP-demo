# u+I^μ[u]=1
include("testmodule.jl")
using Plots, LinearAlgebra, ClassicalOrthogonalPolynomials, SpecialFunctions, LaTeXStrings
pgfplotsx()
α=0.0;β=0.0;b=0;N=50;
xgrid=-1:0.01:1;
f=ones(Inclusion(-1..1));
setprecision(80); # for Mittag-Leffler function

# μ=1/p for different p
P=[1 2 4 31]; # different p
S=Jacobi(α,β);
fv=(S \ ((1 .+ Inclusion(-1..1)).^(-b) .* f))[1:N+1];
err=zeros(N,length(P));
for k=1:length(P)
    p=P[k];
    μ=1/p;
    ygrid=testmodule.x2y.(xgrid,p);
    Q=(1 .+ ygrid).^b .* S[ygrid,1:N+1];
    print("computing iop for p=$p:")
    @time iop=testmodule.OpI22_stable(α,β,b,p,μ,N);
    uv=(I+iop[1:N+1,1:N+1])\fv
    print("computing mittleff for p=$p:")
    @time sol=Float64.(testmodule.mittleff_matlab.(μ,-(xgrid.+1).^BigFloat(μ)));
    for n in 1:N
        err[n,k]=maximum(abs.(sol-Q[:,1:n]*uv[1:n]))
    end
end
plot(err, xlims=(0,30), ylims=(1e-16,1), yaxis=:log, size=(300,250), linestyle=[:dashdot :dash :solid :dot], labels=string.("\$\\mu=1/",P,", p=",P,"\$"), yticks=10. .^(0:-4:-16), legend=:topright)

# irrational μ, p=1/μ
U=[exp(1)/π π/2 1/π^2]; # different μ
S=Jacobi(α,β);
fv=(S \ ((1 .+ Inclusion(-1..1)).^(-b) .* f))[1:N+1];
err=zeros(N,length(U));
for k=1:length(U)
    μ=U[k];
    p=1/μ;
    ygrid=testmodule.x2y.(xgrid,p);
    Q=(1 .+ ygrid).^b .* S[ygrid,1:N+1];
    print("computing iop for p=$p:")
    @time iop=testmodule.OpI11_stable(α,β,b,p,μ,N);
    uv=(I+iop[1:N+1,1:N+1])\fv
    print("computing mittleff for p=$p:")
    @time sol=Float64.(testmodule.mittleff_matlab.(μ,-(xgrid.+1).^BigFloat(μ)));
    for n in 1:N
        err[n,k]=maximum(abs.(sol-Q[:,1:n]*uv[1:n]))
    end
end
plot(err, xlims=(0,30), ylims=(1e-16,1), yaxis=:log, size=(300,250), linestyle=[:solid :dash :dot], yticks=10. .^(0:-4:-16), labels=["\$\\mu=e/\\pi\$, \$p=\\pi/e\$" "\$\\mu=\\pi/2\$, \$p=2/\\pi\$" "\$\\mu=1/\\pi^2\$, \$p=\\pi^2\$"], legend=:topright)

# fix p, pμ=1:2*p
p=3;
U=(1:2*p)./p; # different μ
ygrid=testmodule.x2y.(xgrid,p);
S=Jacobi(α,β);
Q=(1 .+ ygrid).^b .* S[ygrid,1:N+1];
fv=(S \ ((1 .+ Inclusion(-1..1)).^(-b) .* f))[1:N+1];
err=zeros(N,length(U));
for k=1:length(U)
    μ=U[k];
    print("computing iop for μ=$μ:")
    @time iop=testmodule.OpI22_stable(α,β,b,p,μ,N);
    uv=(I+iop[1:N+1,1:N+1])\fv
    print("computing mittleff for μ=$μ:")
    @time sol=Float64.(testmodule.mittleff_matlab.(μ,-(xgrid.+1).^BigFloat(μ)));
    for n in 1:N
        err[n,k]=maximum(abs.(sol-Q[:,1:n]*uv[1:n]))
    end
end
plot(err, xlims=(0,30), ylims=(1e-16,1), yaxis=:log, size=(300,250), legend=:topright, 
    linewidth=[2 2 2 1 1 1],
    linestyle=[:solid :dash :dot :solid :dash :dot],
    labels=string.("\$\\mu=",["1/3" "2/3" "1" "4/3" "5/3" "2"],"\$"),
    yticks=(10. .^(0:-4:-16))
)

# fix μ, different p
μ=0.5;
P=2:2:10; # different p
S=Jacobi(α,β);
fv=(S \ ((1 .+ Inclusion(-1..1)).^(-b) .* f))[1:N+1];
print("computing mittleff for μ=$μ:")
@time sol=Float64.(testmodule.mittleff_matlab.(μ,-(xgrid.+1).^BigFloat(μ)));
err=zeros(N,length(P));
for k=1:length(P)
    p=P[k];
    ygrid=testmodule.x2y.(xgrid,p);
    Q=(1 .+ ygrid).^b .* S[ygrid,1:N+1];
    print("computing iop for p=$p:")
    @time iop=testmodule.OpI22_stable(α,β,b,p,μ,N);
    uv=(I+iop[1:N+1,1:N+1])\fv
    for n in 1:N
        err[n,k]=maximum(abs.(sol-Q[:,1:n]*uv[1:n]))
    end
end
plot(err, ylims=(1e-16,1), yaxis=:log, xlims=(1,40), size=(300,250), legend=:topright, yticks=10. .^(0:-4:-16), linestyle=[:solid :dash :dashdot :dashdotdot :dot], labels=string.("\$p=",P',"\$"))