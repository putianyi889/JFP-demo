include("testmodule.jl")
using Plots, ClassicalOrthogonalPolynomials, LinearAlgebra, LaTeXStrings
pgfplotsx()

α=-1.0;β=1.0;b=0;p=2;μ=0.5;N=500; # parameters
xgrid=-1:0.01:1;
ygrid=testmodule.x2y.(xgrid,p);
S=Jacobi(α,β);
@time iop=testmodule.OpI22_stable(α,β,b,p,μ,N); # FI operator
f=ones(Inclusion(-1..1)); # RHS function
fv=(S \ ((1 .+ Inclusion(-1..1)).^(-b) .* f))[1:N+1]; # coefficients of f
Qgrid=(1 .+ ygrid).^b .* S[ygrid,1:N+1]; # evaluation of Q on grid
λ=[1 2 4 8 16 32];
err=zeros(N,length(λ)); # max abs error w.r.t. truncation size and λ
uv=zeros(N+1,length(λ)); # coefficients of solution w.r.t. λ
truncsize=zeros(length(λ)); # truncation size needed for convergence
condnumb=zeros(N,length(λ)); # condition numbers w.r.t. truncation size and λ

for k in 1:length(λ)
    uv[:,k]=(I+λ[k]^2*iop[1:N+1,1:N+1])\fv
    u=testmodule.mittleff_matlab.(μ,-λ[k]^2*(0:0.01:2).^μ)
    for n in 1:N
        print(n," ")
        err[n,k]=maximum(abs.(u-Qgrid[:,1:n]*uv[1:n,k]))
        condnumb[n,k]=cond(I+λ[k]^2*iop[1:n,1:n])
    end
    for n in 1:N
        if err[n,k]<1e-14
            truncsize[k]=n
            break
        end
    end
end

plot(abs.(uv), ylims=(1e-17,1), yaxis=:log, legend=:topright, size=(300,250), yticks=10.0.^(-16:4:0),xlabel="coefficient index", ylabel="coefficient value", linestyle=[:solid :dash :dot :solid :dash :dot], linewidth=[1 1 1 2 2 2], labels=latexstring.("\\lambda=",λ))
plot(abs.(err), ylims=(1e-17,1), yaxis=:log, legend=:topright, size=(300,250), yticks=10.0.^(-16:4:0),xlabel="truncation size", ylabel="maximum error", linestyle=[:solid :dash :dot :solid :dash :dot], linewidth=[1 1 1 2 2 2], labels=latexstring.("\\lambda=",λ))

plot(condnumb,yaxis=:log2, yticks=2.0.^(0:2:10), size=(300,250), xlabel="truncation size", ylabel="condition number", linestyle=[:solid :dash :dot :solid :dash :dot], linewidth=[1 1 1 2 2 2], labels=latexstring.("\\lambda=",λ), legend=:bottomright)

plot(λ',truncsize,xaxis=:log2, yaxis=:log2, xticks=λ', legend=false, size=(300,250), xlabel=L"\lambda", yticks=16*λ', markers=true, markersize=2, ylabel="truncation size")
plot(λ',condnumb[end-1,:],xaxis=:log2, yaxis=:log2, xticks=λ', legend=false, size=(300,250), xlabel=L"\lambda", yticks=2 .^(1:2:11), markers=true, markersize=2, ylabel="condition number")