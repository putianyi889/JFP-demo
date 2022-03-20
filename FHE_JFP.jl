include("testmodule.jl")
using Plots, ClassicalOrthogonalPolynomials, LinearAlgebra

α=-1.0;β=1.0;b=0;p=2;μ=0.5;N=500; # parameters
xgrid=-1:0.01:1;
ygrid=testmodule.x2y.(xgrid,p);
S=Jacobi(α,β);
iop=testmodule.OpI22_stable(α,β,b,p,μ,N); # FI operator
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

plot(condnumb,yaxis=:log)
plot(λ',condnumb[end-1,:],xaxis=:log, yaxis=:log, markers=true)