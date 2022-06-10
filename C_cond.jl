include("testmodule.jl")
using Plots, LinearAlgebra, GenericLinearAlgebra, ThreadPools, LaTeXStrings
pgfplotsx()

α=BigFloat(-0.5); β=BigFloat(-0.5)
N=100
condition=zeros(N)
op=testmodule.OpC(α,β,N)
@qbthreads for n in 1:N
    print(n," ")
    condition[n]=cond(op[1:n,1:n])
end

plot(condition,yaxis=:log, xlabel=L"N", ylabel=L"cond$\left(C^{(\alpha,\beta)}_{1:N,1:N}\right)$", size=(300,250), yticks=10. .^(0:20:100), legend=false)