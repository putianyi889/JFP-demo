include("testmodule.jl")
using Plots, ApproxFun, LinearAlgebra, GenericLinearAlgebra, BandedMatrices, ThreadPools, JLD

setprecision(4096)
Sb = Jacobi(BigFloat(0.0),1.0) ⊕ JacobiWeight(BigFloat(0.5),0.,Jacobi(BigFloat(0.5),0.5));
Qb = LeftIntegral(Sb,BigFloat(0.5));
condnumb=zeros(100,10)
for n=1:10
    @qbthreads for N=1:100
        println("N=$N, n=$n")
        condnumb[N,n]=cond(Matrix((I+n*Qb)[1:N,1:N]))
    end
end
plot(condnumb,yaxis=:log)

n=25
NN=1:2500
# testcondnumb=zeros(length(NN))
op=(I+n*Qb)[1:NN[end],1:NN[end]];
# @time B=testmodule.bidiagonalize!(op);
# B=BandedMatrix(0=>B.dv,1=>B.ev);
# C=B'*B;
@qbthreads for k in 1056:length(NN)
    N=NN[k]
    # @time v=eigvals(SymTridiagonal(C[1:N,1:N]))
    # testcondnumb[k]=sqrt(v[end]/v[1])
    # @time testcondnumb[k]=cond(Matrix(op[1:N,1:N]))
    @time begin
        BB=testmodule.bidiagonalize!(op[1:N,1:N])
        BB=BandedMatrix(0=>BB.dv,1=>BB.ev);
        CC=SymTridiagonal(BB'*BB);
        v=eigvals(CC);
        testcondnumb[k]=sqrt(v[end]/v[1])
        print("N=$N, n=$n ")
    end
end
plot(NN,testcondnumb, yaxis=:log, legend=:topleft)

test=cond(Matrix(op))

@load "cond_sumspace_1.jld"
@load "cond_sumspace_2.jld"
@load "cond_sumspace_3.jld"
@load "cond_sumspace_4.jld"
@load "cond_sumspace_5.jld"
plot(cond4,yaxis=:log,legend=false)