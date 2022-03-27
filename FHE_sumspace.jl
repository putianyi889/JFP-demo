# Fractional heat equation (sumspace)
include("testmodule.jl")
using Plots, ApproxFun, LinearAlgebra, SpecialFunctions, Beep, GenericSVD, ThreadPools
using Plots.PlotMeasures
pgfplotsx()

# Float64
N=600;
S = Jacobi(0.0,0.0) ⊕ JacobiWeight(-0.5,0.,Jacobi(-0.5,0.5));
Q = LeftIntegral(S,0.5);
y=zeros(N,5);
for n=1:5
    @time yy=((I+n^2*Q) \ 1).coefficients
    y[1:min(length(yy),N),n]=yy[1:min(length(yy),N)]
end
plot(abs.(y), size=(400,250), xlabel="coefficient index", ylabel="coefficient value (abs)", ylims=(1e-17,1e14), yaxis=:log, legend=:topright, label=false)

# BigFloat
Nb=5000;
setprecision(3072);
Sb = Jacobi(BigFloat(0.0),1.0) ⊕ JacobiWeight(BigFloat(0.5),0.,Jacobi(BigFloat(0.5),0.5));
Qb = LeftIntegral(Sb,BigFloat(0.5));
yb=zeros(BigFloat,Nb,5);
for n=1:5
    @time yy=((I+n^2*Qb) \ BigFloat(1)).coefficients
    yb[1:min(length(yy),Nb),n]=yy[1:min(length(yy),Nb)]
end
yt=-800:200:600;
plot(log10.(abs.(yb)), size=(400,250), xlabel="coefficient index", ylabel="coefficient value (abs)", legend=:topright, label=false, yticks=(yt,string.("10^{",yt,"}")))

# Sketch of solution
xgrid=vcat(-0.99999,-0.99:0.01:1);
f1=Fun(S,y[:,1]).(xgrid));
f2=Fun(S,y[:,2]).(xgrid);
f3=Fun(S,y[:,3]).(xgrid);
@time fb1=Fun(Sb,yb[:,1]).(BigFloat.(xgrid));
@time fb2=Fun(Sb,yb[:,2]).(BigFloat.(xgrid));
@time fb3=Fun(Sb,yb[:,3]).(BigFloat.(xgrid));
plot(-1:0.01:1,[f1 f2 f3 fb1 fb2 fb3], xlims=(-1,1), ylims=(0,1), label=false, size=(300,250), linestyle=[:solid :solid :solid :dash :dash :dash], color=[:blue :red :green :blue :red :green], legend=:topright);
annotate!(0.5,0.42,"\$\\lambda=1\$");
annotate!(0.0,0.2,"\$\\lambda=2\$");
annotate!(-0.8,0.06,"\$\\lambda=3\$");
plot!([0 0],[-1 -1], color=:black, label=["numerical" "exact"], linestyle=[:solid :dash])

# Condition
using Plots, ApproxFun, LinearAlgebra
setprecision(896)
Sb = Jacobi(BigFloat(0.0),1.0) ⊕ JacobiWeight(BigFloat(0.5),0.,Jacobi(BigFloat(0.5),0.5));
Qb = LeftIntegral(Sb,BigFloat(0.5));
condnum=zeros(200,16)
for λ in 1:4
    op=(I+λ^2*Q)[1:800,1:800]
    for N=1:800
        condnum[N,λ]=cond((I+λ^2*Q)[1:N,1:N])
    end
end
plot(condnum,yaxis=:log)
plot(condnum[end,:],yaxis=:log)