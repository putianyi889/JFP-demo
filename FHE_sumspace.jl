# Fractional heat equation (sumspace)
include("testmodule.jl")
using Plots, ApproxFun, LinearAlgebra, SpecialFunctions, GenericLinearAlgebra, ThreadPools, LaTeXStrings
using Plots.PlotMeasures
pgfplotsx()

# Float64 (7a)
N=600;
S = Jacobi(0.0,0.0) ⊕ JacobiWeight(-0.5,0.,Jacobi(-0.5,0.5));
Q = LeftIntegral(S,0.5);
y=zeros(N,5);
for n=1:5
    @time yy=((I+n^2*Q) \ 1).coefficients
    y[1:min(length(yy),N),n]=yy[1:min(length(yy),N)]
end
plot(abs.(y), size=(300,250), xlabel="coefficient index", ylabel="coefficient value (abs)", ylims=(1e-17,1e14), yaxis=:log, legend=:topright, label=false)
plot!([25],[abs(y[25,1])],markers=:xcross, label=L"\lambda=1", color=1, markersize=5)
plot!([100],[abs(y[100,2])],markers=:circle, label=L"\lambda=2", color=2, markersize=3)
plot!([190],[abs(y[190,3])],markers=:star5, label=L"\lambda=3", color=3)
plot!([300],[abs(y[300,4])],markers=:diamond, label=L"\lambda=4", color=4)
plot!([400],[abs(y[400,5])],markers=:utriangle, label=L"\lambda=5", color=5)

# BigFloat (7b)
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
plot(log10.(abs.(yb)), size=(300,250), xlabel="coefficient index", ylabel="coefficient value (abs)", legend=:topright, label=false, yticks=(yt,latexstring.("10^{",yt,"}")))
plot!([630],[log10(abs(yb[630,1]))],markers=:xcross, label=L"\lambda=1", color=1, markersize=5)
plot!([1000],[log10(abs(yb[1000,2]))],markers=:circle, label=L"\lambda=2", color=2, markersize=3)
plot!([1500],[log10(abs(yb[1500,3]))],markers=:star5, label=L"\lambda=3", color=3)
plot!([2200],[log10(abs(yb[2200,4]))],markers=:diamond, label=L"\lambda=4", color=4)
plot!([3300],[log10(abs(yb[3300,5]))],markers=:utriangle, label=L"\lambda=5", color=5)

# Sketch of solution (8)
xgrid=vcat(-0.99999,-0.99:0.01:1);
f1=Fun(S,y[:,1]).(xgrid);
f2=Fun(S,y[:,2]).(xgrid);
f3=Fun(S,y[:,3]).(xgrid);
@time fb1=Fun(Sb,yb[:,1]).(BigFloat.(xgrid));
@time fb2=Fun(Sb,yb[:,2]).(BigFloat.(xgrid));
@time fb3=Fun(Sb,yb[:,3]).(BigFloat.(xgrid));
plot(-1:0.01:1,[f1 f2 f3 fb1 fb2 fb3], xlims=(-1,1), ylims=(0,1), label=false, size=(300,250), linestyle=[:solid :solid :solid :dash :dash :dash], color=[1 2 3 1 2 3], legend=:topright, xlabel=L"x", ylabel=L"u(x)");
annotate!(0.5,0.42,text(L"\lambda=1",10));
annotate!(0.0,0.2,text(L"\lambda=2",10));
annotate!(-0.8,0.06,text(L"\lambda=3",10));
plot!([0 0],[-1 -1], color=:black, label=["numerical" "exact"], linestyle=[:solid :dash])

# Condition numbers (7e)
include("testmodule.jl")
using Plots, ApproxFun, LinearAlgebra, GenericLinearAlgebra, ThreadPools
pgfplotsx()

setprecision(2048)
Sb = Jacobi(BigFloat(0.0),1.0) ⊕ JacobiWeight(BigFloat(0.5),0.,Jacobi(BigFloat(0.5),0.5));
Qb = LeftIntegral(Sb,BigFloat(0.5));

# Compute condition numbers for every truncation size. 
# This process can take a long time (O(N^2) in time and O(N) in space for each condition number). Consider using a large step to get a brief view of the overall complexity before going for a full run.
# Consider MultiFloats.jl for faster high-precision operations, but keep in mind that the package can be neither fully reliable nor actively maintained (which was the reason we didn't use it in our paper).
condnum=zeros(800,4)
for λ in 1:4
    op=(I+λ^2*Qb)[1:800,1:800]
    @qthreads for N=1:800
        println("λ=$λ, N=$N: ")
        condnum[N,λ]=cond((I+λ^2*Qb)[1:N,1:N])
    end
end

# Or load our results directly:
using JLD
condnum=load("condnum.jld", "condnum") # download the file from https://github.com/putianyi889/JFP-demo/issues/4

# Plot the results
plot(condnum,yaxis=:log, legend=:topleft, labels=string.("\$\\lambda=",(1:4)',"\$"), xlabel="truncation size", ylabel="condition number", linestyle=[:solid :dash :dashdot :dot], size=(300,250), framestyle=:box, xlims=(0,800))
plot(log10.(log10.(condnum[end,:])), xaxis=:log, xlabel="\$\\lambda\$", ylabel="condition number", framestyle=:box, markers=true, grid=false, size=(300,250), legend=false, yticks=(0:1:3, latexstring.("\$10^{10^{",0:1:3,"}}\$"))) # cond=10^(O(λ^4))
