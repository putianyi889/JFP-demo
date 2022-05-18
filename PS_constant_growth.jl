include("testmodule.jl")
using Plots, LaTeXStrings
pgfplotsx()
default(size=(300,250), foreground_color_grid=:gray90, gridalpha=1)
α=0;β=0;b=0;p=2;μ=0.5;N=1000;
@time iopref=testmodule.OpI22_stable(α,β,b,p,μ,N);
ioplow=zeros(BigFloat,N+2,N+1);
ioplow[1:3,1:2]=iopref[1:3,1:2];
iop=testmodule.OpI(α,β,b,p)[1:N+10,1:N+10];
@time testmodule.bandedSylvester_halfforward!(ioplow,iop,iop);
err=maximum(abs.(iopref-ioplow), dims=1);
errtest=zeros(BigFloat,N+1,5)
for k=1:5
    ioptest=zeros(BigFloat,N+2,N+1);
    ioptest[(k<3 ? (k,1) : (k-2,2))...]=1;
    testmodule.bandedSylvester_halfforward!(ioptest,iop,iop);
    errtest[:,k]=maximum(abs.(ioptest), dims=1)';
end
plt=plot((err'./errtest)[3:end,2:5], yaxis=:log, legend=:topright, linestyle=[:dash :solid :dot :dashdot], xlims=(0,1000), xlabel="\$n\$", labels=["[1,1] or [2,2]=1" "[2,1]=1" "[1,2]=1" "[3,2]=1"], legendfonthalign=:right, ylims=(4e-18,4e-15))
plt=plot(testmodule.finitequot(errtest))
plot([plot((err'./errtest)[3:end,k], yaxis=:log) for k in 1:5]...)
plot((err'./errtest)[3:end,:], yaxis=:log)