include("testmodule.jl")
using Plots
pgfplotsx()
α=0.0;β=0.0;b=0;p=2;μ=0.5;N=500;

# compute the operator in different ways
@time iopref=testmodule.OpI22_stable(BigFloat.((α,β,b,p,μ))...,N); # accurate result
@time ioplow=testmodule.OpI22_unstable(α,β,b,p,μ,N); # double precision
@time iophigh=testmodule.OpI22_unstable(BigFloat.((α,β,b,p,μ))...,N); # 256-bit precision
iophybrid=zeros(BigFloat,N+2,N+1); # hybrid precision
iophybrid[1:3,1:2]=ioplow[1:3,1:2];
iop=testmodule.OpI(BigFloat.((α,β,b,p))...)[1:N+10,1:N+10];
@time testmodule.bandedSylvester_halfforward!(iophybrid,iop,iop);

# errors
errlow=maximum(abs.(iopref-ioplow),dims=1);
errhigh=maximum(abs.(iopref-iophigh),dims=1);
errhybrid=maximum(abs.(iopref-iophybrid),dims=1);

# error overview
plt=plot([errlow' errhigh' errhybrid'], ylims=(1e-90,1e300), yaxis=:log, size=(400,300), legend=:topleft, linestyle=[:dash :solid :dot], linecolor=[:blue :black :red], labels=["double precision" "256-bit precision" "hybrid precision"], xlabel="column index", ylabel="maximum error", yticks=10. .^(-90:45:270))

# growth rate
plt=plot((x->x[2:end,:]./x[1:end-1,:])([errlow' errhigh' errhybrid']), xlabel="column index", ylabel="growth rate", labels=["double precision" "256-bit precision" "hybrid precision"], linestyle=[:solid :dash :dot], linecolor=[:black :blue :red], legend=:topleft, size=(400,300))