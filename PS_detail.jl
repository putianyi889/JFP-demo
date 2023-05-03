include("testmodule.jl")
using Plots, LaTeXStrings
pgfplotsx()

# a
α=0.0;β=0.0;b=0;p=2;μ=0.5;N=10000;
setprecision(128);
grow128=testmodule.OpI22_precShifts(BigFloat.((α,β,b,p,μ))...,N);
setprecision(256);
grow256=testmodule.OpI22_precShifts(BigFloat.((α,β,b,p,μ))...,N);
setprecision(384);
grow384=testmodule.OpI22_precShifts(BigFloat.((α,β,b,p,μ))...,N);
data=2. .^ testmodule.finitediff([grow128 grow256 grow384],step=1);
plt=Plots.plot(data, xlabel=L"n", legend=false, xlims=(0,10000), size=(400,300))
annotate!(6000,12,(L"r(384,n)",10))
annotate!(9000,23,(L"r(256,n)",10))
annotate!(2000,18,(L"r(128,n)",10))

# b and c
α=0.0;β=0.0;b=0;p=3;μ=2//3;N=10000;
setprecision(128);
grow128=testmodule.OpI22_precShifts(BigFloat.((α,β,b,p,μ))...,N);
setprecision(256);
grow256=testmodule.OpI22_precShifts(BigFloat.((α,β,b,p,μ))...,N);
setprecision(384);
grow384=testmodule.OpI22_precShifts(BigFloat.((α,β,b,p,μ))...,N);

# b
data=2. .^ testmodule.finitediff([grow128 grow256 grow384],step=128);
plt=Plots.plot(data, ylims=(5,15), xlabel=L"n", legend=false, xlims=(0,10000), size=(400,300))
annotate!(6000,7.5,(L"q=384",10))
annotate!(9000,12,(L"q=256",10))
annotate!(2000,10,(L"q=128",10))

# c
data=2. .^ testmodule.finitediff([grow128 grow256 grow384],step=1);
plt=Plots.plot(data, xlabel=L"n", legend=:topleft, xlims=(0,10000), size=(400,300), linestyle=[:dot :dash :solid], labels=latexstring.("\$r(",[128 256 384],",n)\$"))