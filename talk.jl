plot(1:N+1, abs.([err1 err2 err3 err4]), yaxis=:log, xlabel="coefficient index", labels=["multi-order" "mul operator" "RL type" "Caputo type"], ylabel="absolute value of coefficients", legend=:bottomleft)

plot(condnumb,yaxis=:log2, yticks=2.0.^(0:2:10), size=(300,250), xlabel="truncation size", ylabel="condition number", labels=latexstring.("\\lambda=",λ), legend=:bottomright)

savefig("jfpcond.pdf")

plot(condnum,yaxis=:log, legend=:topleft, labels=latexstring.("\$\\lambda=",(1:4)',"\$"), xlabel="truncation size", ylabel="condition number", size=(300,250), framestyle=:box, xlims=(0,800))
plot(log10.(log10.(condnum[end,:])), xaxis=:log, xlabel=L"\lambda", ylabel="condition number", framestyle=:box, markers=true, grid=true, size=(300,250), legend=false, yticks=(0:1:3, latexstring.("\$10^{10^{",0:1:3,"}}\$")))


# FHE_sumspace
using SpecialFunctions, Plots
cT(λ, n, μ) = (-λ)^n / gamma(1+n*μ)
cT(λ, n, μ, x) = cT(λ, n, μ) * (1+x) ^ (n*μ)

setprecision(256)
cT1 = Float64.(cT.(BigFloat(1), (0:50)', BigFloat(0.5), xgrid))
cT2 = Float64.(cT.(BigFloat(2), (0:100)', BigFloat(0.5), xgrid))
cT3 = Float64.(cT.(BigFloat(3), (0:160)', BigFloat(0.5), xgrid))

mT1 = cT1[end, :]
mT2 = cT2[end, :]
mT3 = cT3[end, :]

plot(yaxis=:log10, ylims=(1e-17, 1e10), yticks=10.0 .^(-15:5:20), xlabel="term index", ylabel="term value (abs)", size=(300,250), legend=:topright)
plot!(abs.(mT1), label=L"\lambda=1")
plot!(abs.(mT2), label=L"\lambda=2")
plot!(abs.(mT3), label=L"\lambda=3")

sT1 = sum(cT1, dims=2)
sT2 = sum(cT2, dims=2)
sT3 = sum(cT3, dims=2)

plot(size=(300,250), ylims=(0,1), xlabel=L"x", ylabel=L"u(x)", legend=:topright)
plot!(-1:0.01:1, sT1, label=L"\lambda=3")
plot!(-1:0.01:1, sT2, label=L"\lambda=4")
plot!(-1:0.01:1, sT3, label=L"\lambda=5", ylims=:auto, legend=:topleft)
plot!(ylims=(0,1))
plot(-1:0.01:1,[sT1 sT2 sT3 fb1 fb2 fb3], xlims=(-1,1), ylims=(0,1), label=false, size=(300,250), linestyle=[:solid :solid :solid :dash :dash :dash], color=[1 2 3 1 2 3], legend=:topright, xlabel=L"x", ylabel=L"u(x)");

plot(abs.(y[:, 1:3]), yaxis=:log, size=(300,250), xlims=(0,250), ylims=(1e-17, 1e15), legend=:topright, labels=[L"\lambda=1" L"\lambda=2" L"\lambda=3"], xlabel="coefficient index", ylabel="coefficient value (abs)", yticks=[1e-15, 1e-10, 1e-5, 1, 1e5, 1e10, 1e15])

plot(abs.(yb[:, 1:3]), yaxis=:log, size=(300,250), xlims=(0,500), ylims=(1e-17, 1e75), legend=:topright, labels=[L"\lambda=1" L"\lambda=2" L"\lambda=3"], xlabel="coefficient index", ylabel="coefficient value (abs)", yticks = 10.0 .^(-15:15:75))

plot(abs.(uv), ylims=(1e-17,1), yaxis=:log, legend=:topright, size=(300,250), yticks=10.0.^(-16:4:0),xlabel="coefficient index", ylabel="coefficient value (abs)", labels=latexstring.("\\lambda=",λ))

plot(abs.(err), ylims=(1e-17,1), yaxis=:log, legend=:topright, size=(300,250), yticks=10.0.^(-16:4:0),xlabel="truncation size", ylabel="maximum error", labels=latexstring.("\\lambda=",λ))

x = 0:0.01:2
λ = [1 2 3]
u=testmodule.mittleff_matlab.(μ,-λ.^2 .*(0:0.01:2).^μ)
plot(x, u, labels=latexstring.("\\lambda=",λ), legend=:topright, xlabel=L"x", ylabel=L"u(x)", size=(400,250))