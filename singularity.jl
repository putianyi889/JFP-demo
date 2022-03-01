using ClassicalOrthogonalPolynomials, Plots, SpecialFunctions, LaTeXStrings
pgfplotsx()
N=2 .^ (1:20); U=0.1:0.1:0.5
S=Chebyshev();
xgrid=-1:0.01:1;
f=[1/gamma(1+u)*(Inclusion(-1..1).+1).^u for u in U]
err=zeros(length(N),length(U))
for n=1:length(N)
    for u=1:length(U)
        fv=S[:,1:N[n]]\f[u]
        err[n,u]=maximum(abs.(S[xgrid,1:N[n]]*fv-f[u][xgrid]))
    end
end

# Plots of fractional integrals
plot([x->f[u][x] for u in 1:length(U)], xlims=(-1,1), labels=string.("\$\\mu=",U',"\$"), linestyle=[:solid :dash :dashdot :dashdotdot :dot], legend=:bottomright, size=(250,250), xlabel="\$x\$", ylabel=L"$\mathcal{I}^\mu_{-1+}[1](x)$")

# Errors from polynomial approximation
plot(N,err, xaxis=:log, yaxis=:log, xlabel="truncation size", ylabel="error", labels=string.("\$\\mu=",U',"\$"), yticks=10.0.^(0:-1:-7), xticks=10.0.^(0:6), legend=:bottomleft, size=(250,250), linestyle=[:solid :dash :dashdot :dashdotdot :dot])