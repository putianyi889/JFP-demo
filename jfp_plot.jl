include("testmodule.jl")
using Plots, ClassicalOrthogonalPolynomials, LaTeXStrings
pgfplotsx()
S=Jacobi(0,0)
f=[x->S[x,n] for n in 2:5]
g=[x->S[testmodule.x2y(x,2.0),n] for n in 2:5]

# connection between x and y
plot(x->sqrt(2*(1+x))-1, xlims=(-1,1), ylims=(-1,1), size=(250,250), legend=false, framestyle=:box, xlabel="\$x\$", ylabel="\$y\$")

# Legendre polynomials
plot(f, xlims=(-1,1), legend=false, framestyle=:box, size=(250,250), xlabel=L"$y$", linestyle=[:solid :dash :solid :dash], ylabel=L"$P_n(y)$", ylims=(-1,1))

# Legendre fractional polynomials
plot(g, xlims=(-1,1), ylims=(-1,1), legend=false, framestyle=:box, size=(250,250), xlabel="\$x\$", ylabel=L"$Q_n^{(0,0,0,2)}(x)$", linestyle=[:solid :dash :solid :dash])