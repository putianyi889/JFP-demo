using Pkg
Pkg.activate("JFP-demo")

# Functions for operators, coefficients, etc.
module testmodule
using ClassicalOrthogonalPolynomials
using ClassicalOrthogonalPolynomials: (..), Inclusion, massmatrix
using KahanSummation, SparseArrays, QuasiArrays, LaTeXStrings, Plots, CircularArrays, LinearAlgebra, FillArrays, FastTransforms, LazyArrays, BandedMatrices, ThreadPools, SpecialFunctions

include("testpkg.jl")
include("testplot.jl")
include("testmisc.jl")
include("testop.jl")
include("testfiop.jl")
include("testsyl.jl")
include("testml.jl")

end
