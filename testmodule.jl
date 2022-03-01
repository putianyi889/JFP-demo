# Functions for operators, coefficients, etc.
# Basis and variables are defined on Overleaf subsection 6.2 and 6.5
module testmodule
using ClassicalOrthogonalPolynomials
using ClassicalOrthogonalPolynomials: (..), Inclusion, massmatrix
using KahanSummation, BlockArrays, SparseArrays, QuasiArrays, LaTeXStrings, Plots, CircularArrays, LinearAlgebra, FillArrays, FastTransforms, LazyArrays, BandedMatrices, ThreadPools, SpecialFunctions

include("testpkg.jl")
include("testplot.jl")
include("testmisc.jl")
include("testop.jl")
include("testfiop.jl")
include("testsyl.jl")
include("testml.jl")

# Fractional inner product
function Fip(a::AbstractVector,b::AbstractVector)
    T=promote_type(eltype(a),eltype(b));
    w=beta.(T(0.5) * (axes(a,1).+axes(b,1)'), 0.5);
    return sum(w.*(a*b'));
end

# Fractional Legendre basis
function FLb(N::Integer, T::Type=Float64)
    coef=zeros(T,N,N);
    for k=1:N
        coef[k,k]=1;
        for j=1:k-1
            coef[1:j,k] -= Fip(coef[1:j,j],coef[1:k,k]) * coef[1:j,j];
        end
        coef[1:k,k] /= sqrt(Fip(coef[1:k,k],coef[1:k,k]));
    end
    return coef;
end

# Fractional integral operator on underlying basis
FIOpUb(T::Type=Float64) = Vcat(Zeros{T}(1,∞),Diagonal(gamma.(T(1):0.5:∞)./gamma.(T(1.5):0.5:∞)))

# Fractional integral operator on fractional Legendre basis
function FIOpFLb(N::Integer,T::Type=Float64)
    C = FLb(N+1,T);
    return C \ (FIOpUb(T)[1:N+1,1:N] * C[1:N,1:N]);
end

function guessCoef(α,β,n)
    k=0:n;
    return gamma((n+β)/2+1)./gamma.((k.+β)./2 .+1).*gamma.((n+α+β.+k)./2 .+1)./gamma((n+α+β)/2+1)./gamma.((n.-k)./2 .+1)./gamma.(k./2 .+1).*(-1).^(n.-k);
end

end
