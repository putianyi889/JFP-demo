using QuasiArrays
using BandedMatrices: BandedColumnMajor, bidiagonalize!

export mul_kbn

function mul_kbn(A,B)
    C=zeros(promote_type(eltype(A),eltype(B)),size(A,1),size(B,2))
    for m=1:size(A,1)
        for n=1:size(B,2)
            C[m,n]=sum_kbn(A[m,:] .* B[:,n])
        end
    end
    return C
end

Base.issubset(A::AbstractQuasiArray,B::AbstractQuasiArray) = true

function BandedMatrices._bidiagonalize!(A::AbstractMatrix{BigFloat}, M::BandedColumnMajor)
    m,n=size(A)
    if m!=n
        ErrorException("Not implemented for m=$m and n=$n.")
    end
    α = Vector{BigFloat}(undef, m)
    β = Vector{BigFloat}(undef, m-1)
    V1=zeros(BigFloat,n); V1[1]=1;
    U1=A[:,1]; α[1]=norm(U1); U1=U1/α[1];
    V2=Vector{BigFloat}(undef,n);
    U2=Vector{BigFloat}(undef,n);
    for k=2:m
        V2=A'*U1-α[k-1]*V1; β[k-1]=norm(V2); V2=V2/β[k-1];
        U2=A*V2-β[k-1]*U1; α[k]=norm(U2); U2=U2/α[k];
        U1,U2=U2,U1; V1,V2=V2,V1;
    end
    Bidiagonal(α,β,:U)
end

function BandedMatrices.bidiagonalize!(A::AbstractMatrix{BigFloat}, M::BandedColumnMajor)
    m,n=size(A)
    if m!=n
        ErrorException("Not implemented for m=$m and n=$n.")
    end
    α = Vector{BigFloat}(undef, m)
    β = Vector{BigFloat}(undef, m-1)
    V1=zeros(BigFloat,n); V1[1]=1;
    U1=A[:,1]; α[1]=norm(U1); U1=U1/α[1];
    V2=Vector{BigFloat}(undef,n);
    U2=Vector{BigFloat}(undef,n);
    for k=2:m
        V2=A'*U1-α[k-1]*V1; β[k-1]=norm(V2); V2=V2/β[k-1];
        U2=A*V2-β[k-1]*U1; α[k]=norm(U2); U2=U2/α[k];
        U1,U2=U2,U1; V1,V2=V2,V1;
    end
    Bidiagonal(α,β,:U)
end

function LinearAlgebra.cond(A::BandedMatrix{BigFloat})
    B=bidiagonalize!(A)
    B=BandedMatrix(0=>B.dv,1=>B.ev)
    v=eigvals(SymTridiagonal(B'*B))
    sqrt(v[end]/v[1])
end

# New package
export FunArray,FunVector,FunMatrix

struct FunArray{T,N} <: AbstractArray{T,N}
    fun::Function
    sz::NTuple{N,TT} where TT <: Integer
end
Base.size(A::FunArray)=A.sz
Base.getindex(A::FunArray{T,N},I::Vararg{Int,N}) where {T,N}=T(A.fun(I...))

const FunVector{T}=FunArray{T,1}
const FunMatrix{T}=FunArray{T,2}

FunArray(T::Type,fun::Function,sz...)=FunArray{T,length(sz)}(fun,sz)
FunVector(T::Type,fun::Function,sz)=FunArray{T,1}(fun,(sz,))
FunMatrix(T::Type,fun::Function,sz1,sz2)=FunArray{T,2}(fun,(sz1,sz2))
