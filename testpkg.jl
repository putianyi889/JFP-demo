using BlockArrays, QuasiArrays
using BandedMatrices: BandedColumnMajor, bidiagonalize!

export mul_kbn

# Package patches
blocksizes(A) = blocklengths.(axes(A))
blockvcat(a::BlockedUnitRange{T}...) where T = blockedrange(vcat(map(blocklengths,a)...))
function cat_shape(a::BlockedUnitRange{T}...) where T
    ret=deepcopy(a[1])
    for k=2:length(a)
        l=pop!(ret.lasts)
        append!(ret.lasts,l.+a[k].lasts)
    end
    return ret
end
function blockcat_shape(a::BlockedUnitRange{T}...) where T
    ret=deepcopy(a[1])
    for k=2:length(a)
        append!(ret.lasts,ret.lasts[end].+a[k].lasts)
    end
    return ret
end

function blockcat_check(A::AbstractBlockArray...; dims)
    for a in A
        for k=1:length(size(A[1]))
            if k!=dims && axes(a)[k]!=axes(A[1])[k]
                throw(DimensionMismatch("mismatch in dimension $k (expected $(A[1].axes[k]) got $(a.axes[k]))"))
            end
        end
    end
end
function cat_shape(A::AbstractBlockArray...; dims)
    blockcat_check(A...,dims=dims)
    d=vcat([k for k in A[1].axes])
    d[dims]=cat_shape(map(a->axes(a)[dims],A)...)
    return tuple(d...)
end
function blockcat_shape(A::AbstractBlockArray...; dims)
    blockcat_check(A...,dims=dims)
    d=[k for k in A[1].axes]
    d[dims]=blockcat(map(a->axes(a)[dims],A)...)
    return tuple(d...)
end
blockvcat(A::AbstractBlockArray...)=blockcat(A...,dims=Val(1))
blockhcat(A::AbstractBlockArray...)=blockcat(A...,dims=Val(2))
blockcat(A::BlockArray...; dims) = _BlockArray(cat(map(a->a.blocks,A)...,dims=dims),blockcat_shape(A...,dims=dims))
blockcat(A::PseudoBlockArray...; dims) = PseudoBlockArray(cat(map(a->a.blocks,A)...,dims=dims),blockcat_shape(A...,dims=dims))

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
    ?? = Vector{BigFloat}(undef, m)
    ?? = Vector{BigFloat}(undef, m-1)
    V1=zeros(BigFloat,n); V1[1]=1;
    U1=A[:,1]; ??[1]=norm(U1); U1=U1/??[1];
    V2=Vector{BigFloat}(undef,n);
    U2=Vector{BigFloat}(undef,n);
    for k=2:m
        V2=A'*U1-??[k-1]*V1; ??[k-1]=norm(V2); V2=V2/??[k-1];
        U2=A*V2-??[k-1]*U1; ??[k]=norm(U2); U2=U2/??[k];
        U1,U2=U2,U1; V1,V2=V2,V1;
    end
    Bidiagonal(??,??,:U)
end

function BandedMatrices.bidiagonalize!(A::AbstractMatrix{BigFloat}, M::BandedColumnMajor)
    m,n=size(A)
    if m!=n
        ErrorException("Not implemented for m=$m and n=$n.")
    end
    ?? = Vector{BigFloat}(undef, m)
    ?? = Vector{BigFloat}(undef, m-1)
    V1=zeros(BigFloat,n); V1[1]=1;
    U1=A[:,1]; ??[1]=norm(U1); U1=U1/??[1];
    V2=Vector{BigFloat}(undef,n);
    U2=Vector{BigFloat}(undef,n);
    for k=2:m
        V2=A'*U1-??[k-1]*V1; ??[k-1]=norm(V2); V2=V2/??[k-1];
        U2=A*V2-??[k-1]*U1; ??[k]=norm(U2); U2=U2/??[k];
        U1,U2=U2,U1; V1,V2=V2,V1;
    end
    Bidiagonal(??,??,:U)
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
