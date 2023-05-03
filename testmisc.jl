export fracpochhammer, logheatmap, PEval, QEval, x2y, y2x

# Number functions
fracpochhammer(a,b,n) = prod(range(a,length=n)./range(b,length=n));
fracpochhammer(a,b,stepa,stepb,n) = prod(range(a,step=stepa,length=n)./range(b,step=stepb,length=n));

# Change variables
function x2y(x,p)
    (x,p)=promote(x,p);
    return (2^(p-1)*(1+x))^(1/p) - 1;
end
function y2x(y,p)
    (y,p)=promote(y,p);
    return (1+y)^p * 2^(1-p) - 1;
end

function x2y(p)
    x=Inclusion(-one(p)..one(p))
    return 2^((p-one(p))/p) * (one(p) .+x).^(inv(p)) .- one(p);
end
function y2x(p)
    y=Inclusion(-one(p)..one(p))
    return (one(p) .+y).^p .* 2^(one(p)-p) .- one(p);
end

# for argument k, see document for "triu"
function triu2vec(M,k::Integer)
    s=size(M)
    if length(s) != 2
        @error "Not implemented for $s-dimensional array!"
    end
    if k>0
        return triu2vec(M[:,1+k:end],0)
    end
    if s[2]>s[1]+k
        v=Vector{eltype(M)}(undef,s[1]*s[2]+((-k-s[1]+1)*(s[1]+k))÷2)
        pivot=1
        for n=1:s[1]+k
            v[pivot:pivot+n-k-1]=M[1:n-k,n]
            pivot += n-k
        end
        for n=s[1]+k+1:s[2]
            v[pivot:pivot+s[1]-1]=M[:,n]
            pivot += s[1]
        end
    else
        v=Vector{eltype(M)}(undef,(s[2]*(s[2]-2*k+1))÷2)
        pivot=1
        for n=1:s[2]
            v[pivot:pivot+n-k-1]=M[1:n-k,n]
            pivot += n-k
        end
    end
    return v
end

# conversion between two indexing rules
# diagonal rule:
    # d is the index of diagonal (main->0, upper->positive, lower->negative)
    # ind is the index within that diagonal, starting from 1
# rectangle rule:
    # r is the row index and c is the column index
@inline diag2rec(d,ind) = d>0 ? (ind,ind+d) : (ind-d,ind)
@inline rec2diag(r,c) = (c-r,min(c,r))
@inline rec2vec(r,c,k) = (2*k+c) * (c-1) ÷ 2 + r

# Polynomials
PEval(α,β,f,y) = Jacobi(α,β)[y,1:size(f,2)]*f;

function QEval(α,β,b,p,f,x)
    y=x2y.(x,p)
    return (1 .+ y).^b .* PEval(α,β,f,y)
end

finitediff(v::AbstractVector;step=1)=(v[step+1:end]-v[1:end-step])./step
finitediff(v::AbstractMatrix;step=1)=(v[step+1:end,:]-v[1:end-step,:])./step
finitequot(v::AbstractVector)=v[2:end]./v[1:end-1]
finitequot(v::AbstractMatrix)=v[2:end,:]./v[1:end-1,:]