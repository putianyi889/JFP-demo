export OpΛ, OpC, OpI, OpInvW, OpL, OpM, OpM0, OpR, OpW

OpR(γ,δ,α,β) = Jacobi(γ,δ) \ Jacobi(α,β); # R_{(α,β)}^{(γ,δ)}
OpL(α,β,k,j) = Jacobi(α,β) \ (JacobiWeight(k,j).*Jacobi(α+k,β+j)); # L_{(α+k,β+j)}^{(α,β)}
OpW(β) = Diagonal(β+1:∞); # W^{(β)}
OpInvW(β) = Diagonal(inv.(β+1:∞)); # W^{(β)}^{-1}
OpI(α,β,b,p) = p/2^(p-1)*OpL(α,β,0,p)*OpR(α,β+p,α-1,b+p)*OpInvW(b+p-1)*OpR(α,b+p-1,α,β); # I^{(α,β)}_{b,p}
function OpI(α,β,b,p,N) 
    if p>1
        L=*(map(x->x[1:N,1:N],OpL(α,β,0,p).args)...)
    else
        L=OpL(α,β,0,p)[1:N,1:N]
    end
    if b!=β
        R1=*(map(x->x[1:N,1:N],OpR(α,β+p,α-1,b+p).args)...)
    else
        R1=OpR(α,β+p,α-1,b+p)[1:N,1:N]
    end
    if b+p-1-β>1
        R2=*(map(x->x[1:N,1:N],OpR(α,b+p-1,α,β).args)...)
    else
        R2=OpR(α,b+p-1,α,β)[1:N,1:N]
    end
    return p/2^(p-1)*L*R1*BandedMatrix(OpInvW(b+p-1)[1:N,1:N])*R2
end
OpM(α,β) = massmatrix(Jacobi(α,β))
OpM(α,β,N) = Diagonal(Vector(diag(OpM(α,β)[1:N+1,1:N+1])))
OpX(α,β,b,p) = (2.0^(1-p)*(OpL(α,β,0,1)*OpR(α,β+1,α,β))^p-Eye(∞))

# C^{(α,β)}
OpCsingle(α,β,k,n) = (-1)^(n-k)*fracpochhammer(k+β+1,one(β),n-k)*fracpochhammer(n+α+β+1,2*one(β),1,2,k); # [k,n]
function OpCcolumn(α,β,n) # [:,n]
    (α,β)=promote(α,β)
    ret=zeros(typeof(α),n+1)
    ret[1]=(-1)^n*fracpochhammer(β+1,1,n);
    for k=1:n
        ret[k+1]=ret[k]*(n+α+β+k)*(n-k+1)/(-2*k*(k+β));
    end
    return ret;
end
function OpC(α,β,N)
    (α,β)=promote(α,β)
    ret=zeros(typeof(α),N+1,N+1)
    for n=0:N
        ret[1:n+1,n+1] = OpCcolumn(α,β,n)
    end
    return ret
end

OpD(b,p) = Diagonal((2.0^(1-1/p)) .^ (b.+(0:∞)))
OpD(b,p,N) = Diagonal((2.0^(1-1/p)) .^ (b.+(0:N)))
OpC(α,β,b,p,N) = OpD(b,p,N)*OpC(α,β,N)

OpΛ(r,μ) = BandedMatrix(-1=>(gamma.((r+1):μ:∞)./gamma.((r+μ+1):μ:∞)))
OpΛ(r,μ,k) = BandedMatrix(-k=>(gamma.((r+1):μ:∞)./gamma.((r+k*μ+1):μ:∞)))

OpM0(α,β) = FunMatrix(
    promote_type(typeof(α),typeof(β)),
    (m,n)->2^(α+β+m+n-1)*beta(β+m+n-1,α+1),
    ∞,∞)
OpM0(α,β,N) = OpM0(α,β)[1:N+1,1:N+1]
