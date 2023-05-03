using CircularArrays
using ThreadPools

export OpI11, OpI12, OpI13, OpI22, OpI223, OpI21

function OpI1_check(p,μ)
    if !isinteger(p*μ)
        @warn "pμ must be integer. Currently: p=$p, μ=$μ, pμ=$(p*μ)"
    end
end

function OpI2_check(β,b,p,μ)
    if !isinteger(p*μ)
        @warn "pμ must be integer. Currently: p=$p, μ=$μ, pμ=$(p*μ)"
    end
    if !isinteger(p) || p<0
        @warn "p must be non-negative integer. Currently: p=$p"
    end
    if !isinteger(β-b) || β-b<0
        @warn "β-b must be non-negative integer. Currently: β-b=$(β-b)"
    end
    if b+p-1-β<0
        @warn "b+p-1-β must be non-negative. Currently: b+p-1-β=$(b+p-1-β)"
    end
end

function OpI11(α,β,b,p,μ,N)
    #OpI1_check(p,μ)
    k=Int(round(p*μ))
    A = OpC(α,β,N+k)
    (α,β,b,p,μ) = promote(α,β,b,p,μ)
    return A\(2^(μ-k)*OpΛ(b/p,1/p,k)[1:N+k+1,1:N+1]*A[1:N+1,1:N+1])
end

function OpI11_stable(α,β,b,p,μ,N)
    k=Int(p*μ)
    (α,β,b,p,μ) = promote(α,β,b,p,μ)
    envprec=precision(BigFloat) # save environment precision
    varprec=precision(α) # save variable precision
    precup=Int(max(ceil(OpI22_precShift(Float64.((α,β,b))...,2,0.5,N)),0)+ceil(3*p)) # test and compute precision offset
    setprecision(varprec+Int(floor(1.1*precup))) # set underlying precision
    unstableret=OpI11(BigFloat.((α,β,b,p,μ))...,N) # compute operator
    setprecision(varprec)
    ret=typeof(α).(unstableret) # truncate to variable precision
    setprecision(envprec) # recover environment precision
    return ret
end

function OpI12(α,β,b,p,μ,N)
    OpI1_check(p,μ)
    k=Int(p*μ)
    A = OpC(α,β,N+k)
    (α,β,b,p,μ) = promote(α,β,b,p,μ)
    return inv(A)*(2^(μ-k)*OpΛ(b/p,1/p,k)[1:N+k+1,1:N+1]*A[1:N+1,1:N+1])
end

function OpI13(α,β,b,p,μ,N)
    OpI1_check(p,μ)
    k=Int(p*μ)
    A = OpC(α,β,N+k)
    (α,β,b,p,μ) = promote(α,β,b,p,μ)
    return 2^(μ-k) * inv(OpM(α,β,N+k)) *
            A' * (OpM0(α,β,N+k) *
            OpΛ(b/p,1/p,k)[1:N+k+1,1:N+1] *
            A[1:N+1,1:N+1])
end

function OpI13_kbn(α,β,b,p,μ,N)
    OpI1_check(p,μ)
    k=Int(p*μ)
    A = OpC(α,β,N+k)
    (α,β,b,p,μ) = promote(α,β,b,p,μ)
    return 2^(μ-k) * inv(OpM(α,β,N+k)) *
            mul_kbn(triu!(mul_kbn(A', OpM0(α,β,N+k))) ,
            OpΛ(b/p,1/p,k)[1:N+k+1,1:N+1] *
            A[1:N+1,1:N+1])
end

function OpI21(α,β,b,p,μ,N)
    OpI2_check(β,b,p,μ)
    if N<p
        return OpI11(α,β,b,p,μ,N)
    end
    k=Int(p*μ)
    (α,β,b,p,μ) = promote(α,β,b,p,μ)
    ret=zeros(typeof(α),N+k+1,N+1)
    xop=OpX(α,β,b,p)[1:N+k+1,1:N+k+1]
    xiop=xop+μ*OpI(α,β,b,p,N+k+1)
    p=Int(p)
    ret[1:p+k,1:p]=OpI11(α,β,b,p,μ,p-1);
    for n=1:N-p+1
        left=max(n-p,1)
        for m=1:n+p+k
            top=max(m-p,1); bottom=min(m+p,N+k+1)
            ret[m,n+p]=(xop[m,top:bottom]⋅ret[top:bottom,n]-ret[m,left:n+p-1]⋅xiop[left:n+p-1,n])/xiop[n+p,n]
        end
    end
    return ret;
end
OpI22(α,β,b,p,μ,N)=OpI22_stable(α,β,b,p,μ,N)

function OpI22_unstable(α,β,b,p,μ,N)
    #OpI2_check(β,b,p,μ)
    if N<p
        return OpI11(α,β,b,p,μ,N)
    end
    k=Int(round(p*μ))
    (α,β,b,p,μ) = promote(α,β,b,p,μ)
    ret=zeros(typeof(α),N+k+1,N+1)
    iop=OpI(α,β,b,p,N+k+1)
    p=Int(p)
    ret[1:p+k,1:p]=OpI11(α,β,b,p,μ,p-1);
    for n=1:N-p+1
        left=max(n-p,1)
        for m=1:n+p+k
            top=max(m-p,1); bottom=min(m+p,N+k+1)
            ret[m,n+p]=(iop[m,top:bottom]⋅ret[top:bottom,n]-ret[m,left:n+p-1]⋅iop[left:n+p-1,n])/iop[n+p,n]
        end
    end
    return ret;
end

function OpI22_stable(α,β,b,p,μ,N)
    #OpI2_check(β,b,p,μ)
    k=Int(p*μ)
    (α,β,b,p,μ) = promote(α,β,b,p,μ)
    envprec=precision(BigFloat) # save environment precision
    varprec=precision(α) # save variable precision
    n=min(N,100);
    precup=Int(max(ceil(OpI22_precShift(Float64.((α,β,b,p,μ))...,n)/n*N),0)+3*p) # test and compute precision offset
    setprecision(varprec+Int(floor(1.1*precup))) # set underlying precision
    unstableret=OpI22_unstable(BigFloat.((α,β,b,p,μ))...,N) # compute operator
    setprecision(varprec)
    ret=typeof(α).(unstableret) # truncate to variable precision
    setprecision(envprec) # recover environment precision
    return ret
end

function OpI22_precShift(α,β,b,p,μ,N)
    p=Int(p)
    k=Int(p*μ)
    iop=OpI(α,β,b,p,N+k+1)
    test=CircularArray(0.0,(N+2,2*p+1)); test[1,1]=1;
    shift=0;
    for n=1:N-p+1
        for m=1:n+p+k
            left=max(n-p,m-k,1)
            top=max(m-p,1); bottom=min(m+p,n+k)
            test[m,n+p]=(iop[m,top:bottom]⋅test[top:bottom,n]-test[m,left:n+p-1]⋅iop[left:n+p-1,n])/iop[n+p,n]
            if abs(test[m,n+p]) > 1e300
                shift += log2(abs(test[m,n+p]));
                test = CircularArray(test ./ test[m,n+p]);
            end
        end
    end
    shift += log2(maximum(abs.(test[:,N+1])))
    return shift
end

function OpI22_precShifts(α,β,b,p,μ,N;iop=0)
    p=Int(p)
    k=Int(p*μ)
    if iop==0
        @time iop=OpI(α,β,b,p,N+k+1)
    end
    (α,β,b,μ) = promote(α,β,b,μ)
    test=CircularArray(zero(α),(N+k+1,2*p+1)); test[1,1]=1;
    shift=0;
    shifts=zeros(N+1);
    for n=1:N-p+1
        print(n," ")
        Threads.@threads for m=1:n+p+k
            left=max(n-p,m-k,1)
            top=max(m-p,1); bottom=min(m+p,n+k)
            test[m,n+p]=(iop[m,top:bottom]⋅test[top:bottom,n]-test[m,left:n+p-1]⋅iop[left:n+p-1,n])/iop[n+p,n]
        end
        testmax=maximum(abs.(test[:,n+p]))
        shifts[n+p]=shift+log2(testmax)
        if testmax>1e300
            shift += log2(testmax);
            test ./= testmax;
        end
    end
    return shifts
end
