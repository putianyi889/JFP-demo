export bandedSylvester_halfforward!, bandedSylvester_colBack!, bandedSylvester_doublesweep, bandedSylvester_systemSolve!

function bandedSylvester_init(A,B,C)
    if !isbanded(B) || !isbanded(C)
        throw(ArgumentError("B and C must be banded."))
    end
    p=max(bandwidths(B)...,bandwidths(C)...)
    N=size(A,2)-1
    k=size(A,1)-size(A,2)
    return (p,k,N)
end

function bandedSylvester_halfforward!(A,B,C;start=1)
    (p,k,N)=bandedSylvester_init(A,B,C)
    for n=start:N-p+1
        Threads.@threads for m=1:n+p+k
            left=max(n-p,m-k,1)
            top=max(m-p,1); bottom=min(m+p,n+k)
            A[m,n+p]=(C[m,top:bottom]⋅A[top:bottom,n]-A[m,left:n+p-1]⋅B[left:n+p-1,n])/B[n+p,n]
        end
    end
end

function bandedSylvester_rowForward!(A,B,C)
    (p,k,N)=bandedSylvester_init(A,B,C)
    for m=1:N+k-p+1
        for n=1:N+1
            left=max(n-p,1)
            right=min(n+p,N+1)
            top=max(m-p,1)
            A[m+p,n]=-(C[m,top:m+p-1]⋅A[top:m+p-1,n]-A[m,left:right]⋅B[left:right,n])/C[m,m+p]
        end
    end
end

function bandedSylvester_halfbackward!(A,B,C)
    (p,k,N)=bandedSylvester_init(A,B,C)
    for n=N-p+1:-1:p+1
        for m=1:n-p+k
            top=max(m-p,1); bottom=min(m+p,n+k)
            A[m,n-p]=(C[m,top:bottom]⋅A[top:bottom,n]-A[m,n-p+1:n+p]⋅B[n-p+1:n+p,n])/B[n-p,n]
        end
    end
end

function bandedSylvester_colBack!(A,B,C)
    (p,k,N)=bandedSylvester_init(A,B,C)
    bandedSylvester_halfforward!(A,B,C)
    guesscount=2*p*(N+k+1)
    guesses=zeros(eltype(A),size(A)...,guesscount)
    guessfinal=Matrix{eltype(A)}(undef,guesscount,guesscount)
    for m=1:N+k+1
        for n=N+2-2*p:N+1
            count=m*2*p+n-N-1
            guesses[m,n,count]=1
            bandedSylvester_backward!(view(guesses,:,:,count),B,C)
            guessfinal[:,count]=reshape(guesses[:,1:2*p,count],:,1)
        end
    end
    guesses=reshape(guesses,:,guesscount)
    guesscoef=guessfinal\reshape(A[:,1:2*p],:,1)
    return reshape(guesses*guesscoef,size(A)...)
end

function bandedSylvester_fullforward!(A,B,C)
    (p,k,N)=bandedSylvester_init(A,B,C)
    for n=1:N+1-p
        left=max(1,n-p)
        for m=1:N+k+1
            top=max(m-p,1); bottom=min(m+p,N+k)
            A[m,n+p]=(C[m,top:bottom]⋅A[top:bottom,n]-A[m,left:n+p-1]⋅B[left:n+p-1,n])/B[n+p,n]
        end
    end
end

function bandedSylvester_backward!(A,B,C)
    (p,k,N)=bandedSylvester_init(A,B,C)
    for n=N+1-p:-1:p+1
        for m=1:N+k+1
            top=max(m-p,1); bottom=min(m+p,N+k)
            A[m,n-p]=(C[m,top:bottom]⋅A[top:bottom,n]-A[m,n-p+1:n+p]⋅B[n-p+1:n+p,n])/B[n-p,n]
        end
    end
end

function bandedSylvester_doublesweep(A,B,C;tol=1e-15)
    r=copy(A); AA=zero(A)
    for k=1:100
        bandedSylvester_halfforward!(r,B,C)
        AA=AA+r
        bandedSylvester_backward!(AA,B,C)
        r=(A-AA)
        println("Round $k: error=$(maximum(abs.(r[:,1])))")
        if maximum(abs.(r[:,1]))<=tol
            break
        end
    end
    return AA;
end

function bandedSylvester_systemCol(A,B,C) # AB = CA
    (p,k,N)=bandedSylvester_init(A,B,C)
    l=rec2vec(size(A)...,k)
    S=spzeros(0,l)
    for n=1:N+1-p
        for m=1:n+k+p
            newrow=spzeros(1,l)
            for mm=max(1,m-p):min(m+p,n+k)
                newrow[rec2vec(mm,n,k)] -= C[m,mm]
            end
            for nn=max(1,m-k,n-p):n+p
                newrow[rec2vec(m,nn,k)] += B[nn,n]
            end
            S = vcat(S,newrow)
        end
    end
    return S
end

function bandedSylvester_systemSolve!(A,B,C)
    (p,k,N)=bandedSylvester_init(A,B,C)
    S=bandedSylvester_systemCol(A,B,C)
    nz=findnz(A)
    l=length(nz[1])
    b=vcat(zeros(eltype(S),size(S,1)),nz[3])
    SS=spzeros(eltype(S),length(nz[1]),size(S,2))
    for (ll,m,n) in zip(1:l,nz[1],nz[2])
        SS[ll,rec2vec(m,n,k)]=one(eltype(S))
    end
    S=vcat(S,SS)
    AA=S\b
    for n=1:N+1
        A[1:n+k,n]=AA[rec2vec(1,n,k):rec2vec(1,n,k)+n+k-1]
    end
end
