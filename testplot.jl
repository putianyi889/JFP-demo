using LaTeXStrings
using Plots

export              signheatmap,    logheatmap
export colgif,      signcolgif,     logcolgif
export rowgif,      signrowgif,     logrowgif
export triucolgif,  signtriucolgif, logtriucolgif
export triurowgif,  signtriurowgif, logtriurowgif

colgif(A;B...)=@gif for n=axes(A,2)
    plot(axes(A,1),A[:,n];B...)
end
rowgif(A;B...)=@gif for n=axes(A,1)
    plot(axes(A,2),A[n,:];B...)
end
function triucolgif(A;B...)
    k=size(A,1)-size(A,2)
    @gif for n=axes(A,2)
        plot(1:n+k,A[1:n+k,n];B...)
    end
end
function triurowgif(A;B...)
    k=size(A,1)-size(A,2)
    @gif for n=axes(A,1)
        left=max(n-k,1)
        plot(left:size(A,2),A[n,left:end];B...)
    end
end

for plotmethod in (:heatmap,:contour,:colgif,:rowgif,:triucolgif,:triurowgif)
    @eval $(Symbol(string("sign",plotmethod)))(A;B...)=$plotmethod(sign.(A);B...)
    @eval $(Symbol(string("log",plotmethod)))(A;B...)=$plotmethod(log10.(abs.(A));B...)
end

logticks(t)=(t,string.("10^{",t,"}"))
logticklabel(t)=(10. .^t,string.("\$10^{",t,"}\$"))

@recipe function f(::Type{Val{:markervery}}, x, y, z; markerrange=0)
    sx, sy = x[markerrange], y[markerrange]
    # add an empty series with the correct type for legend markers
    @series begin
        seriestype := :path
        markershape --> :auto
        x := [Inf]
        y := [Inf]
    end
    # add a series for the line
    @series begin
        primary := false # no legend entry
        markershape := :none # ensure no markers
        seriestype := :path
        seriescolor := get(plotattributes, :seriescolor, :auto)
        x := x
        y := y
    end
    # return  a series for the sampled markers
    primary := false
    seriestype := :scatter
    markershape --> :auto
    x := sx
    y := sy
end