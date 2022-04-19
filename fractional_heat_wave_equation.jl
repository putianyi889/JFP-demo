# Fractional heat/wave equation (example)
include("testmodule.jl")
using ContinuumArrays, ClassicalOrthogonalPolynomials, Plots, LinearAlgebra, PlotlyJS
pgfplotsx()

# x dim
x=Inclusion(0..2π);
F=Fourier()[x,:];
xgrid=0:0.005π:2π;
K=40;
xksol=F[xgrid,1:2*K+1]';

# initial condition
f=exp.(-cos.(2*x)+0.5*sin.(x))-2*sin.(sin.(x));
Plots.plot(xgrid,f[xgrid], xlabel="\$x\$", ylabel="\$f(x)\$", label=false, xticks=((0:0.5π:2π),["0" "\$\\pi/2\$" "\$\\pi\$" "\$3\\pi/2\$" "\$2\\pi\$"]), xlims=(0,2π))

# Fourier coefficients of the initial condition
fv=(F\f)[1:2*K+1];
Plots.plot(abs.(fv), ylims=(1e-18,1), yaxis=:log)

# t dim
α=0.0;β=0.0;b=0;
p=5;μ=0.2;
T=8; # time scale
N=1200; # increases with T
S=Jacobi(α,β);
taxis=Inclusion(-1..1);
iop=testmodule.OpI22(α,β,b,p,μ,N)[1:N+1,1:N+1]; # the slowest process
tf=ones(Inclusion(-1..1));
tfv=(S \ ((1 .+ Inclusion(-1..1)).^(-b) .* tf))[1:N+1];
tv=(I+(T/2)^μ*K^2*iop)\tfv;
if abs(tv[end])>1e-15
    # check if N is large enough
    @warn "Possible loss of accuracy: the last coefficient is $(tv[end])."
end
tgridshift=0:0.01:2;
tksol=zeros(length(tgridshift),2*K+1);
fill!(view(tksol,:,1),1);
for k=1:K
    tkgridshift=BigFloat(k/K)^(2/μ).*tgridshift;
    tkgrid=tkgridshift.-1;
    tkygrid=testmodule.x2y.(tkgrid,p);
    tkQ=Float64.((1 .+ tkygrid).^b .* S[tkygrid,1:N+1]);
    tksol[:,2*k]=tkQ*tv;
    tksol[:,2*k+1]=tksol[:,2*k];
end

# assemble
sol=tksol*diagm(fv)*xksol;

# interactive surface plot
PlotlyJS.plot(
    PlotlyJS.surface(
        z=sol,y=xgrid,x=T/2*tgridshift,
        hidesurface=true,
        contours=attr(
            x_show=true, x_color="blue", x_size=1, x_start=0, x_end=T,
            y_show=true, y_color="red", y_size=0.5π, y_start=0, y_end=2π
        )
    ),
    Layout(size=(700,700), font_family="Times New Roman",
        scene=attr(
            aspectmode="cube", zaxis_backgroundcolor="black", 
            xaxis_title="t", yaxis_title="x", zaxis_title="u",
            yaxis_tickmode="array",
            yaxis_tickvals=π*(0:0.5:2),
            yaxis_ticktext=["0" "π/2" "π" "3π/2" "2π"]
        )
    )
)

# contour plot
Plots.contour(xgrid, T/2*tgridshift, sol, xlabel="\$x\$", ylabel="\$t\$", xticks=((0:0.5π:2π),["0" "\$\\pi/2\$" "\$\\pi\$" "\$3\\pi/2\$" "\$2\\pi\$"]), xlims=(0,2π), framestyle=:box, size=(400,400))

# wireframe plot
Plots.wireframe(xgrid[1:10:end], T/2*tgridshift[1:10:end], sol[1:10:end,1:10:end], xlabel="\$x\$", ylabel="\$t\$", zlabel="u", xticks=((0:0.5π:2π),["0" "\$\\pi/2\$" "\$\\pi\$" "\$3\\pi/2\$" "\$2\\pi\$"]), xlims=(0,2π), framestyle=:box, size=(400,400), legend=false)
