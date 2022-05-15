# Fractional heat/wave equation (example)
include("testmodule.jl")
using ContinuumArrays, ClassicalOrthogonalPolynomials, Plots, LinearAlgebra, PlotlyJS
pgfplotsx()

# x dim
x=Inclusion(0..2π);
F=Fourier()[x,:];
xgrid=range(0,2π,length=201);
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
tgridshift=range(0,2,length=201);
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
plt=PlotlyJS.plot(
    PlotlyJS.surface(
        z=sol,y=xgrid,x=T/2*tgridshift,
        hidesurface=true,
        showscale=false,
        contours=attr(
            x_show=true, x_color="blue", x_size=0.5, x_start=0, x_end=T,
            y_show=true, y_color="red", y_size=0.25π, y_start=0, y_end=2π
        )
    ),
    Layout(size=(700,700), font_family="Helvetica",
        margin=attr(b=0,l=0,r=0,t=0),
        scene=attr(
            aspectmode="cube",
            camera_eye=attr(autoexpand=false,x=1.5,y=-1.5,z=1.5),
            xaxis=attr(
                title="<i>t</i>",
                range=(0,8)
            ),
            yaxis=attr(
                title="<i>x</i>",
                tickmode="array",
                tickvals=π*(0:0.5:2),
                ticktext=["0" "π/2" "π" "3π/2" "2π"]
            ),
            zaxis=attr(
                title="<i>u</i>",
            )
        )
    )
)

# contour plot (Plots.pgfplotsx and PlotlyJS)
Plots.contour(Float16.(xgrid), Float16.(T/2*tgridshift), Float16.(sol), xlabel="\$x\$", ylabel="\$t\$", xticks=((0:0.5π:2π),["0" "\$\\pi/2\$" "\$\\pi\$" "\$3\\pi/2\$" "\$2\\pi\$"]), xlims=(0,2π), framestyle=:box, size=(300,270), colorbar=false)

layout=Layout(
    xaxis=attr(
        title=attr(
            text="<i>x</i>",
            font_size=15,
            font_family="Helvetica",
            font_color="black",
            standoff=5
        ),
        tickfont_family="Helvetica",
        tickfont_size=10,
        tickfont_color="black",
        tickvals=π*(0:0.5:2),
        ticktext=["0" "<i>π</i>/2" "<i>π</i>" "3<i>π</i>/2" "2<i>π</i>"],
        mirror=true, showline=true,
        linecolor="black",
    ), 
    yaxis=attr(
        title=attr(
            text="<i>t</i>",
            font_size=15,
            font_family="Helvetica",
            font_color="black",
            standoff=10
        ),
        tickfont_family="Helvetica",
        tickfont_size=10,
        tickfont_color="black",
        #autorange="reversed",
        mirror=true,
        showline=true,
        linecolor="black"
    ), 
    margin=attr(l=0,r=8,b=0,t=0),
    plot_bgcolor="white",
    width=300,
    height=300
)
plt1=PlotlyJS.plot(
    PlotlyJS.contour(
        z=sol, x=xgrid, y=T/2*tgridshift, 
        contours_showlabels=true, 
        contours_coloring="lines", 
        contours_labelfont_family="Helvetica",
        contours_labelfont_size=8,
        contours_start=-1, contours_size=0.2, contours_end=4,
        connectgaps=true, 
        line_smoothing=0.5, 
        showscale=false, 
        # colorscale=[[0, "black"], [1, "black"], ]
    ),
    layout
)

# wireframe plot
l=(minimum(sol),maximum(sol))
Plots.plot(size=(300,300), camera=(30,30), legend=false, xlims=(0,2π), ylims=(0,8), zlims=l, xlabel="\$x\$", ylabel="\$t\$", zlabel="\$u\$", xticks=((0:0.5π:2π),["0" "\$\\pi/2\$" "\$\\pi\$" "\$3\\pi/2\$" "\$2\\pi\$"]), foreground_color_grid=:gray90, gridalpha=1, framestyle=:box)
for x in 1:8:length(xgrid)
    Plots.plot!(Float16.(xgrid[x]*ones(length(tgridshift))),Float16.(T/2*tgridshift),Float16.(sol[:,x]), color=:blue, linewidth=0.3)
end
for t in 1:5:length(tgridshift)
    Plots.plot!(Float16.(xgrid), Float16.(T/2*tgridshift[t]*ones(length(xgrid))), Float16.(sol[t,:]),color=:blue, linewidth=0.3)
end
Plots.plot!([0;0], [0;0], [l[1];l[2]], color=:black)
Plots.plot!([0;0], [8;8], [l[1];l[2]], color=:black)
Plots.plot!([2π;2π], [8;8], [l[1];l[2]], color=:black)