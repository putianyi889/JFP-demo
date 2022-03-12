include("testmodule.jl")
using PlotlyJS, PlotlyBase
α=0.0;β=0.0;b=-1;p=2;μ=0.5;N=200;
(α,β,b,μ)=BigFloat.((α,β,b,μ));
iopref=testmodule.OpI22_stable(α,β,b,p,μ,N);
iop21=testmodule.OpI21(α,β,b,p,μ,N);
iop22=testmodule.OpI22_unstable(α,β,b,p,μ,N);

# plot configuration
layout=Layout(
    xaxis=attr(
        title=attr(
            text="column index",
            font_size=20,
            font_family="Times New Roman",
            font_color="black",
        ),
        tickfont_family="Times New Roman",
        tickfont_size=15,
        tickfont_color="black",
        mirror=true, showline=true,
        linecolor="black",
        side="top"
    ), 
    yaxis=attr(
        title=attr(
            text="row index",
            font_size=20,
            font_family="Times New Roman",
            font_color="black"
        ),
        tickfont_family="Times New Roman",
        tickfont_size=15,
        tickfont_color="black",
        autorange="reversed",
        mirror=true,
        showline=true,
        linecolor="black"
    ), 
    plot_bgcolor="white",
    width=400,
    height=400
)

# integral by part
plt1=PlotlyJS.plot(
    PlotlyJS.contour(
        z=log10.(abs.(iopref-iop21)), 
        contours_showlabels=true, 
        contours_coloring="lines", 
        contours_labelfont_family="Times New Roman",
        connectgaps=true, 
        line_smoothing=0.5, 
        showscale=false, 
        colorscale=[[0, "black"], [1, "black"], ]
    ), 
    layout
)

# commutativity
plt2=PlotlyJS.plot(
    PlotlyJS.contour(
        z=log10.(abs.(iopref-iop22)), 
        contours_showlabels=true, 
        contours_coloring="lines", 
        contours_labelfont_family="Times New Roman",
        connectgaps=true, 
        line_smoothing=0.5, 
        showscale=false, 
        colorscale=[[0, "black"], [1, "black"], ]
    ), 
    layout
)