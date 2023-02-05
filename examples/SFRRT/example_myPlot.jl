using Dionysos
UT = Dionysos.Utils
using Plots, Colors
using LinearAlgebra

function example_colormap()
    n_x = 2
    Ellispoids = [UT.Ellipsoid(Matrix{Float64}(I(n_x))*8.0, [-10.0;-10.0]),
                  UT.Ellipsoid(Matrix{Float64}(I(n_x))*5.0, [0.0;-10.0]),
                  UT.Ellipsoid(Matrix{Float64}(I(n_x))*1.0, [-10.0;0.0]),
                  UT.Ellipsoid(Matrix{Float64}(I(n_x))*3.0, [10.0;10.0])]
    vals = [4.0,1.0,10.0,6.0]
    minVal = 0.0 #min(vals...) #0.0
    maxVal = max(vals...)
    colormap = Colors.colormap("Blues")
    mycolorMap = UT.Colormap([minVal,maxVal], colormap)
    plot(aspect_ratio=:equal)
    for i=1:4
        UT.plotE!(Ellispoids[i]; color=UT.get_color(mycolorMap, vals[i]))
    end
    UT.plot_colorBar!(mycolorMap)
end

function example_arrow()
    n_x = 2
    E0 = UT.Ellipsoid(Matrix{Float64}(I(n_x))*8.0, [-10.0;-10.0])
    EF = UT.Ellipsoid(Matrix{Float64}(I(n_x))*3.0, [10.0;10.0])
    plot(aspect_ratio=:equal)
    UT.plotE!(E0)
    UT.plotE!(EF)
    UT.plot_arrow!(E0.c,EF.c;dims=[1,2],color=:black,markeralpha=0.0)
end

example_colormap()
# example_arrow()
