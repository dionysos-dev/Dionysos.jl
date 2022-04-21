module TestMain
using StaticArrays,Plots


function test_periodicity()
    fig = plot(aspect_ratio = 1,legend = false)
    hyperrectangle = AB.HyperRectangle(SVector(-1.0, 0.0), SVector(2*π+1.0, 10.0))
    hx = [0.3, 0.3]
    periodic = [1,2]
    periods = [2.0*pi,10.0]
    Xdom = D.GeneralDomainList(hx;periodic=periodic,periods=periods)
    AB.add_set!(Xdom, hyperrectangle , AB.INNER)
    plot!(U.rectangle(hyperrectangle.lb,hyperrectangle.ub), opacity=.4,color=:blue)
    U.plot_domain!(Xdom)


    h2 = AB.HyperRectangle(SVector(2*pi-0.5, 2*pi+0.5), SVector(9.0, 11.0))
    plot!(U.rectangle(h2.lb,h2.ub), opacity=.4,color=:yellow)


    rectI = AB.get_pos_lims_outer(Xdom.grid, h2)
    ypos_iter = Iterators.product(AB._ranges(rectI)...)
    grid = Xdom.grid
    r = grid.h/2
    for ypos in ypos_iter
        ypos = D.set_in_period_pos(Xdom,ypos)
        center = AB.get_coord_by_pos(grid, ypos)
        plot!(U.rectangle2(center,r), opacity=.2,color=:green)
    end
    display(fig)
end
test_periodicity()
end # end module
