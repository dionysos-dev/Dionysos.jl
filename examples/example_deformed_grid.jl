module Test
using ..Dionysos
using StaticArrays, Plots
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain

function f1(x)
    return x
end
function fi1(x)
    return x
end
function f2(x)
    return SVector(x[2]+sin(x[1]),x[1])
end
function fi2(x)
    return SVector(x[2],x[1]-sin(x[2]))
end
function f3(x)
    return SVector(x[1]*cos(x[2]),x[1]*sin(x[2]))
end
function fi3(x)
    return SVector(sqrt(x[1]*x[1] + x[2]*x[2]),atan(x[2],x[1]))
end
function rotate(x,θ)
    R = @SMatrix [ cos(θ) -sin(θ) ;
                   sin(θ)  cos(θ)]
    return R*x
end
function build_f_rotation(θ; c=SVector(0.0,0.0))
    function f(x)
        return rotate(x-c,θ)+c
    end
    function fi(x)
        return rotate(x-c,-θ)+c
    end
    return f,fi
 end


 function test_with_General_Domain(f,fi)
    X = UT.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 2*π))
    obstacle = UT.HyperRectangle(SVector(10.0, 10.0), SVector(15.0, 15.0))
    hx = [3.0, 0.3]
    d = DO.RectanglularObstacles(X, [obstacle])
    dom = DO.GeneralDomainList(hx;elems=d,f=f,fi=fi,fit=true)
    fig = plot(aspect_ratio = 1,legend = false)
    Plots.plot!(dom)
    display(fig)
end

function test_with_DomainList(f,fi)
    X = UT.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 2*π))
    grid = DO.GridFree(SVector(0.0,0.0),SVector(3.0,0.3))
    Dgrid = DO.DeformedGrid(grid,f,fi)
    dom = DO.DomainList(Dgrid)
    DO.add_set!(dom, X, DO.INNER)
    fig = plot(aspect_ratio = 1,legend = false)
    Plots.plot!(dom)
    display(fig)
end

function test()
    test_with_General_Domain(f1,fi1)

    f, fi = build_f_rotation(π/3.0)
    test_with_General_Domain(f,fi)

    test_with_General_Domain(f2,fi2)

    test_with_General_Domain(f3,fi3)

    test_with_DomainList(f1,fi1)

    f, fi = build_f_rotation(π/3.0)
    test_with_DomainList(f,fi)

    test_with_DomainList(f2,fi2)

    test_with_DomainList(f3,fi3)
end

test()


end
