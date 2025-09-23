struct Colormap
    range::Vector{Float64}
    colormap::Vector{RGB{Float64}}
    function Colormap(range, mycolor)
        v1 = minimum(range)
        v2 = maximum(range)
        if abs(v2 - v1) <= 10e-8
            v2 += 10e-2
        end
        return new(range, mycolor)
    end
end

function Colormap(range)
    return Colormap(range, Colors.colormap("Reds"))
end

function get_color(colorMap::Colormap, val::Float64)
    minV, maxV = colorMap.range
    l = length(colorMap.colormap)
    a = Int(round(((val - minV) ./ (maxV - minV)) * l))
    b = min(l, max(1, a))
    return colorMap.colormap[Int(b)]
end

@recipe function f(colorMap::Colormap)
    marker_z --> colorMap.range
    colorbar_ticks --> (colorMap.range, string.(round.(colorMap.range, digits=2)))
    colorbar --> true
    color --> palette(colorMap.colormap)
    label --> ""
    return [], []
end
