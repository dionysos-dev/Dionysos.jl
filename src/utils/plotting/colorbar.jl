struct Colormap
    range::Vector{Float64}
    colormap::Vector{RGB{Float64}}
end

function Colormap(range)
    return Colormap([min(range), max(range)],Colors.colormap("Blues"))
end

function get_color(colorMap::Colormap, val::Float64)
    minV,maxV = colorMap.range
    l = length(colorMap.colormap)
    a = Int(round(((val-minV)./(maxV-minV))*l))
    b = min(l,max(1,a))
    return colorMap.colormap[Int(b)]
end

function plot_colorBar!(colorMap::Colormap)
    scatter!([], [],marker_z=colorMap.range,legend=false,colorbar=true,color=palette(colorMap.colormap))
end

