

function plot_arrow!(p1,p2;dims=[1,2],color=:black,markeralpha=0.0)
    p1x,p1y = p1[dims]
    p2x,p2y = p2[dims]
    plot!([p1x,p2x],[p1y,p2y], marker=:circle, markeralpha=markeralpha, arrow=(:closed, 2.0),color=color,legend=false)
end

function plot_segment!(p1, p2; linestyle=:dash, color=:black)
    plot!([p1[1], p2[1]], [p1[2], p2[2]], linestyle=linestyle, color=color)
end