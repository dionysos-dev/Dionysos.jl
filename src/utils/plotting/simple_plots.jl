

function plot_arrow!(p1, p2; dims=[1,2], color=:black, markeralpha=0.0)
    p1x,p1y = p1[dims]
    p2x,p2y = p2[dims]
    plot!([p1x,p2x],[p1y,p2y], marker=:circle, markeralpha=markeralpha, arrow=(:closed, 2.0),color=color,legend=false)
end

function plot_segment!(p1, p2; linestyle=:dash, color=:black)
    plot!([p1[1], p2[1]], [p1[2], p2[2]], linestyle=linestyle, color=color)
end

function plot_point!(p; dims=[1,2], color=:black, markeralpha=0.0, marker=:circle)
    pp = p[dims]
    scatter!([pp[1]], [pp[2]], color=color, marker=marker)
end

function plot_traj!(trajCoord; color=:black)
    plot_point!(trajCoord[1])
    for i in 1:length(trajCoord)-1
        plot_point!(trajCoord[i+1], color=color)
        plot_arrow!(trajCoord[i], trajCoord[i+1], color=color)
    end
end