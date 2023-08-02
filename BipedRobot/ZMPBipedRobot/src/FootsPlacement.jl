"""

A Foot Planner as defined in sec. 2 of : 
    1.  Khusainov R, Sagitov A, Klimchik A, Magid E. Arbitrary Trajectory Foot Planner for Bipedal Walking: 
        In: Proceedings of the 14th International Conference on Informatics in Control, Automation and Robotics [Internet]. 
        Madrid, Spain: SCITEPRESS - Science and Technology Publications; 2017 [cité 20 déc 2022]. p. 417‑24. 
        Disponible sur: http://www.scitepress.org/DigitalLibrary/Link.aspx?doi=10.5220/0006442504170424

A `FootPlanner` has 3 differents outputs :
    1.  `right` which describe the global position of the right foot next to the path 
    2.  `left` which describe the global position of the left foot next to the path 
    3.  `center` which describe the global position of the reference point along the path 

The input of of the `FootPlanner` require the robot model which is describe by `BipedRobot` : 
    1.  `Lmax` is the step length 
    2.  `θ_max` is the maximal orientation which is allowed by the robot configuration 
    3.  `offset` is the centre-to-centre distance between both feet
    4.  `xPath` and `yPath` are the coordonates of the path. It has to be sampled at least with 0.01 m of precision
    5.  `isLeftSupport` is a boolean variable which describe the initial step of the motion 

"""
struct FootPlanner
    right::Array
    left::Array
    center::Array
end

"""

Constructor 
"""
function FootPlanner(; br::BipedRobot, check::Bool = false)
    right, left, center = computeFootsPlacement(br)

    if (check)
        xpath = br.xPath
        ypath = br.yPath
        right_plot = reduce(hcat, right)
        left_plot = reduce(hcat, left)

        plt = plot(; title = "Foot Planner", xlabel = "X[m]", ylabel = "Y[m]", dpi = 600)

        plot!(xpath, ypath; label = "Reference path")
        scatter!(left_plot[1, :], left_plot[2, :]; shape = :rect, label = "Left")
        scatter!(right_plot[1, :], right_plot[2, :]; shape = :rect, label = "Right")
        scatter!(
            getindex.(center, 1),
            getindex.(center, 2);
            label = "Center",
            markershape = :cross,
        )
        display(plt)
    end
    return FootPlanner(right, left, center)
end

"""

We compute the foot placement based on the path to follow. 
The result give us an array which each element coorespond to a global coordinate of the point 

This algorithm is describe in the paper mentionned earlier with some modification.

"""
function computeFootsPlacement(br::BipedRobot)
    Lmax = br.Lmax
    θ_max = br.θ_max
    offset = br.d
    initial_position = br.initial_position
    xPath = br.xPath
    yPath = br.yPath
    isLeftSupport = br.isLeftSupport

    k = 1
    x_p = initial_position[1]
    y_p = initial_position[2]
    θ = initial_position[3]

    # Position of the left foot 
    left = Array{Float64}[]
    temp = footPosition(x_p, y_p, θ, offset, true)
    push!(left, temp)

    # Position of the right foot
    right = Array{Float64}[]
    temp = footPosition(x_p, y_p, θ, offset, false)
    push!(right, temp)

    # Position of the reference point on the path
    center = Array{Float64}[]
    push!(center, [x_p; y_p; θ])

    # Initiate the foot to move                         
    left_flag = ~isLeftSupport

    while (k < length(xPath))
        result = nextPoint(k, x_p, y_p, θ, Lmax, θ_max, xPath, yPath)
        push!(center, result[2:4])  # New reference point 
        # Store the computed vector in the right place 
        if (left_flag == true)
            push!(left, footPosition(result[2], result[3], result[4], offset, left_flag))
            push!(right, right[end]) # Support foot -> no movement 
        else
            push!(right, footPosition(result[2], result[3], result[4], offset, left_flag))
            push!(left, left[end]) # Support foot -> no movement   
        end
        left_flag = ~left_flag         # Invert support foot 
        k = Int(result[1])             # Update path index 
        x_p = result[2]                # X component of the next point
        y_p = result[3]                # Y component of the next point
        θ = result[4]                  # Angle component of the next point
    end
    return right, left, center
end

"""

Compute the right and left foot position according to the reference point. 
"""
function footPosition(
    x_p::Float64,
    y_p::Float64,
    θ_p::Float64,
    offset::Float64,
    isleft::Bool,
)
    if (isleft == true)
        x = x_p + offset / 2 * -sin(θ_p)
        y = y_p + offset / 2 * cos(θ_p)
    else
        x = x_p - offset / 2 * -sin(θ_p)
        y = y_p - offset / 2 * cos(θ_p)
    end
    return [x; y; θ_p]
end

"""

This algorithm is explained in the first paper. 
I have ajudted the algorithm to take into account the orientation of the robot with 
the variable `direction` but overall, the algorithm is still the same than the paper 

# """
function nextPoint(
    k_p::Int64,
    x_p::Float64,
    y_p::Float64,
    θ_p::Float64,
    Lmax::Float64,
    θ_max::Float64,
    xPath::T,
    yPath::T,
) where {T <: Union{Vector, StepRangeLen}}
    L0 = Lmax
    k = k_p
    result = []
    e = 1e-12
    direction = 1
    while (L0 >= 0 - e)
        if (L0 <= 0 + e && L0 >= 0 - e)
            θ_new = θ_p + direction * θ_max
            if (θ_new >= pi)
                θ_new = -pi + (θ_new - pi)
            elseif (θ_new <= -pi)
                θ_new = pi + (θ_new + pi)
            end
            return result = [k_p; x_p; y_p; θ_new]
        else
            L = sqrt((xPath[k] - x_p)^2 + (yPath[k] - y_p)^2)
            while (k < length(xPath) && L < L0)
                L = sqrt((xPath[k] - x_p)^2 + (yPath[k] - y_p)^2)
                if (L < L0)
                    k = k + 1
                end
            end
            θ = atan(yPath[k] - yPath[k - 1], xPath[k] - xPath[k - 1])
            Δ = (θ - θ_p)
            if (Δ > pi)
                Δ = 2 * pi - Δ
            elseif (Δ < -pi)
                Δ = 2 * pi + Δ
            end
            if (Δ > 0)
                direction = 1
            elseif (Δ < 0)
                direction = -1
            end
            if (abs(Δ) <= θ_max)
                return result = [k; xPath[k]; yPath[k]; θ]
            else
                L0 = L0 - 0.01
                k = k_p
            end
        end
    end
end
