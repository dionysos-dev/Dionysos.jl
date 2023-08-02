mutable struct BipedRobot
    # Global parameters
    Ts::Union{Int64, Float64}
    zc::Union{Int64, Float64}
    g::Union{Int64, Float64}
    Ad::Matrix
    Bd::Vector
    Cd::Matrix

    # Robot description 
    zmax::Union{Int64, Float64}
    Δz::Union{Int64, Float64}
    L_thigh::Union{Int64, Float64}
    L_leg::Union{Int64, Float64}
    offset_hip_to_motor::Union{Int64, Float64}
    offset_ankle_to_foot::Union{Int64, Float64}

    # Preview controller parameters
    Qx::Matrix
    Qe::Matrix
    R::Matrix
    previewTime::Union{Int64, Float64}

    # Foot pattern parameters
    Lmax::Union{Int64, Float64}
    θ_max::Union{Int64, Float64}
    d::Float64
    initial_position::Vector
    isLeftSupport::Bool
    xPath::Union{Vector, StepRangeLen}
    yPath::Union{Vector, StepRangeLen}

    # ZMP trajectory generator parameters
    Tstep::Union{Int64, Float64}
    Tdelay::Union{Int64, Float64}
    Twait::Union{Int64, Float64}
    δ::Float64

    # CoM trajectory generator parameters
    xinit::Vector
    yinit::Vector
    ΔCoMz::Float64

    # Swing Foot trajectory parameters
    hstep::Float64
    Tver::Union{Int64, Float64}
end

function BipedRobot(;
    readFile::Bool = true,
    URDFfileName::String = "ZMP_2DBipedRobot.urdf",
    paramFileName::String = "param.jl",
)
    if (readFile)
        urdfpath() = joinpath(packagepath(), URDFfileName)
        include(joinpath(packagepath(), paramFileName))

        A, B, C = cartTableModel(zc, g)
        Ad, Bd, Cd = continuous2discrete(A, B, C, Ts)
        p = rank(Cd)       # rank of Cd
        n = rank(Ad)       # rank of Ad

        # Define cost function weights
        Qe = q_e * eye(p)    # tracking error weight matrix
        Qx = 0.0 * eye(n)     # incremental state weight matrix
        R = r_u * eye(p)    # control vector weight matrix

        L1 = 0.0
        L2 = 0.0
        d = 0.0
        offset_hip_to_motor = 0.0
        offset_ankle_to_foot = 0.0

        visual = URDFVisuals(urdfpath())
        e = root(visual.xdoc)

        for (i, joint) in enumerate(e["joint"])
            if (attribute(joint, "name")) == "l_hip_to_motor"
                s = attribute(joint["origin"]..., "xyz")
                numbers = split(s, " ")
                offset_hip_to_motor = abs(parse(Float64, numbers[3]))
                d = 2 * abs(parse(Float64, numbers[2]))
            end
            if (attribute(joint, "name")) == "l_thigh_link_to_motor"
                s = attribute(joint["origin"]..., "xyz")
                numbers = split(s, " ")
                L1 = abs(parse(Float64, numbers[3]))
            end
            if (attribute(joint, "name")) == "l_ankle"
                s = attribute(joint["origin"]..., "xyz")
                numbers = split(s, " ")
                L2 = abs(parse(Float64, numbers[3]))
            end
        end

        for (i, link) in enumerate(e["link"])
            if (attribute(link, "name") == "l_foot_link")
                s = attribute(link["visual"][1]["geometry"][1]["box"]..., "size")
                numbers = split(s, " ")
                offset_ankle_to_foot = abs(parse(Float64, numbers[3]))
            end
        end
        return BipedRobot(
            Ts,
            zc,
            g,
            Ad,
            Bd,
            Cd,
            zmax,
            Δz,
            L1,
            L2,
            offset_hip_to_motor,
            offset_ankle_to_foot,
            Qx,
            Qe,
            R,
            previewTime,
            Lmax,
            θ_max,
            d,
            initial_position,
            isLeftSupport,
            xPath,
            yPath,
            Tstep,
            Tdelay,
            Twait,
            δ,
            xinit,
            yinit,
            ΔCoMz,
            hstep,
            Tver,
        )
    end
end
