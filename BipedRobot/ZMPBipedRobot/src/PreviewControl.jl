"""

A Preview controller as defined in : 
    1. Katayama T, Ohki T, Inoue T, Kato T. Design of an optimal controller for a discrete-time system subject to previewable demand. International Journal of Control. mars 1985;41(3):677‑99. 
    2. Takaba K. A Tutorial on Preview Control Systems. 

A `PreviewController` has 3 differents actions :
    1.  `Gx` which describe the state-feedback action 
    2.  `Gi` which describe the interal action 
    3.  `Gd` which describe the preview action

The implemented controller only support a state-space representation 
according theses equations : 

    * `x(k+1) = Ad x(k) + Bd u(k) `
    * `y(k) = Cd x(k)`

Where Ad is a square matrix, Bd is a vector and Cd is a vector which 
mean that the code only support SISO (Single Input Single Output) 
system  

"""
struct PreviewController{T <: Union{Vector, Matrix}}
    Gx::Matrix
    Gi::Vector
    Gd::T
end

"""

Constructor 
"""
function PreviewController(; br::BipedRobot, check::Bool = false)
    Gx, Gi, Gd = computeGains(br)
    if (check)
        Ts = br.Ts
        plt = plot(;
            title = "Preview Controller",
            xlabel = "Preview time [s]",
            ylabel = "Preview Gain -Gd",
            dpi = 600,
        )
        plot!(0:Ts:((length(Gd) - 1) * Ts), -reduce(vcat, Gd); label = false)
        display(plt)
    end

    return PreviewController(Matrix(Gx), vec(Gi), reduce(vcat, Gd))
end

"""

Compute the gains `Gx`, `Gi` and `Gd` based on the paper 1 in sec. 3. 
"""
function computeGains(br::BipedRobot)
    Ad = br.Ad
    Bd = br.Bd
    Cd = br.Cd
    Qx = br.Qx
    Qe = br.Qe
    R = br.R
    Ts = br.Ts
    previewTime = br.previewTime

    # Preview step 
    h = Int(round(previewTime / Ts) + 1)

    # tolerence on finding the Ricatti equation 
    ϵ = 1e-10

    # Define each dimension
    p = rank(Cd)       # rank of Cd
    r = size(Bd)[1]    # rank of Bd, can't do rank([n][1] vector) 
    n = rank(Ad)       # rank of Ad

    # Define the augmentated state (see Tohru Katayama's article) 
    Ip = eye(p)
    Ba = [
        Cd * Bd
        Bd
    ]
    Ia = [
        Ip
        0
        0
        0
    ]
    Fa = [
        Cd * Ad
        Ad
    ]
    Aa = [Ia Fa]

    Qa = [
        Qe zeros(p, n)
        zeros(n, p) Qx
    ]

    # Ricatti's equation solver with fixed point method 
    Ka = Matrix{Float64}[zeros(p + n, p + n)]
    temp = eye(p + n)
    push!(Ka, temp)
    while (norm(Ka[2] - Ka[1], 2) > ϵ)
        Ka[1] = Ka[2]
        Ka[2] =
            transpose(Aa) * Ka[1] * Aa -
            transpose(Aa) *
            Ka[1] *
            Ba *
            (R[] + transpose(Ba) * Ka[1] * Ba)^-1 *
            transpose(Ba) *
            Ka[1] *
            Aa + Qa
    end
    # Gains of the controller based on the weights matrix Qe, Qx and R 
    Gi = (R[] + transpose(Ba) * Ka[2] * Ba)^-1 * transpose(Ba) * Ka[2] * Ia
    Gx = (R[] + transpose(Ba) * Ka[2] * Ba)^-1 * transpose(Ba) * Ka[2] * Fa

    Ac = Aa - Ba * (R[] + transpose(Ba) * Ka[2] * Ba)^-1 * transpose(Ba) * Ka[2] * Aa

    Xa = Matrix{Float64}[]
    temp = -transpose(Ac) * Ka[2] * Ia
    push!(Xa, temp)

    Gd = Array{Float64}[]
    push!(Gd, -Gi)
    for l in 2:h
        temp = transpose(Ac) * Xa[l - 1]
        push!(Xa, temp)
        temp = (R[] + transpose(Ba) * Ka[2] * Ba)^-1 * (transpose(Ba) * Xa[l - 1])
        push!(Gd, temp)
    end
    return Gx, Gi, Gd
end
