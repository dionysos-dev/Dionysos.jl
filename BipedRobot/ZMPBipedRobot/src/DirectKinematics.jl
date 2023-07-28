@variables a α d θ q1 q2

"""

Get the Homogenous transforamtion matrix from Denavit-Hartenberg's conventioon from the i-th body to the i+1-th body in Symbolic format 
"""
function getTransformationMatrix(a_i::Vector, α_i::Vector, d_i::Vector, θ_i::Vector)
    H_DH = [
        cos(θ) -sin(θ)*cos(α) sin(θ)*sin(α) a*cos(θ)
        sin(θ) cos(θ)*cos(α) -cos(θ)*sin(α) a*sin(θ)
        0 sin(α) cos(α) d
        0 0 0 1
    ]
    H_ij = Array[]
    for i in eachindex(a_i)
        H =
            simplify.(
                substitute.(
                    H_DH,
                    (Dict(a => a_i[i], α => α_i[i], d => d_i[i], θ => θ_i[i]),),
                )
            )
        push!(H_ij, H)
    end
    H_0n = eye(4)
    for i in eachindex(H_ij)
        H_0n = H_0n * H_ij[i]
    end
    H_0n = simplify.(H_0n)
    return H_0n
end

"""

Return the location of the ankle joint w.r.t the hip frame in function format
"""
function foot2hip(br::BipedRobot)
    L1 = br.L_thigh
    L2 = br.L_leg
    d1 = br.offset_hip_to_motor

    a_i = [0, L1, L2, 0]
    α_i = [pi / 2, 0, 0, 0]
    d_i = [-d1, 0, 0, 0]
    θ_i = [0, q1 - pi / 2, q2, q1 + pi / 2 + q2]

    H_foot2hip = getTransformationMatrix(a_i, α_i, d_i, θ_i)
    p_ground_foot = [0; 0; 0; 1]

    p_ground_hip = H_foot2hip * p_ground_foot
    p = build_function.(p_ground_hip, q1, q2)
    return eval.(p)[1:3]
end

"""

Return the location of the hip joint w.r.t the ankle jo in function format 
"""
function hip2foot(br::BipedRobot)
    L1 = br.L_thigh
    L2 = br.L_leg

    a_i = [L2, L1, 0, 0]
    α_i = [0, 0, 0, 0]
    d_i = [0, 0, 0, 0]
    θ_i = [pi / 2 + q2 + q1, -q2, -pi / 2 + q1, 0]

    H_hip2foot = getTransformationMatrix(a_i, α_i, d_i, θ_i)
    p_hip_hip = [0; 0; 0; 1]

    p_hip_foot = *(H_hip2foot, p_hip_hip)
    p = build_function.(p_hip_foot, q1, q2)
    return eval.(p)[1:3]
end
