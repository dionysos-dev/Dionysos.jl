
"""

Construct the cart-table model described by :

1. Kajita S, Kanehiro F, Kaneko K, Fujiwara K, Harada K, Yokoi K, et al. 
Biped walking pattern generation by using preview control of zero-moment point. 
In: 2003 IEEE International Conference on Robotics and Automation (Cat No03CH37422) [Internet]. 
Taipei, Taiwan: IEEE; 2003 [cité 20 déc 2022]. p. 1620‑6. 
Disponible sur: http://ieeexplore.ieee.org/document/1241826/

"""
function cartTableModel(zc::Union{Int64, Float64}, g::Union{Int64, Float64})
    A = [
        0 1.0 0
        0 0 1.0
        0 0 0
    ]
    B = [
        0.0
        0.0
        1.0
    ]
    C = [1 0 -zc / g]
    return A, B, C
end

"""

Convert a continous state-space model to a discrete time state-space model 

Here, the cart-table model has a `rank(A) < 3` which lead to use the 
Taylor definition to find Bd with finite terms 
The higher order terms can be computed but they are equal to 0

"""
function continuous2discrete(
    A::T,
    B::Vector,
    C::T,
    Ts::Union{Int64, Float64},
) where {T <: Union{Vector, Matrix}}
    Ad = exp(A * Ts)
    Bd =
        B * Ts +
        1 / 2 * A * B * Ts^2 +
        1 / 6 * A^2 * B * Ts^3 +
        1 / (4 * 6) * A^3 * B * Ts^4
    Cd = C
    return Ad, Bd, Cd
end
"""

Create an Identity matrix with the dimension `dims`
"""
function eye(dims::Int64)
    return Matrix(I, dims, dims)
end

"""

Here we solve `A x = b` system given by the initials conditions for 
the cubic equation : 

    `x(t) = a_0 + a_1 t + a_2 t^2 + a_3 t^3`

Where `x(t)` is the ZMP (x or y axis), `t` is the time and `a_i` is a coefficent. 
The initial condition is given by : 

    * `x(0) = lastZMP`
    * `x(δ * Tstep) = nextZMP`
    * `x'(0) = 0`
    * `x'(δ * Tstep) = 0 `

"""
function getSplineCoeff(xStart::Float64, xEnd::Float64, yStart::Float64, yEnd::Float64)
    A = [
        1 xStart xStart^2 xStart^3
        1 xEnd xEnd^2 xEnd^3
        0 1 2*xStart 3*xStart^2
        0 1 2*xEnd 3*xEnd^2
    ]
    b = [yStart; yEnd; 0; 0]
    return A \ b
end

"""

Here we construct the cubic equation : 

    `x(t) = a_0 + a_1 t + a_2 t^2 + a_3 t^3 `

for the given coefficents and x vector 
"""
function spline(x::T, coeff::Vector) where {T <: Union{Vector, StepRangeLen, Float64}}
    sum = zeros(length(x))
    for i in 1:length(coeff)
        sum = sum .+ coeff[i] * x .^ (i - 1)
    end
    return sum
end

"""

Open a CSV file in a table format 
"""
function openCSV(filename::String)
    # Define the path to the CSV file
    csvpath() = joinpath("$(filename)")
    # Read the CSV file into a DataFrame
    data = (CSV.read(csvpath(), DataFrame; header = [2], delim = ','))
    return data
end
