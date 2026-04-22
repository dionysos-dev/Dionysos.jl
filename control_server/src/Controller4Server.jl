module Controller4Server

export Controller

mutable struct Controller
    x::Vector{Float64}
    f::Function
    g::Function
end

end #module
