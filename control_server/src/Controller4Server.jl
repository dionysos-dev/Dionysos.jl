module Controller4Server

export Controller

mutable struct Controller{TF, TG}
    x::Vector{Float64}
    f::TF
    g::TG
end

end #module
