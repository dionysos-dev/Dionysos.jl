module Controller4Server

export Controller

mutable struct Controller{TX, TF, TG}
    x::TX
    f::TF
    g::TG
end

end # module
