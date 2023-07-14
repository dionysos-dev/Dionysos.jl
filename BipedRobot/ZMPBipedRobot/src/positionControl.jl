struct PID 
    Kp::Real
    Ki::Real
    Kd::Real
end 

function pid_control!(pid::PID, Δq::Vector, Δq̇::Vector, ∫edt::Vector, Δt::Float64)
    output = pid.Kp .* Δq .+ pid.Ki .* ∫edt .+ pid.Kd .* Δq̇;
    ∫edt = ∫edt .+ Δq * Δt
    # pid.int_error = pid.int_error .+ Δq * Δt; 
    return output
end
