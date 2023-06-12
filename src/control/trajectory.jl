
# Trajectory closed loop
function get_closed_loop_trajectory(contsys, controller, x0, nstep; stopping=(x)->false)
    x = x0
    x_traj = [x]
    u_traj = []
    i = 0
    while !stopping(x) && i≤nstep
        u = controller(x)
        if u === nothing
            break
        end
        x = contsys.sys_map(x, u, contsys.tstep)
        push!(x_traj, x)
        push!(u_traj, u)   
        i = i+1
    end
    return x_traj, u_traj
end

function get_closed_loop_trajectory(f_eval, c_eval, cost_eval, x0, nstep; stopping=(x)->false)
    x = x0
    x_traj = [x]
    u_traj = []
    cost_traj = []
    i = 0
    while !stopping(x) && i≤nstep
          u = c_eval(x)
          if u === nothing
                break
          end
          push!(u_traj, u) 
          push!(cost_traj, cost_eval(x, u)) 
          x = f_eval(x, u)
          push!(x_traj, x)
          i = i+1
    end
    return x_traj, u_traj, cost_traj
end