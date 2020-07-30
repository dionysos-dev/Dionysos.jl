function RungeKutta4(F, x, u, tstep, nsub::Int)
	h = tstep/nsub
	for i = 1:nsub
		Fx1 = F(x, u)
		xrk = x + h*Fx1/2
		Fx2 = F(xrk, u)
		xrk = x + h*Fx2/2
		Fx3 = F(xrk, u)
		xrk = x + h*Fx3
		Fx4 = F(xrk, u)
		x += h*(Fx1 + 2*Fx2 + 2*Fx3 + Fx4)/6
	end
    return x
end

function NewControlSystemRK4(tstep, F_sys, L_bound, sys_noise, meas_noise, n_sys, n_bound)
	@assert length(sys_noise) == length(meas_noise)
	sys_map = (x, u, tstep) -> RungeKutta4(F_sys, x, u, tstep, n_sys)
	F_bound = (r, u) -> L_bound(u)*r + sys_noise
	bound_map = (r, u, tstep) -> RungeKutta4(F_bound, r, u, tstep, n_bound)
	return ControlSystem(tstep, sys_noise, meas_noise, sys_map, bound_map)
end
