function RungeKutta4!(F, x, u, tstep, nsub::Int)
	τ = tstep/nsub
	for i = 1:nsub
		Fx1 = F(x, u)
		xrk = x .+ Fx1.*(τ/2)
		Fx2 = F(xrk, u)
		xrk = x .+ Fx2.*(τ/2)
		Fx3 = F(xrk, u)
		xrk = x .+ Fx3.*τ
		Fx4 = F(xrk, u)
		x = x .+ (Fx1 .+ Fx2.*2 .+ Fx3.*2 .+ Fx4).*(τ/6)
	end
	return x
end

function NewControlSystemRK4(tstep, F_sys, L_bound, sys_noise::NTuple{N, Float64},
		meas_noise::NTuple{N, Float64}, n_sys, n_bound) where N
	sys_map = (x, u, tstep) -> RungeKutta4!(F_sys, x, u, tstep, n_sys)
	F_bound = (r, u) -> L_bound(r, u) .+ sys_noise
	bound_map = (r, u, tstep) -> RungeKutta4!(F_bound, r, u, tstep, n_bound)
	return ControlSystem(tstep, sys_noise, meas_noise, sys_map, bound_map)
end