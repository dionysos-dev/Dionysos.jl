# function NewSymbolicModelHash(X_grid, U_grid, cont_sys)
# 	trans_map = NewTransitionMapHash(X_grid, U_grid, X_grid)
# 	infos = Dict("transitions set" => false)
# 	return SymbolicModel(X_grid, U_grid, cont_sys, trans_map, infos)
# end

function set_transitions_from_controlsystem!(trans_map, cont_sys)
	println("set_transitions_from_controlsystem! sarted")
	remove_all_transitions!(trans_map)
	X_grid = trans_map.X_grid
	U_grid = trans_map.U_grid
	Y_grid = trans_map.Y_grid
	tstep = cont_sys.tstep
	n_trans = 0

	for u_rp in enumerate_gridspace_ref_pos(U_grid)
		u = get_coords_by_pos(U_grid, u_rp[2])
		r = X_grid.h/2 + cont_sys.meas_noise
		cont_sys.bound_map!(r, u, cont_sys.tstep)
		r[:] += cont_sys.meas_noise
		for x_rp in enumerate_gridspace_ref_pos(X_grid)
			x = get_coords_by_pos(X_grid, x_rp[2])
			cont_sys.sys_map!(x, u, tstep)
			y_lbI, y_ubI = get_pos_lim_from_box_outer(X_grid, x - r, x + r)
		    y_pos_coll = _make_iterator_from_lims(y_lbI, y_ubI)
			y_ref_coll = (get_ref_by_pos(X_grid, y_pos) for y_pos in y_pos_coll)
			if all(ref -> ref !== X_grid.overflow_ref, y_ref_coll)
				for y_ref in y_ref_coll
					add_transition_by_ref!(trans_map, x_rp[1], u_rp[1], y_ref)
					n_trans += 1
				end
			end
        end
    end
	println("set_transitions_from_controlsystem! terminated with success: $(n_trans) transitions created")
end

function set_controller_reach!(trans_map_contr, trans_map_sys, X_init, X_reach)
	@assert X_init.grid_space == X_reach.grid_space ==
		trans_map_contr.X_grid == trans_map_contr.Y_grid
		trans_map_sys.X_grid == trans_map_sys.Y_grid
	println("set_controller_reach! started")
	X_grid = trans_map_contr.X_grid
	X_remain = NewSubSpace(X_grid)
	for x_ref in enumerate_gridspace_ref(X_grid)
		if is_xref_controllable(trans_map_sys, x_ref)
			add_to_subspace_by_ref!(X_remain, x_ref)
		end
	end
	X_init2 = NewSubSpace(X_grid)
	for x_ref in enumerate_subspace_ref(X_init)
		add_to_subspace_by_ref!(X_init2, x_ref)
	end
	X_reach2 = NewSubSpace(X_grid)
	for x_ref in enumerate_subspace_ref(X_reach)
		add_to_subspace_by_ref!(X_reach2, x_ref)
		remove_from_subspace_by_ref!(X_remain, x_ref)
		remove_from_subspace_by_ref!(X_init2, x_ref)
	end
	TypeRef = typeof(X_grid.overflow_ref)

	while ~is_subspace_empty(X_init2)
		x_ref_to_add = Set{TypeRef}()
		for x_ref in enumerate_subspace_ref(X_remain)
			uy_ref_coll = get_transition_image(trans_map_sys, x_ref)
			for (u_ref, y_ref_coll) in uy_ref_coll
				if all(is_ref_in_subspace(X_reach2, y_ref) for y_ref in y_ref_coll)
					push!(x_ref_to_add, x_ref)
					for y_ref in y_ref_coll
						add_transition_by_ref!(trans_map_contr, x_ref, u_ref, y_ref)
					end
				end
			end
		end
		if isempty(x_ref_to_add)
			println("set_controller_reach! terminated without covering init set")
			return
		end
		for x_ref in x_ref_to_add
			add_to_subspace_by_ref!(X_reach2, x_ref)
			remove_from_subspace_by_ref!(X_remain, x_ref)
			remove_from_subspace_by_ref!(X_init2, x_ref)
		end
	end
	println("set_controller_reach! terminated with success")
end
