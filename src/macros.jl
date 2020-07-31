# "Set" assumes sym_model is empty initially... May fail if not respected
function set_symmodel_from_controlsystem!(sym_model, cont_sys)
	println("set_symmodel_from_controlsystem! sarted")
	X_grid = sym_model.X_grid
	U_grid = sym_model.U_grid
	Y_grid = sym_model.Y_grid
	tstep = cont_sys.tstep
	n_trans = 0

	for u_rp in enumerate_gridspace_ref_pos(U_grid)
		u = get_coords_by_pos(U_grid, u_rp[2])
		r = X_grid.h./2 .+ cont_sys.meas_noise
		r = cont_sys.bound_map(r, u, cont_sys.tstep)
		r = r .+ cont_sys.meas_noise
		for x_rp in enumerate_gridspace_ref_pos(X_grid)
			x = get_coords_by_pos(X_grid, x_rp[2])
			x = cont_sys.sys_map(x, u, tstep)
            rectI = get_pos_lims_outer(X_grid, HyperRectangle(x .- r, x .+ r))
		    ypos_iter = _make_iterator_from_lims(rectI)
			yref_iter = (get_ref_by_pos(X_grid, y_pos) for y_pos in ypos_iter)
			if any(yref_iter .=== X_grid.overflow_ref)
				continue
			end
			add_to_symmodel_by_new_refs_coll!(sym_model, (x_rp[1], u_rp[1], y_ref) for y_ref in yref_iter)
			n_trans += length(yref_iter)
        end
    end
	println("set_symmodel_from_controlsystem! terminated with success: $(n_trans) transitions created")
end

# "Set" assumes sym_model_contr is empty initially... May fail if not respected
function set_controller_reach!(sym_model_contr, sym_model_sys, X_init, X_target)
	@assert X_init.grid_space == X_target.grid_space ==
		sym_model_contr.X_grid == sym_model_contr.Y_grid
		sym_model_sys.X_grid == sym_model_sys.Y_grid
	println("set_controller_reach! started")
	# Commented sizehints because seems not to improve performances
	# sizehint_symmodel!(sym_model_contr, get_symmodel_size(sym_model_sys))
	X_grid = sym_model_sys.X_grid
	U_grid = sym_model_sys.U_grid
	X_remain = NewSubSet(X_grid)
	for x_ref in enumerate_gridspace_ref(X_grid)
		if is_xref_controllable(sym_model_sys, x_ref)
			add_to_subset_by_new_ref!(X_remain, x_ref)
		end
	end
	# size_remain = get_subset_size(X_remain)
	X_init2 = NewSubSet(X_grid)
	union_subsets!(X_init2, X_init)
	X_target2 = NewSubSet(X_grid)
	# sizehint_subset!(X_target2, size_remain)
	union_subsets!(X_target2, X_target)
	setdiff_subsets!(X_remain, X_target)
	setdiff_subsets!(X_init2, X_target)
	X_target_new = NewSubSet(X_grid)
	U_enabled = NewSubSet(U_grid)

	while ~is_subset_empty(X_init2)
		print(".")
		remove_from_subset_all!(X_target_new)
		for x_ref in enumerate_subset_ref(X_remain)
			remove_from_subset_all!(U_enabled)
			add_inputs_by_xref_ysub!(U_enabled, sym_model_sys, x_ref, X_target2)
			if is_subset_empty(U_enabled)
				continue
			end
			add_to_subset_by_new_ref!(X_target_new, x_ref)
			add_to_symmodel_by_new_refs_coll!(sym_model_contr,
				(x_ref, u_ref, x_ref) for u_ref in enumerate_subset_ref(U_enabled))
		end
		if is_subset_empty(X_target_new)
			println("\nset_controller_reach! terminated without covering init set")
			return
		end
		union_new_subsets!(X_target2, X_target_new)
		setdiff_subsets!(X_remain, X_target_new)
		setdiff_subsets!(X_init2, X_target_new)
	end
	println("\nset_controller_reach! terminated with success")
end
