# "Set" assumes sym_model is empty initially... May fail if not respected
function set_symmodel_from_controlsystem!(sym_model, cont_sys)
	println("set_symmodel_from_controlsystem! started")
	X_grid = sym_model.X_grid
	U_grid = sym_model.U_grid
	Y_grid = sym_model.Y_grid
	tstep = cont_sys.tstep
	n_trans = 0

	# Updates every 1 seconds
	@showprogress 1 "Computing symbolic control system: " for u_rp in enumerate_gridspace_ref_pos(U_grid)
		u = get_coords_by_pos(U_grid, u_rp[2])
		r = X_grid.h./2 .+ cont_sys.meas_noise
		r = cont_sys.bound_map(r, u, cont_sys.tstep)
		r = r .+ cont_sys.meas_noise
		for x_rp in enumerate_gridspace_ref_pos(X_grid)
			x = get_coords_by_pos(X_grid, x_rp[2])
			x = cont_sys.sys_map(x, u, tstep)
            rectI = get_pos_lims_outer(X_grid, HyperRectangle(x .- r, x .+ r))
		    ypos_iter = Iterators.product(_ranges(rectI)...)
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
	@assert X_init.grid_space === X_target.grid_space ===
		sym_model_contr.X_grid === sym_model_contr.Y_grid
		sym_model_sys.X_grid === sym_model_sys.Y_grid
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
	xref_new_coll = get_gridspace_reftype(sym_model_sys.X_grid)[]
	uref_enabled_coll = get_gridspace_reftype(sym_model_sys.U_grid)[]

	prog = ProgressUnknown("# iterations computing controller:")
	while !is_subset_empty(X_init2)
		ProgressMeter.next!(prog)
		empty!(xref_new_coll)
		for x_ref in enumerate_subset_ref(X_remain)
			empty!(uref_enabled_coll)
			add_inputs_by_xref_ysub!(uref_enabled_coll, sym_model_sys, x_ref, X_target2)
			if isempty(uref_enabled_coll)
				continue
			end
			push!(xref_new_coll, x_ref)
			add_to_symmodel_by_new_refs_coll!(sym_model_contr,
				(x_ref, u_ref, x_ref) for u_ref in uref_enabled_coll)
		end
		if isempty(xref_new_coll)
			println("\nset_controller_reach! terminated without covering init set")
			return
		end
		add_to_subset_by_new_ref_coll!(X_target2, xref_new_coll)
		remove_from_subset_by_ref_coll!(X_remain, xref_new_coll)
		remove_from_subset_by_ref_coll!(X_init2, xref_new_coll)
	end
	ProgressMeter.finish!(prog)
	println("\nset_controller_reach! terminated with success")
end

# "Set" assumes sym_model_contr is empty initially... May fail if not respected
function set_controller_safe!(sym_model_contr, sym_model_sys, X_init, X_safe)
	@assert X_init.grid_space === X_safe.grid_space ===
		sym_model_contr.X_grid === sym_model_contr.Y_grid
		sym_model_sys.X_grid === sym_model_sys.Y_grid
	println("set_controller_reach! started")
	# Commented sizehints because seems not to improve performances
	# sizehint_symmodel!(sym_model_contr, get_symmodel_size(sym_model_sys))
	X_grid = sym_model_sys.X_grid
	U_grid = sym_model_sys.U_grid
	X_safe2 = NewSubSet(X_grid)
	for x_ref in enumerate_subset_ref(X_safe)
		if is_xref_controllable(sym_model_sys, x_ref)
			add_to_subset_by_new_ref!(X_safe2, x_ref)
		end
	end
	xref_remove_coll = get_gridspace_reftype(sym_model_sys.X_grid)[]
	uref_enabled_coll = get_gridspace_reftype(sym_model_sys.U_grid)[]

	while true
		print(".")
		display(get_subset_size(X_safe2))
		empty!(xref_remove_coll)
		for x_ref in enumerate_subset_ref(X_safe2)
			empty!(uref_enabled_coll)
			add_inputs_by_xref_ysub!(uref_enabled_coll, sym_model_sys, x_ref, X_safe2)
			if isempty(uref_enabled_coll)
				push!(xref_remove_coll, x_ref)
			end
		end
		if isempty(xref_remove_coll)
			break
		end
		remove_from_subset_by_ref_coll!(X_safe2, xref_remove_coll)
	end

	for x_ref in enumerate_subset_ref(X_safe2)
		empty!(uref_enabled_coll)
		add_inputs_by_xref_ysub!(uref_enabled_coll, sym_model_sys, x_ref, X_safe2)
		add_to_symmodel_by_new_refs_coll!(sym_model_contr,
			(x_ref, u_ref, x_ref) for u_ref in uref_enabled_coll)
	end

	if is_subset1_in_subset2(X_init, X_safe2)
		println("\nset_controller_safe! terminated with success")
	else
		println("\nset_controller_safe! terminated without covering init set")
	end
end
