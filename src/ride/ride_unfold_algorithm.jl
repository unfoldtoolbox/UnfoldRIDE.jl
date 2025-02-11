function ride_algorithm(Modus::Type{RideUnfold}, data, evts, cfg::RideConfig)
    @debug "Running RIDE algorithm with cfg: $cfg"
    ## data_preparation
    data_reshaped = reshape(data, (1, :))
    evts_s = @subset(evts, :event .== 'S')
    evts_r = @subset(evts, :event .== 'R')
    results = RideResults(s_erp = [], r_erp = [], c_erp = [], c_latencies = [])


    #epoch data with the cfg.epoch_range to see how many epochs we have
    #cut evts to match the determined number of epochs
    #the resulting data_epoched is also used for the c latency estimation
    data_epoched, data_epoched_times = Unfold.epoch(
        data = data_reshaped,
        tbl = evts_s,
        τ = cfg.epoch_range,
        sfreq = cfg.sfreq,
    )
    n, data_epoched = Unfold.drop_missing_epochs(evts_s, data_epoched)
    number_epochs = size(data_epoched, 3)
    #@assert size(evts) == (number_epochs * 2) "Size of evts is $(size(evts)) but should be $(number_epochs * 2)"
    evts_s = evts_s[1:number_epochs, :]
    evts_r = evts_r[1:number_epochs, :]

    #reduce evts to the number of epochs
    while size(evts, 1) > number_epochs * 2
        deleteat!(evts, size(evts, 1))
    end
    @assert size(evts, 1) == number_epochs * 2 "Size of evts is $(size(evts,1)) but should be $(number_epochs*2)"
    ##

    ## initial unfold deconvolution
    m = fit(
        UnfoldModel,
        [
            'S' => (@formula(0 ~ 1), firbasis(cfg.s_range, cfg.sfreq, "")),
            'R' => (@formula(0 ~ 1), firbasis(cfg.r_range, cfg.sfreq, "")),
        ],
        evts,
        data_reshaped,
    )
    c_table = coeftable(m)
    s_erp = c_table[c_table.eventname.=='S', :estimate]
    r_erp = c_table[c_table.eventname.=='R', :estimate]
    ##

    ## initial residue calculation (data minus S and R)
    #exclude='C' doesn't seem to do anything. Probably because the model doesn't know what 'C' is. 
    #Evts doesn't contain 'C' at this point.
    yhat = predict(m, exclude_basis = 'C', overlap = true)
    y = data_reshaped
    residuals_without_SR = Unfold._residuals(UnfoldModel, yhat, y)
    ##

    ## initial C latency estimation
    c_latencies_df = initial_peak_estimation(residuals_without_SR, evts_s, cfg)
    evts_c = build_c_evts_table(c_latencies_df, evts_s, cfg)
    ##

    ## calculate first c_erp from initial latencies and residue/data
    data_residuals_c_epoched, times = Unfold.epoch(
        data = residuals_without_SR,
        tbl = evts_c,
        τ = c_range_adjusted(cfg.c_range),
        sfreq = cfg.sfreq,
    )
    n, data_residuals_c_epoched =
        Unfold.drop_missing_epochs(evts_c, data_residuals_c_epoched)
    c_erp = median(data_residuals_c_epoched, dims = 3)
    ##

    ## save interim results
    if cfg.save_interim_results
        save_interim_results!(results, s_erp, r_erp, c_erp[1, :, 1], c_latencies_df)
    end

    ## initial pattern matching with the first calculated c_erp
    c_latencies_df_prev_prev = nothing
    c_latencies_df_prev = deepcopy(c_latencies_df)
    c_latencies_df, xcorr = unfold_pattern_matching(
        c_latencies_df,
        residuals_without_SR,
        c_erp[1, :, 1],
        evts_s,
        cfg,
    )
    evts_c = build_c_evts_table(c_latencies_df, evts_s, cfg)
    ##

    ## save interim results
    if cfg.save_interim_results
        save_interim_results!(results, s_erp, r_erp, c_erp[1, :, 1], c_latencies_df)
    end


    ## iteration start
    for i in range(1, cfg.iteration_limit)
        ## decompose data into S, R and C components using the current C latencies
        evts_with_c = sort(vcat(evts, evts_c), [:latency])
        s_erp, r_erp, c_erp, residue = unfold_decomposition(data_reshaped, evts_with_c, cfg)
        ##

        ## update C latencies via pattern matching
        if cfg.filtering
            residue = dspfilter(residue[1, :], 5, 20)
        end
        c_latencies_df, xcorr, onset =
            unfold_pattern_matching(c_latencies_df, residue, c_erp, evts_s, cfg)
        ##

        ## heuristics
        if cfg.heuristic1 &&
           !isnothing(c_latencies_df_prev) &&
           !isnothing(c_latencies_df_prev_prev)
            heuristic1_monoton_latency_changes!(
                c_latencies_df,
                c_latencies_df_prev,
                c_latencies_df_prev_prev,
            )
        end

        if cfg.heuristic2 && !isnothing(c_latencies_df_prev)
            heuristic2_randommize_latency_on_convex_xcorr!(
                c_latencies_df,
                c_latencies_df_prev,
                xcorr,
                cfg.heuristic2_rng,
            )
        end

        if cfg.heuristic3 && !isnothing(c_latencies_df_prev)
            heuristic3_pick_closest_xcorr_peak!(
                c_latencies_df,
                c_latencies_df_prev,
                xcorr,
                cfg.heuristic3_threshhold,
                onset = onset,
            )
        end

        c_latencies_df_prev_prev = deepcopy(c_latencies_df_prev)
        c_latencies_df_prev = deepcopy(c_latencies_df)
        evts_c = build_c_evts_table(c_latencies_df, evts_s, cfg)
        ##

        ## save interim results
        if cfg.save_interim_results
            save_interim_results!(results, s_erp, r_erp, c_erp, c_latencies_df)
        end
    end

    results.s_erp = s_erp
    results.r_erp = r_erp
    results.c_erp = c_erp
    results.c_latencies = c_latencies_df.latency

    return results
end