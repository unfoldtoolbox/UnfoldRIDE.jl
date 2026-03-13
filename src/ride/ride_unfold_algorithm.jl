ride_algorithm(Modus::Type{UnfoldMode}, data::Array{Float64}, evts, cfg::RideConfig; kwargs...) =
    ride_algorithm(Modus, reshape(data, (1, :)), evts, cfg; kwargs...)

function ride_algorithm(
    Modus::Type{UnfoldMode},
    data::Array{Float64,2},
    evts,
    cfg::RideConfig;
    kwargs...,
)
    @debug "Running RIDE algorithm with cfg: $cfg"
    @assert cfg.s_range[1] >= cfg.epoch_range[1] && cfg.s_range[2] <= cfg.epoch_range[2] "S range must be within the epoch range"
    @assert cfg.c_estimation_range[1] >= cfg.epoch_range[1] &&
            cfg.c_estimation_range[2] <= cfg.epoch_range[2] "C estimation range must be within the epoch range"

    # Check kwargs
    fit_keys = (
        :fit,
        :contrasts,
        :eventcolumn,
        :solver,
        :show_progress,
        :eventfields,
        :show_warnings,
    )
    fit_kwargs = (; (k => v for (k, v) in pairs(kwargs) if k ∈ fit_keys)...)

    ## data_preparation
    interim_results = Vector{Vector}()
    for i in axes(data, 1)
        push!(interim_results, Vector{RideResults}())
    end
    data_epoched, evts_s, evts_r, evts, number_epochs = prepare_epoch_info(data, evts, cfg)
    raw_erp = mean(data_epoched, dims = 3)[:, :, 1]
    ##

    ## initial unfold deconvolution
    m = fit(
        UnfoldModel,
        [
            'S' => (cfg.formulas[1], firbasis(cfg.s_range, cfg.sfreq, "")), # TODO: Let user supply bfdict
            'R' => (cfg.formulas[2], firbasis(cfg.r_range, cfg.sfreq, "")),
        ],
        evts,
        data;
        fit_kwargs...
    )
    c_table = coeftable(m)
    erps = extract_erps_from_coeftable(c_table, size(data, 1), ['S', 'R'])
    s_erp = erps['S']
    r_erp = erps['R']
    ##

    ## initial residue calculation (data minus S and R)
    yhat = predict(m, overlap = true)
    y = data
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
    c_erp = c_erp[:, :, 1]
    ##

    ## save interim results
    if cfg.save_interim_results
        save_interim_results!(
            interim_results,
            evts,
            raw_erp,
            s_erp,
            r_erp,
            c_erp,
            c_latencies_df,
            cfg,
        )
    end

    ## initial pattern matching with the first calculated c_erp
    c_latencies_df_prev_prev = nothing
    c_latencies_df_prev = deepcopy(c_latencies_df)
    for i in range(1, size(c_latencies_df, 1))
        c_latencies_df[i], xcorr = unfold_pattern_matching(
            c_latencies_df[i],
            residuals_without_SR[i, :],
            c_erp[i, :],
            evts_s,
            cfg,
        )
    end

    evts_c = build_c_evts_table(c_latencies_df, evts_s, cfg)
    ##

    ## save interim results
    if cfg.save_interim_results
        save_interim_results!(
            interim_results,
            evts,
            raw_erp,
            s_erp,
            r_erp,
            c_erp,
            c_latencies_df,
            cfg,
        )
    end


    ## iteration start
    model = nothing # init the model because we want the final model later and Julia only has local scope for loops
    for i in range(1, cfg.iteration_limit)
        ## decompose data into S, R and C components using the current C latencies
        evts_with_c = sort(vcat(evts, evts_c), [:latency])
        s_erp, r_erp, c_erp, residue, model = unfold_decomposition(data, evts_with_c, cfg; fit_kwargs = fit_kwargs)
        ##

        ## update C latencies and apply heuristics
        for n in axes(data, 1)
            ## update C latencies via pattern matching
            if cfg.filtering
                residue[n, :] = dspfilter(residue[n, :], 5, cfg.sfreq)
            end
            c_latencies_df[n], xcorr, onset = unfold_pattern_matching(
                c_latencies_df[n],
                residue[n, :],
                c_erp[n, :],
                evts_s,
                cfg,
            )
            ##

            ## heuristics
            if cfg.heuristic1 && i > 1
                heuristic1_monoton_latency_changes!(
                    c_latencies_df[n],
                    c_latencies_df_prev[n],
                    c_latencies_df_prev_prev[n],
                )
            end

            if cfg.heuristic2
                heuristic2_randomize_latency_on_convex_xcorr!(
                    c_latencies_df[n],
                    c_latencies_df_prev[n],
                    xcorr,
                    cfg.heuristic2_rng,
                )
            end

            if cfg.heuristic3
                heuristic3_pick_closest_xcorr_peak!(
                    c_latencies_df[n],
                    c_latencies_df_prev[n],
                    xcorr;
                    equality_threshold = cfg.heuristic3_threshhold,
                    onset = onset,
                )
            end
        end

        c_latencies_df_prev_prev = deepcopy(c_latencies_df_prev)
        c_latencies_df_prev = deepcopy(c_latencies_df)
        evts_c = build_c_evts_table(c_latencies_df, evts_s, cfg)
        ##

        ## save interim results
        if cfg.save_interim_results
            save_interim_results!(
                interim_results,
                evts,
                raw_erp,
                s_erp,
                r_erp,
                c_erp,
                c_latencies_df,
                cfg,
            )
        end
    end

    # Do a last model fit with the final C latencies to get the final S, R and C ERPs; this one should always be done without the robust solver
    evts_with_c = sort(vcat(evts, evts_c), [:latency])
    fit_keys = (
        :fit,
        :contrasts,
        :eventcolumn,
        :show_progress,
        :eventfields,
        :show_warnings,
    )
    fit_kwargs = (; (k => v for (k, v) in pairs(kwargs) if k ∈ fit_keys)...)
    @show "Running final model fit with final C latencies"
    @show fit_kwargs
    s_erp, r_erp, c_erp, residue, model = unfold_decomposition(data, evts_with_c, cfg; fit_kwargs = fit_kwargs)


    results = Vector{RideResults}()
    for i in axes(data, 1)
        r = create_results(
            evts,
            raw_erp[i, :],
            s_erp[i, :],
            r_erp[i, :],
            c_erp[i, :],
            c_latencies_df[i],
            cfg,
        )
        r.interim_results = interim_results[i]
        push!(results, r)
    end

    return results, model # Return both the full results and the final model
end
