ride_algorithm(Modus::Type{ClassicMode}, data::Vector{Float64}, evts; kwargs...) =
    ride_algorithm(Modus, data, evts, RideConfig(; kwargs...))

function ride_algorithm(
    Modus::Type{ClassicMode},
    data::Array{Float64,2},
    evts,
    cfg::RideConfig,
)
    r = Vector()
    for i in range(1, size(data, 1))
        push!(r, ride_algorithm(Modus, data[i, :], evts, cfg)[1])
    end
    return r
end

function ride_algorithm(
    Modus::Type{ClassicMode},
    data::Vector{Float64},
    evts,
    cfg::RideConfig,
)
    @debug "Running RIDE algorithm with cfg: $cfg"

    @assert cfg.s_range[1] >= cfg.epoch_range[1] && cfg.s_range[2] <= cfg.epoch_range[2] "S range must be within the epoch range"
    @assert cfg.c_estimation_range[1] >= cfg.epoch_range[1] &&
            cfg.c_estimation_range[2] <= cfg.epoch_range[2] "C estimation range must be within the epoch range"

    ## data_preparation
    data_reshaped = reshape(data, (1, :))
    data_epoched, evts_s, evts_r, evts, number_epochs =
        prepare_epoch_info(data_reshaped, evts, cfg)
    raw_erp = mean(data_epoched, dims = 3)[1, :, 1]
    interim_results = Vector{RideResults}()
    ##

    ## initial C latency estimation
    c_latencies_df = initial_peak_estimation(data, evts_s, cfg)
    evts_c = build_c_evts_table(c_latencies_df, evts_s, cfg)
    c_range_adj = c_range_adjusted(cfg.c_range)
    ##

    ## initial erp calculation
    #calculate initial  erp of S
    data_epoched_s, data_epoched_s_times =
        Unfold.epoch(data = data_reshaped, tbl = evts_s, τ = cfg.s_range, sfreq = cfg.sfreq)
    n, data_epoched_s = Unfold.drop_missing_epochs(evts_s, data_epoched_s)
    s_erp = median(data_epoched_s, dims = 3)

    #calculate initial erp of R
    data_epoched_r, data_epoched_r_times =
        Unfold.epoch(data = data_reshaped, tbl = evts_r, τ = cfg.r_range, sfreq = cfg.sfreq)
    n, data_epoched_r = Unfold.drop_missing_epochs(evts_r, data_epoched_r)
    r_erp = median(data_epoched_r, dims = 3)

    ## save interim results
    if cfg.save_interim_results
        push!(
            interim_results,
            create_results(
                evts,
                raw_erp,
                s_erp[1, :, 1],
                r_erp[1, :, 1],
                zeros(1),
                c_latencies_df,
                cfg,
            ),
        )
    end

    c_latencies_df_prev_prev = nothing
    c_latencies_df_prev = nothing

    ## iteration start
    for i in range(1, cfg.iteration_limit)
        ## decompose data into S, R and C components using the current C latencies
        for j = 1:25
            #calculate erp of C by subtracting S and R from the data
            data_subtracted_s_and_r = subtract_to_data_epoched(
                data_reshaped,
                evts_c,
                c_range_adj,
                [(evts_s, s_erp, cfg.s_range), (evts_r, r_erp, cfg.r_range)],
                cfg.sfreq,
            )
            c_erp = median(data_subtracted_s_and_r, dims = 3)
            #calculate erp of S
            data_subtracted_c_and_r = subtract_to_data_epoched(
                data_reshaped,
                evts_s,
                cfg.s_range,
                [(evts_c, c_erp, c_range_adj), (evts_r, r_erp, cfg.r_range)],
                cfg.sfreq,
            )
            s_erp = median(data_subtracted_c_and_r, dims = 3)
            #calculate erp of R
            data_subtracted_s_and_c = subtract_to_data_epoched(
                data_reshaped,
                evts_r,
                cfg.r_range,
                [(evts_s, s_erp, cfg.s_range), (evts_c, c_erp, c_range_adj)],
                cfg.sfreq,
            )
            r_erp = median(data_subtracted_s_and_c, dims = 3)
        end
        ##

        ## update C latencies via pattern matching
        #perform the pattern matching on the data with the S and R components subtracted
        data_subtracted_s_and_r = subtract_to_data(
            data_reshaped,
            [(evts_s, s_erp, cfg.s_range), (evts_r, r_erp, cfg.r_range)],
            cfg.sfreq,
        )
        if cfg.filtering # TODO: Check if this is correct; also check for filter artefacts
            data_subtracted_s_and_r = dspfilter(data_subtracted_s_and_r[1, :], 5, cfg.sfreq)
        end
        data_epoched_subtracted_s_and_r, n = Unfold.epoch(
            data = data_subtracted_s_and_r,
            tbl = evts_s,
            τ = cfg.epoch_range,
            sfreq = cfg.sfreq,
        )
        n, data_epoched_subtracted_s_and_r =
            Unfold.drop_missing_epochs(evts_s, data_epoched_subtracted_s_and_r)
        xcorr, m, onset =
            findxcorrpeak(data_epoched_subtracted_s_and_r[1, :, :], c_erp[1, :, 1])
        c_latencies = reshape(m, (1, :))
        for (i, row) in enumerate(eachrow(c_latencies_df))
            if (row.fixed)
                continue
            end
            row.latency = c_latencies[i]
        end
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
            heuristic2_randomize_latency_on_convex_xcorr!(
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
                xcorr;
                equality_threshold = cfg.heuristic3_threshhold,
                onset = onset,
            )
        end

        c_latencies_df_prev_prev = deepcopy(c_latencies_df_prev)
        c_latencies_df_prev = deepcopy(c_latencies_df)
        evts_c = build_c_evts_table(c_latencies_df, evts_s, cfg)
        ##

        ## save interim results
        if cfg.save_interim_results
            push!(
                interim_results,
                create_results(
                    evts,
                    raw_erp,
                    s_erp[1, :, 1],
                    r_erp[1, :, 1],
                    c_erp[1, :, 1],
                    c_latencies_df,
                    cfg,
                ),
            )
        end
    end

    ##last iteration using the mean instead of median
    #calculate erp of C by subtracting S and R from the data
    data_subtracted_s_and_r = subtract_to_data_epoched(
        data_reshaped,
        evts_c,
        c_range_adj,
        [(evts_s, s_erp, cfg.s_range), (evts_r, r_erp, cfg.r_range)],
        cfg.sfreq,
    )
    c_erp = mean(data_subtracted_s_and_r, dims = 3)
    #calculate erp of S
    data_subtracted_c_and_r = subtract_to_data_epoched(
        data_reshaped,
        evts_s,
        cfg.s_range,
        [(evts_c, c_erp, c_range_adj), (evts_r, r_erp, cfg.r_range)],
        cfg.sfreq,
    )
    s_erp = mean(data_subtracted_c_and_r, dims = 3)
    #calculate erp of R
    data_subtracted_s_and_c = subtract_to_data_epoched(
        data_reshaped,
        evts_r,
        cfg.r_range,
        [(evts_s, s_erp, cfg.s_range), (evts_c, c_erp, c_range_adj)],
        cfg.sfreq,
    )
    r_erp = mean(data_subtracted_s_and_c, dims = 3)
    ##

    results = create_results(
        evts,
        raw_erp,
        s_erp[1, :, 1],
        r_erp[1, :, 1],
        c_erp[1, :, 1],
        c_latencies_df,
        cfg,
    )
    results.interim_results = interim_results

    return [results]
end