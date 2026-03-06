function unfold_pattern_matching(latencies_df, data_residuals_continous, c_erp, evts, cfg)
    #epoch residue
    evts_s = @subset(evts, :event .== 'S')
    data_residuals_epoched, times = Unfold.epoch(
        data = data_residuals_continous,
        tbl = evts_s,
        τ = cfg.epoch_range,
        sfreq = cfg.sfreq,
    )
    n, data_residuals_epoched = Unfold.drop_missing_epochs(evts_s, data_residuals_epoched)

    xc, result, onset = findxcorrpeak(data_residuals_epoched[1, :, :], c_erp)

    for (i, row) in enumerate(eachrow(latencies_df))
        if (row.fixed)
            continue
        end
        row.latency = result[i]
    end

    return latencies_df, xc, onset
end

function unfold_decomposition(data, evts_with_c, cfg; fit_kwargs)
    #unfold deconvolution; TODO: make the fit more general, i.e. let the user provide the model structure
    m = fit(
        UnfoldModel,
        [
            'S' => (cfg.formulas[1], firbasis(cfg.s_range, cfg.sfreq, "")), # TODO: Let user supply bfdict
            'R' => (cfg.formulas[2], firbasis(cfg.r_range, cfg.sfreq, "")),
            'C' => (
                cfg.formulas[3],
                firbasis(c_range_adjusted(cfg.c_range), cfg.sfreq, ""),
            ),
        ],
        evts_with_c,
        data;
        fit_kwargs...
    )
    c_table = coeftable(m)
    erps = extract_erps_from_coeftable(c_table, size(data, 1), ['S', 'R', 'C'])
    s_erp = erps['S']
    r_erp = erps['R']
    c_erp = erps['C']

    yhat = predict(m, exclude_basis = 'C', overlap = true)
    y = data
    residuals_without_SR = Unfold._residuals(UnfoldModel, yhat, y)

    return s_erp, r_erp, c_erp, residuals_without_SR, m
end
