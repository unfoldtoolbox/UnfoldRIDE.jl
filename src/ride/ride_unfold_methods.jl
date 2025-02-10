function unfold_pattern_matching(latencies_df, data_residuals_continous, c_erp, evts, cfg)
    #epoch residue
    evts_s = @subset(evts, :event .== 'S')
    data_residuals_epoched, times = Unfold.epoch(data = data_residuals_continous, tbl = evts_s, τ = cfg.epoch_range, sfreq = cfg.sfreq)
    n, data_residuals_epoched = Unfold.drop_missing_epochs(evts_s, data_residuals_epoched)

    xc, result, onset = findxcorrpeak(data_residuals_epoched[1,:,:],c_erp)
    c_latencies = result .- round(Int,  (c_range_adjusted(cfg.c_range)[1] * cfg.sfreq))

    for (i, row) in enumerate(eachrow(latencies_df))
        if(row.fixed) continue end
        row.latency = c_latencies[i]
    end

    return latencies_df, xc, onset
end


function unfold_decomposition(data, evts_with_c, cfg)
    #unfold deconvolution
    m = fit(UnfoldModel,[
        'S' => (@formula(0~1),firbasis(cfg.s_range,cfg.sfreq,"")),
        'R' => (@formula(0~1),firbasis(cfg.r_range,cfg.sfreq,"")),
        'C' => (@formula(0~1),firbasis(c_range_adjusted(cfg.c_range),cfg.sfreq,""))],
        evts_with_c,data)
    c_table = coeftable(m)
    s_erp = c_table[c_table.eventname .== 'S',:estimate]
    r_erp = c_table[c_table.eventname .== 'R',:estimate]
    c_erp = c_table[c_table.eventname .== 'C',:estimate]

    yhat = predict(m,exclude_basis = 'C', overlap = true)
    y = data
    residuals_without_SR = Unfold._residuals(UnfoldModel,yhat,y)

    return s_erp, r_erp, c_erp, residuals_without_SR
end