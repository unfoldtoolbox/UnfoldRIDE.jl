function subtract_to_data(data, others_evts_erp_tuples, sfreq)
    data_subtracted = copy(data)
    for (evts, erp, range) in others_evts_erp_tuples
        for i in evts.latency
            sub_range =
                i+round(Int, range[1] * sfreq):i+round(Int, (range[1] * sfreq))+size(
                    erp[1, :, 1],
                )[1]-1
            if sub_range[end] > length(data_subtracted)
                continue
            end
            data_subtracted[1, sub_range] -= erp[1, :, 1]
        end
    end
    return data_subtracted
end

function subtract_to_data_epoched(
    data,
    target_evts,
    target_range,
    others_evts_erp_tuples,
    sfreq,
)
    data_subtracted = subtract_to_data(data, others_evts_erp_tuples, sfreq)
    data_epoched_subtracted, n = Unfold.epoch(
        data = data_subtracted,
        tbl = target_evts,
        τ = target_range,
        sfreq = sfreq,
    )
    n, data_epoched_subtracted =
        Unfold.drop_missing_epochs(target_evts, data_epoched_subtracted)
    #new_erp = median(data_epoched_subtracted, dims = 3)
    return data_epoched_subtracted
end