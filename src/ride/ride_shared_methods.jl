"""
    findxcorrpeak(d::Matrix, kernel::Vector; window = false)

Calculate the cross correlation between the data and the kernel for each epoch and find the peak.

# Arguments
- `data::Matrix{Float64}`: Input data in shape (samples, epochs)
- `kernel::Vector{Float64}`: Kernel to cross correlate with the data.

# Keyword arguments
- `window::Bool = false`: Apply a Hanning window to the kernel before cross-correlation, 
    favoring central values of the kernel.

# Returns
- `xc::Vector{Vector{Float64}}` : Cross correlation result per epoch.
- `maxima::Vector{Int}` : Maxima of the cross correlation per epoch.
- `onset::Int` : Onset of the kernel.
"""
function findxcorrpeak(
    data::Union{Matrix{Float64},Vector{Float64}},
    kernel::Vector{Float64};
    window::Bool = false,
)
    # Apply optional Hanning window to favor central values of kernel
    weightedkernel = window ? kernel .* hanning(length(kernel)) : kernel
    xc::Vector{Vector{Float64}} =
        xcorr.(eachcol(data), Ref(weightedkernel); padmode = :none)
    onset::Int = length(kernel)
    maxima::Vector{Int} = [findmax(x)[2] for x in xc] .- onset
    return xc, maxima, onset
end

"""
    initial_peak_estimation(data_continous::Matrix{Float64}, evts::DataFrame, cfg::RideConfig)

Multi channel version of initial_peak_estimation.

# Arguments
- `data_continous::Matrix{Float64}`: Input data in shape (channels, samples).
- `evts::DataFrame`: Event table. S events are extracted from this table.
- `cfg::RideConfig`: Configuration for the RIDE algorithm.

# Returns
- `latencies_df_vector::Vector{DataFrame}` : Vector of DataFrames containing the estimated latencies for every channel. Latencies are from the start of the epoch to the beginning of the c_range.
"""
function initial_peak_estimation(
    data_continous::Matrix{Float64},
    evts::DataFrame,
    cfg::RideConfig,
)
    latencies_df_vector = Vector()
    for i = 1:size(data_continous, 1)
        push!(latencies_df_vector, initial_peak_estimation(data_continous[i, :], evts, cfg))
    end
    return latencies_df_vector
end

"""
    initial_peak_estimation(data_continous::Array{Float64}, evts::DataFrame, cfg::RideConfig)

Estimates the C latencies of every epoch using peak picking in the given cfg.c_estimation_range.
The resulting latencies are from the start of the epoch to the beginning of the c_range.

# Arguments
- `data_continous::Array{Float64}`: Continuous data on which the peak picking is performed.
- `evts::DataFrame`: Event table. S events are extracted from this table.
- `cfg::RideConfig`: Configuration for the RIDE algorithm.

# Returns
- `latencies_df::DataFrame` : DataFrame containing the estimated latencies for every epoch and a fixed=false column. Latencies are from the start of the epoch to the beginning of the c_range.
"""
function initial_peak_estimation(
    data_continous::Array{Float64},
    evts::DataFrame,
    cfg::RideConfig,
)
    evts_s = @subset(evts, :event .== 'S')
    data_residuals_epoched, times = Unfold.epoch(
        data = data_continous,
        tbl = evts_s,
        τ = cfg.epoch_range,
        sfreq = cfg.sfreq,
    )
    n, data_residuals_epoched = Unfold.drop_missing_epochs(evts_s, data_residuals_epoched)
    c_latencies = Vector{Float64}(undef, size(data_residuals_epoched, 3))

    # Peak estimation for initial C latencies for every epoch
    for a in (1:size(data_residuals_epoched, 3))
        range_start =
            round(Int, (cfg.c_estimation_range[1] - cfg.epoch_range[1]) * cfg.sfreq)
        range_end = round(Int, (cfg.c_estimation_range[2] - cfg.epoch_range[1]) * cfg.sfreq)

        # Validate range is valid and within bounds
        range_start = max(1, range_start)
        range_end = min(size(data_residuals_epoched, 2), range_end)

        if range_start >= range_end
            @warn "Invalid c_estimation_range for epoch $a: range_start=$range_start >= range_end=$range_end"
            c_latencies[a] = range_start
            continue
        end

        range = range_start:range_end
        # Find maximum absolute value in the given c_estimation range
        maximum = findmax(abs.(data_residuals_epoched[1, range, a]))[2]
        # Format latency to be from epoch start to the start of the c_range window
        c_latencies[a] = maximum + range[1] - 1 + round(Int, cfg.c_range[1] * cfg.sfreq)
    end
    latencies_df = DataFrame(latency = c_latencies, fixed = false)
    return latencies_df
end

"""
    c_range_adjusted(c_range::Vector{Float64})

Helper Function to adjust the c_range to be from 0 to the difference of the original c_range.
This simplifies several of the algorithms functions without the need to create a new variable.
"""
function c_range_adjusted(c_range::Vector{Float64})
    return [0, c_range[2] - c_range[1]]
end

"""
    build_c_evts_table(latencies_df_vector::Vector, evts::DataFrame, cfg::RideConfig)

Create C event table by copying S and adding the estimated latency.

This version is used for multi channel data. It calculates the mean of all channels as the estimated latencies.

# Arguments
- `latencies_df_vector::Vector`: Vector of DataFrames containing the Latency field. Latencies from epoch start are expected.
- `evts::DataFrame`: Event table. S events are extracted from this table.

# Returns
- `evts_c::DataFrame` : New event table containing only the C events.
"""
function build_c_evts_table(latencies_df_vector::Vector, evts::DataFrame, cfg::RideConfig)
    # Calculate mean latencies across all channels
    latencies_mean = deepcopy(latencies_df_vector[1])
    for i in range(2, size(latencies_df_vector, 1))
        latencies_mean.latency .+= latencies_df_vector[i].latency
    end
    latencies_mean.latency ./= size(latencies_df_vector, 1)
    # Use the standard build_c_evts_table function with mean latencies
    return build_c_evts_table(latencies_mean, evts, cfg)
end

"""
    build_c_evts_table(latencies_df::DataFrame, evts, cfg::RideConfig)

Create C event table by copying S and adding the estimated latency.

# Arguments
- `latencies_df::DataFrame`: DataFrame containing the Latency field. Latencies from epoch start are expected.
- `evts::DataFrame`: Event table. S events are extracted from this table.

# Returns
- `evts_c::DataFrame` : New event table containing only the C events.
"""
function build_c_evts_table(latencies_df::DataFrame, evts::DataFrame, cfg::RideConfig)
    evts_s = @subset(evts, :event .== 'S')
    @assert size(latencies_df, 1) == size(evts_s, 1) "latencies_df and evts_s must have the same size"
    evts_c = copy(evts_s)
    evts_c[!, :latency] .=
        round.(
            Int,
            evts_s[!, :latency] + latencies_df[!, :latency] .+
            (cfg.epoch_range[1] * cfg.sfreq),
        )
    evts_c[!, :event] .= 'C'
    return evts_c
end

"""
    heuristic1_monoton_latency_changes!(latencies_df::DataFrame, latencies_df_old::DataFrame, latencies_df_old_old::DataFrame)

Assure that the changes in the latencies are monotonic, i.e. they always change in one direction.
If a non monotonic change is detected, revert the change and set the latency as fixed in the DataFrame.

# Arguments
- `latencies_df::DataFrame`: The current latencies DataFrame.
- `latencies_df_old::DataFrame`: The previous latencies DataFrame, used to calculate the new change in latency.
- `latencies_df_old_old::DataFrame`: The latencies DataFrame before the previous one, used to calculate the previous change in latency.
"""
function heuristic1_monoton_latency_changes!(
    latencies_df::DataFrame,
    latencies_df_old::DataFrame,
    latencies_df_old_old::DataFrame,
)
    for (i, row) in enumerate(eachrow(latencies_df))
        if row.fixed
            continue
        end
        prev_change = latencies_df_old.latency[i] - latencies_df_old_old.latency[i]
        new_change = row.latency - latencies_df_old.latency[i]
        if prev_change > 0 && new_change < 0 || prev_change < 0 && new_change > 0
            row.latency = latencies_df_old.latency[i]
            row.fixed = true
            @debug "heuristic1 reverted non-monoton latency change for epoch $i"
        end
    end
end

"""
    heuristic2_randomize_latency_on_convex_xcorr!(latencies_df, latencies_df_old, xcorr, rng = MersenneTwister(1234),)

Assure that the cross correlation is not convex, randomize the latency if it is.

Check if the cross correlation is convex by searching for peaks in the cross correlation.
When no peak is found, the xcorrelation is considered convex and 
the latency is randomized with a gaussian distribution over the previous latencies.

# Arguments
- `latencies_df::DataFrame`: The current latencies DataFrame.
- `latencies_df_old::DataFrame`: The previous latencies DataFrame. 
    Used to calculate the gaussian distribution for the randomization.
- `xcorr::Vector{Vector{Float64}}`: Cross correlation results for every epoch.

# Keyword arguments
- `rng::AbstractRNG = MersenneTwister(1234)`: RNG used for random latency generation.
"""
function heuristic2_randomize_latency_on_convex_xcorr!(
    latencies_df::DataFrame,
    latencies_df_old::DataFrame,
    xcorr::Vector{Vector{Float64}},
    rng::AbstractRNG = MersenneTwister(1234),
)
    @assert size(latencies_df, 1) == size(xcorr, 1) "latencies_df and xcorr must have the same size"
    ##you cannot calculate a standard deviation with less than 2 values
    standard_deviation = std(latencies_df_old.latency)
    if (isnan(standard_deviation))
        standard_deviation = 1
    end
    normal_distribution = Normal(mean(latencies_df_old.latency), standard_deviation)
    for (i, row) in enumerate(eachrow(latencies_df))
        if row.fixed
            continue
        end
        maxima = findmaxima(xcorr[i])
        maxima_length = length(maxima.indices)
        if maxima_length == 0
            row.latency = round(Int, rand(rng, normal_distribution))
            row.fixed = true
            @debug "heuristic2 randomized latency for epoch $i"
        end
    end
end

"""
    heuristic3_pick_closest_xcorr_peak!(latencies_df::DataFrame, latencies_df_old::DataFrame, xcorr::Vector{Vector{Float64}}; equality_threshold::Float64 = 0.9, onset::Int64 = 0)

Check for multiple "competing" peaks in the cross correlation and pick the closest one to the previous latency.

Any peak with a value > maximum * equality_threshold is considered a competing peak.
Latencies marked as fixed are skipped.

# Arguments
- `latencies_df::DataFrame`: The current latencies DataFrame.
- `latencies_df_old::DataFrame`: The previous latencies DataFrame. 
    Necessary to calculate the closest peak to the previous latency.
- `xcorr::Vector{Vector{Float64}}`: Cross correlation results for every epoch.

# Keyword arguments (if needed)
- `equality_threshold::Float64 = 0.9`: Threshold to determine which peaks are considered as competing.
    Any peak with a value > maximum * equality_threshold is considered a competing peak, with maximum being the highest peak.
    The threshold must be between 0 and 1.
- `onset::Int64 = 0`: Onset of the kernel used in the cross correlation. Used to bring the latencies into the proper format.
"""
function heuristic3_pick_closest_xcorr_peak!(
    latencies_df::DataFrame,
    latencies_df_old::DataFrame,
    xcorr::Vector{Vector{Float64}};
    equality_threshold::Float64 = 0.9,
    onset::Int64 = 0,
)
    @assert size(latencies_df, 1) == size(xcorr, 1) "latencies_df and xcorr must have the same size"
    @assert size(latencies_df_old, 1) == size(latencies_df, 1) "latencies_df and latencies_df_old must have the same size"
    @assert equality_threshold > 0 && equality_threshold <= 1 "equality_threshold must be between 0 and 1"
    for (i, row) in enumerate(eachrow(latencies_df))
        # skip fixed entries
        if row.fixed
            continue
        end
        maxima = findmaxima(xcorr[i])
        if isempty(maxima)
            continue
        end
        maximum_value = findmax(maxima.heights)[1]
        competing_peaks = []
        for (j, peak) in enumerate(maxima.indices)
            if xcorr[i][peak] > maximum_value * equality_threshold
                push!(competing_peaks, peak - onset)
            end
        end
        if isempty(competing_peaks)
            continue
        end
        closest_peak = argmin(abs.(competing_peaks .- latencies_df_old.latency[i]))
        row.latency = competing_peaks[closest_peak]
    end
end


"""
    dspfilter(signal_to_filter::Vector{Float64}, filter_at::Int64, sampling_rate::Int64)

Filter the given signal with a lowpass filter at the given frequency and sampling rate.
"""
function dspfilter(
    signal_to_filter::Vector{Float64},
    filter_at::Int64,
    sampling_rate::Int64,
)
    @assert filter_at * 2 < sampling_rate "Filter frequency must be less than half the sampling rate"
    p = round(
        3.3 / (
            min(max(filter_at * 0.25, 2.0), sampling_rate / 2 - filter_at) / sampling_rate
        ),
    )
    order = Int(p)
    order = Int(order ÷ 2 * 2) # ensure even filter order

    f = DSP.Filters.digitalfilter(
        Lowpass(filter_at / (sampling_rate / 2)),
        FIRWindow(DSP.hanning(order)),
    )
    b::Vector{Float64} = filtfilt(f, signal_to_filter)

    return b
end

"""
    save_interim_results!(results_vector::Vector{Vector}, evts, raw_erp::Matrix, s_erp::Matrix, r_erp::Matrix, c_erp::Matrix, c_latencies_df::Vector, cfg::RideConfig)

Helper function to create a RideResults object for every channel and add it to the given results_vector.
"""
function save_interim_results!(
    results_vector::Vector{Vector},
    evts,
    raw_erp::Matrix,
    s_erp::Matrix,
    r_erp::Matrix,
    c_erp::Matrix,
    c_latencies_df::Vector,
    cfg::RideConfig,
)
    for i in axes(results_vector, 1)
        r = create_results(
            evts,
            raw_erp[i, :],
            s_erp[i, :],
            r_erp[i, :],
            c_erp[i, :],
            c_latencies_df[i],
            cfg,
        )
        push!(results_vector[i], r)
    end
end

"""
    create_results(evts, raw_erp::Vector, s_erp::Vector, r_erp::Vector, c_erp::Vector, c_latencies_df, cfg::RideConfig)

Create a RideResults object from the given parameters.
Pads component erps to be the same length as one epoch. Result also includes unpadded versions.
Changes the c_latencies to be from stimulus onset instead of epoch start.

# Returns
- `result::RideResults` : The created RideResults object, the interim_results field is left empty.
"""
function create_results(
    evts,
    raw_erp::Vector,
    s_erp::Vector,
    r_erp::Vector,
    c_erp::Vector,
    c_latencies_df,
    cfg::RideConfig,
)
    evts_s = @subset(evts, :event .== 'S')
    evts_r = @subset(evts, :event .== 'R')

    # pad erps to have the same size as one epoch
    mean_s_latency = round(Int, ((-cfg.epoch_range[1] + cfg.s_range[1]) * cfg.sfreq))
    evts_r_latencies_from_s = (evts_r.latency - evts_s.latency)
    median_r_latency = round(
        Int,
        median(evts_r_latencies_from_s) + (cfg.r_range[1] * cfg.sfreq) -
        (cfg.epoch_range[1] * cfg.sfreq),
    )
    median_c_latency = round(Int, median(c_latencies_df.latency))
    s_erp_padded = pad_erp_to_epoch_size(s_erp, mean_s_latency, cfg)
    r_erp_padded = pad_erp_to_epoch_size(r_erp, median_r_latency, cfg)
    c_erp_padded = pad_erp_to_epoch_size(c_erp, median_c_latency, cfg)

    #add the epoch range to the c_latencies as the output should be latency 
    #from stimulus onset, not from epoch onset
    c_latencies_from_stimulus_onset =
        round.(c_latencies_df.latency .+ (cfg.epoch_range[1] * cfg.sfreq)) # Rounding cause latency has to be in samples

    result = RideResults(
        raw_erp = raw_erp,
        s_erp = s_erp_padded,
        r_erp = r_erp_padded,
        c_erp = c_erp_padded,
        s_erp_unpadded = s_erp,
        r_erp_unpadded = r_erp,
        c_erp_unpadded = c_erp,
        c_latencies = c_latencies_from_stimulus_onset,
    )
    return result
end

"""
    pad_erp_to_epoch_size(erp::Vector{Float64}, latency_from_epoch_start::Int64, cfg::RideConfig)

Pad the given erp to the epoch size using the given latency.

# Returns
- `padded_erp::Vector{Float64}` : A new padded erp with length exactly equal to epoch_length.
"""
function pad_erp_to_epoch_size(
    erp::Vector{Float64},
    latency_from_epoch_start::Int64,
    cfg::RideConfig,
)
    epoch_length = round(Int, (cfg.epoch_range[2] - cfg.epoch_range[1]) * cfg.sfreq)
    padding_front_length = max(round(Int, latency_from_epoch_start), 0)
    padding_front = zeros(Float64, padding_front_length)

    # Calculate back padding to ensure total length equals epoch_length
    padding_back_length = max(epoch_length - length(padding_front) - length(erp), 0)
    padding_back = zeros(Float64, padding_back_length)

    padded_erp = vcat(padding_front, erp, padding_back)

    # Ensure exact length (trim if necessary)
    if length(padded_erp) > epoch_length
        padded_erp = padded_erp[1:epoch_length]
    elseif length(padded_erp) < epoch_length
        padded_erp = vcat(padded_erp, zeros(Float64, epoch_length - length(padded_erp)))
    end

    return padded_erp
end

"""
    extract_erps_from_coeftable(c_table::DataFrame, data_channels::Int64, eventnames::Vector{Char})

Helper to extract ERPs for specified events from an Unfold coefficient table.

# Arguments
- `c_table::DataFrame`: Coefficient table from `coeftable(model)`.
- `data_channels::Int64`: Number of channels in the data.
- `eventnames::Vector{Char}`: Event names to extract (e.g., `['S', 'R', 'C']`).

# Returns
- `erps::Dict{Char, Matrix{Float64}}`: Dictionary mapping event name to its ERP matrix (channels, samples).
"""
function extract_erps_from_coeftable(
    c_table::DataFrame,
    data_channels::Int64,
    eventnames::Vector{Char},
)
    erps = Dict{Char,Matrix{Float64}}()

    for eventname in eventnames
        erp = Matrix{Float64}(
            undef,
            data_channels,
            size(@subset(c_table, :eventname .== eventname, :channel .== 1), 1),
        )
        for i = 1:data_channels
            erp[i, :] = @subset(c_table, :eventname .== eventname, :channel .== i).estimate
        end
        erps[eventname] = erp
    end

    return erps
end


"""
    prepare_epoch_info(data::Array{Float64,2}, evts::DataFrame, cfg::RideConfig)

Helper to consolidate epoch preparation logic shared by both Classic and Unfold modes.
Performs epoching, drop missing epochs, and trims events to match number of epochs.

Returns a tuple `(data_epoched, evts_s, evts_r, evts_trimmed, number_epochs)` where
- `data_epoched` is the epoched data (channels, samples, epochs),
- `evts_s` and `evts_r` are the S/R event tables trimmed to the number of epochs,
- `evts_trimmed` is the original `evts` trimmed to match epochs,
- `number_epochs` is the number of valid epochs.
"""
function prepare_epoch_info(data::Array{Float64,2}, evts::DataFrame, cfg::RideConfig)
    evts_s = @subset(evts, :event .== 'S')
    evts_r = @subset(evts, :event .== 'R')
    evts_temp = deepcopy(evts)
    
    data_epoched, _times =
        Unfold.epoch(data = data, tbl = evts_s, τ = cfg.epoch_range, sfreq = cfg.sfreq)
    _n, data_epoched = Unfold.drop_missing_epochs(evts_s, data_epoched)
    number_epochs = size(data_epoched, 3)

    evts_s = evts_s[1:number_epochs, :]
    evts_r = evts_r[1:number_epochs, :]

    # trim original evts to match expected number of rows (S+R per epoch)
    while size(evts_temp, 1) > number_epochs * 2
        deleteat!(evts_temp, size(evts_temp, 1))
    end
    @assert size(evts_temp, 1) == number_epochs * 2 "Size of evts is $(size(evts_temp,1)) but should be $(number_epochs*2)"

    return data_epoched, evts_s, evts_r, evts_temp, number_epochs
end
