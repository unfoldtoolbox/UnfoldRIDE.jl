function findxcorrpeak(d, kernel; window = false)
    #the purpose of this method is to find the peak of the cross correlation between the kernel and the data
    #kernel = C component erp. Hanning is applied to factor the center of the C erp more than the edges.
    weightedkernel = window ? kernel .* hanning(length(kernel)) : kernel
    xc = xcorr.(eachcol(d), Ref(weightedkernel); padmode = :none)
    onset = length(kernel)
    m = [findmax(x)[2] for x in xc] .- onset
    return xc, m, onset
end

function filtering20(data, a, b)
    temp = copy(data)
    i = size(data, 2)
    for i in size(data, 3)
        temp[1, :, i] = dspfilter(temp[1, :, i], a, b)
    end
    return temp
end


function initial_peak_estimation(data_residuals_continous::Array{Float64,2}, evts_s, cfg)
    latencies_df_vector = Vector()
    for i in 1:size(data_residuals_continous, 1)
        push!(latencies_df_vector, initial_peak_estimation(data_residuals_continous[i, :], evts_s, cfg))
    end
    return latencies_df_vector
end

function initial_peak_estimation(data_residuals_continous::Array{Float64}, evts_s, cfg)
    ## initial C latency estimation
    data_residuals_epoched, times = Unfold.epoch(
        data = data_residuals_continous,
        tbl = evts_s,
        τ = cfg.epoch_range,
        sfreq = cfg.sfreq,
    )
    n, data_residuals_epoched = Unfold.drop_missing_epochs(evts_s, data_residuals_epoched)
    #Peak estimation/algorithm for initial c latencies
    c_latencies = Matrix{Float64}(undef, 1, size(data_residuals_epoched, 3))
    for a in (1:size(data_residuals_epoched, 3))
        range =
            round.(Int, (cfg.c_estimation_range[1] - cfg.epoch_range[1]) * cfg.sfreq):round(
                Int,
                (cfg.c_estimation_range[2] - cfg.epoch_range[1]) * cfg.sfreq,
            )
        #todo C_range still used here? Does this make sense?
        c_latencies[1, a] =
            (findmax(abs.(data_residuals_epoched[1, range, a])).+range[1].-1)[2] +
            round(Int, cfg.c_range[1] * cfg.sfreq)
    end
    latencies_df = DataFrame(latency = c_latencies[1, :], fixed = false)
    return latencies_df
end

function c_range_adjusted(c_range::Vector{Float64})
    return [0, c_range[2] - c_range[1]]
end

function build_c_evts_table(latencies_df_vector::Vector, evts, cfg)
    #@assert size(latencies_df_vector, 1) > 1 "latencies_df_vector must have at least 2 elements"
    latencies_mean = deepcopy(latencies_df_vector[1])
    for i in range(2,size(latencies_df_vector, 1))
        latencies_mean.latency .+= latencies_df_vector[i].latency
    end
    latencies_mean.latency ./= size(latencies_df_vector, 1)
    return build_c_evts_table(latencies_mean, evts, cfg)
end

#Create C event table by copying S and adding the estimated latency
function build_c_evts_table(latencies_df::DataFrame, evts, cfg)
    evts_s = @subset(evts, :event .== 'S')
    evts_s = first(evts_s, size(latencies_df, 1))
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

function save_interim_results!(results::Vector{}, s_erp, r_erp, c_erp, c_latencies_df)
    for i in 1:size(results, 1)
        save_interim_results!(results[i], s_erp[i, :], r_erp[i, :], c_erp[i, :], c_latencies_df[i])
    end
end

function save_interim_results!(result::RideResults, s_erp, r_erp, c_erp, c_latencies_df)
    temp_result = RideResults(
        s_erp = copy(s_erp),
        r_erp = copy(r_erp),
        c_erp = copy(c_erp),
        c_latencies = copy(c_latencies_df.latency),
    )
    push!(result.interim_results, temp_result)
end

# check for multiple "competing" peaks in the xcorrelation
# any peak with a value > maximum * equality_threshold is considered a competing peak
# the peak closest to the previous latency is chosen
function heuristic3_pick_closest_xcorr_peak!(
    latencies_df,
    latencies_df_old,
    xcorr,
    equality_threshold::Float64;
    onset::Int64 = 0,
)
    @assert size(latencies_df, 1) == size(xcorr, 1) "latencies_df and xcorr must have the same size"
    @assert size(latencies_df_old, 1) == size(latencies_df, 1) "latencies_df and latencies_df_old must have the same size"
    @assert equality_threshold > 0 && equality_threshold <= 1 "equality_threshold must be between 0 and 1"
    for (i, row) in enumerate(eachrow(latencies_df_old))
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
        closest_peak = argmin(abs.(competing_peaks .- row.latency))
        latencies_df.latency[i] = competing_peaks[closest_peak]
    end
end

# check if the xcorrelation is convex by searching for peaks
# when no peak is found, the xcorrelation is considered convex and 
# the latency is randomized with a gaussian distribution over the previous latencies
function heuristic2_randomize_latency_on_convex_xcorr!(
    latencies_df,
    latencies_df_old,
    xcorr,
    rng = MersenneTwister(1234),
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

# make sure the changes in the latencies are monoton
# if a non monoton change is detected, revert the change and set the latency as fixed
function heuristic1_monoton_latency_changes!(
    latencies_df,
    latencies_df_old,
    latencies_df_old_old,
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

function filtering(x::Vector{Float64}, a::Int, b::Int)
    # Ensure x is a column vector
    x = vec(x)
    n = length(x)
    Y = fft(x)  # Compute FFT of input signal

    # Create filter mask H based on the conditions
    H = zeros(ComplexF64, n)

    if b == 0
        H[1] = Y[1]  # Keep only DC component
    elseif a == 0 && b != 0
        H[1:b+1] .= Y[1:b+1]  # Retain low frequencies up to b
        H[end-b+1:end] .= Y[end-b+1:end]  # Retain symmetric frequencies
    else
        H[a+1:b+1] .= Y[a+1:b+1]  # Keep bandpass frequencies
        H[end-b+1:end-a+1] .= Y[end-b+1:end-a+1]  # Retain symmetric part
    end

    # Compute inverse FFT and return the real part
    f = real(ifft(H))
    return f
end

# translated from matlab using gpt
function filtering10(x::Vector{Float64}, a::Int, b::Int)
    # Ensure x is a column vector and remove DC component
    x = vec(x)
    x0 = mean(x)
    x = x .- x0  # Subtract mean to remove DC component

    # Increase FFT resolution by a factor of 10
    n = 10 * length(x)
    x_padded = vcat(x, zeros(n - length(x)))  # Pad with zeros if n is greater than length(x)
    Y = fft(x_padded)

    # Create filter mask H based on conditions
    H = zeros(ComplexF64, n)

    if b == 0
        H[1] = Y[1]  # Keep only DC component
    elseif a == 0 && b != 0
        H[1:b+1] .= Y[1:b+1]  # Retain low frequencies up to b
        H[end-b+1:end] .= Y[end-b+1:end]  # Retain symmetric frequencies
    else
        H[a+1:b+1] .= Y[a+1:b+1]  # Keep bandpass frequencies
        H[end-b+1:end-a+1] .= Y[end-b+1:end-a+1]  # Retain symmetric part
    end

    # Compute inverse FFT and return real part
    f = real(ifft(H))

    # Trim to original signal length
    f = f[1:length(x)]

    # Restore DC component if a == 0
    if a == 0
        f .+= x0
    end

    return f
end

function dspfilter(signal_to_filter, filter_at::Int64, sampling_rate)
    @assert filter_at * 2 < sampling_rate "Filter frequency must be less than half the sampling rate"
    a = signal_to_filter
    #for a in eachcol(signal_to_filter)
    p = round(
        3.3 / (
            min(max(filter_at * 0.25, 2.0), sampling_rate / 2 - filter_at) / sampling_rate
        ),
    )
    order = Int(p)
    order = Int(order ÷ 2 * 2) # we need even filter order

    f = DSP.Filters.digitalfilter(
        Lowpass(filter_at / (sampling_rate / 2)),
        FIRWindow(DSP.hanning(order)),
    )
    b = filtfilt(f, a)
    #end

    return b
end