# The Config Struct

When using UnfoldRIDE you will have to supply the function with a `RideConfig` struct. Additionally, the result of UnfoldRIDE will be given to you in a `RideResult` struct. 

While the Docstrings is already extensive, we will take the time here to explain every option in detail.

```julia
RideConfig(
    sfreq::Int
    s_range::Vector{Float64}
    r_range::Vector{Float64}
    c_range::Vector{Float64}
    c_estimation_range::Vector{Float64}
    tukey_window::Vector{Float64}
    formulas = [@formula(0 ~ 1), @formula(0 ~ 1), @formula(0 ~ 1)]
    epoch_range::Vector{Float64}
    iteration_limit::Int = 4
    heuristic1::Bool = true
    heuristic2::Bool = true
    heuristic2_rng = MersenneTwister(1234)
    heuristic3::Bool = true
    heuristic3_threshhold::Float64 = 0.9
    filtering::Tuple{Bool,Bool} = (true, true)
    save_interim_results::Bool = false
)
```

## `sfreq`
The sample frequency of your data. Should be an Integer.

## Component and epoch ranges
`s_range`, `c_range`, and `r_range` indicate the different time-windows around each component, the same way you would choose a time-window to cut your data into epochs.

All component ranges are relative to their respective latency (or estimated latency) in the event structure.

Opposed to that `epoch_range` indicated the entire epoch relative to the S component. This range should encompass all activity from all components.

## `c_estimation_range`
This time window is used for the initial estimation of the C component through peak picking. The range is relative to the S component.

## `tukey_window`
During the cross-correlation a tukey window is applied to each epoch, i.e. every value outside of this range is set to 0. This ensures that a C component can only be found in the range indicated here.

## `formulas`
A Vector containing formulas for the S, R, and C component. These will be used for the Unfold regression model. The formulas will be ignored in the ClassicRIDE algorithm.

## The three heuristics

The RIDE algorithm implements three heuristics to improve the stability and reliability of C component latency estimation:

### Heuristic 1: Monotonic Latency Changes
- **Parameter**: `heuristic1::Bool = true`
- **Purpose**: Ensures that C component latencies change monotonically across iterations, preventing oscillatory behavior.
- **How it works**: If a latency change reverses direction compared to the previous two iterations, the change is reverted and the latency is marked as fixed.
- **When to adjust**: Disable if you observe overly conservative latency tracking that prevents convergence.

### Heuristic 2: Convex Cross-Correlation Fix
- **Parameter**: `heuristic2::Bool = true`
- **Parameter**: `heuristic2_rng = MersenneTwister(1234)`
- **Purpose**: Detects and handles convex cross-correlation results that indicate poor peak detection.
- **How it works**: When no peaks are found in the cross-correlation (indicating a convex shape), the latency is randomized using a Gaussian distribution based on previous latencies.
- **When to adjust**: Disable only if you're confident your data always has clear peaks and want to avoid randomization.

### Heuristic 3: Competing Peaks Selection
- **Parameter**: `heuristic3::Bool = true`
- **Parameter**: `heuristic3_threshhold::Float64 = 0.9`
- **Purpose**: Handles cases where multiple peaks of similar magnitude exist in the cross-correlation.
- **How it works**: When multiple peaks exceed `threshold × maximum_peak`, the algorithm selects the peak closest to the previous latency.
- **When to adjust**: Increase the threshold (closer to 1.0) to be more selective about competing peaks, or decrease it (closer to 0.0) to be more inclusive.

## Filtering
During the algorithm two filters are applied to the data. By default, both filters are enabled (because they are present in the original Matlab implementation), but both can be disabled.

### `filtering[1]`: 3Hz Low-Pass Filter
- **When applied**: Only during the cross-correlation step
- **Purpose**: Deals with the intrinsic oscillatory nature of EEG noise (especially alpha waves)
- **When to disable**: If you suspect the filter is affecting your specific frequency bands of interest

### `filtering[2]`: 20Hz Low-Pass Filter
- **When applied**: During the (outer) iteration process for latency estimation
- **Purpose**: Attenuates high frequency noise during latency estimation
- **Important note**: The final component estimation uses the original (i.e., unfiltered) data
- **When to disable**: If you're working with high-frequency EEG components that you want to preserve

# Results Struct

The `ride_algorithm` function returns a `RideResults` struct containing the following fields:

## ERP Components
- `raw_erp::Array{Float64}`: The raw ERP computed from the epoched data
- `s_erp::Array{Float64}`: The ERP of the S component, padded to match epoch length
- `r_erp::Array{Float64}`: The ERP of the R component, padded to match epoch length
- `c_erp::Array{Float64}`: The ERP of the C component, padded to match epoch length

## Unpadded ERP Components
These fields contain the ERPs without padding, preserving their original length:
- `s_erp_unpadded::Array{Float64}`: The ERP of the S component without padding
- `r_erp_unpadded::Array{Float64}`: The ERP of the R component without padding
- `c_erp_unpadded::Array{Float64}`: The ERP of the C component without padding

## Latency Information
- `c_latencies::Array{Int64}`: The latencies of the C component onset (in samples from stimulus onset)

## Interim Results
- `interim_results::Vector{RideResults}`: A vector of `RideResults` structs containing the results from each iteration. This field is only populated when the algorithm is run with the `save_interim_results = true` flag.

## Usage Example
```julia
results = ride_algorithm(UnfoldMode, data, evts, cfg)

# Access ERP components
s_component = results.s_erp
c_component = results.c_erp

# Access C component latencies
c_latencies = results.c_latencies

# Access interim results (if save_interim_results was enabled)
if !isempty(results.interim_results)
    first_iteration = results.interim_results[1]
end
``` 