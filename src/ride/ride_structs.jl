"""
    @with_kw struct RideConfig

A struct holding the configuration values of the RIDE algorithm.

# Fields
- `sfreq::Int`: The sampling frequency of the data.
- `s_range::Vector{Float64}`: The range of the S component. Usually determined
through manual inspection of the data.
- `r_range::Vector{Float64}`: The range of the R component. Usually determined
through manual inspection of the data.
- `c_range::Vector{Float64}`: The range of the C component. Usually determined
through manual inspection of the data.
- `c_estimation_range::Vector{Float64}`: The range used for the intial C 
component latency estimation through peak picking.
- `epoch_range::Vector{Float64}`: The range of one epoch centered around the stimulus onset.
- `iteration_limit::Int = 4`: The maximum number of iterations of the RIDE algorithm. This 
is for the outer decomposition-latency estimation loop.
- `heuristic1::Bool = true`: A flag to enable/disable heuristic 1. This heursitic ensures a 
monoton latency evolution.
- `heuristic2::Bool = true`: A flag to enable/disable heuristic 2. This heuristic randomizes
and fixes a latency on encountering a convex xcorrelation result.
- `heuristic2_rng = MersenneTwister(1234)`: The random number generator used for heuristic 2.
- `heuristic3::Bool = true`: A flag to enable/disable heuristic 3. This heuristic searches for
competing peaks in the xcorrelation results. The peak closest to the previous latency is chosen.
- `heuristic3_threshhold::Float64 = 0.9`: The threshold used for heuristic 3. If the peak is
below this threshold * the maximum peak, it is considered a competing peak.
- `filtering::Bool = true`: A flag to enable/disable filtering of the data before performing the 
cross correlation.
- `save_interim_results::Bool = false`: A flag to enable/disable saving the interim results of
each iteration of the RIDE algorithm.


# Examples
```julia-repl
cfg = RideConfig(
    sfreq = 100,
    s_range = [-0.2, 0.4],
    r_range = [0, 0.8],
    c_range = [-0.4, 0.4],
    c_estimation_range = [-0.1, 0.9],
    epoch_range = [-0.3, 1.6],
    iteration_limit = 5,
    heuristic1 = true,
    heuristic2 = true,
    heuristic3 = true,
    save_interim_results = true,
)
```
"""
@with_kw struct RideConfig
    sfreq::Int
    s_range::Vector{Float64}
    r_range::Vector{Float64}
    c_range::Vector{Float64}
    c_estimation_range::Vector{Float64}
    epoch_range::Vector{Float64}
    iteration_limit::Int = 4
    heuristic1::Bool = true
    heuristic2::Bool = true
    heuristic2_rng = MersenneTwister(1234)
    heuristic3::Bool = true
    heuristic3_threshhold::Float64 = 0.9
    filtering::Bool = true
    save_interim_results::Bool = false
end

"""
    abstract type AbstractMode

An abstract type that defines different modes of the RIDE algorithm.
"""
abstract type AbstractMode end

"""
    struct ClassicMode <: AbstractMode

A struct representing the ClassicMode of RIDE. Pass it to the RIDE algorithm to 
    run in Classic mode.
"""
struct ClassicMode <: AbstractMode end

"""
    struct UnfoldMode <: AbstractMode

A struct representing the UnfoldMode of RIDE. Pass it to the RIDE algorithm to 
    run in Unfold mode.
"""
struct UnfoldMode <: AbstractMode end

"""
    @with_kw mutable struct RideResults

A struct holding the results of a single run of the RIDE algorithm.

# Fields
- `interim_results::Vector{RideResults}`: A vector of `RideResults` structs that 
hold the interim results of each iteration of the RIDE algorithm. This field is 
only filled when the algorithm is run with the `save_interim_results` flag.
- `raw_erp::Array{Float64}`: The raw ERP.
- `s_erp::Array{Float64}`: The ERP of the S component.
- `r_erp::Array{Float64}`: The ERP of the R component.
- `c_erp::Array{Float64}`: The ERP of the C component.
- `s_erp_unpadded::Array{Float64}`: The ERP of the S component, unpadded.
- `r_erp_unpadded::Array{Float64}`: The ERP of the R component, unpadded.
- `c_erp_unpadded::Array{Float64}`: The ERP of the C component, unpadded.
- `c_latencies::Array{Int64}`: The latencies of the C component from the stimulus 
onset.
"""
@with_kw mutable struct RideResults
    interim_results::Vector{RideResults} = []
    raw_erp::Array{Float64} = []
    s_erp::Array{Float64} = []
    r_erp::Array{Float64} = []
    c_erp::Array{Float64} = []
    s_erp_unpadded::Array{Float64} = []
    r_erp_unpadded::Array{Float64} = []
    c_erp_unpadded::Array{Float64} = []
    c_latencies::Array{Int64} = []
end