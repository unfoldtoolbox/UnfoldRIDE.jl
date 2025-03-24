@with_kw struct RideConfig
    sfreq::Int
    s_range::Vector{Float64}
    r_range::Vector{Float64}
    c_range::Vector{Float64}
    c_estimation_range::Vector{Float64}
    epoch_range::Vector{Float64}
    iteration_limit::Int = 5
    heuristic1::Bool = true
    heuristic2::Bool = true
    heuristic2_rng = MersenneTwister(1234)
    heuristic3::Bool = true
    heuristic3_threshhold::Float64 = 0.9
    filtering::Bool = true
    save_interim_results::Bool = false
end

abstract type AbstractRIDE end
struct ClassicMode <: AbstractRIDE end
struct UnfoldMode <: AbstractRIDE end

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