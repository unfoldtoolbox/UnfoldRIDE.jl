#module ride_structs
#
#    export ride_config, ride_modus, ride_original, ride_unfold
#
@with_kw struct RideConfig
    sfreq::Int
    s_range::Vector{Float64}
    r_range::Vector{Float64}
    c_range::Vector{Float64}
    c_estimation_range::Vector{Float64}
    epoch_range::Vector{Float64}
    epoch_event_name::Char
    iteration_limit::Int = 5
    heuristic1::Bool = true
    heuristic2::Bool = true
    heuristic2_rng = MersenneTwister(1234)
    heuristic3::Bool = true
    heuristic3_threshhold::Float64 = 0.9
    filtering::Bool = true
    save_interim_results::Bool = false
end

abstract type RideModus end
struct RideOriginal <: RideModus end
struct RideUnfold <: RideModus end

@with_kw mutable struct RideResults
    interim_results::Vector{RideResults} = Vector{RideResults}()
    s_erp::Vector{Float64}
    r_erp::Vector{Float64}
    c_erp::Vector{Float64}
    c_latencies::Vector{Int64}
end

#end