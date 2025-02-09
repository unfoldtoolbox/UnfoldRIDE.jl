#module ride_structs
#
#    export ride_config, ride_modus, ride_original, ride_unfold
#
    @with_kw struct ride_config
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
    end

    abstract type ride_modus end
    struct ride_original <: ride_modus end
    struct ride_unfold <: ride_modus end

#end