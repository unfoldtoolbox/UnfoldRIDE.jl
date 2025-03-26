using Pkg
Pkg.activate("./dev")

using Revise
includet("../src/UnfoldRIDE.jl")
includet("../test/simulate_test_data.jl")
includet("./plotting_methods.jl")
using .UnfoldRIDE
using CairoMakie
using StableRNGs
using BenchmarkTools
#simulate data
begin
    sim_inputs = simulation_inputs()
    #sim_inputs.noise = PinkNoise()
    data, evts, data_clean, evts_clean, data_clean_s, data_clean_r, data_clean_c =
        simulate_default_plus_clean(sim_inputs)
end

begin
    #ENV["JULIA_DEBUG"] = "UnfoldRIDE"
    #config for ride algorithm
    cfg = RideConfig(
        sfreq = 100,
        s_range = [-0.2, 0.4],
        r_range = [0, 0.8],
        c_range = [-0.4, 0.4], # change to -0.4 , 0.4 or something because it's attached to the latency of C
        c_estimation_range = [-0.1, 0.9],
        epoch_range = [-0.3, 1.6],
        iteration_limit = 4,
        heuristic1 = true,
        heuristic2 = true,
        heuristic3 = true,
        save_interim_results = false,
    )

    channels = 5
    data_channels_vector = Vector()
    noise = PinkNoise(; noiselevel = 1)
    for i in 1:channels
        data_channel1 = reshape(deepcopy(data), (1,:))
        UnfoldSim.add_noise!(MersenneTwister(i), noise, data_channel1)
        push!(data_channels_vector, data_channel1)
    end
    data_channels = reduce(vcat, data_channels_vector)



    for i in axes(data_channels, 1)
        plot_first_three_epochs_of_raw_data(reshape(data_channels[i,:], (1,:)), evts);
    end


    #remove the C events from the evts table, these will be estimated by the ride algorithm
    evts_without_c = @subset(evts, :event .!= 'C')

    #run the ride algorithm
    @profview results = ride_algorithm(UnfoldMode, data_channels, evts_without_c, cfg)
    for i in axes(results, 1)
        plot_interim_results(reshape(data_channels[i,:], (1,:)), evts, results[i], cfg)
    end
end

if true == false
    @benchmark ride_algorithm(ClassicMode, data_channels, evts_without_c, cfg)
    
    @benchmark ride_algorithm(UnfoldMode, data_channels, evts_without_c, cfg)
end