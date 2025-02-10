using Revise
includet("./UnfoldRIDE.jl")
includet("./simulation/simulate_test_data.jl")
includet("./plotting_methods.jl")
using .UnfoldRIDE

#simulate data
begin 
    sim_inputs = simulation_inputs()
    sim_inputs.noise = PinkNoise(; noiselevel = 1)
    data, evts, data_clean, evts_clean, data_clean_s, data_clean_r, data_clean_c = simulate_default_plus_clean(sim_inputs)
    #plot_first_three_epochs_of_raw_data(data_clean_s, evts);
    #plot_first_three_epochs_of_raw_data(data, evts);
end

#run the ride algorithm on the simulated data
begin
    #ENV["JULIA_DEBUG"] = "UnfoldRIDE"
    #config for ride algorithm
    cfg = RideConfig(
        sfreq = 100,
        s_range = [-0.2, 0.4],
        r_range = [0, 0.8],
        c_range = [-0.4, 0.4], # change to -0.4 , 0.4 or something because it's attached to the latency of C
        c_estimation_range = [-0.1, 0.9],
        epoch_range = [-0.3,1.6],
        epoch_event_name = 'S',
        iteration_limit = 5,
        heuristic1 = true,
        heuristic2 = true,
        heuristic3 = true,
        save_interim_results = true
    )

    save_to_hdf5_ride_format(data, evts, cfg.epoch_range, cfg.epoch_event_name, 'R', cfg.sfreq)

    #remove the C events from the evts table, these will be estimated by the ride algorithm
    evts_without_c = @subset(evts, :event .!= 'C')

    #run the ride algorithm
    results = ride_algorithm(RideOriginal, data, evts_without_c, cfg)
    s_erp = results.s_erp
    r_erp = results.r_erp
    c_erp = results.c_erp
    c_latencies = results.c_latencies

    plot_interim_results(data,  evts, results, cfg)
end

# calculate and plot clean erps from the simulated data
# these represent the optimal result from the algorithm
begin
    evts_clean_s = @subset(evts_clean, :event .== 'S')
    data_epoched_clean = Unfold.epoch(data = data_clean, tbl = evts_clean_s, τ = cfg.epoch_range, sfreq = cfg.sfreq)[1]
    n, data_epoched_clean = Unfold.drop_missing_epochs(evts_clean_s, data_epoched_clean)
    erp_clean = mean(data_epoched_clean, dims = 3)

    data_epoched_clean = Unfold.epoch(data = data_clean_s, tbl = evts_clean_s, τ = cfg.epoch_range, sfreq = cfg.sfreq)[1]
    n, data_epoched_clean = Unfold.drop_missing_epochs(evts_clean_s, data_epoched_clean)
    erp_clean_s = mean(data_epoched_clean, dims = 3)

    data_epoched_clean = Unfold.epoch(data = data_clean_c, tbl = evts_clean_s, τ = cfg.epoch_range, sfreq = cfg.sfreq)[1]
    n, data_epoched_clean = Unfold.drop_missing_epochs(evts_clean_s, data_epoched_clean)
    erp_clean_c = mean(data_epoched_clean, dims = 3)

    data_epoched_clean = Unfold.epoch(data = data_clean_r, tbl = evts_clean_s, τ = cfg.epoch_range, sfreq = cfg.sfreq)[1]
    n, data_epoched_clean = Unfold.drop_missing_epochs(evts_clean_s, data_epoched_clean)
    erp_clean_r = mean(data_epoched_clean, dims = 3)

    #plot the results
    f = Figure()
    Axis(f[1,1], yticks = -100:100)
    raw = lines!(f[1,1],  erp_clean[1,:,1], color = "black")
    s = lines!(f[1,1],  erp_clean_s[1,:,1], color = "blue")
    c = lines!(f[1,1],  erp_clean_c[1,:,1], color = "red")
    r = lines!(f[1,1],  erp_clean_r[1,:,1], color = "green")
    Legend(f[1,2]
        , [raw, s, r, c]
        , ["Raw ERP", "S ERP", "R ERP", "C ERP"]
    )
    Label(f[0, :], text = "Simulated clean ERPs")
    display(f)
    save("actual_erps.png", f)
end