function plot_c_latency_estimation_four_epochs(data_epoched, c_latencies, c_erp)
    #plotting the estimated c latencies on a couple of epochs
    f = Figure()
    c_latency_from_epoch_start = c_latencies .- round(Int, cfg.epoch_range[1] * cfg.sfreq)
    for a = 1:2, b = 1:2
        i = (a - 1) * 2 + b
        Axis(f[a, b], title = "Estimated C latency epoch $i")
        lines!(f[a, b], data_epoched[1, :, i]; color = "black")
        y = c_latency_from_epoch_start[i]:(c_latency_from_epoch_start[i]+length(c_erp)-1)
        lines!(f[a, b], y, c_erp; color = "red")
    end
    return f
end

function plot_first_three_epochs_of_raw_data(data, evts)
    f = Figure()
    Axis(f[1, 1], title = "First three epochs")
    evts_s = @subset(evts, :event .== 'S')
    evts_r = @subset(evts, :event .== 'R')
    graph = lines!(first(vec(data), evts_s.latency[4]); color = "black")
    graph_r = vlines!(first(evts_r.latency, 3), color = "blue")
    graph_s = vlines!(first(evts_s.latency, 3), color = "red")

    Legend(
        f[1, 2],
        [graph, graph_r, graph_s],
        ["Data", "Reaction Times", "Stimulus Onsets"],
    )
    display(f)
end

function plot_first_three_epochs_of_raw_data(data)
    f = Figure()
    Axis(f[1, 1], title = "First 600 data points")
    graph = lines!(first(vec(data), 600); color = "black")
    Legend(f[1, 2], [graph], ["Data"])
    display(f)
end

function plot_interim_results(data, evts, results, cfg)
    evts_s = @subset(evts, :event .== 'S')
    data_epoched =
        Unfold.epoch(data = data, tbl = evts_s, τ = cfg.epoch_range, sfreq = cfg.sfreq)[1]
    data_epoched = Unfold.drop_missing_epochs(evts_s, data_epoched)[2]
    raw_erp = mean(data_epoched, dims = 3)[1, :, 1]
    for (i, r) in enumerate(vcat(results.interim_results, results))
        f = plot_c_latency_estimation_four_epochs(
            data_epoched,
            r.c_latencies,
            r.c_erp_unpadded,
        )
        Label(f[0, :], text = "Estimated C latency, Iteration $(i-1)", halign = :center)
        display(f)
    end
    for (i, r) in enumerate(vcat(results.interim_results, results))
        f = plot_simple_results(
            raw_erp,
            r.s_erp,
            r.r_erp,
            r.c_erp;
            legend = true,
            title = "Calculated Components, Iteration $(i-1)",
        )
        display(f)
    end
end

function plot_simple_results(raw_erp, s, r, c; legend = false, f = Figure(), title = "")
    ax = Axis(f[1, 1], yticks = -100:100, title = title)
    raw = lines!(raw_erp; color = "black", linewidth = 3, label = "ERP")
    s = lines!(s; color = "blue", label = "S")
    c = lines!(c; color = "red", label = "C")
    r = lines!(r; color = "green", label = "R")
    if legend
        axislegend(ax)
    end
    return f
end