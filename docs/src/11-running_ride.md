# Running Classic and Unfold RIDE

## Data preparation

If you want to use Classic/UnfoldRIDE we assume that your data has already been properly preprocessed.

> [!WARNING]
> Your data must absolutely be preprocessed before using UnfoldRIDE; noisy data can have a drastic influence on variable component estimation.

Additionaly, as opposed to the original [Matlab implementation of RIDE](https://github.com/guangouyang/RIDE), UnfoldRIDE works on continuous data (instead of epoched data). Thus your input data should typically be of `size() = channel x timepoints`. 

## Running Classic/ UnfoldRIDE

If you have your ```data``` and your ```evts_without_c::DataFrame``` we can define a suitable configuration for ride and run the algorithm. The ranges for the individual components have to be determined through manual observation of the data.

```julia
#config for ride algorithm
cfg = RideConfig(
    sfreq = 100,
    s_range = [-0.1, 0.3],
    r_range = [0, 0.4],
    c_range = [-0.4, 0.4],
    tukey_window = (0.3, 0.8),
	formulas = [@formula(0 ~ 1), @formula(0 ~ 1), @formula(0 ~ 1)], c_estimation_range = [0.25, 0.55],
    epoch_range = [-0.1, 1]
)

#run the ride algorithm
#We only have one channel, so we only need the first entry from the results vector.
resultsClassic = ride_algorithm(ClassicMode, data_noisy, evts_without_c, cfg)[1]
resultsUnfold = ride_algorithm(UnfoldMode, data_noisy, evts_without_c, cfg)[1]
```

And that's it. Now, we can plot the results of both algorithm modes.

![Results for Classic and Unfold RIDE](https://github.com/unfoldtoolbox/UnfoldRIDE.jl/blob/initial_commit/docs/images/simulated_EEG_tutorial.png "Results of running Classic and Unfold RIDE on the simulated dataset.")
<details>
<summary>Code used for Graph Creation</summary>

```julia
#plot the results
begin
    f = Figure(size = (1000, 400))

    #plot classic results
    ax = Axis(f[1, 1], yticks = -100:100, title="Classic RIDE")
    raw = lines!(resultsClassic.raw_erp; color = "black", linewidth = 3, label="ERP")
    s = lines!(resultsClassic.s_erp; color = "blue", label="S")
    c = lines!(resultsClassic.c_erp; color = "red", label="C")
    r = lines!(resultsClassic.r_erp; color = "green", label="R")
    
    #plot unfold results
    ax = Axis(f[1, 2], yticks = -100:100, title="Unfold RIDE")
    raw = lines!(resultsUnfold.raw_erp; color = "black", linewidth = 3, label="ERP")
    s = lines!(resultsUnfold.s_erp; color = "blue", label="S")
    c = lines!(resultsUnfold.c_erp; color = "red", label="C")
    r = lines!(resultsUnfold.r_erp; color = "green", label="R")
    axislegend(ax)

    display(f)
end
```
</details>