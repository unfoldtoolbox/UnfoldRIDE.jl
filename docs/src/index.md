```@meta
CurrentModule = UnfoldRIDE
```

# UnfoldRIDE.jl Documentation 

Welcome to the documentation for [UnfoldRIDE](https://github.com/unfoldtoolbox/UnfoldRIDE.jl), a re-implementation of [RIDE](https://cns.hkbu.edu.hk/RIDE.htm) algorithm in Julia with an extension to replace the RIDEs iterative decomposition with an [Unfold](https://github.com/unfoldtoolbox/Unfold.jl) deconvolution.

If you need more information on BIDS, a quick overview and further reading can be found at [Reference/Brain Imaging Data Structure](./generated/reference/BIDS.md)


```@raw html
<div style="width:60%; margin: auto;">
</div>
```

## Key features
- Component Decomposition using the RIDE method
- C-Latency Estimation
- RIDE with an Unfold deconvolution

## Installation
```julia-repl
julia> using Pkg; Pkg.add("UnfoldRIDE")
```
For more detailed instructions please refer to [Installing Julia & Unfold Packages](https://unfoldtoolbox.github.io/UnfoldDocs/main/installation/).


## Usage example

```Julia
cfg = RideConfig(
    sfreq = 100,
    s_range = [-0.2, 0.4],
    r_range = [0, 0.8],
    c_range = [-0.4, 0.4],
    c_estimation_range = [-0.1, 0.9],
    epoch_range = [-0.3, 1.6],
)

resultsClassic = ride_algorithm(ClassicMode, data, evts, cfg)
resultsUnfold = ride_algorithm(ClassicMode, data, evts, cfg)
```

## Where to start: Learning roadmap
### 1. First steps
🔗 [Data Simulation](10-data_simulation.md)

🔗 [Running UnfoldRIDE](11-running_ride.md)

### 2. Intermediate topics
📌 Goal: 
🔗

### 3. Advanced topics
📌 Goal: 
🔗


## Statement of need
TBD