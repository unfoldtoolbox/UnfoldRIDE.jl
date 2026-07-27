# Simulation with Variable Latency Components

If you want to test or develop a method, it's usually a good idea to simulate some data. This is because in simulated data you know how the results should look like in the end (and in EEG data you never know the "Ground Truth"). Luckily for us, another [package of the Unfold family](https://unfoldtoolbox.github.io/UnfoldDocs/UnfoldSim.jl/stable/) makes this easy for us.

To properly test the RIDE algorithm, we need a dataset including at least one component with a variable latency. [UnfoldSim](https://github.com/unfoldtoolbox/UnfoldSim.jl/tree/main) can be used to generate the EEG data and the ```SequenceDesign``` allows us to define a sequence of components with one shared onset. 

First we will load UnfoldSim and define the sample frequency and noise parameters
```julia
using UnfoldSim # NOTE: You will need the rev="feature-sequentialSamplingModels" branch (07/2026)

# Parameters
sfreq = 100
noise = PinkNoise(; noiselevel = 3)
```


In RIDE, we generally differentiate between three different component clusters: S,C and R:
- S represents the Stimulus
- R represents a response to the Stimulus with a known variable latency
- C represents a response to the Stimulus with an uknown variable latency

We define components for each of these component clusters and use [UnfoldSim](https://github.com/unfoldtoolbox/UnfoldSim.jl/tree/main) to simulate them. If there are questions about the next code sections, we recommend checking out the [UnfoldSim documentation](https://unfoldtoolbox.github.io/UnfoldSim.jl/stable/generated/tutorials/quickstart/#Specify-the-simulation-ingredients).

```julia

p1 = LinearModelComponent(;
basis = p100(; sfreq = sfreq),
formula = @formula(0 ~ 1),
β = [s_amp],
)

n1 = LinearModelComponent(;
    basis = n170(; sfreq = sfreq),
    formula = @formula(0 ~ 1),
    β = [s_amp],
)

p3 = LinearModelComponent(;
    basis = UnfoldSim.hanning(0.2, 0.1, sfreq), # sfreq = 100 for the other bases
    formula = @formula(0 ~ 1 + cond),
    β = [c_amp, 0],
)

resp = LinearModelComponent(;
    #basis = UnfoldSim.hanning(Int(0.5 * sfreq)), # sfreq = 100 for the other bases
    basis = UnfoldSim.hanning(0.2, 0.1, sfreq),
    formula = @formula(0 ~ 1),
    β = [r_amp],
    offset = -10*2.5,
)

components = Dict('S' => [p1, n1], 'C' => [p3], 'R' => [resp])
```

Additionally, we will have to overload the ```SequenceOnset``` function of UnfoldSim so that each component cluster gets a variable onset
```julia
seq = UnfoldSim.SequenceOnset(Dict('S'=>LogNormalOnset(μ=4,σ=0.25),
									   'C'=>UniformOnset(0.1*sfreq,0.1*sfreq),
									   'R'=>UniformOnset(1*sfreq,1*sfreq)))
```


```julia
sequencedesign = RepeatDesign(
    SequenceDesign(
        SingleSubjectDesign(conditions=Dict(:cond=>[:A,:B])),"SCR"
        ),
        20)
```

And lastly we can simulate our data.
```julia
data, events = simulate(MersenneTwister(42), sequencedesign, components,seq, noise)
```

Finally, to run the RIDE algorithm, we need to remove the C events from the evts dataframe. They will be estimated during the algorithm and can then be compared to the actual latencies.

```julia
using DataFrames, DataFramesMeta
#only keep the S and R events, the C events will be calculated by the RIDE algorithm
evts_without_c = @subset(events, :event .== 'S' .|| :event .== 'R')
```