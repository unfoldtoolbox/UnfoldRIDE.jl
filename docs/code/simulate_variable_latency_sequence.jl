using UnfoldSim
using Parameters
using Random
using CairoMakie
using DataFramesMeta

@with_kw struct SequenceOnset <: AbstractOnset
    stimulus_onset::AbstractOnset
    components_onset::Vector{AbstractOnset}
end

function UnfoldSim.simulate_onsets(rng, onset::SequenceOnset, simulation::Simulation)
    #calculate stimulus onsets
    stimulus_onsets =
        simulate_interonset_distances(rng, onset.stimulus_onset, simulation.design)
    stimulus_offset_accumulated = accumulate(+, stimulus_onsets, dims = 1, init = 0)

    #calculate component offsets
    components_onsets = Vector{Vector{Int}}()
    for obj in onset.components_onset
        Random.seed!(rng, rand(rng, 1:10000))
        push!(components_onsets, simulate_interonset_distances(rng, obj, simulation.design))
    end

    #combine the stimulus offsets and component offsets into one vector
    result = Vector{Int}()
    for i in axes(stimulus_offset_accumulated, 1)
        current_offset = stimulus_offset_accumulated[i]
        push!(result, current_offset)
        for component_onsets in components_onsets
            push!(result, current_offset + component_onsets[i])
        end
    end

    #cut result to the design size
    result = result[1:size(simulation.design)]
    return result
end

#Define the design
design =
    SingleSubjectDesign(; conditions = Dict(:cond => ["car", "face"])) |>
    x -> RepeatDesign(x, 40)

#Create a sequence design with three components (S, C, R)
sequence_design = SequenceDesign(design, "SCR")

# Define the components
s_component_1 =
    LinearModelComponent(; basis = vcat(p100()), formula = @formula(0 ~ 1), β = [1])
s_component_2 =
    LinearModelComponent(; basis = n170(), formula = @formula(0 ~ 1 + cond), β = [1, 0])
c_component =
    LinearModelComponent(; basis = p100(), formula = @formula(0 ~ 1 + cond), β = [-4, 2])
r_component =
    LinearModelComponent(; basis = p100(), formula = @formula(0 ~ 1 + cond), β = [6, 0])

#offset of the stimulus defines the distance between two epochs
#offsets for the components are applied from the stimulus due to the custom simulate_onsets method
#the latencies for the first three epochs would be:
# stimulus  = 100,      200,        300
# c         = 110:140,  210:240,    310:340
# r         = 120:160,  220:260,    320:360
onset_stimulus = UniformOnset(width = 0, offset = 100)
onset_c = UniformOnset(width = 30, offset = 10)
onset_r = UniformOnset(width = 40, offset = 20)
sequence_onset = SequenceOnset(onset_stimulus, [onset_c, onset_r])
#the components dict has to be consistent with the sequence design, i.e. contain S, C, R
components =
    Dict('S' => [s_component_1, s_component_2], 'C' => [c_component], 'R' => [r_component])

#simulate the data
data, evts = simulate(
    MersenneTwister(7),
    sequence_design,
    components,
    sequence_onset,
    PinkNoise(noiselevel = 0.1),
)

#only keep the S and R events, the C events will be calculated by the RIDE algorithm
evts_without_c = @subset(evts, :event .== 'S' .|| :event .== 'R')

#plotting
begin
    f = Figure(size = (1000, 400))
    ax = Axis(
        f[1, 1],
        title = "Simulated EEG data",
        titlesize = 18,
        xlabel = "Time [samples]",
        ylabel = "Amplitude [µV]",
        xlabelsize = 16,
        ylabelsize = 16,
        xgridvisible = false,
        ygridvisible = false,
        limits = ((190, 700), nothing),
    )

    lines!(data; color = "black")

    #plot the event onsets
    evts_s = @subset(evts, :event .== 'S')
    evts_c = @subset(evts, :event .== 'C')
    evts_r = @subset(evts, :event .== 'R')

    vlines!(
        ax,
        evts_s.latency,
        color = "red",
        linestyle = :dash,
        linewidth = 2,
        label = "Stimulus",
    )
    vlines!(
        ax,
        evts_c.latency,
        color = "green",
        linestyle = :dash,
        linewidth = 2,
        label = "C",
    )
    vlines!(
        ax,
        evts_r.latency,
        color = "blue",
        linestyle = :dash,
        linewidth = 2,
        label = "R",
    )
    axislegend("Event onset"; unique = true)
    display(f)
end