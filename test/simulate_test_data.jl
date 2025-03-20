using UnfoldSim
using Random
using Unfold
using StableRNGs
using Parameters
using HDF5
using DataFrames
using DataFramesMeta
using Statistics
using DSP

@with_kw mutable struct simulation_inputs
    rng::AbstractRNG = MersenneTwister(1234)
    noise::AbstractNoise = NoNoise()
    s_width::Int = 0
    s_offset::Int = 200
    s_beta::Int = 1
    s_continous::Int = 1
    c_width::Int = 30
    c_offset::Int = 30
    c_beta::Int = 5
    r_width::Int = 60
    r_offset::Int = 15
    r_beta::Int = 5
    r_continous::Int = 1
end

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

function simulate_default_plus_clean(simulation_inputs = simulation_inputs())
    data, evts = default_sequence_design(simulation_inputs)

    clean_inputs = deepcopy(simulation_inputs)
    clean_inputs.s_offset = clean_inputs.s_offset + Int(round(clean_inputs.s_width / 2))
    clean_inputs.s_width = 0
    clean_inputs.r_offset = clean_inputs.r_offset + Int(round(clean_inputs.r_width / 2))
    clean_inputs.r_width = 0
    clean_inputs.c_offset = clean_inputs.c_offset + Int(round(clean_inputs.c_width / 2))
    clean_inputs.c_width = 0
    data_clean, evts_clean = default_sequence_design(clean_inputs)

    s_clean_inputs = deepcopy(clean_inputs)
    s_clean_inputs.r_beta = 0
    s_clean_inputs.r_continous = 0
    s_clean_inputs.c_beta = 0
    data_clean_s, n = default_sequence_design(s_clean_inputs)

    r_clean_inputs = deepcopy(clean_inputs)
    r_clean_inputs.s_beta = 0
    r_clean_inputs.s_continous = 0
    r_clean_inputs.c_beta = 0
    data_clean_r, n = default_sequence_design(r_clean_inputs)

    c_clean_inputs = deepcopy(clean_inputs)
    c_clean_inputs.s_beta = 0
    c_clean_inputs.s_continous = 0
    c_clean_inputs.r_beta = 0
    c_clean_inputs.r_continous = 0
    data_clean_c, n = default_sequence_design(c_clean_inputs)

    return data, evts, data_clean, evts_clean, data_clean_s, data_clean_r, data_clean_c
end

function simple_example_design()
    design =
        SingleSubjectDesign(;
            conditions = Dict(
                :condition => ["car", "face"],
                :continuous => range(0, 1, length = 10),
            ),
        ) |> x -> RepeatDesign(x, 4)

    # Define the components
    s_component_p = LinearModelComponent(;
        basis = p100(),
        formula = @formula(0 ~ 1),
        β = [simulation_inputs.s_beta],
    )

    s_component_n = LinearModelComponent(;
        basis = n170(),
        formula = @formula(0 ~ 1 + condition),
        β = [simulation_inputs.s_beta, simulation_inputs.s_continous],
    )

    r_component = LinearModelComponent(;
        basis = p300(),
        formula = @formula(0 ~ 1 + continuous),
        β = [simulation_inputs.r_beta, simulation_inputs.r_continous],
    )

    c_component = LinearModelComponent(;
        basis = p100(), #vcat(p100().*0.3,zeros(7))+vcat(zeros(7),p100()), #
        formula = @formula(0 ~ 1),
        β = [simulation_inputs.c_beta],
    )

    onsetStimulus = UniformOnset(width = 0, offset = 200)

    #onsetC =
    #    UniformOnset(width = simulation_inputs.c_width, offset = simulation_inputs.c_offset)

    #onsetR =
    #    UniformOnset(width = simulation_inputs.r_width, offset = simulation_inputs.r_offset)

    #multi_onset = MultiOnset(onsetStimulus, [onsetC, onsetR])

    components = Dict('S' => [s_component_p, s_component_n])

    data_s, evts_s =
        simulate(MersenneTwister(1234), design, components, onsetStimulus, NoNoise())




end

function default_sequence_design(simulation_inputs = simulation_inputs())
    # Define the design
    design =
        SingleSubjectDesign(;
            conditions = Dict(
                :condition => ["car", "face"],
                :continuous => range(0, 1, length = 10),
            ),
        ) |> x -> RepeatDesign(x, 4)

    sequence_design = SequenceDesign(design, "SCR")

    # Define the components
    s_component_p = LinearModelComponent(;
        basis = p100(),
        formula = @formula(0 ~ 1),
        β = [simulation_inputs.s_beta],
    )

    s_component_n = LinearModelComponent(;
        basis = n170(),
        formula = @formula(0 ~ 1 + condition),
        β = [simulation_inputs.s_beta, simulation_inputs.s_continous],
    )

    r_component = LinearModelComponent(;
        basis = p300(),
        formula = @formula(0 ~ 1 + continuous),
        β = [simulation_inputs.r_beta, simulation_inputs.r_continous],
    )

    c_component = LinearModelComponent(;
        basis = p100(), #vcat(p100().*0.3,zeros(7))+vcat(zeros(7),p100()), #
        formula = @formula(0 ~ 1),
        β = [simulation_inputs.c_beta],
    )

    onsetStimulus =
        UniformOnset(width = simulation_inputs.s_width, offset = simulation_inputs.s_offset)

    onsetC =
        UniformOnset(width = simulation_inputs.c_width, offset = simulation_inputs.c_offset)

    onsetR =
        UniformOnset(width = simulation_inputs.r_width, offset = simulation_inputs.r_offset)

    sequence_onset = SequenceOnset(onsetStimulus, [onsetC, onsetR])

    components = Dict(
        'S' => [s_component_p, s_component_n],
        'C' => [c_component],
        'R' => [r_component],
    )

    data, evts = simulate(
        simulation_inputs.rng,
        sequence_design,
        components,
        sequence_onset,
        simulation_inputs.noise,
    )

    return data, evts
end

function save_to_hdf5_ride_format(data, evts, epoch_range, epoch_char, reaction_char, sfreq)
    evts_epoch_temp = @subset(evts, :event .== epoch_char)
    data_epoched_temp, times =
        Unfold.epoch(data = data, tbl = evts_epoch_temp, τ = epoch_range, sfreq = sfreq)
    evts_epoch, data_epoched =
        Unfold.drop_missing_epochs(evts_epoch_temp, data_epoched_temp)

    #grab the reaction times from the epoched data
    evts_r = @subset(evts, :event .== reaction_char)[!, :latency][1:size(data_epoched, 3)]
    reaction_times = evts_r - evts_epoch[!, :latency]

    #matlab_ride format: Matrix{Float64} with TimeStep:Channel:Trial
    ride_matrix =
        zeros(Float64, size(data_epoched, 2), size(data_epoched, 1), size(data_epoched, 3))

    #fill the ride_matrix with the data
    for x in axes(ride_matrix, 1)
        for y in axes(ride_matrix, 2)
            for z in axes(ride_matrix, 3)
                ride_matrix[x, y, z] = data_epoched[y, x, z]
            end
        end
    end

    #todo maybe create chanlocs and save them here

    h5open("simulated_data.h5", "w") do file
        # Save the 3D array
        write(file, "dataset_data", ride_matrix)

        # Save the vector rt
        write(file, "dataset_rt", Float64.(reaction_times))
    end

end