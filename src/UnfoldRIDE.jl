module UnfoldRIDE
#using Revise

using UnfoldSim
using CairoMakie
using Random
using Unfold
using UnfoldMakie
using StableRNGs
using Parameters
using HDF5
using DataFrames
using DataFramesMeta
using Statistics
using DSP
using Peaks
using Distributions
using FFTW


include("./ride/ride_structs.jl")
#using .ride_structs
include("./ride/ride_original_algorithm.jl")
include("./ride/ride_original_methods.jl")
include("./ride/ride_unfold_algorithm.jl")
include("./ride/ride_unfold_methods.jl")
include("./ride/ride_shared_methods.jl")
include("./ride/plotting_methods.jl")

export ride_algorithm, ride_algorithm_unfold, ride_config, ride_original, ride_unfold


end
