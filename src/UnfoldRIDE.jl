module UnfoldRIDE

using Random
using Unfold
using Parameters
using DataFrames
using DataFramesMeta
using Statistics
using DSP
using Peaks
using Distributions
using FFTW


include("./ride/ride_structs.jl")
include("./ride/ride_original_algorithm.jl")
include("./ride/ride_original_methods.jl")
include("./ride/ride_unfold_algorithm.jl")
include("./ride/ride_unfold_methods.jl")
include("./ride/ride_shared_methods.jl")

export ride_algorithm, ride_algorithm_unfold, RideConfig, RideOriginal, RideUnfold


end