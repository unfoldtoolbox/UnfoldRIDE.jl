### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# ╔═╡ e0b8127c-c319-11f1-9ef8-d5957f31f505
begin
	import Pkg
    # activate a temporary environment
    Pkg.activate(mktempdir())
    Pkg.add([
       # Pkg.PackageSpec(name="UnfoldSim"),
        Pkg.PackageSpec(name="PlutoUI"),
        Pkg.PackageSpec(name="Unfold"),
        Pkg.PackageSpec(name="UnfoldMakie"),
        Pkg.PackageSpec(name="UnfoldBIDS"),
        Pkg.PackageSpec(name="DataFrames"),
        Pkg.PackageSpec(name="DataFramesMeta"),
        Pkg.PackageSpec(name="CSV"),
        Pkg.PackageSpec(name="HDF5"),
        Pkg.PackageSpec(name="Random"),
        Pkg.PackageSpec(name="CairoMakie"),
        Pkg.PackageSpec(name="Statistics"),
        Pkg.PackageSpec(name="DSP"),
        Pkg.PackageSpec(name="Colors"),
        Pkg.PackageSpec(name="ColorSchemes"),
    ])
end

# ╔═╡ 6ae058bf-264c-4044-a1e8-f5ad20a95196
Pkg.add([
        Pkg.PackageSpec(url="https://github.com/unfoldtoolbox/UnfoldRIDE.jl", rev="tukey+filtering"),
        ])

# ╔═╡ 03e87633-ac93-47a2-aeaa-88f78076d8ac
begin
	using Unfold, UnfoldMakie
	using UnfoldBIDS
	using UnfoldRIDE
	using DataFrames, DataFramesMeta, CSV
	using Random
	using CairoMakie # for plotting
	using PlutoUI
	using Statistics
	using DSP
	using HDF5
end

# ╔═╡ d22da4b2-fc81-4c7f-9a09-e9e262220ecf
PlutoUI.TableOfContents()

# ╔═╡ 546546a4-dd14-469e-8e2f-6073f6613717
md"""
# Load data
"""

# ╔═╡ e3ca3d26-0cc5-40b3-8be6-b9ce9482d9f7
#import data from hdf5. Created for use with the sample data exported from matlab
function import_data_from_hdf5(file_path, channel)
    #imported data is in matlab format: timestep, channel, epoch
    import_data = h5read(file_path, "/dataset_data")
    rt = h5read(file_path, "/dataset_rt")

    data = Vector{Float64}()
    #fill the ride_matrix with the data
    offset = 800
    for i = 1:offset
        push!(data, 0)
    end
    for x in axes(import_data, 3)
        for y in axes(import_data, 1)
            push!(data, import_data[y, channel, x])
        end
        for i = 1:offset
            push!(data, 0)
        end
    end
    s =
        Vector(1:length(rt)) .* (offset + size(import_data, 1)) .- size(import_data, 1) .+
        50
    rt = Int.(round.(rt ./ 2))
    rt = rt .+ s
    evts = DataFrame()
    for i in eachindex(rt)
        push!(evts, (event = 'S', latency = s[i]))
        push!(evts, (event = 'R', latency = rt[i]))
    end

    return data, evts
end

# ╔═╡ 7b0204d8-2f8d-4b6d-af29-af31ab43f9a3
begin
    path = "./data/matlab_ride_samp_face.h5"
    data_channels_vector = Vector()
    for i = 1:65
        data_channel1 = reshape(import_data_from_hdf5(path, i)[1], (1, :))
        push!(data_channels_vector, data_channel1)
    end
    data = reduce(vcat, data_channels_vector)
    evts = import_data_from_hdf5(path, 1)[2]
    @assert(evts == import_data_from_hdf5(path, 2)[2])
end

# ╔═╡ e08d5086-9c45-41d0-b8b0-ac37ebadf066
evts

# ╔═╡ 104a4457-7264-4d29-8d4e-404a524b7df5
md"""
# Run RIDE
"""

# ╔═╡ 11c235fb-229a-493a-bc9e-b8ab53a82853
#config for ride algorithm
cfg = RideConfig(
        sfreq = 500,
        s_range = [0, 0.5],
        r_range = [-0.3, 0.4],
        c_range = [-0.300, 0.700],
        c_estimation_range = [0.1, 0.900],
        tukey_window = [0.200, 0.800],
        epoch_range = [-0.100, 1.000],
        iteration_limit = 4,
        heuristic1 = true,
        heuristic2 = true,
        heuristic3 = true,
        save_interim_results = true,
    )

# ╔═╡ 34f499b5-56d7-417d-a6d0-9b477d3f5023
cfg.filtering

# ╔═╡ ccb18bc3-84ec-481e-abf0-584fbeea6b89
begin
    d = reshape(data[44, :], (1, :))
    
    #run the ride algorithm
    resultsUnfold, UnfoldEvts, model = ride_algorithm(UnfoldMode, d, evts, cfg)

    resultsClassic = ride_algorithm(ClassicMode, d, evts, cfg)
    #plot_interim_results(data[44, :], evts, results[1], cfg)
end

# ╔═╡ 446bb549-9acf-440e-81bf-ebaf13aa5cc4
resultsClassic[1]

# ╔═╡ 146362ee-d655-4220-95bc-5b1d1a258c07
md"""
# Plotting
"""

# ╔═╡ 8a3c7d91-cee9-466e-9b7a-41dd4f790c60
let
	colors = Makie.wong_colors()

	
	f = Figure(;size = (600, 600))
	ax1 = Axis(f[1, 1], ylabel = "S-comp",
			   xticksvisible = false, 
			   xticklabelsvisible = false, title = "Real example data")
    ax2 = Axis(f[2, 1], ylabel = "Amplitude [μV]\nC-comp", 
			   xticksvisible = false, 
			   xticklabelsvisible = false)
	ax3 = Axis(f[3, 1], ylabel = "R-comp", xlabel = "Time [sample point]")

	all_ax = [ax1 ax2 ax3]
	

	for i in 1:3

		if i == 1
			lines!(all_ax[i], resultsClassic[1].s_erp_unpadded, color = colors[1], label = "ClassicRIDE")
			lines!(all_ax[i], resultsUnfold[1].s_erp_unpadded, color = colors[2], label = "UnfoldRIDE")
		elseif i == 2
			lines!(all_ax[i], resultsClassic[1].c_erp_unpadded, color = colors[1], label = "ClassicRIDE")
			lines!(all_ax[i], resultsUnfold[1].c_erp_unpadded, color = colors[2], label = "UnfoldRIDE")
		else
			lines!(all_ax[i], resultsClassic[1].r_erp_unpadded, color = colors[1], label = "ClassicRIDE")
			lines!(all_ax[i], resultsUnfold[1].r_erp_unpadded, color = colors[2], label = "UnfoldRIDE")
		end
	end

	axislegend(all_ax[1]; position = :lt)
	#linkyaxes!(all_ax...)
	hidespines!.(all_ax, (:r, :t)...)
	yspace = maximum(tight_yticklabel_spacing!, all_ax)
	
	all_ax[1].yticklabelspace = yspace
	all_ax[2].yticklabelspace = yspace
	all_ax[3].yticklabelspace = yspace
	f
end

# ╔═╡ 518371cd-4133-4b31-9520-b7308ab69eb0
let
	f = Figure()
	ax1 = Axis(f[1,1], title="Results Unfold")
	ax2 = Axis(f[2,1], title="Results Classic")
	lines!(ax1, resultsUnfold[1].s_erp, label="S")
	lines!(ax1, resultsUnfold[1].c_erp, label="C")
	lines!(ax1, resultsUnfold[1].r_erp, label="R")
	
	lines!(ax2, resultsClassic[1].s_erp, label="S")
	lines!(ax2, resultsClassic[1].c_erp, label="C")
	lines!(ax2, resultsClassic[1].r_erp, label="R")

	linkaxes!(ax1, ax2)
	axislegend(ax1)
	current_figure()
end

# ╔═╡ Cell order:
# ╠═e0b8127c-c319-11f1-9ef8-d5957f31f505
# ╠═6ae058bf-264c-4044-a1e8-f5ad20a95196
# ╠═03e87633-ac93-47a2-aeaa-88f78076d8ac
# ╠═d22da4b2-fc81-4c7f-9a09-e9e262220ecf
# ╟─546546a4-dd14-469e-8e2f-6073f6613717
# ╠═e3ca3d26-0cc5-40b3-8be6-b9ce9482d9f7
# ╠═7b0204d8-2f8d-4b6d-af29-af31ab43f9a3
# ╠═e08d5086-9c45-41d0-b8b0-ac37ebadf066
# ╟─104a4457-7264-4d29-8d4e-404a524b7df5
# ╠═11c235fb-229a-493a-bc9e-b8ab53a82853
# ╠═34f499b5-56d7-417d-a6d0-9b477d3f5023
# ╠═ccb18bc3-84ec-481e-abf0-584fbeea6b89
# ╠═446bb549-9acf-440e-81bf-ebaf13aa5cc4
# ╠═146362ee-d655-4220-95bc-5b1d1a258c07
# ╠═8a3c7d91-cee9-466e-9b7a-41dd4f790c60
# ╠═518371cd-4133-4b31-9520-b7308ab69eb0
