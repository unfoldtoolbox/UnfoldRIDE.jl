### A Pluto.jl notebook ###
# v1.0.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ b43c1614-3d50-42ae-8b98-fec5980b8d66
begin
	import Pkg
    # activate a temporary environment
    Pkg.activate(mktempdir())
    Pkg.add([
       # Pkg.PackageSpec(name="UnfoldSim"),
        Pkg.PackageSpec(name="PlutoUI"),
        Pkg.PackageSpec(name="Unfold"),
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

# ╔═╡ 5c55af8d-7baa-4dea-9a1c-2564fd5645ca
Pkg.add(url="https://github.com/unfoldtoolbox/UnfoldSim.jl",rev="feature-sequentialSamplingModels")

# ╔═╡ 6aee3f78-4c4d-4147-bf1b-1bec7859cbb7
Pkg.add([
        Pkg.PackageSpec(url="https://github.com/unfoldtoolbox/UnfoldRIDE.jl", rev="tukey+filtering"),
        ])

# ╔═╡ 94194cfa-720f-11f1-ae6e-1f0c5bc3f500
begin
	using UnfoldSim
	using Unfold
	using UnfoldBIDS
	using UnfoldRIDE
	using DataFrames, DataFramesMeta, CSV
	using Random
	using CairoMakie # for plotting
	using PlutoUI
	using Statistics
	using DSP
end

# ╔═╡ 7eb581c3-2c2c-4400-86b3-f4696ee0fe1c
using Colors, ColorSchemes

# ╔═╡ 51d323ce-c439-4ffb-a472-978c438fb832
using HDF5

# ╔═╡ 9d1c647f-6c31-4063-a950-b3c4023110c3


# ╔═╡ 59987f84-1a4c-418b-9be3-137ab9720415
PlutoUI.TableOfContents()

# ╔═╡ fb5e13d7-bd31-42e7-b09d-a3a93f5cb46b
md"""
# Parameters
"""

# ╔═╡ 4e22effa-b136-4b60-b3f2-09f3b66d2970
@bind sub PlutoUI.Select([1, 4, 5, 14, 15, 20])

# ╔═╡ 26c27b46-9b8a-4536-843c-597701ec9e47
@bind r_amp PlutoUI.Select([0, 5, 10, 20, 50, 100])

# ╔═╡ 968c8b50-2a29-4359-b5c4-6347402ee055
@bind c_amp PlutoUI.Select([5, 10, 20, 50, 100, -10, -5])

# ╔═╡ 9379625a-ed23-44e2-a11d-4d2e8920bc1f
nl = 5

# ╔═╡ 23865862-aeb3-42f9-9119-e3d9991a76e6
begin
	noise = PinkNoise(; noiselevel = nl)
	#noise = NoNoise()
	noise = ExponentialNoise(; noiselevel = nl, τ = 4)
	#noise = IntrinsicNoise(10, sfreq, 10) # broken
end

# ╔═╡ adf9ea1e-7120-49f0-bf36-9787c7ebc80f
@bind rep PlutoUI.Select([20, 50, 100])

# ╔═╡ 4107d55b-a92e-4be3-90d6-b579e67eb6b1
#@bind nl PlutoUI.Select([10, 5, 15, 2])

# ╔═╡ 58be6593-096f-4371-b7dc-982691de7b4e
@bind hpf PlutoUI.Select([false, true])

# ╔═╡ 4c258133-d406-4db7-9c6f-426bbc63070a
@bind algo_filt PlutoUI.Select([(true, true),(true, false), (false, false), (false, true)])

# ╔═╡ 4715997f-2e4a-44fc-9dfd-6497a01e61b2
#(true, false), (false, false), (false, true),(true, true)

# ╔═╡ f5a019fa-00ef-4398-8ce1-979e1fc509ee
md"""
## Solver
"""

# ╔═╡ 232dfd69-4877-438c-b4c5-3e3e8b3ed9ff
begin
	_my_prepare = (x, y) -> Unfold.prepare(collect(x), y) |> Unfold.prepare_XTX

    _my_solver!(Ĥ, data, Xt, R_xx, R_xy) = Ĥ .= coef(rlm(R_xx, Xt*data, MEstimator{TukeyLoss}(); initial_scale = :mad,method=:cg))

    se_solver = (x, y) -> Unfold.solver_main(x, y; prepare_fun = 		_my_prepare,solver_fun! = _my_solver!, show_time = true,)
end

# ╔═╡ 758ff988-68f9-4e28-b077-0a5c84b53230
md"""
# Load data and simulate
"""

# ╔═╡ 9ba775f0-fcc7-44ef-b8ec-e95832b87afd
begin
	# Find real events
	
	path = "/store/projects/unfold_duration/local/P3_clean/"
	specific_folder = nothing
	task = nothing
	csv_events = "/store/projects/unfold_duration/local/P3_clean/csv/"
end

# ╔═╡ 36c369e8-1d81-4943-89ce-d243bff8174e
sfreq = 256

# ╔═╡ 8a397871-e9ad-4c12-a1bf-457d59e052d5
md"""
## Components
"""

# ╔═╡ 39161940-9458-44e0-83d9-89219adac04d
#=╠═╡
begin
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

	cresp = LinearModelComponent(;
	    #basis = UnfoldSim.hanning(Int(0.5 * sfreq)), # sfreq = 100 for the other bases
		basis = UnfoldSim.hanning(0.2, 0.1, sfreq),
	    formula = @formula(0 ~ 1),
	    β = [r_amp],
	    offset = -10*2.5,
	)
	
	resp = LinearModelComponent(;
	    #basis = UnfoldSim.hanning(Int(0.5 * sfreq)), # sfreq = 100 for the other bases
		basis = UnfoldSim.hanning(0.2, 0.1, sfreq),
	    formula = @formula(0 ~ 1),
	    β = [r_amp],
	    offset = -10*2.5,
	)
end
  ╠═╡ =#

# ╔═╡ 13ea4788-6800-4f3b-95ef-12125c7c6f46
#=╠═╡
components = Dict('S' => [p1, n1], 'C' => [p3], 'R' => [resp])
  ╠═╡ =#

# ╔═╡ 28ace3a7-136b-43ac-ac95-748f6a17b1c5
md"""
## Continuous data
"""

# ╔═╡ b017f0d0-b28b-4f03-916d-b0b4d763da82
md"""
# Fit RIDE
"""

# ╔═╡ 440fefab-69bf-4ae5-8f5c-c1cfa880387c
md"""
## Config
"""

# ╔═╡ 5f4758d0-6a09-4719-964c-385d19dcad42
cfg = RideConfig(
    #sfreq is the sampling frequency of the data
    sfreq = 256,
    #ranges for the individual components have to be determined through manual inspection of the data
    s_range = [-0.1, 0.2], # s locked
    r_range = [-0.1, 0.4], # r locked
    c_range = [-0.2, 0.4], # mean-c_n-1 latency
    tukey_window = [0.3, 0.4],
	formulas = [@formula(0 ~ 1), @formula(0 ~ 1), @formula(0 ~ 1)],
    #the range in which the initial peak estimation for the C component is performed
    c_estimation_range = [0.25, 0.4],
    #the range for one epoch
    epoch_range = [-0.1, 1.2],
	filtering = algo_filt,
	save_interim_results = true
)

# ╔═╡ a93faa72-819c-457b-9a51-bc9d430e89c4
md"""
## HP function
"""

# ╔═╡ 04152ac0-69ae-4a98-900c-9a7c82c635c9
function hp_filter(    signal_to_filter::Vector{Float64},
    filter_at,
    sampling_rate::Int64,
)
    responsetype = Highpass(0.1)
	designmethod = Butterworth(4)
	
	b = filtfilt(digitalfilter(
		responsetype, 
		designmethod; 
		fs=sampling_rate), 
			 signal_to_filter)

    return b
	
end

# ╔═╡ c229cfcd-1522-4fcc-93d2-9f2e04eaba92
md"""
## Run
"""

# ╔═╡ 247d06f3-2389-4ab4-9fcd-de7dd537f78d
md"""
## Component results
"""

# ╔═╡ 75f0eb62-4ac7-4452-9d65-a651ef982618


# ╔═╡ 60a569c6-72ec-42a5-8a05-0616bcf4868e
seq = UnfoldSim.SequenceOnset(Dict('S'=>LogNormalOnset(μ=4,σ=0.25),
									   'C'=>UniformOnset(0.01*sfreq,0.5*sfreq),
									   'R'=>UniformOnset(1*sfreq,1*sfreq)))

# ╔═╡ f16378da-9d1a-47fa-b3b3-468d812f6dee
@bind seed PlutoUI.Slider(1:42;show_value=true)

# ╔═╡ d25e53d4-459c-49f1-aa1b-837da10e302a
sequencedesign = RepeatDesign(SequenceDesign(SingleSubjectDesign(conditions=Dict(:cond=>[:A,:B])),"SCR"),rep)

# ╔═╡ eb4b96bb-897a-48e9-85cf-c666c17a3386
#=╠═╡
begin

	
	data, events =
	    simulate(MersenneTwister(seed), sequencedesign, components,seq, noise)
	idx = 3000:5000
	
	lines(idx, data[idx]) # plotting
	lats = filter(:latency => <(idx[end]), filter(:latency => >(idx[1]), events))
	@show size(lats)
	vlines!(@rsubset(lats, :event == 'S')[:, :latency]; linestyle = :dash, label="S")
	vlines!(@rsubset(lats, :event == 'C')[:, :latency]; linestyle = :dash, label="C")
	vlines!(@rsubset(lats, :event == 'R')[:, :latency]; linestyle = :dash, label="R")
	axislegend(current_axis())
	current_figure()
end
  ╠═╡ =#

# ╔═╡ c8439f1e-89d5-43de-bf28-5ff960aeb68a
#=╠═╡
begin
	sel_cols = [:event,:cond,:latency] # due to a (now fixed) unfold bug, we have to select our columns and remove all Unioin{Missing,...}
	evts = @rsubset(events, :event == 'S' || :event == 'R')[:,sel_cols]

	if hpf
		d = hp_filter(data, 0.5, sfreq)
	else
		d = data
	end
	
	disallowmissing!(evts,sel_cols)
	#run the ride algorithm
	#resultsClassic = ride_algorithm(ClassicMode, data_noisy, evts_without_c, cfg)
	resultsUnfold, UnfoldEvts, model = ride_algorithm(UnfoldMode, d, evts, cfg)

	resultsClassic, ClassicEvts = ride_algorithm(ClassicMode, d, evts, cfg)
end;
  ╠═╡ =#

# ╔═╡ edfce652-a04c-4424-ab4e-b3dc34afe008
#=╠═╡
let
	colors = Makie.wong_colors()

	
	f = Figure()
	ax1 = Axis(f[1, 1], ylabel = "S-comp",
			   xticksvisible = false, 
			   xticklabelsvisible = false, title = "Sub-"*string(sub))
    ax2 = Axis(f[2, 1], ylabel = "Amplitude [μV]\nC-comp", 
			   xticksvisible = false, 
			   xticklabelsvisible = false)
	ax3 = Axis(f[3, 1], ylabel = "R-comp", xlabel = "Time [sample point]")

	all_ax = [ax1 ax2 ax3]
	

	for i in 1:3

		if i == 1
			lines!(all_ax[i], resultsClassic.s_erp_unpadded, color = colors[1], label = "ClassicRIDE")
			lines!(all_ax[i], resultsUnfold[1].s_erp_unpadded, color = colors[2], label = "UnfoldRIDE")
		elseif i == 2
			lines!(all_ax[i], resultsClassic.c_erp_unpadded, color = colors[1], label = "ClassicRIDE")
			lines!(all_ax[i], resultsUnfold[1].c_erp_unpadded, color = colors[2], label = "UnfoldRIDE")
		else
			lines!(all_ax[i], resultsClassic.r_erp_unpadded, color = colors[1], label = "ClassicRIDE")
			lines!(all_ax[i], resultsUnfold[1].r_erp_unpadded, color = colors[2], label = "UnfoldRIDE")
		end
	end

	axislegend(all_ax[1]; position = :lt)
	linkyaxes!(all_ax...)
	hidespines!.(all_ax, (:r, :t)...)
	yspace = maximum(tight_yticklabel_spacing!, all_ax)
	
	all_ax[1].yticklabelspace = yspace
	all_ax[2].yticklabelspace = yspace
	all_ax[3].yticklabelspace = yspace
	f
end
  ╠═╡ =#

# ╔═╡ b605baa1-b93a-45fb-9fd6-294eb017db23
generate_events(sequencedesign)

# ╔═╡ 510e2294-504d-4b8c-8081-6775a96df76a
#=╠═╡
let
	f = Figure()
	ax1 = Axis(f[1,1], title="Results Unfold")
	ax2 = Axis(f[2,1], title="Results Classic")
	lines!(ax1, resultsUnfold[1].s_erp, label="S")
	lines!(ax1, resultsUnfold[1].c_erp, label="C")
	lines!(ax1, resultsUnfold[1].r_erp, label="R")
	
	lines!(ax2, resultsClassic.s_erp, label="S")
	lines!(ax2, resultsClassic.c_erp, label="C")
	lines!(ax2, resultsClassic.r_erp, label="R")

	linkaxes!(ax1, ax2)
	axislegend(ax1)
	current_figure()
end
  ╠═╡ =#

# ╔═╡ d04674d8-4fec-4000-a732-77be82a4f378
md"""
## C latencies
"""

# ╔═╡ 16976dff-7617-40a6-953d-eeba514ee5f9
md"""
### C Latencies minus GT
"""

# ╔═╡ 96d2448d-4491-41e5-9e28-fd3a4ee4a388
let
	f = Figure(size = (800, 600))
	ax1 = Axis(f[1,1], title = "C Latencies minus ground truth")

	UF = @rsubset(UnfoldEvts, :event == 'C').latency - @rsubset(tmp, :event == 'C').latency[1:end-1]
	CL = @rsubset(ClassicEvts, :event == 'C').latency - @rsubset(tmp, :event == 'C').latency[1:end-1]
	
	lines!(ax1, UF, label = "UnfoldRIDE")
	lines!(ax1, CL, label = "ClassicRIDE")
	hlines!(ax1, [sqrt(mean(UF.^2))], label = "UnfoldRIDE RMSE")
	hlines!(ax1, [sqrt(mean(CL.^2))], label = "ClassicRIDE RMSE")
	#scatter!(ax1, @rsubset(ClassicEvts, :event == 'C').latency, label = "ClassicRIDE")
	#scatter!(ax1, @rsubset(tmp, :event == 'C').latency, label = "GT")
	axislegend(ax1)
	f
end

# ╔═╡ df855033-9015-4c49-b54e-e22d7e9798a0
let
	subs = [1 4 15]#[5, 9, 10, 11, 14, 19, 23]
	f = Figure()
	ax = Axis(f[1,1])

	for s in subs
		density!(ax, @rsubset(all_events[s, :events], :event == 'S').rt, alpha = 0.7, label = string(s))
	end

	axislegend(ax)
	current_figure()
end

# ╔═╡ d299fbc4-f927-4bc4-89ff-60c9f207a0f8
let
	bins = 45
	f = Figure()
	ax = Axis(f[1,1])
	hist!(ax, @rsubset(all_events[14, :events], :event == 'S').rt, bins = bins, color = :teal, alpha = 0.7)
	hist!(ax, @rsubset(all_events[5, :events], :event == 'S').rt, bins = bins, color = :orange, alpha = 0.7)
	hist!(ax, @rsubset(all_events[9, :events], :event == 'S').rt, bins = bins, color = :red, alpha = 0.7)
	hist!(ax, @rsubset(all_events[10, :events], :event == 'S').rt, bins = bins, color = :green, alpha = 0.7)
	hist!(ax, @rsubset(all_events[6, :events], :event == 'S').rt, bins = bins, alpha = 0.7)
	current_figure()
end

# ╔═╡ 325c2900-6607-408e-bb3b-a63615ae75c7
md"""
### C Latencies over iterations
"""

# ╔═╡ 170a45bf-dc86-4ae6-9df4-f414f2834669
#=╠═╡
let
	f = Figure()
	ax = Axis(f[1,1])
	col1 = ColorScheme(range(colorant"blue", colorant"orange", length=length(resultsClassic.interim_results)+1))

	col2 = ColorScheme(range(colorant"green", colorant"red", length=length(resultsClassic.interim_results)+1))
	
	
	for (i, it) in enumerate(resultsClassic.interim_results)
		scatter!(ax, it.c_latencies, label="Iteration "* string(i), color=col1[i])
	end
	scatter!(ax, resultsClassic.c_latencies, label="Final latency", color=col1[end])

	#ax2 = Axis(f[2,1])
	
	for (i, it) in enumerate(resultsUnfold[1].interim_results)
		scatter!(ax, it.c_latencies, label="Iteration "* string(i), color=col2[i])
	end
	scatter!(ax, resultsUnfold[1].c_latencies, label="Final latency", color=col2[end])

	#linkaxes!(ax, ax2)
	axislegend(ax)
	#axislegend(ax2)
	f
end
  ╠═╡ =#

# ╔═╡ 73407145-9fea-454f-b03d-0651c1b2883f
#=╠═╡
let
	f = Figure()
	ax = Axis(f[1,1])
	
	for (i, it) in enumerate(resultsUnfold[1].interim_results)
		scatter!(ax, it.c_latencies, label="Iteration "* string(i))
	end
	scatter!(ax, resultsUnfold[1].c_latencies, label="Final latency")

	axislegend(ax)
	f
end
  ╠═╡ =#

# ╔═╡ 1531239e-63db-4dbb-bcb5-1e44a0ce4592


# ╔═╡ d440e9b8-3a04-4344-8619-686d900d311d
md"""
# Last ClassicRIDE with UnfoldRIDE latencies
"""

# ╔═╡ ab9f6a8f-40d7-44ca-b14d-148f3fe4cbe1
function RIDE_decompo(data, evts, cfg)

    # Prep
    c_range_adj = UnfoldRIDE.c_range_adjusted(cfg.c_range)
    evts_c = @rsubset(evts, :event == 'C')
    evts_s = @rsubset(evts, :event == 'S')
    evts_r = @rsubset(evts, :event == 'R')

    
	data_reshaped = reshape(data, (1, :)) # Use original unfiltered data for the last iteration

    # Need S and R erps first
    ## initial erp calculation
    #calculate initial  erp of S
    data_epoched_s, data_epoched_s_times =
        Unfold.epoch(data = data_reshaped, tbl = evts_s, τ = cfg.s_range, sfreq = cfg.sfreq)
    n, data_epoched_s = Unfold.drop_missing_epochs(evts_s, data_epoched_s)
    s_erp = median(data_epoched_s, dims = 3)

    #calculate initial erp of R
    data_epoched_r, data_epoched_r_times =
        Unfold.epoch(data = data_reshaped, tbl = evts_r, τ = cfg.r_range, sfreq = cfg.sfreq)
    n, data_epoched_r = Unfold.drop_missing_epochs(evts_r, data_epoched_r)
    r_erp = median(data_epoched_r, dims = 3)



    # Calculate actual ERPs
    
    data_subtracted_s_and_r = UnfoldRIDE.subtract_to_data_epoched(
        data_reshaped,
        evts_c,
        c_range_adj,
        [(evts_s, s_erp, cfg.s_range), (evts_r, r_erp, cfg.r_range)],
        cfg.sfreq,
    )
    c_erp = mean(data_subtracted_s_and_r, dims = 3)
    #calculate erp of S
    data_subtracted_c_and_r = UnfoldRIDE.subtract_to_data_epoched(
        data_reshaped,
        evts_s,
        cfg.s_range,
        [(evts_c, c_erp, c_range_adj), (evts_r, r_erp, cfg.r_range)],
        cfg.sfreq,
    )
    s_erp = mean(data_subtracted_c_and_r, dims = 3)
    #calculate erp of R
    data_subtracted_s_and_c = UnfoldRIDE.subtract_to_data_epoched(
        data_reshaped,
        evts_r,
        cfg.r_range,
        [(evts_s, s_erp, cfg.s_range), (evts_c, c_erp, c_range_adj)],
        cfg.sfreq,
    )
    r_erp = mean(data_subtracted_s_and_c, dims = 3)

    return c_erp, s_erp, r_erp
end

# ╔═╡ 651482d2-43b0-4a9b-a873-3060ad9e483d
#=╠═╡
c_gt, s_gt, r_gt = RIDE_decompo(d, events, cfg)
  ╠═╡ =#

# ╔═╡ 82024cf1-56d1-4fc0-bcbe-4be3d6da0293
#=╠═╡
c_cl, s_cl, r_cl = RIDE_decompo(d, ClassicEvts, cfg)
  ╠═╡ =#

# ╔═╡ 23b68cc8-5647-4be6-adf1-bd8ec29fee3c
#=╠═╡
c, s, r = RIDE_decompo(d, UnfoldEvts, cfg)
  ╠═╡ =#

# ╔═╡ b4303d71-0c42-42d5-b144-82a7166488ed
#=╠═╡
let
	lines(c[1,:], label="Unfold Latencies")
	lines!(current_axis(), c_gt[1,:], label="GT")
	lines!(current_axis(), c_cl[1,:], label="Classic Latencies")
	axislegend(current_axis())
	current_figure()

end
  ╠═╡ =#

# ╔═╡ 228d2f12-74b4-4fbc-9f5a-ae55e45ceaaf
#=╠═╡
let
	lines(r[1,:], label="Results")
	lines!(current_axis(), r_gt[1,:], label="GT")
	lines!(current_axis(), r_cl[1,:], label="Classic")
	axislegend(current_axis())
	current_figure()
end
  ╠═╡ =#

# ╔═╡ b399931b-7743-475f-a136-278db788e33a
#=╠═╡
let
	lines(s[1,:], label="Results")
	lines!(current_axis(), s_gt[1,:], label="GT")
	axislegend(current_axis())
	current_figure()

end
  ╠═╡ =#

# ╔═╡ df8974ff-8e84-49fe-93b3-d2221ceb3183
md"""
# Run for multiple seeds
"""

# ╔═╡ 657c79a5-38bb-41b6-9fc5-de3646bd906e
#seeds = 1:2

# ╔═╡ c7c56a5d-8a24-47e2-88da-897619f463a9
#=╠═╡
begin
	rmse_cl = []
	rmse_uf = []
	for s in seeds
	
		subject = @rsubset(all_events, :subject .== sub)[1,:events]
	
		tmpl = subject[!, [:latency, :trialtype, :rt, :event]]
		tmpl = insert_C_rows!(tmp)
		design = MyManualDesign(tmpl)

		dat, evnts =
	    simulate(MersenneTwister(s), design, components, MyManualOnset(), noise)
		
		evnts = @rsubset(evnts, :event == 'S' || :event == 'R')[:,sel_cols]


		if hpf
			d2 = hp_filter(dat, 0.5, sfreq)
		else
			d2 = dat
		end
		
		disallowmissing!(evnts,sel_cols)
		#run the ride algorithm
		#resultsClassic = ride_algorithm(ClassicMode, data_noisy, evts_without_c, cfg)
		resultsUF, UFEvts, m = ride_algorithm(UnfoldMode, d2, evts, cfg)
	
		resultsCL, CLEvts = ride_algorithm(ClassicMode, d2, evts, cfg)

		UF = @rsubset(UFEvts, :event == 'C').latency - @rsubset(tmpl, :event == 'C').latency[1:end-1]
		CL = @rsubset(CLEvts, :event == 'C').latency - @rsubset(tmpl, :event == 'C').latency[1:end-1]

		push!(rmse_uf, sqrt(mean(UF.^2)))
		push!(rmse_cl, sqrt(mean(CL.^2)))
	end

end
  ╠═╡ =#

# ╔═╡ a475da9e-149f-4897-9eb6-f83bc4398654
md"""
## Variance over seeds
"""

# ╔═╡ 181c458e-baa3-4978-81a7-6919d39aaaaa
#=╠═╡
begin
    categories = ["Classic", "Unfold"]
    colors = ColorScheme(range(colorant"red", colorant"green", length=length(categories)))
    
    fig = Figure(size = (600, 400))
    ax = Axis(fig[1, 1], xticks = (1:length(categories), categories))
    all =    [rmse_cl, rmse_uf] 
    for (indx, f) in enumerate(categories)
        datam = all[indx]
        a = fill(indx, length(datam))
        boxplot!(ax, a, datam; whiskerwidth = 1, width = 0.35,
            color = (colors[indx], 0.45), whiskercolor = (colors[indx], 1),
            mediancolor = :black)
    end
    fig

end
  ╠═╡ =#

# ╔═╡ 703de8c2-1205-487c-90fc-bd9086c63a74
md"""
# Save to HDF
"""

# ╔═╡ 99c8067a-2f80-43f5-9255-f7b40e48d269
function save_to_hdf5_ride_format(
    filepath,
	basename, # Base filename
    data,
    evts,
    epoch_range,
    epoch_char,
    reaction_char,
    sfreq,
)
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

    h5open(joinpath(filepath, basename * ".h5"), "w") do file
        # Save the 3D array
        write(file, "dataset_data", ride_matrix)

        # Save the vector rt
        write(file, "dataset_rt", Float64.(reaction_times))
    end

end

# ╔═╡ 50f1ea96-36a8-4429-87ea-216a73a2d241
savepath = "/store/users/skukies/Projects/2026UnfoldRIDE/hdf-data"

# ╔═╡ d238394e-edc8-486c-bff2-9047582258bb
#=╠═╡
bname = "2026-07-23-Simulations-S_amp_" * string(s_amp) * "C_amp_" * string(c_amp) * "R_amp_" * string(r_amp) * "-NL_"* string(nl) * "-TN_" * string(rep*2)
  ╠═╡ =#

# ╔═╡ cb90256d-330d-4e7a-9a71-e765e593af8c
@bind safe PlutoUI.Select([false, true])

# ╔═╡ d3865d74-9da3-4046-81c5-3e2e5f90b698
#=╠═╡
if safe
	save_to_hdf5_ride_format(
	    savepath,
		bname,
	    data,
	    evts,
	    cfg.epoch_range,
	    'S',
	    'R',
	    sfreq,
	)
end
  ╠═╡ =#

# ╔═╡ cfe43446-547d-45ef-a8b2-32695ea60299
#=╠═╡
let
    evts_epoch_temp = @subset(evts, :event .== 'S')
    data_epoched_temp, times =
        Unfold.epoch(data = data, tbl = evts_epoch_temp, τ = cfg.epoch_range, sfreq = sfreq)
    evts_epoch, data_epoched =
        Unfold.drop_missing_epochs(evts_epoch_temp, data_epoched_temp)
    
    evts_r = @subset(evts, :event .== 'R')[!, :latency][1:size(data_epoched, 3)]
    reaction_times = evts_r - evts_epoch[!, :latency]
end
  ╠═╡ =#

# ╔═╡ cfa763c4-5522-4693-9836-bdb137d758bd
#=╠═╡
evts
  ╠═╡ =#

# ╔═╡ ea904b2b-29e4-455e-b172-42f329c3841c
md"""
# Misc
"""

# ╔═╡ c07fd76f-1335-44da-9aa9-b44867baf7a8
let
	#w = Int(round(0.2, 0.12, sfreq))
	h = UnfoldSim.hanning(0.2, 0.1, sfreq)
lines(h)

end

# ╔═╡ 339c13d9-52b3-40bf-abbf-50332a355b41
lines(p100(; sfreq = sfreq))

# ╔═╡ 7ef6ad3b-d101-4e9a-bf8b-4c5f20203d2b
# ╠═╡ disabled = true
#=╠═╡
@bind s_amp PlutoUI.Select([0, 5, 10, 20, 50, 100])
  ╠═╡ =#

# ╔═╡ 035ccc9b-b1fa-4cc9-9ae0-b79af2bd4910
#=╠═╡
@bind s_amp PlutoUI.Select([0, 5, 10, 20, 50, 100]) 
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═b43c1614-3d50-42ae-8b98-fec5980b8d66
# ╠═5c55af8d-7baa-4dea-9a1c-2564fd5645ca
# ╠═6aee3f78-4c4d-4147-bf1b-1bec7859cbb7
# ╠═94194cfa-720f-11f1-ae6e-1f0c5bc3f500
# ╠═9d1c647f-6c31-4063-a950-b3c4023110c3
# ╠═59987f84-1a4c-418b-9be3-137ab9720415
# ╟─fb5e13d7-bd31-42e7-b09d-a3a93f5cb46b
# ╠═23865862-aeb3-42f9-9119-e3d9991a76e6
# ╠═4e22effa-b136-4b60-b3f2-09f3b66d2970
# ╠═26c27b46-9b8a-4536-843c-597701ec9e47
# ╠═968c8b50-2a29-4359-b5c4-6347402ee055
# ╠═035ccc9b-b1fa-4cc9-9ae0-b79af2bd4910
# ╠═9379625a-ed23-44e2-a11d-4d2e8920bc1f
# ╠═adf9ea1e-7120-49f0-bf36-9787c7ebc80f
# ╠═7ef6ad3b-d101-4e9a-bf8b-4c5f20203d2b
# ╠═4107d55b-a92e-4be3-90d6-b579e67eb6b1
# ╠═58be6593-096f-4371-b7dc-982691de7b4e
# ╠═4c258133-d406-4db7-9c6f-426bbc63070a
# ╠═4715997f-2e4a-44fc-9dfd-6497a01e61b2
# ╟─f5a019fa-00ef-4398-8ce1-979e1fc509ee
# ╠═232dfd69-4877-438c-b4c5-3e3e8b3ed9ff
# ╟─758ff988-68f9-4e28-b077-0a5c84b53230
# ╟─9ba775f0-fcc7-44ef-b8ec-e95832b87afd
# ╠═36c369e8-1d81-4943-89ce-d243bff8174e
# ╠═8a397871-e9ad-4c12-a1bf-457d59e052d5
# ╠═39161940-9458-44e0-83d9-89219adac04d
# ╠═13ea4788-6800-4f3b-95ef-12125c7c6f46
# ╟─28ace3a7-136b-43ac-ac95-748f6a17b1c5
# ╠═eb4b96bb-897a-48e9-85cf-c666c17a3386
# ╠═b605baa1-b93a-45fb-9fd6-294eb017db23
# ╟─b017f0d0-b28b-4f03-916d-b0b4d763da82
# ╟─440fefab-69bf-4ae5-8f5c-c1cfa880387c
# ╠═5f4758d0-6a09-4719-964c-385d19dcad42
# ╠═a93faa72-819c-457b-9a51-bc9d430e89c4
# ╠═04152ac0-69ae-4a98-900c-9a7c82c635c9
# ╠═c229cfcd-1522-4fcc-93d2-9f2e04eaba92
# ╠═c8439f1e-89d5-43de-bf28-5ff960aeb68a
# ╠═247d06f3-2389-4ab4-9fcd-de7dd537f78d
# ╠═75f0eb62-4ac7-4452-9d65-a651ef982618
# ╠═60a569c6-72ec-42a5-8a05-0616bcf4868e
# ╟─edfce652-a04c-4424-ab4e-b3dc34afe008
# ╠═f16378da-9d1a-47fa-b3b3-468d812f6dee
# ╠═d25e53d4-459c-49f1-aa1b-837da10e302a
# ╠═510e2294-504d-4b8c-8081-6775a96df76a
# ╟─d04674d8-4fec-4000-a732-77be82a4f378
# ╟─16976dff-7617-40a6-953d-eeba514ee5f9
# ╠═96d2448d-4491-41e5-9e28-fd3a4ee4a388
# ╠═df855033-9015-4c49-b54e-e22d7e9798a0
# ╠═d299fbc4-f927-4bc4-89ff-60c9f207a0f8
# ╠═325c2900-6607-408e-bb3b-a63615ae75c7
# ╟─170a45bf-dc86-4ae6-9df4-f414f2834669
# ╠═73407145-9fea-454f-b03d-0651c1b2883f
# ╠═1531239e-63db-4dbb-bcb5-1e44a0ce4592
# ╠═d440e9b8-3a04-4344-8619-686d900d311d
# ╠═ab9f6a8f-40d7-44ca-b14d-148f3fe4cbe1
# ╠═651482d2-43b0-4a9b-a873-3060ad9e483d
# ╠═82024cf1-56d1-4fc0-bcbe-4be3d6da0293
# ╠═23b68cc8-5647-4be6-adf1-bd8ec29fee3c
# ╠═b4303d71-0c42-42d5-b144-82a7166488ed
# ╠═228d2f12-74b4-4fbc-9f5a-ae55e45ceaaf
# ╠═b399931b-7743-475f-a136-278db788e33a
# ╟─df8974ff-8e84-49fe-93b3-d2221ceb3183
# ╠═657c79a5-38bb-41b6-9fc5-de3646bd906e
# ╠═c7c56a5d-8a24-47e2-88da-897619f463a9
# ╠═7eb581c3-2c2c-4400-86b3-f4696ee0fe1c
# ╟─a475da9e-149f-4897-9eb6-f83bc4398654
# ╠═181c458e-baa3-4978-81a7-6919d39aaaaa
# ╠═703de8c2-1205-487c-90fc-bd9086c63a74
# ╠═51d323ce-c439-4ffb-a472-978c438fb832
# ╠═99c8067a-2f80-43f5-9255-f7b40e48d269
# ╠═50f1ea96-36a8-4429-87ea-216a73a2d241
# ╠═d238394e-edc8-486c-bff2-9047582258bb
# ╠═cb90256d-330d-4e7a-9a71-e765e593af8c
# ╠═d3865d74-9da3-4046-81c5-3e2e5f90b698
# ╠═cfe43446-547d-45ef-a8b2-32695ea60299
# ╠═cfa763c4-5522-4693-9836-bdb137d758bd
# ╠═ea904b2b-29e4-455e-b172-42f329c3841c
# ╠═c07fd76f-1335-44da-9aa9-b44867baf7a8
# ╠═339c13d9-52b3-40bf-abbf-50332a355b41
