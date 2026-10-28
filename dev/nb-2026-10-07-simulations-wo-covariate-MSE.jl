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

# ╔═╡ 5c55af8d-7baa-4dea-9a1c-2564fd5645ca
Pkg.add(url="https://github.com/unfoldtoolbox/UnfoldSim.jl",rev="feature-sequentialSamplingModels")

# ╔═╡ 6aee3f78-4c4d-4147-bf1b-1bec7859cbb7
Pkg.add([
        Pkg.PackageSpec(url="https://github.com/unfoldtoolbox/UnfoldRIDE.jl", rev="tukey+filtering"),
        ])

# ╔═╡ 94194cfa-720f-11f1-ae6e-1f0c5bc3f500
begin
	using UnfoldSim
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

# ╔═╡ 7eb581c3-2c2c-4400-86b3-f4696ee0fe1c
using Colors, ColorSchemes

# ╔═╡ 7b63ce13-93eb-4c5b-aabe-9a7d935f39b4
using Dates

# ╔═╡ e809e05c-3e18-4057-ae36-cc78331a5120
md"""
# TODO

- Check if latency ground truth is correct
   - UnfoldSim puts event marker at beginning (depending on offset), but RIDE looks for peak
- Get EffectsDesign \checkmark (sort of, still weird with diff. windows)
   - Why does padding lead to dim mismatch error depending on window size (rounding, but how to solve)
- Calculate RMSE for ERPs
"""

# ╔═╡ b1449cc6-f26f-44a3-8be9-99db8411e575
md"""
# Installation
"""

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

# ╔═╡ 035ccc9b-b1fa-4cc9-9ae0-b79af2bd4910
@bind s_amp PlutoUI.Select([0, 5, 10, 20, 50, 100]) 

# ╔═╡ 9379625a-ed23-44e2-a11d-4d2e8920bc1f
nl = 8 #@bind nl PlutoUI.Select([10, 5, 15, 2])

# ╔═╡ 23865862-aeb3-42f9-9119-e3d9991a76e6
begin
	noise = PinkNoise(; noiselevel = nl)
	#noise = NoNoise()
	noise = ExponentialNoise(; noiselevel = nl, τ = 4)
	#noise = IntrinsicNoise(10, sfreq, 10) # broken
end

# ╔═╡ adf9ea1e-7120-49f0-bf36-9787c7ebc80f
@bind rep PlutoUI.Select([5, 10, 20, 30, 50, 100])

# ╔═╡ 58be6593-096f-4371-b7dc-982691de7b4e
@bind hpf PlutoUI.Select([false, true])

# ╔═╡ 4c258133-d406-4db7-9c6f-426bbc63070a
@bind algo_filt PlutoUI.Select([(true, true),(true, false), (false, false), (false, true)])

# ╔═╡ 4715997f-2e4a-44fc-9dfd-6497a01e61b2
#(true, false), (false, false), (false, true),(true, true)

# ╔═╡ 36c369e8-1d81-4943-89ce-d243bff8174e
sfreq = 256

# ╔═╡ 3ea0e47a-ff79-4ba7-a5a2-66380ac1f69f
_, peak_offset = findmax(UnfoldSim.hanning(0.2, 0.1, sfreq)) # true peak of the C component; needed for GT latency calculation; set manually after inspection so far 

# ╔═╡ 657c79a5-38bb-41b6-9fc5-de3646bd906e
seeds = 1:50

# ╔═╡ 03d36b61-4d5a-489d-aacf-ca9594e417e6
md"""
## Design
"""

# ╔═╡ d25e53d4-459c-49f1-aa1b-837da10e302a
sequencedesign = RepeatDesign(SequenceDesign(SingleSubjectDesign(conditions=Dict( :condition=>["bike","face"])),"SCR"),rep)

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
# Simulate
"""

# ╔═╡ 8a397871-e9ad-4c12-a1bf-457d59e052d5
md"""
## Components
"""

# ╔═╡ 39161940-9458-44e0-83d9-89219adac04d
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
	    formula = @formula(0 ~ 1),
	    β = [c_amp],
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
	   # offset = -10*2.5,
	)
end

# ╔═╡ 13ea4788-6800-4f3b-95ef-12125c7c6f46
components = Dict('S' => [p1, n1], 'C' => [p3], 'R' => [resp])

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
    s_range = [-0.1, 0.25], # s locked
    r_range = [-0.1, 0.3], # r locked
    c_range = [-0.4, 0.4], # mean-c_n-1 latency
    tukey_window = [0.2, 0.5],
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

# ╔═╡ f16378da-9d1a-47fa-b3b3-468d812f6dee
@bind seed PlutoUI.Slider(1:42;show_value=true)

# ╔═╡ 60a569c6-72ec-42a5-8a05-0616bcf4868e
seq = UnfoldSim.SequenceOnset(Dict('S'=>LogNormalOnset(μ=4,σ=0.25),
									   'C'=>UniformOnset(0.1*sfreq,0.2*sfreq),
									   'R'=>UniformOnset(0.15*sfreq,0.3*sfreq)))

# ╔═╡ eb4b96bb-897a-48e9-85cf-c666c17a3386
	data, events =
	    simulate(MersenneTwister(seed), sequencedesign, components,seq, noise);

# ╔═╡ eb337f29-56c7-46a6-921e-78096c40fae5
begin
	# Plot continuous data
	idx = 3000:4000
	
	lines(idx, data[idx]) # plotting
	lats = filter(:latency => <(idx[end]), filter(:latency => >(idx[1]), events))
	#@show size(lats)
	vlines!(@rsubset(lats, :event == 'S')[:, :latency]; linestyle = :dash, label="S")
	vlines!(@rsubset(lats, :event == 'C')[:, :latency]; linestyle = :dash, label="C")
	vlines!(@rsubset(lats, :event == 'R')[:, :latency]; linestyle = :dash, label="R")
	axislegend(current_axis())
	current_figure()
end

# ╔═╡ c8439f1e-89d5-43de-bf28-5ff960aeb68a
begin
	evts = @rsubset(events, :event == 'S' || :event == 'R')
	N = size(@rsubset(events, :event == 'S')) 
	if hpf
		d = hp_filter(data, 0.5, sfreq)
	else
		d = data
	end
	
	disallowmissing!(evts)
	#run the ride algorithm
	#resultsClassic = ride_algorithm(ClassicMode, data_noisy, evts_without_c, cfg)
	resultsUnfold, UnfoldEvts, model = ride_algorithm(UnfoldMode, d, evts, cfg)

	resultsClassic, ClassicEvts = ride_algorithm(ClassicMode, d, evts, cfg)
end;

# ╔═╡ edfce652-a04c-4424-ab4e-b3dc34afe008
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

# ╔═╡ ea348fb8-1346-4885-b52d-669244954b49
let
	N = size(@rsubset(events, :event == 'S')) 
	print("Final number of epochs $N")
end

# ╔═╡ 510e2294-504d-4b8c-8081-6775a96df76a
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

# ╔═╡ d125abe0-a056-442a-939d-5aba8ee7d9c0
md"""
# Inspect marginalized effects
"""

# ╔═╡ f9ce46d4-076a-40fb-a299-6c7046514c2e
begin
	efdict = Dict()

	ef = effects(efdict, model)
end

# ╔═╡ be4b3fb1-bc20-451b-a21b-803d4ebd8bef
let
	plot_erp(@rsubset(ef, :eventname == 'C'), mapping = (; row = :eventname))
	lines!(current_axis(), range(start=0, stop=0.8, length=length(resultsClassic.c_erp_unpadded)),resultsClassic.c_erp_unpadded, label="C Classic RIDE", color = :green, linestyle = :dash)
	axislegend(current_axis())
	current_figure()
end

# ╔═╡ 6bf6730e-3a8d-47c7-a3f5-892c80b23a12
size(@rsubset(ef, :eventname == 'C'))

# ╔═╡ df8974ff-8e84-49fe-93b3-d2221ceb3183
md"""
# Get ground truth
Use EffectsDesign to get ground truth data
"""

# ╔═╡ 2e0f4230-f55b-48f5-be2b-63753e409f35
begin
    
    # Function to zero-pad ground_truth and get into correct format
    function format_gt(data, τ, sfreq, times; centering=false)

        #@show size(data)
        sp = length(range(start=τ[1], stop=τ[2], step=1/sfreq))
        @assert sp >= size(data, 1) "padded data shorter than ground truth data. Use bigger time window"
        padded_data = zeros(sp, size(data, 2))
        
    
        for col in 1:size(data, 2)
            N = size(data[:, col], 1)  # number of samples in signal

            before = length(range(start=τ[1], stop=0, step=1/sfreq))
            # Number of samples after the signal (from last sample to τ[2])
            after = length(range(start=0+1/sfreq, stop=τ[2], step=1/sfreq)) - N
            
            if centering # Center the components on half of the signal
                
                # Find peak (i.e. max) index
                _, idx = findmax(data[:, col])

         #       @show idx

                # Cut to center of peak
                before = before - idx # only works without offset
                
                # Add Überschuss
                after = after + idx

          #      @show abs(before)
           #     @show after
            end
           
            tmp = pad_array(
                reshape(data[:,col], size(data[:,col], 1)), 
                (-before, after),
                0) # pad ground truth to be same length as estimates

          #  @show sp
          #  @show size(tmp)
            padded_data[:, col] = tmp
        end
        padded_data = reshape(padded_data, 1, size(padded_data)...) # reshape to be channel x samplepoints x event
        
        return padded_data
    end
    
    # Changes effects design event DataFrame to be used with Unfold.result_to_table
    #df_to_vec(df) = [@rsubset(df, :event == e) for e in unique(df.event)] 
    
    function data_to_effects_table(data, events::DataFrame, τ::Tuple, sfreq::Int; centering=false)
          ## change gt to effects style DataFrame
          times = range(τ[1], stop = τ[2], step = 1 / sfreq) # get correct times vector
          g = format_gt(data, τ, sfreq, times, centering=centering) # zero-pad ground_truth and get into correct format
          effects = Unfold.result_to_table([g], [events], [times], ["effects"])
        return effects
    end
end

# ╔═╡ 4a7bd5c1-439a-4b8d-86d5-fba4fc845422
τ_vec = [
	'S' => Tuple(cfg.s_range), 
	'C' => Tuple(cfg.c_range), 
	'R' => Tuple(cfg.r_range)
];

# ╔═╡ 563ac9fa-cdbe-4761-87c1-8229bdc8b1ac
effects_design = EffectsDesign(sequencedesign, efdict)

# ╔═╡ 8d32b5e6-a2ee-47fb-af0c-90030c476150
gt_data, gt_events = simulate(
    MersenneTwister(1),
    effects_design,
    components,
    NoOnset(),
    NoNoise(),
    return_epoched = true,
);

# ╔═╡ 7c2b087a-206a-48c9-8667-2485a5cc1d8c
begin
	all_ef = DataFrame()
	for p in τ_vec[1:2]
		tmp_idx = findall(x -> x == p[1], gt_events.event)

		if p[1] == 'C'
			tmp_ef = data_to_effects_table(gt_data[:, tmp_idx], gt_events[tmp_idx, :], p[2], sfreq, centering=true)
		else
			tmp_ef = data_to_effects_table(gt_data[:, tmp_idx], gt_events[tmp_idx, :], p[2], sfreq)
		end
		all_ef = vcat(all_ef, tmp_ef, cols = :union)
	end
end

# ╔═╡ 6578832d-3110-4f10-8a54-8a68cedf2676
all_ef

# ╔═╡ eadb5c72-7acb-495b-9473-44dfc5129a73
begin
	#gt_effects.time = gt_effects.time ./ 100 .- 1 / 100
	plot_erp(@rsubset(all_ef, :event == 'C', :condition == "bike"), mapping=(;row=:event))
	plot_erp!(current_figure()[2,1],@rsubset(ef, :eventname == 'C'), mapping=(;row=:eventname))
	current_figure()
end

# ╔═╡ 1087f8d7-4df1-467c-98c1-2366cd4c6455
#c_gt_vec = @rsubset(all_ef, :event == 'C').yhat

# ╔═╡ f8fb7123-ab67-4001-abae-464e21136d66
#c_vec = @rsubset(ef, :eventname == 'C').yhat

# ╔═╡ b4a17a71-5f18-45b6-9ed8-24601d5ee211
begin # Calculate RMSE over all cont predictors
	# Get classic vector once
	c_vec_cl = (resultsClassic.c_erp_unpadded)[1:end-1]

	# Init rmse vectors
	all_rmse_uf = []
	all_rmse_cl = []
	
	# Loop over predictors
	# Get ground truth vector
	c_gt_vec = @rsubset(all_ef, :event == 'C', :condition == "bike").yhat

	# Get unfold results vector
	c_vec = @rsubset(ef, :eventname == 'C').yhat[1:end-1]

	# Calc. RMSE
	all_rmse_c_uf = sqrt(mean((c_gt_vec - c_vec).^2))
	all_rmse_c_cl = sqrt(mean((c_gt_vec - c_vec_cl).^2))
	
end

# ╔═╡ b7ff4859-ad0c-4d09-8406-c99f96351c2a
all_rmse_c_uf

# ╔═╡ 71cd87e2-663b-46a0-aef8-8b56bfb42f50
md"""
## Over seeds
"""

# ╔═╡ 1c383a6e-ac72-48e2-862c-fc7dd9d5c3a0
begin
	rmse_cl = []
	rmse_uf = []
	for s in seeds

		dat, gt_evnts =
	    simulate(MersenneTwister(s), sequencedesign, components,seq, noise);
		
		evnts = @rsubset(gt_evnts, :event == 'S' || :event == 'R')


		if hpf
			d2 = hp_filter(dat, 0.5, sfreq)
		else
			d2 = dat
		end
		
		disallowmissing!(evnts)
		#run the ride algorithm
		#resultsClassic = ride_algorithm(ClassicMode, data_noisy, evts_without_c, cfg)
		resultsUF, UFEvts, m = ride_algorithm(UnfoldMode, d2, evnts, cfg)
	
		resultsCL, CLEvts = ride_algorithm(ClassicMode, d2, evnts, cfg)

		# Get effects for uf results
		tmpefdict = Dict()

		tmp_ef = effects(efdict, m)
		
		#tmp_unique_cont = unique(all_ef.continuous) # get all unique predictors

		# Get classic vector once
		tmp_c_vec_cl = (resultsCL.c_erp_unpadded)[1:end-1]
	
		# Init rmse vectors
		tmp_all_rmse_uf = []
		tmp_all_rmse_cl = []
		
			# Loop over predictors
		# Get ground truth vector
		tmp_c_gt_vec = @rsubset(all_ef, :event == 'C', :condition == "bike").yhat
	
		# Get unfold results vector
		tmp_c_vec = @rsubset(tmp_ef, :eventname == 'C').yhat[1:end-1]
	
		# Calc. RMSE
		tmp_all_rmse_c_uf = sqrt(mean((tmp_c_gt_vec - tmp_c_vec).^2))
		tmp_all_rmse_c_cl = sqrt(mean((tmp_c_gt_vec - tmp_c_vec_cl).^2))
		
		push!(rmse_uf, tmp_all_rmse_c_uf)
		push!(rmse_cl, tmp_all_rmse_c_cl)
	end
	
end

# ╔═╡ c1388558-ed4b-4107-a91d-8dce734135db
md"""
## Plotting
"""

# ╔═╡ d164d5ed-48fe-4908-9a73-9ebd614987e8
let
	scatter(rmse_uf, label="Unfold")
	scatter!(current_axis(), rmse_cl, label="Classic")
	#limits!(current_axis(), nothing, nothing, 2.5, 2.6)
	axislegend(current_axis())
	current_figure()

end

# ╔═╡ c4a2db16-d3ce-4a38-a898-988333fab443
let
    categories = ["Classic", "Unfold"]
    colors = ColorScheme(range(colorant"red", colorant"blue", length=length(categories)))
    
    fig = Figure(size = (600, 400); title="RMSE for both algorithms over $(length(seeds)) seeds")
    ax = Axis(fig[1, 1], xticks = (1:length(categories), categories), 
             # limits = (nothing, (2, 3)),
             title="RMSE for both algorithms over $(length(seeds)) seeds",
              ylabel = "Mean RMSE over entire epoch (μV)"
             )
    all =    [rmse_cl, rmse_uf] 
    for (indx, f) in enumerate(categories)
        datam = all[indx]
        a = fill(indx, length(datam))
        boxplot!(ax, a, datam; whiskerwidth = 1, width = .45,
            color = (colors[indx], 0.45), whiskercolor = (colors[indx], 1),
            mediancolor = :black)
    end
    fig

end

# ╔═╡ cd296375-a0cb-4c44-89d8-4384e1e0af5b
md"""
## Misc
"""

# ╔═╡ 715f2c9f-4719-45a8-9813-7fbe17a0b5d6
resultsClassic.c_latencies

# ╔═╡ 7262b007-e52a-4de8-8d73-47cc9f7c03b7
resultsUnfold[1].c_latencies

# ╔═╡ a951850d-c126-4268-a8b9-5d11f0eb93e7
let
	lines(resultsClassic.c_latencies, label="Classic")
	lines!(current_axis(), resultsUnfold[1].c_latencies, label="Unfold")
	axislegend(current_axis())
	current_figure()
end

# ╔═╡ 30401e4a-092b-4bdf-847d-408ec39bf3ab
md"""
# Latency ground truth
"""

# ╔═╡ 58de58cf-26ea-476a-bc7e-b503dc9e6a19
md"""
## C latencies
"""

# ╔═╡ c7c56a5d-8a24-47e2-88da-897619f463a9
begin
	lat_rmse_cl = []
	lat_rmse_uf = []
	for s in seeds

		dat, gt_evnts =
	    simulate(MersenneTwister(s), sequencedesign, components,seq, noise);
		
		evnts = @rsubset(gt_evnts, :event == 'S' || :event == 'R')


		if hpf
			d2 = hp_filter(dat, 0.5, sfreq)
		else
			d2 = dat
		end
		
		disallowmissing!(evnts)
		#run the ride algorithm
		#resultsClassic = ride_algorithm(ClassicMode, data_noisy, evts_without_c, cfg)
		resultsUF, UFEvts, m = ride_algorithm(UnfoldMode, d2, evnts, cfg)
	
		resultsCL, CLEvts = ride_algorithm(ClassicMode, d2, evnts, cfg)

		
		UF = (@rsubset(gt_evnts, :event == 'C').latency[1:end-1] .+ peak_offset) - @rsubset(UFEvts, :event == 'C').latency
		CL = (@rsubset(gt_evnts, :event == 'C').latency[1:end-1] .+ peak_offset) - @rsubset(CLEvts, :event == 'C').latency

		push!(lat_rmse_uf, sqrt(mean(UF.^2)))
		push!(lat_rmse_cl, sqrt(mean(CL.^2)))
	end

end;

# ╔═╡ a475da9e-149f-4897-9eb6-f83bc4398654
md"""
## Variance over seeds
"""

# ╔═╡ 181c458e-baa3-4978-81a7-6919d39aaaaa
let
    categories = ["Classic", "Unfold"]
    colors = ColorScheme(range(colorant"red", colorant"blue", length=length(categories)))
    
    fig = Figure(size = (600, 400))
    ax = Axis(fig[1, 1], xticks = (1:length(categories), categories),
             title="RSME difference in latencies over $(length(seeds)) seeds",
              ylabel = "RMSE of latencies (samples)"
             )
    all =    [lat_rmse_cl, lat_rmse_uf] 
    for (indx, f) in enumerate(categories)
        datam = all[indx]
        a = fill(indx, length(datam))
        boxplot!(ax, a, datam; whiskerwidth = 1, width = .45,
            color = (colors[indx], 0.45), whiskercolor = (colors[indx], 1),
            mediancolor = :black)
    end
    fig

end

# ╔═╡ 4c54e883-d54b-48c1-bc7b-4e1e9120db11
md"""
# Save
"""

# ╔═╡ 3635ae13-69ba-4b0f-a9f0-00dfb20150a6
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

# ╔═╡ f447147c-9b6c-4a07-8163-f18ba5c24cb0
savepath = "/store/users/skukies/Projects/2026UnfoldRIDE/hdf-data"

# ╔═╡ 809b2c9f-26b7-493c-a7ea-cb5dc4140a72
bname = Dates.format(today(), "yyyy-mm-dd") * "-Simulations-S_amp_" * string(s_amp) * "C_amp_" * string(c_amp) * "R_amp_" * string(r_amp) * "-NL_"* string(nl) * "-TN_" * string(rep*2)

# ╔═╡ dc05150b-ca4a-4501-901f-7043c2883afb
@bind safe PlutoUI.Select([false, true])

# ╔═╡ 2631afb9-2370-4988-a6ab-eb8e3d6d2003
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

# ╔═╡ Cell order:
# ╟─e809e05c-3e18-4057-ae36-cc78331a5120
# ╟─b1449cc6-f26f-44a3-8be9-99db8411e575
# ╠═b43c1614-3d50-42ae-8b98-fec5980b8d66
# ╠═5c55af8d-7baa-4dea-9a1c-2564fd5645ca
# ╠═6aee3f78-4c4d-4147-bf1b-1bec7859cbb7
# ╠═94194cfa-720f-11f1-ae6e-1f0c5bc3f500
# ╠═59987f84-1a4c-418b-9be3-137ab9720415
# ╟─fb5e13d7-bd31-42e7-b09d-a3a93f5cb46b
# ╠═23865862-aeb3-42f9-9119-e3d9991a76e6
# ╠═4e22effa-b136-4b60-b3f2-09f3b66d2970
# ╠═26c27b46-9b8a-4536-843c-597701ec9e47
# ╠═968c8b50-2a29-4359-b5c4-6347402ee055
# ╠═035ccc9b-b1fa-4cc9-9ae0-b79af2bd4910
# ╠═9379625a-ed23-44e2-a11d-4d2e8920bc1f
# ╠═adf9ea1e-7120-49f0-bf36-9787c7ebc80f
# ╠═58be6593-096f-4371-b7dc-982691de7b4e
# ╠═4c258133-d406-4db7-9c6f-426bbc63070a
# ╠═4715997f-2e4a-44fc-9dfd-6497a01e61b2
# ╠═36c369e8-1d81-4943-89ce-d243bff8174e
# ╠═3ea0e47a-ff79-4ba7-a5a2-66380ac1f69f
# ╠═657c79a5-38bb-41b6-9fc5-de3646bd906e
# ╟─03d36b61-4d5a-489d-aacf-ca9594e417e6
# ╠═d25e53d4-459c-49f1-aa1b-837da10e302a
# ╟─f5a019fa-00ef-4398-8ce1-979e1fc509ee
# ╠═232dfd69-4877-438c-b4c5-3e3e8b3ed9ff
# ╟─758ff988-68f9-4e28-b077-0a5c84b53230
# ╟─8a397871-e9ad-4c12-a1bf-457d59e052d5
# ╠═39161940-9458-44e0-83d9-89219adac04d
# ╠═13ea4788-6800-4f3b-95ef-12125c7c6f46
# ╟─28ace3a7-136b-43ac-ac95-748f6a17b1c5
# ╠═eb4b96bb-897a-48e9-85cf-c666c17a3386
# ╠═eb337f29-56c7-46a6-921e-78096c40fae5
# ╟─b017f0d0-b28b-4f03-916d-b0b4d763da82
# ╟─440fefab-69bf-4ae5-8f5c-c1cfa880387c
# ╠═5f4758d0-6a09-4719-964c-385d19dcad42
# ╟─a93faa72-819c-457b-9a51-bc9d430e89c4
# ╟─04152ac0-69ae-4a98-900c-9a7c82c635c9
# ╟─c229cfcd-1522-4fcc-93d2-9f2e04eaba92
# ╠═c8439f1e-89d5-43de-bf28-5ff960aeb68a
# ╠═ea348fb8-1346-4885-b52d-669244954b49
# ╟─247d06f3-2389-4ab4-9fcd-de7dd537f78d
# ╠═edfce652-a04c-4424-ab4e-b3dc34afe008
# ╠═f16378da-9d1a-47fa-b3b3-468d812f6dee
# ╠═60a569c6-72ec-42a5-8a05-0616bcf4868e
# ╠═510e2294-504d-4b8c-8081-6775a96df76a
# ╟─d125abe0-a056-442a-939d-5aba8ee7d9c0
# ╠═f9ce46d4-076a-40fb-a299-6c7046514c2e
# ╠═be4b3fb1-bc20-451b-a21b-803d4ebd8bef
# ╠═6bf6730e-3a8d-47c7-a3f5-892c80b23a12
# ╟─df8974ff-8e84-49fe-93b3-d2221ceb3183
# ╟─2e0f4230-f55b-48f5-be2b-63753e409f35
# ╠═8d32b5e6-a2ee-47fb-af0c-90030c476150
# ╠═4a7bd5c1-439a-4b8d-86d5-fba4fc845422
# ╠═7c2b087a-206a-48c9-8667-2485a5cc1d8c
# ╠═6578832d-3110-4f10-8a54-8a68cedf2676
# ╠═563ac9fa-cdbe-4761-87c1-8229bdc8b1ac
# ╠═eadb5c72-7acb-495b-9473-44dfc5129a73
# ╠═1087f8d7-4df1-467c-98c1-2366cd4c6455
# ╠═f8fb7123-ab67-4001-abae-464e21136d66
# ╠═b4a17a71-5f18-45b6-9ed8-24601d5ee211
# ╠═b7ff4859-ad0c-4d09-8406-c99f96351c2a
# ╟─71cd87e2-663b-46a0-aef8-8b56bfb42f50
# ╠═1c383a6e-ac72-48e2-862c-fc7dd9d5c3a0
# ╟─c1388558-ed4b-4107-a91d-8dce734135db
# ╠═d164d5ed-48fe-4908-9a73-9ebd614987e8
# ╠═c4a2db16-d3ce-4a38-a898-988333fab443
# ╟─cd296375-a0cb-4c44-89d8-4384e1e0af5b
# ╠═715f2c9f-4719-45a8-9813-7fbe17a0b5d6
# ╠═7262b007-e52a-4de8-8d73-47cc9f7c03b7
# ╠═a951850d-c126-4268-a8b9-5d11f0eb93e7
# ╟─30401e4a-092b-4bdf-847d-408ec39bf3ab
# ╟─58de58cf-26ea-476a-bc7e-b503dc9e6a19
# ╠═c7c56a5d-8a24-47e2-88da-897619f463a9
# ╠═7eb581c3-2c2c-4400-86b3-f4696ee0fe1c
# ╟─a475da9e-149f-4897-9eb6-f83bc4398654
# ╠═181c458e-baa3-4978-81a7-6919d39aaaaa
# ╟─4c54e883-d54b-48c1-bc7b-4e1e9120db11
# ╠═7b63ce13-93eb-4c5b-aabe-9a7d935f39b4
# ╠═3635ae13-69ba-4b0f-a9f0-00dfb20150a6
# ╠═f447147c-9b6c-4a07-8163-f18ba5c24cb0
# ╠═809b2c9f-26b7-493c-a7ea-cb5dc4140a72
# ╠═dc05150b-ca4a-4501-901f-7043c2883afb
# ╠═2631afb9-2370-4988-a6ab-eb8e3d6d2003
