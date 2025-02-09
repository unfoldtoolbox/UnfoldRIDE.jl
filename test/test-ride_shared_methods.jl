using UnfoldRIDE
using Test
using Unfold
using UnfoldSim
using Statistics
using SignalAnalysis
using DataFrames
using DataFramesMeta
using Distributions
using Random
using Parameters
using CairoMakie
using DSP
using FFTW
using Peaks

include("../src/ride/ride_shared_methods.jl")

@testset "ride_shared_methods.jl" begin

    @testset "findxcorrpeak" begin
        d = UnfoldSim.pad_array(hanning(10),-35,0)
        kernel = hanning(20)

        f = Figure()
        lines(f[1,1],d)
        lines!(kernel)
        xc, m, onset = findxcorrpeak(d,kernel)
        lines(f[1,2],xc[1])

        lines(f[2,1],d)
        vlines!([m[1]])
        lines!(m[1].+(1:length(kernel)),kernel)

        display(f)
        using Test
        @test findxcorrpeak(d,kernel)[2] == [30]
        @test findxcorrpeak(d,kernel;window=true)[2] == [30]
    end

    ## test heuristic1
    @testset "heuristic1" begin
        latencies_df = DataFrame(latency = [75,33], fixed = [false, false])
        latencies_df_old = DataFrame(latency = [70,30], fixed = [false, false])
        latencies_df_old_old = DataFrame(latency = [50,40], fixed = [false, false])

        heuristic1_monoton_latency_changes!(latencies_df, latencies_df_old, latencies_df_old_old)

        @test latencies_df.latency == [75,30]
        @test latencies_df.fixed == [false, true]
    end

    ## test heuristic2
    @testset "heuristic2" begin
        epoch1 = reshape(vcat(zeros(100), hanning(100) .* -1), (1,:,1))
        epoch2 = reshape(vcat(zeros(100), hanning(100)), (1,:,1))
        data_epoched = cat(epoch1, epoch2, dims=3)
        
        latencies_df = DataFrame(latency = [0,100], fixed = [false, false])
        latencies_df_old = DataFrame(latency = [50,60], fixed = [false, false])

        xc, xc_values = findxcorrpeak(data_epoched[1,:,:], hanning(100))
        latencies_df_old = copy(latencies_df)
        latencies_df.latency = xc_values

        xc = [hanning(100) .* -1, hanning(100)]

        rng = MersenneTwister(12345)
        heuristic2_randommize_latency_on_convex_xcorr!(latencies_df, latencies_df_old, xc, rng)

        @test latencies_df.latency == [111,100]
        @test latencies_df.fixed == [true, false]
    end

    ## test heuristic3
    @testset "heuristic3" begin
        cfg = ride_config(
            sfreq = 100,
            c_range = [-0.4, 0.4],
            s_range = [0, 0],
            r_range = [0, 0],
            c_estimation_range = [0.2, 1.2],
            epoch_range = [-0.3,1.6],
            epoch_event_name = 'S',
        )
        #identical epochs with perfect match at 100 and subpar match at 300
        epoch1 = reshape(vcat(zeros(100), hanning(100), zeros(100), hanning(100).*0.9, zeros(100)), (1,:,1))
        epoch2 = reshape(vcat(zeros(100), hanning(100), zeros(100), hanning(100).*0.9, zeros(100)), (1,:,1))
        data_epoched = cat(epoch1, epoch2, dims=3)

        #same latency for both epochs
        latencies_df = DataFrame(latency = [100,100], fixed = [false, false])
        #201 is closer to the subpar 300 peak in the previous latencies, should trigger heuristic
        latencies_df_old = DataFrame(latency = [201,200], fixed = [false, false])
        xc, xc_values, onset = findxcorrpeak(data_epoched[1,:,:], hanning(100))

        heuristic3_pick_closest_xcorr_peak!(latencies_df, latencies_df_old, xc, 0.8, onset=onset)

        @test latencies_df.latency == [300, 100]
        @test latencies_df.fixed == [false, false]
    end

    #test filtering10
    @testset "filtering10" begin
        #sfreq = 100
        #data, evts = createTestData()
        #range_test = [0.0, 1.0]
        #data = reshape(data, (1,:))
        data = reshape(vcat(zeros(100), hanning(100) .* -1), (1,:))

        data_filtered = filtering10(data[1,:], 100, 0)

        f = Figure()
        lines(f[1,1], data_filtered[1,:], linewidth = 3)
        lines!(f[1,1], data[1,:], color = :red)
        display(f)
        
        #@test result_zero[1,:] == zeros(length(result_zero[1,:]))
    end

    #test dspfilter
    # this doesn't seem to work as expected. Does xcorr even make sense here?
    # how can I turn the xcorr result into a score 0-100 ?
    @testset "dspfilter" begin
        data = vcat(zeros(100), hanning(100), hanning(100) .* -1)

        noise = PinkNoise(;noiselevel = 0.1)
        data_noisy = copy(data)
        UnfoldSim.add_noise!(MersenneTwister(1234), noise, data_noisy)

        data_filtered = dspfilter(data_noisy, 5, 100)

        xcorr_result = xcorr(data, data_filtered; padmode = :none)  
        xcorr_max = findmax(xcorr_result)
        @test xcorr_max[1] > 50
        @test (xcorr_max[2] > 295 && xcorr_max[2] < 305)

        xcorr_terrible = xcorr(data, data_noisy; padmode = :none)
        xcorr_max_terrible = findmax(xcorr_terrible)
    end
end