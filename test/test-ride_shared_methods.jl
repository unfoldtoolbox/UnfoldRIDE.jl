include("../src/ride/ride_shared_methods.jl")

@testset "ride_shared_methods.jl" begin

    @testset "findxcorrpeak" begin
        d = UnfoldSim.pad_array(hanning(10), -35, 0)
        kernel = hanning(20)


        f = lineplot(d, title = "data + kernel")
        lineplot!(f, kernel)
        display(f)
        xc, m, onset = findxcorrpeak(d, kernel)
        display(lineplot(xc[1], title = "xcorr"))

        f = lineplot(d, title = "data + kernel corrected")
        lineplot!(f, vcat(zeros(m[1]), kernel))
        display(f)

        using Test
        @test findxcorrpeak(d, kernel)[2] == [30]
        @test findxcorrpeak(d, kernel; window = true)[2] == [30]
    end

    ## test heuristic1
    @testset "heuristic1" begin
        latencies_df = DataFrame(latency = [75, 33], fixed = [false, false])
        latencies_df_old = DataFrame(latency = [70, 30], fixed = [false, false])
        latencies_df_old_old = DataFrame(latency = [50, 40], fixed = [false, false])

        heuristic1_monoton_latency_changes!(
            latencies_df,
            latencies_df_old,
            latencies_df_old_old,
        )

        @test latencies_df.latency == [75, 30]
        @test latencies_df.fixed == [false, true]
    end

    ## test heuristic2
    @testset "heuristic2" begin
        epoch1 = reshape(vcat(zeros(100), hanning(100) .* -1), (1, :, 1))
        epoch2 = reshape(vcat(zeros(100), hanning(100)), (1, :, 1))
        data_epoched = cat(epoch1, epoch2, dims = 3)

        latencies_df = DataFrame(latency = [0, 100], fixed = [false, false])
        latencies_df_old = DataFrame(latency = [50, 60], fixed = [false, false])

        xc, xc_values = findxcorrpeak(data_epoched[1, :, :], hanning(100))
        latencies_df_old = copy(latencies_df)
        latencies_df.latency = copy(xc_values)

        xc = [hanning(100) .* -1, hanning(100)]

        rng = MersenneTwister(12345)
        heuristic2_randomize_latency_on_convex_xcorr!(
            latencies_df,
            latencies_df_old,
            xc,
            rng,
        )

        #check that latencies_df.latency was randomized during the heuristic
        @test latencies_df.latency != xc_values
        #after randomization, the latency should be fixed
        @test latencies_df.fixed == [true, false]
    end

    ## test heuristic3
    @testset "heuristic3" begin
        cfg = RideConfig(
            sfreq = 100,
            c_range = [-0.4, 0.4],
            s_range = [0, 0],
            r_range = [0, 0],
            c_estimation_range = [0.2, 1.2],
            epoch_range = [-0.3, 1.6],
        )
        #identical epochs with perfect match at 100 and subpar match at 300
        epoch1 = reshape(
            vcat(zeros(100), hanning(100), zeros(100), hanning(100) .* 0.9, zeros(100)),
            (1, :, 1),
        )
        epoch2 = reshape(
            vcat(zeros(100), hanning(100), zeros(100), hanning(100) .* 0.9, zeros(100)),
            (1, :, 1),
        )
        data_epoched = cat(epoch1, epoch2, dims = 3)

        #same latency for both epochs
        latencies_df = DataFrame(latency = [100, 100], fixed = [false, false])
        #201 is closer to the subpar 300 peak in the previous latencies, should trigger heuristic
        latencies_df_old = DataFrame(latency = [201, 200], fixed = [false, false])
        xc, xc_values, onset = findxcorrpeak(data_epoched[1, :, :], hanning(100))

        heuristic3_pick_closest_xcorr_peak!(
            latencies_df,
            latencies_df_old,
            xc;
            equality_threshold = 0.8,
            onset = onset,
        )

        @test latencies_df.latency == [300, 100]
        @test latencies_df.fixed == [false, false]
    end

    #test dspfilter
    @testset "dspfilter" begin
        xcorr_max = []
        xcorr_max_noisy = []
        for i = 1:10000
            data = vcat(zeros(100), hanning(100), hanning(100) .* -1)

            noise = PinkNoise(; noiselevel = 0.1)
            data_noisy = copy(data)
            UnfoldSim.add_noise!(MersenneTwister(i), noise, data_noisy)

            data_filtered = dspfilter(data_noisy, 5, 100)

            xcorr_result = xcorr(normalize(data), normalize(data_filtered); padmode = :none)
            push!(xcorr_max, findmax(xcorr_result)[1])

            xcorr_noisy = xcorr(normalize(data), normalize(data_noisy); padmode = :none)
            push!(xcorr_max_noisy, findmax(xcorr_noisy)[1])
        end
        @test !(false in (xcorr_max .> 0.9))
        @test !(0 in (xcorr_max_noisy .< xcorr_max))
    end
end