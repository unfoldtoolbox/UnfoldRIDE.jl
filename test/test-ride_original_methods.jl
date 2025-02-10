include("../src/ride/ride_original_methods.jl")
using UnicodePlots

function createTestData()
    design = SingleSubjectDesign(;
        conditions = Dict(
            :condA => ["LevelA"],
        ),
    ) |> x -> RepeatDesign(x, 5);
    p1 = LinearModelComponent(;
        basis = hanning(100), 
        formula = @formula(0 ~ 1), 
        β = [1]
    );
    onset = UniformOnset(
        width = 0,
        offset = 200,
    );
    data, evts = simulate(
        MersenneTwister(1),
        design,
        [p1],
        onset,
    );
    return data, evts
end

@testset "ride_original_methods.jl" begin
    @testset "subtract_to_data1" begin
        sfreq = 100

        data = reshape(vcat(zeros(100), hanning(100), zeros(100), hanning(100), zeros(100)), (1,:))
        evts = DataFrame(:event => ['B','B'], :latency => [101,301])
        range_test = [0.0, 1.0]

        erp_to_subtract = reshape(hanning(100),(1,:,1))

        result_zero = subtract_to_data(data, [(evts, erp_to_subtract, range_test)], sfreq)

        display(lineplot(data[1,:],title = "data" ))
        display(lineplot(erp_to_subtract[1,:], title = "erp_to_subtract"))
        display(lineplot(result_zero[1,:], title = "result after subtraction (should be zeros)"))

        @test result_zero[1,:] == zeros(length(result_zero[1,:]))
    end

    @testset "subtract_to_data2" begin
        sfreq = 100
        data, evts = createTestData()
        range_test = [0.0, 1.0]
        data = reshape(data, (1,:))

        erp_to_subtract = reshape(hanning(100),(1,:,1))

        result_zero = subtract_to_data(data, [(evts, erp_to_subtract, range_test)], sfreq)

        f = Figure()
        lines(f[1,1], data[1,:])
        lines(f[1,2], erp_to_subtract[1,:])
        lines(f[2,1], result_zero[1,:])
        display(f)
        
        @test result_zero[1,:] == zeros(length(result_zero[1,:]))
    end
end