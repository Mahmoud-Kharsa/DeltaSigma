using MAT

@testset "simulateHBF" begin
    # MATLAB output
    test = matread("resources/simulateHBF.mat")

    x = [1; zeros(2^11-1)]
    for mode = 0:3
        for i = 1:4
            f1, f2 = exampleHBF(i)
            y = simulateHBF(x, f1, f2, mode)
            @test isapprox(test["y_$mode"][:,i], y)
        end
    end
end
