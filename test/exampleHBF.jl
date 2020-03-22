using MAT

@testset "exampleHBF" begin
    # MATLAB output
    test = matread("resources/exampleHBF.mat")

    # HBFex1
    F1, F2 = exampleHBF(1)
    for i = 1:length(F1)
        @test test["F1_1"]["val"][i] == F1[1][i]
        @test test["F1_1"]["csd"][i] == F1[2][i]
    end
    for i = 1:length(F2)
        @test test["F2_1"]["val"][i] == F2[1][i]
        @test test["F2_1"]["csd"][i] == F2[2][i]
    end

    # HBFex2
    F1, F2 = exampleHBF(2)
    for i = 1:length(F1)
        @test test["F1_2"]["val"][i] == F1[1][i]
        @test test["F1_2"]["csd"][i] == F1[2][i]
    end
    for i = 1:length(F2)
        @test test["F2_2"]["val"][i] == F2[1][i]
        @test test["F2_2"]["csd"][i] == F2[2][i]
    end

    # Saramaki88
    F1, F2 = exampleHBF(3)
    for i = 1:length(F1)
        @test test["F1_3"]["val"][i] == F1[1][i]
        @test test["F1_3"]["csd"][i] == F1[2][i]
    end
    for i = 1:length(F2)
        @test test["F2_3"]["val"][i] == F2[1][i]
        @test test["F2_3"]["csd"][i] == F2[2][i]
    end

    # Saramaki90
    F1, F2 = exampleHBF(4)
    for i = 1:length(F1)
        @test test["F1_4"]["val"][i] == F1[1][i]
        @test test["F1_4"]["csd"][i] == F1[2][i]
    end
    for i = 1:length(F2)
        @test test["F2_4"]["val"][i] == F2[1][i]
        @test test["F2_4"]["csd"][i] == F2[2][i]
    end
end
