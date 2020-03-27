using MAT

@testset "find2dPIS" begin
    # MATLAB output
    test = matread("resources/find2dPIS.mat")

    A = [1 0; 1 1]
    B = [1 -1; 1 -2]
    C = [0 1]
    D = [0 0]
    ABCD = [A B; C D]
    u = 1 / pi

    s = find2dPIS(u, ABCD, 0)
    @test test["s"] == s
end
