using ControlSystems
using DeltaSigma
using MAT

@testset "simulateDSM" begin
    vars = matread("simulateDSM.mat")
    u = vars["u"]
    H = zpk(vars["H_z"][:,1], vars["H_p"][:,1], 1, 1)
    v_test = vars["v"]
    y_test = vars["y"]
    
    v, _, _, y = simulateDSM(u, H)
    @test v == v_test
    @test isapprox(y, y_test)
end