using ControlSystems
using DeltaSigma
using MAT

@testset "simulateSNR" begin
    vars = matread("simulateSNR.mat")
    OSR = 32
    H = zpk(vars["H_z"][:,1], vars["H_p"][:,1], 1, 1)
    snr_test = vars["snr"]
    amp_test = vars["amp"]
    
    snr, amp = simulateSNR(H, OSR)

end