module DeltaSigma

export
    calculateSNR,
    calculateTF,
    dbp,
    dbv,
    designHBF,
    ds_hann,
    dsdemo1,
    dsdemo2,
    dsdemo6,
    dsdemo7,
    exampleHBF,
    find2dPIS,
    mapABCD,
    partitionABCD,
    peakSNR,
    predictSNR,
    realizeNTF,
    rmsGain,
    scaleABCD,
    simulateDSM,
    simulateHBF,
    simulateMS,
    simulateSNR,
    stuffABCD,
    synthesizeNTF,
    undbp,
    undbv

# Key functions
include("calculateTF.jl")
include("mapABCD.jl")
include("realizeNTF.jl")
include("scaleABCD.jl")
include("simulateDSM.jl")
include("simulateMS.jl")
include("simulateSNR.jl")
include("stuffABCD.jl")
include("synthesizeNTF.jl")

# Speciality Functions
include("designHBF.jl")
include("exampleHBF.jl")
include("predictSNR.jl")
include("simulateHBF.jl")
include("PosInvSet/find2dPIS.jl")

# Utility functions
include("calculateSNR.jl")
include("dbp.jl")
include("dbv.jl")
include("ds_hann.jl")
include("partitionABCD.jl")
include("peakSNR.jl")
include("rmsGain.jl")
include("undbp.jl")
include("undbv.jl")

# Demos
include("demos/dsdemo1.jl")
include("demos/dsdemo2.jl")
include("demos/dsdemo6.jl")
include("demos/dsdemo7.jl")

# Internal helper functions
include("bquantize.jl")
include("bunquantize.jl")
include("delay.jl")
include("ds_optzeros.jl")
include("evalF0.jl")
include("evalF1.jl")
include("evalRPoly.jl")
include("evalTF.jl")
include("frespHBF.jl")
include("selectElement.jl")
include("PosInvSet/dscut.jl")
include("PosInvSet/dsmap.jl")
include("PosInvSet/dssplit2d.jl")
include("PosInvSet/hull2d.jl")
include("PosInvSet/leftof.jl")
include("PosInvSet/outconvex2d.jl")

# Internal graphing functions
include("axisLabels.jl")
include("figureMagic.jl")
include("plotPZ.jl")
include("PosInvSet/dotplot.jl")
include("PosInvSet/polyplot.jl")

end # module
