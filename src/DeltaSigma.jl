module DeltaSigma

export
    calculateSNR,
    dbv,
    simulateDSM,
    simulateMS,
    simulateSNR,
    synthesizeNTF

#Main functions
include("calculateSNR.jl")
include("dbv.jl")
include("simulateDSM.jl")
include("simulateMS.jl")
include("simulateSNR.jl")
include("synthesizeNTF.jl")

#Helper functions
include("createZPK.jl")
include("ds_optzeros.jl")

end # module
