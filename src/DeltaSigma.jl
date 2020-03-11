module DeltaSigma

export
    calculateSNR,
    dbv,
    simulateDSM,
    simulateSNR,
    synthesizeNTF

#Main functions
include("calculateSNR.jl")
include("dbv.jl")
include("simulateDSM.jl")
include("simulateSNR.jl")
include("synthesizeNTF.jl")

#Helper functions
include("createZPK.jl")
include("ds_optzeros.jl")

end # module
