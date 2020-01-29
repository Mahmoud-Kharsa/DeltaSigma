module DeltaSigma

export
    simulateDSM,
    simulateSNR,
    synthesizeNTF

#Main functions
include("simulateDSM.jl")
include("simulateSNR.jl")
include("synthesizeNTF.jl")

#Helper functions
include("createZPK.jl")
include("ds_optzeros.jl")

end # module
