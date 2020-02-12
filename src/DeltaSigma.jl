module DeltaSigma

export 
    simulateSNR,
    synthesizeNTF

#Main functions
include("simulateSNR.jl")
include("synthesizeNTF.jl")

#Helper functions
include("createZPK.jl")
include("ds_optzeros.jl")

end # module
