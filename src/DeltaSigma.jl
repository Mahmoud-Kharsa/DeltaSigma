module DeltaSigma

__precompile__(false)

export simulateSNR, synthesizeNTF

using ControlSystems

#Main functions
include("simulateSNR.jl")
include("synthesizeNTF.jl")

#Helper functions
include("createZPK.jl")
include("ds_optzeros.jl")

end # module
