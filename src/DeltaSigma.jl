module DeltaSigma

export
    calculateSNR,
    simulateDSM,
    simulateMS,
    simulateSNR,
    synthesizeNTF

# Main functions
include("calculateSNR.jl")
include("simulateDSM.jl")
include("simulateMS.jl")
include("simulateSNR.jl")
include("synthesizeNTF.jl")

# Helper functions
include("ds_optzeros.jl")
include("selectElement.jl")

end # module
