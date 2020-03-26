module DeltaSigma

export
    calculateSNR,
    mapABCD,
    peakSNR,
    predictSNR,
    realizeNTF,
    simulateDSM,
    simulateMS,
    simulateSNR,
    stuffABCD,
    synthesizeNTF

# Main functions
include("calculateSNR.jl")
include("mapABCD.jl")
include("peakSNR.jl")
include("predictSNR.jl")
include("realizeNTF.jl")
include("simulateDSM.jl")
include("simulateMS.jl")
include("simulateSNR.jl")
include("stuffABCD.jl")
include("synthesizeNTF.jl")

# Helper functions
include("calculateTF.jl")
include("ds_optzeros.jl")
include("partitionABCD.jl")
include("selectElement.jl")

end # module
