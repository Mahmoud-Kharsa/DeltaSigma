module DeltaSigma

export
    calculateSNR,
    peakSNR,
    predictSNR,
    simulateDSM,
    simulateMS,
    simulateSNR,
    stuffABCD,
    synthesizeNTF

# Main functions
include("calculateSNR.jl")
include("peakSNR.jl")
include("predictSNR.jl")
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
