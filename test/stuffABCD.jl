using ControlSystems
using MAT

@testset "stuffABCD" begin
    # MATLAB output
    test = matread("resources/stuffABCD.mat")

    ntf = synthesizeNTF(5, 42, 1)
    forms = ["CRFB", "CRFF", "CRFBD", "CRFFD", "Stratos", "DSFB"]
    for i = 1:length(forms)
        form = forms[i]
        a, g, b, c = realizeNTF(ntf, form)
        ABCD = stuffABCD(a, g, b, c, form)
        @test isapprox(test["ABCD_$form"], ABCD)
    end

    # "CIFB" and "CIFF" expect the real part of zeros to be 1
    ntf = zpk([1, 1+2im, 1-2im, 1+4im, 1-4im], [2, 2+1im, 2-1im, 2+3im, 2-3im], 1, 1)
    forms = ["CIFB", "CIFF"]
    for i = 1:length(forms)
        form = forms[i]
        a, g, b, c = realizeNTF(ntf, form)
        ABCD = stuffABCD(a, g, b, c, form)
        @test isapprox(test["ABCD_$form"], ABCD)
    end

    # "PFF" expects zeros to be split into two groups separated by the 0.5
    # radian line in the complex plane
    ntf = zpk([0.5, 0.7+0.2im, 0.7-0.2im, 0.3+0.4im, 0.3-0.4im], [2, 2+1im, 2-1im, 2+3im, 2-3im], 1, 1)
    a, g, b, c = realizeNTF(ntf, "PFF")
    ABCD = stuffABCD(a, g, b, c, "PFF")
    @test isapprox(test["ABCD_PFF"], ABCD)
end
