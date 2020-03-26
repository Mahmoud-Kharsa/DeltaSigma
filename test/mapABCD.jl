using ControlSystems
using MAT

@testset "mapABCD" begin
    ntf = synthesizeNTF(5, 42, 1)
    forms = ["CRFB", "CRFF", "CRFBD", "CRFFD", "Stratos"]
    for i = 1:length(forms)
        form = forms[i]
        a1, g1, b1, c1 = realizeNTF(ntf, form)
        ABCD = stuffABCD(a1, g1, b1, c1, form)
        a2, g2, b2, c2 = mapABCD(ABCD, form)
        # test that mapABCD undoes stuffABCD, allow error in last digit
        @test isapprox(a1, a2, rtol=eps(Float64))
        @test isapprox(g1, g2, rtol=eps(Float64))
        @test isapprox(b1, b2, rtol=eps(Float64))
        @test isapprox(c1, c2, rtol=eps(Float64))
    end

    # "CIFB" and "CIFF" expect the real part of zeros to be 1
    ntf = zpk([1, 1+2im, 1-2im, 1+4im, 1-4im], [2, 2+1im, 2-1im, 2+3im, 2-3im], 1, 1)
    forms = ["CIFB", "CIFF"]
    for i = 1:length(forms)
        form = forms[i]
        a1, g1, b1, c1 = realizeNTF(ntf, form)
        ABCD = stuffABCD(a1, g1, b1, c1, form)
        a2, g2, b2, c2 = mapABCD(ABCD, form)
        # test that mapABCD undoes stuffABCD, allow error in last digit
        @test isapprox(a1, a2, rtol=eps(Float64))
        @test isapprox(g1, g2, rtol=eps(Float64))
        @test isapprox(b1, b2, rtol=eps(Float64))
        @test isapprox(c1, c2, rtol=eps(Float64))
    end

    # "PFF" and "DSFB" aren't supported by MATLAB toolbox yet
end
