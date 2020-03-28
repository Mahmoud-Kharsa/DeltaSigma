using ControlSystems
using MAT

@testset "realizeNTF" begin
    # MATLAB output
    test = matread("resources/realizeNTF.mat")

    ntf = synthesizeNTF(5, 42, 1)
    forms = ["CRFB", "CRFF", "CRFBD", "CRFFD", "Stratos", "DSFB"]
    for i = 1:length(forms)
        form = forms[i]
        a, g, b, c = realizeNTF(ntf, form)
        @test isapprox(vec(test["a_$form"]), a)
        @test isapprox(vec(test["g_$form"]), g)
        @test isapprox(vec(test["b_$form"]), b)
        @test isapprox(vec(test["c_$form"]), c)
    end

    # "CIFB" and "CIFF" expect the real part of zeros to be 1
    ntf = zpk([1, 1+2im, 1-2im, 1+4im, 1-4im], [2, 2+1im, 2-1im, 2+3im, 2-3im], 1, 1)
    forms = ["CIFB", "CIFF"]
    for i = 1:length(forms)
        form = forms[i]
        a, g, b, c = realizeNTF(ntf, form)
        @test isapprox(vec(test["a_$form"]), a)
        @test isapprox(vec(test["g_$form"]), g)
        @test isapprox(vec(test["b_$form"]), b)
        @test isapprox(vec(test["c_$form"]), c)
    end

    # "PFF" expects zeros to be split into two groups separated by the 0.5
    # radian line in the complex plane
    ntf = zpk([0.5, 0.7+0.2im, 0.7-0.2im, 0.3+0.4im, 0.3-0.4im], [2, 2+1im, 2-1im, 2+3im, 2-3im], 1, 1)
    a, g, b, c = realizeNTF(ntf, "PFF")
    @test isapprox(vec(test["a_PFF"]), a)
    @test isapprox(vec(test["g_PFF"]), g)
    @test isapprox(vec(test["b_PFF"]), b)
    @test isapprox(vec(test["c_PFF"]), c)

    # test with stf
    ntf = synthesizeNTF(5, 42, 1)
    stf = tf(1, [1, 0, 1], 1)
    a, g, b, c = realizeNTF(ntf, "CRFB", stf)
    @test isapprox(vec(test["a_stf"]), a)
    @test isapprox(vec(test["g_stf"]), g)
    @test isapprox(vec(test["b_stf"]), b)
    @test isapprox(vec(test["c_stf"]), c)
end
