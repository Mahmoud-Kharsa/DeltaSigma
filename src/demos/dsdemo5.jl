using ControlSystems
using FFTW
using LaTeXStrings
using PyPlot
using Printf

function dsdemo5()
    show_usage = true
    osr = 25
    ntf = synthesizeNTF(6, osr, 1, 4)
    M = 16
    N = 2^14
    sigma_d = 0.01  # 1% mismatch
    window = ds_hann(N) / (M*N/8)
    windowM = repeat(window', outer=(M,1))

    mtf1 = zpk([1], [0], 1, 1)                  # First-order shaping
    mtf2 = zpk([1, 1], [0.3, 0.3], 1, 1)        # Second-order shaping
    mtf2b = synthesizeNTF(2, osr, 1, 2)         # Second-order shaping
    mtf4 = synthesizeNTF(4, osr*0.9, 1, 1.3)    # Fourth-order shaping

    f  = 0.01
    dw = []
    cases = [
        ( undbv(-3), nothing, 0.0, "Thermometer"        , "Thermometer"),
        ( undbv(-3),    mtf1, 0.0, "Rotation"           , "Rotation"),
        ( undbv(-3),    mtf2, 0.0, "2nd-order"          , string(L"2^{nd}", "-order")),
        (undbv(-30),    mtf1, 0.0, "Rotation"           , "Rotation"),
        (undbv(-30),    mtf1, 0.5, "Rot + dither"       , "Rot + dither"),
        ( undbv(-3),   mtf2b, 0.0, "2nd-order with zero", string(L"2^{nd}", "-order with zero")),
        ( undbv(-3),    mtf4, 0.0, "4th-order"          , string(L"4^{th}", "-order")),
    ]
    comparisons = [[1,2], [2,3], [4,5], [3,6], [6,7]]

    sv = Array{Matrix{Float64}}(undef, length(cases))
    leg = Array{String}(undef, length(cases))
    ltx_leg = Array{String}(undef, length(cases))
    Svv = Array{Vector{Float64}}(undef, length(cases))
    Sdd = Array{Vector{Float64}}(undef, length(cases))
    fin = zeros(Int, length(cases))
    for i = 1:length(cases)
        A = cases[i][1]
        mtf = cases[i][2]
        dither = cases[i][3]
        leg_i = cases[i][4]
        ltx_leg_i = cases[i][5]

        fin[i] = round(Int, f*N)
        inband = setdiff(2:(1+ceil(Int, 0.5*N/osr)), 1 .+ [0; 1; fin[i] .+ [-1, 0, 1]])
        w = (2*pi/N) * fin[i]
        u = M * A * sin.(w * (0:N-1))
        v, = simulateDSM(u, ntf, M+1)    # M unit elements requires an (M+1)-level quant.
        Svv[i] = abs.(fft(v .* window)).^2
        if isnothing(mtf)
            sv[i] = ds_therm(v, M)
        else
            sv[i], = simulateMS(v, M, mtf, dither, dw)
        end
        Sdd[i] = sigma_d^2 * sum(abs.(fft(sv[i] .* windowM, 2)).^2, dims=1)[:]
        mnp = sum(Sdd[i][inband]) / 1.5
        leg[i] = @sprintf("%s (MNP = %.0f dBFS)", leg_i, dbp(mnp))
        ltx_leg[i] = @sprintf("%s (MNP = %.0f dBFS)", ltx_leg_i, dbp(mnp))
    end

    # Plot results
    for comp_i = 1:length(comparisons)
        case_nums = comparisons[comp_i]
        println("Mismatch-Shaping Unit-Element DAC")
        @printf("Comparing %s vs. %s.\n", leg[case_nums[1]], leg[case_nums[2]])
        nc = length(case_nums)

        if show_usage
            fig = figure(1)
            fig.clf()
            fig.canvas.set_window_title("Element Usage")
            T = 25
            for i = 1:nc
                subplot(nc, 1, i)
                ci = case_nums[i]
                plotUsage(sv[ci][:,1:T])
            end
            fig.canvas.manager.window.raise_()
        end

        fig = figure(2)
        fig.clf()
        fig.canvas.set_window_title("Error Spectra")
        cols = ["b", "m", "r"]
        cleg = Array{String}(undef, nc)
        for i = 1:nc
            ci = case_nums[i]
            plotSpectrum(sqrt.(Sdd[ci]), fin[ci], cols[i], 8, 4)
            cleg[i] = ltx_leg[ci]
        end
        ci = case_nums[nc]
        plotSpectrum(sqrt.(Svv[ci]), fin[ci], "g", 8, 5)
        axis([1e-3, 0.5, -140, -50])
        plot([1e-3, 0.5/osr], -140*[1, 1], color="k", linewidth=4)
        text(0.5, -140, @sprintf("NBW=%.1e", 1.5/N), horizontalalignment="right", verticalalignment="bottom")
        grid()
        grid(which="minor", linestyle="dotted")
        ylabel("Error PSD")
        xlabel("Normalized Frequency")
        title(@sprintf("A = %.0fdBFS", dbp(Svv[ci][fin[ci]])))
        legend(cleg, loc="upper right")

        fig.canvas.manager.window.raise_()

        println("paused")
        readline()
    end

end
