using PSSFSS, Test
using LinearAlgebra: I
@testset "sympixels" begin
    P = 5
    units = mm
    nrim = 0
    halfnint = 2
    pdiv = 2

    veclen = halfnint * (halfnint + 1) ÷ 2
    patternvec = ones(Bool, veclen)

    sheet = sympixels(; P, nrim, halfnint, patternvec, units, pdiv, class='J')
    flist = 1
    steering = (θ=0, ϕ=0)
    logfile = resultfile = devnull
    showprogress = false
    res = analyze([Layer(), sheet, Layer()], flist, steering; showprogress, logfile, resultfile) |> only
    gsm4x4 = [res.gsm[1,1] res.gsm[1,2]; res.gsm[2,1] res.gsm[2,2]]
    @test maximum(abs, gsm4x4 + I) < 1e-5
end