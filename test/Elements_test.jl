using PSSFSS
using Test
using LinearAlgebra: norm

sh1 = rectstrip(Lx=1, Ly=1.0, Nx=1, Ny=1, Px=1, Py=1, units=inch)
sh2 = deepcopy(sh1)
@testset "recstrip" begin
    @test length(sh1.ρ) == 4
    @test sh1.ρ[2] - sh1.ρ[1] == [1, 0]
    @test sh1.ρ[3] - sh1.ρ[1] == [0, 1]
    @test sh1.ρ[4] - sh1.ρ[1] == [1, 1]
    @test sh1.e1 == [1, 3, 1, 2, 1]
    @test sh1.e2 == [2, 4, 3, 4, 4]
    @test sh1.fe == reshape([2, 3, 5, 4, 5, 1], 3, 2)
    @test sh1.fv == reshape([1, 4, 3, 1, 2, 4], 3, 2)
    @test sh1.fufp
    @test sh1.dx == sh1.dy == 0
    @test sh1.rot == 0
    @test sh1.units == inch
    @test sh1.class == 'J'
    @test sh1.s₁ == [1, 0]
    @test sh1.s₂ == [0, 1]
    @test sh1.β₁ ≈ 2π .* [1, 0]
    @test sh1.β₂ ≈ 2π .* [0, 1]
    @test sh1.Zs == 0
end

@testset "orient!" begin
    PSSFSS.Sheets.orient!(sh2, 45, [0,0])
    PSSFSS.Sheets.orient!(sh2, -45, [0,0])
    @test sh1.ρ ≈ sh2.ρ
end

Z = 0.2 + 0.3im
sh2 = rectstrip(Lx=1, Ly=1.0, Nx=1, Ny=1, Px=1, Py=1, units=inch, Zsheet=Z)
@testset "Zsheet" begin
    @test Z == sh2.Zs
end

R = 0.25
sh3 = rectstrip(Lx=1, Ly=1.0, Nx=1, Ny=1, Px=1, Py=1, units=inch, Rsheet=R)
@testset "Rsheet" begin
    @test R == sh3.Zs
end

@testset "badZsheet" begin
    @test_throws ErrorException rectstrip(Lx=1, Ly=1.0, Nx=1, Ny=1, Px=1, Py=1, units=inch, Zsheet=-0.2+0.2im)
end

@testset "StructuredPolyring" begin
    s = 0.4; ntri = 100
    s1 = [s, 0]; s2 = [0, s]; units = cm; orient = 45; a=[0.0]; b=[0.5*s/sqrt(2)]; sides = 4
    showprogress = false; logfile = resultfile = devnull

    sheet_structured = polyring(; ntri, units, s1, s2, sides, a, b, orient, structuredtri=true)
    sheet_unstructured = polyring(; ntri, units, s1, s2, sides, a, b, orient, structuredtri=false)

    FGHz = 30
    gsm_s = analyze([Layer(), sheet_structured, Layer()], FGHz, (ϕ=0, θ=0); showprogress, logfile, resultfile)[1].gsm
    gsm_u = analyze([Layer(), sheet_unstructured, Layer()], FGHz, (ϕ=0, θ=0); showprogress, logfile, resultfile)[1].gsm
    for i in 1:2, j in 1:2
        @test norm(gsm_s[i,j] - gsm_u[i,j], Inf) < 5e-3
    end
end

@testset "StructuredLoadedcross" begin
    s = 0.4; ntri = 400
    s1 = [s, 0]; s2 = [0, s]; units = cm; L1 = 0.8s; L2 = 0.3L1; w = 0.6L2
    showprogress = false; logfile = resultfile = devnull

    sheet_structured = loadedcross(; ntri, units, s1, s2, L1, L2, w, structuredtri=true)
    sheet_unstructured = loadedcross(; ntri, units, s1, s2, L1, L2, w, structuredtri=false)

    FGHz = 30
    gsm_s = analyze([Layer(), sheet_structured, Layer()], FGHz, (ϕ=0, θ=0); showprogress, logfile, resultfile)[1].gsm
    gsm_u = analyze([Layer(), sheet_unstructured, Layer()], FGHz, (ϕ=0, θ=0); showprogress, logfile, resultfile)[1].gsm
    for i in 1:2, j in 1:2
        @test norm(gsm_s[i,j] - gsm_u[i,j], Inf) < 5e-3
    end
end

@testset "StructuredJerusalemcross" begin
    s = 0.4; ntri = 600
    P = s; units = cm; L1 = 0.6s; L2 = 0.1L1; w = 0.6L2; A = 0.5L1; B = 0.12L1
    showprogress = false; logfile = resultfile = devnull

    sheet_structured = jerusalemcross(; ntri, units, A, B, P, L1, L2, w, structuredtri=true)
    sheet_unstructured = jerusalemcross(; ntri, units, A, B, P, L1, L2, w, structuredtri=false)

    FGHz = 30
    gsm_s = analyze([Layer(), sheet_structured, Layer()], FGHz, (ϕ=0, θ=0); showprogress, logfile, resultfile)[1].gsm
    gsm_u = analyze([Layer(), sheet_unstructured, Layer()], FGHz, (ϕ=0, θ=0); showprogress, logfile, resultfile)[1].gsm
    for i in 1:2, j in 1:2
        @test norm(gsm_s[i,j] - gsm_u[i,j], Inf) < 5e-3
    end
end

@testset "Manji" begin
    L1 = 0.481
    L2 = 0.22
    w = 0.06
    L3 = 0.1
    a = 0.15
    w2 = 0.04
    sheet = manji(; s1=[1.1, 0], s2=[0, 1.1], units=cm, L1, L2, L3, w, w2, a, ntri=1000)
    @test facecount(sheet) == 1634

    L3 = 0.16
    sheet = manji(; s1=[1.1, 0], s2=[0, 1.1], units=cm, L1, L2, L3, w, w2, a, ntri=1000)
    @test facecount(sheet) == 1448
end

@testset "Sinuous" begin
    P = 0.55; L2 = 0.95P
    s1 = P * [1, 0]; s2 = P * [0, 1]; sides = 45; units = cm
    ntri = 1400; w = w2 = 0.03; g = 0.02; rc = 0.05
    b=[0.1, 0.15, 0.2]
    sheet = sinuous(; arms=2, b, w, rc, g, w2, L2, sides, ntri, units, s1, s2)
    @test facecount(sheet) == 1525


    s1 = [1, 0]; s2 = [0, 1]
    b = [0.12, 0.2, 0.3]
    sides = 50; ntri = 2800; units = cm
    sheet = sinuous(; arms=4, b, w=0.03, rc=0.05, s1, s2,
                   L2=0.95, w2=0.03, c2=0.12, g=0.04, sides, ntri, units)
    @test facecount(sheet) > 3300 # Different on Apple, and depends on Julia version
end

@testset "ManjiFull" begin
    #=
    Angular Selective Surface
    From paper Zhenting Chen, Chao Du, Jie Liu, Di Zhou, and Zhongxiang Shen,
    "Design Methodology of Dual-Polarized Angle-Selective Surface Based
    on Three-Layer Frequency-Selective Surfaces", IEEE Trans. Antennas Propagat., Vol. 71, No. 11, November 2023,
    pp. 8704--8713.
    =#
    smat_true, smat_false = map([true, false]) do fufp
        P_paper = 50 # Unit cell dimension
        L_paper = 115.0
        L1_paper = 24.7
        L2_paper = 42.0
        L3_paper = 21.5
        L4_paper = 45.0
        W1_paper = 7.4
        W2_paper = 3.0
        kws = (; units = mm, class = 'J', w2 = 0.0, a = 0.0, s1 = [P_paper, 0], s2 = [0, P_paper], fufp, ntri=10)
        outer = manji(; L3=W1_paper, w=W1_paper, L1=L1_paper, L2=L2_paper/2-W1_paper, kws...)
        strata = [Layer(), outer, Layer()]
        flist = 3.0652
        steering = (θ=0, ϕ=0)
        result = analyze(strata, flist, steering; showprogress=false, logfile=devnull, resultfile=devnull)
        gsm = result[1].gsm
        smat = [gsm[1,1] gsm[1,2]; gsm[2,1] gsm[2,2]]
        return smat
    end
    @test norm(smat_true - smat_false) < 1e-12

end
