using PSSFSS
using LinearAlgebra: norm, diagind
using PSSFSS.GSMs: GSM, cascade, cascade!, initialize_gsm_file, 
                  append_gsm_data, read_gsm_file, choose_gblocks
using PSSFSS.Substrate: Layer, TEorTM, TE, TM
using PSSFSS.PSSFSSLen
using PSSFSS.Substrate: Layer
using PSSFSS.Elements: rectstrip
using PSSFSS.Constants: c₀, twopi
using Test


atol = 1e-16 # For comparison to zero


@testset "cascade1" begin
    na1 = 10
    na2 = nb1 = 20
    nb2 = 15
    
    a = GSM(na1, na2)
    b = GSM(nb1, nb2)
    c = cascade(a, b)

    @test norm(c.s11) == norm(c.s22) == 0
    @test all(c.s12[diagind(c.s12)] .== 1)
    @test all(c.s12[setdiff(LinearIndices(c.s12), diagind(c.s12))] .== 0)
    @test all(c.s21[diagind(c.s21)] .== 1)
    @test all(c.s21[setdiff(LinearIndices(c.s21), diagind(c.s21))] .== 0)
end

@testset "cascade2" begin
    na1 = 1
    na2 = nb1 = 1
    nb2 = 1
    α = 0.35
    β = sqrt(1 - α^2)
    a = GSM(na1, na2)
    a.s11 .= a.s22 .= α
    a.s12 .= a.s21 .= β

    b = GSM(nb1, nb2)
    γ = 0.23
    δ = sqrt(1 - γ^2)
    a = GSM(na1, na2)
    b.s11 .= b.s22 .= γ
    b.s12 .= b.s21 .= δ
    c = cascade(a, b)
    Δ = 1 / (1 - a.s22[1,1] * b.s11[1,1])
    s11 = a.s11 + a.s12 .* a.s21 .* b.s11 ./ Δ
    s12 = a.s12 .* b.s12 ./ Δ
    s21 = a.s21 .* b.s21 ./ Δ
    s22 = b.s22 + b.s12 .* b.s21 .* a.s22 ./ Δ
    @test norm(c.s11 - s11) ≤ atol
    @test norm(c.s12 - s12) ≤ atol
    @test norm(c.s21 - s21) ≤ atol
    @test norm(c.s22 - s22) ≤ atol
end

@testset "cascade!" begin
    na1 = 1
    na2 = nb1 = 1
    nb2 = 1
    α = 0.35
    β = sqrt(1 - α^2)
    a = GSM(na1, na2)
    a.s11 .= a.s22 .= α
    a.s12 .= a.s21 .= β
    b = deepcopy(a)

    width = 25mm
    layer = Layer(width=width)
    layer.γ = [1.0im]

    prop = exp(-layer.γ[1] * ustrip(Float64, u"m", width))
    s11 = a.s11 
    s12 = a.s12 * prop
    s21 = a.s21 * prop
    s22 = b.s22 * prop^2


    cascade!(a, layer, ustrip(Float64, mm, width))
    cascade!(b, layer)

    @test norm(a.s11 - b.s11) ≤ atol
    @test norm(a.s12 - b.s12) ≤ atol
    @test norm(a.s21 - b.s21) ≤ atol
    @test norm(a.s22 - b.s22) ≤ atol

    @test norm(a.s11 - s11) ≤ atol
    @test norm(a.s12 - s12) ≤ atol
    @test norm(a.s21 - s21) ≤ atol
    @test norm(a.s22 - s22) ≤ atol
end

@testset "GSM file" begin
    na1 = 1
    na2 = nb1 = 1
    nb2 = 1
    α = 0.35
    β = sqrt(1 - α^2)
    a = GSM(na1, na2)
    a.s11 .= a.s22 .= α
    a.s12 .= a.s21 .= β

    l1 = Layer(ϵᵣ = 2.0)
    l1.P = [TE, TM]
    l1.M = [0, 0]
    l1.N = [0, 0]
    
    l2 = deepcopy(l1)
    l2.ϵᵣ = 1
    
    sh = rectstrip(Lx=1, Ly=1, Px=1, Py=1, Nx=2, Ny=1, units=mm, fufp=true)

    fname = joinpath(tempdir(), "gsmfile.gsm")

    initialize_gsm_file(fname, l1, l2, sh, sh)
    case = Dict("θ"=>0.0, "ϕ"=>0.0, "FGHz"=>2.0)
    append_gsm_data(fname, "1", a, l1, l2, case)
    dat = read_gsm_file(fname)
    @test dat["ModelistIn"] == (TEorTM[TE, TM], [0, 0], [0, 0])
    @test dat["ModelistOut"] == (TEorTM[TE, TM], [0, 0], [0, 0])
    @test dat["s1in"] == [1.0, 0.0]
    @test dat["s1out"] == [1.0, 0.0]
    @test dat["s2in"] == [0.0, 1.0]
    @test dat["s2out"] == [0.0, 1.0]
    @test dat["1/gsm"].s11 == a.s11
    @test dat["1/gsm"].s12 == a.s12
    @test dat["1/gsm"].s21 == a.s21
    @test dat["1/gsm"].s22 == a.s22
    @test dat["1/epsmu1"] == (2.0 - 0.0im, 1.0 - 0.0im)
    @test dat["1/epsmu2"] == (1.0 - 0.0im, 1.0 - 0.0im)
    @test dat["1/case"] == case
    @test dat["unitsin"] == mm
    @test dat["unitsout"] == mm
end

@testset "choose_gblocks" begin
    strata = 
        [Layer()
         rectstrip(units=mm, Px = 1, Py = 1, Lx=1, Ly=1, Nx=2, Ny=2)
         Layer(width=20mm)
         Layer(width=1mil)
         Layer(width=10mm)
         Layer(width=1mil, ϵᵣ=2.2)
         rectstrip(units=mm, Px = 1, Py = 1, Lx=1, Ly=1, Nx=2, Ny=2)
         Layer()
        ]
    FGHz = 2.0; λ = c₀ / (FGHz*1e9); k0 = twopi / λ
    gbl = choose_gblocks(strata, k0)
    @test length(gbl) == 4
    @test gbl[1] == Gblock(1:1,1)
    @test gbl[2] == Gblock(2:2,0)
    @test gbl[3] == Gblock(3:3,0)
    @test gbl[4] == Gblock(4:5,5)
end