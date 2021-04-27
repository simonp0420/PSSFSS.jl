using PSSFSS
using LinearAlgebra: norm, diagind
using PSSFSS.GSMs: GSM, cascade, cascade!, initialize_gsm_file, 
                  append_gsm_data, read_gsm_file, choose_gblocks,
                  Gblock
using PSSFSS.Layers: Layer, TEorTM, TE, TM
using PSSFSS.PSSFSSLen
using PSSFSS.Elements: rectstrip
using PSSFSS.Sheets: RWGSheet
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

@testset "choose_gblocks1" begin
    strata = 
    [Layer()
     rectstrip(units=cm, Px = 1, Py = 1, Lx=0.5, Ly=0.5, Nx=10, Ny=10)
     Layer(width=1mil, ϵᵣ=2.2)
     Layer(width=1mil)
     Layer(width=3mm)
     rectstrip(units=cm, Px = 1, Py = 1, Lx=0.5, Ly=0.5, Nx=10, Ny=10)
     Layer(width=1mil, ϵᵣ=2.2)
     Layer()
    ]
    islayer = map(x -> x isa Layer, strata)
    issheet = map(x -> x isa RWGSheet, strata)
    layers = convert(Vector{Layer}, strata[islayer])
    sheets = convert(Vector{RWGSheet}, strata[issheet])
    nl = length(layers)
    nj = nl - 1
    ns = length(sheets)
    sint = cumsum(islayer)[issheet] # sint[k] contains dielectric interface number of k'th sheet 
    junc = zeros(Int, nj)
    junc[sint] = 1:ns #  junc[i] is the sheet number present at interface i, or 0 if no sheet is there
    
    FGHz = 10.0; k0 = twopi * FGHz*1e9 / c₀
    gbl = choose_gblocks(layers, sheets, junc, k0)
    @test length(gbl) == 2
    @test gbl[1] == Gblock(1:3,1)
    @test gbl[2] == Gblock(4:5,4)
end

@testset "choose_gblocks2" begin
    duroid = Layer(name="Duroid 6035", ϵᵣ=3.6, tanδ=0.0013, width=20mil)
    strata = 
    [   
        Layer(name="Vacuum")
        duroid
        polyring(sides=4, orient=45, a=[0.1902], b=[0.2169], units=inch, 
                    s1=[0.5117, 0], s2=[0, 0.5117], ntri=100)
        duroid
        Layer(name="Vacuum", width=0.05inch)
        duroid
        polyring(sides=4, orient=45, a=[0.1902], b=[0.2169], units=inch, 
                    s1=[0.5117, 0], s2=[0, 0.5117], ntri=100)
        duroid
        Layer(name="Vacuum")
    ]

    islayer = map(x -> x isa Layer, strata)
    issheet = map(x -> x isa RWGSheet, strata)
    layers = convert(Vector{Layer}, strata[islayer])
    sheets = convert(Vector{RWGSheet}, strata[issheet])
    nl = length(layers)
    nj = nl - 1
    ns = length(sheets)
    sint = cumsum(islayer)[issheet] # sint[k] contains dielectric interface number of k'th sheet 
    junc = zeros(Int, nj)
    junc[sint] = 1:ns #  junc[i] is the sheet number present at interface i, or 0 if no sheet is there
    fmin = 2.0 * 1e9  # minimum frequency [Hz]
    k0min = twopi * fmin / c₀
    gbl = choose_gblocks(layers, sheets, junc, k0min)
    @test gbl[1] == Gblock(1:3,2)
    @test gbl[2] == Gblock(4:6,5)
end

@testset "choose_gblocks3" begin
    duroid = Layer(name="Duroid 6035", ϵᵣ=3.6, tanδ=0.0013, width=20mil)
    strata = 
    [   
        Layer(name="Vacuum")
        duroid
        polyring(sides=4, orient=45, a=[0.1902], b=[0.2169], units=inch, 
                    s1=[0.5117, 0], s2=[0, 0.5117], ntri=100)
        duroid
        duroid
        polyring(sides=4, orient=45, a=[0.1902], b=[0.2169], units=inch, 
                    s1=[0.5117, 0], s2=[0, 0.5117], ntri=100)
        duroid
        Layer(name="Vacuum")
    ]

    islayer = map(x -> x isa Layer, strata)
    issheet = map(x -> x isa RWGSheet, strata)
    layers = convert(Vector{Layer}, strata[islayer])
    sheets = convert(Vector{RWGSheet}, strata[issheet])
    nl = length(layers)
    nj = nl - 1
    ns = length(sheets)
    sint = cumsum(islayer)[issheet] # sint[k] contains dielectric interface number of k'th sheet 
    junc = zeros(Int, nj)
    junc[sint] = 1:ns #  junc[i] is the sheet number present at interface i, or 0 if no sheet is there
    fmin = 2.0 * 1e9  # minimum frequency [Hz]
    k0min = twopi * fmin / c₀
    @test choose_gblocks(layers, sheets, junc, k0min) == [Gblock(1:2,2), Gblock(3:5,4)]
end

@testset "choose_gblocks4" begin
    duroid = Layer(name="Duroid 6035", ϵᵣ=3.6, tanδ=0.0013, width=20mil)
    strata = 
    [   
        Layer(name="Vacuum")
        duroid
        polyring(sides=4, orient=45, a=[0.1902], b=[0.2169], units=inch, 
                    s1=[0.5117, 0], s2=[0, 0.5117], ntri=100)
        duroid
        Layer(name="Vacuum", width=5mil)
        duroid
        polyring(sides=4, orient=45, a=[0.1902], b=[0.2169], units=inch, 
                    s1=[0.5117, 0], s2=[0, 0.5117], ntri=100)
        duroid
        Layer(name="Vacuum", width=5mil)
        duroid
        polyring(sides=4, orient=45, a=[0.1902], b=[0.2169], units=inch, 
                    s1=[0.5117, 0], s2=[0, 0.5117], ntri=100)
        duroid
        Layer(name="Vacuum")
    ]

    islayer = map(x -> x isa Layer, strata)
    issheet = map(x -> x isa RWGSheet, strata)
    layers = convert(Vector{Layer}, strata[islayer])
    sheets = convert(Vector{RWGSheet}, strata[issheet])
    nl = length(layers)
    nj = nl - 1
    ns = length(sheets)
    sint = cumsum(islayer)[issheet] # sint[k] contains dielectric interface number of k'th sheet 
    junc = zeros(Int, nj)
    junc[sint] = 1:ns #  junc[i] is the sheet number present at interface i, or 0 if no sheet is there
    fmin = 2.0 * 1e9  # minimum frequency [Hz]
    k0min = twopi * fmin / c₀

    @test choose_gblocks(layers, sheets, junc, k0min) == [Gblock(1:2,2), Gblock(3:5,5),Gblock(6:9,8)]

end

