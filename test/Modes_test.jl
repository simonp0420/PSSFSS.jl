using PSSFSS
using PSSFSS.Modes: choose_layer_modes!
#using PSSFSS.Elements: rectstrip, polyring
using PSSFSS.Layers: Layer
using PSSFSS.Sheets: RWGSheet
#using PSSFSS.PSSFSSLen: mm, cm, inch, mil
using PSSFSS.GSMs: Gblock, choose_gblocks
using PSSFSS.Constants: c₀, twopi
using Test
using Logging: Error, ConsoleLogger, default_metafmt, global_logger

testlogger = ConsoleLogger(stderr, Error,
                       meta_formatter=default_metafmt, show_limited=true,
                       right_justify=0)
oldlogger = global_logger(testlogger)

@testset "choose_layer_modes!1" begin
    strata = 
    [   
        Layer(name="Vacuum")
        Layer(name="Duroid 6035", ϵᵣ=3.6, tanδ=0.0013, width=20mil)
        polyring(sides=4, orient=45, a=[0.1902], b=[0.2169], units=inch, 
                    s1=[0.5117, 0], s2=[0, 0.5117], ntri=100)
        Layer(name="Duroid 6035", ϵᵣ=3.6, tanδ=0.0013, width=20mil)
        Layer(name="Vacuum", width=0.2inch)
        Layer(name="Duroid 6035", ϵᵣ=3.6, tanδ=0.0013, width=20mil)
        polyring(sides=4, orient=45, a=[0.1902], b=[0.2169], units=inch, 
                    s1=[0.5117, 0], s2=[0, 0.5117], ntri=100)
        Layer(name="Duroid 6035", ϵᵣ=3.6, tanδ=0.0013, width=20mil)
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
    fmax = 8.0 * 1e9  # minimum frequency [Hz]
    k0max = twopi * fmax / c₀
    dbmin = 25.0

    gbl = choose_gblocks(layers, sheets, junc, k0min)
    choose_layer_modes!(layers, sheets, junc, gbl, k0max, dbmin)
    @test [length(l.P) for l in strata if l isa Layer] == [2,0,0,18,0,0,2]
end

@testset "choose_layer_modes!2" begin
    strata = 
    [   
        Layer(name="Vacuum")
        Layer(name="Duroid 6035", ϵᵣ=3.6, tanδ=0.0013, width=20mil)
        polyring(sides=4, orient=45, a=[0.1902], b=[0.2169], units=inch, 
                    s1=[0.5117, 0], s2=[0, 0.5117], ntri=100)
        Layer(name="Duroid 6035", ϵᵣ=3.6, tanδ=0.0013, width=20mil)
        Layer(name="Vacuum", width=0.1inch)
        Layer(name="Vacuum", width=0.1inch)
        Layer(name="Duroid 6035", ϵᵣ=3.6, tanδ=0.0013, width=20mil)
        polyring(sides=4, orient=45, a=[0.1902], b=[0.2169], units=inch, 
                    s1=[0.9005, 0], s2=[0, 0.9005], ntri=100)
        Layer(name="Duroid 6035", ϵᵣ=3.6, tanδ=0.0013, width=20mil)
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
    fmax = 8.0 * 1e9  # minimum frequency [Hz]
    k0max = twopi * fmax / c₀
    dbmin = 25.0

    gbl = choose_gblocks(layers, sheets, junc, k0min)
    choose_layer_modes!(layers, sheets, junc, gbl, k0max, dbmin)
    @test [length(l.P) for l in strata if l isa Layer] == [2,0,0,2,2,0,0,2]
end

@testset "choose_layer_modes!3" begin
    strata = 
    [   
        Layer(name="Vacuum")
        Layer(name="Duroid 6035", ϵᵣ=3.6, tanδ=0.0013, width=20mil)
        polyring(sides=4, orient=45, a=[0.1902], b=[0.2169], units=inch, 
                    s1=[0.5117, 0], s2=[0, 0.5117], ntri=100)
        Layer(name="Duroid 6035", ϵᵣ=3.6, tanδ=0.0013, width=20mil)
        Layer(name="Vacuum", width=0.1inch)
        Layer(name="Vacuum", width=0.1inch)
        Layer(name="Duroid 6035", ϵᵣ=3.6, tanδ=0.0013, width=20mil)
        polyring(sides=4, orient=45, a=[0.1902], b=[0.2169], units=inch, 
                    s1=[0.5117, 0], s2=[0, 0.5117], ntri=100)
        Layer(name="Duroid 6035", ϵᵣ=3.6, tanδ=0.0013, width=20mil)
        Layer(name="Vacuum")
    ]

    fmin = 2.0 * 1e9  # minimum frequency [Hz]
    k0min = twopi * fmin / c₀
    fmax = 8.0 * 1e9  # minimum frequency [Hz]
    k0max = twopi * fmax / c₀
    dbmin = 25.0
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

    gbl = choose_gblocks(layers, sheets, junc, k0min)
    choose_layer_modes!(layers, sheets, junc, gbl, k0max, dbmin)
    #@test_logs [(:info,),(:info,),(:info,)] choose_layer_modes!(strata, gbl, k0max, dbmin)
    @test [length(l.P) for l in strata if l isa Layer] == [2,0,0,42,0,0,0,2]

end

global_logger(oldlogger)
nothing
