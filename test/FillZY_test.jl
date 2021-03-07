using LinearAlgebra: norm
using PSSFSS
using PSSFSS.FillZY: fillz, filly
using PSSFSS.Elements: s₁s₂2β₁β₂
using PSSFSS.Sheets: Sheet
using PSSFSS.RWG: setup_rwg
using Unitful: ustrip
using Logging: Error, ConsoleLogger, default_metafmt, global_logger
using Test



testlogger = ConsoleLogger(stderr, Error,
                       meta_formatter=default_metafmt, show_limited=true,
                       right_justify=0)
oldlogger = global_logger(testlogger)


@testset "fillz" begin

    FGHz = 2.0; wl_inch = 11.80285 / FGHz
    inch2m = 2.54/100
    wl_m = wl_inch * inch2m
    k0 = 2π / wl_m
    strata = [ Layer()
               Layer(ϵᵣ=3.4, tanδ=0.003, width=1mil)
               rectstrip(Nx=20, Ny=20, Lx=20, Ly=50, Px=100, Py=100, units=mil)
               Layer(ϵᵣ=3.4, tanδ=0.003, width=1mil)
             ]
         
    layers = [s for s in strata if s isa Layer]

    s = findfirst(x -> x isa Sheet, strata) - 1 # location of Sheet
    metal = [x for x in strata if x isa Sheet][1]
    β₁, β₂ = metal.β₁, metal.β₂
    ufact = 0.5
    units_per_meter = ustrip(Float64, metal.units, 1u"m")
    u = ufact * max(norm(β₁), norm(β₂)) * units_per_meter

    rwgdat = setup_rwg(metal)
    ψ₁ = ψ₂ = 0.0

    zmat = fillz(k0,u,layers,s,ψ₁,ψ₂,metal,rwgdat)

    @test abs(zmat[20,30] - (0.8990238712012797 - 354.1682770720877im)) < 0.005
    @test abs(zmat[100,100] - (110.95156005861463 - 37529.820828933574im)) < 0.01

end


@testset "filly" begin

    FGHz = 2.0; wl_inch = 11.80285 / FGHz
    inch2m = 2.54/100
    wl_m = wl_inch * inch2m
    k0 = 2π / wl_m
    strata = [ Layer()
              Layer(ϵᵣ=3.4, tanδ=0.003, width=1mil)
              rectstrip(class='M', Nx=20, Ny=20, Lx=20, Ly=50, Px=100, Py=100, units=mil)
              Layer(ϵᵣ=3.4, tanδ=0.003, width=1mil)
            ]
        
    layers = [s for s in strata if s isa Layer]

    s = findfirst(x -> x isa Sheet, strata) - 1 # location of Sheet
    apert = [x for x in strata if x isa Sheet][1]
    β₁, β₂ = apert.β₁, apert.β₂
    ufact = 0.5
    units_per_meter = ustrip(Float64, apert.units, 1u"m")
    u = ufact * max(norm(β₁), norm(β₂)) * units_per_meter

    rwgdat = setup_rwg(apert)
    ψ₁ = ψ₂ = 0.0

    ymat = filly(k0,u,layers,s,ψ₁,ψ₂,apert,rwgdat)

    @test abs(ymat[20,30] - (1.2168396602563192e-6 - 0.02476973884018213im)) < 0.0001
    @test abs(ymat[100,100] - (2.2127855040521283e-6 - 3.5131557304457406im)) < 0.0001 

end



global_logger(oldlogger)


nothing