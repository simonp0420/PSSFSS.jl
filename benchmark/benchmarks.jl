using BenchmarkTools
using PSSFSS
using LinearAlgebra: ×, norm
using PSSFSS.Elements: s₁s₂2β₁β₂
using PSSFSS.Sheets: SV2, Sheet
using PSSFSS.PGF: direct_electric_modal_series, direct_magnetic_modal_series, jksums, c3_calc,
    d3_calc, electric_modal_sum_funcs, magnetic_modal_sum_funcs
using PSSFSS.FillZY: fillz, filly
using PSSFSS.RWG: setup_rwg
using PSSFSS.GSMs: GSM, cascade
    

const SUITE = @benchmarkset "PSSFSS" begin
    @benchmarkset "PGF" begin 
        FGHz = 2.0
        wl_inch = 11.80285 / FGHz
        inch2m = 2.54 / 100
        wl_m = wl_inch * inch2m
        k0 = 2π / wl_m
        layers = Layer[Layer()
            Layer(width=20mil, ϵᵣ=3.6, tanδ=0.0013, μᵣ=5.0)
            Layer(width=10mil, ϵᵣ=8.0, tanδ=0.0013, μᵣ=2.0)
            Layer(width=50mil, ϵᵣ=2.6, tanδ=0.0003)
            Layer()]
        s = 2 # Source location
        ψ₁, ψ₂ = 0.4, -0.6  # Incremental phase shifts (rad)
        s₁ = SV2([1.14 * inch2m, 0.0]) # Lattice vector
        s₂ = SV2([0.5707, 0.9885] * inch2m) # Lattice vector
        β₁, β₂ = s₁s₂2β₁β₂(s₁, s₂)
        ufact = 0.5
        u = ufact * max(norm(β₁), norm(β₂))
        β₀₀ = (ψ₁ * β₁ + ψ₂ * β₂) / (2π)
        extract = true
        ρdifs_mil = [SV2([2, 6]), SV2([200, 0]), SV2([0, 800]), SV2([570, 980])] 
        ρdifs = ρdifs_mil .* (inv(1000) * 0.0254)
        
        @benchmarkset "jksums" begin
            @case "jksums1" jksums($u * $ρdifs[1], $ψ₁, $ψ₂, $u * $s₁, $u * $s₂, $extract)
            @case "jksums2" jksums($u * $ρdifs[2], $ψ₁, $ψ₂, $u * $s₁, $u * $s₂, $extract)
            @case "jksums3" jksums($u * $ρdifs[3], $ψ₁, $ψ₂, $u * $s₁, $u * $s₂, $extract)
            @case "jksums4" jksums($u * $ρdifs[4], $ψ₁, $ψ₂, $u * $s₁, $u * $s₂, $extract)
        end
        @benchmarkset "modal sum funcs" begin
            @case "electric" electric_modal_sum_funcs($k0, $u, $ψ₁, $ψ₂, $layers, $s, $β₁, $β₂, $β₀₀)
            @case "magnetic" magnetic_modal_sum_funcs($k0, $u, $ψ₁, $ψ₂, $layers, $s, $β₁, $β₂, $β₀₀)
        end
    end
    @benchmarkset "ZY Fill" begin
        FGHz = 2.0
        wl_inch = 11.80285 / FGHz
        inch2m = 2.54 / 100
        wl_m = wl_inch * inch2m
        k0 = 2π / wl_m
        @benchmarkset "Z" begin
            strata = [Layer()
            Layer(ϵᵣ=3.4, tanδ=0.003, width=1mil)
            rectstrip(Nx=20, Ny=20, Lx=20, Ly=50, Px=100, Py=100, units=mil)
            Layer(ϵᵣ=3.4, tanδ=0.003, width=1mil)]
    
            layers = [s for s in strata if s isa Layer]
    
            s = findfirst(x -> x isa Sheet, strata) - 1 # location of Sheet
            metal = [x for x in strata if x isa Sheet][1]
            β₁, β₂ = metal.β₁, metal.β₂
            ufact = 0.5
            units_per_meter = ustrip(Float64, metal.units, 1u"m")
            u = ufact * max(norm(β₁), norm(β₂)) * units_per_meter
            rwgdat = setup_rwg(metal)
            ψ₁ = ψ₂ = 0.0
            @case "normalinc fufp" fillz($k0, $u, $layers, $s, $ψ₁, $ψ₂, $metal, $rwgdat) setup=GC.gc()
            ψ₁ = ψ₂ = 0.05
            @case "obliqueinc fufp" fillz($k0, $u, $layers, $s, $ψ₁, $ψ₂, $metal, $rwgdat) setup=GC.gc()
            metal = rectstrip(Nx=20, Ny=20, Lx=20, Ly=50, Px=100, Py=100, units=mil,fufp=false)
            rwgdat = setup_rwg(metal)
            ψ₁ = ψ₂ = 0.0
            @case "normalinc !fufp" fillz($k0, $u, $layers, $s, $ψ₁, $ψ₂, $metal, $rwgdat) setup=GC.gc()
            ψ₁ = ψ₂ = 0.05
            @case "obliqueinc !fufp" fillz($k0, $u, $layers, $s, $ψ₁, $ψ₂, $metal, $rwgdat) setup=GC.gc()
        end

        @benchmarkset "Y" begin
            strata = [Layer()
            Layer(ϵᵣ=3.4, tanδ=0.003, width=1mil)
            rectstrip(Nx=20, Ny=20, Lx=20, Ly=50, Px=100, Py=100, units=mil, class='M')
            Layer(ϵᵣ=3.4, tanδ=0.003, width=1mil)]
    
            layers = [s for s in strata if s isa Layer]
    
            s = findfirst(x -> x isa Sheet, strata) - 1 # location of Sheet
            metal = [x for x in strata if x isa Sheet][1]
            β₁, β₂ = metal.β₁, metal.β₂
            ufact = 0.5
            units_per_meter = ustrip(Float64, metal.units, 1u"m")
            u = ufact * max(norm(β₁), norm(β₂)) * units_per_meter
            rwgdat = setup_rwg(metal)
            ψ₁ = ψ₂ = 0.0
            @case "normalinc fufp" filly($k0, $u, $layers, $s, $ψ₁, $ψ₂, $metal, $rwgdat) setup=GC.gc()
            ψ₁ = ψ₂ = 0.05
            @case "obliqueinc fufp" filly($k0, $u, $layers, $s, $ψ₁, $ψ₂, $metal, $rwgdat) setup=GC.gc()
            metal = rectstrip(Nx=20, Ny=20, Lx=20, Ly=50, Px=100, Py=100, units=mil, class='M', fufp=false)
            ψ₁ = ψ₂ = 0.0
            @case "normalinc !fufp" filly($k0, $u, $layers, $s, $ψ₁, $ψ₂, $metal, $rwgdat) setup=GC.gc()
            ψ₁ = ψ₂ = 0.05
            @case "obliqueinc !fufp" filly($k0, $u, $layers, $s, $ψ₁, $ψ₂, $metal, $rwgdat) setup=GC.gc()
        end

    

    end

    @benchmarkset "Cascade" begin
        n1, n2, n3, n4 = 20, 100, 1000, 2000
        a11 = 0.3 * rand(ComplexF64, n1, n1)
        a12 = rand(ComplexF64, n1, n2)
        a22 = 0.3 * rand(ComplexF64, n2, n2)
        a = GSM(a11, a12, transpose(a12), a22)
        b11 = 0.3 * rand(ComplexF64, n2, n2)
        b12 = rand(ComplexF64, n2, n3)
        b22 = 0.3 * rand(ComplexF64, n3, n3)
        b = GSM(b11, b12, transpose(b12), b22)
        ab = cascade(a,b)
        c11 = 0.3 * rand(ComplexF64, n3, n3)
        c12 = rand(ComplexF64, n3, n4)
        c22 = 0.3 * rand(ComplexF64, n4, n4)
        c = GSM(c11, c12, transpose(c12), c22)
        bc = cascade(b,c)
        d11 = 0.3 * rand(ComplexF64, n4, n4)
        d12 = rand(ComplexF64, n4, n1)
        d22 = 0.3 * rand(ComplexF64, n1, n1)
        d = GSM(d11, d12, transpose(d12), d22)
        @case "a ⋆ b" cascade($a, $b);
        @case "b ⋆ c" cascade($b, $c);
        @case "(a ⋆ b) ⋆ c" cascade($ab, $c);
        @case "a ⋆ (b ⋆ c)" cascade($a, $bc);
        @case "d ⋆ a" cascade($d, $a);
    end
end

