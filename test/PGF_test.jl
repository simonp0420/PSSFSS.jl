using PSSFSS
using LinearAlgebra: ×, norm
using PSSFSS.Elements: s₁s₂2β₁β₂
using PSSFSS.PGF: direct_electric_modal_series, direct_magnetic_modal_series, jksums, c3_calc, d3_calc
using Logging: Error, ConsoleLogger, default_metafmt, global_logger
using Test



testlogger = ConsoleLogger(stderr, Error,
                       meta_formatter=default_metafmt, show_limited=true,
                       right_justify=0)
oldlogger = global_logger(testlogger)


compute_exact_from_scratch = false

FGHz = 2.0; wl_inch = 11.80285 / FGHz
inch2m = 2.54/100
wl_m = wl_inch * inch2m
k0 = 2π / wl_m
layers = Layer[Layer()
         Layer(width=20mil, ϵᵣ=3.6, tanδ=0.0013, μᵣ=5.0)
         Layer(width=10mil, ϵᵣ=8.0, tanδ=0.0013, μᵣ=2.0)
         Layer(width=50mil, ϵᵣ=2.6, tanδ=0.0003)
         Layer()]
s = 2 # Source location
ψ₁, ψ₂ = 0.4, -0.6  # Incremental phase shifts (rad)
s₁ = [1.14*inch2m, 0.0] # Lattice vector
s₂ = [0.5707, 0.9885] * inch2m # Lattice vector
β₁, β₂ = s₁s₂2β₁β₂(s₁,s₂)
β₀₀ = (ψ₁ * β₁ + ψ₂ * β₂) / (2π)
extract = true
ρdif_mil = [2, 6]; ρdif = ρdif_mil / 1000 * 2.54/100

"Precomputed high-accuracy modal series for electric sources"
function Σm_exact(k0, u, ψ₁, ψ₂, layers, s, β₁, β₂, β₀₀, ρdif)
    if u == 0.1 * max(norm(β₁), norm(β₂))
        Σm1_exact, Σm2_exact = (-232.15311689927114 - 10.286495939161526im,
                                -11.778556545589899 - 160.87504453495308im)
    elseif u == 0.9 * max(norm(β₁), norm(β₂))
        Σm1_exact, Σm2_exact = (-86.47040978143367 - 9.887172408455307im,
                                133.51945406766538 - 160.47675058617088im)
    else
        # Warning: This takes about 15 seconds!
        (Σm1,Σm2) = direct_electric_modal_series(k0, u, ψ₁, ψ₂, layers, s, β₁, β₂, β₀₀, ρdif)
        Σm1_exact, Σm2_exact = sum(Σm1), sum(Σm2)
    end
    return (Σm1_exact, Σm2_exact)
end

"Precomputed high-accuracy modal series for magnetic sources"
function Σpm_exact(k0, u, ψ₁, ψ₂, layers, s, β₁, β₂, β₀₀, ρdif)
    if u == 0.1 * max(norm(β₁), norm(β₂))
        Σpm1_exact, Σpm2_exact = (-5905.287097073966 - 116.24899221632135im,
                                  -65.290894751642 - 115.6027231717079im)
    elseif u == 0.9 * max(norm(β₁), norm(β₂))
        Σpm1_exact, Σpm2_exact = (-2516.511258383936 - 111.36619289656556im,
                                  138.66489521333057 - 115.0436702287188im)
    else
        # Warning: This takes about 15 seconds!
        (Σpm1,Σpm2) = direct_electric_modal_series(k0, u, ψ₁, ψ₂, layers, s, β₁, β₂, β₀₀, ρdif)
        Σpm1_exact, Σpm2_exact = sum(Σpm1), sum(Σpm2)
    end
    return (Σpm1_exact, Σpm2_exact)
end

conv = 2e-4

@testset "Modal Series Test" begin
    for ufact in [0.1, 0.9]
        u = ufact * max(norm(β₁), norm(β₂))
        (Σm1_exact, Σm2_exact) = Σm_exact(k0, u, ψ₁, ψ₂, layers, s, β₁, β₂, β₀₀, ρdif)
        (Σm1_func, Σm2_func) = electric_modal_sum_funcs(k0, u, ψ₁, ψ₂, layers, s, β₁, β₂, β₀₀)
        (Σm1_fast, Σm2_fast) = (Σm1_func(ρdif), Σm2_func(ρdif))
        @test abs(Σm1_fast - Σm1_exact)/abs(Σm1_exact) < conv
        @test abs(Σm2_fast - Σm2_exact)/abs(Σm2_exact) < conv
        (Σpm1_exact, Σpm2_exact) = Σpm_exact(k0, u, ψ₁, ψ₂, layers, s, β₁, β₂, β₀₀, ρdif)
        (Σpm1_func, Σpm2_func) = magnetic_modal_sum_funcs(k0, u, ψ₁, ψ₂, layers, s, β₁, β₂, β₀₀)
        (Σpm1_fast, Σpm2_fast) = (Σpm1_func(ρdif), Σpm2_func(ρdif))
        @test abs(Σpm1_fast - Σpm1_exact)/abs(Σpm1_exact) < conv
        @test abs(Σpm2_fast - Σpm2_exact)/abs(Σpm2_exact) < conv
    end
end


function full_pgfs(k0, u, ψ₁, ψ₂, layers, s, β₁, β₂, β₀₀, ρdif, Σm1, Σm2, Σpm1, Σpm2, Σs1, Σs2)
    μ̃ = 2 * layers[s].μᵣ * layers[s+1].μᵣ / (layers[s].μᵣ + layers[s+1].μᵣ) # Normalized to μ0.
    ϵ̄ = (layers[s].ϵᵣ + layers[s+1].ϵᵣ) / 2 # Normalized to ϵ0.
    # Obtain the appropriate Green's function expansion coefficients:
    c3 = c3_calc(k0, u, layers[s].μᵣ, layers[s].ϵᵣ, layers[s+1].μᵣ, layers[s+1].ϵᵣ)
    d3 = d3_calc(k0, u, layers[s].μᵣ, layers[s].ϵᵣ, layers[s+1].μᵣ, layers[s+1].ϵᵣ)
    GAxx = μ̃ * (Σm1 + u/(4π) * (Σs1 + c3/(u^2) * Σs2)) # Divided by μ₀
    GΦ = 1/ϵ̄ * (Σm2 + u/(4π) * (Σs1 + d3/(u^2) * Σs2)) # Multiplied by ϵ₀

    c3s = c3_calc(k0, u, layers[s].ϵᵣ, layers[s].μᵣ, layers[s].ϵᵣ, layers[s].μᵣ)
    c3sp1 = c3_calc(k0, u, layers[s+1].ϵᵣ, layers[s+1].μᵣ, layers[s+1].ϵᵣ, layers[s+1].μᵣ)
    d3s = d3_calc(k0, u, layers[s].ϵᵣ, layers[s].μᵣ, layers[s].ϵᵣ, layers[s].μᵣ)
    d3sp1 = d3_calc(k0, u, layers[s+1].ϵᵣ, layers[s+1].μᵣ, layers[s+1].ϵᵣ, layers[s+1].μᵣ)
    GFxx = -(Σpm1 + u*ϵ̄/π * Σs1 + (c3s * layers[s].ϵᵣ + c3sp1 * layers[s+1].ϵᵣ)/(2π*u) * Σs2)
    GΨ = Σpm2 + u/(π*μ̃) * Σs1 + (d3s/layers[s].μᵣ + d3sp1/layers[s+1].μᵣ)/(2π*u) * Σs2 # Multiplied by ϵ₀
    return (GAxx, GΦ, GFxx, GΨ)
end


@testset "Full PGFs" begin
    GAxx = zeros(ComplexF64,2)
    GΦ = zeros(ComplexF64,2)
    GFxx = zeros(ComplexF64,2)
    GΨ = zeros(ComplexF64,2)
    for (i,ufact) in enumerate([0.1, 0.9])
        u = ufact*max(norm(β₁), norm(β₂))
        (Σs1, Σs2) = jksums(u*ρdif, ψ₁, ψ₂, u*s₁, u*s₂, extract)
        (Σm1_func, Σm2_func) = electric_modal_sum_funcs(k0, u, ψ₁, ψ₂, layers, s, β₁, β₂, β₀₀)
        (Σm1, Σm2) = (Σm1_func(ρdif), Σm2_func(ρdif))
        (Σpm1_func, Σpm2_func) = magnetic_modal_sum_funcs(k0, u, ψ₁, ψ₂, layers, s, β₁, β₂, β₀₀)
        (Σpm1, Σpm2) = (Σpm1_func(ρdif), Σpm2_func(ρdif))
        (GAxx[i], GΦ[i], GFxx[i], GΨ[i]) = full_pgfs(k0, u, ψ₁, ψ₂, layers, s, β₁, β₂, β₀₀, ρdif,
                                                     Σm1, Σm2, Σpm1, Σpm2, Σs1, Σs2)
    end
    @test 2*abs(GAxx[2]-GAxx[1])/abs(GAxx[2]+GAxx[1]) < conv
    @test 2*abs(GΦ[2]-GΦ[1])/abs(GΦ[2]+GΦ[1]) < conv
    @test 2*abs(GFxx[2]-GFxx[1])/abs(GFxx[2]+GFxx[1]) < conv
    @test 2*abs(GΨ[2]-GΨ[1])/abs(GΨ[2]+GΨ[1]) < conv
end

oldlogger = global_logger(testlogger)
nothing