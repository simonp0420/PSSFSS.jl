using BenchmarkTools
using PSSFSS
using LinearAlgebra: ×, norm
using PSSFSS.Elements: s₁s₂2β₁β₂
using PSSFSS.Sheets: SV2
using PSSFSS.PGF: direct_electric_modal_series, direct_magnetic_modal_series, jksums, c3_calc,
    d3_calc, electric_modal_sum_funcs, magnetic_modal_sum_funcs

#=
using Logging: Error, ConsoleLogger, default_metafmt, global_logger

testlogger = ConsoleLogger(stderr, Error,
    meta_formatter=default_metafmt, show_limited=true,
    right_justify=0)
oldlogger = global_logger(testlogger)
=#


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



const SUITE = @benchmarkset "PSSFSS" begin  # continue here
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

