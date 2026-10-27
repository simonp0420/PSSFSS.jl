# amr_test.jl -- Unit tests for adaptive mesh refinement

using PSSFSS: rectstrip, facecount, @u_str, mm
using PSSFSS.RWG: setup_rwg
using PSSFSS.AutoMeshRefine: tricharge!, tricharge2!

#using StaticArrays: SArray
using Random: seed!
using Test

seed!(1)

@testset "tricharge_test" begin
    omega = 23.2
    sheet = rectstrip(Lx=1, Ly=1, Px=1, Py=1, Nx=3, Ny=3, units=mm)
    sheet.ψ₁ = π / 3
    sheet.ψ₂ = π / 4
    rwg = setup_rwg(sheet)
    charges1, charges2 = (zeros(ComplexF64, facecount(sheet)) for _ in 1:2)

    nbf = size(rwg.bfe, 2)
    currents = rand(ComplexF64, nbf)
    tricharge!(charges1, currents, rwg, sheet, omega)
    tricharge2!(charges2, currents, rwg, sheet, omega)
    @test charges1 ≈ charges2
end
using PSSFSS: RWGSheet, cm
using Unitful: ustrip
using PSSFSS.AutoMeshRefine: tricurrent!
using StaticArrays: SVector

# Two triangles sharing edge (1,2). Local edge k is opposite vertex k.
# Explicit topology keeps the analytic expectations independent of meshing code.
function two_triangle_sheet(points; units=cm, periodic=false)
    sheet = RWGSheet()
    sheet.units = units
    scale = ustrip(Float64, units, 1u"m")
    sheet.ρ = [SVector{2,Float64}(scale .* p) for p in points]
    sheet.e1 = [1, 2, 3, 1, 2]
    sheet.e2 = [2, 3, 1, 4, 4]
    sheet.fv = [3 4; 1 1; 2 2]
    sheet.fe = [1 1; 2 5; 3 4]
    sheet.s₁ = SVector(scale, 0.0)
    sheet.s₂ = SVector(0.0, scale)
    sheet.β₁ = SVector(2pi / scale, 0.0)
    sheet.β₂ = SVector(0.0, 2pi / scale)
    sheet.ξη_check = periodic
    return sheet
end

@testset "Analytic single-basis charge and current" begin
    # Shared edge has length 2 m; areas are 1 and 2 m². PSSFSS's
    # modified RWG basis has unit integrated flux, not unit normal density:
    # J+ = I(r-r+)/2A+, J- = -I(r-r-)/2A-, rho = i div(J)/omega.
    points = [(0.0, 0.0), (2.0, 0.0), (0.0, 1.0), (0.0, -2.0)]
    amplitude = 2.0 - 3.0im
    omega = 7.0
    for units in (cm, mm)
        sheet = two_triangle_sheet(points; units)
        rwg = setup_rwg(sheet)
        @test size(rwg.bfe, 2) == 1
        currents = [amplitude]
        expected_charge = (im * amplitude / omega) .* [1.0, -0.5]
        expected_current = amplitude .* [SVector(1/3, -1/3), SVector(-1/6, -1/3)]
        for charge! in (tricharge!, tricharge2!)
            charges = fill(99.0 + 12.0im, 2)
            charge!(charges, currents, rwg, sheet, omega)
            @test charges ≈ expected_charge
            @test charges[1] + 2charges[2] ≈ 0 atol=1e-14 # integrated charge
            charge!(charges, currents, rwg, sheet, 2omega)
            @test charges ≈ expected_charge / 2
            charge!(charges, zero.(currents), rwg, sheet, omega)
            @test all(iszero, charges) # output must be overwritten
        end
        current_density = fill(SVector(99.0 + 12im, -23.0 + 4im), 2)
        tricurrent!(current_density, currents, rwg, sheet)
        @test current_density ≈ expected_current
        tricurrent!(current_density, zero.(currents), rwg, sheet)
        @test all(iszero, current_density)
    end

    # Geometric dilation: J scales as 1/length and rho as 1/length².
    sheet = two_triangle_sheet([3 .* p for p in points])
    rwg = setup_rwg(sheet)
    charges = zeros(ComplexF64, 2)
    tricharge!(charges, [amplitude], rwg, sheet, omega)
    @test charges ≈ (im * amplitude / omega) .* [1.0, -0.5] / 9
    current_density = zeros(SVector{2,ComplexF64}, 2)
    tricurrent!(current_density, [amplitude], rwg, sheet)
    @test current_density ≈ amplitude .* [SVector(1/3, -1/3), SVector(-1/6, -1/3)] / 3
end

@testset "Analytic periodic half-basis phase factors" begin
    sheet = two_triangle_sheet([(0.0, 0.0), (1.0, 1.0), (0.0, 1.0), (1.0, 0.0)]; periodic=true)
    sheet.ψ₁ = pi / 3
    sheet.ψ₂ = -pi / 4
    rwg = setup_rwg(sheet)
    @test size(rwg.bfe, 2) == 3
    amplitude = 2.0 - 3.0im
    omega = 7.0
    # Areas are 1/2. xi+ is on face 1; eta+ is on face 2.
    # The minus half-basis acquires exp(-i psi) at the opposite boundary.
    for (code, phase, charge_pattern, direction) in
        ((1, cis(-sheet.ψ₁), [1.0, -cis(-sheet.ψ₁)], SVector(-2/3, -1/3)),
         (3, cis(-sheet.ψ₂), [-cis(-sheet.ψ₂), 1.0], SVector(-1/3, -2/3)))
        bf = only(findall(b -> rwg.eci[rwg.bfe[1, b]] == code, axes(rwg.bfe, 2)))
        currents = zeros(ComplexF64, 3)
        currents[bf] = amplitude
        for charge! in (tricharge!, tricharge2!)
            charges = zeros(ComplexF64, 2)
            charge!(charges, currents, rwg, sheet, omega)
            @test charges ≈ (2im * amplitude / omega) .* charge_pattern
        end
        expected = code == 1 ? [amplitude * direction, phase * amplitude * direction] :
                               [phase * amplitude * direction, amplitude * direction]
        current_density = zeros(SVector{2,ComplexF64}, 2)
        tricurrent!(current_density, currents, rwg, sheet)
        @test current_density ≈ expected
    end
end

@testset "Uniform periodic current has zero charge" begin
    sheet = two_triangle_sheet([(0.0, 0.0), (1.0, 1.0), (0.0, 1.0), (1.0, 0.0)]; periodic=true)
    rwg = setup_rwg(sheet)
    target = SVector(2.0 - 3im, -1.0 + 4im)
    # Coefficients are integrated outward fluxes on the plus-face edges:
    # diagonal normal points (1,-1), xi=0 normal (-1,0), eta=0 (0,-1).
    currents = zeros(ComplexF64, size(rwg.bfe, 2))
    for bf in axes(rwg.bfe, 2)
        code = rwg.eci[rwg.bfe[1, bf]]
        currents[bf] = code == 0 ? target[1] - target[2] : code == 1 ? -target[1] : -target[2]
    end
    for charge! in (tricharge!, tricharge2!)
        charges = zeros(ComplexF64, 2)
        charge!(charges, currents, rwg, sheet, 7.0)
        @test charges ≈ zeros(ComplexF64, 2) atol=1e-14
    end
    current_density = zeros(SVector{2,ComplexF64}, 2)
    tricurrent!(current_density, currents, rwg, sheet)
    @test current_density ≈ [target, target]
end
