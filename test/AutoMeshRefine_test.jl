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