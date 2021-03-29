using PSSFSS
using PSSFSS.RWG
using StaticArrays: SArray
using Test


@testset "fufptrue" begin
    sheet = rectstrip(Lx=1, Ly=1, Px=1, Py=1, Nx=2, Ny=1, units=mm, fufp=true)
    rwgdata = setup_rwg(sheet)
    @test rwgdata.bfe == [8  6  9  1  2  5
                          8  6  9  3  4  7]
    @test rwgdata.bff == [1  2  3  2  4  1
                          2  3  4  1  3  4]
    @test rwgdata.ebf == [4, 5, 4, 5, 6, 2, 6, 1, 3]
    @test rwgdata.eci == [3, 3, 4, 4, 1, 0, 2, 0, 0]
    @test rwgdata.ufpm == [1  5   9  11
                           2  6  10  12
                           3  7   1   5
                           4  8   2   6]
    @test rwgdata.ufp2fp == [[1, 11], [2, 12], [3], [4], [5, 15], [6, 16], [7],
                             [8], [9], [10], [13], [14]]
    @test rwgdata.nufp == 12
end

@testset "fufpfalse" begin
    sheet = rectstrip(Lx=1, Ly=1, Px=1, Py=1, Nx=2, Ny=1, units=mm, fufp=false)
    rwgdata = setup_rwg(sheet)
    @test rwgdata.bfe == [8  6  9  1  2  5
                          8  6  9  3  4  7]
    @test rwgdata.bff == [1  2  3  2  4  1
                          2  3  4  1  3  4]
    @test rwgdata.ebf == [4, 5, 4, 5, 6, 2, 6, 1, 3]
    @test rwgdata.eci == [3, 3, 4, 4, 1, 0, 2, 0, 0]
    @test rwgdata.ufpm == reshape(1:16, (4,4))
    @test rwgdata.ufp2fp == [[i] for i in 1:16]
    @test rwgdata.nufp == 16
end

@testset "rwgbfft!" begin
    sheet = rectstrip(Lx=1, Ly=1, Px=1, Py=1, Nx=2, Ny=1, units=mm, fufp=false)
    rwgdat = setup_rwg(sheet)
    nbf = size(rwgdat.bff,2)
    ft = zeros(SArray{Tuple{2},ComplexF64,1,2}, (nbf))
    k = [0.,0.]; ψ₁ = 0.0; ψ₂ = 0.0
    rwgbfft!(ft, rwgdat, sheet, k, ψ₁, ψ₂)
    @test ft[1] ≈ [0.0001666666666666667 + 0.0im, -0.0003333333333333334 + 0.0im]
    @test ft[2] ≈ [0.0003333333333333334 + 0.0im, 0.0003333333333333334 + 0.0im]
    @test ft[3] ≈ [0.00016666666666666663 + 0.0im, -0.0003333333333333334 + 0.0im]

    sheet = rectstrip(Lx=1, Ly=1, Px=1, Py=1, Nx=2, Ny=1, units=inch, fufp=false)
    rwgdat = setup_rwg(sheet)
    nbf = size(rwgdat.bff,2)
    ft = zeros(SArray{Tuple{2},ComplexF64,1,2}, (nbf))
    k = [20.0,-30.3]; ψ₁ = 0.8; ψ₂ = -0.7
    rwgbfft!(ft, rwgdat, sheet, k, ψ₁, ψ₂)
    @test ft[1] ≈ [0.003986416685586114 - 0.0010511300702857503im, -0.008079847138694855 + 0.0021304773084817956im]
    @test ft[2] ≈ [0.008263557134402535 - 0.0010871639073400372im, 0.008284527320353664 - 0.0010899227712197804im]
    @test ft[3] ≈ [0.004122638196063312 - 1.5707327530113332e-5im, -0.008355946971820058 + 3.1836312009155855e-5im]
end