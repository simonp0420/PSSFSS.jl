using PSSFSS
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

