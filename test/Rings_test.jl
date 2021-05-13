using PSSFSS
using PSSFSS.Rings: Ring
using Test

@testset "Rings" begin
    @test collect(Ring(0)) == [(0,0)]
    @test collect(Ring(1)) == [(1, -1), (1, 0), (1, 1), (-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1)]
    @test collect(Ring(2)) == [(2, -2), (2, -1), (2, 0), (2, 1), (2, 2), (-2, -2), (-2, -1), (-2, 0),
                               (-2, 1), (-2, 2), (-1, -2), (0, -2), (1, -2), (-1, 2), (0, 2), (1, 2)]
    @test length(collect(Ring(10))) == length(Ring(10)) == 80
    @test all(map(x -> max(abs(x[1]),abs(x[2])), collect(Ring(100))) .== 100)
end



