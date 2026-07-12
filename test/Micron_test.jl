# Tests that μm functions as a unit

using PSSFSS
using Test


@testset "micron" begin
    @test 1mm == 1000μm
    @test 1mm == 1000micron
end
