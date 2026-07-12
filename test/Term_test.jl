# Tests that μm functions as a unit

using PSSFSS: PSSFSS
using Test


@testset "term" begin
    tests = (PSSFSS._is_ijulia(), PSSFSS._is_vsc_integrated_term(), PSSFSS._is_vsc_notebook(),
             PSSFSS._is_pluto_notebook())
    @test count(tests) < 2
    @test PSSFSS.run_environment() isa String
end
