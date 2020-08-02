using Test
using SafeTestsets

@safetestset "RWGSheet Tests" begin include("RWGSheet_test.jl") end
@safetestset "Elements Tests" begin include("Elements_test.jl") end
@safetestset "RWGData Tests" begin include("RWGData_test.jl") end
@safetestset "PGF Tests" begin include("PGF_test.jl") end
