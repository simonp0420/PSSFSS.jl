using Test
using SafeTestsets

@safetestset "Term Tests" begin
    include("Term_test.jl")
end
@safetestset "Rings Tests" begin
    include("Rings_test.jl")
end
@safetestset "RWGSheet Tests" begin
    include("RWGSheet_test.jl")
end
@safetestset "Elements Tests" begin
    include("Elements_test.jl")
end
@safetestset "RWGData Tests" begin
    include("RWGData_test.jl")
end
@safetestset "PGF Tests" begin
    include("PGF_test.jl")
end
@safetestset "Zint Tests" begin
    include("Zint_test.jl")
end
@safetestset "GSMs Tests" begin
    include("GSMs_test.jl")
end
@safetestset "Modes Tests" begin
    include("Modes_test.jl")
end
@safetestset "FillZY Tests" begin
    include("FillZY_test.jl")
end
@safetestset "FastSweep Tests" begin
    include("FastSweep_test.jl")
end
@safetestset "Full Tests" begin
    include("full_test.jl")
end
@safetestset "TEP test" begin
    include("tep_test.jl")
end
@safetestset "Fresnel test" begin
    include("fresnel_test.jl")
end

@safetestset "sympixels test" begin
    include("sympixels_test.jl")
end

@safetestset "Micron test" begin
    include("Micron_test.jl")
end

@safetestset "Translation Test" begin
    include("Translate_test.jl")
end

if get(ENV, "BENCHMARK", "false") == "true"
    include("benchmark.jl")
end