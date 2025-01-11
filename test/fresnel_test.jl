using PSSFSS
using Test
@testset "TEP Creation" begin
dwidth = 3mm
duroid = Layer(epsr=2.2, tandel=0.0009, width=dwidth)
strata = [Layer(), duroid, Layer(width=-dwidth)]
steering = (θ=0:10:80, ϕ=0)
FGHz = 1:2
rttbl = tempname()
resultfile = tempname()
logfile = devnull
showprogress = false
analyze(strata, FGHz, steering; logfile, resultfile, showprogress)
testfile_fresnel = tempname()
res2fresnel(resultfile, testfile_fresnel)
testdat = readlines(testfile_fresnel)
testdat = testdat[10:end]
gooddat = readlines(joinpath(@__DIR__, "pssfss_diel_fresnel_table.rttbl"))
gooddat = gooddat[10:end]
@test all(x -> isequal(x[1],x[2]), zip(testdat, gooddat))
end
