using PSSFSS
using Test
@testset "TEP Creation" begin
    using PSSFSS
    using TicraUtilities: TicraUtilities as TU
    FGHz = 3.0
    Px = Py = 5.0
    Lx = 4.0
    Ly = 0.2
    Nx = 130
    Ny = 6
    sheet = rectstrip(; Px, Py, Lx, Ly, Nx, Ny, units=cm, sigma=5.7e8)
    steering = (θ = 0:30:60, ϕ = 0:60:300)
    strata = [Layer(), sheet, Layer()]
    resultfile = tempname()
    logfile = devnull
    showprogress = false
    results = analyze(strata, FGHz, steering; resultfile, logfile, showprogress)
    teppssfss_file = tempname()
    res2tep(resultfile, teppssfss_file)
    tep_pssfss = TU.read_tepfile(teppssfss_file)
    tep_ticra = TU.read_tepfile(joinpath(pkgdir(PSSFSS), "test", "ticra_tools_dipole.tep"))
    @test maximum(abs, TU.get_sff(tep_pssfss) - TU.get_sff(tep_ticra)) < 0.001
    @test maximum(abs, TU.get_sfr(tep_pssfss) - TU.get_sfr(tep_ticra)) < 0.001
    @test maximum(abs, TU.get_srf(tep_pssfss) - TU.get_srf(tep_ticra)) < 0.001
    @test maximum(abs, TU.get_srr(tep_pssfss) - TU.get_srr(tep_ticra)) < 0.001
end
