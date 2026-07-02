import PSSFSS.Outputs.read_result_file
import PSSFSS.FastSweep.interpolate_band
using StaticArrays: MArray
using LinearAlgebra: norm
using Test

@testset "Fast Sweep" begin
    direct_results = read_result_file(joinpath(@__DIR__, "wernercross_direct.res"))
    freqs = [dr.FGHz for dr in direct_results]
    smat4x4s = interpolate_band(freqs) do FGHz
        i = findfirst(x -> x.FGHz == FGHz, direct_results)
        gsm = direct_results[i].gsm
        smat4x4 =  MArray{Tuple{4,4}}([gsm[1,1] gsm[1,2]; gsm[2,1] gsm[2,2]])
        return smat4x4
    end
    err_limit_db = -80
    err_limit = 10^(err_limit_db / 20)
    for (s4x4interp, dr) in zip(smat4x4s, direct_results)
        gsm = dr.gsm
        s4x4direct =  MArray{Tuple{4,4}}([gsm[1,1] gsm[1,2]; gsm[2,1] gsm[2,2]])
        @test norm(s4x4interp - s4x4direct) ≤ err_limit
    end
end
nothing
