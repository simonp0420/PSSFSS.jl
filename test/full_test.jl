using PSSFSS
using PSSFSS.Constants: c₀
using PSSFSS.GSMs: GSM
using LinearAlgebra: norm
using Test
using Logging: Error, ConsoleLogger, default_metafmt, global_logger


testlogger = ConsoleLogger(stderr, Error,
                       meta_formatter=default_metafmt, show_limited=true,
                       right_justify=0)
oldlogger = global_logger(testlogger)


"""
    grating(kP, nterms=30) -> (Γ, T)

Compute the normal incidence reflecton and transmission coefficients of a symmetric grid of 
zero-thickness conducting strips.  The product of the period of the strips and the incident
electric field wavenumber is `kP` (dimensionless).  The incident electric field is perpendicular
to the direction along the axis of the strips.  The formula is from Problem 10.6 of Collin,
*Field Theory of Guided Waves, Second Ed.*, IEEE Press, 1991.  The series have been accelerated
by appying a Kummer transformation, using the first two terms in the Maclaurin series for the
inverse sin function.  `kP` must be in the half-open interval [0,1). The default number of summed
terms `nterms` yeilds better than 10 digits of accuracy over the interval [0.01,0.99].
"""
function grating(kP;nterms=30)
    sum1 = 1.3862943611198906 # \sum_{n=1}^{\infty} 1/(n-1/2) - 1/n = log(4) 
    sum3 = 7.2123414189575710 # \sum_{n=1}^{\infty} (n-1/2)^{-3} - n^{-3} = 6 * \zeta(3)
    x = kP/(4π)
    θ = x*sum1 + x^3/6 * sum3
    for n = 1:nterms
        xonmhalf = x/(n - 0.5)
        xon = x/n
        term = asin(xonmhalf) - (xonmhalf + (xonmhalf)^3/6) - 
              (asin(xon) - (xon + xon^3/6))
        θ += term
    end
    Γ = sin(θ) * cis(-π/2 - θ)
    T = 1+Γ
    return (Γ, T)
end


FGHz = 0.8
steering = (θ=0, ϕ=0)
## Compute analytic answer
(r,t) = grating(2π*FGHz)
gsmexact = GSM(2,2)
gsmexact.s11 = [r 0.0; 0.0 -t]
gsmexact.s22 = gsmexact.s11
gsmexact.s21 = [t 0; 0 -r]
gsmexact.s12 = gsmexact.s21

c_inchghz = c₀ * 100/2.54 * 1e-9
period = c_inchghz # so the period/wavelength = freq in GHz
Py = period
Ly = period/2
Px = Lx = Ly/10
Ny = 40
Nx = round(Int, Ny*Lx/Ly)



@testset "JtypeSymmetricStrip" begin

    sheet = rectstrip(;Px, Py, Lx, Ly, Nx, Ny, units=inch)
    strata = [Layer()
              sheet
              Layer()]

    results = analyze(strata, FGHz, steering, logfile=devnull, 
                      resultfile=devnull, showprogress=false)
    gsmj = results[1].gsm
    for m in 1:2, n in 1:2
        @test norm(gsmj[m,n] - gsmexact[m,n], Inf) < 0.01
    end
end    

@testset "MtypeSymmetricStrip" begin

    sheet = rectstrip(;Px, Py, Lx, Ly, Nx, Ny, units=inch, class='M')
    strata = [Layer()
              sheet
              Layer() ]     

    results = analyze(strata, FGHz, steering, logfile=devnull, 
                      resultfile=devnull, showprogress=false)
    gsmm = results[1].gsm
    for m in 1:2, n in 1:2
        @test  norm(gsmm[m,n] - gsmexact[m,n], Inf) < 0.01
    end
end

global_logger(oldlogger)
nothing
