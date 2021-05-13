```@meta
EditURL = "<unknown>/symmetric_strip.jl"
```

## Symmetric Strip Example
This example consists of a symmetric strip grating, i.e. a grating where the strip width
is half the unit cell width.  For normal incidence there is a closed-form solution due to Weinstein,
but a more recent reference is Problem 10.6 of Collin, *Field Theory of Guided Waves, Second Ed.*,
IEEE Press, 1991.  Here is the code for computing the exact solution:

```@repl 1
"""
    grating(kP, nterms=30) -> (Γ, T)

Compute the normal incidence reflecton and transmission coefficients of a symmetric grid of
zero-thickness conducting strips.  The product of the period of the strips and the incident
electric field wavenumber is `kP` (dimensionless).  The incident electric field is perpendicular
to the direction along the axis of the strips.  The series have been accelerated
by applying a Kummer's transformation, using the first two terms in the Maclaurin series for the
inverse sin function.  `kP` must be in the half-open interval [0,1). The default number of summed
terms `nterms` yields better than 10 digits of accuracy over the interval [0.01,0.99].
"""
function grating(kP; nterms=30)
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
    T = 1 + Γ
    return (Γ, T)
end
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

