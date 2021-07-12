```@meta
EditURL = "https://github.com/simonp0420/PSSFSS.jl/tree/main/docs/literate/examples.jl"
```

# [PSSFSS](https://github.com/simonp0420/PSSFSS) Examples

```@meta
EditURL = "https://github.com/simonp0420/PSSFSS.jl/tree/main/docs/literate/symmetric_strip.jl"
```

## Symmetric Strip Grating
This example consists of a symmetric strip grating, i.e. a grating where the strip width
is half the unit cell period ``P``:

![diagram](./assets/symmetric_strip_diagram.png)

Only three of the infinite number of strips in the grating are shown, and they extend infinitely to the left and right.
The grating lies in the ``z=0`` plane with free space on both sides. The shaded areas represent metalization.
The dashed lines show two possible choices for the unit cell location: "J" for a formulation in terms of electric
surface currents, and "M" for magnetic surface currents.

For normal incidence there is a closed-form solution due to Weinstein,
but for a more recent reference one can find the solution in Problem 10.6 of R. E. Collin,
*Field Theory of Guided Waves, Second Ed.*,
IEEE Press, 1991.  Here is the code for computing the exact solution:

````@example symmetric_strip
"""
    grating(kP, nterms=30) -> (Γ, T)

Compute the normal incidence reflecton and transmission coefficients of a symmetric grid of
zero-thickness conducting strips.  The product of the period of the strips and the incident
electric field wavenumber is `kP` (dimensionless).  The incident electric field is
perpendicular to the direction along the axis of the strips.  The series have been
accelerated by applying a Kummer's transformation, using the first two terms in the Maclaurin
series for the inverse sin function.  `kP` must be in the half-open interval [0,1). The
default number of summed terms `nterms` yields better than 10 digits of accuracy over the
interval [0.01,0.99].
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
````

 Note that using the extension of
[Babinet's Principle for electromagnetic fields](http://kirkmcd.princeton.edu/examples/babinet.pdf)
this also provides the solution (upon appropriate interchange and sign change of the coefficients) for
the case where the incident wave polarization is parallel to the direction of the strips.

Here is the PSSFSS code to analyze this structure using electric currents as the unknowns.  We will
scale the geometry so that the frequency in GHz is numerically equal to the period of the strips
measured in wavelengths.

````@example symmetric_strip
using Plots, PSSFSS
c = 11.802852677165355 # light speed [inch*GHz]
period = c  # so the period/wavelength = freq in GHz
Py = period
Ly = period/2
Px = Lx = Ly/10 # Infinite in x direction so this can be anything
Ny = 60
Nx = round(Int, Ny*Lx/Ly)
sheet = rectstrip(;Px, Py, Lx, Ly, Nx, Ny, units=inch)
flist = 0.02:0.02:0.98
steering = (θ=0, ϕ=0)
strata = [Layer()
          sheet
          Layer()]
results_j = analyze(strata, flist, steering, showprogress=false,
                    resultfile=devnull, logfile=devnull);
p1 = plot(sheet)
p2 = plot(sheet, unitcell=true)
ptitle = plot(title = "Symmetric Strip Triangulation",
             grid = false, showaxis = false, xtick=[], ytick=[],
             bottom_margin = -50Plots.px)
plot(ptitle, p1, p2, layout = @layout([A{0.09h}; [B C]]))
savefig("symstrip1.png"); nothing  # hide
````

![](symstrip1.png)
