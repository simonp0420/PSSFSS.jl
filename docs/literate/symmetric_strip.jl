#nb # %% A slide [markdown] {"slideshow": {"slide_type": "subslide"}}
# ## Symmetric Strip Grating
# This example consists of a symmetric strip grating, i.e. a grating where the strip width
# is half the unit cell period ``P``:
#
# ![diagram](../src/assets/symmetric_strip_diagram.png)
#
# Only three of the infinite number of strips in the grating are shown, and they extend infinitely to the left and right.
# The grating lies in the ``z=0`` plane with free space on both sides. The shaded areas represent metalization.
# The dashed lines show two possible choices for the unit cell location: "J" for a formulation in terms of electric 
# surface currents, and "M" for magnetic surface currents.
#
# For normal incidence there is a closed-form solution due to Weinstein, 
# but for a more recent reference one can find the solution in Problem 10.6 of R. E. Collin, 
# *Field Theory of Guided Waves, Second Ed.*, 
# IEEE Press, 1991.  Here is the code for computing the exact solution:
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

#  Note that using the extension of 
# [Babinet's Principle for electromagnetic fields](http://kirkmcd.princeton.edu/examples/babinet.pdf)
# this also provides the solution (upon appropriate interchange and sign change of the coefficients) for 
# the case where the incident wave polarization is parallel to the direction of the strips.

# Here is the PSSFSS code to analyze this structure using electric currents as the unknowns.  We will 
# scale the geometry so that the frequency in GHz is numerically equal to the period of the strips
# measured in wavelengths.
using PSSFSS, Plots
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
results_j = analyze(strata, flist, steering, showprogress=false)
p1 = plot(sheet)
p2 = plot(sheet, unitcell=true)
title = plot(title = "Symmetric Strip Triangulation", 
             grid = false, showaxis = false, xtick=[], ytick=[],
             bottom_margin = -50Plots.px)
plot(title, p1, p2, layout = @layout([A{0.09h}; [B C]]))
# Note that setting `Lx = Px` causes the strip to fully occupy the x-extent
# of the unit cell.  PSSFSS automatically ensures that the triangle edges at these unit
# cell boundaries define basis functions that satisfy the Floquet (phase shift) boundary
# conditions, so that currents are free to flow across these unit cell boundaries.

# We can also analyze the same structure using magnetic currents in the areas free of 
# metallization as the unknowns:
sheet = rectstrip(;class='M', Px, Py, Lx, Ly, Nx, Ny, units=inch)
strata = [Layer()
          sheet
          Layer()]
results_m = analyze(strata, flist, steering, showprogress=false);
# 
# Each 50-frequency run of `analyze` takes about 14 seconds
# for this geometry of 720 triangles on my machine.  
# More detailed timing information is available in the log file.

# We will compare the PSSFSS results to the analytic solution:
## Generate exact results:
rt = grating.(2π*flist)
rperp_exact = first.(rt)
tperp_exact = last.(rt)
rpar_exact = -tperp_exact
tpar_exact = -rperp_exact;

# Obtain PSSFSS results for electric and magnetic currents:
outrequest = @outputs s11(v,v) s21(v,v) s11(h,h) s21(h,h)
rperp_j, tperp_j, rpar_j, tpar_j = 
      collect.(eachcol(extract_result(results_j, outrequest)))
rperp_m, tperp_m, rpar_m, tpar_m = 
      collect.(eachcol(extract_result(results_m, outrequest)));

# Generate the comparison plots:
angdeg(z) = rad2deg(angle(z)) # Convenience function

p1 = plot(title = "Perpendicular Reflection Magnitude",
          xlabel = "Period (wavelengths)",
          ylabel = "Coefficient Magnitude",
          legend=:topleft)
plot!(p1, flist, abs.(rperp_exact), ls=:dash, label="Exact")
plot!(p1, flist, abs.(rperp_j), label="PSSFSS J")
plot!(p1, flist, abs.(rperp_m), label="PSSFSS M")
#-
p2 = plot(title = "Perpendicular Reflection Phase",
          xlabel = "Period (wavelengths)",
          ylabel = "Phase (deg)")
plot!(p2, flist, angdeg.(rperp_exact), ls=:dash, label="Exact")
plot!(p2, flist, angdeg.(rperp_j), label="PSSFSS J")
plot!(p2, flist, angdeg.(rperp_m), label="PSSFSS M")
#-
p1 = plot(title = "Parallel Reflection Magnitude",
          xlabel = "Period (wavelengths)",
          ylabel = "Coefficient Magnitude")
plot!(p1, flist, abs.(rpar_exact), ls=:dash, label="Exact")
plot!(p1, flist, abs.(rpar_j), label="PSSFSS J")
plot!(p1, flist, abs.(rpar_m), label="PSSFSS M")
#-
p2 = plot(title = "Parallel Reflection Phase",
          xlabel = "Period (wavelengths)",
          ylabel = "Phase (deg)")
plot!(p2, flist, angdeg.(rpar_exact), ls=:dash, label="Exact")
plot!(p2, flist, angdeg.(rpar_j), label="PSSFSS J")
plot!(p2, flist, angdeg.(rpar_m), label="PSSFSS M")

# Now look at the transmission coefficients:
p1 = plot(title = "Perpendicular Transmission Magnitude",
          xlabel = "Period (wavelengths)",
          ylabel = "Coefficient Magnitude")
plot!(p1, flist, abs.(tperp_exact), ls=:dash, label="Exact")
plot!(p1, flist, abs.(tperp_j), label="PSSFSS J")
plot!(p1, flist, abs.(tperp_m), label="PSSFSS M")
#-
p2 = plot(title = "Perpendicular Transmission Phase",
          xlabel = "Period (wavelengths)",
          ylabel = "Phase (deg)")
plot!(p2, flist, angdeg.(tperp_exact), ls=:dash, label="Exact")
plot!(p2, flist, angdeg.(tperp_j), label="PSSFSS J")
plot!(p2, flist, angdeg.(tperp_m), label="PSSFSS M")
#-
p1 = plot(title = "Parallel Transmission Magnitude",
          xlabel = "Period (wavelengths)",
          ylabel = "Coefficient Magnitude", legend=:topleft)
plot!(p1, flist, abs.(tpar_exact), ls=:dash, label="Exact")
plot!(p1, flist, abs.(tpar_j), label="PSSFSS J")
plot!(p1, flist, abs.(tpar_m), label="PSSFSS M")
#-
p2 = plot(title = "Parallel Transmission Phase",
          xlabel = "Period (wavelengths)",
          ylabel = "Phase (deg)")
plot!(p2, flist, angdeg.(tpar_exact), ls=:dash, label="Exact")
plot!(p2, flist, angdeg.(tpar_j), label="PSSFSS J")
plot!(p2, flist, angdeg.(tpar_m), label="PSSFSS M")

# ### Conclusion
# Although good agreement is obtained, as expected the best agreement between 
# all three results occurs for the lowest frequencies, where the triangles are
# smallest in terms of wavelength.  This emphasizes the fact that it is necessary for the 
# user to check that enough triangles have been requested for good convergence
# over the frequency band of interest.  This example is an extremely demanding case
# in terms of bandwidth, as the ratio of maximum to minimum frequency here 
# is ``0.98/0.02 = 49:1``
