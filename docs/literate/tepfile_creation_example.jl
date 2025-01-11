#nb # %% A slide [markdown] {"slideshow": {"slide_type": "Slide"}}
# ## TEP File Creation
# Here we show how to create a TICRA-compatible TEP (tabulated electrical properties) file using the [`res2tep`](@ref)
# function included with PSSFSS. The geometry for this example is a rectangular copper strip measuring 
# 4 cm × 0.2 cm in a 5 cm square unit cell. The code for analyzing this geometry and creating the TEP file is shown below:
# ```Julia
# using PSSFSS
# FGHz = 3.0
# Px = Py = 5
# Lx = 4.0
# Ly = 0.2
# Nx = 130
# Ny = 6
# sheet = rectstrip(; Px, Py, Lx, Ly, Nx, Ny, units=cm, sigma=5.7e8)
# steering = (θ=0:5:70, ϕ=0:15:345)
# strata = [Layer(), sheet, Layer()]
# resultfile = "strip.res"
# results = analyze(strata, FGHz, steering; resultfile)
# res2tep(results, "dipole_pssfss.tep"; name = "dipole", class = "pssfss")
# # Alternatively: res2tep(resultfile, "dipole_pssfss.tep"; name = "dipole", class = "pssfss")
# ```
# Note that a very fine discretization has been specified, resulting in 1560 triangles.  There are also a large
# number of steering angles requested (15 × 24 = 360).  This analysis required about 155 seconds on my machine.
# Please see the documentation for [`res2tep`](@ref) for requirements on setting up analysis scan angles.
#
# For comparison, the same geometry was analyzed using the QUPES program of Ticra Tools Student Edition
# 2024 (hence the specification for `sigma` above, which is the default conductivity used by QUPES). For the 
# QUPES analysis, basis function expansion accuracy was set to "Extreme" with the other analysis settings at their 
# default values.  The maximum magnitude of the difference between QUPES and PSSFSS complex scattering coefficients
# for any of these 360 steering angles was approximately 0.0021.  Convergence studies showed that for both codes, the 
# complex scattering coefficients were still varying slightly in the third decimal place for the settings used in this
# example.

