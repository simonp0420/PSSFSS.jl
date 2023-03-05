#nb # %% A slide [markdown] {"slideshow": {"slide_type": "Slide"}}
# ## Reflectarray Element
# This example is taken from Figure 6 of 
# Li, Jiao, and Zhao: "A Novel Microstrip Rectangular-Patch/Ring-
# Combination Reflectarray Element and Its Application", **IEEE Antennas and Wireless Propagation Letters**, 
# VOL. 8, 2009, pp. 1119-1112.
#
# It generates the so-called "S-curve" for reflection phase of a reflectarray element. The element 
# consists of a square patch in a square ring, separated from a ground plane by two dielectric layers.
# Reflection phase is plotted versus the `L2` parameter as defined by Li, et al., which characterizes the 
# overall size of the element.
#

# We start by defining a convenience function to generate a `RWGSheet` for a given value of `L2` in mm and 
# desired number of triangles `ntri`:

function element(L2, ntri::Int)
    L = 17 # Period of square lattice (mm)
    k1 = 0.5
    k2 = 0.125
    L1 = k1 * L2
    w = k2 * L2
    a = [-1.0, (L2-2w)/√2]
    b = [L1/√2, L2/√2]
    return polyring(; sides=4, s1=[L, 0], s2=[0, L], ntri, orient=45, a, b, units=mm)
end

# Here is a plot of the two elements at the size extremes to be examined:
using PSSFSS, Plots
p1 = plot(element(2, 675), title="L2 = 2mm", unitcell=true, linecolor=:red)
p2 = plot(element(16.5, 3416), title="L2 = 16.5mm", unitcell=true, linecolor=:red)
p = plot(p1,p2)
#md savefig(p, "reflectarray_elements.png"); nothing  # hide
#-
#md # ![](reflectarray_elements.png)

# Here is the script that generates the S-curve.  To ensure accurate phases, a convergence check
# is performed for each distinct `L2` value.  The number of triangles is increased by a factor of 1.5
# if the previous increase resulted in a reflection phase change of more than one degree.  After
# computing the phases, comparison is made to the same data computed by Li, et al from Ansoft HFSS and 
# CST Microwave Studio, and displayed in their Figure 6.  

# Note that the `width` for the first `Layer` is set to `-7mm`.  This has the effect of referring
# reflection phases to the location of the ground plane, as was apparently done by Li, et al.  Also
# note that the `unwrap!` function of the `DSP` package is used to "unwrap" phases to remove any
# possible discontinuous jumps of 360 degrees.

# ```Julia
# using PSSFSS
# using Plots, DelimitedFiles
# using DSP: unwrap!
#
# p = plot(title="Reflection Phase Convergence", xlim=(2,18),ylim=(-350,150),xtick=2:2:18,ytick=-350:50:150,
#          xlabel="Large Square Dimension L2 (mm)", ylabel="Reflecton Phase (deg)", legend=:topright)
#
# FGHz = 10.0
# L2s = range(2.0, 16.5, length=30)
# nfactor = 1.5 # Growth factor for ntri
# phasetol = 1.0 # Phase convergence tolerance
# record = Tuple{Float64, Int, Int}[] # Storage for (L2,ntri,facecount(sheet))
# phases = zeros(length(L2s))
# ntri = 300
# for (il2, L2) in pairs(L2s)
#     print("L2 = $L2; ntri =")
#     if il2 > 1
#         global ntri = round(Int, ntri / (nfactor*nfactor))
#     end
#     oldphase = Inf
#     firstrun = true
#     while true
#         ntri = round(Int, nfactor * ntri)
#         if firstrun
#             print(" $ntri")
#             firstrun = false
#         else
#             print(", $ntri")
#         end
#         sheet = element(L2, ntri)
#         strata = [Layer(width=-7mm)
#                   sheet
#                   Layer(ϵᵣ=2.65, width=1mm)
#                   Layer(ϵᵣ=1.07, width=6mm)
#                   pecsheet() # ground plane
#                   Layer()]
#         results = analyze(strata, FGHz, (ϕ=0, θ=0), showprogress=false, logfile="li2009.log", resultfile="li2009.res")
#         phase = extract_result(results, @outputs s11ang(h,h))[1,1]
#         phase < 0 && (phase += 360)
#         Δphase = abs(phase - oldphase)
#         if Δphase > phasetol
#             oldphase = phase
#         else
#             phases[il2] = phase
#             push!(record, (L2, ntri, facecount(sheet)))
#             println("; Δphase = $(round(Δphase, digits=2))°")
#             break
#         end
#     end
# end
# 
# unwrap!(phases, range=360)
# p = plot(title="Li 2009 Reflectarray", xlim=(2,18),ylim=(-350,150),xtick=2:2:18,ytick=-350:50:150,
#     xminorticks=2, yminorticks=2,
#     xlabel="Large Square Dimension L2 (mm)", ylabel="Reflecton Phase (deg)", legend=:topright)
# plot!(p, L2s, phases, label="PSSFSS", color=:red, mscolor=:red, ms=3, shape=:circ)
# dat = readdlm("li2009_ansoft_digitized.csv",  ',')
# scatter!(p, dat[:,1], dat[:,2], label="Li Ansoft", mc=:blue, msc=:blue, markershape=:square)
# dat = readdlm("li2009_cst_digitized.csv",  ',')
# scatter!(p, dat[:,1], dat[:,2], label="Li CST", mc=:green, msc=:green, markershape=:circle)
# display(p)
# println()
# record
# ```

# ![](./assets/li2009_comparison.svg)

# It can be seen that the PSSFSS phases compare well with the HFSS phases, better than the CST and HFSS phases compare
# to each other.  The authors of the paper do not discuss checking the convergence of their results or even any details of
# how they set up their HFSS and CST models.

# The console output from the above script is shown below:
# ```Julia
# julia> include("li2009_convergence.jl")
# 
# L2 = 2.0; ntri = 450, 675; Δphase = 0.02°
# L2 = 2.5; ntri = 450, 675; Δphase = 0.04°
# L2 = 3.0; ntri = 450, 675; Δphase = 0.06°
# L2 = 3.5; ntri = 450, 675; Δphase = 0.1°
# L2 = 4.0; ntri = 450, 675; Δphase = 0.14°
# L2 = 4.5; ntri = 450, 675; Δphase = 0.18°
# L2 = 5.0; ntri = 450, 675; Δphase = 0.19°
# L2 = 5.5; ntri = 450, 675; Δphase = 0.15°
# L2 = 6.0; ntri = 450, 675; Δphase = 0.02°
# L2 = 6.5; ntri = 450, 675; Δphase = 0.41°
# L2 = 7.0; ntri = 450, 675, 1012; Δphase = 0.12°
# L2 = 7.5; ntri = 675, 1012; Δphase = 0.23°
# L2 = 8.0; ntri = 675, 1012; Δphase = 0.38°
# L2 = 8.5; ntri = 675, 1012; Δphase = 0.55°
# L2 = 9.0; ntri = 675, 1012; Δphase = 0.73°
# L2 = 9.5; ntri = 675, 1012; Δphase = 0.93°
# L2 = 10.0; ntri = 675, 1012, 1518, 2277, 3416; Δphase = 0.67°
# L2 = 10.5; ntri = 2277, 3416; Δphase = 0.54°
# L2 = 11.0; ntri = 2277, 3416; Δphase = 0.35°
# L2 = 11.5; ntri = 2277, 3416; Δphase = 0.1°
# L2 = 12.0; ntri = 2277, 3416; Δphase = 0.17°
# L2 = 12.5; ntri = 2277, 3416; Δphase = 0.41°
# L2 = 13.0; ntri = 2277, 3416; Δphase = 0.59°
# L2 = 13.5; ntri = 2277, 3416; Δphase = 0.69°
# L2 = 14.0; ntri = 2277, 3416; Δphase = 0.73°
# L2 = 14.5; ntri = 2277, 3416; Δphase = 0.72°
# L2 = 15.0; ntri = 2277, 3416; Δphase = 0.69°
# L2 = 15.5; ntri = 2277, 3416; Δphase = 0.65°
# L2 = 16.0; ntri = 2277, 3416; Δphase = 0.61°
# L2 = 16.5; ntri = 2277, 3416; Δphase = 0.59°
# 
# 30-element Vector{Tuple{Float64, Int64, Int64}}:
#  (2.0, 675, 722)
#  (2.5, 675, 722)
#  (3.0, 675, 722)
#  (3.5, 675, 722)
#  (4.0, 675, 722)
#  (4.5, 675, 722)
#  (5.0, 675, 722)
#  (5.5, 675, 722)
#  (6.0, 675, 722)
#  (6.5, 675, 722)
#  (7.0, 1012, 944)
#  (7.5, 1012, 944)
#  (8.0, 1012, 944)
#  (8.5, 1012, 944)
#  (9.0, 1012, 944)
#  (9.5, 1012, 944)
#  (10.0, 3416, 3314)
#  (10.5, 3416, 3314)
#  (11.0, 3416, 3314)
#  (11.5, 3416, 3314)
#  (12.0, 3416, 3314)
#  (12.5, 3416, 3314)
#  (13.0, 3416, 3314)
#  (13.5, 3416, 3314)
#  (14.0, 3416, 3314)
#  (14.5, 3416, 3314)
#  (15.0, 3416, 3314)
#  (15.5, 3416, 3314)
#  (16.0, 3416, 3314)
#  (16.5, 3416, 3314)
# ```

# The first list shows how the script increases the number of triangles requested for each geometry
# until the reflection phase is sufficiently converged.  The printout of the `record` array shows
# both the final requested value of `ntri` and actual number of triangle faces generated by the mesher for each
# `L2` value.  The final cases with `ntri=3416` required about 21 seconds each of execution time.
