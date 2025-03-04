#nb # %% A slide [markdown] {"slideshow": {"slide_type": "subslide"}}
# ## Square Loop Absorber
# This example, from Figure 7 of [costa2012frequency](@cite),
# shows how one can use the [`polyring`](@ref)
# function to model square loop elements.  Three different designs are examined
# that employ different loop thicknesses and different values of sheet resistance.
# We compare the reflection coefficient magnitudes computed by PSSFSS with those digitized
# from the cited figure when the sheet is suspended
# 5 mm above a ground plane, hence we will also make use of the [`pecsheet`](@ref) function.

using Plots, PSSFSS, DelimitedFiles
D = 11 # Period of square lattice (mm)
r_outer = √2/2 * D/8 * [5,6,7] # radii of square outer vertices
thickness = D/16 * [1,2,3]
r_inner = r_outer - √2 * thickness  # radii of square inner vertices
Rs = [15,40,70] # Sheet resistance (Ω/□)
labels = ["Thin", "Medium", "Thick"]
colors = [:green, :blue, :red]
p = plot(title="Costa Absorber", xlim=(0,25),ylim=(-35,0),xtick=0:5:25,ytick=-35:5:0,
         xlabel="Frequency (GHz)", ylabel="Reflection Magnitude (dB)", legend=:bottomleft)
ps = []
for (ri, ro, label, color, R) in zip(r_inner, r_outer, labels, colors, Rs)
    sheet = polyring(sides=4, s1=[D, 0], s2=[0, D], ntri=750, orient=45, 
                     a=[ri], b=[ro], Zsheet=R, units=mm)
    push!(ps, plot(sheet, unitcell=true, title=label, lc=color))
    strata = [Layer()
              sheet
              Layer(width=5mm)
              pecsheet() # Perfectly conducting ground plane
              Layer()]
    results = analyze(strata, 1:0.2:25, (ϕ=0, θ=0), showprogress=false,
                      resultfile=devnull, logfile=devnull)
    data = extract_result(results, @outputs FGHz s11dB(h,h))
    plot!(p, data[:,1], data[:,2], label="PSSFSS "*label, lc=color)
    dat = readdlm("../src/assets/costa_2014_" * lowercase(label) * "_reflection.csv", ',')
    plot!(p, dat[:,1], dat[:,2], label="Costa "*label, ls=:dash, lc=color)
end
plot(ps..., layout=(1,3), size=(600,220), margin=3Plots.mm)
#md savefig("sqloop1.png"); nothing  # hide
#-
#md # ![](sqloop1.png)
#-
# This PSSFSS run of three geometries takes about 15 seconds on my machine.
p
#md savefig(p,"sqloop2.png"); nothing  # hide
#-
#md # ![](sqloop2.png)

# It is useful to take a look at the log file created by PSSFSS for the last run above
# (from a previous run where the log file was not discarded):
# ```
# Starting PSSFSS 1.2.2 analysis on 2023-03-04 at 19:44:10.800
# Julia Version 1.8.5
# Commit 17cfb8e65e (2023-01-08 06:45 UTC)
# Platform Info:
#   OS: Windows (x86_64-w64-mingw32)
#   CPU: 8 × Intel(R) Core(TM) i7-9700 CPU @ 3.00GHz
#   WORD_SIZE: 64
#   LIBM: libopenlibm
#   LLVM: libLLVM-13.0.1 (ORCJIT, skylake)
#   Threads: 8 on 8 virtual cores
#   BLAS: LBTConfig([ILP64] libopenblas64_.dll)
# 
# 
# 
# Dielectric layer information... 
# 
#  Layer  Width  units  epsr   tandel   mur  mtandel modes  beta1x  beta1y  beta2x  beta2y
#  ----- ------------- ------- ------ ------- ------ ----- ------- ------- ------- -------
#      1    0.0000  mm    1.00 0.0000    1.00 0.0000     2   571.2    -0.0    -0.0   571.2
#  ==================  Sheet   1  ========================   571.2    -0.0    -0.0   571.2
#      2    5.0000  mm    1.00 0.0000    1.00 0.0000    42   571.2    -0.0    -0.0   571.2
#  ==================  Sheet   2  ========================     0.0     0.0     0.0     0.0
#      3    0.0000  mm    1.00 0.0000    1.00 0.0000     2   571.2    -0.0    -0.0   571.2
# 
# 
# 
# PSS/FSS sheet information...
# 
# Sheet  Loc         Style      Rot  J/M Faces Edges Nodes Unknowns  NUFP
# -----  ---  ---------------- ----- --- ----- ----- ----- -------- ------
#    1     1          polyring   0.0  J    720  1152   432    1008  199676
#    2     2              NULL   0.0  E      0     0     0       0       0
# ⋮
# ```
#-
# Note from the dielectric layer report that there are 42 modes defined in the region between the 
# ground plane and the FSS sheet.  This is the number of modes selected by the code to include
# in the generalized scattering matrix formulation to properly account for electromagnetic coupling
# between the two surfaces. If the 5 mm spacing were increased to, say, 7 mm then fewer modes
# would be needed.  Also note in the FSS sheet information that `NUFP` (the number of unique face pairs)
# 199676, is less than the number of faces squared (``567009 = 753^2``), a consequence of the structured
# triangulation used for a 4-sided `polyring`.  

#-
# ### Conclusion
# PSSFSS results agree very well with those of the paper, except for the medium
# width loop, where the agreement is not quite as good.  It was found empirically that
# using a slightly different value of `Rs = 37` for this ring results in nearly perfect agreement
# with the digitized results.
