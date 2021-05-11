#nb # %% A slide [markdown] {"slideshow": {"slide_type": "subslide"}}
# ## Square Loop Absorber
# This example is from Figure 7 of Costa and Monorchio: "A frequency selective 
# radome with wideband absorbing properties", *IEEE Trans. AP-S*, 
# Vol. 60, no. 6, June 2012, pp. 2740--2747.  It shows how one can use the `polyring`
# function to model square loop elements.  Three different designs are examined
# that employ different loop thicknesses and different values of sheet resistance.
# We compare the reflection coefficient magnitudes computed by PSSFSS with those digitized
# from the cited figure when the sheet is suspended
# 5 mm above a ground plane, hence we will also make use of the `pecsheet` function.

using PSSFSS, Plots, DelimitedFiles
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
for (i,(ri, ro, label, color, R)) in enumerate(zip(r_inner, r_outer, labels, colors, Rs))
    sheet = polyring(sides=4, s1=[D, 0], s2=[0, D], ntri=700, orient=45, 
                     a=[ri], b=[ro], Rsheet=R, units=mm)
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
plot(ps..., layout=(1,3))
# This run takes about 85 seconds on my machine.
p

# It is useful to take a look at the log file created by PSSFSS for the last run above:
# ```
# Starting PSSFSS analysis on 2021-05-05 at 15:22:09.955
# 
# 
# Dielectric layer information... 
# 
# Layer  Width  units  epsr   tandel   mur  mtandel modes  beta1x  beta1y  beta2x  beta2y
# ----- ------------- ------- ------ ------- ------ ----- ------- ------- ------- -------
#     1    0.0000  mm    1.00 0.0000    1.00 0.0000     2   571.2     0.0     0.0   571.2
# ==================  Sheet   1  ========================   571.2     0.0     0.0   571.2
#     2    5.0000  mm    1.00 0.0000    1.00 0.0000    42   571.2     0.0     0.0   571.2
# ==================  Sheet   2  ========================     0.0     0.0     0.0     0.0
#     3    0.0000  mm    1.00 0.0000    1.00 0.0000     2   571.2     0.0     0.0   571.2
# 
# 
# 
# PSS/FSS sheet information...
# 
# Sheet  Loc         Style      Rot  J/M Faces Edges Nodes Unknowns  NUFP
# -----  ---  ---------------- ----- --- ----- ----- ----- -------- ------
#   1     1          polyring   0.0  J    753  1201   448    1058  567009
#   2     2              NULL   0.0  E      0     0     0       0       0
# 
# ...
# ```
#-
# Note from the dielectric layer report that there are 42 modes defined in the region between the 
# ground plane and the FSS sheet.  This is the number of modes selected by the code to include
# in the generalized scattering matrix formulation to properly account for electromagnetic coupling
# between the two surfaces. If the 5 mm spacing were increased to, say, 7 mm then fewer modes
# would be needed.  Also note in the FSS sheet information that `NUFP` (the number of unique face pairs)
# is exactly equal to the number of faces squared (``567009 = 753^2``), a consequence of the unstructured
# triangulation used for a `polyring`.

#-
# ### Conclusion
# PSSFSS results agree very well with those of the paper, except for the medium
# width loop, where the agreement is not quite as good.  The reason for this is
# not known.