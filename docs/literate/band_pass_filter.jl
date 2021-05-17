#nb # %% A slide [markdown] {"slideshow": {"slide_type": "subslide"}}
# ## Loaded Cross Band Pass Filter
# This example is originally from Fig. 7.9 of B. Munk, *Frequency Selective Surfaces,
# Theory and Design,* John Wiley and Sons, 2000.  The same case was analyzed in L. Li,
# D. H. Werner et al, "A Model-Based Parameter Estimation Technique for
# Wide-Band Interpolation of Periodic Moment Method Impedance Matrices With Application to
# Genetic Algorithm Optimization of Frequency Selective Surfaces", *IEEE Trans. AP-S*,
# vol. 54, no. 3, March 2006, pp. 908-924, Fig. 6.  Unfortunately, in neither reference 
# are the dimensions of the loaded cross stated, except for the square unit cell
# period of 8.61 mm.  I estimated the dimensions from the sketch in Fig. 6 of the second
# reference.  To provide a reliable comparison, I enlisted my colleague 
# [Mike Maybell](https://www.linkedin.com/in/mike-maybell-308b77ba),
# principal of Planet Earth Communications, who generously offered to
# analyze the filter using 
# [CST Microwave Studio](https://www.3ds.com/products-services/simulia/products/cst-studio-suite/), 
# a rigorous commercial finite volume electromagnetic solver.

# Two identical loaded cross slot-type elements are separated by a 6 mm layer of dielectric
# constant 1.9.  Outboard of each sheet is a 1.1 cm layer of dielectric constant 1.3.
# The closely spaced sheets are a good test of the generalized scattering formulation
# implemented in PSSFSS.  The sheet geometry is shown below.  Remember that the entire
# sheet is metalized *except* for the region of the triangulation.

#md ENV["GKSwstype"] = "100" # hide
using Plots, PSSFSS
sheet = loadedcross(class='M', w=0.023, L1=0.8, L2=0.14, 
            s1=[0.861,0.0], s2=[0.0,0.861], ntri=600, units=cm)
plot(sheet, unitcell=true)

#-

steering = (ϕ=0, θ=0)
strata = [  Layer()
            Layer(ϵᵣ=1.3, width=1.1cm)
            sheet
            Layer(ϵᵣ=1.9, width=0.6cm)
            sheet
            Layer(ϵᵣ=1.3, width=1.1cm)
            Layer()  ]
flist = 1:0.1:20
results = analyze(strata, flist, steering, resultfile=devnull, 
                  logfile=devnull, showprogress=false)
data = extract_result(results, @outputs FGHz s21db(v,v) s11db(v,v))
using DelimitedFiles
dat = readdlm("../src/assets/MaybellLoadedCrossResults.csv", ',', skipstart=1)
p = plot(xlabel="Frequency (GHz)", ylabel="Reflection Coefficient (dB)",
         legend=:left, title="Loaded Cross Band-Pass Filter", xtick=0:2:20, ytick=-30:5:0,
         xlim=(-0.1,20.1), ylim=(-35,0.1))
plot!(p, data[:,1], data[:,3], label="PSSFSS", color=:red)
plot!(p, dat[:,1], dat[:,2], label="CST", color=:blue)
p
#-
p2 = plot(xlabel="Frequency (GHz)", ylabel="Transmission Coefficient (dB)",
          legend=:bottom, title="Loaded Cross Band-Pass Filter", xtick=0:2:20, ytick=-80:10:0,
         xlim=(-0.1,20.1), ylim=(-80,0.1))
plot!(p2, data[:,1], data[:,2], label="PSSFSS", color=:red)
plot!(p2, dat[:,1], dat[:,4], label="CST", color=:blue)
p2

# This analysis takes about 90 seconds for 191 frequencies on my machine.  Note that
# rather than including two separate invocations of the `loadedcross` function when
# defining the strata, I referenced the same sheet object in the two different locations.
# This allows PSSFSS to recognize that the triangulations are identical, and to exploit
# this fact in making the analysis more efficient.  In fact, if both sheets had been embedded
# in similar dielectric claddings (in the same order), then the GSM (generalized scattering matrix)
# computed for the sheet in its first location could be reused without additional computation for its
# second location.  In this case, though, only the spatial integrals are re-used.  For a oblique
# incidence case, computing the spatial integrals is often the most expensive part of the analysis,
# so the savings from reusing the same sheet definition can be substantial.

# ### Conclusion
# Very good agreement is obtained versus CST over a large dynamic range.