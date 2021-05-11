#nb # %% A slide [markdown] {"slideshow": {"slide_type": "subslide"}}
# ## Resistive Square Patch
# This example will demonstrate the ability of PSSFSS to accurately model finite
# conductivity of FSS metallization.  It consists of a square finitely conducting 
# patch in a square lattice.  It is taken from a paper by Alon S. Barlevy and 
# Yahya Rahmat-Samii, 
# "Fundamental Constraints on the Electrical Characteristics of Frequency Selective 
# Surfaces", **Electromagnetics**, vol. 17, 1997, pp. 41-68. This particular example 
# is from Section 3.2, Figures 7 and 8.  We will compare PSSFSS results to those digitized
# from the cited figures.
#
# We start by defining a function that creates a patch of the desired sheet resistance:

using PSSFSS, Plots
patch(R) = rectstrip(Nx=10, Ny=10, Px=1, Py=1, Lx=0.5, Ly=0.5, units=cm, Rsheet=R)
plot(patch(0), unitcell=true)

# The patches measure 0.5 cm on a side and lie in a square lattice of period 1 cm.
# Now we perform the analysis, looping over the desired values of sheet resistance.

steering = (ϕ=0, θ=0)
flist = 1:0.5:60
Rs = [0, 10, 30, 100]
calculated = zeros(length(flist), length(Rs)) # preallocate storage
outputs = @outputs s11mag(v,v)
for (i,R) in pairs(Rs)
    strata = [Layer(), patch(R), Layer()]
    results = analyze(strata, flist, steering, showprogress=false, 
                      logfile=devnull, resultfile=devnull)
    calculated[:,i] = extract_result(results, outputs)
end

# Looping over the four sheet resistance values, each evaluated at 119 frequencies
# required approximately 20 seconds on my machine.
#
# We plot the results, including those digitized from the paper for comparison:

using DelimitedFiles
markers = (:diamond, :utriangle, :square, :xcross)
colors = (:blue, :red, :green, :black)
p = plot(xlim=(-0.01,60.01), xtick = 0:10:60, ylim=(-0.01,1.01), ytick=0:0.1:1, 
         xlabel="Frequency (GHz)", ylabel="Reflection Coefficient Magnitude",
         title = "Resistive Square Patch",
         legend=:topright)
for (i,R) in pairs(Rs)
    scatter!(p, flist, calculated[:,i], label="PSSFSS $R Ω", ms=2, shape=markers[i], color=colors[i])
    data = readdlm("../src/assets/barlevy_patch_digitized_$(R)ohm.csv", ',')
    plot!(p, data[:,1], data[:,2], label="Barlevy $R Ω", ls=:dash, color=colors[i])
end
p

# ### Conclusion
# PSSFSS results are indistinguishable from those reported in the cited paper.
# 