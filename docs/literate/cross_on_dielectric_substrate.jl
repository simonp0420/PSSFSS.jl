#nb # %% A slide [markdown] {"slideshow": {"slide_type": "subslide"}}
# ## Cross on Dielectric Substrate
# This example is also taken from the paper by Alon S. Barlevy and 
# Yahya Rahmat-Samii, "Fundamental Constraints on the Electrical Characteristics 
# of Frequency Selective Surfaces", **Electromagnetics**, vol. 17, 1997, pp. 41-68. 
# This particular example is from Section 3.2, Figures 7 and 8.  It also appeared at 
# higher resolution in Barlevy's PhD dissertation from which the comparison curves 
# were digitized.
#
# We use the `loadedcross` element where we choose `w > L2/2`, so that the Cross
# is "unloaded", i.e. the center section is filled in with metallization:

#md ENV["GKSwstype"] = "100" # hide
using Plots, PSSFSS, DelimitedFiles
sheet = loadedcross(w=1.0, L1=0.6875, L2=0.0625, s1=[1.0,0.0], 
                    s2=[0.0,1.0], ntri=600, units=cm)
plot(sheet, unitcell=true)

# A few things to note. First, the mesh is **unstructured**.  So there are no redundant 
# triangle face-pairs that PSSFSS can exploit to reduce execution time.  Second, the 
# number of triangle faces generated is only approximately equal to the requested value
# of 600.  This can be verified by entering the Julia variable `sheet` at the 
# [REPL](https://docs.julialang.org/en/v1/manual/getting-started/#man-getting-started) 
# (i.e. the Julia prompt):

sheet

# The cross FSS is etched on a dielectric sheet of thickness 3 mm.  The dielectric 
# constant is varied over the values 1, 2, and 4 to observe the effect on the resonant 
# frequency.  Following the reference, the list of analysis frequencies is varied slightly
# depending on the value of dielectric constant:

resultsstack = Any[]
steering = (ϕ=0, θ=0)
for eps in [1, 2, 4]
    strata = [  Layer()
                sheet
                Layer(ϵᵣ=eps, width=3mm)
                Layer() 
             ]
    if eps == 1
        flist = 1:0.2:30
    elseif eps == 2
        flist = 1:0.2:26
    else
        flist = 1:0.2:20
    end
    results = analyze(strata, flist, steering, showprogress=false, 
                      resultfile=devnull, logfile=devnull)
    push!(resultsstack, results)
end

# The above loop requires about 80 seconds of execution time on my machine.
# Compare PSSFSS results to those digitized from the dissertation figure:

col=[:red,:blue,:green]
p = plot(xlim=(0.,30), xtick = 0:5:30, ylim=(0,1), ytick=0:0.1:1, 
         xlabel="Frequency (GHz)", ylabel="Reflection Coefficient Magnitude",
         legend=:topleft, lw=2)
for (i,eps) in enumerate([1,2,4])
    data = extract_result(resultsstack[i],  @outputs FGHz s11mag(v,v))
    plot!(p, data[:,1], data[:,2], label="PSSFSS ϵᵣ = $eps", lc=col[i])
    data = readdlm("../src/assets/barlevy_diss_eps$(eps).csv", ',')
    plot!(p, data[:,1], data[:,2], label="Barlevy ϵᵣ = $eps", lc=col[i], ls=:dot)
end
p


# ### Conclusion
# PSSFSS results agree very well with those of the cited reference, especially when
# accounting for the fact that the reference results are obtained by digitizing a 
# scanned figure.
# 