#nb # %% A slide [markdown] {"slideshow": {"slide_type": "subslide"}}
# ## Flexible Absorber
# This example is from [li2023ultra](@cite).
# It uses square and circular resistive FSS elements sandwiched between layers of flexible dielectrics to 
# realize a reflective absorber (i.e. a "rabsorber").
# We compare the reflection coefficient magnitude computed by PSSFSS to that digitized
# from the Figure 2(a) of the paper.  The latter was obtained by the authors using CST Microwave Studio.

using PSSFSS, Plots, DelimitedFiles

p = 20 # Period of square unit cell (mm)
a, b, Rsquare = (0.5, 19, 120) # Parameters of squares
c, d, Rcircle = (1, 9, 400) # Parameters of circles

units = mm; Px = Py = p; Lx = Ly = p - 2a; Nx = Ny = 21
squares = rectstrip(; Px, Py, Lx, Ly, Nx, Ny, units, Zsheet=Rsquare)
s1 = [p,0]; s2 = [0, p]; ntri = 800; sides = 40
disks = polyring(; units, sides, a=[0], b=[d], s1, s2, ntri, Zsheet=Rcircle)

psq = plot(squares, unitcell=true, linecolor=:red)
pdi = plot(disks, unitcell=true, linecolor=:red)
psheet = plot(psq, pdi)
plot!(psheet, title=["Resistive Squares" "Resistive Disks"], size=(600,300))
#md savefig(psheet,"flexibleabsorbersheets.png"); nothing  # hide
#-
#md # ![](flexibleabsorbersheets.png)
#-
silicone = Layer(width = 4mm, ϵᵣ = 2.9, tanδ = 0.1) # using unicode keywords
foam = Layer(width = 6mm, epsr = 1.05, tandel = 0.02) # using ASCII keywords

strata = [Layer(), silicone, disks, silicone, squares, foam, pecsheet(), Layer()]

FGHz = 0.5:0.1:20
scan = (θ=0, ϕ=0)

results = analyze(strata, FGHz, scan, showprogress=false, resultfile=devnull, logfile=devnull)

s11db = extract_result(results, @outputs s11db(te,te))
li = readdlm("../src/assets/li2023_fig2a_s11db_digitized.csv", ',', Float64, '\n')

ps11 = plot(title="Normal Incidence Reflection Magnitude", 
            xlabel="Frequency (GHz)", ylabel="20log₁₀|s₁₁|", xlim=(0,20),
            ylim=(-18,0), xtick=0:2:20, ytick=-20:2:0, framestyle=:box)
plot!(ps11, FGHz, s11db, color=:blue, label="PSSFSS", lw=2)
plot!(ps11, li[:,1], li[:,2], color=:red, label="Li et al.", lw=2)
#md savefig(ps11, "flexibleabsorbers11.png"); nothing  # hide
#-
#md # ![](flexibleabsorbers11.png)
#-
# This PSSFSS run takes about 44 seconds on my machine for 196 frequencies covering a 40:1 bandwidth.

#-
# ### Conclusion
# PSSFSS results agree well with those of the paper.
