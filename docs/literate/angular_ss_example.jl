#nb # %% A slide [markdown] {"slideshow": {"slide_type": "Slide"}}
# ## Angle Selective Surface
# This example is taken [chen2023design](@cite).
# A three-sheet FSS is designed that is transparent to normally incident plane waves, but strongly
# attenuates obliquely incident waves.  All three sheets are swastika-shaped, with the outer two
# sheets being identical.  The shapes are generated in PSSFSS using the manji element.
#-

# We begin by setting up and plotting the geometries of the outer and inner sheets:
# ```Julia
# using PSSFSS, Plots
# 
# # Dimensions in mm using the paper's notation:
# Pₚ=50; Lₚ=115; L1ₚ=24.7; L2ₚ=42; L3ₚ=21.5; L4ₚ=45; W1ₚ=7.4; W2ₚ=3
# 
# ntri = 2500
# common_keywords = (units=mm, class='M', w2=0, a=0, s1=[Pₚ, 0], s2=[0, Pₚ], ntri=ntri)
# 
# outer = manji(; L3=W1ₚ, w=W1ₚ, L1=L1ₚ, L2=L2ₚ/2-W1ₚ, common_keywords...)
# pouter = plot(outer, unitcell=true, title="Outer", linecolor=:red, linewidth=0.5)
# 
# inner = manji(; L3=W2ₚ, w=W2ₚ, L1=L3ₚ, L2=L4ₚ/2-W2ₚ, common_keywords...)
# pinner = plot(inner, unitcell=true, title="Inner", linecolor=:red, linewidth=0.5)
# 
# p = plot(pouter, pinner, size=(600,300))
# ```
#-
#md # ![](./assets/ass_geometry.svg)
#-

# Here is the rest of the script that performs the normal incidence analysis and compares it against
# the same analysis performed in HFSS:

# ```Julia
# strata = [Layer()
#           outer 
#           Layer(width=L_paper/2 * 1mm)
#           inner
#           Layer(width=L_paper/2 * 1mm)
#           outer
#           Layer()]
# 
# flist = range(2.5, 3.5, 300)
# steering = (θ=0, ϕ=0)
# result = analyze(strata, flist, steering; logfile="ass1_ntri$ntri.log", resultfile="ass1_ntri$ntri.res")
# 
# s21db = extract_result(result, @outputs s21db(te,te))
# 
# p = plot(xlim=(2.5,3.5),ylim=(-50,0),xticks=2.5:0.1:3.5, minorticks=2,framestyle=:box,
#          title="Normal Incidence Transmission", xlabel="Frequency (GHz)", ylabel="20log₁₀|S₂₁| (dB)")
# plot!(p, flist, s21db, label="PSSFSS", linewidth=2)
# 
# using DelimitedFiles
# hfssdat = readdlm("assets/hfss_ass_full.csv", ',', Float64, skipstart=1)
# freqhfss = hfssdat[:,1]
# s21hfss = hfssdat[:,4]
# plot!(p, freqhfss, s21hfss, label="HFSS", linewidth=2)
# display(p)
# ```

#-
#md # ![](./assets/ass1.png)
#-

# The above run required about 126 seconds on my machine.
#-

# Below is the script that analyzes the angle selective surface at 3 GHz while varying the incidence angle:

# ```Julia
# using PSSFSS, Plots
# 
# # Dimensions in mm using the paper's notation:
# Pₚ=50; Lₚ=115; L1ₚ=24.7; L2ₚ=42; L3ₚ=21.5; L4ₚ=45; W1ₚ=7.4; W2ₚ=3
# 
# common_keywords = (units=mm, class='M', w2=0, a=0, s1=[Pₚ, 0], s2=[0, Pₚ], ntri=2500)
# outer = manji(; L3=W1ₚ, w=W1ₚ, L1=L1ₚ, L2=L2ₚ/2-W1ₚ, common_keywords...)
# inner = manji(; L3=W2ₚ, w=W2ₚ, L1=L3ₚ, L2=L4ₚ/2-W2ₚ, common_keywords...)
# 
# strata = [Layer()
#           outer 
#           Layer(width=Lₚ/2 * 1mm)
#           inner
#           Layer(width=Lₚ/2 * 1mm)
#           outer
#           Layer()]
# 
# flist = 3.0
# steering = (θ=0:2:60, ϕ=0)
# result = analyze(strata, flist, steering; logfile="ass2.log", resultfile="ass2.res")
# 
# s21te, s21tm = eachcol(extract_result(result, @outputs s21db(te,te) s21db(tm,tm)))
# 
# (; ϕ, θ) = steering
# p = plot(xlim=(0,60),ylim=(-40,0),xticks=0:10:60, yticks=-40:10:0, minorticks=2, framestyle=:box,
#          xlabel="Incidence Angle θ (deg)", ylabel="20log₁₀|S₂₁| (dB)", 
#          title="ϕ = $(ϕ)° Transmission at $flist GHz")
# plot!(p, θ, s21te, label="PSSFSS TE", linewidth=2)
# plot!(p, θ, s21tm, label="PSSFSS TM", linewidth=2)
# display(p)
# ```

#-
#md # ![](./assets/ass2.png)
#-

# Note how the transmission amplitude rapidly rolls off beyond about 15°.  This run required 
# about 13.5 minutes to complete, much longer than for normal incidence. This is due to fact that
# the incremental phase shift ψ₁ is not held constant during the analysis, requiring the spatial
# integrals to be recomputed for each new incidence angle.
