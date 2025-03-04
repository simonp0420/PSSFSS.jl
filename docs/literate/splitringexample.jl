#nb # %% A slide [markdown] {"slideshow": {"slide_type": "Slide"}}
# ## Split-Ring Resonator
# This example is taken from [fabian2021independently; Figure 3](@cite).
# It consists of three concentric split rings with gaps sequentially rotated by 180° situated
# on a thin dielectric slab.
# 
#

# Here is a script that analyzes this geometry:

# ```Julia
# using PSSFSS, Plots, DelimitedFiles
# 
# r123 = [3.7, 4.25, 4.8]
# d = 10.8
# w = 0.3
# g = 0.3
# a = r123 .- w/2
# b = a .+ w
# gapwidth = [g, g, g]
# gapcenter = [-90, 90, -90]
# s1 = [d, 0]
# s2 = [0, d]
# 
# sheet = splitring(;class='M', units=mm, sides=42, ntri=1500, a, b, s1, s2, gapwidth, gapcenter)
# display(plot(sheet, unitcell=true))
# strata = [
#     Layer()
#     sheet
#     Layer(ϵᵣ=10.2, tanδ=0.0023, width=0.13mm)
#     Layer()
#     ]
# 
# freqs = 1:0.02:14)
# steering = (θ = 0, ϕ = 0)
# results = analyze(strata, freqs, steering)
# s11dbvv = extract_result(results, @outputs s11db(v,v))
# p = plot(xlabel="Frequency (GHz)", ylabel="S₁₁ Amplitude (dB)",
#     xtick=1:14, xlim=(1, 14), ytick = -40:5:0, ylim=(-30,0),
#     legend=:bottom)
# plot!(p, freqs, s11dbvv, label="PSSFSS")
# dat = readdlm("../src/assets/fabian2021_fig3_digitized.csv", ',')
# plot!(p, dat[:,1], dat[:,2], label="Fabian (CST)")
# ```

# ![](./assets/fabian2021_element.png)

# ![](./assets/fabian2021_comparison.png)

# This run of 651 frequencies requires about 40 seconds on my machine.
# Generally good agreement is seen between the PSSFSS predicted reflection amplitude and that
# digitized from the paper. (The latter was obtained from a CST frequency domain analysis,
# according to the paper's authors.) However, there is a small discrepancy in the predicted resonant 
# frequencies that increases 
# with frequency, likely because both results are less well converged at higher frequencies.
# Also, the reflection amplitudes of the higher-frequency peaks are less than unity for the CST
# results, possibly because the authors may have included the finite conductivity of the metal traces.  
# This detail was not reported in the paper.

