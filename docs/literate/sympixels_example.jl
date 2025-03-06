#nb # %% A slide [markdown] {"slideshow": {"slide_type": "Slide"}}
# ## Pixelated Band Pass Filter
# This example of from [hong2021design; Fig. 4a](@cite).  It consists of a single-layer
# band-pass filter consisting of many square metallic "pixels", optimized using a special binary
# optimizer described in the cited paper. 
#
# The geometry for the pixelated filter is created using the [`sympixels`](@ref) function as 
# show below.
using PSSFSS, Plots
P = 5.4
units = mm
nrim = 1
halfnint = 26
patternvec = Bool[
    0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1,
    0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0,
    1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1,
    0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1,
    1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1,
    0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0,
    0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0,
    1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1,
    1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1,
    1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0,
    1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1,
    1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1,
    0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1,
    1, 0, 1, 1, 1, 1, 0, 1, 1, 1,
    0, 1, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 1, 1, 0, 0, 1, 1,
    0, 0, 0, 0, 1, 1, 0,
    1, 0, 0, 1, 0, 0,
    0, 0, 1, 1, 1,
    1, 0, 0, 0,
    0, 0, 1,
    0, 0,
    0]


sheet1 = sympixels(; P, units, nrim, halfnint, patternvec, pdiv=1, class='M')
plot(sheet1, size=(450,450), unitcell=true)
#-
# As described in the documentation for [`sympixels`](@ref), the `1` entries in `patternvec` above
# denote the locations of the metallized pixels.  However, by setting `class = 'M'` above we choose to 
# triangulate the empty pixel locations, i.e. those that are not occupied by metal.  As discussed in
# [Checkerboard Metasurface](@ref), this choice is necessary to avoid spurious results when analyzing 
# a structure of this type.  
#
# Note that the above triangulation is created by forming a square for each triangulated pixel and then
# adding a single diagonal.  A finer triangulation is required for an accurate analysis result.  This can 
# be accomplished by increasing `pdiv` to a value greater than `1`.  Here is the triangulation that results
# from `pdiv = 2`:
sheet2 = sympixels(; P, units, nrim, halfnint, patternvec, pdiv=2, class='M')
plot(sheet2, lw=0.5, size=(450,450), unitcell=true)
#-
# Now each pixel has been divided into an array of 2×2 squares, each of which receives a diagonal to 
# form triangles. The code for analyzing the `pdiv = 1` case would look like this:
# ```Julia
# strata = [Layer(), sheet1, Layer(epsr=3.28, tandel=0.007, width=0.05mm), Layer()]
# steering = (θ=0, ϕ=0)
# flist = range(start=20, stop=36, length=401)
# logfile = "bp_filter_pdiv1.log"
# resultfile = "bp_filter_pdiv1.res"
# results = analyze(strata, flist, steering; logfile, resultfile)
# s21db = extract_result(results, @outputs s21db(te,te))
# ```
#-
# Accurately analyzing a structure like this consisting of a huge number of pixels can be become
# very expensive (in terms of memory usage and CPU time) very quickly.  For example, the `pdiv = 1` 
# case generated 3056 basis functions and took about 30 seconds to analyze on a 32 GByte machine.
# But the `pdiv = 2` case generated 15,088 unknowns and required 683 seconds of execution time.  It 
# also failed on the 32 GByte machine due to lack of memory, but was run successfully on a 64 GByte
# machine.
#-
# Because the analysis in [hong2021design](@cite) used nonzero metal thickness, I ran the zero-thickness
# geometry in HFSS to observe the significance of the thickness.  Here is a plot comparing the digitized
# results from [hong2021design; Fig. 4a](@cite), with the zero-thickness HFSS analysis and the two
# PSSFSS analyses:
#-
# ![](./assets/sympixels_fig4a_comparison.svg)
#-
# As seen above, the finite metallization thickness does not have a large effect on the transmission trace.
# Using `pdiv = 2` produces much better agreement with HFSS, but it appears that an even larger value would
# be required to achieve full convergence.

