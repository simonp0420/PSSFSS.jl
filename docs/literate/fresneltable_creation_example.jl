#nb # %% A slide [markdown] {"slideshow": {"slide_type": "Slide"}}
# ## Fresnel Table Creation
# Here we show how to create an HFSS-compatible Fresnel table using the [`res2fresnel`](@ref)
# function included with PSSFSS. The geometry for this example is a double square loop
# in a 5 cm square unit cell. The code for analyzing this geometry and creating the 
# Fresnel table is shown below:
# ```Julia
# using PSSFSS
# dwidth = 3mm
# duroid = Layer(epsr=2.2, tandel=0.0009, width=dwidth)
# a = √2 * [1, 2.125]
# b = √2 * [1.5, 3.125]
# units = mm
# sides = 4
# orient = 45
# sheet = polyring(;a, b, units, sides, orient, s1=[8,0], s2=[0,8], ntri=2000)
# strata = [Layer(), duroid, sheet, duroid, Layer(width=-2*dwidth)]
# steering = (θ=0:5:45, ϕ=0)
# FGHz = 10:2:20
# results = analyze(strata, FGHz, steering; resultfile="double_square_loop.res")
# res2fresnel(results, "double_square_loop.rttbl")
# # Alternatively: res2fresnel("double_square_loop.res", "double_square_loop.rttbl")
# ```
# Note that the final layer specifies a width that is the negative of the sum of the remaining layers' widths.
# This has the effect of moving the output phase reference plane to coincide with the input phase reference plane.
# Please see the documentation for [`res2fresnel`](@ref) for discussion of the setup requirements and restrictions
# necessary for creating a valid Fresnel table.

