# ---
# title: Pixelated Square Patch
# cover: "assets/demo_pixelated2.png"
# description: "Pixelated element created by the pixels function"
# ---

# Specifying `'J'` for the element class as in this example will generate a warning in the Julia REPL.
# In this case the warning can be ignored because there are no pixels that intersect with their neighbors
# at only a single point. See the discussion in the "checkerboard" usage example. 

using PSSFSS, Plots
units = cm
P = 1
units = cm
pdiv = 10
patternmat = Bool[1 0 0 1
                  0 0 0 0
                  0 0 0 0
                  1 0 0 1]

sheet = pixels(; P, pdiv, patternmat, units, class='J')
p1 = plot(sheet, linecolor = :red, unitcell = true)
p2 = plot(sheet, linecolor = :blue, rep=(3,3))
plot(sheet, axis=false, xlabel="", ylabel="", xtick=[], ytick=[], linecolor=:green, size=(400,400), rep=(3,3)) #src
savefig("assets/demo_pixelated2.png") #src
plot(p1, p2, layout = (1,2), size=(800,400))

