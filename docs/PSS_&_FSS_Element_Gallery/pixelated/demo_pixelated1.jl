# ---
# title: Symmetric Pixelated Element
# cover: "assets/demo_pixelated1.png"
# description: "Pixelated element created by the sympixels function"
# ---

# For a pixelated element created with the `sympixels` function, one should always specify
# `'M'` for the element class.  This is discussed in the "checkerboard" usage example. 

using PSSFSS, Plots
P = 5.4
units = mm
nrim = 1
halfnint = 26
pdiv = 1
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

sheet = sympixels(; P, units, nrim, halfnint, patternvec, pdiv, class='M')
plot(sheet, axis=false, xlabel="", ylabel="", xtick=[], ytick=[], linecolor=:orange, size=(400,400), rep=(3,3)) #src
savefig("assets/demo_pixelated1.png") #src
p1 = plot(sheet, linecolor = "red", unitcell=true)
p2 = plot(sheet, linecolor = "blue", lw=0.35, rep = (3,3))
plot(p1, p2, layout = (1,2), size=(800,400))

