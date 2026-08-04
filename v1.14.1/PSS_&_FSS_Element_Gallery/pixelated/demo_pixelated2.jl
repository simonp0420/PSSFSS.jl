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
plot(p1, p2, layout = (1,2), size=(800,400))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
