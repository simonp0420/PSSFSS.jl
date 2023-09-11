using PSSFSS, Plots
sheet = diagstrip(P=5.2, w=0.21, units=mm, Nl=60, Nw=4, orient=45)
p1 = plot(sheet, linecolor=:red, unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(4,3))
plot(p1, p2, layout = (1,2), size=(800,400))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
