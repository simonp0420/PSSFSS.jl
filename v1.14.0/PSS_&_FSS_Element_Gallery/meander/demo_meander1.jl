using PSSFSS, Plots
sheet = meander(a=0.3535, b=0.707, h=0.28, w1=0.018, w2=0.018, ntri=600, units=inch,rot=45)
p1 = plot(sheet, linecolor=:green, unitcell=true)
p2 = plot(sheet, linecolor=:green, rep=(4,3))
plot(p1, p2, layout = (1,2), size=(800,400))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
