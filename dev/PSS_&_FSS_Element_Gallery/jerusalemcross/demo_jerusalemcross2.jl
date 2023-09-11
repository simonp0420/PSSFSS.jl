using PSSFSS, Plots
sheet = jerusalemcross(P=1, L1=0.9, L2=0.12, w=0.04, A = 0.4, B = 0.12, units=cm, ntri=600)
p1 = plot(sheet, linecolor=:red, unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(4,3))
plot(p1, p2, layout = (1,2), size=(800,400))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
