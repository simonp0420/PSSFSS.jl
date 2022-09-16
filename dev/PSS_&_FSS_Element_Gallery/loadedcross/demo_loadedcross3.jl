using PSSFSS, Plots
sheet = loadedcross(s1=[1,0], s2=[0.5,sqrt(3)/2], L1=0.9, L2=0.2, w=0.35, units=cm, ntri=300, orient=45)
p1 = plot(sheet, linecolor=:red, unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(4,3))
plot(p1, p2, layout = (1,2), size=(800,400))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

