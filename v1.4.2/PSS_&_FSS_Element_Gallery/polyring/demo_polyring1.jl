using PSSFSS, Plots
sheet = polyring(s1=[1,0], s2=[0,1], a=[0.2, 0.45], b=[0.35, 0.55], sides=4, orient = 45, units=cm, ntri=800)
p1 = plot(sheet, linecolor=:red, unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(4,3))
plot(p1, p2, layout = (1,2), size=(800,400))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

