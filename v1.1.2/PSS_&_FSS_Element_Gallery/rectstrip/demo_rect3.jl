using Plots, PSSFSS
strip = rectstrip(Nx=8, Ny=4, Px=1, Py=0.2, Lx=0.15, Ly=0.2, units=cm)
p1 = plot(strip, unitcell=true, linecolor=:red)
p2 = plot(strip, rep=(3,5), linecolor=:blue)
plot(p1, p2, layout=(1,2), size=(800,400))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

