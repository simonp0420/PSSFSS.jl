using Plots, PSSFSS
strip = rectstrip(Nx=4, Ny=8, Px=0.2, Py=1, Lx=0.2, Ly=0.15, units=cm)
p1 = plot(strip, unitcell=true, linecolor=:red)
p2 = plot(strip, rep=(5,3), linecolor=:blue)
plot(p1, p2, layout=(1,2))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

