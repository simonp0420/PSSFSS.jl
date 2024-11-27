using Plots, PSSFSS
strip = rectstrip(Nx=10, Ny=2, Px=1, Py=1, Lx=0.5, Ly=0.1, orient=45, units=cm)
p1 = plot(strip, unitcell=true, linecolor=:red)
p2 = plot(strip, rep=(3,3), linecolor=:blue)
plot(p1, p2, layout=(1,2), size=(600,300))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
