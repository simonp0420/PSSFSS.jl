# Tests that translating geometry does not break the simulation

using PSSFSS
using Test
using Unitful
using PSSFSS.Sheets:translate!

ntri=100
period = 20
sheet_list = []
name_list = []
loaded_sheet = loadedcross(s1=[period,0.0], s2=[0.0,period], L1=0.6*period, L2=0.2*period, w=0.05*period, units=μm, ntri=ntri, class='M')
push!(sheet_list, loaded_sheet)
push!(name_list, "Loaded Cross")
jerusalem_sheet = jerusalemcross(P=period, L1=0.8*period, L2=0.1*period, A=0.4*period, B=0.1*period, w=0.2*period, units=μm, ntri=ntri)
push!(sheet_list, jerusalem_sheet)
push!(name_list, "Jerusalem Cross")
rectstrip_sheet = rectstrip(Lx=period*0.5, Ly=period*0.5, Nx=10, Ny=10, Px=period, Py=period, units=μm)
push!(sheet_list, rectstrip_sheet)
push!(name_list, "Rectstrip")
polyring_sheet = polyring(units=μm, s1=[period,0.0], s2=[0.0,period], a=[0.2*period], b=[0.4*period], sides=8, ntri=ntri)
push!(sheet_list, polyring_sheet)
push!(name_list, "Polyring")
splitring_sheet = splitring(s1=[period,0.0], s2=[0.0,period], units=μm, a=[0.2*period], b=[0.4*period], sides=16, ntri=ntri, gapcenter=0, gapwidth=period*0.1)
push!(sheet_list, splitring_sheet)
push!(name_list, "Split ring")
sinuous_sheet = sinuous(s1=[period,0.0], s2=[0.0,period], units=μm, L2=period*0.9, w=0.075*period, b=[period*0.3, period*0.45], rc=period*0.15, arms=4, sides=20, ntri=ntri, g=period*0.1)
push!(sheet_list, sinuous_sheet)
push!(name_list, "Sinuous")



function compare_results(sheet1, dx, dy, translate_sheet, debug)
    flist1 = 5u"THz"
    tol = 1e-8
    steering1 = (phi=45, theta=45)

    # Set values for dx and dy for phase shift
    sheet_ps = deepcopy(sheet1)
    sheet_ps.dx = dx
    sheet_ps.dy = dy

    # Translate sheet
    sheet_t = deepcopy(sheet1)
    if translate_sheet
        translate!(sheet_t, dx, dy)
    end

    strata_ps = [Layer(), sheet1, Layer(width=10μm), sheet_ps, Layer()]
    strata_t = [Layer(),sheet1, Layer(width=10μm),  sheet_t, Layer()]
    if debug
        println("Sheet type: ", typeof(sheet1), "  nodes: ", nodecount(sheet1), "  faces: ", facecount(sheet1))
        println("dx, dy: ", dx, ", ", dy)
        println("sheet_ps.dx, sheet_ps.dy: ", sheet_ps.dx, ", ", sheet_ps.dy)
    end

    results_ps = analyze(strata_ps, flist1, steering1, showprogress=false,resultfile=devnull, logfile=devnull)
    results_t = analyze(strata_t, flist1, steering1, showprogress=false,resultfile=devnull, logfile=devnull)
    below_tol_flag = true
    for m in 1:2, n in 1:2
        gsmdiff = results_ps[1].gsm[m, n] - results_t[1].gsm[m, n]
        maxabserror = maximum(abs, gsmdiff)
        if debug
            println("(m, n) = ($m, $n): Max abs. error = ", maxabserror)
        end
        if maxabserror > tol
            below_tol_flag = false
        end
    end
    return below_tol_flag
    
end

dx = 0.05*period
dy = 0.05*period
for (i, sheet) in enumerate(sheet_list)
    @testset "$(name_list[i])" begin
        # With both sheets translated, the results should match
        @test compare_results(sheet, dx, dy, true, false)
        # With only one sheet translated, the results should not match
        @test !compare_results(sheet, dx, dy, false, false)
    end
end
