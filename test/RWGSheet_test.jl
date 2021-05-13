using PSSFSS
import PSSFSS.Sheets: SV2, recttri, combine, read_sheet_data, write_sheet_data
using Test

bl = SV2([0.,0.])
tr = SV2([1.,1.])
nx = ny = 1
sh1 = recttri(bl, tr, nx, ny)   
@testset "recttri" begin
    @test length(sh1.ρ) == 4
    @test sh1.ρ[1] == [0.,0.]
    @test sh1.ρ[2] == [1.,0.]
    @test sh1.ρ[3] == [0.,1.]
    @test sh1.ρ[4] == [1.,1.]
    @test sh1.e1 == [1,3,1,2,1]
    @test sh1.e2 == [2,4,3,4,4]
    @test sh1.fe == reshape([2,3,5, 4,5,1], 3,2)
    @test sh1.fv == reshape([1,4,3, 1,2,4], 3,2)
end

bl = SV2([1.,0.])
tr = SV2([2.,1.])
nx = ny = 1
sh2 = recttri(bl, tr, nx, ny)
dup_coor = 'x'
dup_coor_value = 1
sh3 = combine(sh1, sh2, dup_coor, dup_coor_value)
@testset "combine right" begin
    @test sh3.ρ == [[0,0], [1,0], [0,1], [1,1], [2,0], [2,1]]
    @test sh3.e1 == [1, 3, 1, 2, 1, 2, 4, 5, 2]
    @test sh3.e2 == [2, 4, 3, 4, 4, 5, 6, 6, 6]
    @test sh3.fv == reshape([1,4,3,  1,2,4,  2,6,4,  2,5,6], 3,4)
    @test sh3.fe == reshape([2,3,5,  4,5,1,  7,4,9,  8,9,6], 3,4)
end


bl = SV2([0.,1.])
tr = SV2([1.,2.])
nx = ny = 1
sh2 = recttri(bl, tr, nx, ny)
dup_coor = 'y'
dup_coor_value = 1.0
sh3 = combine(sh1, sh2, dup_coor, dup_coor_value)
@testset "combine top" begin
    @test sh3.ρ == [[0,0], [1,0], [0,1], [1,1], [0,2], [1,2]]
    @test sh3.e1 == [1, 3, 1, 2, 1, 5, 3, 4, 3]
    @test sh3.e2 == [2, 4, 3, 4, 4, 6, 5, 6, 6]
    @test sh3.fe == [2 4 6 8; 3 5 7 9; 5 1 9 2]
    @test sh3.fv == [1 1 3 3; 4 2 6 4; 3 4 5 6]
end

bl = SV2([1.,1.])
tr = SV2([2.,2.])
nx = ny = 1
sh2 = recttri(bl, tr, nx, ny)
dup_coor = 'y'
dup_coor_value = 1.0
sh3 = combine(sh1, sh2, dup_coor, dup_coor_value)
@testset "combine topright" begin
    @test sh3.ρ == [[0, 0], [1, 0], [0, 1], 
                     [1, 1], [2, 1], [1, 2], [2, 2]]
    @test sh3.e1 == [1, 3, 1, 2, 1, 4, 6, 4, 5, 4]
    @test sh3.e2 == [2, 4, 3, 4, 4, 5, 7, 6, 7, 7]
    @test sh3.fe == [2 4 7 9; 3 5 8 10; 5 1 10 6]
    @test sh3.fv == [1 1 4 4; 4 2 7 5; 3 4 6 7]
end

bl = SV2([2.,0.])
tr = SV2([3.,1.])
nx = ny = 1
sh2 = recttri(bl, tr, nx, ny)
dup_coor = ' '
dup_coor_value = NaN
sh3 = combine(sh1, sh2, dup_coor, dup_coor_value)
@testset "combine disjoint" begin
    @test sh3.ρ == [[0,0], [1,0], [0,1], [1,1], [2,0], [3,0], [2,1], [3,1]]
    @test sh3.e1 == [1, 3, 1, 2, 1, 5, 7, 5, 6, 5]
    @test sh3.e2 == [2, 4, 3, 4, 4, 6, 8, 7, 8, 8]
    @test sh3.fe == [2 4 7 9; 3 5 8 10; 5 1 10 6]
    @test sh3.fv == [1 1 5 5; 4 2 8 6; 3 4 7 8]
end

@testset "read and write" begin
    fname = tempname()
    write_sheet_data(fname, sh3)
    sh4 = read_sheet_data(fname)
    @test sh3 == sh4
end