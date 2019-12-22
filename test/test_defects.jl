using Graphene
using Test

@test isequal([Tuple(test_xy[:,1])], index2xy(1, test_indexed_atoms))
@test isequal([Tuple(test_xy[:,12])], index2xy(12, test_indexed_atoms))
@test isequal(map(x -> Tuple(x), eachcol(test_xy[:,1:2:9])), index2xy(1:2:9, test_indexed_atoms))

test_gatom = GAtom(1, 2., 3., [4, 5, 6, 7], "GATOM_Signature", 8, "dataset_9")
test_gbond = GBond(11, 22., 33., [44, 55, 66, 77], "GBOND_Signature", 88, "dataset_99")
test_gpolygon = GPolygon(111, 222., 333., [444, 555, 666, 777], "GPOLYGON_Signature", 888, "dataset_999", 6)
test_gdefect = GDefect(1111, 2222., 3333., [4444, 5555, 6666, 7777], "GDEFECT_Signature", 8888, "dataset_9999", 14, "Divacancy")

@testset "GAtom" begin
    @test typeof(test_gatom) == GAtom
    @test get_id(test_gatom) == 1
    @test get_x(test_gatom) == 2.
    @test get_y(test_gatom) == 3.
    @test get_relatives(test_gatom) == [4, 5, 6, 7]
    @test get_signature(test_gatom) == "GATOM_Signature"
    @test get_frame(test_gatom) == 8
    @test get_dataset(test_gatom) == "dataset_9"
    @test get_noa(test_gatom) == 1
    @test get_type(test_gatom) == "Atom"
end

@testset "GBond" begin
    @test typeof(test_gbond) == GBond
    @test get_id(test_gbond) == 11
    @test get_x(test_gbond) == 22.
    @test get_y(test_gbond) == 33.
    @test get_relatives(test_gbond) == [44, 55, 66, 77]
    @test get_signature(test_gbond) == "GBOND_Signature"
    @test get_frame(test_gbond) == 88
    @test get_dataset(test_gbond) == "dataset_99"
    @test get_noa(test_gbond) == 2
    @test get_type(test_gbond) == "Bond"
end

@testset "GPolygon" begin
    @test typeof(test_gpolygon) == GPolygon
    @test get_id(test_gpolygon) == 111
    @test get_x(test_gpolygon) == 222.
    @test get_y(test_gpolygon) == 333.
    @test get_relatives(test_gpolygon) == [444, 555, 666, 777]
    @test get_signature(test_gpolygon) == "GPOLYGON_Signature"
    @test get_frame(test_gpolygon) == 888
    @test get_dataset(test_gpolygon) == "dataset_999"
    @test get_noa(test_gpolygon) == 6
    @test get_type(test_gpolygon) == "Polygon"
end

@testset "GDefect" begin
    @test typeof(test_gdefect) == GDefect
    @test get_id(test_gdefect) == 1111
    @test get_x(test_gdefect) == 2222.
    @test get_y(test_gdefect) == 3333.
    @test get_relatives(test_gdefect) == [4444, 5555, 6666, 7777]
    @test get_signature(test_gdefect) == "GDEFECT_Signature"
    @test get_frame(test_gdefect) == 8888
    @test get_dataset(test_gdefect) == "dataset_9999"
    @test get_noa(test_gdefect) == 14
    @test get_type(test_gdefect) == "Divacancy"
end
