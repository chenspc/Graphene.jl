using Graphene
using Test
using GeometricalPredicates: Point2D
using Statistics: mean

test_gatom = GAtom(1, 2., 3., Set([4, 5, 6, 7]), "GATOM_Signature", 8, "dataset_9")
test_gbond = GBond(11, 22., 33., Set([44, 55, 66, 77]), "GBOND_Signature", 88, "dataset_99")
test_gpolygon = GPolygon(111, 222., 333., Set([444, 555, 666, 777]), "GPOLYGON_Signature", 888, "dataset_999", 6)
test_gdefect = GDefect(1111, 2222., 3333., Set([4444, 5555, 6666, 7777]), "GDEFECT_Signature", 8888, "dataset_9999", 14, Set([1, 11, 111]), "Divacancy")

test_simple_gatom = GAtom(1, 2., 3., Set([]), "", 0, "dataset")
test_simple_gbond = GBond(11, 22., 33., Set([]), "", 0, "dataset")
test_simple_gpolygon = GPolygon(111, 222., 333., Set([]), "", 0, "dataset", 6)

@testset "GAtom" begin
    @test isa(test_gatom, GAtom)
    @test get_id(test_gatom) == 1
    @test get_x(test_gatom) == 2.
    @test get_y(test_gatom) == 3.
    @test get_relatives(test_gatom) == Set([4, 5, 6, 7])
    @test get_signature(test_gatom) == "GATOM_Signature"
    @test get_frame(test_gatom) == 8
    @test get_dataset(test_gatom) == "dataset_9"
    @test get_noa(test_gatom) == 1
    @test get_type(test_gatom) == "Atom"
    @test isatom(test_gatom)
    @test get_members(test_gatom) == Set([1])

    @test GAtom(1, 2., 3., 0, "dataset") == test_simple_gatom
    @test GAtom(1, 2., 3., 0) == test_simple_gatom
    @test GAtom(1, 2., 3., Set([])) == test_simple_gatom
    @test GAtom(1, 2., 3.) == test_simple_gatom
    @test GAtom(1, Point2D(2., 3.)) == test_simple_gatom
end

@testset "GBond" begin
    @test isa(test_gbond, GBond)
    @test get_id(test_gbond) == 11
    @test get_x(test_gbond) == 22.
    @test get_y(test_gbond) == 33.
    @test get_relatives(test_gbond) == Set([44, 55, 66, 77])
    @test get_signature(test_gbond) == "GBOND_Signature"
    @test get_frame(test_gbond) == 88
    @test get_dataset(test_gbond) == "dataset_99"
    @test get_noa(test_gbond) == 2
    @test get_type(test_gbond) == "Bond"
    @test isbond(test_gbond)
    @test get_members(test_gbond) == Set([11])

    @test GBond(11, 22., 33., 0, "dataset") == test_simple_gbond
    @test GBond(11, 22., 33., 0) == test_simple_gbond
    @test GBond(11, 22., 33., Set([])) == test_simple_gbond
    @test GBond(11, 22., 33.) == test_simple_gbond
end

@testset "GPolygon" begin
    @test isa(test_gpolygon, GPolygon)
    @test get_id(test_gpolygon) == 111
    @test get_x(test_gpolygon) == 222.
    @test get_y(test_gpolygon) == 333.
    @test get_relatives(test_gpolygon) == Set([444, 555, 666, 777])
    @test get_signature(test_gpolygon) == "GPOLYGON_Signature"
    @test get_frame(test_gpolygon) == 888
    @test get_dataset(test_gpolygon) == "dataset_999"
    @test get_noa(test_gpolygon) == 6
    @test get_type(test_gpolygon) == "Polygon"
    @test ispolygon(test_gpolygon)
    @test get_members(test_gpolygon) == Set([111])

    @test GPolygon(111, 222., 333., 0, "dataset", 6) == test_simple_gpolygon
    @test GPolygon(111, 222., 333., 0, 6) == test_simple_gpolygon
    @test GPolygon(111, 222., 333., 6) == test_simple_gpolygon
end

@testset "GDefect" begin
    @test isa(test_gdefect, GDefect)
    @test get_id(test_gdefect) == 1111
    @test get_x(test_gdefect) == 2222.
    @test get_y(test_gdefect) == 3333.
    @test get_relatives(test_gdefect) == Set([4444, 5555, 6666, 7777])
    @test get_signature(test_gdefect) == "GDEFECT_Signature"
    @test get_frame(test_gdefect) == 8888
    @test get_dataset(test_gdefect) == "dataset_9999"
    @test get_noa(test_gdefect) == 14
    @test get_members(test_gdefect) == Set([1, 11, 111])
    @test get_type(test_gdefect) == "Divacancy"
end

test_gatom1 = GAtom(1, 12., 13., Set([4, 5, 6, 7]), "GATOM_Signature", 8, "dataset_9")
test_gatom2 = GAtom(2, 22., 23., Set([4, 5, 6, 7]), "GATOM_Signature", 8, "dataset_9")
test_gbond1 = GBond(11, 122., 133., Set([44, 55, 66, 77]), "GBOND_Signature", 8, "dataset_9")
test_gbond2 = GBond(22, 222., 233., Set([44, 55, 66, 77]), "GBOND_Signature", 8, "dataset_9")
test_gpolygon1 = GPolygon(111, 1222., 1333., Set([444, 555, 666, 777]), "GPOLYGON_Signature", 8, "dataset_9", 5)
test_gpolygon2 = GPolygon(222, 2222., 2333., Set([444, 555, 666, 777]), "GPOLYGON_Signature", 8, "dataset_9", 7)
test_g_for_sort = [test_gatom1, test_gatom2, test_gpolygon2, test_gbond1, test_gpolygon1, test_gbond2]

@testset "isless and sort" begin
    @test isless(test_gatom1, test_gpolygon1)
    @test isless(test_gbond2, test_gpolygon2)
    @test isless(test_gpolygon1, test_gpolygon2)
end

link_relatives!(test_gatom1, test_gatom2)
@test get_relatives(test_gatom1) == Set([4, 5, 6, 7, 2])
@test get_relatives(test_gatom2) == Set([4, 5, 6, 7, 1])
link_relatives!(test_gbond1, test_gbond2)
@test get_relatives(test_gbond1) == Set([44, 55, 66, 77, 22])
@test get_relatives(test_gbond2) == Set([44, 55, 66, 77, 11])
link_relatives!(test_gpolygon1, test_gpolygon2)
@test get_relatives(test_gpolygon1) == Set([444, 555, 666, 777, 222])
@test get_relatives(test_gpolygon2) == Set([444, 555, 666, 777, 111])

link_relatives!(test_gatom1, test_gbond1)
@test get_relatives(test_gatom1) == Set([4, 5, 6, 7, 2, 11])
@test get_relatives(test_gbond1) == Set([44, 55, 66, 77, 22, 1])
link_relatives!(test_gbond2, test_gatom2)
@test get_relatives(test_gatom2) == Set([4, 5, 6, 7, 1, 22])
@test get_relatives(test_gbond2) == Set([44, 55, 66, 77, 11, 2])
link_relatives!(test_gatom1, test_gpolygon1)
@test get_relatives(test_gatom1) == Set([4, 5, 6, 7, 2, 11, 111])
@test get_relatives(test_gpolygon1) == Set([444, 555, 666, 777, 222, 1])
link_relatives!(test_gpolygon2, test_gatom2)
@test get_relatives(test_gatom2) == Set([4, 5, 6, 7, 1, 22, 222])
@test get_relatives(test_gpolygon2) == Set([444, 555, 666, 777, 111, 2])
link_relatives!(test_gbond1, test_gpolygon1)
@test get_relatives(test_gbond1) == Set([44, 55, 66, 77, 22, 1, 111])
@test get_relatives(test_gpolygon1) == Set([444, 555, 666, 777, 222, 1, 11])
link_relatives!(test_gpolygon2, test_gbond2)
@test get_relatives(test_gbond2) == Set([44, 55, 66, 77, 11, 2, 222])
@test get_relatives(test_gpolygon2) == Set([444, 555, 666, 777, 111, 2, 22])

test_gatom1 = GAtom(1, 12., 13., Set([4, 5, 6, 7]), "GATOM_Signature", 8, "dataset_9")
test_gbond1 = GBond(11, 122., 133., Set([44, 55, 66, 77]), "GBOND_Signature", 8, "dataset_9")
test_gpolygon1 = GPolygon(111, 1222., 1333., Set([444, 555, 666, 777]), "GPOLYGON_Signature", 8, "dataset_9", 6)
test_gpolygon2 = GPolygon(222, 2222., 2333., Set([444, 555, 666, 777]), "GPOLYGON_Signature", 8, "dataset_9", 6)

link_relatives!(test_gpolygon2, [test_gatom1, test_gbond1, test_gpolygon1])
@test get_relatives(test_gatom1) == Set([4, 5, 6, 7, 222])
@test get_relatives(test_gbond1) == Set([44, 55, 66, 77, 222])
@test get_relatives(test_gpolygon1) == Set([444, 555, 666, 777, 222])
@test get_relatives(test_gpolygon2) == Set([444, 555, 666, 777, 1, 11, 111])

test_gatom1 = GAtom(1, 12., 13., Set([4, 5, 6, 7]), "GATOM_Signature", 8, "dataset_9")
test_gbond1 = GBond(11, 122., 133., Set([44, 55, 66, 77]), "GBOND_Signature", 8, "dataset_9")
test_gbond2 = GBond(22, 222., 233., Set([44, 55, 66, 77]), "GBOND_Signature", 8, "dataset_9")
test_gpolygon1 = GPolygon(111, 1222., 1333., Set([444, 555, 666, 777]), "GPOLYGON_Signature", 8, "dataset_9", 6)

link_relatives!(test_gbond2, [test_gatom1, test_gbond1, test_gpolygon1])
@test get_relatives(test_gatom1) == Set([4, 5, 6, 7, 22])
@test get_relatives(test_gbond1) == Set([44, 55, 66, 77, 22])
@test get_relatives(test_gpolygon1) == Set([444, 555, 666, 777, 22])

test_gatom1 = GAtom(1, 12., 13., Set([4, 5, 6, 7]), "GATOM_Signature", 8, "dataset_9")
test_gpolygon1 = GPolygon(111, 1222., 1333., Set([444, 555, 666, 777]), "GPOLYGON_Signature", 8, "dataset_9", 6)
test_gpolygon2 = GPolygon(222, 2222., 2333., Set([444, 555, 666, 777]), "GPOLYGON_Signature", 8, "dataset_9", 6)

link_relatives!(test_gatom1, [test_gpolygon1, test_gpolygon2])
@test get_relatives(test_gatom1) == Set([4, 5, 6, 7, 111, 222])
@test get_relatives(test_gpolygon1) == Set([444, 555, 666, 777, 1])
@test get_relatives(test_gpolygon2) == Set([444, 555, 666, 777, 1])

@test make_gatom(collect(test_indexed_atoms)[5]) == GAtom(5, 15., 26.1, Set([]), "", 0, "dataset")
gatoms_vector = make_gatom.(collect(test_indexed_atoms))
@test length(gatoms_vector) == length(test_indexed_atoms)
@test gatoms_vector[5] == GAtom(5, 15., 26.1, Set([]), "", 0, "dataset")

test_gatom1 = GAtom(1, 12., 13., Set([4, 5, 6, 7]), "GATOM_Signature", 8, "dataset_9")
test_gatom2 = GAtom(2, 22., 23., Set([4, 5, 6, 7]), "GATOM_Signature", 8, "dataset_9")
@test make_gbond([test_gatom1, test_gatom2], 12) == GBond(12, 17., 18., Set([1, 2]), "", 0, "dataset")
@test make_gbond!([test_gatom1, test_gatom2], 12) == GBond(12, 17., 18., Set([1, 2]), "", 0, "dataset")
@test test_gatom1 == GAtom(1, 12., 13., Set([4, 5, 6, 7, 12]), "GATOM_Signature", 8, "dataset_9")
@test test_gatom2 == GAtom(2, 22., 23., Set([4, 5, 6, 7, 12]), "GATOM_Signature", 8, "dataset_9")

test_bondmatrix = make_bondmatrix(test_bonds)
@test count(!iszero, test_bondmatrix) == 30
@test sort(unique([x for x in test_bondmatrix if !iszero(x)])) == collect(1:15)

@test make_bondmatrix([(1,2)]) == [0 1; 1 0]
@test make_bondmatrix([(1,2), (1,4), (3,4)]) == [0 1 0 2;
                                                 1 0 0 0;
                                                 0 0 0 3;
                                                 2 0 3 0]

test_gpolygon1 = make_gpolygon([test_gatom1, test_gatom2], [test_gbond1, test_gbond2], 111)
test_gpolygon2 = make_gpolygon!([test_gatom1, test_gatom2], [test_gbond1, test_gbond2], 111)
@test get_x(test_gpolygon1) == get_x(test_gpolygon2) == mean(get_x.([test_gatom1, test_gatom2]))
@test get_y(test_gpolygon1) == get_y(test_gpolygon2) == mean(get_y.([test_gatom1, test_gatom2]))
@test get_id(test_gpolygon1) == get_id(test_gpolygon2)
@test get_relatives(test_gpolygon1) == Set([])
@test get_relatives(test_gpolygon2) == Set([1, 11, 2, 22])

test_gatom1 = GAtom(1, 12., 13., Set([66, 11, 111, 222]), "", 8, "dataset")
test_gatom2 = GAtom(2, 22., 23., Set([11, 22, 111, 222]), "", 8, "dataset")
test_gatom3 = GAtom(3, 32., 33., Set([22, 33, 111, 222]), "", 8, "dataset")
test_gatom4 = GAtom(4, 42., 43., Set([33, 44, 111, 222]), "", 8, "dataset")
test_gatom5 = GAtom(5, 52., 53., Set([44, 55, 111, 222]), "", 8, "dataset")
test_gatom6 = GAtom(6, 62., 63., Set([55, 66, 111, 222]), "", 8, "dataset")
test_gbond1 = GBond(11, 122., 133., Set([1, 2, 111, 222]), "", 8, "dataset")
test_gbond2 = GBond(22, 222., 233., Set([2, 3, 111, 222]), "", 8, "dataset")
test_gbond3 = GBond(33, 322., 333., Set([3, 4, 111, 222]), "", 8, "dataset")
test_gbond4 = GBond(44, 422., 433., Set([4, 5, 111, 222]), "", 8, "dataset")
test_gbond5 = GBond(55, 522., 533., Set([5, 6, 111, 222]), "", 8, "dataset")
test_gbond6 = GBond(66, 622., 633., Set([6, 1, 111, 222]), "", 8, "dataset")
test_gpolygon1 = GPolygon(111, 1222., 1333., Set([1, 2, 3, 4, 5, 6, 11, 22, 33, 44, 55, 66, 111]), "", 8, "dataset", 6)
test_gpolygon2 = GPolygon(222, 2222., 2333., Set([1, 2, 3, 4, 5, 6, 11, 22, 33, 44, 55, 66, 222]), "", 8, "dataset", 6)

test_graphene = vcat([test_gatom1, test_gatom2, test_gatom3, test_gatom4, test_gatom5, test_gatom6],
                      [test_gbond1, test_gbond2, test_gbond3, test_gbond4, test_gbond5, test_gbond6],
                      [test_gpolygon1, test_gpolygon2])
test_noa_dict = Dict([Pair(x._id, get_noa(x)) for x in test_graphene])
map(x -> make_signature!(x, test_noa_dict), test_graphene)

@test test_gatom1 == GAtom(1, 12., 13., Set([66, 11, 111, 222]), "6-2|", 8, "dataset")
@test test_gatom2 == GAtom(2, 22., 23., Set([11, 22, 111, 222]), "6-2|", 8, "dataset")
@test test_gatom3 == GAtom(3, 32., 33., Set([22, 33, 111, 222]), "6-2|", 8, "dataset")
@test test_gatom4 == GAtom(4, 42., 43., Set([33, 44, 111, 222]), "6-2|", 8, "dataset")
@test test_gatom5 == GAtom(5, 52., 53., Set([44, 55, 111, 222]), "6-2|", 8, "dataset")
@test test_gatom6 == GAtom(6, 62., 63., Set([55, 66, 111, 222]), "6-2|", 8, "dataset")
@test test_gbond1 == GBond(11, 122., 133., Set([1, 2, 111, 222]), "6-2|", 8, "dataset")
@test test_gbond2 == GBond(22, 222., 233., Set([2, 3, 111, 222]), "6-2|", 8, "dataset")
@test test_gbond3 == GBond(33, 322., 333., Set([3, 4, 111, 222]), "6-2|", 8, "dataset")
@test test_gbond4 == GBond(44, 422., 433., Set([4, 5, 111, 222]), "6-2|", 8, "dataset")
@test test_gbond5 == GBond(55, 522., 533., Set([5, 6, 111, 222]), "6-2|", 8, "dataset")
@test test_gbond6 == GBond(66, 622., 633., Set([6, 1, 111, 222]), "6-2|", 8, "dataset")
@test test_gpolygon1 == GPolygon(111, 1222., 1333., Set([1, 2, 3, 4, 5, 6, 11, 22, 33, 44, 55, 66, 111]), "6-1|", 8, "dataset", 6)
@test test_gpolygon2 == GPolygon(222, 2222., 2333., Set([1, 2, 3, 4, 5, 6, 11, 22, 33, 44, 55, 66, 222]), "6-1|", 8, "dataset", 6)

@testset "Make Graphene" begin
    test_graphene1 = make_graphene(test_xy, max_bondlength=8, frame=12, dataset="test_dataset")
    test_graphene2 = make_graphene(test_xy, max_bondlength=12, frame=12, dataset="test_dataset")
    test_graphene3 = make_graphene(test_xy[:, 2], max_bondlength=8, frame=12, dataset="test_dataset")
    test_graphene4 = make_graphene(test_xy[:, 2:4], max_bondlength=12, frame=12, dataset="test_dataset")
    test_graphene5 = make_graphene(test_xy[:, [2,4]], max_bondlength=8, frame=12, dataset="test_dataset")
    test_graphene6 = make_graphene(test_xy[:, [2,4]], max_bondlength=12, frame=12, dataset="test_dataset")
    @test length(test_graphene1) == 15
    @test length(test_graphene2) == 33
    @test length(test_graphene3) == 1
    @test length(test_graphene4) == 5
    @test length(test_graphene5) == 2
    @test length(test_graphene6) == 4

    @test get_id.(test_graphene1) == collect(1:length(test_graphene1))
    @test get_id.(test_graphene2) == collect(1:length(test_graphene2))
end

test_defect_image = import_stack("test_files/test_defect_6_bw.h5")
test_graphene = make_graphene(centroid2xy(make_centroids(test_defect_image)), max_bondlength=10)

test_flower = first(find_defect(test_graphene, "Flower"))
test_divacancy = first(find_defect(test_graphene, "V2(585)"))
test_divacancy14 = first(find_defect(test_graphene, "V2(14)"))
test_butterfly = first(find_defect(test_graphene, "Butterfly"))
test_stonewales = first(find_defect(test_graphene, "Stone-Wales"))
test_5775 = first(find_defect(test_graphene, "5775"))
test_defects = find_defect(test_graphene, "Flower", "Divacancy", "Butterfly", "Stone-Wales", "5775")
test_defects[1]._id = length(test_graphene) + 1
test_defects[2]._id = length(test_graphene) + 1
test_defects[3]._id = length(test_graphene) + 1
test_defects[4]._id = length(test_graphene) + 1
test_defects[5]._id = length(test_graphene) + 1
test_defects[6]._id = length(test_graphene) + 1
@test test_defects[1] == test_flower
@test test_defects[2] == test_divacancy
@test test_defects[3] == test_divacancy14
@test test_defects[4] == test_butterfly
@test test_defects[5] == test_stonewales
@test test_defects[6] == test_5775

@testset "Stone-Wales" begin
    @test get_id(test_stonewales) == length(test_graphene) + 1
    @test get_frame(test_stonewales) == 1
    @test get_dataset(test_stonewales) == "dataset"
    @test get_noa(test_stonewales) == 16
    @test length(get_relatives(test_stonewales)) == 10
    @test length(get_members(test_stonewales)) == 39
    @test get_type(test_stonewales) == "Stone-Wales"
end

@testset "5775" begin
    @test get_id(test_5775) == length(test_graphene) + 1
    @test get_frame(test_5775) == 1
    @test get_dataset(test_5775) == "dataset"
    @test get_noa(test_5775) == 18
    @test length(get_relatives(test_5775)) == 12
    @test length(get_members(test_5775)) == 43
    @test get_type(test_5775) == "5775"
end

@testset "Divacancy" begin
    @test get_id(test_divacancy) == length(test_graphene) + 1
    @test get_frame(test_divacancy) == 1
    @test get_dataset(test_divacancy) == "dataset"
    @test get_noa(test_divacancy) == 14
    @test length(get_relatives(test_divacancy)) == 10
    @test length(get_members(test_divacancy)) == 33
    @test get_type(test_divacancy) == "Divacancy"
end

@testset "Divacancy-14" begin
    @test get_id(test_divacancy14) == length(test_graphene) + 1
    @test get_frame(test_divacancy14) == 1
    @test get_dataset(test_divacancy14) == "dataset"
    @test get_noa(test_divacancy14) == 14
    @test length(get_relatives(test_divacancy14)) == 10
    @test length(get_members(test_divacancy14)) == 29
    @test get_type(test_divacancy14) == "Divacancy-14"
end

@testset "Flower" begin
    @test get_id(test_flower) == length(test_graphene) + 1
    @test get_frame(test_flower) == 1
    @test get_dataset(test_flower) == "dataset"
    @test get_noa(test_flower) == 22
    @test length(get_relatives(test_flower)) == 12
    @test length(get_members(test_flower)) == 55
    @test get_type(test_flower) == "Flower"
end

@testset "Butterfly" begin
    @test get_id(test_butterfly) == length(test_graphene) + 1
    @test get_frame(test_butterfly) == 1
    @test get_dataset(test_butterfly) == "dataset"
    @test get_noa(test_butterfly) == 30
    @test length(get_relatives(test_butterfly)) == 14
    @test length(get_members(test_butterfly)) == 77
    @test get_type(test_butterfly) == "Butterfly"
end
