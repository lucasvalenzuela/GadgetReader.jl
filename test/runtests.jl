using GadgetGalaxies
using Test

GG = GadgetGalaxies

@testset "GadgetGalaxies.jl" begin

    @testset "Snapshot" begin
        @test GG.get_snapbase("box", 32) == joinpath("box", "snapdir_032", "snap_032")
        @test GG.get_subbase("box", 2) == joinpath("box", "groups_002", "sub_002")

        snapshot = Snapshot("snapbase", "subbase")
        @test snapshot.snapbase == "snapbase"
        @test snapshot.subbase == "subbase"

        snapshot = Snapshot("snapbase", nothing)
        @test snapshot.snapbase == "snapbase"
        @test isnothing(snapshot.subbase)

        snapshot = Snapshot(nothing, "subbase")
        @test isnothing(snapshot.snapbase)
        @test snapshot.subbase == "subbase"

        @test Snapshot(snapbase="snapbase") == Snapshot("snapbase", nothing)
        @test Snapshot(subbase="subbase") == Snapshot(nothing, "subbase")

        snapbase = joinpath("box", "snapdir_010", "snap_010")
        subbase = joinpath("box", "snapdir_010", "snap_010")

        Snapshot("box", 10) == Snapshot(snapbase, subbase)
        Snapshot("box", 10; snapbase=false) == Snapshot(nothing, subbase)
        Snapshot("box", 10; subbase=false) == Snapshot(snapbase, nothing)
    end
end
