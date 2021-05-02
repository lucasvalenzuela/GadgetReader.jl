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

    @testset "Galaxy" begin
        p = GG.Particles(:stars, Dict("ID"=>[1, 2], "MASS"=>[1e11, 2e11]))
        p.pos = [1 2; 3 4; 5 6]

        @test p.type === :stars
        @test p.id == [1, 2]
        @test p.mass == [1e11, 2e11]
        @test p.pos == [1 2; 3 4; 5 6]
        @test p.pos == p[:pos]
        @test p.properties["POS"] === p.pos
        @test p.properties["MASS"] === p.mass
        @test issetequal(GG.particleproperties(p), [:id, :mass, :pos])
        @test issetequal(propertynames(p), [:type, :id, :mass, :pos])

        p[:pos] = [6 5; 4 3; 2 1]
        @test p.pos == [6 5; 4 3; 2 1]

        @test_throws ErrorException p[:type] = :dm


        snapshot = Snapshot("box", 10)
        dm = GG.Particles(:dm, Dict())
        g = Galaxy(snapshot, 1234, nothing, Dict(:stars=>p))
        g.dm = dm

        @test g.snapshot === snapshot
        @test g.isub === 1234
        @test isnothing(g.subid)
        @test g.stars === g[:stars] === p
        @test g.dm === g[:dm] === dm
        @test g.particles[:stars] === g.stars
        @test g.particles[:dm] === g.dm
        @test issetequal(GG.particletypes(g), [:stars, :dm])
        @test issetequal(propertynames(g), [:snapshot, :isub, :subid, :stars, :dm])

        g[:bh] = dm
        @test g.bh === dm

        g = Galaxy(snapshot, 1234, false)
        @test g.snapshot === snapshot
        @test g.isub === 1234
        @test isnothing(g.subid)
        @test isempty(GG.particletypes(g))
        @test g.isub === 1234

        g.stars = p
        @test g[:stars] === p

        @test_throws ErrorException g[:isub] = 1235
    end
end
