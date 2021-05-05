using Downloads
using GadgetGalaxies
using GadgetIO
using Test
using Unitful
using UnitfulAstro

D = Downloads
GG = GadgetGalaxies

@testset "GadgetGalaxies.jl" begin

    @info "downloading test data..."
    url = "http://www.usm.uni-muenchen.de/~lboess/GadgetIO/"
    for i in 0:3
        D.download(joinpath(url, "sub_002.$i"), "./sub_002.$i")
        D.download(joinpath(url, "snap_002.$i"), "./snap_002.$i")
        D.download(joinpath(url, "snap_002.$i.key"), "./snap_002.$i.key")
    end
    D.download(joinpath(url, "snap_002.key.index"), "./snap_002.key.index")
    @info "done!"

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
        # particles
        p = GG.Particles(:stars, Dict("ID" => [1, 2], "MASS" => [1e11, 2e11]))
        p.zs = [1 2; 3 4; 5 6]
        p.iM = [3, 4]

        @test p.type === :stars
        @test p.id == [1, 2]
        @test p.mass == [1e11, 2e11]
        @test p.zs == [1 2; 3 4; 5 6]
        @test p.iM == [3, 4]
        @test p.zs == p[:zs]
        @test p.properties["Zs"] === p.zs
        @test p.properties["MASS"] === p.mass
        @test p.properties["Zs"] === p["Zs"]
        @test p.properties["iM"] === p.iM
        @test issetequal(keys(p), [:id, :mass, :zs, :im])
        @test values(p) === values(p.properties)
        @test issetequal(propertynames(p), [:type, :id, :mass, :zs, :im])

        p[:zs] = [6 5; 4 3; 2 1]
        @test p.Zs == [6 5; 4 3; 2 1]
        @test p.zs === p[:ZS]

        @test_throws ErrorException p[:type] = :dm

        @test particle_type_id(:gas) == 0


        # galaxy
        snapshot = Snapshot("box", 10)
        dm = GG.Particles(:dm, Dict())
        g = Galaxy(snapshot, 1234, nothing, Dict(:stars => p))
        g.dm = dm

        @test g.snapshot === snapshot
        @test g.isub === 1234
        @test isnothing(g.subid)
        @test g.stars === g[:stars] === p
        @test g.dm === g[:dm] === dm
        @test g.particles[:stars] === g.stars
        @test g.particles[:dm] === g.dm
        @test issetequal(keys(g), [:stars, :dm])
        @test values(g) === values(g.particles)
        @test issetequal(propertynames(g), [:snapshot, :isub, :subid, :stars, :dm])

        g[:bh] = dm
        @test g.bh === dm

        g = Galaxy(snapshot, 1234, false)
        @test g.snapshot === snapshot
        @test g.isub === 1234
        @test isnothing(g.subid)
        @test isempty(keys(g))
        @test g.isub === 1234

        g.stars = p
        @test g[:stars] === p

        @test_throws ErrorException g[:isub] = 1235

        io = IOBuffer()
        show(io, "text/plain", g)
        @test String(take!(io)) == "Galaxy isub 1234\nstars: 2 Particles\n ID MASS Zs iM\n"

        snapshot = Snapshot("snap_002", "sub_002")
        g = Galaxy(snapshot, 0)
        @test g.subid == HaloID(0, 1)
    end


    # create real snapshot for rest of tests
    snapshot = Snapshot("snap_002", "sub_002")
    g = Galaxy(snapshot, 0)
    h = read_header(GadgetIO.select_file(snapshot.subbase, 0))

    @testset "Units" begin
        v, w = 1.0, 1.0f0
        vs = ones(Float64, 3)
        ws = ones(Float32, 3)

        a = 0.277777
        b = a * u"kpc"
        @test convert_units_pos(v, h, :sim) == v
        @test convert_units_pos(v, h, :physical) ≈ a rtol=1e-5
        @test convert_units_physical_pos(v, h) ≈ a rtol=1e-5
        @test convert_units_physical_pos(v, h) |> typeof === Float64
        @test convert_units_physical_pos(w, h) ≈ a rtol=1e-5
        @test convert_units_physical_pos(w, h) |> typeof === Float32
        @test convert_units_physical_pos(vs, h)[1] ≈ a rtol=1e-5
        @test convert_units_physical_pos(vs, h) |> eltype === Float64
        @test convert_units_physical_pos(ws, h)[1] .≈ a rtol=1e-5
        @test convert_units_physical_pos(ws, h) |> eltype === Float32
        @test convert_units_pos(v, h, :full) ≈ b rtol=1e-5
        @test convert_units_full_pos(v, h) ≈ b rtol=1e-5
        @test convert_units_full_pos(v, h) |> typeof <: Quantity{Float64}
        @test convert_units_full_pos(w, h) ≈ b rtol=1e-5
        @test convert_units_full_pos(w, h) |> typeof <: Quantity{Float32}
        @test convert_units_full_pos(vs, h)[1] ≈ b rtol=1e-5
        @test convert_units_full_pos(vs, h) |> eltype <: Quantity{Float64}
        @test convert_units_full_pos(ws, h)[1] .≈ b rtol=1e-5
        @test convert_units_full_pos(ws, h) |> eltype <: Quantity{Float32}

        vs2, ws2 = copy(vs), copy(ws)
        convert_units_physical_pos!(vs2, h)
        @test vs2[1] ≈ a rtol=1e-5
        convert_units_physical_pos!(ws2, h)
        @test ws2[1] ≈ a rtol=1e-5

        # TODO: perhaps also test vel, temp, mass

        @test convert_units_physical(vs, :pos, h) == convert_units_physical_pos(vs, h)
        @test convert_units_physical(vs, :vel, h) == convert_units_physical_vel(vs, h)
        @test convert_units_physical(vs, :temp, h) == convert_units_physical_temp(vs, h)
        @test convert_units_physical(vs, :zs, h) == convert_units_physical_mass(vs, h)

        v, w = 0.1, 0.1
        vs = ones(Float64, 3) / 10
        ws = ones(Float32, 3) / 10

        a = 1.06078
        b = a * Unitful.Gyr
        @test convert_units_age(v, h, :sim) == v
        @test convert_units_age(v, h, :physical) ≈ a rtol=1e-5
        @test convert_units_physical_age(v, h) ≈ a rtol=1e-5
        @test convert_units_physical_age(v, h) |> typeof === Float64
        @test convert_units_physical_age(w, h) ≈ a rtol=1e-5
        @test convert_units_physical_age(w, h) |> typeof === Float64
        @test convert_units_physical_age(vs, h)[1] ≈ a rtol=1e-5
        @test convert_units_physical_age(vs, h) |> eltype === Float64
        @test convert_units_physical_age(ws, h)[1] ≈ a rtol=1e-5
        @test convert_units_physical_age(ws, h) |> eltype === Float64
        @test convert_units_age(v, h, :full) ≈ b rtol=1e-5
        @test convert_units_full_age(v, h) ≈ b rtol=1e-5
        @test convert_units_full_age(v, h) |> typeof <: Quantity{Float64}
        @test convert_units_full_age(w, h) ≈ b rtol=1e-5
        @test convert_units_full_age(w, h) |> typeof <: Quantity{Float64}
        @test convert_units_full_age(vs, h)[1] ≈ b rtol=1e-5
        @test convert_units_full_age(vs, h) |> eltype <: Quantity{Float64}
        @test convert_units_full_age(ws, h)[1] ≈ b rtol=1e-5
        @test convert_units_full_age(ws, h) |> eltype <: Quantity{Float64}
        @test convert_units_physical(vs, :age, h) == convert_units_physical_age(vs, h)
        @test convert_units_physical!(copy(vs), :age, h) == convert_units_physical_age(vs, h)
        @test convert_units_full(vs, :age, h) == convert_units_full_age(vs, h)

        vs2, ws2 = copy(vs), copy(ws)
        convert_units_physical_age!(vs2, h)
        @test vs2[1] ≈ a rtol=1e-5
        convert_units_physical_age!(ws2, h)
        @test ws2[1] ≈ a rtol=1e-5

        zs, mass = ones(11, 2), [42, 34]
        a = 24.0731
        @test convert_units_solar_metallicity(zs, mass)[1] ≈ a rtol=1e-5


        p = GG.Particles(:stars, Dict())
        p.pos, p.vel, p.zs = rand(3, 2), rand(3, 2), rand(11, 2)
        p.iM, p.mass, p.age, p.temp = rand(2), rand(2), rand(2), rand(2)
        convert_units!(p, h, :physical)
        @test eltype(p.pos) === Float64
        @test eltype(p.age) === Float64
        @test eltype(p.zs) === Float64
        @test size(p.zs) == (2,)

        p.pos, p.vel, p.zs = rand(3, 2), rand(3, 2), rand(11, 2)
        p.iM, p.mass, p.age, p.temp = rand(2), rand(2), rand(2), rand(2)
        convert_units!(p, h, :full)
        @test eltype(p.pos) <: Unitful.Length
        @test eltype(p.vel) <: Unitful.Velocity
        @test eltype(p.age) <: Unitful.Time
        @test eltype(p.temp) <: Unitful.Temperature
        @test eltype(p.mass) <: Unitful.Mass
        @test eltype(p.zs) === Float64
        @test size(p.zs) == (2,)
    end

    @testset "Utils" begin
        boxsize, boxsize_half = 10, 5
        @test GG.shift_across_box_border(1, 2, boxsize, boxsize_half) == 1
        @test GG.shift_across_box_border(8, 2, boxsize, boxsize_half) == -2
        @test GG.shift_across_box_border(2, 8, boxsize, boxsize_half) == 12
    end

    @testset "Read" begin
        @test read_redshift(snapshot) ≈ 4 rtol=1e-5
    end

    @testset "Read Halo" begin
        props = ((:dm, ["POS"]),(:gas, ["MASS", "VEL"]))
        read_halo!(g; use_keys=false, props, units=:full)

        @test typeof(g) <: Galaxy
        @test typeof(g.dm) <: GG.Particles
        @test typeof(g.gas) <: GG.Particles

        io = IOBuffer()
        show(io, "text/plain", g)
        @test String(take!(io)) ==
               "Galaxy at z=4.0\n isub 0\ndm: 7 Particles\n ID POS\ngas: 20 Particles\n ID MASS VEL\n"
    end

    @info "deleting test data..."
    for i in 0:3
        rm("./sub_002.$i")
        rm("./snap_002.$i")
        rm("./snap_002.$i.key")
    end
    rm("./snap_002.key.index")
    @info "done!"
end
