using Downloads
using GadgetGalaxies
using GadgetIO
using LinearAlgebra
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
        @test haskey(p, :id)
        @test haskey(p, "ID")
        @test !haskey(p, :pos)
        @test !haskey(p, "POS")
        @test values(p) === values(p.properties)
        @test issetequal(propertynames(p), [:type, :id, :mass, :zs, :im])

        p[:zs] = [6 5; 4 3; 2 1]
        @test p.Zs == [6 5; 4 3; 2 1]
        @test p.zs === p[:ZS]
        p["Zs"] = [1 2; 3 4; 5 6]
        @test p.Zs == [1 2; 3 4; 5 6]
        @test p.zs === p[:ZS]

        @test_throws ErrorException p[:type] = :dm

        @test particle_type_id(:gas) == 0

        q = copy(p)
        @test q !== p
        @test q.id === p.id
        q.id = p.id .+ 1
        @test q.id !== p.id
        @test q.id[1] != p.id[1]


        # galaxy
        snapshot = Snapshot("box", 10)
        dm = GG.Particles(:dm, Dict())
        g = Galaxy(snapshot, 1234, nothing, Dict(:stars => p))
        g.dm = dm

        @test g.snapshot === snapshot
        @test GG.getind(g) === g.isub === 1234
        @test isnothing(g.subid) && isnothing(GG.getid(g))
        @test g.stars === g[:stars] === p
        @test g.dm === g[:dm] === dm
        @test g.particles[:stars] === g.stars
        @test g.particles[:dm] === g.dm
        @test issetequal(keys(g), [:stars, :dm])
        @test haskey(g, :stars)
        @test !haskey(g, :gas)
        @test values(g) === values(g.particles)
        @test issetequal(propertynames(g), [:snapshot, :isub, :subid, :stars, :dm])

        g[:bh] = dm
        @test g.bh === dm

        h = copy(g)
        @test h !== g
        @test h.stars !== g.stars
        @test h.stars.id === g.stars.id
        h.stars.id = g.stars.id .+ 1
        @test h.stars.id !== g.stars.id

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


        # galaxy group
        snapshot = Snapshot("box", 10)

        gr = GalaxyGroup(snapshot, 1234, nothing, Dict(:stars => p))
        gr.dm = dm

        @test gr.snapshot === snapshot
        @test GG.getind(gr) === gr.igroup === 1234
        @test isnothing(gr.groupid) && isnothing(GG.getid(gr))
        @test gr.stars === gr[:stars] === p
        @test gr.dm === gr[:dm] === dm
        @test gr.particles[:stars] === gr.stars
        @test gr.particles[:dm] === gr.dm
        @test issetequal(keys(gr), [:stars, :dm])
        @test values(gr) === values(gr.particles)
        @test issetequal(propertynames(gr), [:snapshot, :igroup, :groupid, :stars, :dm])

        gr[:bh] = dm
        @test gr.bh === dm

        gr = GalaxyGroup(snapshot, 1234, false)
        @test gr.snapshot === snapshot
        @test gr.igroup === 1234
        @test isnothing(gr.groupid)
        @test isempty(keys(gr))
        @test gr.igroup === 1234

        gr.stars = p
        @test gr[:stars] === p

        @test_throws ErrorException gr[:igroup] = 1235

        io = IOBuffer()
        show(io, "text/plain", gr)
        @test String(take!(io)) == "GalaxyGroup igroup 1234\nstars: 2 Particles\n ID MASS Zs iM\n"

        snapshot = Snapshot("snap_002", "sub_002")
        gr = GalaxyGroup(snapshot, 0)
        @test gr.groupid == HaloID(0, 1)
    end


    # create real snapshot for rest of tests
    snapshot = Snapshot("snap_002", "sub_002")
    g = Galaxy(snapshot, 0)
    gr = GalaxyGroup(snapshot, 0)
    h = read_header(snapshot.subbase)

    @testset "Units" begin
        v, w = 1.0, 1.0f0
        vs = ones(Float64, 3)
        ws = ones(Float32, 3)

        a = 0.277777
        b = a * u"kpc"
        @test convert_units_pos(v, h, :sim) == v
        @test convert_units_pos(v, h, :physical) ≈ a rtol = 1e-5
        @test convert_units_physical_pos(v, h) ≈ a rtol = 1e-5
        @test convert_units_physical_pos(v, h) |> typeof === Float64
        @test convert_units_physical_pos(w, h) ≈ a rtol = 1e-5
        @test convert_units_physical_pos(w, h) |> typeof === Float32
        @test convert_units_physical_pos(vs, h)[1] ≈ a rtol = 1e-5
        @test convert_units_physical_pos(vs, h) |> eltype === Float64
        @test convert_units_physical_pos(ws, h)[1] .≈ a rtol = 1e-5
        @test convert_units_physical_pos(ws, h) |> eltype === Float32
        @test convert_units_pos(v, h, :full) ≈ b rtol = 1e-5
        @test convert_units_full_pos(v, h) ≈ b rtol = 1e-5
        @test convert_units_full_pos(v, h) |> typeof <: Quantity{Float64}
        @test convert_units_full_pos(w, h) ≈ b rtol = 1e-5
        @test convert_units_full_pos(w, h) |> typeof <: Quantity{Float32}
        @test convert_units_full_pos(vs, h)[1] ≈ b rtol = 1e-5
        @test convert_units_full_pos(vs, h) |> eltype <: Quantity{Float64}
        @test convert_units_full_pos(ws, h)[1] .≈ b rtol = 1e-5
        @test convert_units_full_pos(ws, h) |> eltype <: Quantity{Float32}

        af = convert(Float32, a)
        bf = af * u"kpc"
        as = [a, a, a]
        asf = [af, af, af]
        bs = [b, b, b]
        bsf = [bf, bf, bf]
        @test simulation_units_pos(a, h) |> typeof === Float64
        @test simulation_units_pos(af, h) |> typeof === Float32
        @test simulation_units_pos(b, h) |> typeof === Float64
        @test simulation_units_pos(bf, h) |> typeof === Float32
        @test simulation_units_pos(a, h) ≈ v rtol = 1e-5
        @test simulation_units_pos(b, h) ≈ v rtol = 1e-5
        @test simulation_units_pos(b |> u"m", h) ≈ v rtol = 1e-5
        @test simulation_units_pos(as, h) |> eltype === Float64
        @test simulation_units_pos(asf, h) |> eltype === Float32
        @test simulation_units_pos(bs, h) |> eltype === Float64
        @test simulation_units_pos(bsf, h) |> eltype === Float32
        @test simulation_units_pos(as, h)[1] ≈ v rtol = 1e-5
        @test simulation_units_pos(asf, h)[1] ≈ v rtol = 1e-5
        @test simulation_units_pos(bs, h)[1] ≈ v rtol = 1e-5
        @test simulation_units_pos(bsf, h)[1] ≈ v rtol = 1e-5
        @test simulation_units_pos(bs .|> u"m", h)[1] ≈ v rtol = 1e-5
        simulation_units_pos!(as, h)
        @test as[1] ≈ v rtol = 1e-5

        vs2, ws2 = copy(vs), copy(ws)
        convert_units_physical_pos!(vs2, h)
        @test vs2[1] ≈ a rtol = 1e-5
        convert_units_physical_pos!(ws2, h)
        @test ws2[1] ≈ a rtol = 1e-5

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
        @test convert_units_age(v, h, :physical) ≈ a rtol = 1e-5
        @test convert_units_physical_age(v, h) ≈ a rtol = 1e-5
        @test convert_units_physical_age(v, h) |> typeof === Float64
        @test convert_units_physical_age(w, h) ≈ a rtol = 1e-5
        @test convert_units_physical_age(w, h) |> typeof === Float64
        @test convert_units_physical_age(vs, h)[1] ≈ a rtol = 1e-5
        @test convert_units_physical_age(vs, h) |> eltype === Float64
        @test convert_units_physical_age(ws, h)[1] ≈ a rtol = 1e-5
        @test convert_units_physical_age(ws, h) |> eltype === Float64
        @test convert_units_age(v, h, :full) ≈ b rtol = 1e-5
        @test convert_units_full_age(v, h) ≈ b rtol = 1e-5
        @test convert_units_full_age(v, h) |> typeof <: Quantity{Float64}
        @test convert_units_full_age(w, h) ≈ b rtol = 1e-5
        @test convert_units_full_age(w, h) |> typeof <: Quantity{Float64}
        @test convert_units_full_age(vs, h)[1] ≈ b rtol = 1e-5
        @test convert_units_full_age(vs, h) |> eltype <: Quantity{Float64}
        @test convert_units_full_age(ws, h)[1] ≈ b rtol = 1e-5
        @test convert_units_full_age(ws, h) |> eltype <: Quantity{Float64}
        @test convert_units_physical(vs, :age, h) == convert_units_physical_age(vs, h)
        @test convert_units_physical!(copy(vs), :age, h) == convert_units_physical_age(vs, h)
        @test convert_units_full(vs, :age, h) == convert_units_full_age(vs, h)

        vs2, ws2 = copy(vs), copy(ws)
        convert_units_physical_age!(vs2, h)
        @test vs2[1] ≈ a rtol = 1e-5
        convert_units_physical_age!(ws2, h)
        @test ws2[1] ≈ a rtol = 1e-5

        zs, mass = ones(11, 2), [42, 34]
        a = 24.0731
        @test convert_units_solar_metallicity(zs, mass)[1] ≈ a rtol = 1e-5


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
        @test read_redshift(snapshot) ≈ 4 rtol = 1e-5
        @test read_header_particle_mass(snapshot, :dm, :physical) == 0
        @test read_header_particle_mass(snapshot, :dm) == 0u"Msun"
        @test read_header_particle_mass(snapshot, :gas) == 0u"Msun"
        @test read_galaxy_prop(g, "MSUB", :sim) ≈ 9.17776 rtol = 1e-5
        @test read_galaxy_prop(g, "MSUB", :physical) ≈ 1.27468f11 rtol = 1e-5
        @test read_galaxy_prop(g, "MSUB", :full) ≈ 1.27468f11u"Msun" rtol = 1e-5
        @test read_galaxy_prop(g, "SPOS") |> length == 3
        @test read_galaxy_prop(g, "SPOS")[1] ≈ -413.854u"kpc" rtol = 1e-5
        @test read_galaxy_prop(g, "SVEL")[1] ≈ 83.1962u"km/s" rtol = 1e-5
        @test read_galaxy_prop(g, "GRNR") == 0
        info = "Reading property GRNR of halo HaloID(0, 1)"
        warn = "The quantity GRNR is returned in simulation units"
        @test_logs (:info, info) (:warn, warn) read_galaxy_prop(g, "GRNR"; verbose=true)
        @test GG.convert_units_subfind_prop(0.1, "SAGE", h) ≈ 1.06078Unitful.Gyr rtol = 1e-5
        @test GG.convert_units_subfind_prop(0.1, "SSFR", h) == 0.1u"Msun" / Unitful.yr

        @test read_galaxy_pos(g) == read_galaxy_pos(gr)
        @test read_galaxy_vel(g) == read_galaxy_vel(gr)
        @test is_main_halo(g)
    end

    @testset "Read Halo" begin
        props = ((:dm, ["POS"]), (:gas, ["MASS", "POS", "VEL"]))
        read_halo!(g; use_keys=false, props, units=:full)

        @test typeof(g) <: Galaxy
        @test typeof(g.dm) <: GG.Particles
        @test typeof(g.gas) <: GG.Particles

        io = IOBuffer()
        show(io, "text/plain", g)
        @test String(take!(io)) ==
              "Galaxy at z=4.0\n isub 0\ndm: 7 Particles\n ID POS\ngas: 20 Particles\n ID MASS POS VEL\n"

        read_halo!(gr; use_keys=false, props, units=:full)

        @test typeof(gr) <: GalaxyGroup
        @test typeof(gr.dm) <: GG.Particles
        @test typeof(gr.gas) <: GG.Particles

        io = IOBuffer()
        show(io, "text/plain", gr)
        @test String(take!(io)) ==
              "GalaxyGroup at z=4.0\n igroup 0\ndm: 70 Particles\n ID POS\ngas: 76 Particles\n ID MASS POS VEL\n"
    end

    @testset "Transformation" begin
        a = [1 2 3 4; 5 6 7 8; 9 0 1 1]
        b = [9 0 1 1; 5 6 7 8; -1 -2 -3 -4]
        R = [0 0 1; 0 1 0; -1 0 0]
        @test rotate(a, R) == b
        rotate!(a, R)
        @test a == b

        p11 = 2.45564u"kpc"
        v12 = 90.39079u"km/s"
        @test g.dm.pos[1, 1] ≈ p11 rtol = 1e-5
        @test g.gas.vel[1, 2] ≈ v12 rtol = 1e-5
        gc = rotate(g, R)
        dmc = rotate(g.dm, R)
        @test gc.dm.pos[3, 1] ≈ -p11 rtol = 1e-5
        @test dmc.pos[3, 1] ≈ -p11 rtol = 1e-5
        @test gc.gas.vel[3, 2] ≈ -v12 rtol = 1e-5
        @test g.dm.pos != gc.dm.pos
        @test g.gas.vel != gc.gas.vel
        rotate!(g, R)
        @test g.dm.pos[3, 1] ≈ -p11 rtol = 1e-5
        @test g.gas.vel[3, 2] ≈ -v12 rtol = 1e-5
        rotate!(dmc, R)
        @test dmc.pos[3, 1] ≈ 3.05406u"kpc" rtol = 1e-5

        @test rotate_edgeon(gr.gas; algorithm=:unw_i_vol).pos != gr.gas.pos
        @test rotate_edgeon(gr, :gas; radius=20u"kpc").dm.pos != gr.dm.pos

        gr2 = deepcopy(gr)
        rotate_edgeon!(gr2.gas).pos
        @test gr2 != gr.gas.pos
        gr2 = deepcopy(gr)
        rotate_edgeon!(gr2, :gas).dm.pos
        @test gr2.dm.pos != gr.dm.pos

        @test rotate_edgeon(gr.gas; axis_ratios=true) isa Tuple
        @test rotate_edgeon!(gr2.gas; axis_ratios=true) isa Tuple
        @test rotate_edgeon(gr, :gas; axis_ratios=true) isa Tuple
        @test rotate_edgeon!(gr2, :gas; axis_ratios=true) isa Tuple
    end

    @testset "Shapes" begin
        let err = nothing
            try
                GG.get_algorithm_variables(:abc)
            catch err
            end
            @test err isa ErrorException
        end

        # 2d to 3d rotation matrix
        θ = π / 3
        R2d = [cos(θ) -sin(θ); sin(θ) cos(θ)]
        R3dx = GG.rotation_matrix_to_3d(R2d, [2, 3])
        R3dy = GG.rotation_matrix_to_3d(R2d, [3, 1])
        R3dz = GG.rotation_matrix_to_3d(R2d, [1, 2])
        @test R3dx == [1 0 0; 0 cos(θ) -sin(θ); 0 sin(θ) cos(θ)]
        @test R3dy == [cos(θ) 0 sin(θ); 0 1 0; -sin(θ) 0 cos(θ)]
        @test R3dz == [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]

        # parameters
        @test eccentricity(0.3) ≈ sqrt(1 - 0.3^2)
        @test ellipticity(0.3) ≈ 1 - 0.3
        @test triaxiality(0.7, 0.6) ≈ (1 - 0.7^2) / (1 - 0.6^2)

        # rotation matrix
        𝐈 = [2 -1 0; -1 2 -1; 0 -1 2]
        Q⁻¹, q, s = rotation_matrix_axis_ratios(𝐈)
        @test Q⁻¹[1, 1] ≈ Q⁻¹[1, 3] ≈ Q⁻¹[3, 1] ≈ Q⁻¹[3, 3] ≈ 0.5
        @test -Q⁻¹[1, 2] ≈ Q⁻¹[2, 1] ≈ -Q⁻¹[2, 3] ≈ Q⁻¹[3, 2] ≈ √0.5
        @test Q⁻¹[2, 2] ≈ 0 atol = 1e-14
        @test q ≈ 0.765366 rtol = 1e-5
        @test s ≈ 0.414213 rtol = 1e-5

        Q⁻¹, q, s, axes = rotation_matrix_axis_ratios(𝐈; return_axes=true)
        @test axes[1][1] ≈ axes[1][3] ≈ axes[3][1] ≈ axes[3][3] ≈ 0.5
        @test axes[2][1] ≈ -axes[2][3] ≈ -axes[1][2] ≈ axes[3][2] ≈ √0.5
        @test axes[2][2] ≈ 0 atol = 1e-14

        # inertia matrix: upper triangle
        𝐈 = zeros(Float32, 3, 3)
        𝐈_unw = [3103.08 -240.712 24.434; 0.0 19302.6 3955.77; 0.0 0.0 5671.47]
        𝐈_red = [16.1944 -2.31039 0.309565; 0.0 39.9156 3.32378; 0.0 0.0 19.89]
        GG.fill_𝐈_up!(𝐈, gr.gas.pos, gr.gas.mass, false, :unweighted)
        @test isapprox.(𝐈, 𝐈_unw) |> all
        GG.fill_𝐈_up!(𝐈, gr.gas.pos, gr.gas.mass, false, :reduced)
        @test isapprox.(𝐈, 𝐈_red) |> all
        GG.fill_𝐈_up!(𝐈, gr.gas.pos, gr.gas.mass, true, :unweighted)
        @test isapprox.(𝐈 ./ ustrip(gr.gas.mass[1]), 𝐈_unw) |> all # only works because of same gas mass of all
        GG.fill_𝐈_up!(𝐈, gr.gas.pos, gr.gas.mass, true, :reduced)
        @test isapprox.(𝐈 ./ ustrip(gr.gas.mass[1]), 𝐈_red) |> all # only works because of same gas mass of all

        # distances
        @test r²_ellipsoid(gr.gas.pos, 0.7, 0.5)[1] ≈ 573.1933u"kpc^2" rtol = 1e-5
        @test r²_ellipsoid(gr.gas.pos, 0.7, nothing)[1] ≈ 514.1716u"kpc^2" rtol = 1e-5
        @test r²_sphere(gr.gas.pos)[1] ≈ r²_ellipsoid(gr.gas.pos, 1, 1)[1] rtol = 1e-5
        @test r²_circle([3 2 3; 4 3 1])[1] == 25
        @test r²_sphere(gr.gas.pos[[1, 2], :])[1] ≈ r²_circle(gr.gas.pos[[1, 2], :])[1] rtol = 1e-5

        # ellipsoidal masking
        sph = Sphere(20u"kpc")
        @test GG.ellipsoidal_mask(gr.gas.pos, sph) |> count == 47
        ell = Ellipsoid(20u"kpc", 0.7, 0.5, :volume)
        @test GG.ellipsoidal_mask(gr.gas.pos, ell) |> count == 45
        @test GG.ellipsoidal_mask(gr.gas.pos, ell; return_r²_ellipsoid=true)[2][1] ≈ 573.1933u"kpc^2" rtol =
            1e-5
        ell = Ellipsoid(20u"kpc", 0.7, 0.5, :axis)
        @test GG.ellipsoidal_mask(gr.gas.pos, ell) |> count == 35

        # inertia matrices
        𝐈_unw = Symmetric(𝐈_unw, :U)
        𝐈_red = Symmetric(𝐈_red, :U)
        @test isapprox.(GG.inertia_matrix(gr.gas.pos, nothing; mass_weighted=false), 𝐈_unw) |> all
        @test isapprox.(GG.inertia_matrix(gr.gas.pos, gr.gas.mass) ./ ustrip(gr.gas.mass[1]), 𝐈_unw) |> all
        @test isapprox.(
            GG.inertia_matrix(gr.gas.pos, nothing; mass_weighted=false, inertia_matrix_type=:reduced),
            𝐈_red,
        ) |> all
        @test isapprox.(
            GG.inertia_matrix(gr.gas.pos, gr.gas.mass; inertia_matrix_type=:reduced) ./
            ustrip(gr.gas.mass[1]),
            𝐈_red,
        ) |> all

        @test isapprox.(
            GG.inertia_matrix_iterative(gr.gas.pos, nothing, :unweighted; mass_weighted=false),
            𝐈_unw,
        ) |> all
        @test isapprox.(
            GG.inertia_matrix_iterative(gr.gas.pos, gr.gas.mass, :unweighted) ./ ustrip(gr.gas.mass[1]),
            𝐈_unw,
        ) |> all
        @test isapprox.(
            GG.inertia_matrix_iterative(gr.gas.pos, nothing, :reduced; mass_weighted=false),
            𝐈_red,
        ) |> all
        @test isapprox.(
            GG.inertia_matrix_iterative(gr.gas.pos, gr.gas.mass, :reduced) ./ ustrip(gr.gas.mass[1]),
            𝐈_red,
        ) |> all

        # iterative inertia matrices
        radius = 20u"kpc"
        𝐈_unwi = [1749.14 -1014.96 453.757; -1014.96 2522.14 -459.417; 453.757 -459.417 1034.89]
        𝐈_redi = [13.8342 -2.85969 1.23256; -2.85969 18.0241 -1.74392; 1.23256 -1.74392 14.1417]
        𝐈_redelli = [7.99452 -4.20033 1.454; -4.20033 12.72 -1.5627; 1.454 -1.5627 6.82152]
        𝐈_unwi_ax = [285.81 -420.516 164.467; -420.516 667.925 -272.656; 164.467 -272.656 145.668]
        𝐈_redi_ax = [13.1609 -2.25873 1.22684; -2.25873 16.9281 -2.09812; 1.22684 -2.09812 13.911]
        @test isapprox.(
            GG.inertia_matrix_iterative(gr.gas.pos, nothing, :unweighted; radius, mass_weighted=false),
            𝐈_unwi,
        ) |> all
        @test isapprox.(
            GG.inertia_matrix_iterative(gr.gas.pos, gr.gas.mass, :unweighted; radius) ./
            ustrip(gr.gas.mass[1]),
            𝐈_unwi,
        ) |> all
        @test isapprox.(
            GG.inertia_matrix_iterative(gr.gas.pos, nothing, :reduced; radius, mass_weighted=false),
            𝐈_redi,
        ) |> all
        @test isapprox.(
            GG.inertia_matrix_iterative(
                gr.gas.pos,
                nothing,
                :reduced;
                radius,
                mass_weighted=false,
                ellipsoidal_distance=true,
            ),
            𝐈_redelli,
        ) |> all
        @test isapprox.(
            GG.inertia_matrix_iterative(
                gr.gas.pos,
                nothing,
                :unweighted;
                radius,
                mass_weighted=false,
                constant_volume=false,
            ),
            𝐈_unwi_ax,
        ) |> all
        @test isapprox.(
            GG.inertia_matrix_iterative(
                gr.gas.pos,
                nothing,
                :reduced;
                radius,
                mass_weighted=false,
                constant_volume=false,
            ),
            𝐈_redi_ax,
        ) |> all

        # full rotation matrix
        R_unw = [
            0.013092664 -0.96556526 -0.25983113
            0.05690142 -0.25871283 0.96427697
            0.99829394 0.027409703 -0.05155479
        ]
        R_unw_r = [
            0.5127014 -0.8163861 0.2658031
            0.7872874 0.57054496 0.23378807
            0.34251392 -0.08939998 -0.9352497
        ]
        R_unw_r_i = [
            0.5581614 -0.78767705 0.2608086
            0.67917794 0.61428696 0.40170747
            0.4766269 0.047082156 -0.8778439
        ]
        R, q, s = rotation_matrix_edgeon(gr.gas; algorithm=:unw)
        @test isapprox.(R, R_unw; rtol=1e-5) |> all
        R, q, s = rotation_matrix_edgeon(gr.gas; algorithm=:unw, radius)
        @test isapprox.(R, R_unw_r; rtol=1e-5) |> all
        R, q, s = rotation_matrix_edgeon(gr.gas; algorithm=:unw_i_vol, radius)
        @test isapprox.(R, R_unw_r_i; rtol=1e-5) |> all

        # 2D rotation matrix
        R_unw = [0.009512025 0.0 0.9999547; 0.0 1.0 0.0; 0.9999547 0.0 -0.009512025]
        R_unw_r = [1.0 0.0 0.0; 0.0 0.9836683 -0.17999081; 0.0 0.17999081 0.98366827]
        R_unw_r_i = [0.43801472 -0.89896786 0.0; 0.898968 0.43801472 0.0; 0.0 0.0 1.0]
        R, q = rotation_matrix_2D(gr.gas; algorithm=:unw, perspective=:edgeon)
        @test isapprox.(R, R_unw; rtol=1e-5) |> all
        R, q = rotation_matrix_2D(gr.gas; algorithm=:unw, radius, perspective=:sideon)
        @test isapprox.(R, R_unw_r; rtol=1e-5) |> all
        R, q = rotation_matrix_2D(gr.gas; algorithm=:unw_i_vol, radius, perspective=:faceon)
        @test isapprox.(R, R_unw_r_i; rtol=1e-5) |> all
        @test rotation_matrix_2D(gr.gas; algorithm=:unw, matrix_2d=true)[1] |> size == (2, 2)

        @test_nowarn rotation_matrix_edgeon(gr.gas; algorithm=:red, radius)
        @test_nowarn rotation_matrix_edgeon(gr.gas; algorithm=:unw_i_ax, radius)
        @test_nowarn rotation_matrix_edgeon(gr.gas; algorithm=:red_i_vol, radius)
        @test_nowarn rotation_matrix_edgeon(gr.gas; algorithm=:red_i_ax, radius)
        @test_nowarn rotation_matrix_edgeon(gr.gas; algorithm=:red_ell_i_vol, radius)
        @test_nowarn rotation_matrix_edgeon(gr.gas; algorithm=:red_ell_i_ax, radius)

        let err = nothing
            try
                GG.rotation_matrix_2D(gr.gas; algorithm=:unw, matrix_2d=false, perspective=:abc)
            catch err
            end
            @test err isa ErrorException
        end

        # translation
        Δx = [1, 2, 3] .* u"kpc"
        gr2 = deepcopy(gr)
        x = gr2.gas.pos[:, 1]
        @test translate(gr2.gas, Δx).pos[:, 1] == x .+ Δx
        translate!(gr2.gas, Δx)
        @test gr2.gas.pos[:, 1] == x .+ Δx

        x = gr2.gas.pos[:, 1]
        @test translate(gr2, Δx).gas.pos[:, 1] == x .+ Δx
        translate!(gr2, Δx)
        @test gr2.gas.pos[:, 1] == x .+ Δx

        # center of mass
        gr2 = deepcopy(gr)
        a = [-0.036875404, -0.22094725, 0.5704074] .* u"kpc"
        rstart = 40u"kpc"
        @test isapprox.(center_of_mass_iterative(gr2, rstart, :gas; limit_fraction=0.01), a; rtol=1e-5) |> all
        @test isapprox.(
            translate_to_center_of_mass_iterative(gr2, rstart, :gas).gas.pos,
            gr2.gas.pos .- a;
            rtol=1e-5,
        ) |> all
        @test gr2.gas.pos == gr.gas.pos
        translate_to_center_of_mass_iterative!(gr2, rstart, :gas)
        @test isapprox.(gr2.gas.pos, gr.gas.pos .- a; rtol=1e-5) |> all
        gr2 = deepcopy(gr)
        @test isapprox.(
            translate_to_center_of_mass_iterative(gr2.gas, rstart).pos,
            gr2.gas.pos .- a;
            rtol=1e-5,
        ) |> all
        translate_to_center_of_mass_iterative!(gr2.gas, rstart)
        @test isapprox.(gr2.gas.pos, gr.gas.pos .- a; rtol=1e-5) |> all

        # center of velocity
        gr2 = deepcopy(gr)
        a = [99.92409, -50.980724, -64.156525] .* u"km/s"
        radius = 40u"kpc"
        @test isapprox.(center_of_velocity(gr2, radius, :gas; p=0.9), a; rtol=1e-5) |> all
        @test isapprox.(
            translate_to_center_of_velocity(gr2, radius, :gas).gas.vel,
            gr2.gas.vel .- a;
            rtol=1e-5,
        ) |> all
        @test gr2.gas.vel == gr.gas.vel
        translate_to_center_of_velocity!(gr2, radius, :gas)
        @test isapprox.(gr2.gas.vel, gr.gas.vel .- a; rtol=1e-5) |> all
        gr2 = deepcopy(gr)
        @test isapprox.(translate_to_center_of_velocity(gr2.gas, radius).vel, gr2.gas.vel .- a; rtol=1e-5) |>
              all
        translate_to_center_of_velocity!(gr2.gas, radius)
        @test isapprox.(gr2.gas.vel, gr.gas.vel .- a; rtol=1e-5) |> all
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
