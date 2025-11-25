using Test
using FrontCutTracking

@testset "FrontTracker" begin
    tracker = FrontTracker()

    @test num_points(tracker) == 0
    idx1 = add_point!(tracker, (0.0, 0, 0))
    idx2 = add_point!(tracker, (1, 0, 0))
    idx3 = add_point!(tracker, (0, 1, 0))
    idx4 = add_point!(tracker, (0, 0, 1))

    @test num_points(tracker) == 4
    @test tracker.points[idx1] == (0.0, 0.0, 0.0)

    line_idx = add_line!(tracker, (idx1, idx2))
    @test num_lines(tracker) == 1
    @test lines_touching(tracker, idx1) == [tracker.lines[line_idx]]

    surf_idx_a = add_surface!(tracker, (idx1, idx2, idx3))
    surf_idx_b = add_surface!(tracker, (idx1, idx3, idx4))
    @test num_surfaces(tracker) == 2
    @test surfaces_touching(tracker, idx3) == [tracker.surfaces[surf_idx_a], tracker.surfaces[surf_idx_b]]

    @test_throws ArgumentError add_line!(tracker, (0, 10))
    @test_throws ArgumentError lines_touching(tracker, 10)

    remove_line!(tracker, line_idx)
    @test num_lines(tracker) == 0

    remove_point!(tracker, idx2)
    @test num_points(tracker) == 3
    @test all(p -> p != (1, 0, 0), tracker.points)
    @test tracker.surfaces == [(idx1, idx3 - 1, idx4 - 1)]

    reset!(tracker)
    @test num_points(tracker) == 0
    @test num_lines(tracker) == 0
    @test num_surfaces(tracker) == 0
end

@testset "Initializers" begin
    pt = point_front(position=(1, 2, 3))
    @test num_points(pt) == 1
    @test num_lines(pt) == 0
    @test num_surfaces(pt) == 0
    @test pt.points[1] == (1.0, 2.0, 3.0)

    line = line_front(length=2.0, segments=2, direction=(0, 0, 1), metadata=Dict(:tag => :line))
    @test num_points(line) == 3
    @test num_lines(line) == 2
    @test line.points[1] == (0.0, 0.0, -1.0)
    @test line.points[end] == (0.0, 0.0, 1.0)
    @test line.metadata[:tag] == :line

    circle = circle_front(radius=1.5, segments=8, axis=:z)
    @test num_points(circle) == 8
    @test num_lines(circle) == 8
    radii = map(p -> sqrt(p[1]^2 + p[2]^2), circle.points)
    @test all(isapprox(r, 1.5; atol=1e-8) for r in radii)

    patch = planar_patch_front(size=(2.0, 1.0), divisions=(2, 1), axis=:y)
    @test num_points(patch) == (2 + 1) * (1 + 1)
    @test num_lines(patch) == (1 + 1) * 2 + (2 + 1) * 1
    @test num_surfaces(patch) == 2 * 2 * 1
    @test all(isapprox(p[2], 0.0; atol=1e-10) for p in patch.points)

    sphere = sphere_front(radius=2.0, rings=3, segments=6, metadata=Dict(:shape => :sphere))
    @test num_points(sphere) == 2 + 3 * 6
    @test num_lines(sphere) == 6 * (2 * 3 + 1)
    @test num_surfaces(sphere) == 2 * 6 * 3
    @test sphere.metadata[:shape] == :sphere
    dists = map(p -> sqrt(p[1]^2 + p[2]^2 + p[3]^2), sphere.points)
    @test all(isapprox(d, 2.0; atol=1e-8) for d in dists)
    @test_throws ArgumentError sphere_front(radius=1.0, rings=0)
    @test_throws ArgumentError sphere_front(radius=1.0, segments=2)

    star = star_front(radius=2.0, spikes=5, dent=0.5, axis=:x)
    @test num_points(star) == 2 * 5 + 1
    @test num_lines(star) == 2 * 5
    @test num_surfaces(star) == 2 * 5
    @test all(isapprox(p[1], 0.0; atol=1e-10) for p in star.points)
    radii = map(p -> sqrt(p[2]^2 + p[3]^2), star.points[1:end-1])
    @test maximum(radii) ≈ 2.0
    @test minimum(radii) ≈ 1.0

    star_wire = star_front(radius=1.0, spikes=6, dent=0.3, filled=false)
    @test num_surfaces(star_wire) == 0
    @test num_points(star_wire) == 2 * 6

    @test_throws ArgumentError star_front(radius=1.0, spikes=2)
    @test_throws ArgumentError star_front(radius=1.0, dent=0.0)

    rings = 3
    segments = 10
    star3d = star3d_front(radius=1.0, spike_length=0.4, petals=7,
                          rings=rings, segments=segments,
                          metadata=Dict(:shape => :star3d))
    expected_points = 2 + rings * segments
    expected_lines = rings * segments + (rings - 1) * segments + 2 * segments
    expected_surfaces = 2 * segments * rings
    @test num_points(star3d) == expected_points
    @test num_lines(star3d) == expected_lines
    @test num_surfaces(star3d) == expected_surfaces
    radii3d = map(p -> sqrt(p[1]^2 + p[2]^2 + p[3]^2), star3d.points)
    min_radius = minimum(radii3d)
    max_radius = maximum(radii3d)
    @test min_radius ≥ 1.0 * (1 - 0.4) - 1e-8
    @test max_radius ≤ 1.0 * (1 + 0.4) + 1e-8
    @test star3d.metadata[:shape] == :star3d
    @test_throws ArgumentError star3d_front(radius=0.0)
    @test_throws ArgumentError star3d_front(spike_length=1.1)
    @test_throws ArgumentError star3d_front(petals=2)
end

@testset "VTK export" begin
    tracker = sphere_front(radius=0.5, rings=2, segments=4)
    mktempdir() do dir
        outfile = write_vtk_front(tracker, joinpath(dir, "sphere"), compress=false)
        @test endswith(outfile, ".vtu")
        @test isfile(outfile)
    end
    empty_tracker = FrontTracker()
    @test_throws ArgumentError write_vtk_front(empty_tracker, "dummy")
end
