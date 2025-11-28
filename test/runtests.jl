using Test
using LinearAlgebra
using LibGEOS
using FrontCutTracking

@testset "FrontTracker3D Basic Construction" begin
    # Create an empty 3D interface
    front = FrontTracker3D()
    @test isempty(front.markers)
    @test isempty(front.faces)
    @test front.is_closed == true
    
    # Create interface with markers and faces
    markers = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, 1.0, 0.0), (0.5, 0.5, 1.0)]
    faces = [(1, 2, 3), (1, 2, 4), (2, 3, 4), (3, 1, 4)]
    front = FrontTracker3D(markers, faces)
    @test length(front.markers) == 4
    @test length(front.faces) == 4
    @test front.is_closed == true
end

@testset "FrontTracker3D Marker Operations" begin
    front = FrontTracker3D()
    
    # Test add_marker!
    add_marker!(front, 0.0, 0.0, 0.0)
    add_marker!(front, 1.0, 0.0, 0.0)
    add_marker!(front, 0.5, 1.0, 0.0)
    @test length(get_markers(front)) == 3
    
    # Test set_markers!
    new_markers = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (1.0, 2.0, 0.0), (1.0, 1.0, 2.0)]
    set_markers!(front, new_markers)
    @test length(get_markers(front)) == 4
    
    # Test add_face!
    add_face!(front, 1, 2, 3)
    add_face!(front, 1, 2, 4)
    @test length(get_faces(front)) == 2
    
    # Test set_faces!
    new_faces = [(1, 2, 3), (1, 2, 4), (2, 3, 4), (3, 1, 4)]
    set_faces!(front, new_faces)
    @test length(get_faces(front)) == 4
end

@testset "FrontTracker3D Shape Creation" begin
    @testset "Sphere" begin
        front = FrontTracker3D()
        create_sphere!(front, 0.5, 0.5, 0.5, 0.3, 10, 20)
        markers = get_markers(front)
        faces = get_faces(front)
        
        # Check that we have markers and faces
        @test length(markers) > 0
        @test length(faces) > 0
        
        # Check that markers are approximately on a sphere
        for (x, y, z) in markers
            distance = sqrt((x - 0.5)^2 + (y - 0.5)^2 + (z - 0.5)^2)
            @test isapprox(distance, 0.3, atol=1e-10)
        end
    end
    
    @testset "Ellipsoid" begin
        front = FrontTracker3D()
        create_ellipsoid!(front, 0.5, 0.5, 0.5, 0.3, 0.2, 0.4, 10, 20)
        markers = get_markers(front)
        faces = get_faces(front)
        
        @test length(markers) > 0
        @test length(faces) > 0
        
        # Check that markers are approximately on an ellipsoid
        for (x, y, z) in markers
            normalized_distance = ((x - 0.5)/0.3)^2 + ((y - 0.5)/0.2)^2 + ((z - 0.5)/0.4)^2
            @test isapprox(normalized_distance, 1.0, atol=1e-10)
        end
    end
    
    @testset "Box" begin
        front = FrontTracker3D()
        create_box!(front, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5)
        markers = get_markers(front)
        faces = get_faces(front)
        
        @test length(markers) > 0
        @test length(faces) > 0
        
        # Check that all markers are on the surface of the box
        for (x, y, z) in markers
            on_surface = (isapprox(x, 0.0, atol=1e-10) || isapprox(x, 1.0, atol=1e-10) ||
                         isapprox(y, 0.0, atol=1e-10) || isapprox(y, 1.0, atol=1e-10) ||
                         isapprox(z, 0.0, atol=1e-10) || isapprox(z, 1.0, atol=1e-10))
            @test on_surface
        end
    end
    
    @testset "Cylinder" begin
        front = FrontTracker3D()
        create_cylinder!(front, 0.5, 0.5, 0.0, 1.0, 0.3, 16, 5, true)
        markers = get_markers(front)
        faces = get_faces(front)
        
        @test length(markers) > 0
        @test length(faces) > 0
        @test front.is_closed == true
    end
end

@testset "FrontTracker3D Geometric Operations" begin
    @testset "Point Inside Sphere" begin
        front = FrontTracker3D()
        create_sphere!(front, 0.5, 0.5, 0.5, 0.3, 20, 40)
        
        # Test point inside (center)
        @test is_point_inside(front, 0.5, 0.5, 0.5) == true
        
        # Test points outside
        @test is_point_inside(front, 0.0, 0.0, 0.0) == false
        @test is_point_inside(front, 1.0, 1.0, 1.0) == false
        @test is_point_inside(front, 0.5, 0.5, 1.0) == false
    end
    
    @testset "SDF Calculation" begin
        front = FrontTracker3D()
        create_sphere!(front, 0.0, 0.0, 0.0, 1.0, 20, 40)
        
        # Test center (should be -1.0 approximately)
        center_sdf = sdf(front, 0.0, 0.0, 0.0)
        @test center_sdf < 0
        @test isapprox(center_sdf, -1.0, atol=0.1)
        
        # Test outside (should be positive)
        outside_sdf = sdf(front, 2.0, 0.0, 0.0)
        @test outside_sdf > 0
        @test isapprox(outside_sdf, 1.0, atol=0.1)
    end
    
    @testset "Surface Area" begin
        front = FrontTracker3D()
        radius = 0.5
        create_sphere!(front, 0.0, 0.0, 0.0, radius, 40, 80)
        
        computed_area = compute_total_surface_area(front)
        expected_area = 4π * radius^2
        
        @test isapprox(computed_area, expected_area, rtol=0.05)  # 5% tolerance for discretization
    end
    
    @testset "Enclosed Volume" begin
        front = FrontTracker3D()
        radius = 0.5
        create_sphere!(front, 0.0, 0.0, 0.0, radius, 40, 80)
        
        computed_volume = compute_enclosed_volume(front)
        expected_volume = (4/3) * π * radius^3
        
        @test isapprox(computed_volume, expected_volume, rtol=0.05)  # 5% tolerance for discretization
    end
    
    @testset "Centroid" begin
        center_x, center_y, center_z = 0.3, 0.4, 0.5
        front = FrontTracker3D()
        create_sphere!(front, center_x, center_y, center_z, 0.2, 20, 40)
        
        cx, cy, cz = compute_centroid(front)
        
        @test isapprox(cx, center_x, atol=0.01)
        @test isapprox(cy, center_y, atol=0.01)
        @test isapprox(cz, center_z, atol=0.01)
    end
end

@testset "FrontTracker3D Face and Edge Operations" begin
    front = FrontTracker3D()
    create_sphere!(front, 0.0, 0.0, 0.0, 1.0, 10, 20)
    
    @testset "Face Parameters" begin
        face_normals, face_areas, face_centroids = compute_face_parameters(front)
        
        @test length(face_normals) == length(front.faces)
        @test length(face_areas) == length(front.faces)
        @test length(face_centroids) == length(front.faces)
        
        # All face areas should be positive
        @test all(a -> a > 0, face_areas)
        
        # All normals should be unit vectors
        for n in face_normals
            len = sqrt(n[1]^2 + n[2]^2 + n[3]^2)
            @test isapprox(len, 1.0, atol=1e-6)
        end
    end
    
    @testset "Edge Parameters" begin
        edges, edge_normals, edge_lengths, edge_midpoints = compute_edge_parameters(front)
        
        @test length(edges) > 0
        @test length(edge_normals) == length(edges)
        @test length(edge_lengths) == length(edges)
        @test length(edge_midpoints) == length(edges)
        
        # All edge lengths should be positive
        @test all(l -> l > 0, edge_lengths)
    end
    
    @testset "Marker Normals" begin
        normals = compute_marker_normals(front)
        
        @test length(normals) == length(front.markers)
        
        # All normals should be unit vectors
        for n in normals
            len = sqrt(n[1]^2 + n[2]^2 + n[3]^2)
            @test isapprox(len, 1.0, atol=1e-6)
        end
        
        # For a sphere centered at origin, normals should point radially outward
        for (i, (x, y, z)) in enumerate(front.markers)
            expected_normal = (x, y, z)
            len = sqrt(x^2 + y^2 + z^2)
            if len > 0
                expected_normal = (x/len, y/len, z/len)
                computed_normal = normals[i]
                
                # Dot product should be close to 1 (same direction)
                dot = expected_normal[1]*computed_normal[1] + 
                      expected_normal[2]*computed_normal[2] + 
                      expected_normal[3]*computed_normal[3]
                @test isapprox(abs(dot), 1.0, atol=0.1)
            end
        end
    end
end

@testset "FrontTracker3D Cell Intersection" begin
    front = FrontTracker3D()
    create_sphere!(front, 0.5, 0.5, 0.5, 0.3, 10, 20)
    
    # Create a simple mesh
    x_nodes = collect(0.0:0.2:1.0)
    y_nodes = collect(0.0:0.2:1.0)
    z_nodes = collect(0.0:0.2:1.0)
    
    intersections = compute_face_cell_intersections((x_nodes, y_nodes, z_nodes), front)
    
    # At least some cells should have intersections
    has_intersections = any(!isempty(v) for v in values(intersections))
    @test has_intersections
    
    # Cells far from sphere should have no intersections
    @test isempty(intersections[(1, 1, 1)])  # Corner far from center
end

@testset "FrontTracker3D Volume Jacobian" begin
    front = FrontTracker3D()
    create_sphere!(front, 0.5, 0.5, 0.5, 0.3, 6, 12)
    
    # Create a finer mesh that the sphere will cross
    x_faces = collect(0.0:0.25:1.0)
    y_faces = collect(0.0:0.25:1.0)
    z_faces = collect(0.0:0.25:1.0)
    
    # Compute the volume Jacobian
    jacobian = compute_volume_jacobian_3d(front, x_faces, y_faces, z_faces)
    
    # Basic checks
    @test isa(jacobian, Dict{Tuple{Int, Int, Int}, Vector{Tuple{Int, Float64}}})
    
    # All cells should exist in the dictionary
    nx, ny, nz = length(x_faces)-1, length(y_faces)-1, length(z_faces)-1
    @test length(jacobian) == nx * ny * nz
    
    # At least some cells should have non-zero Jacobian entries
    # (cells where the interface crosses)
    has_nonzero_entries = any(!isempty(v) for v in values(jacobian))
    @test has_nonzero_entries
end

@testset "Julia FrontTracking Tests" begin
    
    @testset "Basic Construction" begin
        # Create an empty interface
        front = FrontTracker()
        @test isempty(front.markers)
        @test front.is_closed == true
        @test front.interface === nothing
        @test front.interface_poly === nothing
        
        # Create interface with markers
        markers = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        front = FrontTracker(markers)
        # The implementation adds a closing point
        @test length(front.markers) == 5  
        @test front.is_closed == true
        @test front.interface !== nothing
        @test front.interface_poly !== nothing
    end
    
    @testset "Shape Creation" begin
        # Test circle creation
        front = FrontTracker()
        create_circle!(front, 0.5, 0.5, 0.3, 20)
        markers = get_markers(front)
        @test length(markers) == 21  # With closing point
        
        # Check that markers are approximately on a circle
        for (x, y) in markers
            distance = sqrt((x - 0.5)^2 + (y - 0.5)^2)
            @test isapprox(distance, 0.3, atol=1e-10)
        end
        
        # Test rectangle creation
        front = FrontTracker()
        create_rectangle!(front, 0.1, 0.2, 0.8, 0.9)
        markers = get_markers(front)
        @test length(markers) == 97

        
        # Test ellipse creation
        front = FrontTracker()
        create_ellipse!(front, 0.5, 0.5, 0.3, 0.2, 20)
        markers = get_markers(front)
        @test length(markers) == 21  # With closing point
        
        # Check that markers are approximately on an ellipse
        for (x, y) in markers
            normalized_distance = ((x - 0.5)/0.3)^2 + ((y - 0.5)/0.2)^2
            @test isapprox(normalized_distance, 1.0, atol=1e-10)
        end
    end
    
    @testset "Point Inside Tests" begin
        # Create a square
        front = FrontTracker()
        create_rectangle!(front, 0.0, 0.0, 1.0, 1.0)
        
        # Test points inside
        @test is_point_inside(front, 0.5, 0.5) == true
        @test is_point_inside(front, 0.1, 0.1) == true
        @test is_point_inside(front, 0.9, 0.9) == true
        
        # Test points outside
        @test is_point_inside(front, -0.5, 0.5) == false
        @test is_point_inside(front, 1.5, 0.5) == false
        @test is_point_inside(front, 0.5, -0.5) == false
        @test is_point_inside(front, 0.5, 1.5) == false
    end
    
    @testset "SDF Calculation" begin
        # Create a square centered at origin
        front = FrontTracker()
        create_rectangle!(front, -1.0, -1.0, 1.0, 1.0)
        
        # Test points inside (should be negative)
        @test sdf(front, 0.0, 0.0) < 0
        @test isapprox(sdf(front, 0.0, 0.0), -1.0, atol=1e-2)
        
        # Test points outside (should be positive)
        @test sdf(front, 2.0, 0.0) > 0
        @test isapprox(sdf(front, 2.0, 0.0), 1.0, atol=1e-2)
        
        # Test points on the boundary (should be approximately zero)
        @test isapprox(sdf(front, 1.0, 0.0), 0.0, atol=1e-2)
        @test isapprox(sdf(front, 0.0, 1.0), 0.0, atol=1e-2)
        
        # Create a circle for additional tests
        front = FrontTracker()
        create_circle!(front, 0.0, 0.0, 1.0)
        
        # Test points at various distances from circle
        @test isapprox(sdf(front, 0.0, 0.0), -1.0, atol=1e-2)  # Center should be -radius
        @test isapprox(sdf(front, 2.0, 0.0), 1.0, atol=1e-2)   # Outside, distance = 1
        @test isapprox(sdf(front, 0.5, 0.0), -0.5, atol=1e-2)  # Inside, distance = 0.5
    end
    

    @testset "Normals and Curvature" begin
        @testset "Circle Normals" begin
            # Create a circle centered at origin
            radius = 0.5
            center_x, center_y = 0.0, 0.0
            front = FrontTracker()
            create_circle!(front, center_x, center_y, radius, 32)
            
            # Calculate normals
            markers = get_markers(front)
            normals = compute_marker_normals(front, markers)
            
            # For a circle, normals should point away from center with unit length
            for i in 1:length(markers)-1 # Skip last point (duplicate of first for closed shape)
                x, y = markers[i]
                nx, ny = normals[i]
                
                # Expected normal: unit vector from center to point
                dist = sqrt((x - center_x)^2 + (y - center_y)^2)
                expected_nx = (x - center_x) / dist
                expected_ny = (y - center_y) / dist
                
                # Check that normal is a unit vector
                @test isapprox(nx^2 + ny^2, 1.0, atol=1e-6)
                
                # Check normal direction
                @test isapprox(nx, expected_nx, atol=1e-6)
                @test isapprox(ny, expected_ny, atol=1e-6)
            end
        end
    end
    
    @testset "Volume Jacobian" begin
        # Create a circular interface
        front = FrontTracker()
        create_circle!(front, 0.5, 0.5, 0.3)
        
        # Create a simple mesh grid
        x_faces = 0.0:0.1:1.0
        y_faces = 0.0:0.1:1.0
        
        # Compute the volume Jacobian
        jacobian = compute_volume_jacobian(front, x_faces, y_faces)
        
        # Basic checks
        @test isa(jacobian, Dict{Tuple{Int, Int}, Vector{Tuple{Int, Float64}}})
        
        # At least some cells should have non-zero Jacobian entries
        @test any(length(values) > 0 for (_, values) in jacobian)
        
        # Cells far from the interface should have zero Jacobian entries
        # Cells at corner (1,1) and (10,10) should be outside the circle
        @test length(get(jacobian, (1, 1), [])) == 0
        @test length(get(jacobian, (10, 10), [])) == 0
        
        # Cells near the interface should have non-zero Jacobian entries
        # Find where the interface crosses cells
        has_nonzero_entries = false
        for i in 1:length(x_faces)-1
            for j in 1:length(y_faces)-1
                if haskey(jacobian, (i, j)) && !isempty(jacobian[(i, j)])
                    has_nonzero_entries = true
                    break
                end
            end
            if has_nonzero_entries
                break
            end
        end
        @test has_nonzero_entries

    end
end

@testset "FrontTracker1D Basic Functionality" begin
    # Test constructors
    ft_empty = FrontTracker1D()
    @test isempty(get_markers(ft_empty))

    ft = FrontTracker1D([0.3, 0.7])
    @test get_markers(ft) == [0.3, 0.7]

    # Test add_marker!
    add_marker!(ft, 0.5)
    @test get_markers(ft) == sort([0.3, 0.5, 0.7])

    # Test set_markers!
    set_markers!(ft, [0.1, 0.9])
    @test get_markers(ft) == [0.1, 0.9]

    # Test is_point_inside
    @test !is_point_inside(ft, 0.0)
    @test is_point_inside(ft, 0.5)
    @test !is_point_inside(ft, 1.0)

    # Test sdf
    @test sdf(ft, 0.0) > 0
    @test sdf(ft, 0.5) < 0
    @test sdf(ft, 1.0) > 0
end

@testset "FrontTracker1D Capacity Calculation" begin
    # Simple mesh and front
    nx = 10
    lx = 1.0
    x_nodes = collect(range(0.0, stop=lx, length=nx+1))
    ft = FrontTracker1D([0.3, 0.7])

    caps = compute_capacities_1d((x_nodes,), ft)
    @test length(caps[:fractions]) == nx+1
    @test length(caps[:volumes]) == nx+1
    @test length(caps[:centroids_x]) == nx+1
    @test length(caps[:cell_types]) == nx+1
    @test length(caps[:Ax]) == nx+1
    @test length(caps[:Wx]) == nx+1
    @test length(caps[:Bx]) == nx+1

    # Check that fluid fractions are between 0 and 1
    @test all(0 .<= caps[:fractions] .<= 1)
end

@testset "FrontTracker1D Space-Time Capacity Calculation" begin
    nx = 10
    lx = 1.0
    x_nodes = collect(range(0.0, stop=lx, length=nx+1))
    ft_n = FrontTracker1D([0.3, 0.7])
    ft_np1 = FrontTracker1D([0.35, 0.75])
    dt = 0.1

    st_caps = compute_spacetime_capacities_1d((x_nodes,), ft_n, ft_np1, dt)
    @test length(st_caps[:Ax_spacetime]) == nx+1
    @test length(st_caps[:V_spacetime]) == nx+1
    @test length(st_caps[:Bx_spacetime]) == nx+1
    @test length(st_caps[:Wx_spacetime]) == nx+1
    @test length(st_caps[:ST_centroids]) == nx
    @test length(st_caps[:ms_cases]) == nx
    @test length(st_caps[:edge_types]) == nx+1
    @test length(st_caps[:t_crosses]) == nx+1
    @test length(st_caps[:center_types]) == nx
    @test length(st_caps[:t_crosses_center]) == nx
end

@testset "Geometric Metrics" begin
    center_x, center_y = 0.45, 0.55
    radius = 0.27
    front = FrontTracker()
    create_circle!(front, center_x, center_y, radius, 720)

    poly = get_fluid_polygon(front)
    @test !isnothing(poly)

    area = LibGEOS.area(poly)
    expected_area = π * radius^2
    @test isapprox(area, expected_area; atol=1e-3)

    centroid_geom = LibGEOS.centroid(poly)
    centroid = (
        LibGEOS.getcoord(centroid_geom, 1),
        LibGEOS.getcoord(centroid_geom, 2),
    )
    @test isapprox(centroid[1], center_x; atol=1e-3)
    @test isapprox(centroid[2], center_y; atol=1e-3)

    expected_perimeter = 2π * radius

    segments, _, _, segment_lengths, _ = compute_segment_parameters(front)
    discrete_perimeter = sum(segment_lengths)
    @test isapprox(discrete_perimeter, expected_perimeter; atol=1e-3)

    nx = 64
    ny = 64
    x_nodes = collect(range(0.0, stop=1.0, length=nx+1))
    y_nodes = collect(range(0.0, stop=1.0, length=ny+1))
    capacities = compute_capacities((x_nodes, y_nodes), front)
    volume_sum = sum(capacities[:volumes])
    @test isapprox(volume_sum, expected_area; atol=5e-3)

    vol = capacities[:volumes]
    cx = capacities[:centroids_x]
    cy = capacities[:centroids_y]
    total_mass = sum(vol)
    centroid_x = sum(vol .* cx) / total_mass
    centroid_y = sum(vol .* cy) / total_mass
    @test isapprox(centroid_x, center_x; atol=5e-3)
    @test isapprox(centroid_y, center_y; atol=5e-3)
end
