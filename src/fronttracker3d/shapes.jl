"""
    create_sphere!(ft::FrontTracker3D, center_x::Float64, center_y::Float64, center_z::Float64, 
                   radius::Float64, n_lat::Int=20, n_lon::Int=40)

Creates a spherical interface using a UV sphere parameterization.
- center_x, center_y, center_z: Center coordinates of the sphere
- radius: Radius of the sphere
- n_lat: Number of latitude subdivisions
- n_lon: Number of longitude subdivisions
"""
function create_sphere!(ft::FrontTracker3D, center_x::Float64, center_y::Float64, center_z::Float64, 
                        radius::Float64, n_lat::Int=20, n_lon::Int=40)
    markers = Vector{Tuple{Float64, Float64, Float64}}()
    faces = Vector{Tuple{Int, Int, Int}}()
    
    # Add north pole
    push!(markers, (center_x, center_y, center_z + radius))
    north_pole_idx = 1
    
    # Generate latitude rings (excluding poles)
    for i in 1:n_lat-1
        theta = π * i / n_lat  # From 0 to π
        z = center_z + radius * cos(theta)
        ring_radius = radius * sin(theta)
        
        for j in 0:n_lon-1
            phi = 2π * j / n_lon  # From 0 to 2π
            x = center_x + ring_radius * cos(phi)
            y = center_y + ring_radius * sin(phi)
            push!(markers, (x, y, z))
        end
    end
    
    # Add south pole
    push!(markers, (center_x, center_y, center_z - radius))
    south_pole_idx = length(markers)
    
    # Create faces for north pole cap
    for j in 0:n_lon-1
        j_next = (j + 1) % n_lon
        # First ring starts at index 2
        push!(faces, (north_pole_idx, 2 + j, 2 + j_next))
    end
    
    # Create faces for middle rings
    for i in 0:n_lat-3
        ring_start = 2 + i * n_lon
        next_ring_start = ring_start + n_lon
        
        for j in 0:n_lon-1
            j_next = (j + 1) % n_lon
            
            # Two triangles per quad
            v1 = ring_start + j
            v2 = ring_start + j_next
            v3 = next_ring_start + j
            v4 = next_ring_start + j_next
            
            push!(faces, (v1, v3, v2))
            push!(faces, (v2, v3, v4))
        end
    end
    
    # Create faces for south pole cap
    last_ring_start = 2 + (n_lat - 2) * n_lon
    for j in 0:n_lon-1
        j_next = (j + 1) % n_lon
        push!(faces, (last_ring_start + j, south_pole_idx, last_ring_start + j_next))
    end
    
    ft.markers = markers
    ft.faces = faces
    ft.is_closed = true
    
    return ft
end

"""
    create_box!(ft::FrontTracker3D, min_x::Float64, min_y::Float64, min_z::Float64,
                max_x::Float64, max_y::Float64, max_z::Float64, n_per_edge::Int=10)

Creates a box/cuboid interface with markers distributed on the surface.
- min_x, min_y, min_z: Minimum corner coordinates
- max_x, max_y, max_z: Maximum corner coordinates
- n_per_edge: Number of markers per edge (for refinement)
"""
function create_box!(ft::FrontTracker3D, min_x::Float64, min_y::Float64, min_z::Float64,
                     max_x::Float64, max_y::Float64, max_z::Float64, n_per_edge::Int=10)
    markers = Vector{Tuple{Float64, Float64, Float64}}()
    faces = Vector{Tuple{Int, Int, Int}}()
    
    # Calculate step sizes
    dx = (max_x - min_x) / n_per_edge
    dy = (max_y - min_y) / n_per_edge
    dz = (max_z - min_z) / n_per_edge
    
    # Helper function to add a rectangular face with triangulation
    function add_rect_face!(markers, faces, corners, nx, ny, get_point)
        start_idx = length(markers)
        
        # Add vertices in a grid
        for j in 0:ny
            for i in 0:nx
                push!(markers, get_point(i/nx, j/ny))
            end
        end
        
        # Add triangular faces
        for j in 0:ny-1
            for i in 0:nx-1
                v1 = start_idx + j * (nx + 1) + i + 1
                v2 = start_idx + j * (nx + 1) + i + 2
                v3 = start_idx + (j + 1) * (nx + 1) + i + 1
                v4 = start_idx + (j + 1) * (nx + 1) + i + 2
                
                push!(faces, (v1, v3, v2))
                push!(faces, (v2, v3, v4))
            end
        end
    end
    
    # Bottom face (z = min_z) - normal pointing -z
    add_rect_face!(markers, faces, nothing, n_per_edge, n_per_edge, 
        (u, v) -> (min_x + u * (max_x - min_x), min_y + v * (max_y - min_y), min_z))
    
    # Top face (z = max_z) - normal pointing +z
    add_rect_face!(markers, faces, nothing, n_per_edge, n_per_edge, 
        (u, v) -> (min_x + u * (max_x - min_x), max_y - v * (max_y - min_y), max_z))
    
    # Front face (y = min_y) - normal pointing -y
    add_rect_face!(markers, faces, nothing, n_per_edge, n_per_edge, 
        (u, v) -> (min_x + u * (max_x - min_x), min_y, min_z + v * (max_z - min_z)))
    
    # Back face (y = max_y) - normal pointing +y
    add_rect_face!(markers, faces, nothing, n_per_edge, n_per_edge, 
        (u, v) -> (max_x - u * (max_x - min_x), max_y, min_z + v * (max_z - min_z)))
    
    # Left face (x = min_x) - normal pointing -x
    add_rect_face!(markers, faces, nothing, n_per_edge, n_per_edge, 
        (u, v) -> (min_x, max_y - u * (max_y - min_y), min_z + v * (max_z - min_z)))
    
    # Right face (x = max_x) - normal pointing +x
    add_rect_face!(markers, faces, nothing, n_per_edge, n_per_edge, 
        (u, v) -> (max_x, min_y + u * (max_y - min_y), min_z + v * (max_z - min_z)))
    
    ft.markers = markers
    ft.faces = faces
    ft.is_closed = true
    
    return ft
end

"""
    create_ellipsoid!(ft::FrontTracker3D, center_x::Float64, center_y::Float64, center_z::Float64, 
                      radius_x::Float64, radius_y::Float64, radius_z::Float64, 
                      n_lat::Int=20, n_lon::Int=40)

Creates an ellipsoidal interface.
- center_x, center_y, center_z: Center coordinates
- radius_x, radius_y, radius_z: Semi-axes of the ellipsoid
- n_lat: Number of latitude subdivisions
- n_lon: Number of longitude subdivisions
"""
function create_ellipsoid!(ft::FrontTracker3D, center_x::Float64, center_y::Float64, center_z::Float64, 
                           radius_x::Float64, radius_y::Float64, radius_z::Float64, 
                           n_lat::Int=20, n_lon::Int=40)
    markers = Vector{Tuple{Float64, Float64, Float64}}()
    faces = Vector{Tuple{Int, Int, Int}}()
    
    # Add north pole
    push!(markers, (center_x, center_y, center_z + radius_z))
    north_pole_idx = 1
    
    # Generate latitude rings (excluding poles)
    for i in 1:n_lat-1
        theta = π * i / n_lat  # From 0 to π
        z = center_z + radius_z * cos(theta)
        ring_scale = sin(theta)
        
        for j in 0:n_lon-1
            phi = 2π * j / n_lon  # From 0 to 2π
            x = center_x + radius_x * ring_scale * cos(phi)
            y = center_y + radius_y * ring_scale * sin(phi)
            push!(markers, (x, y, z))
        end
    end
    
    # Add south pole
    push!(markers, (center_x, center_y, center_z - radius_z))
    south_pole_idx = length(markers)
    
    # Create faces for north pole cap
    for j in 0:n_lon-1
        j_next = (j + 1) % n_lon
        push!(faces, (north_pole_idx, 2 + j, 2 + j_next))
    end
    
    # Create faces for middle rings
    for i in 0:n_lat-3
        ring_start = 2 + i * n_lon
        next_ring_start = ring_start + n_lon
        
        for j in 0:n_lon-1
            j_next = (j + 1) % n_lon
            
            v1 = ring_start + j
            v2 = ring_start + j_next
            v3 = next_ring_start + j
            v4 = next_ring_start + j_next
            
            push!(faces, (v1, v3, v2))
            push!(faces, (v2, v3, v4))
        end
    end
    
    # Create faces for south pole cap
    last_ring_start = 2 + (n_lat - 2) * n_lon
    for j in 0:n_lon-1
        j_next = (j + 1) % n_lon
        push!(faces, (last_ring_start + j, south_pole_idx, last_ring_start + j_next))
    end
    
    ft.markers = markers
    ft.faces = faces
    ft.is_closed = true
    
    return ft
end

"""
    create_cylinder!(ft::FrontTracker3D, center_x::Float64, center_y::Float64, z_min::Float64, z_max::Float64,
                     radius::Float64, n_radial::Int=32, n_height::Int=10, capped::Bool=true)

Creates a cylindrical interface.
- center_x, center_y: Center coordinates in the xy-plane
- z_min, z_max: Bottom and top z coordinates
- radius: Radius of the cylinder
- n_radial: Number of radial subdivisions
- n_height: Number of height subdivisions
- capped: Whether to include top and bottom caps
"""
function create_cylinder!(ft::FrontTracker3D, center_x::Float64, center_y::Float64, z_min::Float64, z_max::Float64,
                          radius::Float64, n_radial::Int=32, n_height::Int=10, capped::Bool=true)
    markers = Vector{Tuple{Float64, Float64, Float64}}()
    faces = Vector{Tuple{Int, Int, Int}}()
    
    dz = (z_max - z_min) / n_height
    
    # Add lateral surface vertices
    for i in 0:n_height
        z = z_min + i * dz
        for j in 0:n_radial-1
            phi = 2π * j / n_radial
            x = center_x + radius * cos(phi)
            y = center_y + radius * sin(phi)
            push!(markers, (x, y, z))
        end
    end
    
    # Create lateral surface faces
    for i in 0:n_height-1
        ring_start = i * n_radial + 1
        next_ring_start = ring_start + n_radial
        
        for j in 0:n_radial-1
            j_next = (j + 1) % n_radial
            
            v1 = ring_start + j
            v2 = ring_start + j_next
            v3 = next_ring_start + j
            v4 = next_ring_start + j_next
            
            push!(faces, (v1, v3, v2))
            push!(faces, (v2, v3, v4))
        end
    end
    
    if capped
        # Add bottom cap center
        push!(markers, (center_x, center_y, z_min))
        bottom_center_idx = length(markers)
        
        # Bottom cap faces
        for j in 0:n_radial-1
            j_next = (j + 1) % n_radial
            push!(faces, (bottom_center_idx, 1 + j_next, 1 + j))
        end
        
        # Add top cap center
        push!(markers, (center_x, center_y, z_max))
        top_center_idx = length(markers)
        
        # Top cap faces
        top_ring_start = n_height * n_radial + 1
        for j in 0:n_radial-1
            j_next = (j + 1) % n_radial
            push!(faces, (top_center_idx, top_ring_start + j, top_ring_start + j_next))
        end
        
        ft.is_closed = true
    else
        ft.is_closed = false
    end
    
    ft.markers = markers
    ft.faces = faces
    
    return ft
end
