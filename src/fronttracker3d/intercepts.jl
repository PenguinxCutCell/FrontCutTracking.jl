"""
    compute_face_cell_intersections(nodes::NTuple{3, AbstractVector}, ft::FrontTracker3D)

Computes intersections between triangular faces of the interface and cells of a 3D mesh.
Returns a dictionary mapping cell indices (i,j,k) to lists of (face_idx, intersection_area).
"""
function compute_face_cell_intersections(nodes::NTuple{3, AbstractVector}, ft::FrontTracker3D)
    x_nodes = nodes[1]
    y_nodes = nodes[2]
    z_nodes = nodes[3]
    nx = length(x_nodes) - 1
    ny = length(y_nodes) - 1
    nz = length(z_nodes) - 1
    
    # Initialize dictionary for all cells
    cell_face_intersections = Dict{Tuple{Int, Int, Int}, Vector{Tuple{Int, Float64}}}()
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                cell_face_intersections[(i, j, k)] = []
            end
        end
    end
    
    # For each face, find intersecting cells
    for (face_idx, (i1, i2, i3)) in enumerate(ft.faces)
        v1 = ft.markers[i1]
        v2 = ft.markers[i2]
        v3 = ft.markers[i3]
        
        # Get bounding box of triangle
        min_x = min(v1[1], v2[1], v3[1])
        max_x = max(v1[1], v2[1], v3[1])
        min_y = min(v1[2], v2[2], v3[2])
        max_y = max(v1[2], v2[2], v3[2])
        min_z = min(v1[3], v2[3], v3[3])
        max_z = max(v1[3], v2[3], v3[3])
        
        # Find cell range that might intersect
        i_min = max(1, searchsortedlast(x_nodes, min_x))
        i_max = min(nx, searchsortedfirst(x_nodes, max_x))
        j_min = max(1, searchsortedlast(y_nodes, min_y))
        j_max = min(ny, searchsortedfirst(y_nodes, max_y))
        k_min = max(1, searchsortedlast(z_nodes, min_z))
        k_max = min(nz, searchsortedfirst(z_nodes, max_z))
        
        # Check each potentially intersecting cell
        for i in i_min:i_max
            for j in j_min:j_max
                for k in k_min:k_max
                    cell_bounds = (
                        (x_nodes[i], x_nodes[i+1]),
                        (y_nodes[j], y_nodes[j+1]),
                        (z_nodes[k], z_nodes[k+1])
                    )
                    
                    # Compute intersection area
                    intersection_area = triangle_box_intersection_area(v1, v2, v3, cell_bounds)
                    
                    if intersection_area > 1e-10
                        push!(cell_face_intersections[(i, j, k)], (face_idx, intersection_area))
                    end
                end
            end
        end
    end
    
    return cell_face_intersections
end

"""
    triangle_box_intersection_area(v1, v2, v3, cell_bounds)

Computes the area of intersection between a triangle and an axis-aligned box.
Uses Sutherland-Hodgman polygon clipping extended to 3D.
"""
function triangle_box_intersection_area(v1::Tuple{Float64, Float64, Float64},
                                        v2::Tuple{Float64, Float64, Float64},
                                        v3::Tuple{Float64, Float64, Float64},
                                        cell_bounds::Tuple{Tuple{Float64, Float64}, 
                                                         Tuple{Float64, Float64}, 
                                                         Tuple{Float64, Float64}})
    # Start with triangle vertices as polygon
    polygon = [v1, v2, v3]
    
    # Clip against each of the 6 planes of the box
    x_min, x_max = cell_bounds[1]
    y_min, y_max = cell_bounds[2]
    z_min, z_max = cell_bounds[3]
    
    # Clip against x_min plane
    polygon = clip_polygon_plane(polygon, (1.0, 0.0, 0.0), x_min)
    if isempty(polygon) return 0.0 end
    
    # Clip against x_max plane
    polygon = clip_polygon_plane(polygon, (-1.0, 0.0, 0.0), -x_max)
    if isempty(polygon) return 0.0 end
    
    # Clip against y_min plane
    polygon = clip_polygon_plane(polygon, (0.0, 1.0, 0.0), y_min)
    if isempty(polygon) return 0.0 end
    
    # Clip against y_max plane
    polygon = clip_polygon_plane(polygon, (0.0, -1.0, 0.0), -y_max)
    if isempty(polygon) return 0.0 end
    
    # Clip against z_min plane
    polygon = clip_polygon_plane(polygon, (0.0, 0.0, 1.0), z_min)
    if isempty(polygon) return 0.0 end
    
    # Clip against z_max plane
    polygon = clip_polygon_plane(polygon, (0.0, 0.0, -1.0), -z_max)
    if isempty(polygon) return 0.0 end
    
    # Calculate area of clipped polygon
    return polygon_area_3d(polygon)
end

"""
    clip_polygon_plane(polygon, normal, d)

Clips a polygon against a half-plane defined by normal · p >= d.
Uses Sutherland-Hodgman algorithm.
"""
function clip_polygon_plane(polygon::Vector{Tuple{Float64, Float64, Float64}},
                            normal::Tuple{Float64, Float64, Float64},
                            d::Float64)
    if isempty(polygon)
        return polygon
    end
    
    output = Tuple{Float64, Float64, Float64}[]
    n = length(polygon)
    
    for i in 1:n
        current = polygon[i]
        next = polygon[mod1(i + 1, n)]
        
        current_inside = dot3(normal, current) >= d - 1e-10
        next_inside = dot3(normal, next) >= d - 1e-10
        
        if current_inside
            push!(output, current)
            if !next_inside
                # Compute intersection
                intersection = line_plane_intersection(current, next, normal, d)
                push!(output, intersection)
            end
        elseif next_inside
            # Compute intersection
            intersection = line_plane_intersection(current, next, normal, d)
            push!(output, intersection)
        end
    end
    
    return output
end

"""
    line_plane_intersection(p1, p2, normal, d)

Finds the intersection point of a line segment with a plane.
"""
function line_plane_intersection(p1::Tuple{Float64, Float64, Float64},
                                  p2::Tuple{Float64, Float64, Float64},
                                  normal::Tuple{Float64, Float64, Float64},
                                  d::Float64)
    direction = (p2[1] - p1[1], p2[2] - p1[2], p2[3] - p1[3])
    denom = dot3(normal, direction)
    
    if abs(denom) < 1e-10
        return p1  # Line parallel to plane
    end
    
    t = (d - dot3(normal, p1)) / denom
    t = clamp(t, 0.0, 1.0)
    
    return (p1[1] + t * direction[1],
            p1[2] + t * direction[2],
            p1[3] + t * direction[3])
end

"""
    polygon_area_3d(polygon)

Calculates the area of a 3D polygon by projecting to 2D.
"""
function polygon_area_3d(polygon::Vector{Tuple{Float64, Float64, Float64}})
    n = length(polygon)
    if n < 3
        return 0.0
    end
    
    # Compute normal of polygon
    v1 = polygon[1]
    normal = (0.0, 0.0, 0.0)
    
    for i in 2:n-1
        e1 = (polygon[i][1] - v1[1], polygon[i][2] - v1[2], polygon[i][3] - v1[3])
        e2 = (polygon[i+1][1] - v1[1], polygon[i+1][2] - v1[2], polygon[i+1][3] - v1[3])
        cross = cross3(e1, e2)
        normal = (normal[1] + cross[1], normal[2] + cross[2], normal[3] + cross[3])
    end
    
    # The magnitude of the cross product sum / 2 gives the area
    return 0.5 * sqrt(dot3(normal, normal))
end

"""
    compute_intercept_jacobian_3d(nodes::NTuple{3, AbstractVector}, ft::FrontTracker3D; density::Float64=1.0)

Computes the Jacobian of volumes with respect to intercept displacements in 3D.
For each cell (i,j,k) and each face F, J[cell, F] = ∂V_cell/∂δ_F = ρL × A_cell,F,
where A_cell,F is the intersection area of face F with the cell.

Returns:
- intercept_jacobian: Dict mapping cell indices to list of (face_idx, jacobian_value)
- face_normals: Normal vectors for each face
- face_areas: Areas of each face
- face_centroids: Centroids of each face
"""
function compute_intercept_jacobian_3d(nodes::NTuple{3, AbstractVector}, ft::FrontTracker3D; density::Float64=1.0)
    # Compute face-cell intersections
    cell_face_intersections = compute_face_cell_intersections(nodes, ft)
    
    # Compute face parameters
    face_normals, face_areas, face_centroids = compute_face_parameters(ft)
    
    x_nodes = nodes[1]
    y_nodes = nodes[2]
    z_nodes = nodes[3]
    nx = length(x_nodes) - 1
    ny = length(y_nodes) - 1
    nz = length(z_nodes) - 1
    
    # Create Jacobian dictionary
    intercept_jacobian = Dict{Tuple{Int, Int, Int}, Vector{Tuple{Int, Float64}}}()
    
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                intercept_jacobian[(i, j, k)] = []
                
                for (face_idx, intersection_area) in cell_face_intersections[(i, j, k)]
                    # J[cell, F] = ρL × intersection_area
                    jacobian_value = density * intersection_area
                    push!(intercept_jacobian[(i, j, k)], (face_idx, jacobian_value))
                end
            end
        end
    end
    
    return intercept_jacobian, face_normals, face_areas, face_centroids
end

"""
    update_front_with_intercept_displacements_3d!(ft::FrontTracker3D, 
                                                   displacements::AbstractVector{<:Real},
                                                   face_normals::Vector{Tuple{Float64, Float64, Float64}},
                                                   face_areas::Vector{Float64})

Updates the 3D interface by displacing each face according to its intercept displacement.
Applies area-weighted averaging to distribute displacements to shared vertices.
"""
function update_front_with_intercept_displacements_3d!(ft::FrontTracker3D, 
                                                       displacements::AbstractVector{<:Real},
                                                       face_normals::Vector{Tuple{Float64, Float64, Float64}},
                                                       face_areas::Vector{Float64})
    markers = copy(ft.markers)
    n_markers = length(markers)
    n_faces = length(displacements)
    
    # Structure to store weighted contributions to each marker
    marker_contributions = Dict{Int, Vector{Tuple{Float64, Tuple{Float64, Float64, Float64}}}}()
    for i in 1:n_markers
        marker_contributions[i] = []
    end
    
    # Accumulate contributions from each face
    for (face_idx, (i1, i2, i3)) in enumerate(ft.faces)
        if face_idx > n_faces
            break
        end
        
        displacement = displacements[face_idx]
        normal = face_normals[face_idx]
        area = max(face_areas[face_idx], 1e-10)
        
        vector_displacement = (displacement * normal[1],
                              displacement * normal[2],
                              displacement * normal[3])
        
        # Use face area as weight
        weight = area
        
        push!(marker_contributions[i1], (weight, vector_displacement))
        push!(marker_contributions[i2], (weight, vector_displacement))
        push!(marker_contributions[i3], (weight, vector_displacement))
    end
    
    # Apply weighted average displacement to each marker
    for i in 1:n_markers
        contributions = marker_contributions[i]
        if !isempty(contributions)
            total_weight = sum(c[1] for c in contributions)
            
            if total_weight > 0
                avg_dx = sum(c[1] * c[2][1] for c in contributions) / total_weight
                avg_dy = sum(c[1] * c[2][2] for c in contributions) / total_weight
                avg_dz = sum(c[1] * c[2][3] for c in contributions) / total_weight
                
                markers[i] = (markers[i][1] + avg_dx,
                             markers[i][2] + avg_dy,
                             markers[i][3] + avg_dz)
            end
        end
    end
    
    ft.markers = markers
    return ft
end
