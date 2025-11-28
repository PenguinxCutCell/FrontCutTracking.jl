"""
    compute_edge_parameters(ft::FrontTracker3D)

Computes parameters for each edge of the interface mesh:
- edges: List of vertex index pairs defining each edge
- edge_normals: Average normal at each edge midpoint
- edge_lengths: Length of each edge
- edge_midpoints: Midpoint coordinates of each edge

Returns (edges, edge_normals, edge_lengths, edge_midpoints)
"""
function compute_edge_parameters(ft::FrontTracker3D)
    if isempty(ft.faces) || isempty(ft.markers)
        return [], [], [], []
    end
    
    # Build unique edges from faces
    edge_set = Set{Tuple{Int, Int}}()
    edge_to_faces = Dict{Tuple{Int, Int}, Vector{Int}}()
    
    for (face_idx, (i1, i2, i3)) in enumerate(ft.faces)
        # Each triangle has 3 edges
        for (a, b) in [(i1, i2), (i2, i3), (i3, i1)]
            # Store edges with smaller index first
            edge = a < b ? (a, b) : (b, a)
            push!(edge_set, edge)
            
            if !haskey(edge_to_faces, edge)
                edge_to_faces[edge] = Int[]
            end
            push!(edge_to_faces[edge], face_idx)
        end
    end
    
    edges = collect(edge_set)
    n_edges = length(edges)
    
    edge_normals = Vector{Tuple{Float64, Float64, Float64}}(undef, n_edges)
    edge_lengths = Vector{Float64}(undef, n_edges)
    edge_midpoints = Vector{Tuple{Float64, Float64, Float64}}(undef, n_edges)
    
    for (edge_idx, (i1, i2)) in enumerate(edges)
        v1 = ft.markers[i1]
        v2 = ft.markers[i2]
        
        # Edge vector
        edge_vec = (v2[1] - v1[1], v2[2] - v1[2], v2[3] - v1[3])
        
        # Edge length
        edge_lengths[edge_idx] = sqrt(dot3(edge_vec, edge_vec))
        
        # Edge midpoint
        edge_midpoints[edge_idx] = (
            (v1[1] + v2[1]) / 2,
            (v1[2] + v2[2]) / 2,
            (v1[3] + v2[3]) / 2
        )
        
        # Average normal from adjacent faces
        normal = (0.0, 0.0, 0.0)
        adj_faces = get(edge_to_faces, (min(i1, i2), max(i1, i2)), Int[])
        
        for face_idx in adj_faces
            face_normal = compute_face_normal(ft, face_idx)
            face_area = compute_face_area(ft, face_idx)
            normal = (normal[1] + face_normal[1] * face_area,
                     normal[2] + face_normal[2] * face_area,
                     normal[3] + face_normal[3] * face_area)
        end
        
        # Normalize
        len = sqrt(dot3(normal, normal))
        if len > 1e-10
            edge_normals[edge_idx] = (normal[1] / len, normal[2] / len, normal[3] / len)
        else
            edge_normals[edge_idx] = (0.0, 0.0, 1.0)
        end
    end
    
    return edges, edge_normals, edge_lengths, edge_midpoints
end

"""
    compute_face_parameters(ft::FrontTracker3D)

Computes parameters for each triangular face of the interface:
- face_normals: Normal vector for each face
- face_areas: Area of each face
- face_centroids: Centroid of each face

Returns (face_normals, face_areas, face_centroids)
"""
function compute_face_parameters(ft::FrontTracker3D)
    if isempty(ft.faces) || isempty(ft.markers)
        return [], [], []
    end
    
    n_faces = length(ft.faces)
    face_normals = Vector{Tuple{Float64, Float64, Float64}}(undef, n_faces)
    face_areas = Vector{Float64}(undef, n_faces)
    face_centroids = Vector{Tuple{Float64, Float64, Float64}}(undef, n_faces)
    
    for (face_idx, (i1, i2, i3)) in enumerate(ft.faces)
        v1 = ft.markers[i1]
        v2 = ft.markers[i2]
        v3 = ft.markers[i3]
        
        # Face centroid
        face_centroids[face_idx] = (
            (v1[1] + v2[1] + v3[1]) / 3,
            (v1[2] + v2[2] + v3[2]) / 3,
            (v1[3] + v2[3] + v3[3]) / 3
        )
        
        # Face normal and area
        face_normals[face_idx] = compute_face_normal(ft, face_idx)
        face_areas[face_idx] = compute_face_area(ft, face_idx)
    end
    
    return face_normals, face_areas, face_centroids
end

"""
    get_adjacent_faces(ft::FrontTracker3D, marker_idx::Int)

Returns indices of all faces that include the specified marker.
"""
function get_adjacent_faces(ft::FrontTracker3D, marker_idx::Int)
    adjacent = Int[]
    for (face_idx, (i1, i2, i3)) in enumerate(ft.faces)
        if i1 == marker_idx || i2 == marker_idx || i3 == marker_idx
            push!(adjacent, face_idx)
        end
    end
    return adjacent
end

"""
    get_adjacent_markers(ft::FrontTracker3D, marker_idx::Int)

Returns indices of all markers directly connected to the specified marker by an edge.
"""
function get_adjacent_markers(ft::FrontTracker3D, marker_idx::Int)
    adjacent = Set{Int}()
    for (i1, i2, i3) in ft.faces
        if i1 == marker_idx
            push!(adjacent, i2)
            push!(adjacent, i3)
        elseif i2 == marker_idx
            push!(adjacent, i1)
            push!(adjacent, i3)
        elseif i3 == marker_idx
            push!(adjacent, i1)
            push!(adjacent, i2)
        end
    end
    return collect(adjacent)
end

"""
    compute_vertex_area(ft::FrontTracker3D, marker_idx::Int)

Computes the Voronoi area associated with a vertex.
This is the sum of 1/3 of each adjacent triangle's area.
"""
function compute_vertex_area(ft::FrontTracker3D, marker_idx::Int)
    total_area = 0.0
    adjacent_faces = get_adjacent_faces(ft, marker_idx)
    
    for face_idx in adjacent_faces
        total_area += compute_face_area(ft, face_idx) / 3.0
    end
    
    return total_area
end

"""
    compute_mean_curvature(ft::FrontTracker3D, marker_idx::Int)

Computes the mean curvature at a vertex using the Laplace-Beltrami operator.
"""
function compute_mean_curvature(ft::FrontTracker3D, marker_idx::Int)
    vertex = ft.markers[marker_idx]
    adjacent_markers = get_adjacent_markers(ft, marker_idx)
    
    if isempty(adjacent_markers)
        return 0.0
    end
    
    # Compute discrete Laplace-Beltrami operator
    laplacian = (0.0, 0.0, 0.0)
    total_weight = 0.0
    
    for adj_idx in adjacent_markers
        adj_vertex = ft.markers[adj_idx]
        
        # Simple uniform weighting (cotangent weights would be more accurate)
        weight = 1.0
        diff = (adj_vertex[1] - vertex[1], adj_vertex[2] - vertex[2], adj_vertex[3] - vertex[3])
        
        laplacian = (laplacian[1] + weight * diff[1],
                    laplacian[2] + weight * diff[2],
                    laplacian[3] + weight * diff[3])
        total_weight += weight
    end
    
    if total_weight > 0
        laplacian = (laplacian[1] / total_weight,
                    laplacian[2] / total_weight,
                    laplacian[3] / total_weight)
    end
    
    # Mean curvature is half the magnitude of the Laplacian
    return 0.5 * sqrt(dot3(laplacian, laplacian))
end
