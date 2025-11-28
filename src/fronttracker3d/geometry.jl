"""
    is_point_inside(ft::FrontTracker3D, x::Float64, y::Float64, z::Float64)

Checks if a point is inside the 3D interface using the signed solid angle method.
This method computes the total solid angle subtended by the mesh from the point.
For a closed surface, the solid angle is 4π if inside, 0 if outside.
"""
function is_point_inside(ft::FrontTracker3D, x::Float64, y::Float64, z::Float64)
    if isempty(ft.faces) || isempty(ft.markers)
        return false
    end
    
    point = (x, y, z)
    total_solid_angle = 0.0
    
    for (i1, i2, i3) in ft.faces
        v1 = ft.markers[i1]
        v2 = ft.markers[i2]
        v3 = ft.markers[i3]
        
        # Compute vectors from point to vertices
        a = (v1[1] - x, v1[2] - y, v1[3] - z)
        b = (v2[1] - x, v2[2] - y, v2[3] - z)
        c = (v3[1] - x, v3[2] - y, v3[3] - z)
        
        # Compute distances
        la = sqrt(dot3(a, a))
        lb = sqrt(dot3(b, b))
        lc = sqrt(dot3(c, c))
        
        if la < 1e-10 || lb < 1e-10 || lc < 1e-10
            # Point is on or very close to vertex
            return true
        end
        
        # Normalize vectors
        a = (a[1]/la, a[2]/la, a[3]/la)
        b = (b[1]/lb, b[2]/lb, b[3]/lb)
        c = (c[1]/lc, c[2]/lc, c[3]/lc)
        
        # Compute solid angle contribution using the formula:
        # Ω = 2 * atan2(det, trace) where det = a · (b × c) and 
        # trace = 1 + a·b + b·c + c·a
        numerator = dot3(a, cross3(b, c))
        denominator = 1.0 + dot3(a, b) + dot3(b, c) + dot3(c, a)
        
        total_solid_angle += 2.0 * atan(numerator, denominator)
    end
    
    # Inside if solid angle is approximately 4π (or -4π depending on face orientation)
    # Use a threshold of 2π (halfway between 0 and 4π)
    return abs(total_solid_angle) > 2π
end

"""
    ray_intersects_triangle(ox, oy, oz, v1, v2, v3)

Möller–Trumbore ray-triangle intersection algorithm.
Tests if a ray from point (ox, oy, oz) in the +x direction intersects the triangle.
"""
function ray_intersects_triangle(ox::Float64, oy::Float64, oz::Float64,
                                  v1::Tuple{Float64, Float64, Float64},
                                  v2::Tuple{Float64, Float64, Float64},
                                  v3::Tuple{Float64, Float64, Float64})
    epsilon = 1e-10
    
    # Ray direction: +x
    ray_dir = (1.0, 0.0, 0.0)
    
    # Edge vectors
    e1 = (v2[1] - v1[1], v2[2] - v1[2], v2[3] - v1[3])
    e2 = (v3[1] - v1[1], v3[2] - v1[2], v3[3] - v1[3])
    
    # Cross product of ray direction and e2
    h = cross3(ray_dir, e2)
    
    # Dot product of e1 and h
    a = dot3(e1, h)
    
    # Ray is parallel to triangle
    if abs(a) < epsilon
        return false
    end
    
    f = 1.0 / a
    s = (ox - v1[1], oy - v1[2], oz - v1[3])
    u = f * dot3(s, h)
    
    if u < 0.0 || u > 1.0
        return false
    end
    
    q = cross3(s, e1)
    v = f * dot3(ray_dir, q)
    
    if v < 0.0 || u + v > 1.0
        return false
    end
    
    # Compute t to find intersection point
    t = f * dot3(e2, q)
    
    # Ray intersection (t > 0 means intersection is in front of ray origin)
    return t > epsilon
end

"""
Cross product for 3D tuples
"""
function cross3(a::Tuple{Float64, Float64, Float64}, b::Tuple{Float64, Float64, Float64})
    return (a[2]*b[3] - a[3]*b[2],
            a[3]*b[1] - a[1]*b[3],
            a[1]*b[2] - a[2]*b[1])
end

"""
Dot product for 3D tuples
"""
function dot3(a::Tuple{Float64, Float64, Float64}, b::Tuple{Float64, Float64, Float64})
    return a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
end

"""
    sdf(ft::FrontTracker3D, x::Float64, y::Float64, z::Float64)

Calculates the signed distance function for a given point.
Positive outside fluid, negative inside fluid.
"""
function sdf(ft::FrontTracker3D, x::Float64, y::Float64, z::Float64)
    if isempty(ft.markers) || isempty(ft.faces)
        return Inf
    end
    
    # Calculate minimum distance to any triangular face
    min_dist = Inf
    
    for (i1, i2, i3) in ft.faces
        v1 = ft.markers[i1]
        v2 = ft.markers[i2]
        v3 = ft.markers[i3]
        
        dist = point_triangle_distance((x, y, z), v1, v2, v3)
        min_dist = min(min_dist, dist)
    end
    
    # Determine sign (negative inside, positive outside)
    is_inside_val = is_point_inside(ft, x, y, z)
    
    return is_inside_val ? -min_dist : min_dist
end

"""
    point_triangle_distance(p, v1, v2, v3)

Calculates the distance from point p to triangle (v1, v2, v3).
"""
function point_triangle_distance(p::Tuple{Float64, Float64, Float64},
                                  v1::Tuple{Float64, Float64, Float64},
                                  v2::Tuple{Float64, Float64, Float64},
                                  v3::Tuple{Float64, Float64, Float64})
    # Edge vectors
    e0 = (v2[1] - v1[1], v2[2] - v1[2], v2[3] - v1[3])
    e1 = (v3[1] - v1[1], v3[2] - v1[2], v3[3] - v1[3])
    
    # Vector from v1 to p
    d = (v1[1] - p[1], v1[2] - p[2], v1[3] - p[3])
    
    a = dot3(e0, e0)
    b = dot3(e0, e1)
    c = dot3(e1, e1)
    d_dot_e0 = dot3(d, e0)
    d_dot_e1 = dot3(d, e1)
    
    det = a * c - b * b
    s = b * d_dot_e1 - c * d_dot_e0
    t = b * d_dot_e0 - a * d_dot_e1
    
    if det < 1e-10
        # Degenerate triangle
        return sqrt(dot3(d, d))
    end
    
    if s + t <= det
        if s < 0
            if t < 0
                # Region 4
                if d_dot_e0 < 0
                    s = clamp(-d_dot_e0 / a, 0.0, 1.0)
                    t = 0.0
                else
                    s = 0.0
                    t = clamp(-d_dot_e1 / c, 0.0, 1.0)
                end
            else
                # Region 3
                s = 0.0
                t = clamp(-d_dot_e1 / c, 0.0, 1.0)
            end
        elseif t < 0
            # Region 5
            t = 0.0
            s = clamp(-d_dot_e0 / a, 0.0, 1.0)
        else
            # Region 0
            inv_det = 1.0 / det
            s *= inv_det
            t *= inv_det
        end
    else
        if s < 0
            # Region 2
            tmp0 = b + d_dot_e0
            tmp1 = c + d_dot_e1
            if tmp1 > tmp0
                numer = tmp1 - tmp0
                denom = a - 2 * b + c
                s = clamp(numer / denom, 0.0, 1.0)
                t = 1.0 - s
            else
                s = 0.0
                t = clamp(-d_dot_e1 / c, 0.0, 1.0)
            end
        elseif t < 0
            # Region 6
            tmp0 = b + d_dot_e1
            tmp1 = a + d_dot_e0
            if tmp1 > tmp0
                numer = tmp1 - tmp0
                denom = a - 2 * b + c
                t = clamp(numer / denom, 0.0, 1.0)
                s = 1.0 - t
            else
                t = 0.0
                s = clamp(-d_dot_e0 / a, 0.0, 1.0)
            end
        else
            # Region 1
            numer = (c + d_dot_e1) - (b + d_dot_e0)
            if numer <= 0
                s = 0.0
            else
                denom = a - 2 * b + c
                s = clamp(numer / denom, 0.0, 1.0)
            end
            t = 1.0 - s
        end
    end
    
    # Closest point on triangle
    closest = (v1[1] + s * e0[1] + t * e1[1],
               v1[2] + s * e0[2] + t * e1[2],
               v1[3] + s * e0[3] + t * e1[3])
    
    diff = (p[1] - closest[1], p[2] - closest[2], p[3] - closest[3])
    return sqrt(dot3(diff, diff))
end

"""
    compute_face_normal(ft::FrontTracker3D, face_idx::Int)

Computes the outward normal vector for a triangular face.
"""
function compute_face_normal(ft::FrontTracker3D, face_idx::Int)
    if face_idx < 1 || face_idx > length(ft.faces)
        error("Invalid face index: $face_idx")
    end
    
    i1, i2, i3 = ft.faces[face_idx]
    v1 = ft.markers[i1]
    v2 = ft.markers[i2]
    v3 = ft.markers[i3]
    
    # Edge vectors
    e1 = (v2[1] - v1[1], v2[2] - v1[2], v2[3] - v1[3])
    e2 = (v3[1] - v1[1], v3[2] - v1[2], v3[3] - v1[3])
    
    # Cross product gives normal
    normal = cross3(e1, e2)
    
    # Normalize
    len = sqrt(dot3(normal, normal))
    if len > 1e-10
        normal = (normal[1] / len, normal[2] / len, normal[3] / len)
    else
        normal = (0.0, 0.0, 1.0)  # Default normal for degenerate triangles
    end
    
    return normal
end

"""
    compute_face_area(ft::FrontTracker3D, face_idx::Int)

Computes the area of a triangular face.
"""
function compute_face_area(ft::FrontTracker3D, face_idx::Int)
    if face_idx < 1 || face_idx > length(ft.faces)
        error("Invalid face index: $face_idx")
    end
    
    i1, i2, i3 = ft.faces[face_idx]
    v1 = ft.markers[i1]
    v2 = ft.markers[i2]
    v3 = ft.markers[i3]
    
    # Edge vectors
    e1 = (v2[1] - v1[1], v2[2] - v1[2], v2[3] - v1[3])
    e2 = (v3[1] - v1[1], v3[2] - v1[2], v3[3] - v1[3])
    
    # Cross product magnitude / 2 = area
    cross = cross3(e1, e2)
    return 0.5 * sqrt(dot3(cross, cross))
end

"""
    compute_total_surface_area(ft::FrontTracker3D)

Computes the total surface area of the 3D interface.
"""
function compute_total_surface_area(ft::FrontTracker3D)
    total_area = 0.0
    for i in 1:length(ft.faces)
        total_area += compute_face_area(ft, i)
    end
    return total_area
end

"""
    compute_enclosed_volume(ft::FrontTracker3D)

Computes the volume enclosed by the 3D interface using the divergence theorem.
"""
function compute_enclosed_volume(ft::FrontTracker3D)
    if !ft.is_closed || isempty(ft.faces)
        return 0.0
    end
    
    volume = 0.0
    for (i1, i2, i3) in ft.faces
        v1 = ft.markers[i1]
        v2 = ft.markers[i2]
        v3 = ft.markers[i3]
        
        # Signed volume of tetrahedron formed with origin
        # V = (1/6) * dot(v1, cross(v2, v3))
        cross_v2_v3 = cross3(v2, v3)
        volume += dot3(v1, cross_v2_v3)
    end
    
    return abs(volume) / 6.0
end

"""
    compute_centroid(ft::FrontTracker3D)

Computes the centroid of the volume enclosed by the 3D interface.
"""
function compute_centroid(ft::FrontTracker3D)
    if !ft.is_closed || isempty(ft.faces)
        # Return centroid of markers if not closed
        if isempty(ft.markers)
            return (0.0, 0.0, 0.0)
        end
        cx = sum(m[1] for m in ft.markers) / length(ft.markers)
        cy = sum(m[2] for m in ft.markers) / length(ft.markers)
        cz = sum(m[3] for m in ft.markers) / length(ft.markers)
        return (cx, cy, cz)
    end
    
    # For a closed surface, compute volume-weighted centroid
    total_volume = 0.0
    cx, cy, cz = 0.0, 0.0, 0.0
    
    for (i1, i2, i3) in ft.faces
        v1 = ft.markers[i1]
        v2 = ft.markers[i2]
        v3 = ft.markers[i3]
        
        # Signed volume of tetrahedron formed with origin
        cross_v2_v3 = cross3(v2, v3)
        vol = dot3(v1, cross_v2_v3) / 6.0
        
        # Centroid of tetrahedron is average of 4 vertices (origin + 3 triangle vertices)
        tet_cx = (v1[1] + v2[1] + v3[1]) / 4.0
        tet_cy = (v1[2] + v2[2] + v3[2]) / 4.0
        tet_cz = (v1[3] + v2[3] + v3[3]) / 4.0
        
        cx += vol * tet_cx
        cy += vol * tet_cy
        cz += vol * tet_cz
        total_volume += vol
    end
    
    if abs(total_volume) > 1e-10
        cx /= total_volume
        cy /= total_volume
        cz /= total_volume
    end
    
    return (cx, cy, cz)
end

"""
    compute_marker_normals(ft::FrontTracker3D)

Computes normal vectors at each marker by averaging adjacent face normals.
"""
function compute_marker_normals(ft::FrontTracker3D)
    n_markers = length(ft.markers)
    normals = Vector{Tuple{Float64, Float64, Float64}}(undef, n_markers)
    
    # Initialize normals to zero
    for i in 1:n_markers
        normals[i] = (0.0, 0.0, 0.0)
    end
    
    # Accumulate face normals weighted by area
    for (face_idx, (i1, i2, i3)) in enumerate(ft.faces)
        face_normal = compute_face_normal(ft, face_idx)
        face_area = compute_face_area(ft, face_idx)
        
        weighted_normal = (face_normal[1] * face_area, 
                          face_normal[2] * face_area, 
                          face_normal[3] * face_area)
        
        normals[i1] = (normals[i1][1] + weighted_normal[1],
                       normals[i1][2] + weighted_normal[2],
                       normals[i1][3] + weighted_normal[3])
        normals[i2] = (normals[i2][1] + weighted_normal[1],
                       normals[i2][2] + weighted_normal[2],
                       normals[i2][3] + weighted_normal[3])
        normals[i3] = (normals[i3][1] + weighted_normal[1],
                       normals[i3][2] + weighted_normal[2],
                       normals[i3][3] + weighted_normal[3])
    end
    
    # Normalize the accumulated normals
    for i in 1:n_markers
        len = sqrt(dot3(normals[i], normals[i]))
        if len > 1e-10
            normals[i] = (normals[i][1] / len, normals[i][2] / len, normals[i][3] / len)
        else
            normals[i] = (0.0, 0.0, 1.0)  # Default
        end
    end
    
    return normals
end

"""
    compute_fluid_volume_in_cell_3d(ft::FrontTracker3D, cell_bounds::Tuple)

Computes the fluid volume inside a single 3D cell using the divergence theorem.
cell_bounds = ((x_min, x_max), (y_min, y_max), (z_min, z_max))

Uses the divergence theorem: V = (1/3) ∫∫_∂Ω x·n dA
The boundary ∂Ω consists of:
1. Interface triangles clipped to the cell
2. Cell face portions that are inside the fluid

This is a geometric method that computes exact volumes based on edge/face intersections.
"""
function compute_fluid_volume_in_cell_3d(ft::FrontTracker3D, cell_bounds::Tuple)
    x_min, x_max = cell_bounds[1]
    y_min, y_max = cell_bounds[2]
    z_min, z_max = cell_bounds[3]
    
    if isempty(ft.faces) || isempty(ft.markers)
        return 0.0
    end
    
    cell_volume = (x_max - x_min) * (y_max - y_min) * (z_max - z_min)
    
    # Check if cell center is inside
    cx, cy, cz = (x_min + x_max) / 2, (y_min + y_max) / 2, (z_min + z_max) / 2
    center_inside = is_point_inside(ft, cx, cy, cz)
    
    # Check all 8 corners
    corners = [
        (x_min, y_min, z_min), (x_max, y_min, z_min),
        (x_min, y_max, z_min), (x_max, y_max, z_min),
        (x_min, y_min, z_max), (x_max, y_min, z_max),
        (x_min, y_max, z_max), (x_max, y_max, z_max)
    ]
    corners_inside = [is_point_inside(ft, c[1], c[2], c[3]) for c in corners]
    n_corners_inside = sum(corners_inside)
    
    # If all corners are inside, cell is fully inside fluid
    if n_corners_inside == 8
        return cell_volume
    end
    
    # If all corners are outside and center is outside, cell is fully outside
    if n_corners_inside == 0 && !center_inside
        # Double check by testing if any face intersects
        has_intersection = false
        for (i1, i2, i3) in ft.faces
            v1 = ft.markers[i1]
            v2 = ft.markers[i2]
            v3 = ft.markers[i3]
            if triangle_intersects_box(v1, v2, v3, cell_bounds)
                has_intersection = true
                break
            end
        end
        if !has_intersection
            return 0.0
        end
    end
    
    # For cut cells, use divergence theorem: V = (1/3) ∫∫ x·n dA
    # Sum contributions from interface triangles clipped to cell
    volume = 0.0
    
    for (i1, i2, i3) in ft.faces
        v1 = ft.markers[i1]
        v2 = ft.markers[i2]
        v3 = ft.markers[i3]
        
        # Clip triangle to cell bounds
        clipped_polygon = clip_triangle_to_box(v1, v2, v3, cell_bounds)
        
        if length(clipped_polygon) >= 3
            # Compute contribution using divergence theorem
            # For each triangle in the clipped polygon (fan triangulation)
            p0 = clipped_polygon[1]
            for i in 2:length(clipped_polygon)-1
                p1 = clipped_polygon[i]
                p2 = clipped_polygon[i+1]
                
                # Triangle contribution: (1/3) * centroid · normal * area
                # Using signed volume of tetrahedron with origin
                vol_contrib = signed_tetrahedron_volume((0.0, 0.0, 0.0), p0, p1, p2)
                volume += vol_contrib
            end
        end
    end
    
    # Add contributions from cell faces that are inside the fluid
    # For each cell face, compute the area inside the fluid
    volume += compute_cell_face_contributions(ft, cell_bounds, corners_inside)
    
    return abs(volume)
end

"""
    triangle_intersects_box(v1, v2, v3, cell_bounds)

Quick check if a triangle might intersect an axis-aligned box.
"""
function triangle_intersects_box(v1::Tuple{Float64, Float64, Float64},
                                  v2::Tuple{Float64, Float64, Float64},
                                  v3::Tuple{Float64, Float64, Float64},
                                  cell_bounds::Tuple)
    x_min, x_max = cell_bounds[1]
    y_min, y_max = cell_bounds[2]
    z_min, z_max = cell_bounds[3]
    
    # Check bounding box overlap
    tri_x_min = min(v1[1], v2[1], v3[1])
    tri_x_max = max(v1[1], v2[1], v3[1])
    tri_y_min = min(v1[2], v2[2], v3[2])
    tri_y_max = max(v1[2], v2[2], v3[2])
    tri_z_min = min(v1[3], v2[3], v3[3])
    tri_z_max = max(v1[3], v2[3], v3[3])
    
    return !(tri_x_max < x_min || tri_x_min > x_max ||
             tri_y_max < y_min || tri_y_min > y_max ||
             tri_z_max < z_min || tri_z_min > z_max)
end

"""
    clip_triangle_to_box(v1, v2, v3, cell_bounds)

Clips a triangle to an axis-aligned box using Sutherland-Hodgman algorithm.
Returns a polygon (list of vertices) representing the clipped region.
"""
function clip_triangle_to_box(v1::Tuple{Float64, Float64, Float64},
                               v2::Tuple{Float64, Float64, Float64},
                               v3::Tuple{Float64, Float64, Float64},
                               cell_bounds::Tuple)
    x_min, x_max = cell_bounds[1]
    y_min, y_max = cell_bounds[2]
    z_min, z_max = cell_bounds[3]
    
    # Start with triangle vertices
    polygon = [v1, v2, v3]
    
    # Clip against each of the 6 planes
    polygon = clip_polygon_against_plane(polygon, (1.0, 0.0, 0.0), x_min)
    if isempty(polygon) return polygon end
    
    polygon = clip_polygon_against_plane(polygon, (-1.0, 0.0, 0.0), -x_max)
    if isempty(polygon) return polygon end
    
    polygon = clip_polygon_against_plane(polygon, (0.0, 1.0, 0.0), y_min)
    if isempty(polygon) return polygon end
    
    polygon = clip_polygon_against_plane(polygon, (0.0, -1.0, 0.0), -y_max)
    if isempty(polygon) return polygon end
    
    polygon = clip_polygon_against_plane(polygon, (0.0, 0.0, 1.0), z_min)
    if isempty(polygon) return polygon end
    
    polygon = clip_polygon_against_plane(polygon, (0.0, 0.0, -1.0), -z_max)
    
    return polygon
end

"""
    clip_polygon_against_plane(polygon, normal, d)

Clips a polygon against a half-space defined by normal · p >= d.
"""
function clip_polygon_against_plane(polygon::Vector{Tuple{Float64, Float64, Float64}},
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
        
        current_dist = dot3(normal, current) - d
        next_dist = dot3(normal, next) - d
        
        current_inside = current_dist >= -1e-10
        next_inside = next_dist >= -1e-10
        
        if current_inside
            push!(output, current)
            if !next_inside
                # Compute intersection
                t = current_dist / (current_dist - next_dist)
                t = clamp(t, 0.0, 1.0)
                intersection = (
                    current[1] + t * (next[1] - current[1]),
                    current[2] + t * (next[2] - current[2]),
                    current[3] + t * (next[3] - current[3])
                )
                push!(output, intersection)
            end
        elseif next_inside
            # Compute intersection
            t = current_dist / (current_dist - next_dist)
            t = clamp(t, 0.0, 1.0)
            intersection = (
                current[1] + t * (next[1] - current[1]),
                current[2] + t * (next[2] - current[2]),
                current[3] + t * (next[3] - current[3])
            )
            push!(output, intersection)
        end
    end
    
    return output
end

"""
    signed_tetrahedron_volume(p0, p1, p2, p3)

Computes the signed volume of a tetrahedron with vertices p0, p1, p2, p3.
Volume = (1/6) * |det([p1-p0, p2-p0, p3-p0])|
"""
function signed_tetrahedron_volume(p0::Tuple{Float64, Float64, Float64},
                                    p1::Tuple{Float64, Float64, Float64},
                                    p2::Tuple{Float64, Float64, Float64},
                                    p3::Tuple{Float64, Float64, Float64})
    # Vectors from p0 to other points
    a = (p1[1] - p0[1], p1[2] - p0[2], p1[3] - p0[3])
    b = (p2[1] - p0[1], p2[2] - p0[2], p2[3] - p0[3])
    c = (p3[1] - p0[1], p3[2] - p0[2], p3[3] - p0[3])
    
    # Determinant = a · (b × c)
    cross_bc = cross3(b, c)
    det = dot3(a, cross_bc)
    
    return det / 6.0
end

"""
    compute_cell_face_contributions(ft::FrontTracker3D, cell_bounds::Tuple, corners_inside::Vector{Bool})

Computes the volume contributions from cell faces using the divergence theorem.
For each cell face, we need to account for the portion that bounds the fluid region.
"""
function compute_cell_face_contributions(ft::FrontTracker3D, cell_bounds::Tuple, corners_inside::Vector{Bool})
    x_min, x_max = cell_bounds[1]
    y_min, y_max = cell_bounds[2]
    z_min, z_max = cell_bounds[3]
    
    volume = 0.0
    
    # For each face of the cell, if it's entirely inside the fluid or partially inside,
    # we need to add its contribution to the volume integral
    
    # Face at x = x_min (normal pointing -x), corners 1,3,5,7
    face_corners = [corners_inside[1], corners_inside[3], corners_inside[7], corners_inside[5]]
    if all(face_corners)
        # Entire face inside fluid: contribution = x_min * face_area * (-1) / 3 = -x_min * dy * dz / 3
        # Using divergence theorem with n = (-1, 0, 0)
        dy = y_max - y_min
        dz = z_max - z_min
        volume -= x_min * dy * dz / 3.0
    end
    
    # Face at x = x_max (normal pointing +x), corners 2,4,6,8
    face_corners = [corners_inside[2], corners_inside[4], corners_inside[8], corners_inside[6]]
    if all(face_corners)
        dy = y_max - y_min
        dz = z_max - z_min
        volume += x_max * dy * dz / 3.0
    end
    
    # Face at y = y_min (normal pointing -y), corners 1,2,5,6
    face_corners = [corners_inside[1], corners_inside[2], corners_inside[6], corners_inside[5]]
    if all(face_corners)
        dx = x_max - x_min
        dz = z_max - z_min
        volume -= y_min * dx * dz / 3.0
    end
    
    # Face at y = y_max (normal pointing +y), corners 3,4,7,8
    face_corners = [corners_inside[3], corners_inside[4], corners_inside[8], corners_inside[7]]
    if all(face_corners)
        dx = x_max - x_min
        dz = z_max - z_min
        volume += y_max * dx * dz / 3.0
    end
    
    # Face at z = z_min (normal pointing -z), corners 1,2,3,4
    face_corners = [corners_inside[1], corners_inside[2], corners_inside[4], corners_inside[3]]
    if all(face_corners)
        dx = x_max - x_min
        dy = y_max - y_min
        volume -= z_min * dx * dy / 3.0
    end
    
    # Face at z = z_max (normal pointing +z), corners 5,6,7,8
    face_corners = [corners_inside[5], corners_inside[6], corners_inside[8], corners_inside[7]]
    if all(face_corners)
        dx = x_max - x_min
        dy = y_max - y_min
        volume += z_max * dx * dy / 3.0
    end
    
    return volume
end

"""
    compute_volume_jacobian_3d(ft::FrontTracker3D, x_faces::AbstractVector{<:Real}, 
                               y_faces::AbstractVector{<:Real}, z_faces::AbstractVector{<:Real}, 
                               epsilon::Float64=1e-3)

Calculates the volume Jacobian matrix for a given 3D mesh and interface.
The Jacobian represents ∂V_cell/∂δ_marker, the sensitivity of cell fluid volumes
to marker displacements along their normal directions.

Uses central differencing: J[i,j,k,m] = (V⁺ - V⁻) / (2ε)

The default epsilon is 1e-3, which should be large enough for the Gauss-Legendre 
quadrature to detect changes in volume from boundary perturbations.

Returns a dictionary mapping cell indices (i,j,k) to lists of (marker_idx, jacobian_value).
"""
function compute_volume_jacobian_3d(ft::FrontTracker3D, x_faces::AbstractVector{<:Real}, 
                                    y_faces::AbstractVector{<:Real}, z_faces::AbstractVector{<:Real}, 
                                    epsilon::Float64=1e-3)
    # Convert to vectors if needed
    x_faces_vec = collect(x_faces)
    y_faces_vec = collect(y_faces)
    z_faces_vec = collect(z_faces)
    
    # Get mesh dimensions
    nx = length(x_faces_vec) - 1
    ny = length(y_faces_vec) - 1
    nz = length(z_faces_vec) - 1
    
    # Get markers and compute their normals
    markers = get_markers(ft)
    normals = compute_marker_normals(ft)
    n_markers = length(markers)
    
    # Calculate original cell volumes
    original_volumes = Dict{Tuple{Int, Int, Int}, Float64}()
    
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                cell_bounds = (
                    (x_faces_vec[i], x_faces_vec[i+1]),
                    (y_faces_vec[j], y_faces_vec[j+1]),
                    (z_faces_vec[k], z_faces_vec[k+1])
                )
                original_volumes[(i, j, k)] = compute_fluid_volume_in_cell_3d(ft, cell_bounds)
            end
        end
    end
    
    # Initialize dictionary for storing the Jacobian
    volume_jacobian = Dict{Tuple{Int, Int, Int}, Vector{Tuple{Int, Float64}}}()
    for key in keys(original_volumes)
        volume_jacobian[key] = []
    end
    
    # Track which markers have entries
    markers_with_entries = Set{Int}()
    
    for marker_idx in 1:n_markers
        # Original marker position
        original_marker = markers[marker_idx]
        normal = normals[marker_idx]
        
        # Positive perturbation
        pos_perturbed_marker = (
            original_marker[1] + epsilon * normal[1],
            original_marker[2] + epsilon * normal[2],
            original_marker[3] + epsilon * normal[3]
        )
        
        # Negative perturbation
        neg_perturbed_marker = (
            original_marker[1] - epsilon * normal[1],
            original_marker[2] - epsilon * normal[2],
            original_marker[3] - epsilon * normal[3]
        )
        
        # Create copies of markers with perturbations
        pos_perturbed_markers = copy(markers)
        pos_perturbed_markers[marker_idx] = pos_perturbed_marker
        
        neg_perturbed_markers = copy(markers)
        neg_perturbed_markers[marker_idx] = neg_perturbed_marker
        
        # Create new front trackers with perturbed markers
        pos_perturbed_tracker = FrontTracker3D(pos_perturbed_markers, ft.faces, ft.is_closed)
        neg_perturbed_tracker = FrontTracker3D(neg_perturbed_markers, ft.faces, ft.is_closed)
        
        # Track the max Jacobian value for this marker
        max_jac_value = 0.0
        max_jac_cell = nothing
        
        # Calculate perturbed volumes using central differencing
        for ((i, j, k), _) in original_volumes
            cell_bounds = (
                (x_faces_vec[i], x_faces_vec[i+1]),
                (y_faces_vec[j], y_faces_vec[j+1]),
                (z_faces_vec[k], z_faces_vec[k+1])
            )
            
            # Calculate volumes with perturbed fronts
            pos_volume = compute_fluid_volume_in_cell_3d(pos_perturbed_tracker, cell_bounds)
            neg_volume = compute_fluid_volume_in_cell_3d(neg_perturbed_tracker, cell_bounds)
            
            # Calculate Jacobian value using central differencing
            jacobian_value = (pos_volume - neg_volume) / (2.0 * epsilon)
            
            # Store significant changes
            if abs(jacobian_value) > 1e-10
                push!(volume_jacobian[(i, j, k)], (marker_idx, jacobian_value))
                push!(markers_with_entries, marker_idx)
                
                if abs(jacobian_value) > abs(max_jac_value)
                    max_jac_value = jacobian_value
                    max_jac_cell = (i, j, k)
                end
            elseif abs(jacobian_value) > abs(max_jac_value)
                max_jac_value = jacobian_value
                max_jac_cell = (i, j, k)
            end
        end
        
        # If this marker has no entries, add its maximum value entry
        if marker_idx ∉ markers_with_entries && max_jac_cell !== nothing
            push!(volume_jacobian[max_jac_cell], (marker_idx, max_jac_value))
            push!(markers_with_entries, marker_idx)
        end
    end
    
    return volume_jacobian
end
