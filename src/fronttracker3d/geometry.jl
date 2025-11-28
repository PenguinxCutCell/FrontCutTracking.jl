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
