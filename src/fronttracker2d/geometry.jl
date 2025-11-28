"""
    get_fluid_polygon(ft::FrontTracker)

Returns a polygon representing the fluid domain bounded by the interface.
"""
function get_fluid_polygon(ft::FrontTracker)
    if length(ft.markers) < 3
        return nothing
    end
    
    if ft.is_closed
        # For a closed interface, directly return the polygon
        if !LibGEOS.isValid(ft.interface_poly)
            # Try to fix invalid polygon
            return LibGEOS.buffer(ft.interface_poly, 0.0)
        end
        return ft.interface_poly
    else
        # For an open interface, we need to define domain closure
        # Using convex hull as an approximation
        return LibGEOS.convexhull(ft.interface_poly)
    end
end

"""
    is_point_inside(ft::FrontTracker, x::Float64, y::Float64)

Checks if a point is inside the interface.
"""
function is_point_inside(ft::FrontTracker, x::Float64, y::Float64)
    if isnothing(ft.interface_poly)
        return false
    end
    
    point = LibGEOS.Point(x, y)
    return LibGEOS.contains(ft.interface_poly, point)
end

"""
    get_intersection(ft::FrontTracker, other_geometry)

Calculates the intersection with another geometry.
"""
function get_intersection(ft::FrontTracker, other_geometry)
    if isnothing(ft.interface_poly)
        return nothing
    end
    
    return LibGEOS.intersection(other_geometry, ft.interface_poly)
end


"""
    sdf(ft::FrontTracker, x::Float64, y::Float64)

Calculates the signed distance function for a given point.
"""
function sdf(ft::FrontTracker, x, y)
    if isnothing(ft.interface)
        return Inf
    end
    
    # Create a point using LibGEOS
    point = LibGEOS.Point(x, y)
    
    # Calculate the distance to the interface
    distance = LibGEOS.distance(point, ft.interface)
    
    # Determine the sign (negative inside, positive outside)
    is_inside_val = is_point_inside(ft, x, y)
    
    return is_inside_val ? -distance : distance
end

"""
    compute_marker_normals(ft::FrontTracker, markers=nothing)

Calculates the normal vectors for each marker of the interface.
"""
function compute_marker_normals(ft::FrontTracker, markers=nothing)
    if isnothing(markers)
        markers = ft.markers
    end
    
    if length(markers) < 3
        return [(0.0, 1.0) for _ in markers]
    end
    
    normals = []
    n_markers = length(markers)
    is_closed = ft.is_closed
    
    # For maintaining orientation consistency
    prev_normal = nothing
    
    for i in 1:n_markers
        # Handle indices for previous and next points with boundary conditions
        prev_idx = is_closed ? mod1(i-1, n_markers) : max(1, i-1)
        next_idx = is_closed ? mod1(i+1, n_markers) : min(n_markers, i+1)
        
        # MODIFICATION IMPORTANTE: Utiliser la méthode des tangentes pour tous les points du premier marqueur
        # Pour les contours fermés, traiter le premier marqueur comme un point spécial
        if is_closed && (i == 1 || i == n_markers && markers[1] == markers[end])
            # Utiliser les points adjacents pour calculer les tangentes
            prev_idx = n_markers > 2 ? n_markers - 1 : 1
            next_idx = 2
            
            t1_x = markers[i][1] - markers[prev_idx][1]
            t1_y = markers[i][2] - markers[prev_idx][2]
            
            t2_x = markers[next_idx][1] - markers[i][1]
            t2_y = markers[next_idx][2] - markers[i][2]
            
            # Normaliser les vecteurs tangents
            t1_len = sqrt(t1_x^2 + t1_y^2)
            t2_len = sqrt(t2_x^2 + t2_y^2)
            
            if t1_len > 0 && t2_len > 0
                t1_x, t1_y = t1_x/t1_len, t1_y/t1_len
                t2_x, t2_y = t2_x/t2_len, t2_y/t2_len
                
                # Moyenner les tangentes
                tx = (t1_x + t2_x)
                ty = (t1_y + t2_y)
                
                # Normaliser
                t_len = sqrt(tx^2 + ty^2)
                if t_len > 0
                    tx, ty = tx/t_len, ty/t_len
                    
                    # La normale est perpendiculaire à la tangente
                    n_x, n_y = -ty, tx
                    
                    # Vérifier l'orientation
                    test_x = markers[i][1] + 1e-3 * n_x
                    test_y = markers[i][2] + 1e-3 * n_y
                    if is_point_inside(ft, test_x, test_y)
                        n_x, n_y = -n_x, -n_y
                    end
                    
                    push!(normals, (n_x, n_y))
                    prev_normal = [n_x, n_y]
                    
                    # Si c'est le dernier point qui est dupliqué, assurez-vous de sauter l'itération
                    if i == n_markers && markers[1] == markers[end]
                        continue
                    end
                else
                    push!(normals, (0.0, 1.0))
                    prev_normal = [0.0, 1.0]
                end
                
                continue
            end
        end
        # Special handling for endpoints of open curves
        if !is_closed && (i == 1 || i == n_markers)
            # For first point in open curve, use forward difference
            if i == 1
                # Use the vector from first to second marker as tangent
                t_x = markers[2][1] - markers[1][1]
                t_y = markers[2][2] - markers[1][2]
                
                # Normalize tangent
                t_len = sqrt(t_x^2 + t_y^2)
                if t_len > 0
                    t_x /= t_len
                    t_y /= t_len
                    
                    # Rotate 90° to get normal (outward pointing)
                    n_x, n_y = -t_y, t_x
                else
                    n_x, n_y = 0.0, 1.0
                end
                
                push!(normals, (n_x, n_y))
                prev_normal = [n_x, n_y]
                continue
            end
            
            # For last point in open curve, use backward difference
            if i == n_markers
                # Use the vector from second-to-last to last marker as tangent
                t_x = markers[n_markers][1] - markers[n_markers-1][1]
                t_y = markers[n_markers][2] - markers[n_markers-1][2]
                
                # Normalize tangent
                t_len = sqrt(t_x^2 + t_y^2)
                if t_len > 0
                    t_x /= t_len
                    t_y /= t_len
                    
                    # Rotate 90° to get normal (outward pointing)
                    n_x, n_y = -t_y, t_x
                    
                    # Ensure consistency with previous normal
                    if !isnothing(prev_normal)
                        dot_product = n_x * prev_normal[1] + n_y * prev_normal[2]
                        if dot_product < 0
                            n_x, n_y = -n_x, -n_y
                        end
                    end
                else
                    n_x, n_y = prev_normal[1], prev_normal[2]
                end
                
                push!(normals, (n_x, n_y))
                continue
            end
        end
        
        # Regular points use the osculating circle method
        # Points P1, P2, P3 for calculating the osculating circle
        p1 = collect(markers[prev_idx])
        p2 = collect(markers[i])
        p3 = collect(markers[next_idx])
        
        # Rest of the function remains the same...
        
        try
            # Check if points are distinct
            if (norm(p1 - p2) < 1e-10 || 
                norm(p2 - p3) < 1e-10 || 
                norm(p3 - p1) < 1e-10)
                error("Points too close")
            end
            
            # Equations to find the center of the circle through three points
            # Linear system: Ax = b where x = [center_x, center_y, -r²]
            A = [2*p1[1] 2*p1[2] 1;
                 2*p2[1] 2*p2[2] 1;
                 2*p3[1] 2*p3[2] 1]
            
            b = [p1[1]^2 + p1[2]^2,
                 p2[1]^2 + p2[2]^2,
                 p3[1]^2 + p3[2]^2]
            
            # Solve for the circle center
            x = A \ b
            center = x[1:2]
            
            # Calculate the normal vector (from point to center)
            normal_vector = center - p2
            
            # Normalize the vector
            norm_val = norm(normal_vector)
            if norm_val < 1e-10
                error("Norm too small")
            end
            normal_vector = normal_vector / norm_val
            
            # Check orientation
            if is_closed
                # For closed interface, check if normal points outward
                test_point = p2 + 1e-3 * normal_vector
                if is_point_inside(ft, test_point[1], test_point[2])
                    normal_vector = -normal_vector  # Invert normal
                end
            else
                # For open interface, ensure consistency with previous normal
                if !isnothing(prev_normal) && i > 1
                    # Calculate dot product to check if normals point in similar directions
                    dot_product = normal_vector[1] * prev_normal[1] + normal_vector[2] * prev_normal[2]
                    if dot_product < 0
                        normal_vector = -normal_vector
                    end
                elseif i == 1
                    # For first point, prefer upward-pointing normal
                    if normal_vector[2] < 0
                        normal_vector = -normal_vector
                    end
                end
            end
            
            prev_normal = normal_vector
            push!(normals, (normal_vector[1], normal_vector[2]))
            
        catch e
            # Fallback: use tangent method if osculating circle fails
            if i > 1 || is_closed
                t1_x = markers[i][1] - markers[prev_idx][1]
                t1_y = markers[i][2] - markers[prev_idx][2]
            else
                t1_x = markers[next_idx][1] - markers[i][1]
                t1_y = markers[next_idx][2] - markers[i][2]
            end
            
            if i < n_markers || is_closed
                t2_x = markers[next_idx][1] - markers[i][1]
                t2_y = markers[next_idx][2] - markers[i][2]
            else
                t2_x = markers[i][1] - markers[prev_idx][1]
                t2_y = markers[i][2] - markers[prev_idx][2]
            end
            
            # Normalize tangent vectors
            t1_len = sqrt(t1_x^2 + t1_y^2)
            t2_len = sqrt(t2_x^2 + t2_y^2)
            
            if t1_len > 0 && t2_len > 0
                t1_x, t1_y = t1_x/t1_len, t1_y/t1_len
                t2_x, t2_y = t2_x/t2_len, t2_y/t2_len
                
                # Average tangents
                tx = (t1_x + t2_x)
                ty = (t1_y + t2_y)
                
                # Normalize
                t_len = sqrt(tx^2 + ty^2)
                if t_len > 0
                    tx, ty = tx/t_len, ty/t_len
                    
                    # Normal is perpendicular to tangent
                    nx, ny = -ty, tx
                    
                    # Check orientation
                    if is_closed
                        test_x = markers[i][1] + 1e-3 * nx
                        test_y = markers[i][2] + 1e-3 * ny
                        if is_point_inside(ft, test_x, test_y)
                            nx, ny = -nx, -ny
                        end
                    else
                        # Ensure consistency
                        normal_vector = [nx, ny]
                        if !isnothing(prev_normal) && i > 1
                            dot_product = normal_vector[1] * prev_normal[1] + normal_vector[2] * prev_normal[2]
                            if dot_product < 0
                                nx, ny = -nx, -ny
                            end
                        elseif i == 1
                            if ny < 0
                                nx, ny = -nx, -ny
                            end
                        end
                    end
                    
                    prev_normal = [nx, ny]
                    push!(normals, (nx, ny))
                else
                    # Degenerate case
                    if !isnothing(prev_normal)
                        push!(normals, (prev_normal[1], prev_normal[2]))
                    else
                        push!(normals, (0.0, 1.0))
                    end
                end
            else
                # Degenerate case
                if !isnothing(prev_normal)
                    push!(normals, (prev_normal[1], prev_normal[2]))
                else
                    push!(normals, (0.0, 1.0))
                end
            end
        end
    end
    
    return normals
end
