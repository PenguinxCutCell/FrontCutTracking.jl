function compute_spacetime_capacities(nodes::NTuple{2, AbstractVector}, front_n::FrontTracker, front_np1::FrontTracker, dt::Float64)
    # Extract mesh information
    x_nodes = nodes[1]
    y_nodes = nodes[2]
    nx = length(x_nodes) - 1
    ny = length(y_nodes) - 1
    
    # Initialize capacity arrays
    Ax_spacetime = zeros(nx+1, ny+1)
    Ay_spacetime = zeros(nx+1, ny+1)
    V_spacetime = zeros(nx+1, ny+1)  # Space-time volumes for each cell

    # Initialize space-time centroid arrays
    centroid_x_spacetime = zeros(nx+1, ny+1)
    centroid_y_spacetime = zeros(nx+1, ny+1)
    
    # Initialize additional arrays in the main function
    Bx_spacetime = zeros(nx+1, ny+1)
    By_spacetime = zeros(nx+1, ny+1)
    Wx_spacetime = zeros(nx+1, ny+1)
    Wy_spacetime = zeros(nx+1, ny+1)

    # Arrays to store crossing times
    t_crosses_x = zeros(nx+1, ny)
    t_crosses_y = zeros(nx, ny+1)
    t_crosses_v = zeros(nx, ny)  # For volume crossing times

    # Pre-compute fluid cell properties at initial and final times
    _, volumes_n, _, _, _ = fluid_cell_properties(nodes, front_n)
    _, volumes_np1, _, _, _ = fluid_cell_properties(nodes, front_np1)
    
    # 1. Calculate Ax capacities (vertical faces)
    for i in 1:nx+1
        for j in 1:ny
            x_face = x_nodes[i]
            y_min, y_max = y_nodes[j], y_nodes[j+1]
            face_height = y_max - y_min
            
            # Check if vertices are inside fluid at time n and n+1
            bottom_n = is_point_inside(front_n, x_face, y_min)
            bottom_np1 = is_point_inside(front_np1, x_face, y_min)
            top_n = is_point_inside(front_n, x_face, y_max)
            top_np1 = is_point_inside(front_np1, x_face, y_max)
            
            # Find all crossing times
            crossings = Float64[]
            push!(crossings, 0.0)  # Add t_n
            push!(crossings, dt)   # Add t_{n+1}
            
            # Add bottom and top crossing times if the state changes
            if bottom_n != bottom_np1
                t_cross_bottom = find_crossing_time(front_n, front_np1, x_face, y_min, dt)
                push!(crossings, t_cross_bottom)
            end
            if top_n != top_np1
                t_cross_top = find_crossing_time(front_n, front_np1, x_face, y_max, dt)
                push!(crossings, t_cross_top)
            end
            
            # Sort crossing times
            sort!(crossings)
            
            # Integrate capacity over each subinterval using trapezoidal rule
            for k in 1:(length(crossings)-1)
                t_start = crossings[k]
                t_end = crossings[k+1]
                
                # Calculate τ values (normalized time)
                τ_start = t_start / dt
                τ_end = t_end / dt
                
                # Create intermediate fronts at t_start and t_end by linear interpolation
                front_start = interpolate_front(front_n, front_np1, τ_start)
                front_end = interpolate_front(front_n, front_np1, τ_end)
                
                # Calculate face capacity at t_start and t_end
                Ax_start = calculate_face_capacity_x(front_start, x_face, y_min, y_max)
                Ax_end = calculate_face_capacity_x(front_end, x_face, y_min, y_max)
                
                # Trapezoidal rule integration: ∫f(t)dt ≈ (t_end - t_start) * (f(t_start) + f(t_end)) / 2
                subinterval = (t_end - t_start) * (Ax_start + Ax_end) / 2
                
                # Add contribution to total capacity
                Ax_spacetime[i, j] += subinterval
            end
            
            # Store representative crossing time for visualization
            if bottom_n != bottom_np1
                t_crosses_x[i, j] = find_crossing_time(front_n, front_np1, x_face, y_min, dt)
            elseif top_n != top_np1
                t_crosses_x[i, j] = find_crossing_time(front_n, front_np1, x_face, y_max, dt)
            else
                t_crosses_x[i, j] = 0.5 * dt  # Default if no crossing
            end
            
            # Safety bounds check
            Ax_spacetime[i, j] = clamp(Ax_spacetime[i, j], 0.0, face_height * dt)
        end
    end
    
    # 2. Calculate Ay capacities (horizontal faces) - similar approach
    for i in 1:nx
        for j in 1:ny+1
            y_face = y_nodes[j]
            x_min, x_max = x_nodes[i], x_nodes[i+1]
            face_width = x_max - x_min
            
            # Check if vertices are inside fluid at time n and n+1
            left_n = is_point_inside(front_n, x_min, y_face)
            left_np1 = is_point_inside(front_np1, x_min, y_face)
            right_n = is_point_inside(front_n, x_max, y_face)
            right_np1 = is_point_inside(front_np1, x_max, y_face)
            
            # Find all crossing times
            crossings = Float64[]
            push!(crossings, 0.0)  # Add t_n
            push!(crossings, dt)   # Add t_{n+1}
            
            # Add left and right crossing times if the state changes
            if left_n != left_np1
                t_cross_left = find_crossing_time(front_n, front_np1, x_min, y_face, dt)
                push!(crossings, t_cross_left)
            end
            if right_n != right_np1
                t_cross_right = find_crossing_time(front_n, front_np1, x_max, y_face, dt)
                push!(crossings, t_cross_right)
            end
            
            # Sort crossing times
            sort!(crossings)
            
            # Integrate capacity over each subinterval using trapezoidal rule
            for k in 1:(length(crossings)-1)
                t_start = crossings[k]
                t_end = crossings[k+1]
                
                # Calculate τ values
                τ_start = t_start / dt
                τ_end = t_end / dt
                
                # Create intermediate fronts
                front_start = interpolate_front(front_n, front_np1, τ_start)
                front_end = interpolate_front(front_n, front_np1, τ_end)
                
                # Calculate face capacity
                Ay_start = calculate_face_capacity_y(front_start, y_face, x_min, x_max)
                Ay_end = calculate_face_capacity_y(front_end, y_face, x_min, x_max)
                
                # Trapezoidal rule integration
                subinterval = (t_end - t_start) * (Ay_start + Ay_end) / 2
                
                # Add contribution
                Ay_spacetime[i, j] += subinterval
            end
            
            # Store representative crossing time
            if left_n != left_np1
                t_crosses_y[i, j] = find_crossing_time(front_n, front_np1, x_min, y_face, dt)
            elseif right_n != right_np1
                t_crosses_y[i, j] = find_crossing_time(front_n, front_np1, x_max, y_face, dt)
            else
                t_crosses_y[i, j] = 0.5 * dt
            end
            
            # Safety bounds check
            Ay_spacetime[i, j] = clamp(Ay_spacetime[i, j], 0.0, face_width * dt)
        end
    end

    # 3. Calculate V (volume) space-time capacities
    for i in 1:nx
        for j in 1:ny
            # Get cell dimensions
            cell_width = x_nodes[i+1] - x_nodes[i]
            cell_height = y_nodes[j+1] - y_nodes[j]
            cell_area = cell_width * cell_height
            
            # Initial and final volumes
            vol_n = volumes_n[i, j]
            vol_np1 = volumes_np1[i, j]
            
            """
            # If volume didn't change much, use simple trapezoidal rule
            if abs(vol_n - vol_np1) < 1e-10 * cell_area
                V_spacetime[i, j] = 0.5 * (vol_n + vol_np1) * dt
                t_crosses_v[i, j] = 0.5 * dt
                continue
            end
            """

            # Get cell vertices
            vertices = [
                (x_nodes[i], y_nodes[j]),        # Bottom-left
                (x_nodes[i+1], y_nodes[j]),      # Bottom-right
                (x_nodes[i+1], y_nodes[j+1]),    # Top-right
                (x_nodes[i], y_nodes[j+1])       # Top-left
            ]
            
            # Get cell center
            cell_center_x = (x_nodes[i] + x_nodes[i+1]) / 2
            cell_center_y = (y_nodes[j] + y_nodes[j+1]) / 2
            
            # Collect ALL relevant crossing times from the face crossings we already computed
            crossings = Float64[]
            push!(crossings, 0.0)  # Add t_n
            push!(crossings, dt)   # Add t_{n+1}
            
            # Add crossing times from x-faces (left and right face of this cell)
            if t_crosses_x[i, j] > 0
                push!(crossings, t_crosses_x[i, j])
            end
            if t_crosses_x[i+1, j] > 0
                push!(crossings, t_crosses_x[i+1, j])
            end
            
            # Add crossing times from y-faces (bottom and top face of this cell)
            if t_crosses_y[i, j] > 0
                push!(crossings, t_crosses_y[i, j])
            end
            if t_crosses_y[i, j+1] > 0
                push!(crossings, t_crosses_y[i, j+1])
            end
            
            """
            # Check if cell center changes state
            center_n = is_point_inside(front_n, cell_center_x, cell_center_y)
            center_np1 = is_point_inside(front_np1, cell_center_x, cell_center_y)
            
            if center_n != center_np1
                t_cross_center = find_crossing_time(front_n, front_np1, cell_center_x, cell_center_y, dt)
                push!(crossings, t_cross_center)
                t_crosses_v[i, j] = t_cross_center
            else
                t_crosses_v[i, j] = center_n ? 0.0 : dt
            end
            """
            
            # Sort and remove duplicates
            sort!(unique!(crossings))
            
            # For cells with significant volume change, add additional midpoints for better sampling
            if length(crossings) > 2 && abs(vol_np1 - vol_n) > 0.01 * cell_area
                refined_crossings = Float64[]
                for k in 1:(length(crossings)-1)
                    push!(refined_crossings, crossings[k])
                    # Add midpoint between each pair of crossing times
                    midpoint = 0.5 * (crossings[k] + crossings[k+1])
                    push!(refined_crossings, midpoint)
                end
                push!(refined_crossings, crossings[end])
                crossings = sort(refined_crossings)
            end
            
            # Initialize space-time volume for this cell
            V_spacetime[i, j] = 0.0
            weighted_centroid_x = 0.0
            weighted_centroid_y = 0.0
            
            # Perform quadrature integration with improved accuracy
            if length(crossings) >= 3  # Simpson's rule requires at least 3 points
                for k in 1:2:length(crossings)-2
                    t_start = crossings[k]
                    t_mid = crossings[k+1]
                    t_end = crossings[k+2]
                    
                    # Calculate τ values (normalized time)
                    τ_start = t_start / dt
                    τ_mid = t_mid / dt
                    τ_end = t_end / dt
                    
                    # Create intermediate fronts
                    front_start = interpolate_front(front_n, front_np1, τ_start)
                    front_mid = interpolate_front(front_n, front_np1, τ_mid)
                    front_end = interpolate_front(front_n, front_np1, τ_end)
                    
                    # Calculate volumes at these times
                    _, vol_start, cx_start, cy_start, _ = fluid_cell_properties(nodes, front_start)
                    _, vol_mid, cx_mid, cy_mid, _ = fluid_cell_properties(nodes, front_mid)
                    _, vol_end, cx_end, cy_end, _ = fluid_cell_properties(nodes, front_end)
                    
                    # Extract volumes for this cell
                    v_start = vol_start[i, j]
                    v_mid = vol_mid[i, j]
                    v_end = vol_end[i, j]

                    # Weight centroids by volumes
                    wx_start = v_start > 0 ? cx_start[i, j] * v_start : 0.0
                    wy_start = v_start > 0 ? cy_start[i, j] * v_start : 0.0
                    
                    wx_mid = v_mid > 0 ? cx_mid[i, j] * v_mid : 0.0
                    wy_mid = v_mid > 0 ? cy_mid[i, j] * v_mid : 0.0
                    
                    wx_end = v_end > 0 ? cx_end[i, j] * v_end : 0.0
                    wy_end = v_end > 0 ? cy_end[i, j] * v_end : 0.0
                    
                    
                    # Apply Simpson's rule for this interval: (h/6)*(f₁ + 4f₂ + f₃)
                    h = t_end - t_start
                    sub_integral = (h/6) * (v_start + 4*v_mid + v_end)

                    # Weighted centroid integration
                    sub_integral_wx = (h/6) * (wx_start + 4*wx_mid + wx_end)
                    sub_integral_wy = (h/6) * (wy_start + 4*wy_mid + wy_end)
                    
                    # Add contribution to space-time volume
                    V_spacetime[i, j] += sub_integral
                    weighted_centroid_x += sub_integral_wx
                    weighted_centroid_y += sub_integral_wy
                end
                
                # Handle any remaining interval with trapezoidal rule
                if length(crossings) % 2 == 0
                    last_idx = length(crossings)
                    t_start = crossings[last_idx-1]
                    t_end = crossings[last_idx]
                    
                    # Calculate τ values
                    τ_start = t_start / dt
                    τ_end = t_end / dt
                    
                    # Create intermediate fronts
                    front_start = interpolate_front(front_n, front_np1, τ_start)
                    front_end = interpolate_front(front_n, front_np1, τ_end)
                    
                    # Calculate volumes
                    _, vol_start, cx_start, cy_start, _ = fluid_cell_properties(nodes, front_start)
                    _, vol_end,  cx_end, cy_end, _ = fluid_cell_properties(nodes, front_end)
                    
                    # Extract volumes
                    v_start = vol_start[i, j]
                    v_end = vol_end[i, j]

                    # Weight centroids by volumes
                    wx_start = v_start > 0 ? cx_start[i, j] * v_start : 0.0
                    wy_start = v_start > 0 ? cy_start[i, j] * v_start : 0.0
                    
                    wx_end = v_end > 0 ? cx_end[i, j] * v_end : 0.0
                    wy_end = v_end > 0 ? cy_end[i, j] * v_end : 0.0
                    
                     # Apply trapezoidal rule
                    h = t_end - t_start
                    sub_integral_v = h * (v_start + v_end) / 2
                    sub_integral_wx = h * (wx_start + wx_end) / 2
                    sub_integral_wy = h * (wy_start + wy_end) / 2
                    
                    # Add contributions
                    V_spacetime[i, j] += sub_integral_v
                    weighted_centroid_x += sub_integral_wx
                    weighted_centroid_y += sub_integral_wy
                end
            else
                # Use standard trapezoidal rule for simple cases
                for k in 1:(length(crossings)-1)
                    t_start = crossings[k]
                    t_end = crossings[k+1]
                    
                    # Calculate τ values
                    τ_start = t_start / dt
                    τ_end = t_end / dt
                    
                    # Create intermediate fronts
                    front_start = interpolate_front(front_n, front_np1, τ_start)
                    front_end = interpolate_front(front_n, front_np1, τ_end)
                    
                    # Calculate volumes
                    _, vol_start, cx_start, cy_start, _ = fluid_cell_properties(nodes, front_start)
                    _, vol_end, cx_end, cy_end, _ = fluid_cell_properties(nodes, front_end)
                    
                    # Extract volumes
                    v_start = vol_start[i, j]
                    v_end = vol_end[i, j]
                    
                    # Weight centroids by volumes
                    wx_start = v_start > 0 ? cx_start[i, j] * v_start : 0.0
                    wy_start = v_start > 0 ? cy_start[i, j] * v_start : 0.0
                    
                    wx_end = v_end > 0 ? cx_end[i, j] * v_end : 0.0
                    wy_end = v_end > 0 ? cy_end[i, j] * v_end : 0.0
                    
                    # Apply trapezoidal rule
                    h = t_end - t_start
                    sub_integral_v = h * (v_start + v_end) / 2
                    sub_integral_wx = h * (wx_start + wx_end) / 2
                    sub_integral_wy = h * (wy_start + wy_end) / 2
                    
                    # Add contributions
                    V_spacetime[i, j] += sub_integral_v
                    weighted_centroid_x += sub_integral_wx
                    weighted_centroid_y += sub_integral_wy
                end
            end
            
            # Safety bounds check
            V_spacetime[i, j] = clamp(V_spacetime[i, j], 0.0, cell_area * dt)
        end
    end

        # 4. Calculate Bx and By space-time capacities (wet lengths through cell centroids)
    for i in 1:nx
        for j in 1:ny
            # Get cell dimensions
            cell_width = x_nodes[i+1] - x_nodes[i]
            cell_height = y_nodes[j+1] - y_nodes[j]
            
            # Reuse the crossing times calculated for volumes
            crossings = Float64[]
            push!(crossings, 0.0)  # Add t_n
            push!(crossings, dt)   # Add t_{n+1}
            
            # Add crossing times from surrounding faces
            if t_crosses_x[i, j] > 0
                push!(crossings, t_crosses_x[i, j])
            end
            if t_crosses_x[i+1, j] > 0
                push!(crossings, t_crosses_x[i+1, j])
            end
            if t_crosses_y[i, j] > 0
                push!(crossings, t_crosses_y[i, j])
            end
            if t_crosses_y[i, j+1] > 0
                push!(crossings, t_crosses_y[i, j+1])
            end
            
            # Sort and remove duplicates
            sort!(unique!(crossings))
            
            # Initialize space-time capacities
            Bx_spacetime[i, j] = 0.0
            By_spacetime[i, j] = 0.0
            
            # For each time interval, integrate using trapezoidal rule
            for k in 1:(length(crossings)-1)
                t_start = crossings[k]
                t_end = crossings[k+1]
                
                # Calculate τ values
                τ_start = t_start / dt
                τ_end = t_end / dt
                
                # Create intermediate fronts
                front_start = interpolate_front(front_n, front_np1, τ_start)
                front_end = interpolate_front(front_n, front_np1, τ_end)
                
                # Calculate fluid properties at these times
                _, _, cx_start, cy_start, _ = fluid_cell_properties(nodes, front_start)
                _, _, cx_end, cy_end, _ = fluid_cell_properties(nodes, front_end)
                
                # Calculate Bx, By at start and end times
                _, _, Bx_start, By_start = compute_second_type_capacities(
                    nodes, front_start, cx_start, cy_start
                )
                _, _, Bx_end, By_end = compute_second_type_capacities(
                    nodes, front_end, cx_end, cy_end
                )
                
                # Trapezoidal rule integration
                h = t_end - t_start
                Bx_subinterval = h * (Bx_start[i, j] + Bx_end[i, j]) / 2
                By_subinterval = h * (By_start[i, j] + By_end[i, j]) / 2
                
                # Add contributions
                Bx_spacetime[i, j] += Bx_subinterval
                By_spacetime[i, j] += By_subinterval
            end
            
            # Safety bounds check
            Bx_spacetime[i, j] = clamp(Bx_spacetime[i, j], 0.0, cell_height * dt)
            By_spacetime[i, j] = clamp(By_spacetime[i, j], 0.0, cell_width * dt)
        end
    end
    
    # 5. Calculate Wx and Wy space-time capacities (fluid volumes between centroids)
    # Use Simpson's rule integration like for volumes
    for i in 1:nx
        for j in 1:ny
            # Calculate Wx_spacetime at i+1/2, j
            if i < nx
                # Get relevant cell dimensions
                y_min = y_nodes[j]
                y_max = y_nodes[j+1]
                
                # Collect crossing times for this interface
                wx_crossings = Float64[]
                push!(wx_crossings, 0.0)  # Add t_n
                push!(wx_crossings, dt)   # Add t_{n+1}
                
                # Add crossing times from surrounding faces
                if t_crosses_x[i, j] > 0
                    push!(wx_crossings, t_crosses_x[i, j])
                end
                if t_crosses_x[i+1, j] > 0
                    push!(wx_crossings, t_crosses_x[i+1, j])
                end
                if t_crosses_y[i, j] > 0
                    push!(wx_crossings, t_crosses_y[i, j])
                end
                if t_crosses_y[i, j+1] > 0
                    push!(wx_crossings, t_crosses_y[i, j+1])
                end
                if t_crosses_y[i+1, j] > 0
                    push!(wx_crossings, t_crosses_y[i+1, j])
                end
                if t_crosses_y[i+1, j+1] > 0
                    push!(wx_crossings, t_crosses_y[i+1, j+1])
                end
                
                # Sort and remove duplicates
                sort!(unique!(wx_crossings))
                
                # For significant volume changes, add additional midpoints for better sampling
                if length(wx_crossings) > 2
                    refined_crossings = Float64[]
                    for k in 1:(length(wx_crossings)-1)
                        push!(refined_crossings, wx_crossings[k])
                        # Add midpoint between each pair of crossing times
                        midpoint = 0.5 * (wx_crossings[k] + wx_crossings[k+1])
                        push!(refined_crossings, midpoint)
                    end
                    push!(refined_crossings, wx_crossings[end])
                    wx_crossings = sort(refined_crossings)
                end
                
                # Initialize space-time capacity
                Wx_spacetime[i+1, j] = 0.0
                
                # Perform quadrature integration with improved accuracy
                if length(wx_crossings) >= 3  # Simpson's rule requires at least 3 points
                    for k in 1:2:length(wx_crossings)-2
                        t_start = wx_crossings[k]
                        t_mid = wx_crossings[k+1]
                        t_end = wx_crossings[k+2]
                        
                        # Calculate τ values (normalized time)
                        τ_start = t_start / dt
                        τ_mid = t_mid / dt
                        τ_end = t_end / dt
                        
                        # Create intermediate fronts
                        front_start = interpolate_front(front_n, front_np1, τ_start)
                        front_mid = interpolate_front(front_n, front_np1, τ_mid)
                        front_end = interpolate_front(front_n, front_np1, τ_end)
                        
                        # Calculate fluid properties at these times
                        _, _, cx_start, cy_start, _ = fluid_cell_properties(nodes, front_start)
                        _, _, cx_mid, cy_mid, _ = fluid_cell_properties(nodes, front_mid)
                        _, _, cx_end, cy_end, _ = fluid_cell_properties(nodes, front_end)
                        
                        # Calculate Wx at start, mid, and end times
                        Wx_start, _, _, _ = compute_second_type_capacities(
                            nodes, front_start, cx_start, cy_start
                        )
                        Wx_mid, _, _, _ = compute_second_type_capacities(
                            nodes, front_mid, cx_mid, cy_mid
                        )
                        Wx_end, _, _, _ = compute_second_type_capacities(
                            nodes, front_end, cx_end, cy_end
                        )
                        
                        # Apply Simpson's rule: (h/6)*(f₁ + 4*f₂ + f₃)
                        h = t_end - t_start
                        Wx_subinterval = (h/6) * (Wx_start[i+1, j] + 4*Wx_mid[i+1, j] + Wx_end[i+1, j])
                        
                        # Add contribution
                        Wx_spacetime[i+1, j] += Wx_subinterval
                    end
                    
                    # Handle any remaining interval with trapezoidal rule
                    if length(wx_crossings) % 2 == 0
                        last_idx = length(wx_crossings)
                        t_start = wx_crossings[last_idx-1]
                        t_end = wx_crossings[last_idx]
                        
                        # Calculate τ values
                        τ_start = t_start / dt
                        τ_end = t_end / dt
                        
                        # Create intermediate fronts
                        front_start = interpolate_front(front_n, front_np1, τ_start)
                        front_end = interpolate_front(front_n, front_np1, τ_end)
                        
                        # Calculate fluid properties
                        _, _, cx_start, cy_start, _ = fluid_cell_properties(nodes, front_start)
                        _, _, cx_end, cy_end, _ = fluid_cell_properties(nodes, front_end)
                        
                        # Calculate Wx at start and end times
                        Wx_start, _, _, _ = compute_second_type_capacities(
                            nodes, front_start, cx_start, cy_start
                        )
                        Wx_end, _, _, _ = compute_second_type_capacities(
                            nodes, front_end, cx_end, cy_end
                        )
                        
                        # Apply trapezoidal rule
                        h = t_end - t_start
                        Wx_subinterval = h * (Wx_start[i+1, j] + Wx_end[i+1, j]) / 2
                        
                        # Add contribution
                        Wx_spacetime[i+1, j] += Wx_subinterval
                    end
                else
                    # Use standard trapezoidal rule for simple cases
                    for k in 1:(length(wx_crossings)-1)
                        t_start = wx_crossings[k]
                        t_end = wx_crossings[k+1]
                        
                        # Calculate τ values
                        τ_start = t_start / dt
                        τ_end = t_end / dt
                        
                        # Create intermediate fronts
                        front_start = interpolate_front(front_n, front_np1, τ_start)
                        front_end = interpolate_front(front_n, front_np1, τ_end)
                        
                        # Calculate fluid properties
                        _, _, cx_start, cy_start, _ = fluid_cell_properties(nodes, front_start)
                        _, _, cx_end, cy_end, _ = fluid_cell_properties(nodes, front_end)
                        
                        # Calculate Wx at start and end times
                        Wx_start, _, _, _ = compute_second_type_capacities(
                            nodes, front_start, cx_start, cy_start
                        )
                        Wx_end, _, _, _ = compute_second_type_capacities(
                            nodes, front_end, cx_end, cy_end
                        )
                        
                        # Apply trapezoidal rule
                        h = t_end - t_start
                        Wx_subinterval = h * (Wx_start[i+1, j] + Wx_end[i+1, j]) / 2
                        
                        # Add contribution
                        Wx_spacetime[i+1, j] += Wx_subinterval
                    end
                end
            end
            
            # Calculate Wy_spacetime at i, j+1/2 (following the same approach)
            if j < ny
                # Same structure as Wx but for Wy
                # [Similar code structure with Wy calculations]
                # Collect crossing times for this interface
                wy_crossings = Float64[]
                push!(wy_crossings, 0.0)  # Add t_n
                push!(wy_crossings, dt)   # Add t_{n+1}
                
                # Add crossing times from surrounding faces
                if t_crosses_x[i, j] > 0
                    push!(wy_crossings, t_crosses_x[i, j])
                end
                if t_crosses_x[i, j+1] > 0
                    push!(wy_crossings, t_crosses_x[i, j+1])
                end
                if t_crosses_y[i, j] > 0
                    push!(wy_crossings, t_crosses_y[i, j])
                end
                if t_crosses_y[i, j+1] > 0
                    push!(wy_crossings, t_crosses_y[i, j+1])
                end
                if t_crosses_x[i+1, j] > 0
                    push!(wy_crossings, t_crosses_x[i+1, j])
                end
                if t_crosses_x[i+1, j+1] > 0
                    push!(wy_crossings, t_crosses_x[i+1, j+1])
                end
                
                # Sort and remove duplicates
                sort!(unique!(wy_crossings))
                
                # For significant volume changes, add additional midpoints for better sampling
                if length(wy_crossings) > 2
                    refined_crossings = Float64[]
                    for k in 1:(length(wy_crossings)-1)
                        push!(refined_crossings, wy_crossings[k])
                        # Add midpoint between each pair of crossing times
                        midpoint = 0.5 * (wy_crossings[k] + wy_crossings[k+1])
                        push!(refined_crossings, midpoint)
                    end
                    push!(refined_crossings, wy_crossings[end])
                    wy_crossings = sort(refined_crossings)
                end
                
                # Initialize space-time capacity
                Wy_spacetime[i, j+1] = 0.0
                
                # Perform quadrature integration with improved accuracy
                if length(wy_crossings) >= 3  # Simpson's rule requires at least 3 points
                    for k in 1:2:length(wy_crossings)-2
                        t_start = wy_crossings[k]
                        t_mid = wy_crossings[k+1]
                        t_end = wy_crossings[k+2]
                        
                        # Calculate τ values (normalized time)
                        τ_start = t_start / dt
                        τ_mid = t_mid / dt
                        τ_end = t_end / dt
                        
                        # Create intermediate fronts
                        front_start = interpolate_front(front_n, front_np1, τ_start)
                        front_mid = interpolate_front(front_n, front_np1, τ_mid)
                        front_end = interpolate_front(front_n, front_np1, τ_end)
                        
                        # Calculate fluid properties at these times
                        _, _, cx_start, cy_start, _ = fluid_cell_properties(nodes, front_start)
                        _, _, cx_mid, cy_mid, _ = fluid_cell_properties(nodes, front_mid)
                        _, _, cx_end, cy_end, _ = fluid_cell_properties(nodes, front_end)
                        
                        # Calculate Wy at start, mid, and end times
                        _, Wy_start, _, _ = compute_second_type_capacities(
                            nodes, front_start, cx_start, cy_start
                        )
                        _, Wy_mid, _, _ = compute_second_type_capacities(
                            nodes, front_mid, cx_mid, cy_mid
                        )
                        _, Wy_end, _, _ = compute_second_type_capacities(
                            nodes, front_end, cx_end, cy_end
                        )
                        
                        # Apply Simpson's rule: (h/6)*(f₁ + 4*f₂ + f₃)
                        h = t_end - t_start
                        Wy_subinterval = (h/6) * (Wy_start[i, j+1] + 4*Wy_mid[i, j+1] + Wy_end[i, j+1])
                        
                        # Add contribution
                        Wy_spacetime[i, j+1] += Wy_subinterval
                    end
                    
                    # Handle any remaining interval with trapezoidal rule
                    if length(wy_crossings) % 2 == 0
                        last_idx = length(wy_crossings)
                        t_start = wy_crossings[last_idx-1]
                        t_end = wy_crossings[last_idx]
                        
                        # Calculate τ values
                        τ_start = t_start / dt
                        τ_end = t_end / dt
                        
                        # Create intermediate fronts
                        front_start = interpolate_front(front_n, front_np1, τ_start)
                        front_end = interpolate_front(front_n, front_np1, τ_end)
                        
                        # Calculate fluid properties
                        _, _, cx_start, cy_start, _ = fluid_cell_properties(nodes, front_start)
                        _, _, cx_end, cy_end, _ = fluid_cell_properties(nodes, front_end)
                        
                        # Calculate Wy at start and end times
                        _, Wy_start, _, _ = compute_second_type_capacities(
                            nodes, front_start, cx_start, cy_start
                        )
                        _, Wy_end, _, _ = compute_second_type_capacities(
                            nodes, front_end, cx_end, cy_end
                        )
                        
                        # Apply trapezoidal rule
                        h = t_end - t_start
                        Wy_subinterval = h * (Wy_start[i, j+1] + Wy_end[i, j+1]) / 2
                        
                        # Add contribution
                        Wy_spacetime[i, j+1] += Wy_subinterval
                    end
                else
                    # Use standard trapezoidal rule for simple cases
                    for k in 1:(length(wy_crossings)-1)
                        t_start = wy_crossings[k]
                        t_end = wy_crossings[k+1]
                        
                        # Calculate τ values
                        τ_start = t_start / dt
                        τ_end = t_end / dt
                        
                        # Create intermediate fronts
                        front_start = interpolate_front(front_n, front_np1, τ_start)
                        front_end = interpolate_front(front_n, front_np1, τ_end)
                        
                        # Calculate fluid properties
                        _, _, cx_start, cy_start, _ = fluid_cell_properties(nodes, front_start)
                        _, _, cx_end, cy_end, _ = fluid_cell_properties(nodes, front_end)
                        
                        # Calculate Wy at start and end times
                        _, Wy_start, _, _ = compute_second_type_capacities(
                            nodes, front_start, cx_start, cy_start
                        )
                        _, Wy_end, _, _ = compute_second_type_capacities(
                            nodes, front_end, cx_end, cy_end
                        )
                        
                        # Apply trapezoidal rule
                        h = t_end - t_start
                        Wy_subinterval = h * (Wy_start[i, j+1] + Wy_end[i, j+1]) / 2
                        
                        # Add contribution
                        Wy_spacetime[i, j+1] += Wy_subinterval
                    end
                end
            end
        end
    end

    return Dict(
        :Ax_spacetime => Ax_spacetime,
        :Ay_spacetime => Ay_spacetime,
        :V_spacetime => V_spacetime,
        :Bx_spacetime => Bx_spacetime,
        :By_spacetime => By_spacetime,
        :Wx_spacetime => Wx_spacetime,
        :Wy_spacetime => Wy_spacetime,
        :centroid_x_spacetime => centroid_x_spacetime,
        :centroid_y_spacetime => centroid_y_spacetime,
        :t_crosses_x => t_crosses_x,
        :t_crosses_y => t_crosses_y,
        :t_crosses_v => t_crosses_v   
    )
end

"""
    interpolate_front(front_n::FrontTracker, front_np1::FrontTracker, τ::Float64)

Create an intermediate front by linear interpolation between two fronts.
τ should be between 0 and 1 (normalized time).
"""
function interpolate_front(front_n::FrontTracker, front_np1::FrontTracker, τ::Float64)
    # Get markers from both fronts
    markers_n = get_markers(front_n)
    markers_np1 = get_markers(front_np1)
    
    # Ensure we have the same number of markers
    if length(markers_n) != length(markers_np1)
        error("Fronts must have the same number of markers for interpolation")
    end
    
    # Interpolate marker positions
    new_markers = [(
        (1-τ) * markers_n[i][1] + τ * markers_np1[i][1],
        (1-τ) * markers_n[i][2] + τ * markers_np1[i][2]
    ) for i in 1:length(markers_n)]
    
    # Create new front
    return FrontTracker(new_markers, front_n.is_closed)
end

"""
    calculate_face_capacity_x(front::FrontTracker, x::Float64, y_min::Float64, y_max::Float64)

Calculate the capacity (wet length) of a vertical face at position x from y_min to y_max.
"""
function calculate_face_capacity_x(front::FrontTracker, x::Float64, y_min::Float64, y_max::Float64)
    # Create a line representing the face
    face_line = LibGEOS.LineString([[x, y_min], [x, y_max]])
    
    # Get fluid polygon
    fluid_poly = get_fluid_polygon(front)
    
    # Calculate intersection
    if LibGEOS.intersects(face_line, fluid_poly)
        intersection = LibGEOS.intersection(face_line, fluid_poly)
        
        if isa(intersection, LibGEOS.LineString)
            return LibGEOS.geomLength(intersection)
        elseif isa(intersection, LibGEOS.MultiLineString)
            total_length = 0.0
            for k in 1:LibGEOS.getNumGeometries(intersection)
                line = LibGEOS.getGeometry(intersection, k-1)
                total_length += LibGEOS.geomLength(line)
            end
            return total_length
        elseif isa(intersection, LibGEOS.Point) || isa(intersection, LibGEOS.MultiPoint)
            return 0.0
        end
    else
        # Check if midpoint is inside
        mid_y = (y_min + y_max) / 2
        if is_point_inside(front, x, mid_y)
            return y_max - y_min
        else
            return 0.0
        end
    end
    
    return 0.0
end

"""
    calculate_face_capacity_y(front::FrontTracker, y::Float64, x_min::Float64, x_max::Float64)

Calculate the capacity (wet length) of a horizontal face at position y from x_min to x_max.
"""
function calculate_face_capacity_y(front::FrontTracker, y::Float64, x_min::Float64, x_max::Float64)
    # Create a line representing the face
    face_line = LibGEOS.LineString([[x_min, y], [x_max, y]])
    
    # Get fluid polygon
    fluid_poly = get_fluid_polygon(front)
    
    # Calculate intersection
    if LibGEOS.intersects(face_line, fluid_poly)
        intersection = LibGEOS.intersection(face_line, fluid_poly)
        
        if isa(intersection, LibGEOS.LineString)
            return LibGEOS.geomLength(intersection)
        elseif isa(intersection, LibGEOS.MultiLineString)
            total_length = 0.0
            for k in 1:LibGEOS.getNumGeometries(intersection)
                line = LibGEOS.getGeometry(intersection, k-1)
                total_length += LibGEOS.geomLength(line)
            end
            return total_length
        elseif isa(intersection, LibGEOS.Point) || isa(intersection, LibGEOS.MultiPoint)
            return 0.0
        end
    else
        # Check if midpoint is inside
        mid_x = (x_min + x_max) / 2
        if is_point_inside(front, mid_x, y)
            return x_max - x_min
        else
            return 0.0
        end
    end
    
    return 0.0
end

"""
    find_crossing_time(front_n::FrontTracker, front_np1::FrontTracker, x::Float64, y::Float64, dt::Float64)

Estimates the time when the interface crosses a specific point (x,y) between two time steps.
Uses linear interpolation of the signed distance function between time steps.

Parameters:
- front_n: Front tracker at time t_n
- front_np1: Front tracker at time t_n+1
- x, y: Coordinates of the point to check
- dt: Time step size

Returns:
- t_cross: Estimated crossing time within [0, dt] interval
"""
function find_crossing_time(front_n::FrontTracker, front_np1::FrontTracker, x::Float64, y::Float64, dt::Float64)
    # Get the signed distance at both time steps
    sdf_n = sdf(front_n, x, y)
    sdf_np1 = sdf(front_np1, x, y)
    
    # Check if the interface actually crosses this point
    if (sdf_n * sdf_np1 > 0)
        # No crossing (same sign at both time steps)
        if sdf_n < 0
            # Inside fluid at both times
            return dt/2
        else
            # Outside fluid at both times
            return dt/2
        end
    end
    
    # If one SDF is zero, return that time
    if abs(sdf_n) < 1e-10
        return 0.0
    elseif abs(sdf_np1) < 1e-10
        return dt
    end
    
    # Linear interpolation to find crossing time:
    t_cross = dt * abs(sdf_n) / (abs(sdf_n) + abs(sdf_np1))
    
    # Ensure the result is within [0, dt]
    return clamp(t_cross, 0.0, dt)
end

