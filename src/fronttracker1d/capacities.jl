function compute_capacities_1d(nodes::NTuple{1, AbstractVector}, front::FrontTracker1D)
    # Extract mesh information
    x_nodes = nodes[1]
    nx = length(x_nodes) - 1
    
    # Initialize capacity arrays
    fractions = zeros(nx+1)    # Fluid fractions
    volumes = zeros(nx+1)      # Fluid volumes (lengths in 1D)
    centroids_x = zeros(nx+1)  # Fluid centroids
    cell_types = zeros(Int, nx+1)  # Cell types (0: solid, 1: fluid, -1: cut)
    Ax = zeros(nx+1)           # Face capacities
    Wx = zeros(nx+1)           # Staggered volumes 
    Bx = zeros(nx+1)           # Center-line capacities
    
    # Interface information
    interface_positions = Dict{Int, Float64}()
    
    # Process each cell
    for i in 1:nx
        x_min, x_max = x_nodes[i], x_nodes[i+1]
        cell_width = x_max - x_min
        
        # Find interface points inside this cell
        interface_points_in_cell = filter(m -> x_min <= m <= x_max, front.markers)
        
        # Check if the cell contains interface points
        if isempty(interface_points_in_cell)
            # Cell is either fully fluid or fully solid
            if is_point_inside(front, (x_min + x_max)/2)
                # Fully fluid cell
                fractions[i] = 1.0
                volumes[i] = cell_width
                centroids_x[i] = (x_min + x_max) / 2
                cell_types[i] = 1
                
                # For fully fluid cells:
                Bx[i] = 1.0  # Full capacity at center
                
                # Wx values will be computed after processing all cells
            else
                # Fully solid cell
                fractions[i] = 0.0
                volumes[i] = 0.0
                centroids_x[i] = (x_min + x_max) / 2
                cell_types[i] = 0
                
                # For fully solid cells:
                Bx[i] = 0.0  # No capacity at center
            end
        else
            # Cut cell
            cell_types[i] = -1
            
            # Store interface positions for later use
            for position in interface_points_in_cell
                interface_positions[i] = position
            end
            
            # Sort interface points within the cell
            sort!(interface_points_in_cell)
            
            # Calculate fluid fraction and volume
            fluid_segments = []
            
            # Start with cell left boundary
            current_x = x_min
            is_fluid = is_point_inside(front, current_x + 1e-10)
            
            # Process each interface point and collect fluid segments
            for point in interface_points_in_cell
                if is_fluid
                    # Add fluid segment from current_x to interface point
                    push!(fluid_segments, (current_x, point))
                end
                
                # Toggle fluid state
                is_fluid = !is_fluid
                current_x = point
            end
            
            # Process the last segment to the right boundary
            if is_fluid
                push!(fluid_segments, (current_x, x_max))
            end
            
            # Calculate total fluid volume and centroid
            total_fluid_volume = 0.0
            weighted_centroid_x = 0.0
            
            for (start_x, end_x) in fluid_segments
                segment_length = end_x - start_x
                total_fluid_volume += segment_length
                
                # Weighted contribution to centroid
                segment_centroid = (start_x + end_x) / 2
                weighted_centroid_x += segment_length * segment_centroid
            end
            
            # Store computed values
            volumes[i] = total_fluid_volume
            fractions[i] = total_fluid_volume / cell_width
            
            if total_fluid_volume > 0
                centroids_x[i] = weighted_centroid_x / total_fluid_volume
            else
                centroids_x[i] = (x_min + x_max) / 2
            end
            
            # Calculate Bx for cut cells (depends on position of centroid relative to interface)
            if fractions[i] > 0
                mid_point = centroids_x[i]
                Bx[i] = is_point_inside(front, mid_point) ? 1.0 : 0.0
            else
                Bx[i] = 0.0
            end
        end
    end
    
    # Calculate face capacities (Ax) at cell boundaries
    for i in 1:nx+1
        x_face = x_nodes[i]
        Ax[i] = is_point_inside(front, x_face) ? 1.0 : 0.0
    end
    
    # Calculate staggered volumes (Wx) between cell centers
    for i in 2:nx
        # Position between cell centroids
        x_left = centroids_x[i-1] 
        x_right = centroids_x[i]
        
        # If either adjacent cell has zero volume, Wx is zero
        if volumes[i-1] == 0.0 && volumes[i] == 0.0
            Wx[i] = 0.0
            continue
        end
        
        # Check for interfaces between the centroids
        interfaces_between = Float64[]
        for marker in front.markers
            if x_left < marker < x_right
                push!(interfaces_between, marker)
            end
        end
        
        if isempty(interfaces_between)
            # No interface between centroids, full connectivity
            Wx[i] = x_right - x_left
        else
            # For cells with interface, calculate connection based on fluid side
            # Find the interface closest to either centroid
            closest_interface = interfaces_between[1]
            min_distance = min(abs(closest_interface - x_left), abs(closest_interface - x_right))
            
            for marker in interfaces_between
                dist_left = abs(marker - x_left)
                dist_right = abs(marker - x_right)
                if min(dist_left, dist_right) < min_distance
                    closest_interface = marker
                    min_distance = min(dist_left, dist_right)
                end
            end
            
            # Determine which side is fluid
            if is_point_inside(front, x_left)
                # Left centroid is in fluid, Wx is distance from left centroid to interface
                Wx[i] = closest_interface - x_left
            elseif is_point_inside(front, x_right)
                # Right centroid is in fluid, Wx is distance from interface to right centroid
                Wx[i] = x_right - closest_interface
            else
                # Neither centroid is in fluid, no connection
                Wx[i] = 0.0
            end
        end
    end
    
    # Return all capacities in a dictionary
    return Dict(
        :fractions => fractions,          # Fluid fractions
        :volumes => volumes,              # Fluid volumes (lengths in 1D)
        :centroids_x => centroids_x,      # Fluid centroids
        :cell_types => cell_types,        # Cell types (0: solid, 1: fluid, -1: cut)
        :Ax => Ax,                        # Face capacities
        :Wx => Wx,                        # Staggered volumes
        :Bx => Bx,                        # Center line capacities
        :interface_positions => interface_positions  # Interface positions
    )
end

"""
    compute_spacetime_capacities_1d(nodes::NTuple{1, AbstractVector}, front_n::FrontTracker1D, front_np1::FrontTracker1D, dt::Float64)

Calculates all space-time geometric capacities for a 1D mesh between two time steps.
Now includes Wx time-integrated connection capacities using centroids.
"""
