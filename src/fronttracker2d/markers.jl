
function update_geometry!(ft::FrontTracker)
    if length(ft.markers) < 3
        # Not enough points to create a valid geometry
        ft.interface = nothing
        ft.interface_poly = nothing
        return
    end
    
    # Ensure interface is closed for polygon creation
    markers_to_use = ft.markers
    if ft.is_closed && length(markers_to_use) > 0 && markers_to_use[1] != markers_to_use[end]
        # Add closing point
        push!(markers_to_use, markers_to_use[1])
    end
    
    # For LibGEOS, we need to create the geometry directly from coordinates
    coords = [collect(point) for point in markers_to_use]
    
    # Create LineString
    ft.interface = LibGEOS.LineString(coords)
    
    # Create Polygon (if closed and valid)
    if ft.is_closed && length(ft.markers) >= 3
        ft.interface_poly = LibGEOS.Polygon([coords])
    else
        ft.interface_poly = nothing
    end
end

"""
Gets all markers of the interface.
"""
function get_markers(ft::FrontTracker)
    return ft.markers
end

"""
Adds a marker to the interface.
"""
function add_marker!(ft::FrontTracker, x::Float64, y::Float64)
    push!(ft.markers, (x, y))
    if length(ft.markers) >= 3  # Only update geometry when we have enough points
        update_geometry!(ft)
    end
    return ft
end

"""
Sets all markers of the interface.
"""
function set_markers!(ft::FrontTracker, markers::AbstractVector, is_closed=nothing)
    # Convert to proper type if needed
    typed_markers = convert(Vector{Tuple{Float64, Float64}}, markers)
    ft.markers = typed_markers
    if is_closed !== nothing
        ft.is_closed = is_closed
    end
    update_geometry!(ft)
    return ft
end

"""
Creates a circular interface.
"""
