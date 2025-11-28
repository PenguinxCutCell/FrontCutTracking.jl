function create_circle!(ft::FrontTracker, center_x::Float64, center_y::Float64, radius::Float64, n_markers::Int=100)
    # Create properly typed vector
    markers = Vector{Tuple{Float64, Float64}}(undef, n_markers)
    
    for i in 0:n_markers-1
        angle = 2.0 * π * (i / n_markers)
        x = center_x + radius * cos(angle)
        y = center_y + radius * sin(angle)
        markers[i+1] = (x, y)
    end
    
    set_markers!(ft, markers, true)
    return ft
end

"""
    create_rectangle!(ft::FrontTracker, min_x::Float64, min_y::Float64, max_x::Float64, max_y::Float64, n_markers::Int=100)

Creates a rectangular interface with markers distributed evenly around the perimeter.
Avoids placing markers at corners for numerical stability.
Starts in the middle of the bottom edge and traces the contour.
- min_x, min_y: Coordinates of bottom-left corner
- max_x, max_y: Coordinates of top-right corner
- n_markers: Number of points to place around the perimeter
"""
function create_rectangle!(ft::FrontTracker, min_x::Float64, min_y::Float64, max_x::Float64, max_y::Float64, n_markers::Int=100)
    # Calculate side lengths
    width = max_x - min_x
    height = max_y - min_y
    
    # Calculate perimeter
    perimeter = 2 * (width + height)
    
    # Calculate markers per side (proportional to side length)
    total_markers = n_markers - 4  # Reserve space to skip corners
    
    # Points per unit length (used for distribution)
    points_per_unit = total_markers / perimeter
    
    # Calculate number of points on each edge (excluding corners)
    n_bottom = max(1, round(Int, width * points_per_unit))
    n_right = max(1, round(Int, height * points_per_unit))
    n_top = max(1, round(Int, width * points_per_unit))
    n_left = max(1, round(Int, height * points_per_unit))
    
    # Adjust total to ensure we use exactly n_markers
    actual_total = n_bottom + n_right + n_top + n_left
    if actual_total < total_markers
        # Add extra points proportionally to the longest sides
        extra_points = total_markers - actual_total
        if width >= height
            extra_bottom = round(Int, extra_points / 2)
            extra_top = extra_points - extra_bottom
            n_bottom += extra_bottom
            n_top += extra_top
        else
            extra_right = round(Int, extra_points / 2)
            extra_left = extra_points - extra_right
            n_right += extra_right
            n_left += extra_left
        end
    elseif actual_total > total_markers
        # Remove points proportionally
        excess_points = actual_total - total_markers
        if width >= height
            excess_bottom = round(Int, excess_points / 2)
            excess_top = excess_points - excess_bottom
            n_bottom = max(1, n_bottom - excess_bottom)
            n_top = max(1, n_top - excess_top)
        else
            excess_right = round(Int, excess_points / 2)
            excess_left = excess_points - excess_right
            n_right = max(1, n_right - excess_right)
            n_left = max(1, n_left - excess_left)
        end
    end
    
    # Create properly typed vector (add +1 for closing point)
    markers = Vector{Tuple{Float64, Float64}}(undef, n_markers+1)
    
    # Start placing markers from the middle of bottom edge
    current_marker = 1
    
    # Bottom edge (left to right)
    for i in 1:n_bottom
        t = i / (n_bottom + 1)  # Skip 0 and 1 to avoid corners
        markers[current_marker] = (min_x + t * width, min_y)
        current_marker += 1
    end
    
    # Right edge (bottom to top)
    for i in 1:n_right
        t = i / (n_right + 1)  # Skip 0 and 1 to avoid corners
        markers[current_marker] = (max_x, min_y + t * height)
        current_marker += 1
    end
    
    # Top edge (right to left)
    for i in 1:n_top
        t = i / (n_top + 1)  # Skip 0 and 1 to avoid corners
        markers[current_marker] = (max_x - t * width, max_y)
        current_marker += 1
    end
    
    # Left edge (top to bottom)
    for i in 1:n_left
        t = i / (n_left + 1)  # Skip 0 and 1 to avoid corners
        markers[current_marker] = (min_x, max_y - t * height)
        current_marker += 1
    end
    
    # Ensure we have the right number of markers
    actual_markers = current_marker - 1
    if actual_markers != n_markers
        println("Warning: Requested $n_markers markers but placed $actual_markers markers")
    end
    
    # Close the curve by duplicating the first marker
    markers[actual_markers+1] = markers[1]
    
    # Set markers in the front tracker
    set_markers!(ft, markers[1:actual_markers+1], true)
    return ft
end

"""
    create_ellipse!(ft::FrontTracker, center_x::Float64, center_y::Float64, radius_x::Float64, radius_y::Float64, n_markers::Int=100)

Creates an elliptical interface.
"""
function create_ellipse!(ft::FrontTracker, center_x::Float64, center_y::Float64, radius_x::Float64, radius_y::Float64, n_markers::Int=100)
    theta = range(0, 2π, length=n_markers+1)[1:end-1] # Avoid duplicating the first point
    markers = [(center_x + radius_x*cos(t), center_y + radius_y*sin(t)) for t in theta]
    set_markers!(ft, markers, true)
    return ft
end

"""
    create_crystal!(ft::FrontTracker, center_x::Float64, center_y::Float64, 
                   base_radius::Float64, n_lobes::Int=6, amplitude::Float64=0.2, 
                   n_markers::Int=100)

Crée une interface en forme de cristal (un cercle avec perturbation angulaire).
- center_x, center_y: coordonnées du centre
- base_radius: rayon de base
- n_lobes: nombre de lobes du cristal (6 pour symétrie hexagonale)
- amplitude: amplitude de la perturbation (0-1)
- n_markers: nombre de points sur l'interface
"""
function create_crystal!(ft::FrontTracker, center_x::Float64, center_y::Float64, 
                        base_radius::Float64, n_lobes::Int=6, amplitude::Float64=0.2, 
                        n_markers::Int=100)
    # Créer un vecteur avec le type approprié
    markers = Vector{Tuple{Float64, Float64}}(undef, n_markers+1)
    
    for i in 1:n_markers
        # Angle pour ce marqueur
        θ = 2.0 * π * (i-1) / n_markers
        
        # Rayon perturbé pour un effet cristallin
        r = base_radius * (1.0 + amplitude * cos(n_lobes * θ))
        
        # Coordonnées cartésiennes
        x = center_x + r * cos(θ)
        y = center_y + r * sin(θ)
        
        markers[i] = (x, y)
    end
    
    # Fermer la courbe en répétant le premier point
    markers[n_markers+1] = markers[1]
    
    # Définir les marqueurs dans le FrontTracker
    set_markers!(ft, markers, true)
    return ft
end

"""
    get_fluid_polygon(ft::FrontTracker)

Returns a polygon representing the fluid domain bounded by the interface.
"""
