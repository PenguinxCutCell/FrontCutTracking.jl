function compute_segment_parameters(ft::FrontTracker)
    markers = ft.markers
    n_markers = length(markers)
    
    if n_markers < 2
        return [], [], [], [], []
    end
    
    # Nombre de segments (n_markers pour un contour fermé, n_markers-1 pour un contour ouvert)
    n_segments = if ft.is_closed
        if n_markers > 0 && markers[1] == markers[end]
            n_markers - 1
        else
            n_markers
        end
    else
        n_markers - 1
    end
    
    # Initialiser les tableaux de résultats
    segments = Vector{Tuple{Int, Int}}(undef, n_segments)
    segment_normals = Vector{Tuple{Float64, Float64}}(undef, n_segments)
    segment_intercepts = Vector{Float64}(undef, n_segments)
    segment_lengths = Vector{Float64}(undef, n_segments)
    segment_midpoints = Vector{Tuple{Float64, Float64}}(undef, n_segments)
    
    # Parcourir tous les segments
    for i in 1:n_segments
        # Indice du marqueur suivant (avec gestion de la boucle pour contour fermé)
        next_i = i < n_markers ? i + 1 : 1
        
        # Points de début et de fin du segment
        p1 = markers[i]
        p2 = markers[next_i]
        
        # Vecteur du segment (de p1 à p2)
        segment_vector = (p2[1] - p1[1], p2[2] - p1[2])
        
        # Longueur du segment
        segment_length = sqrt(segment_vector[1]^2 + segment_vector[2]^2)
        
        if segment_length < 1e-15
            # Éviter la division par zéro pour les segments très courts
            segment_normals[i] = (0.0, 1.0)  # Normal arbitraire
            segment_intercepts[i] = p1[1] * 0.0 + p1[2] * 1.0  # α_I = n_I ⋅ p1
            segment_lengths[i] = 0.0
            segment_midpoints[i] = p1  # Point de milieu = point de début pour segments courts
        else
            # Normale unitaire (rotation de 90° dans le sens trigonométrique)
            normal = (-segment_vector[2] / segment_length, segment_vector[1] / segment_length)
            
            # Vérifier l'orientation de la normale (doit pointer à l'extérieur)
            if ft.is_closed
                # Pour un contour fermé, vérifier si la normale pointe vers l'extérieur
                test_point = (p1[1] + 1e-3 * normal[1], p1[2] + 1e-3 * normal[2])
                if is_point_inside(ft, test_point[1], test_point[2])
                    # Si le point test est à l'intérieur, inverser la normale
                    normal = (-normal[1], -normal[2])
                end
            end
            
            # Calcul de l'intercept α_I = n_I ⋅ p1
            intercept = normal[1] * p1[1] + normal[2] * p1[2]
            
            # Stockage des résultats
            segments[i] = (i, next_i)
            segment_normals[i] = normal
            segment_intercepts[i] = intercept
            segment_lengths[i] = segment_length
            segment_midpoints[i] = ((p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2)
        end
    end
    
    return segments, segment_normals, segment_intercepts, segment_lengths, segment_midpoints
end

"""
    create_segment_line(ft::FrontTracker, segment_idx::Int)

Crée une LineString représentant un segment de l'interface à partir des marqueurs.
"""
function create_segment_line(ft::FrontTracker, segment_idx::Int)
    markers = ft.markers
    n_markers = length(markers)
    
    if segment_idx < 1 || segment_idx > (ft.is_closed ? n_markers : n_markers - 1)
        error("Indice de segment invalide: $segment_idx")
    end
    
    # Récupérer les indices des marqueurs qui définissent le segment
    next_idx = segment_idx < n_markers ? segment_idx + 1 : 1
    
    # Récupérer les coordonnées des marqueurs
    start_point = markers[segment_idx]
    end_point = markers[next_idx]
    
    # Créer la LineString
    return LibGEOS.LineString([[start_point[1], start_point[2]], [end_point[1], end_point[2]]])
end

"""
    compute_segment_cell_intersections(nodes::NTuple{2, AbstractVector}, ft::FrontTracker)

Calcule les intersections entre les segments de l'interface et les cellules du maillage.
Retourne un dictionnaire où les clés sont les indices de cellules (i,j) et les valeurs
sont des listes de tuples (segment_idx, intersection_length).
"""
function compute_segment_cell_intersections(nodes::NTuple{2, AbstractVector}, ft::FrontTracker)
    # Calculer les paramètres des segments
    segments, segment_normals, segment_intercepts, segment_lengths, segment_midpoints = 
        compute_segment_parameters(ft)
    
    x_nodes = nodes[1]
    y_nodes = nodes[2]
    nx = length(x_nodes) - 1
    ny = length(y_nodes) - 1
    
    # Dictionnaire pour stocker les intersections segment-cellule
    cell_segment_intersections = Dict{Tuple{Int, Int}, Vector{Tuple{Int, Float64}}}()
    
    # Initialiser le dictionnaire pour toutes les cellules
    for i in 1:nx
        for j in 1:ny
            cell_segment_intersections[(i,j)] = []
        end
    end
    
    # Pour chaque segment, calculer les intersections avec les cellules
    n_segments = length(segments)
    for segment_idx in 1:n_segments
        # Créer une ligne représentant le segment
        segment_line = create_segment_line(ft, segment_idx)
        
        # Calculer les intersections avec toutes les cellules
        for i in 1:nx
            for j in 1:ny
                # Créer le polygone de la cellule
                cell_coords = [
                    [x_nodes[i], y_nodes[j]],
                    [x_nodes[i+1], y_nodes[j]],
                    [x_nodes[i+1], y_nodes[j+1]],
                    [x_nodes[i], y_nodes[j+1]],
                    [x_nodes[i], y_nodes[j]]
                ]
                cell_poly = LibGEOS.Polygon([cell_coords])
                
                # Vérifier l'intersection
                if LibGEOS.intersects(cell_poly, segment_line)
                    intersection = LibGEOS.intersection(cell_poly, segment_line)
                    
                    # Calculer la longueur d'intersection
                    if isa(intersection, LibGEOS.LineString)
                        intersection_length = LibGEOS.geomLength(intersection)
                        if intersection_length > 1e-10
                            push!(cell_segment_intersections[(i,j)], (segment_idx, intersection_length))
                        end
                    elseif isa(intersection, LibGEOS.MultiLineString)
                        total_length = 0.0
                        for k in 1:LibGEOS.getNumGeometries(intersection)
                            line = LibGEOS.getGeometry(intersection, k-1)
                            total_length += LibGEOS.geomLength(line)
                        end
                        if total_length > 1e-10
                            push!(cell_segment_intersections[(i,j)], (segment_idx, total_length))
                        end
                    end
                end
            end
        end
    end
    
    return cell_segment_intersections, segments, segment_normals, segment_intercepts, segment_lengths
end

"""
    compute_intercept_jacobian(nodes::NTuple{2, AbstractVector}, ft::FrontTracker; density::Float64=1.0)

Calcule la jacobienne des volumes par rapport aux déplacements d'intercept.
Pour chaque cellule k=(i,j) et chaque segment I, J[k,I] = ∂V_k/∂δ_I = ρL × A_k,I,
où A_k,I est la longueur d'intersection du segment I avec la cellule k,
et ρL est un facteur physique (densité × latent heat)

Retourne:
- intercept_jacobian: Dict{Tuple{Int, Int}, Vector{Tuple{Int, Float64}}} - la jacobienne
- segments: vecteur des segments (i, j) où i, j sont les indices des marqueurs
- segment_normals: vecteur des normales unitaires pour chaque segment
- segment_intercepts: vecteur des intercepts initiaux pour chaque segment
- segment_lengths: vecteur des longueurs de chaque segment
"""
