function compute_intercept_jacobian(nodes::NTuple{2, AbstractVector}, ft::FrontTracker; density::Float64=1.0)
    # Calculer les intersections segment-cellule
    cell_segment_intersections, segments, segment_normals, segment_intercepts, segment_lengths = 
        compute_segment_cell_intersections(nodes, ft)
    
    # Récupérer les dimensions du maillage
    x_nodes = nodes[1]
    y_nodes = nodes[2]
    nx = length(x_nodes) - 1
    ny = length(y_nodes) - 1
    
    # Créer un dictionnaire pour stocker la jacobienne
    # Format: Dict{Tuple{Int, Int}, Vector{Tuple{Int, Float64}}}
    # Clé: (i,j) = indice de la cellule, Valeur: liste de (segment_idx, jacobian_value)
    intercept_jacobian = Dict{Tuple{Int, Int}, Vector{Tuple{Int, Float64}}}()
    
    # Pour chaque cellule, calculer sa contribution à la jacobienne
    for i in 1:nx
        for j in 1:ny
            intercept_jacobian[(i,j)] = []
            
            # Parcourir tous les segments qui intersectent cette cellule
            for (segment_idx, intersection_length) in cell_segment_intersections[(i,j)]
                # J[k,I] = ρL × A_k,I
                jacobian_value = density * intersection_length
                
                # Stocker la valeur dans la jacobienne
                push!(intercept_jacobian[(i,j)], (segment_idx, jacobian_value))
            end
        end
    end
    
    return intercept_jacobian, segments, segment_normals, segment_intercepts, segment_lengths
end

"""
    update_front_with_intercept_displacements!(ft::FrontTracker, displacements::AbstractVector{<:Real}, 
                                            segment_normals::Vector{Tuple{Float64, Float64}},
                                            segment_lengths::Vector{Float64})

Met à jour l'interface en déplaçant chaque segment selon ses déplacements d'intercept.
Applique une pondération basée sur la longueur des segments pour distribuer les déplacements
aux marqueurs partagés entre plusieurs segments.

Paramètres:
- ft: l'objet FrontTracker à mettre à jour
- displacements: vecteur des déplacements δ_I pour chaque segment
- segment_normals: vecteur des normales unitaires pour chaque segment
- segment_lengths: vecteur des longueurs de chaque segment

Retourne le FrontTracker mis à jour.
"""
function update_front_with_intercept_displacements!(ft::FrontTracker, displacements::AbstractVector{<:Real}, 
                                                  segment_normals::Vector{Tuple{Float64, Float64}},
                                                  segment_lengths::Vector{Float64})
    markers = copy(ft.markers)
    n_markers = length(markers)
    n_segments = length(displacements)
    
  
    
    # Structure pour stocker les contributions pondérées de chaque segment à chaque marqueur
    segment_contributions = Dict{Int, Vector{Tuple{Float64, Tuple{Float64, Float64}}}}()
    for i in 1:n_markers
        segment_contributions[i] = []
    end
    
    # Pour chaque segment, calculer sa contribution aux marqueurs
    for (s_idx, displacement) in enumerate(displacements)
        # Récupérer les indices des marqueurs aux extrémités du segment
        start_idx = s_idx
        end_idx = s_idx < n_markers ? s_idx + 1 : 1
        
        # Calculer le vecteur de déplacement
        normal = segment_normals[s_idx]
        vector_displacement = (displacement * normal[1], displacement * normal[2])
        
        # Utiliser la longueur du segment comme poids
        segment_weight = max(segment_lengths[s_idx], 1e-10)  # Éviter division par zéro
        
        # Enregistrer la contribution pondérée de ce segment pour chaque marqueur
        push!(segment_contributions[start_idx], (segment_weight, vector_displacement))
        push!(segment_contributions[end_idx], (segment_weight, vector_displacement))
    end
    
    # Calculer et appliquer le déplacement moyen pondéré pour chaque marqueur
    for i in 1:n_markers
        contributions = segment_contributions[i]
        if !isempty(contributions)
            # Calculer la somme des poids
            total_weight = sum(contrib[1] for contrib in contributions)
            
            if total_weight > 0
                # Calculer le déplacement moyen pondéré
                avg_dx = sum(contrib[1] * contrib[2][1] for contrib in contributions) / total_weight
                avg_dy = sum(contrib[1] * contrib[2][2] for contrib in contributions) / total_weight
                
                # Appliquer le déplacement
                markers[i] = (markers[i][1] + avg_dx, markers[i][2] + avg_dy)
            end
        end
    end
    
    # Mettre à jour l'interface avec les nouveaux marqueurs
    set_markers!(ft, markers, ft.is_closed)
    
    return ft
end
