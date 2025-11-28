"""
A Julia implementation of 3D front tracking for fluid-solid interfaces.
In 3D, the interface is a closed surface (mesh) made of triangular faces.
"""
mutable struct FrontTracker3D
    # Markers are 3D coordinates (x, y, z)
    markers::Vector{Tuple{Float64, Float64, Float64}}
    # Triangular faces defined by indices into the markers array
    faces::Vector{Tuple{Int, Int, Int}}
    # Whether the surface is closed
    is_closed::Bool
    
    # Constructor for empty front
    function FrontTracker3D()
        return new([], [], true)
    end
    
    # Constructor with markers and faces
    function FrontTracker3D(markers::Vector{Tuple{Float64, Float64, Float64}}, 
                            faces::Vector{Tuple{Int, Int, Int}}, 
                            is_closed::Bool=true)
        return new(markers, faces, is_closed)
    end
    
    # Constructor with markers only (no faces)
    function FrontTracker3D(markers::Vector{Tuple{Float64, Float64, Float64}}, is_closed::Bool=true)
        return new(markers, [], is_closed)
    end
end
