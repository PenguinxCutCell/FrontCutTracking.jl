"""
    get_markers(ft::FrontTracker3D)

Gets all markers of the 3D interface.
"""
function get_markers(ft::FrontTracker3D)
    return ft.markers
end

"""
    get_faces(ft::FrontTracker3D)

Gets all triangular faces of the 3D interface.
"""
function get_faces(ft::FrontTracker3D)
    return ft.faces
end

"""
    add_marker!(ft::FrontTracker3D, x::Float64, y::Float64, z::Float64)

Adds a marker to the 3D interface.
"""
function add_marker!(ft::FrontTracker3D, x::Float64, y::Float64, z::Float64)
    push!(ft.markers, (x, y, z))
    return ft
end

"""
    set_markers!(ft::FrontTracker3D, markers::AbstractVector, is_closed=nothing)

Sets all markers of the 3D interface.
"""
function set_markers!(ft::FrontTracker3D, markers::AbstractVector, is_closed=nothing)
    # Convert to proper type if needed
    typed_markers = convert(Vector{Tuple{Float64, Float64, Float64}}, markers)
    ft.markers = typed_markers
    if is_closed !== nothing
        ft.is_closed = is_closed
    end
    return ft
end

"""
    set_faces!(ft::FrontTracker3D, faces::Vector{Tuple{Int, Int, Int}})

Sets all triangular faces of the 3D interface.
"""
function set_faces!(ft::FrontTracker3D, faces::Vector{Tuple{Int, Int, Int}})
    ft.faces = faces
    return ft
end

"""
    add_face!(ft::FrontTracker3D, i1::Int, i2::Int, i3::Int)

Adds a triangular face to the 3D interface.
"""
function add_face!(ft::FrontTracker3D, i1::Int, i2::Int, i3::Int)
    push!(ft.faces, (i1, i2, i3))
    return ft
end
