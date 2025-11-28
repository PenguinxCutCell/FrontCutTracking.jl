"""
A Julia implementation of 1D front tracking for fluid-solid interfaces.
"""
mutable struct FrontTracker1D
    # In 1D, markers are just x-coordinates of interface points
    markers::Vector{Float64}
    
    # Constructor for empty front
    function FrontTracker1D()
        return new([])
    end
    
    # Constructor with markers
    function FrontTracker1D(markers::Vector{Float64})
        return new(sort(markers))  # Keep markers sorted for easier processing
    end
end

"""
    get_markers(ft::FrontTracker1D)

Gets all markers of the 1D interface.
"""
function get_markers(ft::FrontTracker1D)
    return ft.markers
end

"""
    add_marker!(ft::FrontTracker1D, x::Float64)

Adds a marker to the 1D interface and keeps the collection sorted.
"""
function add_marker!(ft::FrontTracker1D, x::Float64)
    push!(ft.markers, x)
    sort!(ft.markers)
    return ft
end

"""
    set_markers!(ft::FrontTracker1D, markers::AbstractVector{<:Real})

Sets all markers of the 1D interface.
"""
function set_markers!(ft::FrontTracker1D, markers::AbstractVector{<:Real})
    ft.markers = convert(Vector{Float64}, markers)
    sort!(ft.markers)
    return ft
end

"""
    is_point_inside(ft::FrontTracker1D, x::Float64)

Checks if a point is inside the fluid region.
For 1D, we define "inside" as being to the left of an odd-indexed interface point
or to the right of an even-indexed interface point.
"""
