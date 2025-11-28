function is_point_inside(ft::FrontTracker1D, x::Float64)
    if isempty(ft.markers)
        return false
    end
    
    # Find number of interface points to the left of x
    count = sum(m < x for m in ft.markers)
    
    # If odd number of interface points to the left, point is inside
    return count % 2 == 1
end

"""
    sdf(ft::FrontTracker1D, x::Float64)

Calculates the signed distance function for a given point.
Positive outside fluid, negative inside fluid.
"""
function sdf(ft::FrontTracker1D, x::Float64)
    if isempty(ft.markers)
        return Inf
    end
    
    # Calculate distance to nearest interface point
    distance = minimum(abs.(ft.markers .- x))
    
    # Determine sign (negative inside, positive outside)
    sign_val = is_point_inside(ft, x) ? -1.0 : 1.0
    
    return sign_val * distance
end

"""
    compute_capacities_1d(nodes::NTuple{1, AbstractVector}, front::FrontTracker1D)

Calculates all geometric capacities for a 1D mesh described by its node coordinates.
Returns a dictionary with all results.
"""
