module FrontCutTracking

using GeoInterface
using LibGEOS

include("fronttracker/types.jl")
include("fronttracker/markers.jl")
include("fronttracker/shapes.jl")
include("fronttracker/geometry.jl")
include("fronttracker/capacities.jl")
include("fronttracker/spacetime.jl")
include("fronttracker/segments.jl")
include("fronttracker/intercepts.jl")

include("fronttracker1d/types.jl")
include("fronttracker1d/geometry.jl")
include("fronttracker1d/capacities.jl")
include("fronttracker1d/spacetime.jl")

export FrontTracker,
	   create_circle!,
	   create_rectangle!,
	   create_ellipse!,
	   update_geometry!,
	   create_crystal!,
	   set_markers!,
	   get_markers,
	   add_marker!,
	   get_fluid_polygon,
	   is_point_inside,
	   get_intersection,
	   sdf,
	   compute_marker_normals,
	   compute_volume_jacobian,
	   compute_capacities,
	   fluid_cell_properties,
	   compute_surface_capacities,
	   compute_second_type_capacities,
	   compute_intercept_jacobian,
	   compute_segment_cell_intersections,
	   create_segment_line,
	   compute_segment_parameters,
	   update_front_with_intercept_displacements!,
	   compute_spacetime_capacities,
	   FrontTracker1D,
	   compute_capacities_1d,
	   compute_spacetime_capacities_1d

end # module
