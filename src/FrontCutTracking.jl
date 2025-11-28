module FrontCutTracking

using GeoInterface
using LibGEOS

include("fronttracker2d/types.jl")
include("fronttracker2d/markers.jl")
include("fronttracker2d/shapes.jl")
include("fronttracker2d/geometry.jl")
include("fronttracker2d/capacities.jl")
include("fronttracker2d/spacetime.jl")
include("fronttracker2d/segments.jl")
include("fronttracker2d/intercepts.jl")

include("fronttracker1d/types.jl")
include("fronttracker1d/geometry.jl")
include("fronttracker1d/capacities.jl")
include("fronttracker1d/spacetime.jl")

include("fronttracker3d/types.jl")
include("fronttracker3d/markers.jl")
include("fronttracker3d/shapes.jl")
include("fronttracker3d/geometry.jl")
include("fronttracker3d/segments.jl")
include("fronttracker3d/intercepts.jl")

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
	   compute_spacetime_capacities_1d,
	   # 3D FrontTracker exports
	   FrontTracker3D,
	   get_faces,
	   set_faces!,
	   add_face!,
	   create_sphere!,
	   create_box!,
	   create_ellipsoid!,
	   create_cylinder!,
	   compute_face_normal,
	   compute_face_area,
	   compute_total_surface_area,
	   compute_enclosed_volume,
	   compute_centroid,
	   compute_edge_parameters,
	   compute_face_parameters,
	   get_adjacent_faces,
	   get_adjacent_markers,
	   compute_vertex_area,
	   compute_mean_curvature,
	   compute_face_cell_intersections,
	   compute_intercept_jacobian_3d,
	   update_front_with_intercept_displacements_3d!,
	   compute_volume_jacobian_3d

end # module
