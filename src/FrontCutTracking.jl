module FrontCutTracking

using WriteVTK: vtk_grid, vtk_save, MeshCell, VTKCellTypes

export FrontTracker,
       add_point!,
       add_line!,
       add_surface!,
       remove_point!,
       remove_line!,
       remove_surface!,
       num_points,
       num_lines,
       num_surfaces,
       lines_touching,
       surfaces_touching,
    reset!,
    point_front,
    line_front,
    circle_front,
    planar_patch_front,
    sphere_front,
    star_front,
    star3d_front,
    write_vtk_front

const Point3D{T} = NTuple{3,T}
const Line = NTuple{2,Int}
const Surface = NTuple{3,Int}

"""
    FrontTracker(; T::Type{<:Real}=Float64, metadata=Dict{Symbol,Any}())

Mutable container that stores 0D points, 1D line segments, and 2D triangular
surfaces for front-tracking style algorithms. All points share the same floating
point element type `T`.
"""
mutable struct FrontTracker{T<:Real}
    points::Vector{Point3D{T}}
    lines::Vector{Line}
    surfaces::Vector{Surface}
    metadata::Dict{Symbol,Any}

    function FrontTracker{T}(points::Vector{Point3D{T}},
                             lines::Vector{Line},
                             surfaces::Vector{Surface},
                             metadata::Dict{Symbol,Any}) where {T<:Real}
        new{T}(points, lines, surfaces, metadata)
    end
end
FrontTracker(; T::Type{<:Real}=Float64, metadata=Dict{Symbol,Any}()) =
    FrontTracker{T}(Point3D{T}[], Line[], Surface[], _metadata_dict(metadata))

Base.copy(tracker::FrontTracker{T}) where {T} =
    FrontTracker{T}(copy(tracker.points), copy(tracker.lines), copy(tracker.surfaces), copy(tracker.metadata))

num_points(tracker::FrontTracker) = length(tracker.points)
num_lines(tracker::FrontTracker) = length(tracker.lines)
num_surfaces(tracker::FrontTracker) = length(tracker.surfaces)

"""
    reset!(tracker)

Remove all entities and metadata from `tracker` in-place.
"""
function reset!(tracker::FrontTracker)
    empty!(tracker.points)
    empty!(tracker.lines)
    empty!(tracker.surfaces)
    empty!(tracker.metadata)
    tracker
end

# -- Point helpers -----------------------------------------------------------

function _point_tuple(point, ::Type{T}) where {T}
    length(point) == 3 || throw(ArgumentError("Points must have three coordinates."))
    ntuple(i -> convert(T, getindex(point, i)), 3)
end

_point_tuple(point::Point3D{Ti}, ::Type{T}) where {Ti,T} =
    ntuple(i -> convert(T, point[i]), 3)

_tuple_add(a, b) = (a[1] + b[1], a[2] + b[2], a[3] + b[3])
_tuple_scale(a, s) = (a[1] * s, a[2] * s, a[3] * s)

function _normalize_direction(direction)
    vec = _point_tuple(direction, Float64)
    mag = sqrt(vec[1]^2 + vec[2]^2 + vec[3]^2)
    mag > 0 || throw(ArgumentError("Direction vector must be non-zero."))
    _tuple_scale(vec, 1 / mag)
end

function _plane_axes(axis::Symbol)
    axis === :x && return ((0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    axis === :y && return ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0))
    axis === :z && return ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
    throw(ArgumentError("Axis must be one of :x, :y, or :z."))
end

_metadata_dict(metadata::Dict{Symbol,Any}) = metadata
_metadata_dict(metadata) = Dict{Symbol,Any}(metadata)

function _points_matrix(points::Vector{Point3D{T}}) where {T}
    if isempty(points)
        return Matrix{T}(undef, 3, 0)
    end
    mat = Matrix{T}(undef, 3, length(points))
    @inbounds for (col, point) in enumerate(points)
        mat[1, col] = point[1]
        mat[2, col] = point[2]
        mat[3, col] = point[3]
    end
    mat
end

"""
    add_point!(tracker, point) -> Int

Append a 3D `point` (any indexable of length 3) to the tracker and return its
1-based index.
"""
function add_point!(tracker::FrontTracker{T}, point) where {T}
    push!(tracker.points, _point_tuple(point, T))
    num_points(tracker)
end

# -- Line helpers ------------------------------------------------------------

function _check_indices(tracker::FrontTracker, idxs)
    npts = num_points(tracker)
    any(i -> i < 1 || i > npts, idxs) && throw(ArgumentError("Point index out of bounds."))
    nothing
end

"""
    add_line!(tracker, line) -> Int

Register a 1D line segment defined by two point indices. Returns the new line
count.
"""
function add_line!(tracker::FrontTracker, line::Line)
    _check_indices(tracker, line)
    push!(tracker.lines, line)
    num_lines(tracker)
end

# -- Surface helpers ---------------------------------------------------------

"""
    add_surface!(tracker, surface) -> Int

Register a 2D triangular surface defined by three point indices. Returns the new
surface count.
"""
function add_surface!(tracker::FrontTracker, surface::Surface)
    _check_indices(tracker, surface)
    push!(tracker.surfaces, surface)
    num_surfaces(tracker)
end

# -- Removal utilities -------------------------------------------------------

"""
    remove_line!(tracker, idx)

Remove the `idx`-th stored line.
"""
function remove_line!(tracker::FrontTracker, idx::Int)
    1 <= idx <= num_lines(tracker) || throw(ArgumentError("Line index out of bounds."))
    deleteat!(tracker.lines, idx)
    tracker
end

"""
    remove_surface!(tracker, idx)

Remove the `idx`-th stored surface.
"""
function remove_surface!(tracker::FrontTracker, idx::Int)
    1 <= idx <= num_surfaces(tracker) || throw(ArgumentError("Surface index out of bounds."))
    deleteat!(tracker.surfaces, idx)
    tracker
end

function _reindex(val::Int, removed::Int)
    val > removed ? val - 1 : val
end

"""
    remove_point!(tracker, idx)

Remove the `idx`-th point and reindex any remaining connectivity.
"""
function remove_point!(tracker::FrontTracker, idx::Int)
    1 <= idx <= num_points(tracker) || throw(ArgumentError("Point index out of bounds."))
    deleteat!(tracker.points, idx)

    tracker.lines = Line[( _reindex(a, idx), _reindex(b, idx) )
                         for (a, b) in tracker.lines if a != idx && b != idx]

    tracker.surfaces = Surface[( _reindex(a, idx), _reindex(b, idx), _reindex(c, idx) )
                               for (a, b, c) in tracker.surfaces if a != idx && b != idx && c != idx]
    tracker
end

# -- Queries -----------------------------------------------------------------

"""
    lines_touching(tracker, point_idx)

Return all line segments that reference `point_idx`.
"""
function lines_touching(tracker::FrontTracker, point_idx::Int)
    _check_indices(tracker, [point_idx])
    [line for line in tracker.lines if point_idx in line]
end

"""
    surfaces_touching(tracker, point_idx)

Return all surfaces that reference `point_idx`.
"""
function surfaces_touching(tracker::FrontTracker, point_idx::Int)
    _check_indices(tracker, [point_idx])
    [face for face in tracker.surfaces if point_idx in face]
end

function Base.show(io::IO, tracker::FrontTracker{T}) where {T}
    print(io, "FrontTracker{$T}(points=$(num_points(tracker)), lines=$(num_lines(tracker)), surfaces=$(num_surfaces(tracker)))")
end

# -- Classical initializers -------------------------------------------------

"""
    point_front(; position=(0.0, 0.0, 0.0), T=Float64, metadata=Dict{Symbol,Any}())

Create a tracker that only stores a single 0D point.
"""
function point_front(; position=(0.0, 0.0, 0.0),
                      T::Type{<:Real}=Float64,
                      metadata=Dict{Symbol,Any}())
    tracker = FrontTracker(T=T, metadata=_metadata_dict(metadata))
    add_point!(tracker, position)
    tracker
end

"""
    line_front(; center=(0.0, 0.0, 0.0), direction=(1.0, 0.0, 0.0),
                length=1.0, segments::Integer=1, T=Float64, metadata=Dict{Symbol,Any}())

Create a polyline aligned with `direction` that spans `length` and is split into
`segments` equal edges.
"""
function line_front(; center=(0.0, 0.0, 0.0),
                    direction=(1.0, 0.0, 0.0),
                    length::Real=1.0,
                    segments::Integer=1,
                    T::Type{<:Real}=Float64,
                    metadata=Dict{Symbol,Any}())
    length > 0 || throw(ArgumentError("Length must be positive."))
    segments >= 1 || throw(ArgumentError("Segments must be >= 1."))

    tracker = FrontTracker(T=T, metadata=_metadata_dict(metadata))
    unit_dir = _normalize_direction(direction)
    center_pt = _point_tuple(center, Float64)
    start = _tuple_add(center_pt, _tuple_scale(unit_dir, -0.5 * length))
    step = length / segments

    prev_idx = nothing
    for i in 0:segments
        pos = _tuple_add(start, _tuple_scale(unit_dir, step * i))
        idx = add_point!(tracker, pos)
        if prev_idx !== nothing
            add_line!(tracker, (prev_idx, idx))
        end
        prev_idx = idx
    end
    tracker
end

"""
    circle_front(; center=(0.0, 0.0, 0.0), radius=1.0, segments::Integer=32,
                  axis::Symbol=:z, T=Float64, metadata=Dict{Symbol,Any}())

Approximate a circular 1D interface lying in the plane orthogonal to `axis`.
`segments` controls the polygon resolution.
"""
function circle_front(; center=(0.0, 0.0, 0.0),
                       radius::Real=1.0,
                       segments::Integer=32,
                       axis::Symbol=:z,
                       T::Type{<:Real}=Float64,
                       metadata=Dict{Symbol,Any}())
    radius > 0 || throw(ArgumentError("Radius must be positive."))
    segments >= 3 || throw(ArgumentError("Segments must be >= 3."))

    tracker = FrontTracker(T=T, metadata=_metadata_dict(metadata))
    basis_u, basis_v = _plane_axes(axis)
    center_pt = _point_tuple(center, Float64)

    first_idx = nothing
    prev_idx = nothing
    for k in 0:(segments - 1)
        angle = 2 * pi * k / segments
        offset = _tuple_add(_tuple_scale(basis_u, radius * cos(angle)),
                            _tuple_scale(basis_v, radius * sin(angle)))
        idx = add_point!(tracker, _tuple_add(center_pt, offset))
        if first_idx === nothing
            first_idx = idx
        end
        if prev_idx !== nothing
            add_line!(tracker, (prev_idx, idx))
        end
        prev_idx = idx
    end
    add_line!(tracker, (prev_idx, first_idx))
    tracker
end

"""
    planar_patch_front(; center=(0.0, 0.0, 0.0), axis=:z,
                        size=(1.0, 1.0), divisions=(1, 1),
                        T=Float64, metadata=Dict{Symbol,Any}())

Create a rectangular grid of points that tiles a plane orthogonal to `axis` with
`divisions` triangles per direction. Both line connectivity and triangular
surfaces are generated.
"""
function planar_patch_front(; center=(0.0, 0.0, 0.0),
                             axis::Symbol=:z,
                             size::Tuple{<:Real,<:Real}=(1.0, 1.0),
                             divisions::Tuple{Int,Int}=(1, 1),
                             T::Type{<:Real}=Float64,
                             metadata=Dict{Symbol,Any}())
    nx, ny = divisions
    nx >= 1 || throw(ArgumentError("divisions[1] must be >= 1."))
    ny >= 1 || throw(ArgumentError("divisions[2] must be >= 1."))
    size_u, size_v = size
    size_u > 0 || throw(ArgumentError("size[1] must be positive."))
    size_v > 0 || throw(ArgumentError("size[2] must be positive."))

    tracker = FrontTracker(T=T, metadata=_metadata_dict(metadata))
    basis_u, basis_v = _plane_axes(axis)
    center_pt = _point_tuple(center, Float64)
    step_u = size_u / nx
    step_v = size_v / ny
    start_u = -size_u / 2
    start_v = -size_v / 2

    grid = Matrix{Int}(undef, nx + 1, ny + 1)
    for j in 0:ny
        for i in 0:nx
            local_u = start_u + step_u * i
            local_v = start_v + step_v * j
            offset = _tuple_add(_tuple_scale(basis_u, local_u), _tuple_scale(basis_v, local_v))
            grid[i + 1, j + 1] = add_point!(tracker, _tuple_add(center_pt, offset))
        end
    end

    for j in 1:(ny + 1)
        for i in 1:nx
            add_line!(tracker, (grid[i, j], grid[i + 1, j]))
        end
    end
    for j in 1:ny
        for i in 1:(nx + 1)
            add_line!(tracker, (grid[i, j], grid[i, j + 1]))
        end
    end

    for j in 1:ny
        for i in 1:nx
            p00 = grid[i, j]
            p10 = grid[i + 1, j]
            p01 = grid[i, j + 1]
            p11 = grid[i + 1, j + 1]
            add_surface!(tracker, (p00, p10, p11))
            add_surface!(tracker, (p00, p11, p01))
        end
    end
    tracker
end

"""
    sphere_front(; center=(0.0, 0.0, 0.0), radius=1.0,
                  rings::Integer=8, segments::Integer=16,
                  T=Float64, metadata=Dict{Symbol,Any}())

Construct a triangulated sphere using a latitude/longitude layout. `rings`
controls the number of interior latitude bands (excluding the poles) and
`segments` controls the azimuthal resolution.
"""
function sphere_front(; center=(0.0, 0.0, 0.0),
                        radius::Real=1.0,
                        rings::Integer=8,
                        segments::Integer=16,
                        T::Type{<:Real}=Float64,
                        metadata=Dict{Symbol,Any}())
    radius > 0 || throw(ArgumentError("Radius must be positive."))
    rings >= 1 || throw(ArgumentError("rings must be >= 1."))
    segments >= 3 || throw(ArgumentError("segments must be >= 3."))

    tracker = FrontTracker(T=T, metadata=_metadata_dict(metadata))
    center_pt = _point_tuple(center, Float64)

    north = add_point!(tracker, _tuple_add(center_pt, (0.0, 0.0, radius)))
    south = add_point!(tracker, _tuple_add(center_pt, (0.0, 0.0, -radius)))

    ring_indices = Matrix{Int}(undef, segments, rings)
    for ring in 1:rings
        θ = pi * ring / (rings + 1)
        sinθ = sin(θ)
        cosθ = cos(θ)
        for seg in 0:(segments - 1)
            φ = 2 * pi * seg / segments
            offset = (radius * sinθ * cos(φ),
                      radius * sinθ * sin(φ),
                      radius * cosθ)
            ring_indices[seg + 1, ring] = add_point!(tracker, _tuple_add(center_pt, offset))
        end
    end

    wrap_next(idx) = idx == segments ? 1 : idx + 1

    for ring in 1:rings
        for seg in 1:segments
            add_line!(tracker, (ring_indices[seg, ring], ring_indices[wrap_next(seg), ring]))
        end
    end

    for seg in 1:segments
        add_line!(tracker, (north, ring_indices[seg, 1]))
        add_line!(tracker, (ring_indices[seg, rings], south))
    end

    for ring in 1:(rings - 1)
        for seg in 1:segments
            add_line!(tracker, (ring_indices[seg, ring], ring_indices[seg, ring + 1]))
        end
    end

    for seg in 1:segments
        next_seg = wrap_next(seg)
        add_surface!(tracker, (north, ring_indices[seg, 1], ring_indices[next_seg, 1]))
        add_surface!(tracker, (ring_indices[seg, rings], ring_indices[next_seg, rings], south))
    end

    for ring in 1:(rings - 1)
        for seg in 1:segments
            next_seg = wrap_next(seg)
            p00 = ring_indices[seg, ring]
            p10 = ring_indices[next_seg, ring]
            p01 = ring_indices[seg, ring + 1]
            p11 = ring_indices[next_seg, ring + 1]
            add_surface!(tracker, (p00, p10, p11))
            add_surface!(tracker, (p00, p11, p01))
        end
    end

    tracker
end

"""
    star_front(; center=(0.0, 0.0, 0.0), radius=1.0,
                spikes::Integer=5, dent::Real=0.4,
                axis::Symbol=:z, filled::Bool=true,
                T=Float64, metadata=Dict{Symbol,Any}())

Create a star-shaped polygon lying on the plane orthogonal to `axis`.
`spikes` controls the number of star tips while `dent` controls the
indentation depth (0 < dent < 1). When `filled=true`, triangular surfaces
fan out from the center.
"""
function star_front(; center=(0.0, 0.0, 0.0),
                     radius::Real=1.0,
                     spikes::Integer=5,
                     dent::Real=0.4,
                     axis::Symbol=:z,
                     filled::Bool=true,
                     T::Type{<:Real}=Float64,
                     metadata=Dict{Symbol,Any}())
    radius > 0 || throw(ArgumentError("Radius must be positive."))
    spikes >= 3 || throw(ArgumentError("spikes must be >= 3."))
    0 < dent < 1 || throw(ArgumentError("dent must lie between 0 and 1."))

    tracker = FrontTracker(T=T, metadata=_metadata_dict(metadata))
    basis_u, basis_v = _plane_axes(axis)
    center_pt = _point_tuple(center, Float64)
    inner_radius = radius * (1 - dent)
    angles = 0:(2 * spikes - 1)

    polygon_indices = Int[]
    prev_idx = nothing
    first_idx = nothing
    for k in angles
        r = iseven(k) ? radius : inner_radius
        angle = pi * k / spikes
        offset = _tuple_add(_tuple_scale(basis_u, r * cos(angle)),
                            _tuple_scale(basis_v, r * sin(angle)))
        idx = add_point!(tracker, _tuple_add(center_pt, offset))
        push!(polygon_indices, idx)
        if first_idx === nothing
            first_idx = idx
        elseif prev_idx !== nothing
            add_line!(tracker, (prev_idx, idx))
        end
        prev_idx = idx
    end
    add_line!(tracker, (prev_idx, first_idx))

    if filled
        center_idx = add_point!(tracker, center_pt)
        for (i, idx) in enumerate(polygon_indices)
            next_idx = polygon_indices[mod1(i + 1, length(polygon_indices))]
            add_surface!(tracker, (center_idx, idx, next_idx))
        end
    end

    tracker
end

"""
    star3d_front(; center=(0.0, 0.0, 0.0), radius=1.0,
                  spike_length=0.25, petals::Integer=6,
                  rings::Integer=6, segments::Integer=24,
                  T=Float64, metadata=Dict{Symbol,Any}())

Generate a smooth, flower-like surface by modulating a spherical radius with
`petals` lobes around the azimuthal direction. The `spike_length` parameter
controls how far the lobes bulge from the base radius (0 ≤ spike_length < 1).
"""
function star3d_front(; center=(0.0, 0.0, 0.0),
                       radius::Real=1.0,
                       spike_length::Real=0.25,
                       petals::Integer=6,
                       rings::Integer=6,
                       segments::Integer=24,
                       T::Type{<:Real}=Float64,
                       metadata=Dict{Symbol,Any}())
    radius > 0 || throw(ArgumentError("Radius must be positive."))
    0 <= spike_length < 1 || throw(ArgumentError("spike_length must satisfy 0 ≤ value < 1."))
    petals >= 3 || throw(ArgumentError("petals must be >= 3."))
    rings >= 1 || throw(ArgumentError("rings must be >= 1."))
    segments >= 3 || throw(ArgumentError("segments must be >= 3."))

    tracker = FrontTracker(T=T, metadata=_metadata_dict(metadata))
    center_pt = _point_tuple(center, Float64)

    north = add_point!(tracker, _tuple_add(center_pt, (0.0, 0.0, radius)))
    south = add_point!(tracker, _tuple_add(center_pt, (0.0, 0.0, -radius)))

    ring_indices = Matrix{Int}(undef, segments, rings)
    for ring in 1:rings
        θ = pi * ring / (rings + 1)
        sinθ = sin(θ)
        cosθ = cos(θ)
        envelope = sinθ^2
        for seg in 0:(segments - 1)
            φ = 2 * pi * seg / segments
            modulation = 1 + spike_length * envelope * cos(petals * φ)
            local_radius = radius * modulation
            offset = (local_radius * sinθ * cos(φ),
                      local_radius * sinθ * sin(φ),
                      local_radius * cosθ)
            ring_indices[seg + 1, ring] = add_point!(tracker, _tuple_add(center_pt, offset))
        end
    end

    wrap_next(idx) = idx == segments ? 1 : idx + 1

    for ring in 1:rings
        for seg in 1:segments
            add_line!(tracker, (ring_indices[seg, ring], ring_indices[wrap_next(seg), ring]))
        end
    end

    for seg in 1:segments
        add_line!(tracker, (north, ring_indices[seg, 1]))
        add_line!(tracker, (ring_indices[seg, rings], south))
    end

    for ring in 1:(rings - 1)
        for seg in 1:segments
            add_line!(tracker, (ring_indices[seg, ring], ring_indices[seg, ring + 1]))
        end
    end

    for seg in 1:segments
        next_seg = wrap_next(seg)
        add_surface!(tracker, (north, ring_indices[seg, 1], ring_indices[next_seg, 1]))
        add_surface!(tracker, (ring_indices[seg, rings], ring_indices[next_seg, rings], south))
    end

    for ring in 1:(rings - 1)
        for seg in 1:segments
            next_seg = wrap_next(seg)
            p00 = ring_indices[seg, ring]
            p10 = ring_indices[next_seg, ring]
            p01 = ring_indices[seg, ring + 1]
            p11 = ring_indices[next_seg, ring + 1]
            add_surface!(tracker, (p00, p10, p11))
            add_surface!(tracker, (p00, p11, p01))
        end
    end

    tracker
end

"""
    write_vtk_front(tracker, filepath; compress=true)

Export the contents of `tracker` as an unstructured VTK `.vtu` file using
WriteVTK. Returns the generated filename, which ParaView can open directly.
"""
function write_vtk_front(tracker::FrontTracker, filepath::AbstractString; compress::Bool=true)
    num_points(tracker) > 0 || throw(ArgumentError("Tracker must contain at least one point to export."))

    points = _points_matrix(tracker.points)
    cells = MeshCell[]

    for idx in 1:num_points(tracker)
        push!(cells, MeshCell(VTKCellTypes.VTK_VERTEX, (idx,)))
    end
    for line in tracker.lines
        push!(cells, MeshCell(VTKCellTypes.VTK_LINE, line))
    end
    for face in tracker.surfaces
        push!(cells, MeshCell(VTKCellTypes.VTK_TRIANGLE, face))
    end

    dataset = vtk_grid(filepath, points, cells; compress=compress)
    vtk_save(dataset)
    dataset.path
end

end # module
