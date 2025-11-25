# FrontCutTracking.jl

Minimal utilities for building and manipulating 0D/1D/2D fronts used in front
tracking style algorithms. The `FrontTracker` container stores points, lines,
and triangular surfaces with simple mutation helpers plus a few ready-made
initializers for classical interface setups.

## Quick start

```julia
julia> using FrontCutTracking

julia> tracker = FrontTracker();

julia> a = add_point!(tracker, (0, 0, 0));
julia> b = add_point!(tracker, (1, 0, 0));
julia> add_line!(tracker, (a, b));
1
```

## Classical initializers

```julia
circle = circle_front(radius = 1.0, segments = 32)
line   = line_front(length = 0.1, segments = 5, direction = (0, 1, 0))
patch  = planar_patch_front(size = (1.0, 1.0), divisions = (8, 8))
sphere = sphere_front(radius = 0.5, rings = 6, segments = 24)
star   = star_front(radius = 2.0, spikes = 5, dent = 0.4)
star3d = star3d_front(radius = 1.0, petals = 8, spike_length = 0.25)
```

`star3d_front` modulates a sphere with sinusoidal lobes, yielding a smooth,
flower-like surface; tweak the `petals` count and `spike_length` amplitude to
match the effect you need.

Each initializer returns a fully populated `FrontTracker`, so you can start
advecting or augmenting them immediately (for instance, to tag metadata or add
additional connectivity).

## Exporting to ParaView

Install `WriteVTK` (already listed as a dependency) and call:

```julia
tracker = sphere_front(radius = 0.5)
write_vtk_front(tracker, "front_sphere")  # creates front_sphere.vtu
```

Open the resulting `.vtu` file inside ParaView to inspect the mesh.

## Running the tests

```bash
julia --project -e 'using Pkg; Pkg.test()'
```