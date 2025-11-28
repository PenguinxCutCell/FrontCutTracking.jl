# FrontCutTracking.jl

[CI](https://github.com/PenguinxCutCell/FrontCutTracking.jl/actions/workflows/CI.yml)

FrontCutTracking.jl provides lightweight 1D/2D front-tracking helpers in pure Julia.
It exposes `FrontTracker` (2D) and `FrontTracker1D` containers along with geometry,
capacity, and space–time utilities used in cut-cell and interface-capturing codes.

## Installing

```julia
julia> ]
pkg> dev https://github.com/PenguinxCutCell/FrontCutTracking.jl
```

The project ships with a `Project.toml`, so activating the repo directory is
enough for running tests or experimenting with the API.

## Quick Start (2D)

```julia
using FrontCutTracking

front = FrontTracker()
create_circle!(front, 0.5, 0.5, 0.25, 120) # populate markers

nodes = ([0.0:0.05:1.0;], [0.0:0.05:1.0;])
capacities = compute_capacities(nodes, front)

segments = compute_segment_parameters(front)
segment_lengths = segments[4]
perimeter = sum(segment_lengths) # discrete perimeter via segment lengths
```

Key helpers:
- `create_circle!`, `create_rectangle!`, `create_ellipse!`, `create_crystal!`
- `fluid_cell_properties`, `compute_surface_capacities`, `compute_capacities`
- `compute_spacetime_capacities`, `compute_segment_cell_intersections`
- `compute_intercept_jacobian`, `update_front_with_intercept_displacements!`

## Quick Start (1D)

```julia
front = FrontTracker1D([0.3, 0.7])
mesh = ([0.0:0.1:1.0;],)

caps = compute_capacities_1d(mesh, front)
st_caps = compute_spacetime_capacities_1d(mesh, front, front, 0.01)
```

## Project Layout

- `src/fronttracker/`: 2D front tracker split into `types`, `markers`, `shapes`,
  `geometry`, `capacities`, `spacetime`, `segments`, and `intercepts`.
- `src/fronttracker1d/`: 1D equivalents (`types`, `geometry`, `capacities`, `spacetime`).
- `test/`: unit tests covering FrontTracker, FrontTracker1D, and analytic geometry checks.

## Running Tests

```bash
julia --project -e "using Pkg; Pkg.test()"
```

