# Salad Mixing Lab

A browser-based dimensionless 2D granular simulation of rods ("leaves") and disks ("lentils"). See the website and source code for model details.

## Running locally
`python3 -m http.server 8000`, then open http://localhost:8000.

## Testing
`node --test tests/*.test.mjs`

## Physics
Rigid capsules and disks interact through repulsive normal spring-dashpots plus regularized dynamic Coulomb friction. Rod–disk forces generate torques. A periodically reversing rotational drive creates motion. This is a preliminary toy model, **not** a validated leaf/lentil simulation. Static friction memory, flexible sheets, gravity and hydrodynamic forces are not yet included.
