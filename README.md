# Salad Mixing Lab

Interactive exploratory granular DEM model of rigid rods (leaves) mixed with small disks (lentils) under **alternating sine flow** on a 2D periodic torus.

**Website:** https://tsufli.github.io/salad-mixing-lab/

## Stirring protocol

Over successive half-periods of length tau=2.5, a prescribed velocity field acts on all particles:

- First: u_x = (W K / tau) sin(2 pi y / H), u_y = 0
- Next: u_x = 0, u_y = (H K / tau) sin(2 pi x / W)

For passive tracers the **exact** one-period map in normalized coordinates X=x/W, Y=y/H is

X' = [ X + K sin(2 pi Y) ] mod 1

Y' = [ Y + K sin(2 pi X') ] mod 1.

This area-preserving map is an analytical reference, **not** the exact discrete particle dynamics.

## Discrete model

Disks and spherocylinders interact via excluded-volume normal springs with normal damping and smooth regularized dynamic Coulomb friction. Rod–disk forces generate torques. All pair contacts use nearest periodic images, and particles wrap through periodic boundaries. Finite mass, velocity relaxation toward the flow, and rotational drag are included. **No physical walls, fluid solver, static friction history, gravity or flexibility.**

The simulation initializes disks preferentially in the left half, while rods are distributed throughout the periodic box. Mixing is reported as the spatial variance of disk counts over a fixed 4x3 grid, normalized by their mean squared. For finite samples uniform random placement has a nonzero noise floor.

## Run and test

Open index.html using a local static HTTP server:

```bash
python3 -m http.server 8000
node --test tests/*.test.mjs
```

GitHub Actions runs tests and deploys Pages.
