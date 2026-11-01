# Salad Mixing Lab — Flexible Fiber Network

Interactive 2D periodic simulation of small disks and **flexible fibers**, driven by the alternating sine flow.

**Live website:** https://tsufli.github.io/salad-mixing-lab/

**3D flexible strips:** https://tsufli.github.io/salad-mixing-lab/strips3d.html

The website now has two tabs: the original 2D flexible-fiber / disk experiment and a new 3D strip-only model.

## Physical model

Each of the default 30 fibers is a five-node bead chain (total 150 nodes) with nearest-neighbor stretching stiffness `ks=400` and bending energy

`E_b = (bending/2) * sum_i |r_(i-1) - 2*r_i + r_(i+1)|^2`.

There are 180 small independent disks by default; the box dimensions are 18 × 12 with periodic boundaries. The stretching/bending energy is not calibrated to real sheet elasticity.

All node/disk pairs on different fibers (or involving disks) interact via repulsive normal spring-dashpot contacts and a **history-dependent tangential spring** limited by Coulomb sliding friction. The tangential spring state is kept only while a particular pair remains in contact. There are no permanent bonds between different fibers. Rigid-body contact moments and continuous rod-segment exclusion are replaced by node-level contacts, a deliberate approximation.

Particles obey

`m_a r_a'' = -d(E_s+E_b)/dr_a + sum_b F_ab(contact) + gamma*m_a*(u(r_a,t)-v_a)`

where the elastic term is absent for disks and the imposed alternating sine field is

- first half period: `u=(W*K/tau*sin(2*pi*y/H),0)`
- second half period: `u=(0,H*K/tau*sin(2*pi*x/W))`

with `tau=2.5` and `gamma=8`. Contact parameters are `kn=650, etaN=10, kt=130, etaT=2`. Numerical integration uses explicit symplectic Euler with `dt=.002` and capped particle velocities.

The exact passive-tracer map is a useful comparison but **not** the interacting particle solution.

## Measurements

- Disk count variance on a 4×3 spatial grid.
- Fraction of fibers in the largest *instantaneous* connected component of the interfiber contact graph.
- Mean contact degree and current uninterrupted contact ages, with a statistic for completed contact durations.
- Radially averaged disk static structure factor `S_dd(k)`, averaged after t >= 10.
- CSV export, new random seeds, and same-seed restarts.

Long-lived contact connections are not permanent bonds. Connectivity is not equivalent to elastic rigidity. The contact duration measurements are node-pair measures, not directly fiber–fiber lifetime survival probabilities.

## Run and tests

```bash
python3 -m http.server 8000
node --test tests/*.test.mjs
```

The GitHub Pages workflow publishes on pushes to main.

## Scope and limitations

This is a **2D flexible bead-chain** model rather than a 3D triangulated flexible sheet model, so it cannot model overlapping leaves or genuine three-dimensional topological interlocking. It uses an imposed flow instead of coupled fluid dynamics. It has not been physically calibrated, validated against experimental mixing rates, or checked systematically for time-step and size convergence. Its numerical stability does not establish predictive realism.


## 3D triangulated-strip experiment

Open `strips3d.html` in the published site. The code is in `strips3d.mjs` (3D physics) and `strips3d-worker.mjs` (background simulation worker). The original 2D model is retained independently.

**Scope:** A periodic cube with deformable finite-width triangulated strip surfaces, no lentils, gravity, adhesive links or fluid back-reaction. Each strip uses a fast 5 × 3 vertex grid by default, with an optional detailed 7 × 3 mesh. Total strip mass is held constant across these mesh resolutions, but a mesh-convergence study of the elastic discretization is still required. Springs on longitudinal, transverse and diagonal mesh edges penalize in-plane stretch/shear. Longitudinal/transverse second differences and mixed differences approximate bending and twist. These are *discrete mesh energies*, not an exact thin-shell constitutive law.

**Drive:** Three successive incompressible sine shears, x depending on y, y depending on z, z depending on x. One phase lasts τ=2.5. Control parameters are strip number and size, periodic box length (density), bending penalty, friction and shear displacement K.

**Contact:** Pairwise distances between triangle surfaces use vertex-face, edge-edge and edge-through-face intersection tests with periodic images. Normal penalty contact and history-dependent Coulomb-limited tangential springs distribute forces to triangle vertices using barycentric weights. Nonlocal self-contact is included. A spatial hash and per-triangle bounding-box filter accelerate candidate discovery. Temporary tangential memory is tracked per triangle pair; continuous inter-strip contact lifetime is tracked per strip pair. The code does *not* implement continuous collision detection or exact nonpenetration; numerical crossing remains possible.

**Measurements:** Largest connected component fraction; number and mean age of inter-strip contact links; mean duration of ended contact links; and a velocity coherence diagnostic after subtracting translation and best-fit angular velocity. Coherence is undefined for clusters of fewer than three strips. Neither connectivity nor coherence by itself proves true topological entanglement or material rigidity.

**Count-preserving packing:** Requested strips are initialized in a larger dilute periodic cube and then compressed to the selected target box. At each preparation step the center of every strip is translated according to a smooth cube-contraction schedule, without scaling the strip's internal vertices or its reference spring lengths. During packing, background stirring is off and the strips relax under elastic and contact forces. This is a prescribed geometrical compaction, not a pressure-controlled barostat. The simulation explicitly displays packing progress and starts stirring time at t=0 when compaction completes. A dropdown allows quick, balanced and gentle packing; fast compression may induce more penetration.

**Performance:** The physics runs in a Web Worker; the main thread renders only when the simulation delivers a frame or the camera changes. The triangle broad phase now indexes periodic, thickness-inflated triangle axis-aligned bounding boxes in a spatial hash, and uses cheap separating-plane rejection before full proximity evaluation. The default mesh is 5 × 3; 7 × 3 remains selectable. The time step was increased from 0.002 to 0.003 after basic numerical checks. In a same-timestep sparse benchmark the new contact algorithm was approximately 2.6× faster for 22–40 strips. Accuracy and performance at dense packing can differ markedly.

The numerical simulation is still exploratory and computationally intensive at high density. It uses soft contacts, no continuous collision detection, and **does not guarantee no tunneling or topologically faithful entanglement**. Monitor the maximum overlap during compression and stirring. Parameter/time-step/mesh convergence and independent validation are required before interpreting a possible transition.

Run all repository tests:

```bash
node --test tests/*.test.mjs
```
