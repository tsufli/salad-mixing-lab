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

**Scope:** A periodic cube with deformable finite-width triangulated strip surfaces, no lentils, gravity, adhesive links or fluid back-reaction. Each strip uses a fast 5 × 2 vertex grid by default, with 5 × 3 and 7 × 3 mesh comparison options. Total strip mass is held constant across these mesh resolutions, but a mesh-convergence study of the elastic discretization is still required. Springs on longitudinal, transverse and diagonal mesh edges penalize in-plane stretch/shear. Longitudinal/transverse second differences and mixed differences approximate bending and twist. These are *discrete mesh energies*, not an exact thin-shell constitutive law.

**Drive:** Three successive incompressible sine shears, x depending on y, y depending on z, z depending on x. One phase lasts τ=2.5. Control parameters are strip number and size, periodic box length (density), bending penalty, friction and shear displacement K.

**Contact:** Pairwise distances between triangle surfaces use vertex-face, edge-edge and edge-through-face intersection tests with periodic images. Normal penalty contact and history-dependent Coulomb-limited tangential springs distribute forces to triangle vertices using barycentric weights. Nonlocal self-contact is included. A spatial hash and per-triangle bounding-box filter accelerate candidate discovery. Temporary tangential memory is tracked per triangle pair; continuous inter-strip contact lifetime is tracked per strip pair. The code does *not* implement continuous collision detection or exact nonpenetration; numerical crossing remains possible.

**Measurements:** Largest connected component fraction; number and mean age of inter-strip contact links; mean duration of ended contact links; and a velocity coherence diagnostic after subtracting translation and best-fit angular velocity. Coherence is undefined for clusters of fewer than three strips. Neither connectivity nor coherence by itself proves true topological entanglement or material rigidity.

**Count-preserving packing:** Requested strips are initialized in a larger dilute periodic cube and then compressed to the selected target box. At each preparation step the center of every strip is translated according to a smooth cube-contraction schedule, without scaling the strip's internal vertices or its reference spring lengths. During packing, background stirring is off and the strips relax under elastic and contact forces. This is a prescribed geometrical compaction, not a pressure-controlled barostat. The simulation displays progress through **three stages**: smooth compression, an optional unforced relaxation at constant box size, and stirred dynamics. Stirring time starts at t=0 after relaxation. It reports the largest cluster at the end of compression and again after unforced relaxation, so transient compression contacts can be distinguished from persistent connectivity. Packing and relaxation duration are separately configurable; fast compression may cause excessive penetration.

**Performance:** The physics runs in a Web Worker; the main thread renders only when the simulation delivers a frame or the camera changes. The triangle broad phase now indexes periodic, thickness-inflated triangle axis-aligned bounding boxes in a spatial hash, and uses cheap separating-plane rejection before full proximity evaluation. The default mesh is 5 × 2; 5 × 3 and 7 × 3 remain selectable. The time step was increased from 0.002 to 0.003 after basic numerical checks. In a same-timestep sparse benchmark the new contact algorithm was approximately 2.6× faster for 22–40 strips. Accuracy and performance at dense packing can differ markedly.

The numerical simulation is still exploratory and computationally intensive at high density. It uses soft contacts, no continuous collision detection, and **does not guarantee no tunneling or topologically faithful entanglement**. Monitor the maximum overlap during compression and stirring. Parameter/time-step/mesh convergence and independent validation are required before interpreting a possible transition.

Run all repository tests:

```bash
node --test tests/*.test.mjs
```

## Dense-packing experiments (November 2026)

The 3D tab now supports **12–100 strips** and a **target effective volume fraction** control (3%–28%), plus manual cube side length and loose/dense/very-dense presets. The estimated density
`phi_eff = N * strip_length * strip_width * contact_thickness / box_side^3`
is a computational control parameter, not a measurement of real solid salad volume. Choosing phi automatically computes the target box side. The default setup is 64 strips, length 3.5, width 1.2, phi_eff 0.18, a 5×2 coarse mesh, 420 packing steps and 240 relaxation steps.

The plots include a **compression trajectory: largest instantaneous contact cluster fraction versus effective phi**, with an additional marker for the cluster surviving unforced relaxation. This trajectory is protocol- and history-dependent and should **not be mistaken for an equilibrium percolation curve**.

Normal contact stiffness/damping were increased to `kn=3200, etaN=35` and drag to `gamma=9` to reduce spurious overlap during high-density compression. A 64-strip, 18%-density, balanced 420-step compression test showed peak penetration of about 24% of effective contact thickness, whereas the earlier weak-contact settings yielded roughly 57% with the same setup. These parameters are not experimentally calibrated. A high-density run can still experience very large penetration, especially under rapid compression or extreme density; the website shows the **maximum observed overlap divided by contact thickness** and warns above 25%.

The browser worker reports solver steps/second. Large dense packings remain expensive; for higher quality or scientific claims check time step, mesh and contact stiffness convergence, and use continuous collision detection or a collision-safe shell method before interpreting topological locking.

## Optional lentils and unified leaf coloring (November 2026)

The 3D tab has a **Lentils** checkbox (off by default) and an optional 20–200 particle count control (default 110). After changing the checkbox, use **New run** or **Restart same seed** to rebuild the simulation. With the checkbox off, no lentil particles are created and the strips-only experiment is preserved. All strip meshes are rendered in the same leafy-green palette, shaded only by orientation; lentils appear as brown spheres.

Lentils are actual **finite-mass 3D spheres**, currently a geometrically simplified proxy rather than lentil-shaped disks, with default radius 0.115 and mass 0.22. They are initialized without overlap in the dilute box, participate in packing, zero-drive relaxation and stirring, and feel the same prescribed sine-flow drag as strip vertices. Contacts between lentils, and between lentils and strip triangles, use the same normal spring/dashpot and Coulomb-limited tangential contact law as strip contacts. Sphere–surface forces are distributed barycentrically to triangle vertices with equal and opposite reaction on the sphere. Their tangential contact memory clears when contact breaks. Short-range sphere–strip candidates use a periodic spatial hash. The leaf-only connectivity statistic still connects strips to other strips, not via an intervening lentil.

This remains an exploratory penalty-contact numerical model. Especially in dense preparations, sphere–surface crossing and overlap are possible because the simulation does not have continuous collision detection. The lentil count is deliberately independently tunable rather than constrained to a calibrated salad composition.

## Constant-volume mixture comparisons (November 2026)

The 3D tab now includes an opt-in **Fixed total ingredient volume** comparison mode (enabled by default). The lentils checkbox still controls whether the second species is present. When it is present, a new **lentils per leaf** ratio slider sets `r = N_lentils / N_leaves` over 0–4. The existing leaf count slider becomes the reference leaf-only count `N0`. The nominal ingredient-volume budget is `V0 = N0 * leafLength * leafWidth * thickness`. A model leaf contributes `vLeaf = length * width * 0.16` and a lentil sphere contributes `vLentil = 4*pi*radius^3/3`.

For each requested ratio, the model selects integer populations near `N_leaves = V0/(vLeaf+r*vLentil)`, `N_lentils = r*N_leaves`. Because species must be integer counts, the realized total volume differs slightly from the target; **the UI displays the actual counts, ratio and relative volume error**. With target effective fraction selected, `boxSide = cbrt(V0/targetPhi)`, independent of mixture ratio, so both target total material volume and periodic box volume remain constant as the mixture changes. The actual effective fraction includes *both* species. The manual-box mode can also be used without changing box side during a ratio scan.

For comparisons, keep the reference leaf count, individual ingredient sizes, friction, bending, stirring, compression protocol, and random seed fixed. Use multiple seeds for ensemble means and uncertainties. The CSV export records the complete composition and parameters preceding the time-series columns. Manual mixture mode remains available by disabling Fixed total ingredient volume.

**Important interpretation:** This geometric volume is a simple modeling surrogate, not experimentally calibrated mass, excluded volume or the volume of salad as prepared in a bowl. With the current model sizes one lentil occupies about 0.95% of one leaf's effective volume, so a number ratio of 2 replaces very little leaf volume. A future comparison may prefer the *lentil share of occupied volume* as a second composition variable.
