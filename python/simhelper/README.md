# LLP event sampling in `gen_llp`

`gen_llp` generates LLP decays inside a fixed rectangular MATHUSLA decay
volume. It combines a supplied set of LLP four-vectors with a separate pool of
decays generated in the LLP rest frame, places the decays inside the detector,
boosts the decay products to the lab frame, and writes the resulting particles
to a ROOT tree.

The sampling is an importance-sampling and rejection-sampling procedure. The
description below follows the implementation exactly; it does not assume that
the input four-vectors are distributed uniformly.

## Inputs and geometry

The principal inputs to the sampling are:

- `filename_llp4vec`: a comma-separated array of LLP four-vectors in the order
  `(E, px, py, pz)`. If the file has five columns, the first is treated as a
  weight column, removed from the kinematics, and the first row's weight is
  saved as `weight0`.
- `filename_products`: a pool of LLP decay records. Each daughter is represented
  by `(E, px, py, pz, PDG ID)` in the LLP rest frame.
- `mX`: the LLP mass, in the same energy units as the input four-vectors.
- `ctau`: the proper decay length, in metres.
- `N_MC`: the number of accepted events to write.
- `rand_seed`: the seed used for both NumPy's random generator and Python's
  `random` module.

The geometry is hard-coded in CMS coordinates, in metres. The decay volume is
39 m wide in `x`, 11 m high in `y`, and 39 m long in `z`. It spans

```text
x = [-19.5, 19.5]
y = [86.47, 97.47]
z = [70, 109]
```

The generated azimuth is restricted to

```text
phi = [-atan(19.5 / 86.47), +atan(19.5 / 86.47)].
```

The polar-angle limits implied by the near-top and far-bottom corners of the
volume are converted to a pseudorapidity interval. In the code, the LLP
pseudorapidity is evaluated as

```text
eta = arctanh(|pz| / |p|).
```

Only four-vectors inside that interval enter the sampling pool. The absolute
value means that the original sign of the longitudinal momentum is discarded.

## Four-vector preprocessing and importance probabilities

For each LLP four-vector that passes the pseudorapidity cut, the function
computes

```text
P       = |p|
P_T     = sqrt(px^2 + py^2)
P_L     = |pz|
betaGamma = P / mX
gamma     = sqrt(mX^2 + P^2) / mX.
```

It then estimates how likely that momentum is to produce a decay in the volume.
Thirty azimuths are drawn uniformly from the restricted `phi` interval. At each
azimuth, the momentum direction is reconstructed as

```text
Px = P_T sin(phi)
Py = P_T cos(phi)
Pz = P_L.
```

The ray from the CMS interaction point along this direction is intersected with
the decay-volume box. If the entry and exit distances from the origin are
`L1` and `L2`, the probability that the LLP decays between them is

```text
P_decay(phi) = exp[-L1 / (betaGamma * ctau)]
             - exp[-L2 / (betaGamma * ctau)].
```

Azimuths that do not intersect the volume contribute zero. The sum is always
divided by 30, including missed intersections. Two quantities are retained for
each four-vector:

- `Pbar_decay`: the 30-draw mean of `P_decay(phi)`;
- `P_decay_max`: the largest `P_decay(phi)` seen in those 30 draws.

The discrete proposal distribution for the LLP four-vector is the normalized
mean probability

```text
q(i) = Pbar_decay(i) / sum_j Pbar_decay(j).
```

Thus, four-vectors estimated to be more likely to decay inside the detector are
proposed more often. These estimates are Monte Carlo estimates made from the
same seeded NumPy generator, rather than exact angular integrals.

## Accepted-event sampling

Candidates are generated until `N_MC` events have been accepted:

1. Choose a four-vector index `i` according to `q(i)`.
2. Draw a new `phi` uniformly in the restricted azimuthal interval and construct
   `(Px, Py, P_L)` as above.
3. Intersect its ray with the volume. If there is no nonzero path through the
   box, reject the candidate and restart.
4. Recompute `P_decay(phi)` from the resulting `L1` and `L2`.
5. Draw the decay distance from an exponential distribution with mean
   `betaGamma * ctau`, truncated to `[L1, L2]`.
6. Accept the candidate with probability

   ```text
   P_decay(phi) / P_decay_max(i).
   ```

   Otherwise, reject it and restart. Here `P_decay_max(i)` is only the maximum
   observed during the earlier 30-draw scan, not an analytic upper bound.
7. Place the accepted vertex at the sampled distance along the ray.

In pseudocode, the core procedure is:

```text
for each LLP four-vector passing the eta cut:
    estimate Pbar_decay and P_decay_max from 30 random phi values

while accepted_events < N_MC:
    i <- sample with probability proportional to Pbar_decay[i]
    phi <- uniform(phi_min, phi_max)
    L1, L2 <- ray-box intersection
    if there is no path through the box: continue

    P_decay <- exp(-L1 / (betaGamma*ctau))
             - exp(-L2 / (betaGamma*ctau))
    L <- truncated_exponential(betaGamma*ctau, L1, L2)
    if uniform(0, 1) > P_decay / P_decay_max[i]: continue

    construct vertex at L
    attach and boost one rest-frame decay
    write event
```

## Decay products, weights, and output

Before event generation, the decay-record indices are randomly permuted. The
`n`th accepted event uses entry `n mod N_decay` from that permutation, so each
decay record is used once per cycle in a randomized order.

Every daughter four-vector is Lorentz-boosted from the LLP rest frame along the
accepted LLP flight direction, using the selected LLP's `gamma`. The decay
vertex and daughter momenta are then converted from CMS coordinates to the
Geant coordinate convention. Positions are converted from metres to
millimetres, and momenta are converted from GeV to MeV.

The output ROOT tree is named `data` and contains vector branches
`Gen_pdgID`, `Gen_x`, `Gen_y`, `Gen_z`, `Gen_t`, `Gen_px`, `Gen_py`, and
`Gen_pz`. All daughters in one event share the sampled vertex, and `Gen_t` is
set to zero.

For every accepted event, the function computes

```text
weight_total = weight0
             * (N_LLP / N_MC)
             * ((phi_max - phi_min) / (2*pi))
             * P_decay_max(i).
```

In the current implementation this value is not written to a ROOT branch or
returned. The function returns `vecs_used_list`, the indices into the
post-pseudorapidity-cut four-vector list for all accepted events.

With a fixed `rand_seed` and unchanged inputs, both the NumPy draws and the
randomized decay ordering are reproducible. `metadata` is accepted as an
argument but is not currently used.

## Summary

LLP four-momenta were sampled from a pre-generated distribution after restricting the pseudorapidity to trajectories compatible with the detector acceptance. For each retained four-momentum, the probability of decay inside the detector was estimated by averaging over 30 azimuthal angles drawn uniformly within the geometrically allowed range. At each angle, the entry and exit distances, $L_1$ and $L_2$, were obtained from the intersection of the LLP trajectory with the detector volume, and the in-volume decay probability was evaluated as $P_{\mathrm{dec}}=\exp[-L_1/(\beta\gamma c\tau)]-\exp[-L_2/(\beta\gamma c\tau)]$, where $\beta\gamma=|\mathbf{p}|/m_X$. LLP four-momenta were then selected with probabilities proportional to their azimuthally averaged decay probabilities, thereby preferentially sampling kinematic configurations with greater detector acceptance.

For each selected four-momentum, a new azimuthal angle was generated and the decay position was drawn from the exponential decay-length distribution, truncated to the segment of the trajectory contained within the detector. The candidate was retained using rejection sampling with an acceptance factor given by the ratio of its in-volume decay probability to the maximum probability found in the corresponding 30-angle scan. A decay final state was drawn from a randomly ordered library of rest-frame decays, Lorentz-boosted along the LLP flight direction, and transformed to the detector coordinate and unit conventions. The nominal event normalization was calculated as the product of the input-sample weight, the ratio of available LLP four-momenta to generated events, the sampled azimuthal fraction, and the maximum in-volume decay probability; in the current implementation, this value is calculated but is not stored in the output tree.
