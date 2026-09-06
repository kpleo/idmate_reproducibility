# IDMate numerical comparisons

## Historical window bounds

6/10 distances exceed the raw bound; 0/10 exceed the 0.050 window criterion.
Every positive excess is smaller than that row's recorded trace deviation.
This empirical comparison is not a reconstruction of a trace-corrected bound.

| Candidate | Raw bound | Window distance | Full distance | Trace deviation | Window / raw bound |
| --- | ---: | ---: | ---: | ---: | ---: |
| Si_2 | 6.405e-10 | 1.342e-09 | 5.401e-01 | 2.328e-09 | 2.095 |
| Si_3 | 1.241e-10 | 1.252e-10 | 9.479e-02 | 2.292e-10 | 1.009 |
| Si_5 | 1.138e-09 | 1.307e-10 | 4.754e-01 | 2.193e-10 | 0.115 |
| Si_6 | 2.891e-10 | 1.240e-10 | 3.119e-02 | 2.207e-10 | 0.429 |
| graphene_5 | 5.014e-11 | 8.193e-12 | 1.202e+00 | 9.921e-12 | 0.163 |
| graphene_6 | 3.624e-03 | 1.644e-10 | 4.032e-01 | 1.898e-10 | 4.535e-08 |
| Al_2 | 3.234e-10 | 1.751e-09 | 3.545e-01 | 4.083e-09 | 5.414 |
| Al_3 | 5.245e-10 | 7.423e-10 | 6.820e-03 | 1.752e-09 | 1.415 |
| Al_5 | 3.412e-10 | 1.771e-09 | 3.574e-01 | 4.140e-09 | 5.189 |
| Al_6 | 5.207e-10 | 7.545e-10 | 4.500e-03 | 1.783e-09 | 1.449 |

## Candidate separation

40 accepted, 18 rejected; 8 abstentions excluded.

| Pair set | Concordant | Tied | Total pairs | AUC |
| --- | ---: | ---: | ---: | ---: |
| All | 500 | 0 | 720 | 0.694 |
| Same material | 156 | 0 | 240 | 0.650 |
| Cross material | 344 | 0 | 480 | 0.717 |

| Accepted material | Rejected material | Concordant / total | AUC |
| --- | --- | ---: | ---: |
| Si | Si | 48 / 96 | 0.500 |
| Si | graphene | 80 / 96 | 0.833 |
| Si | Al | 96 / 96 | 1.000 |
| graphene | Si | 0 / 48 | 0.000 |
| graphene | graphene | 12 / 48 | 0.250 |
| graphene | Al | 24 / 48 | 0.500 |
| Al | Si | 48 / 96 | 0.500 |
| Al | graphene | 96 / 96 | 1.000 |
| Al | Al | 96 / 96 | 1.000 |

| Material | Within-material AUC | AUC after deleting this material |
| --- | ---: | ---: |
| Si | 0.500 | 0.792 |
| graphene | 0.250 | 0.750 |
| Al | 1.000 | 0.486 |

Equal-material mean AUC: 0.583.
Candidate-level bootstrap 95% interval: [0.544, 0.831].
The interval does not use material or trajectory clusters. Deletion is a sensitivity
calculation, not held-out prediction. Decision and proposal family are confounded.

| Proposal family | Accept | Reject | Abstain |
| --- | ---: | ---: | ---: |
| low_cutoff | 22 | 0 | 2 |
| target_subspace | 18 | 0 | 6 |
| coarse_k | 0 | 18 | 0 |

## Silicon shell profiles

The three calculations differ in settings and number of derivative columns.

| Configuration | Columns | Share at S=11 | Share at S=24 | Share at S=27 | 90% radius | 99% radius |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| exploratory | 128 | 0.927 | 0.992 | 0.996 | 2.031 | 3.000 |
| selection | 136 | 0.909 | 0.987 | 0.993 | 2.031 | 3.182 |
| frozen | 168 | 0.904 | 0.982 | 0.993 | 2.031 | 3.182 |

Radii are in inverse bohr. Rounded printed profiles do not reproduce full-matrix closure.
Exploratory log-share slope: -1.529 decades per inverse bohr (34 shells with share above 1.000e-08).
Exploratory mean self-shell fractions span 0.141 to 0.590; they do not decrease monotonically.

## H2 reuse conditions

| Condition | Number of attempts |
| --- | ---: |
| Retained-state spectral margin positive | 54 / 60 |
| Occupied-state spectral margin positive | 60 / 60 |
| Full provisional reuse test passed | 29 / 60 |
| Separate first-order comparison within 10% | 60 / 60 |

Retained windows are open at iterations 8-61; reuse passes occur at 33-61.
The 25 open-window failures all have unconverged candidate residuals.
The perturbative comparison is distinct from the provisional reuse test and is not a decision condition.
