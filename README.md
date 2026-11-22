# Extremal black hole counts

Code and computation logs for the series *Black Hole Counts in Three-Dimensional Gravity*.

## Part I: Positivity for every extremal partition function

Witten's extremal partition function `Z_k` (central charge `c = 24k`) is the unique modular function
equal to the Virasoro vacuum character up to `O(q)`. Its Virasoro primary counts `p_k(n)` (primaries of
weight `k + n`) are the proposed black hole microstate counts. Part I proves that every `p_k(n)` is a positive
integer, increasing in `n` and in `k`, with a Cardy formula for heavy states.

The proof reduces the infinite problem to finite tasks. This repository carries them out.

| File | What it does | Log |
|---|---|---|
| `src/jcoeffs.py` | exact coefficients of `J = j - 744` from `E4^3/Delta` | |
| `src/extremal.py` | `w_k(n)` by the Hecke formula, `p_k(n)` exactly; independent check via polynomials in `J` | |
| `src/check_finite.py` | finite checks F1-F3; positivity, monotonicity, log-concavity for 85,559 values (`k <= 40`, `kn <= 20000`) | `logs/finite_checks.log` |
| `src/bounds.py` | explicit Rademacher-Weil constants used in the analytic estimates | |
| `src/check_light.py` | light regime constants (`n <= 9k`, `k >= 10`) | `logs/light_regime.log` |
| `src/check_heavy.py` | heavy regime error envelope (`n >= 8k`) and monotonicity checks | `logs/heavy_regime.log` |
| `src/tables.py`, `src/appendix.py` | data for the tables and figures | `logs/tables.log`, `logs/appendix.log` |

### Running

Requires Python 3 with `python-flint` and `mpmath`.

```
cd src
python3 check_finite.py
python3 check_light.py
python3 check_heavy.py
python3 tables.py
python3 appendix.py
```

All integer computations are exact. Real constants are evaluated with 40 to 400 significant digits; every
inequality holds with a margin far larger than the rounding error.

## Part II: Supersymmetric extremal partition functions

Witten's N = 1 extremal NS partition function `Phi_k` (central charge `c = 12k`) is invariant under
`Gamma_0(4)+` and equals the NS vacuum character up to `O(q)`. Its superconformal primary counts `N_k(j)`
(weight `(k+j)/2`) are the proposed black hole microstate counts, and `beta_k = (-1)^k Phi_k(1/2)` is
Witten's Ramond index. Part II proves that every `N_k(j)` is positive and increasing in `j` and in `k`,
a Cardy formula with `c_eff = c - 3/2`, and that `beta_k > 0` for odd `k` and `beta_k < 0` for even `k >= 8`.

All files are in `part2/`.

| File | What it does | Log |
|---|---|---|
| `src/superext.py` | exact `kappa(N)` (Hauptmodul of `Gamma_0(4)+`), `d(s)`, Hecke formula, `N_k(j)`, `beta_k` | |
| `src/jcoeffs.py`, `src/bounds2.py` | coefficients of `J`; Bessel main terms and Rademacher-Weil constants | |
| `src/consts2.py` | certified constants `delta_N`, `eta_M` for the coefficients of the Hauptmodul | `logs/constants.log` |
| `src/check_finite.py` | replication formula check, Witten's values, finite checks F1-F4; 128,344 values (`k <= 40`, `kj <= 30000`) | `logs/finite_checks.log` |
| `src/check_light2.py` | light regime constants (`j <= 9k`, `k >= 10`) | `logs/light_regime.log` |
| `src/check_heavy2.py` | a priori bounds, error envelope and monotonicity in the heavy regime | `logs/heavy_regime.log` |
| `src/check_beta.py` | Ramond index: weight 3/2 circle method, exact signs for `k <= 3000` | `logs/beta.log` |
| `src/tables2.py` | data for the tables and figures, Ramond partition functions | `logs/tables.log` |

### Running

```
cd part2
python3 src/consts2.py
python3 src/check_finite.py
python3 src/check_light2.py
python3 src/check_heavy2.py
python3 src/check_beta.py
python3 src/tables2.py
```

## Part III: Ramond primaries

The Ramond sector of Witten's N = 1 extremal theory at `c = 12k` is `H_k = (-1)^k f_k(g_R)`, with
`g_R = -24 - 4096 (eta(2tau)/eta(tau))^24` the Hauptmodul at the Ramond cusp. Part III proves an exact winding formula
`R_k(n) = (-1)^m beta_k [n = m^2] + 1/2 sum_r 4096^r rho_r(n-r) e_r(k-r)` for the number of Ramond primary multiplets,
with positive coefficients `rho_r` and `e_r`. Consequently every `R_k(n)` is a positive integer, strictly increasing
in `n` and `k`, divisible by 4096 away from perfect squares, and obeys a Cardy formula with explicit error for `n >= k`.

All files are in `part3/`. Run everything with `sh part3/run_all.sh`; it takes under a minute.

| File | What it does | Log |
|---|---|---|
| `src/ramond.py` | exact `H_k`, `E_r`, `rho_r` and the winding formula for `R_k(n)` | |
| `src/superext.py`, `src/jcoeffs.py` | Hauptmodul of `Gamma_0(4)+`, NS counts, `beta_k`; coefficients of `J` | |
| `src/taylor_check.py` | Taylor identity for `f_k(-24-u)` (`k <= 15`) | `logs/taylor_check.log` |
| `src/stest.py` | the S-transformation of the Ramond function, to 40 digits | `logs/stest.log` |
| `src/lemma_check.py` | key inequalities K1, K2 exactly for `k <= 3000`; tail bounds | `logs/lemma_check.log` |
| `src/cong.py`, `src/v2check.py` | congruences mod 4096 and 8192; 2-adic data | `logs/cong.log`, `logs/v2check.log` |
| `src/nscompare.py` | `R_k(1) > N_k(1)` for `k <= 3000` | `logs/nscompare.log` |
| `src/consts.py`, `src/apriori.py` | a priori bound on the Hauptmodul; certified constants | `logs/consts.log`, `logs/apriori.log` |
| `src/cardy.py`, `src/envelope.py` | circle method errors, error envelope `E+` and its corner values | `logs/cardy.log`, `logs/envelope.log` |
| `src/unif.py`, `src/lc.py` | data for the uniform Cardy and log-concavity conjectures | `logs/unif.log`, `logs/lc.log` |
| `src/tables3.py`, `src/figdata.py` | data for the tables and figures | `logs/tables3.log`, `logs/figdata.log` |

## Part IV: Extremal N = 2 elliptic genera

An extremal elliptic genus of index `m` (central charge `c = 6m`) is a weak Jacobi form of weight 0 whose polar
coefficients are those of the N = 2 vacuum character and its spectral flows (Gaberdiel-Gukov-Keller-Moore-Ooguri).
Part IV proves that one exists exactly for `m in {1, 2, 3, 4, 5, 7, 8, 11, 13}`, and is then unique.
The proof evaluates the Neveu-Schwarz function at the Ramond cusp of the theta group, which gives one necessary
condition for each residue of `m` mod 4:

* odd `m`: `2 C(2m-1) - rho(2m-1) = 3m - 1`, a lattice-point count in an ellipse, impossible for `m > 3655`;
* `m = 2 mod 4`: `L(m) = -12m`, false for `m >= 6`;
* `m = 0 mod 4`: `X(m) = 0`, false for `m >= 16`.

The even cases use exact series to `m = 40000` and an explicit circle method beyond. Exact linear algebra
(rank certificates mod the prime 33554393) excludes the remaining levels 12, 19, 31, 41 and produces the nine genera.

All files are in `part4/`. Run everything with `sh part4/run_all.sh`.

| File | What it does | Log |
|---|---|---|
| `src/jacobi.py` | two-variable series for `phi_{-2,1}`, `phi_{0,1}`, `E4`, `E6` and the flowed vacuum `W_m` | |
| `src/constraints.py` | the three constraints from theta series in `q^(1/8)`; lattice identity check | `logs/constraints.log` |
| `src/zero4.py` | the `m = 0 mod 4` constraint derived directly from the vacuum character | `logs/zero4.log` |
| `src/sseries.py` | exact `s`-series of `G`, `T`, `D`; `L(m) < -12m` and `X(m) < 0` for all `m <= 40000` | `logs/sseries.log` |
| `src/analytic.py` | certified constants and error/main ratios of the circle method for `N >= 20000` | `logs/analytic.log` |
| `src/transf.py`, `src/identities.py` | transformation laws, theta identities, cusp values of the Hauptmoduln | `logs/transf.log`, `logs/identities.log` |
| `src/polarlemma.py` | the NS polar structure, coefficient by coefficient, for `m <= 60` | `logs/polarlemma.log` |
| `src/cuspcheck.py` | Ramond-cusp evaluations on the actual extremal genera | `logs/cuspcheck.log` |
| `src/polar.py` | rank certificates mod `p` for `m <= 70` | `logs/polar.log` |
| `src/exist.py` | exact rational solutions for the nine levels; no solution for 12 and 19 | `logs/exist.log` |
| `src/tables4.py`, `src/genus_rows.py` | data for the tables and figures; `q^1`, `q^2` rows of the nine genera | `logs/tables4.log`, `logs/genus_rows.log` |
| `src/integral4.py`, `src/signs.py` | integrality and sign pattern of the nine genera through `q^6` | `logs/integral.log`, `logs/signs.log` |

## Further parts

Further parts of the programme will be added here.
