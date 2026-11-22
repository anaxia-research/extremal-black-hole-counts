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

## Further parts

Part III (Ramond primaries) and Part IV (the Maloney-Witten partition function) will be added here.
