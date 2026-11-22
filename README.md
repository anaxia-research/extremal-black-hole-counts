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

## Parts II and III

Part II (N = 1 super-extremal partition functions) and Part III (the Maloney-Witten-Keller sign map) will be added here.
