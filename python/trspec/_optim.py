"""Nelder-Mead simplex search with the interface of MATLAB's fminsearch."""
import numpy as np
from scipy.optimize import minimize


def fminsearch(fun, x0, tol_x=1e-4, tol_fun=1e-4, max_iter=None, max_fun_evals=None):
    """Minimize fun(x) from x0 with the Nelder-Mead simplex method.

    Same start simplex (5 % steps, 0.00025 for zero entries), coefficients and
    stopping rule (TolX and TolFun both met) as fminsearch.

    Returns (x, fval, exitflag, iterations); exitflag 1 = converged,
    0 = stopped at the iteration or evaluation limit.
    """
    x0 = np.atleast_1d(np.asarray(x0, dtype=float)).ravel()
    n = x0.size
    opts = {'xatol': tol_x, 'fatol': tol_fun,
            'maxiter': 200 * n if max_iter is None else int(max_iter),
            'maxfev': 200 * n if max_fun_evals is None else int(max_fun_evals)}
    res = minimize(lambda x: float(fun(x)), x0, method='Nelder-Mead', options=opts)
    exitflag = 1 if res.status == 0 else 0
    return np.asarray(res.x, dtype=float), float(res.fun), exitflag, int(res.nit)
