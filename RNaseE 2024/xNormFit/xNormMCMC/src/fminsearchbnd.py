import numpy as np

from scipy.optimize import fmin


def xtransform(x, params):
    xtrans = np.zeros(params['xsize'])
    k = 0
    for i in range(params['n']):
        if params['BoundClass'][i] == 1:
            xtrans[i] = params['LB'][i] + x[k] ** 2
            k += 1
        elif params['BoundClass'][i] == 2:
            xtrans[i] = params['UB'][i] - x[k] ** 2
            k += 1
        elif params['BoundClass'][i] == 3:
            xtrans[i] = (np.sin(x[k]) + 1) / 2
            xtrans[i] = xtrans[i] * (params['UB'][i] - params['LB'][i]) + params['LB'][i]
            xtrans[i] = max(params['LB'][i], min(params['UB'][i], xtrans[i]))
            k += 1
        elif params['BoundClass'][i] == 4:
            xtrans[i] = params['LB'][i]
        elif params['BoundClass'][i] == 0:
            xtrans[i] = x[k]
            k += 1
    return xtrans


def fminsearchbnd(fun, x0, LB, UB, options=None, *args):
    n = len(x0)
    x0 = x0.flatten()

    if LB is None:
        LB = -np.inf * np.ones(n)
    else:
        LB = LB.flatten()

    if UB is None:
        UB = np.inf * np.ones(n)
    else:
        UB = UB.flatten()

    if (n != len(LB)) or (n != len(UB)):
        raise ValueError('x0 is incompatible in size with either LB or UB.')

    params = {'args': args, 'LB': LB, 'UB': UB, 'fun': fun, 'n': n, 'xsize': x0.shape, 'BoundClass': np.zeros(n)}

    for i in range(n):
        k = int(np.isfinite(LB[i])) + 2 * int(np.isfinite(UB[i]))
        params['BoundClass'][i] = k
        if (k == 3) and (LB[i] == UB[i]):
            params['BoundClass'][i] = 4

    x0u = x0.copy()
    k = 0
    for i in range(n):
        if params['BoundClass'][i] == 1:
            if x0[i] <= LB[i]:
                x0u[k] = 0
            else:
                x0u[k] = np.sqrt(x0[i] - LB[i])
            k += 1
        elif params['BoundClass'][i] == 2:
            if x0[i] >= UB[i]:
                x0u[k] = 0
            else:
                x0u[k] = np.sqrt(UB[i] - x0[i])
            k += 1
        elif params['BoundClass'][i] == 3:
            if x0[i] <= LB[i]:
                x0u[k] = -np.pi / 2
            elif x0[i] >= UB[i]:
                x0u[k] = np.pi / 2
            else:
                x0u[k] = 2 * (x0[i] - LB[i]) / (UB[i] - LB[i]) - 1
                x0u[k] = 2 * np.pi + np.arcsin(np.clip(x0u[k], -1, 1))
            k += 1
        elif params['BoundClass'][i] == 0:
            x0u[k] = x0[i]
            k += 1

    if k <= n:
        x0u = x0u[:k]

    if not options:
        options = {}

    # Use a nested function as the OutputFcn wrapper
    def outfun_wrapper(x, *args):
        xtrans = xtransform(x, params)
        if 'outputfcn' in params:
            return params['outputfcn'](xtrans, *args)
        return False

    # Check for an outputfcn. If there is any, then substitute my
    # own wrapper function.
    if 'outputfcn' in options:
        params['outputfcn'] = options['outputfcn']
        options['outputfcn'] = outfun_wrapper

    def intrafun(x):
        xtrans = xtransform(x, params)
        fval = fun(xtrans, *params['args'])
        return fval

    xu = fmin(intrafun, x0u, xtol=1e-4, ftol=1e-4, maxfun=10000, maxiter=10000)

    x = xtransform(xu, params)
    x = np.reshape(x, params['xsize'])

    if 'outputfcn' in params:
        stop = params['outputfcn'](x, params['args'])
        if stop:
            exitflag = 1
        else:
            exitflag = 0
    else:
        exitflag = 0

    output = {
        'iterations': 10000,
        'funcCount': 10000,
        'algorithm': 'fminsearch',
        'message': 'Optimization terminated successfully.'
    }

    return x, intrafun(xu), exitflag, output
