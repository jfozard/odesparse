# Author: Travis Oliphant

__all__ = ['odeints','odeints_cfunc','dlopen', 'dlclose']



_msgs = {2: "Integration successful.",
         -1: "Excess work done on this call (perhaps wrong Dfun type).",
         -2: "Excess accuracy requested (tolerances too small).",
         -3: "Illegal input detected (internal error).",
         -4: "Repeated error test failures (internal error).",
         -5: "Repeated convergence failures (perhaps bad Jacobian or tolerances).",
         -6: "Error weight became zero during problem.",
         -7: "Fatal error in sparse solver"
         }

import _odesparse
from copy import copy
from scipy.sparse import csr_matrix
import numpy as np

def odeints(func, y0, t, args=(), full_output=0,
            rtol=None, atol=None, tcrit=None, h0=0.0,
            hmax=0.0, hmin=0.0, seth=0.0, mxords=5, mxstep=0, mxhnil=0, lrw=0, nnz=0,
            printmesg=0, JPat=None):
    """
    Integrate a system of ordinary differential equations.

    Solve a system of ordinary differential equations using lsodes from the
    FORTRAN library odepack.

    Solves the initial value problem for stiff or non-stiff systems
    of first order ode-s::

        dy/dt = func(y,t0,...)

    where y can be a vector.

    Parameters
    ----------
    func : callable(y, t0, ...)
        Computes the derivative of y at t0.
    y0 : array
        Initial condition on y (can be a vector).
    t : array
        A sequence of time points for which to solve for y.  The initial
        value point should be the first element of this sequence.
    args : tuple
        Extra arguments to pass to function.
    full_output : boolean
        True if to return a dictionary of optional outputs as the second output
    printmessg : boolean
        Whether to print the convergence message

    Returns
    -------
    y : array, shape (len(y0), len(t))
        Array containing the value of y for each desired time in t,
        with the initial value y0 in the first row.

    infodict : dict, only returned if full_output == True
        Dictionary containing additional output information

        =======  ============================================================
        key      meaning
        =======  ============================================================
        'hu'     vector of step sizes successfully used for each time step.
        'tcur'   vector with the value of t reached for each time step.
                 (will always be at least as large as the input times).
        'tolsf'  vector of tolerance scale factors, greater than 1.0,
                 computed when a request for too much accuracy was detected.
        'nst'    cumulative number of time steps
        'nfe'    cumulative number of function evaluations for each time step
        'nje'    cumulative number of jacobian evaluations for each time step
        'nqu'    a vector of method orders for each successful step.
        'imxer'  index of the component of largest magnitude in the
                 weighted local error vector (e / ewt) on an error return, -1
                 otherwise.
        'lenrw'  the length of the double work array required.
        'leniw'  the length of integer work array required.
        'nnz'    the number of nonzero elements in the Jacobian matrix
        'ngp'    the number of groups of column indices, used in 
                 difference quotient Jacobian aproximations.
        'nzl'    the number of nonzero elements in the strict upper
                 triangle of the LU factorization used in the chord iteration
        'nzu'    the number of nonzero elements in the strict upper
                 triangle of the LU factorization used in the chord iteration

        =======  ============================================================

    Other Parameters
    ----------------

    rtol, atol : float
        The input parameters rtol and atol determine the error
        control performed by the solver.  The solver will control the
        vector, e, of estimated local errors in y, according to an
        inequality of the form ``max-norm of (e / ewt) <= 1``,
        where ewt is a vector of positive error weights computed as:
        ``ewt = rtol * abs(y) + atol``
        rtol and atol can be either vectors the same length as y or scalars.
        Defaults to 1.49012e-8.
    tcrit : array
        Vector of critical points (e.g. singularities) where integration
        care should be taken.
    h0 : float, (0: solver-determined)
        The step size to be attempted on the first step.
    hmax : float, (0: solver-determined)
        The maximum absolute step size allowed.
    hmin : float, (0: solver-determined)
        The minimum absolute step size allowed.
    seth : float
        Threshold for sparsity in jacobian
    mxords : integer, (0: solver-determined)
        Maximum order to be allowed for the stiff (BDF) method.
    mxstep : integer, (0: solver-determined)
        Maximum number of (internally defined) steps allowed for each
        integration point in t.
    mxhnil : integer, (0: solver-determined)
        Maximum number of messages printed.
    JPat : sparse array (nodes x nodes) giving sparsity pattern of Jacobian

    """


    t = copy(t)
    y0 = copy(y0)
    
    if JPat is not None:
        J=csr_matrix(JPat)
        ia=J.indptr+1
        ja=J.indices+1
        nnz=J.nnz
    else:
        ia=None
        ja=None


    output = _odesparse.odeints(func, y0, t, args, 
                             full_output, rtol, atol, tcrit, h0, hmax, hmin,
                             seth, mxords, mxstep, mxhnil, lrw, nnz, ia, ja)
    if output[-1] < 0:
        print _msgs[output[-1]]
        print "Run with full_output = 1 to get quantitative information."
    else:
        if printmesg:
            print _msgs[output[-1]]

    if full_output:
        output[1]['message'] = _msgs[output[-1]]

    output = output[:-1]
    if len(output) == 1:
        return output[0]
    else:
        return output


def odeints_cfunc(cfunc, y0, t, full_output=0, neq=None,
            rtol=None, atol=None, tcrit=None, h0=0.0,
            hmax=0.0, hmin=0.0, seth=0.0, mxords=5, mxstep=0, mxhnil=0, lrw=0, nnz=0,
            printmesg=0, JPat=None):
    """
    Integrate a system of ordinary differential equations.

    Solve a system of ordinary differential equations using lsodes from the
    FORTRAN library odepack.

    Solves the initial value problem for stiff or non-stiff systems
    of first order ode-s::

        dy/dt = func(y,t0,...)

    where y can be a vector.

    Parameters
    ----------
    cfunc : callable(y, t0, ...)
        C function which computes the derivative of y at t0.
    y0 : array
        Initial condition on y (can be a vector).
    t : array
        A sequence of time points for which to solve for y.  The initial
        value point should be the first element of this sequence.
    args : tuple
        Extra arguments to pass to function.
    full_output : boolean
        True if to return a dictionary of optional outputs as the second output
    printmessg : boolean
        Whether to print the convergence message

    Returns
    -------
    y : array, shape (len(y0), len(t))
        Array containing the value of y for each desired time in t,
        with the initial value y0 in the first row.

    infodict : dict, only returned if full_output == True
        Dictionary containing additional output information

        =======  ============================================================
        key      meaning
        =======  ============================================================
        'hu'     vector of step sizes successfully used for each time step.
        'tcur'   vector with the value of t reached for each time step.
                 (will always be at least as large as the input times).
        'tolsf'  vector of tolerance scale factors, greater than 1.0,
                 computed when a request for too much accuracy was detected.
        'nst'    cumulative number of time steps
        'nfe'    cumulative number of function evaluations for each time step
        'nje'    cumulative number of jacobian evaluations for each time step
        'nqu'    a vector of method orders for each successful step.
        'imxer'  index of the component of largest magnitude in the
                 weighted local error vector (e / ewt) on an error return, -1
                 otherwise.
        'lenrw'  the length of the double work array required.
        'leniw'  the length of integer work array required.
        'nnz'    the number of nonzero elements in the Jacobian matrix
        'ngp'    the number of groups of column indices, used in 
                 difference quotient Jacobian aproximations.
        'nzl'    the number of nonzero elements in the strict upper
                 triangle of the LU factorization used in the chord iteration
        'nzu'    the number of nonzero elements in the strict upper
                 triangle of the LU factorization used in the chord iteration

        =======  ============================================================

    Other Parameters
    ----------------

    rtol, atol : float
        The input parameters rtol and atol determine the error
        control performed by the solver.  The solver will control the
        vector, e, of estimated local errors in y, according to an
        inequality of the form ``max-norm of (e / ewt) <= 1``,
        where ewt is a vector of positive error weights computed as:
        ``ewt = rtol * abs(y) + atol``
        rtol and atol can be either vectors the same length as y or scalars.
        Defaults to 1.49012e-8.
    tcrit : array
        Vector of critical points (e.g. singularities) where integration
        care should be taken.
    h0 : float, (0: solver-determined)
        The step size to be attempted on the first step.
    hmax : float, (0: solver-determined)
        The maximum absolute step size allowed.
    hmin : float, (0: solver-determined)
        The minimum absolute step size allowed.
    seth : float
        Threshold for sparsity in jacobian
    mxords : integer, (0: solver-determined)
        Maximum order to be allowed for the stiff (BDF) method.
    mxstep : integer, (0: solver-determined)
        Maximum number of (internally defined) steps allowed for each
        integration point in t.
    mxhnil : integer, (0: solver-determined)
        Maximum number of messages printed.
    JPat : sparse array (nodes x nodes) giving sparsity pattern of Jacobian

    """


    t = copy(t)
    y0 = copy(y0)
    
    if JPat is not None:
        J=csr_matrix(JPat)
        ia=J.indptr+1
        ja=J.indices+1
        nnz=J.nnz
    else:
        ia=None
        ja=None

    if neq is None:
        neq = [len(y0)] 

    print 'neq -> : ', neq

    output = _odesparse.odeints_cfunc(cfunc, y0, t, neq,
                             full_output, rtol, atol, tcrit, h0, hmax, hmin,
                             seth, mxords, mxstep, mxhnil, lrw, nnz, ia, ja)


    if output[-1] < 0:
        print _msgs[output[-1]]
        print "Run with full_output = 1 to get quantitative information."
    else:
        if printmesg:
            print _msgs[output[-1]]

    if full_output:
        output[1]['message'] = _msgs[output[-1]]

    output = output[:-1]
    if len(output) == 1:
        return output[0]
    else:
        return output

def dlopen(dllname, fname):
    return _odesparse.dlopen(dllname, fname)

def dlclose():
    _odesparse.dlclose()
