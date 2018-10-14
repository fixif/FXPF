# coding: utf8

"""
Python wrapper for the Fixed-Point Formats library

Volkova, A.; Hilaire, T.; Lauter, C., “Determining Fixed-Point Formats for a Digital Filter Implementation
 using the Worst-Case Peak Gain Measure,” in 49th Asilomar Conference on Signals, Systems and Computers,
 pages 737-741, 8-11 November 2015

"""

__author__ = "Anastasia Volkova"
__copyright__ = ""
__credits__ = ["Anastasia Volkova"]

__license__ = "GPL v3"
__version__ = "1.0"
__maintainer__ = "Anastasia Volkova"
__email__ = "anastasia.volkova@inria.fr"
__status__ = "Beta"

import ctypes
import ctypes.util
from numpy import empty, float64, uint64, int, mat

if not ctypes.util.find_library('fxpf'):
    raise ValueError("The FXPF library cannot be found (is it installed?")
_FXPFlib = ctypes.CDLL(ctypes.util.find_library('fxpf'))

"""int FXPF                                    (uint64_t *msb, int *lsb, int *additionalSteps, double *error, int *wl,
											double *A, double *B, double *C, double *D,
											int n, int p, int q, double *u_bound);
"""
_FXPFfun = _FXPFlib.FXPF
_FXPFfun.argtypes = ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(
    ctypes.c_int), ctypes.POINTER(ctypes.c_double),
ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(
    ctypes.c_double),
ctypes.POINTER(ctypes.c_double), ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_double)


def FXPF_ABCD(A, B, C, D, u_bound, wl):
    """Compute the WCPG from the matrices A, B, C, D
		A,B,C and D are numpy matrices or array of the right size"""
    # get the sizes
    n = A.shape[1]
    p, q = D.shape
    # get the pointer to the double arrays
    pA = A.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    pB = B.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    pC = C.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    pD = D.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    pwl = wl.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64))
    pu_bound = u_bound.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    # run the function to fill the empty array W
    msb = empty((p + n), dtype=uint64)
    lsb = empty((p + n), dtype=ctypes.c_int)
    additionalSteps = empty((1, 1), dtype=ctypes.c_int)
    error = empty((p + n), dtype=ctypes.c_double)

    ret = _FXPFfun(msb.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
                   lsb.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                   additionalSteps.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
                   error.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                   pwl, pA, pB, pC, pD, n, p, q, pu_bound)
    if ret < 0:
        raise ValueError("Something went wrong during the FXPF computation.")
    else:
        return mat(msb), mat(lsb), mat(error), mat(additionalSteps)


def FXPF_readme():
    print("Python wrapper for the Fixed-Point Formats library")
