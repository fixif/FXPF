# coding=utf8

"""
This file contains tests for the FXPF wrapper
"""


__author__ = "Anastasia Volkova"
__copyright__ = ""
__credits__ = ["Anastasia Volkova"]

__license__ = "GPL v3"
__version__ = "1.0"
__maintainer__ = "Anastasia Volkova"
__email__ = "anastasia.volkova@inria.fr"
__status__ = "Beta"


import pytest
import numpy

from numpy import array, zeros, absolute, eye, dot, diagflat, ones, fliplr, atleast_2d, r_
from numpy import matrix as mat
from numpy.testing import assert_allclose
from numpy.random.mtrand import rand, randn, randint, choice, random_sample
from numpy.core.umath import pi, cos, sin
from numpy.linalg import solve
from numpy.linalg.linalg import LinAlgError
from scipy.signal import butter


from fixif.FXPF import FXPF_ABCD, FXPF_readme


def random_ABCD(n, p, q, pRepeat=0.01, pReal=0.5, pBCmask=0.90, pDmask=0.8, pDzero=0.5):
	"""
	Generate ONE n-th order random  stable state-spaces, with q inputs and p outputs
	copy/adapted from control-python library (Richard Murray): https://sourceforge.net/projects/python-control/
	(thanks guys!)
	possibly already adpated/copied from Mathworks or Octave

	Parameters:
	- n: number of states (default:  random between 5 and 10)
	- p: number of outputs (default: 1)
	- q: number of inputs (default: 1)

	- pRepeat: Probability of repeating a previous root (default: 0.01)
	- pReal: Probability of choosing a real root (default: 0.5). Note that when choosing a complex root,
		the conjugate gets chosen as well. So the expected proportion of real roots is pReal / (pReal + 2 * (1 - pReal))
	- pBCmask: Probability that an element in B or C will not be masked out (default: 0.90)
	- pDmask: Probability that an element in D will not be masked out (default: 0.8)
	- pDzero: Probability that D = 0 (default: 0.5)

	Returns a four numpy matrices A,B,C,D
	"""

	# Make some poles for A.  Preallocate a complex array.
	poles = zeros(n) + zeros(n) * 0.j
	i = 0

	while i < n:

		if rand() < pRepeat and i != 0 and i != n - 1:
			# Small chance of copying poles, if we're not at the first or last  element.
			if poles[i - 1].imag == 0:
				poles[i] = poles[i - 1]     # Copy previous real pole.
				i += 1

			else:
				poles[i:i + 2] = poles[i - 2:i]     # Copy previous complex conjugate pair of poles.
				i += 2

		elif rand() < pReal or i == n - 1:
			poles[i] = 2. * rand() - 1.     # No-oscillation pole.
			i += 1

		else:
			mag = rand()    # Complex conjugate pair of oscillating poles.
			phase = 2. * pi * rand()
			poles[i] = complex(mag * cos(phase), mag * sin(phase))
			poles[i + 1] = complex(poles[i].real, -poles[i].imag)
			i += 2

	# Now put the poles in A as real blocks on the diagonal.

	A = zeros((n, n))
	i = 0

	while i < n:

		if poles[i].imag == 0:
			A[i, i] = poles[i].real
			i += 1

		else:
			A[i, i] = A[i + 1, i + 1] = poles[i].real
			A[i, i + 1] = poles[i].imag
			A[i + 1, i] = -poles[i].imag
			i += 2


	while True:     # Finally, apply a transformation so that A is not block-diagonal.
		T = randn(n, n)

		try:
			A = dot(solve(T, A), T)  # A = T \ A * T
			break

		except LinAlgError:
			# In the unlikely event that T is rank-deficient, iterate again.
			pass

	# Make the remaining matrices.
	B = randn(n, q)
	C = randn(p, n)
	D = randn(p, q)

	# Make masks to zero out some of the elements.
	while True:
		Bmask = rand(n, q) < pBCmask
		if not Bmask.all():  # Retry if we get all zeros.
			break

	while True:
		Cmask = rand(p, n) < pBCmask
		if not Cmask.all():  # Retry if we get all zeros.
			break

	if rand() < pDzero:
		Dmask = zeros((p, q))
	else:
		while True:
			Dmask = rand(p, q) < pDmask
			if not Dmask.all():  # Retry if we get all zeros.
				break


	# Apply masks.
	B *= Bmask
	C *= Cmask
	# D *= Dmask

	return A, B, C, D


def iter_random_ABCD(number, n=(5, 10), p=(1, 5), q=(1, 5), pRepeat=0.01, pReal=0.5, pBCmask=0.90, pDmask=0.8, pDzero=0.5):
	"""
	Generate some n-th order random (stable) state-spaces, with q inputs and p outputs
	Returns:
		- returns a generator of numpy matrices (A,B,C,D)  (to use in a for loop for example)
	"""
	for _ in range(number):
		yield random_ABCD(randint(*n), randint(*p), randint(*q), pRepeat, pReal, pBCmask, pDmask, pDzero)


def test_FXPF_simple():
	wl = numpy.array([12, 12, 12, 12])

	A = numpy.array([[2.076913426318702e+00,-1.534298334880014e+00, 3.908833248284007e-01],
					 [1.0, 0, 0],[0, 1.0, 0]])
	B = numpy.array([[1.0], [0], [0]])
	C = numpy.array([[4.220284791563418e-02, 1.218393525130549e-02,1.156199298609820e-02]])
	D = numpy.array([[8.312697966613880e-03]])
	u_bound = numpy.array([1.0])
	msb, lsb, error, additionalSteps = FXPF_ABCD(A, B, C, D, u_bound, wl)

	print("\n FXPF are determined using %s additional steps." % (mat(additionalSteps)))
	print("The system is:")
	print (mat(A), B, C, D)
	print ("Wordlengths: %s" % (wl))
	print ("MSBs: %s" % (msb))
	print("LSBs: %s" % (lsb))
	print("Errors: %s" % (error))


@pytest.mark.parametrize("S", iter_random_ABCD(1, (5, 10), (1, 2), (1, 2), pBCmask=0.1))
def test_FXPF(S):
	A,B,C,D = S
	order = A.shape[0]
	outputs, inputs = D.shape
	wl = 16 * ones([order + outputs, 1])

	u_bound = numpy.array([1.0])
	msb, lsb, error, additionalSteps = FXPF_ABCD(A, B, C, D, u_bound, wl)
	print("\n FXPF are determined using %s additional steps." % mat(additionalSteps))
	print("The system is: \n A=%s \n B= %s \n C=%s \n D= %s \n" %(A, B, C, D))
	print ("Wordlengths: %s" % (wl))
	print ("MSBs: %s" % (msb))
	print("LSBs: %s" % (lsb))
	print("Errors: %s" % (error))



@pytest.mark.parametrize("S", iter_random_ABCD(1, (3, 10), (1, 2), (1, 2), pBCmask=0.1))
def test_FXPF(S):
	A,B,C,D = S
	order = A.shape[0]
	outputs, inputs = D.shape
	wl = 16 * ones([order + outputs, 1])

	u_bound = numpy.array([1.0])
	msb, lsb, error, additionalSteps = FXPF_ABCD(A, B, C, D, u_bound, wl)
	print("\n FXPF are determined using %s additional steps." % (mat(additionalSteps)))
	print("The system is: \n A=%s \n B= %s \n C=%s \n D= %s \n" %(A, B, C, D))
	print("Wordlengths: %s" % (wl))
	print("MSBs: %s" % (msb))
	print("Errors: %s" % (error))
	print("LSBs: %s" % (lsb))

def test_FXPF_readme():
	FXPF_readme()