from mpmath.libmpf import fzero, fone, finf, fninf, fnan
from mpmath.libmpf import bitcount, to_fixed, from_man_exp
from mpmath.libmpf import MP_ZERO, MP_ONE, MP_BASE
from mpmath.gammazeta import mpf_zeta_int, mpf_euler, euler_fixed

gamma_taylor_cache = {}

from mpmath import *
from mpmath.libmpf import *
from mpmath.libelefun import *
from mpmath.gammazeta import *
from mpmath.libintmath import int_fac

def zeta_array(N, prec):
    """
    zeta(n) = A * pi**n / n! + B

    where A is a rational number (A = Bernoulli number
    for n even) and B is an infinite sum over powers of exp(2*pi).
    (B = 0 for n even).
    """
    extra = 30
    wp = prec+extra
    zeta_values = [MP_ZERO] * (N+2)
    pi = pi_fixed(wp)
    # STEP 1:
    one = MP_ONE << wp
    zeta_values[0] = -one//2
    f_2pi = mpf_shift(mpf_pi(wp),1)
    exp_2pi_k = exp_2pi = mpf_exp(f_2pi, wp)
    # Compute exponential series
    # Store values of 1/(exp(2*pi*k)-1),
    # exp(2*pi*k)/(exp(2*pi*k)-1)**2, 1/(exp(2*pi*k)-1)**2
    # pi*k*exp(2*pi*k)/(exp(2*pi*k)-1)**2
    exps3 = []
    k = 1
    while 1:
        tp = wp - 9*k
        if tp < 1:
            break
        # 1/(exp(2*pi*k-1)
        q1 = mpf_div(fone, mpf_sub(exp_2pi_k, fone, tp), tp)
        # pi*k*exp(2*pi*k)/(exp(2*pi*k)-1)**2
        q2 = mpf_mul(exp_2pi_k, mpf_mul(q1,q1,tp), tp)
        q1 = to_fixed(q1, wp)
        q2 = to_fixed(q2, wp)
        q2 = (k * q2 * pi) >> wp
        exps3.append((q1, q2))
        # Multiply for next round
        exp_2pi_k = mpf_mul(exp_2pi_k, exp_2pi, wp)
        k += 1
    # Exponential sum
    for n in xrange(3, N+1, 2):
        s = MP_ZERO
        k = 1
        for e1, e2 in exps3:
            if n%4 == 3:
                t = e1 // k**n
            else:
                U = (n-1)//4
                t = (e1 + e2//U) // k**n
            if not t:
                break
            s += t
            k += 1
        zeta_values[n] = -2*s
    # Even zeta values
    B = [mpf_abs(mpf_bernoulli(k,wp)) for k in xrange(N+2)]
    pi_pow = fpi = mpf_pow_int(mpf_shift(mpf_pi(wp), 1), 2, wp)
    pi_pow = mpf_div(pi_pow, from_int(4), wp)
    for n in xrange(2,N+2,2):
        z = mpf_mul(B[n], pi_pow, wp)
        zeta_values[n] = to_fixed(z, wp)
        pi_pow = mpf_mul(pi_pow, fpi, wp)
        pi_pow = mpf_div(pi_pow, from_int((n+1)*(n+2)), wp)
    # Zeta sum
    reciprocal_pi = (one << wp) // pi
    for n in xrange(3, N+1, 4):
        U = (n-3)//4
        s = zeta_values[4*U+4]*(4*U+7)//4
        for k in xrange(1, U+1):
            s -= (zeta_values[4*k] * zeta_values[4*U+4-4*k]) >> wp
        zeta_values[n] += (2*s*reciprocal_pi) >> wp
    for n in xrange(5, N+1, 4):
        U = (n-1)//4
        s = zeta_values[4*U+2]*(2*U+1)
        for k in xrange(1, 2*U+1):
            s += ((-1)**k*2*k* zeta_values[2*k] * zeta_values[4*U+2-2*k])>>wp
        zeta_values[n] += ((s*reciprocal_pi)>>wp)//(2*U)
    #return [mpf(from_man_exp(x,-wp,prec,'n')) for x in zeta_values]
    return [x>>extra for x in zeta_values]

def gamma_taylor_coefficients(prec):
    """
    Gives the Taylor coefficients of 1/gamma(x) around x = 1 as
    a list of fixed-point numbers. Enough coefficients are returned
    to ensure that the series converges to the given precision
    when x is in [0.5, 1.5].
    """
    if prec in gamma_taylor_cache:
        return gamma_taylor_cache[prec]
    # The estimate is valid up to a precision of at least 15000 bits
    if prec < 1000:
        N = int(prec**0.76 + 2)
    else:
        N = int(prec**0.787 + 2)
    wp = prec + 20
    A = [0] * N
    A[0] = MP_ZERO
    A[1] = MP_ONE << wp
    A[2] = euler_fixed(wp)
    #zeta_values = [0,0]+[to_fixed(mpf_zeta_int(k,wp),wp) for k in xrange(2,N)]
    zeta_values = zeta_array(N, wp)
    for k in xrange(3, N):
        a = (-A[2]*A[k-1])>>wp
        for j in xrange(2,k):
            a += ((-1)**j * zeta_values[j] * A[k-j]) >> wp
        a //= (1-k)
        A[k] = a
    A = [a>>20 for a in A]
    A = A[::-1]
    gamma_taylor_cache[prec] = A
    return A
