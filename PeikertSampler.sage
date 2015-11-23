# T = t_matrix(n,prec = prec)
# discrete Gaussian rounding.
from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler


def pei_t_matrix(s, prec = 100):
    """
    return the unitary matrix 1/sqrt{2} (1,1,-i, i) \otimes I.
    """
    C = ComplexField(prec)
    eyer = Matrix.identity(s)

    return C(1/sqrt(2))*block_matrix([[eyer, eyer], [-C(I)*eyer, C(I)*eyer]])

def _real_part(mat):
    """
    real part of a matrix.
    """
    return Matrix([[matentry.real_part() for matentry in matrixrow] for matrixrow in list(mat)])


def _index_matrix(p):
    Top = Matrix(ZZ, [[Mod(i*j,p) for i in range(p-1)] for j in range(1, (p+1)//2)])
    Bottom = Matrix(ZZ, [[Mod(i*(p-j),p) for i in range(p-1)] for j in range(1, (p+1)//2)])
    return Top.stack(Bottom)


class PeikertSampler:

    def __init__(self, p, s, r = None, prec = 100):
        """
        r -- a small rounding parameter  1 <= r ~ \sqrt{log(p-1)}.
        s -- the discrete Gaussian width normalized by 1/sqrt(p). We require
        s > r.
        """
        if r is None:
            r = RR(sqrt(log(p-1)))
        if s <= r:
            raise ValueError('s must be large enough')
        self.p = p
        self.r = r
        self.s = s*RR(sqrt(p))
        self.prec = prec
        self.n = p-1
        self.C = ComplexField(self.prec)


        t  = cputime()
        self.Bcomplex = self._Bcomplex()
        print 'Bc. %s'%cputime(t)
        self.T = pei_t_matrix((self.p-1)//2 , self.prec)
        self.B1 = _real_part(self.T*self.Bcomplex)
        print 'B1. %s'%cputime(t)

        self.Sigma1 = self.B1*self.B1.transpose()
        self.B1inv = self.compute_B1inv()
        print 'B1inv %s'%cputime(t)


        self.B2 = self.compute_B2()
        print 'B2 %s'%cputime(t)

        # Precompute the matrix B1^-1*B2.
        self.B1invB2 = self.B1inv*self.B2
        print 'B1invB2 %s'%cputime(t)

        # continuous Gaussian
        self.Dcont = RealDistribution('gaussian', 1/RR(sqrt(2*pi)))


    def __repr__(self):
        pass


    def _Bcomplex(self):
        p = self.p
        C = self.C
        zetap = C.zeta(p)
        Index = _index_matrix(p)
        return Matrix(C, [[zetap**Index[i][j] for j in range(p-1)] for i in range(p-1)])

    def compute_B1inv(self):
        """
        compute the inverse of B1 without directly calling.
        """
        n = self.n
        Bc = self.Bcomplex
        Bcstar = Bc.conjugate_transpose()
        C = self.C
        ones = ones_matrix(C, n)
        eye = identity_matrix(C, n)
        Bcinv =  C(1/self.p)*(ones+eye)*Bcstar
        Tstar = self.T.conjugate_transpose()
        return _real_part(Bcinv*Tstar)

    def compute_B2(self):
        s = self.s
        r = self.r
        prec = self.prec
        n = self.n
        Sigma2 =  s**2*identity_matrix(self.n) - r**2*self.Sigma1
        from mpmath import *
        mp.dps = prec // 2
        Sigma2_mp = mp.matrix([list(ww) for ww in list(Sigma2)])
        B2 = cholesky(Sigma2_mp)
        B2_sage = Matrix([[RealField(prec)(B2[i,j]) for j in range(n)] for i in range(n)])
        return B2_sage


    def __call__(self):
        # first we sample a multivariate Gaussian of width one.
        v = vector([self.Dcont.get_random_element() for _ in range(self.n)])
        x2 = self.B1invB2*v
        return _randomized_rounding(x2, self.r)


def _randomized_rounding(v,r):
    """
    v -- a vector in RR^n.
    r >0  -- rounding width.
    return a vector in ZZ^n close to v.
    """
    result = []
    sigma = RR(r/sqrt(2*pi))
    R = v[0].parent()
    for i in range(len(v)):
        vprime = R(v[i] - floor(v[i]))
        D = DiscreteGaussianDistributionIntegerSampler(sigma = sigma, c = vprime)
        result.append( ZZ(floor(v[i]) + D()))
    return result

