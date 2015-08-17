# The goal is to implement my own lattice sampler.

from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler


def _fpbkz(A, K = 10**30, block = 8, delta = 0.75):
    """
    including a transpose operation.
    """
    print 'blocksize for bkz = %s'%block
    At = A.transpose()
    RF = A[0][0].parent()
    AA = Matrix(ZZ, [[round(K*a) for a in row] for row in list(At)])
    F = FP_LLL(AA)
    F.BKZ(block_size = block, delta= delta)
    B = F._sage_()
    T = B*AA**(-1)
    B1 = Matrix(RF, [[a/RF(K) for a in row] for row in list(B)])
    return T.transpose(), B1.transpose()


class MyLatticeSampler:
    """
    Sampling from discrete Gaussian.
    """

    def __init__(self,A,sigma = 1,dps = 60, method = 'LLL', already_orthogonal = False, gram_schmidt_norms = None):
        self.A = A # we are using column span instead of rowspan
        self.sigma = sigma

        self._degree = A.nrows()
        if method == 'LLL':
            self.T =  self._lll_reduce()
        elif method == 'BKZ':
            self.T = self._bkz_reduce()
        else:
            raise NotImplementedError
        self.B  = self.A*self.T


        if already_orthogonal: # The columns of A are already gram-schmidt.
            self._G = self.A
            if gram_schmidt_norms is None:
                self._gs_norms = [self._G.column(i).norm() for i in range(self._degree)]
            else:
                self._gs_norms = gram_schmidt_norms
        else:
            # Compute the gram-schmidt ourselves. Can be slow.
            self._gs_norms, self._G = self.compute_G(dps =  dps)

        self.final_sigma = sigma*(prod(self._gs_norms))**(1/self._degree)

        self.Ainv = (self.A)**(-1)

    def _bkz_reduce(self,block = None):
        print 'bkz being performed...'
        if block is None:
            block = ZZ(self._degree // 2)
        return _fpbkz(self.A, block = block)[0]

    def _lll_reduce(self):
        print 'lll being performed...'
        A = self.A
        return gp(A).qflll().sage()

    def col_sum(self):
        """
        related to the evaluation attack, return the list a where
                   a[i] = colsum(A^-1,i)
        """
        return vector([1 for _ in range(self._degree)])*self.Ainv
    
    def babai_quality(self):
        """
        inspired by Kim's explanation, I think the quality of a basis
        for babai should be the ratio ||\tilde{bn}||/||\tilde{b1}||
        """
        gs_norms = self._gs_norms
        return float(min(gs_norms)/max(gs_norms))

    def __repr__(self):
        return 'Discrete Gaussian sampler with dimension %s and sigma = %s'%(self._degree, self.final_sigma.n())

    def compute_G(self, dps = 50):
        t = cputime()
        B = self.B
        n = self._degree
        from mpmath import *
        mp.dps = dps
        # print 'dps = %s'%dps
        prec = dps*6
        AA = mp.matrix([list(w) for w in list(B)])
        Q,R = qr(AA) # QR decomposition

        M = mp.matrix([list(Q.column(i)*R[i,i]) for i in range(n)]);
        M_sage = Matrix([[RealField(prec)(M[i,j]) for i in range(n)] for j in range(n)])
        verbose('gram schmidt computation took %s'%cputime(t))
        return [abs(RealField(prec)(R[i,i])) for i in range(n)], M_sage

    def set_sigma(self,newsigma):
        self.final_sigma = newsigma

    def babai(self,c):
        """
        run babai's algorithm and find a lattice vector close to the
        input point c.
        Note this is super similar to the __call__ function

        Returns a tuple (v,z), where v is the actual vector in R^n,
        and z is its coordinate *in terms of a*. So we have
        v = Az.
        """
        # also need matrix. Okay.
        n = self._degree
        try:
            c = vector(c)
        except:
            pass
        G, norms  = self._G, self._gs_norms
        B = self.B
        T = self.T
        zs = []
        v = c

        for i in range(n)[::-1]:
            b_ = G.column(i)
            v_ = v.dot_product(b_) / norms[i]**2
            # sigma_ = sigma/b_.norm()
            z = ZZ(round(v_))
            v = v - z*B.column(i)
            zs.append(z)
        return c - v, T*(vector(zs[::-1]))

    def __call__(self, c = None):
        """
        c -- an n-dimensional vector, so that we are sampling a discrete gaussian
        centered at c.
        """
        v = 0
        # also need matrix. Okay.
        sigma, G = self.final_sigma, self._G
        n = self._degree
        if c is None:
            c = zero_vector(n)
        B = self.B
        T = self.T
        zs = []
        norms = self._gs_norms
        # print 'stds = %s'%stds
        for i in range(n)[::-1]:
            b_ = G.column(i)
            c_ = c.dot_product(b_) / norms[i]**2
            # sigma_ = sigma/b_.norm()
            sigma_ = sigma/norms[i]
            #print 'sigma = %s'%sigma_
            #print 'c_ = %s'%c_
            assert(sigma_ > 0)
            z = DiscreteGaussianDistributionIntegerSampler(sigma=sigma_, c=c_, algorithm="uniform+table")()
            c = c - z*B.column(i)
            v = v + z*B.column(i)
            zs.append(z)
        return v, T*vector(zs[::-1])

