

def _split_primes(m, H, min_prime = 2, max_prime = 10):
    Zm = Integers(m)
    return [p for p in primes(min_prime, max_prime) if Zm(p) in [Zm(h) for h in H]]


from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler


class SubCycSampler:
    """
    We write our own GPV sampler for sub-cyclotomic fields.
    It also has the functionality of simulating an attack.

    Caution: according to GPV, we need to have s >= ||tilde(B)||*log(n)
    for the sampler to approximate discrete lattice Gaussian. So if
    s is smaller than what is required, the __call__() method is not
    guaranteed to output discrete Gaussian.
    """

    def __init__(self,m,H,sigma = 1,prec = 100,  disc = None):
        """
        require: m must be square free and odd.

        disc: the discriminant of K =  Q(zeta_m)^H. We pass it
        as an optional parameter, since when the order of H is
        large, the computation could be very slow.
        """

        self.m = m
        self.H = H
        self.H1 =  H.compute_elements()
        self.cosets = H.cosets()
        self.sigma = sigma
        self.prec = prec

        self._degree = euler_phi(m) // len(self.H1)

        self._is_totally_real = H.is_totally_real()


        if not self._is_totally_real:
            merged_cosets = []
            for c in self.cosets:
                if not any([-c/d in self.H1 for d in merged_cosets]):
                    merged_cosets.append(ZZ(c))
            newcosets = merged_cosets + [-a for a in merged_cosets]
            self.cosets = newcosets

        self.TstarA, self.A = self.embedding_matrix(prec = self.prec)

        self.disc = (self.TstarA).det() #maybe this is faster in computing discriminants.

        self.adj = RR(abs(self.disc)**(1.0/(self._degree)))

        self.final_sigma = self.sigma*self.adj
        self._T = self.lll_transform_matrix()
        self.Ared = self.TstarA*self._T
        # gram-schmidt basis and norms.
        self._G, self.gs_norms = self.compute_G()
        # self._dd_gen = DiscreteGaussianDistributionIntegerSampler(sigma = sigma)

        # # self.D = [ DiscreteGaussianDistributionIntegerSampler(sigma=ss) for ss in self.stds]


    def __repr__(self):
        return 'RLWE error sampler with m = %s,  H = %s and final sigma = %s'%(self.m, self.H, self.final_sigma)

    def minpoly(self):
        K.<z> = CyclotomicField(self.m)
        return sum([z**h for h in self.H1]).minpoly()

    def compute_G(self, prec = 53):
        """
        computing a colum gram-schmidt basis for the embedded lattice O_K.
        return the basis and the length of each vector as a list.

        Modified on 8/2: do this after using LLL to reduce the basis.
        """
        B = self.Ared
        n = self._degree
        from mpmath import *
        mp.dps = prec // 2
        BB = mp.matrix([list(ww) for ww in list(B)])
        Q,R = qr(BB) # QR decomposition
        M = mp.matrix([list(Q.column(i)*R[i,i]) for i in range(n)]);
        M_sage = Matrix([[RealField(prec)(M[i,j]) for i in range(n)] for j in range(n)])
        v = [abs(R[i,i]) for i in range(n)]
        return M_sage,v # vectors are columns

    def basis_lengths(self):
        return [self.Ared.column(i).norm() for i in range(self._degree)]

    def embedding_matrix(self, prec = None):
        """
        to-do: need to separate totally real case from totally complex
        """
        m = self.m
        H1 = self.H1
        if prec is None:
            prec = self.prec
        C = ComplexField(prec)
        zetam = C.zeta(m)
        cosets = self.cosets
        n = self._degree

        if self._is_totally_real:
            A = _real_part(Matrix([[sum([zetam**(ZZ(l*k*h)) for h in H1]) for l in cosets] for k in cosets]))
            return A,A
        else:
            T = t_matrix(n,prec = prec)
            A = Matrix([[sum([zetam**(ZZ(l*k*h)) for h in H1]) for l in cosets] for k in cosets])
            # print A0.det(), T.det()
            return _real_part(T.conjugate_transpose()*A),A

    def lll_transform_matrix(self):
        A = self.TstarA
        return gp(A).qflll().sage()

    def coset_reps(self):
        """
        I need this for representing the basis vectors. Each coset rep c
        represents the element \alpha_c =  \sum_{h \in H} \zeta_m^{ch}.
        """
        return self.cosets

    def _call_dd(self, basis_only = True):
        """
        sample from DD. Most naive sampling
        """
        D = self._dd_gen
        return [D() for _ in range(self._degree)]


    def __call__(self,  method = 'GPV'):
        """
        return an integer vector a = (a_c) indexed by the coset reps of self,
        which represents the vector \sum_c a_c \alpha_c
        Use the algorithm of [GPV].
        http://www.cc.gatech.edu/~cpeikert/pubs/trap_lattice.pdf

        If minkowski = True, return the lattice vector in R^n. Otherwise,
        return the coordinate of the vector in terms of the embedding matrix of self.
        """
        if method == 'DD':
            return self._call_dd()
        elif method != 'GPV':
            raise ValueError
        else:
            v = 0
            sigma, Ared, G, norms = self.final_sigma, self.Ared, self._G, self.gs_norms
            n = Ared.nrows()
            c = zero_vector(n)
            zs = []
            for i in range(n)[::-1]:
                b_ = G.column(i)
                c_ = c.dot_product(b_) / (norms[i]**2)
                sigma_ = sigma/norms[i]
                assert(sigma_ > 0)
                z = DiscreteGaussianDistributionIntegerSampler(sigma=sigma_, c=c_, algorithm="uniform+online")()
                c = c - z*Ared.column(i)
                v = v + z*Ared.column(i)
                zs.append(z)
        return v, zs[::-1]

    def split_primes(self, min_prime, max_prime):
        """
        a bunch of split primes
        """
        return _split_primes(self.m,self.H1,min_prime = min_prime, max_prime = max_prime)
    def _modq_dict(self,q):
        """
        a sanity check of the generators modulo q.
        """
        cc = self.cosets
        vv = self.vecs_modq(q, reduced = False)
        return dict(zip(cc,vv))

    def vecs_modq(self,q, reduced = True):
        """
        the basis elements (normal integral basis) modulo q.

        If reduced is true, return the LLL-reduced basis mod q
        """
        m = self.m
        if Integers(m)(q) not in self.H1:
            raise ValueError('q (= %s) is not a split prime'%q)
        v = finite_cyclo_traces(m,q,self.cosets,self.H1) # could be slow
        if not reduced:
            result = v
        else:
            result =  vector(v)*self._T
        return result

    """
    def min_vecs(self,q, stds = None):
        v = self.vecs_modq(q)
        red_roots =[a.lift() if a.lift() < q//2 else a.lift() - q for a in v]

        min_norm, scale =
        if stds is None:
            # elos case
            stds = [1 for _ in range(self.degree)]
        for b in range(2,q//2):
            b = GF(q)(b)
            scaled_roots  =[b*a for a in v]
            red_roots =[a.lift() if a.lift() < q//2 else a.lift() - q for a in scaled_roots]
            _norm = RR(sqrt(sum([(abs(a)*s)**2 for a,s in zip(red_roots, stds)])))
            if  _norm > min_norm:
                scale = b
                min_norm = _norm
      return [a*scale for a in v], scale
    """

    #def elos_quality(self,q):
    #    """
    ##    See tookit paper, lemma 2.8, and [ELOS].
    ##    """
    #    adj = self.adj
    #   min_norm, scale = self.min_vecs_norm(q)
    #    return 4*RR(sqrt(2*pi*n))*adj*min_norm/ RR(q)


    def ltwo_quality(self,q):
        vecs = self.vecs_modq(q)
        _norm, scale = RR(sqrt(sum([ZZ(a)**2 for a in reduce_roots(vecs,q)]))), 1
        for b in range(2,q//2):
            newnorm =RR(sqrt(sum([ZZ(c)**2 for c in reduce_roots([b*a for a in vecs],q)])))
            if newnorm < _norm:
                _norm, scale = newnorm, b
        return 2*_norm/q, scale

    def loo_quality(self,q):
        vecs = self.vecs_modq(q)
        _norm, scale = RR(max([abs(a) for a in reduce_roots(vecs,q)])), 1
        for b in range(2,q//2):
            newnorm =RR(max([abs(c) for c in reduce_roots([b*a for a in vecs],q)]))
            if newnorm < _norm:
                _norm, scale = newnorm, b
        return 2*_norm/q, scale


    def lone_quality(self,q):
        """
        We use l
        """
        vecs = self.vecs_modq(q)
        _norm, scale = RR(sum([abs(ZZ(a)) for a in reduce_roots(vecs,q)])), 1
        for b in range(2,q//2+1):
            newnorm = RR(sum([abs(ZZ(c)) for c in reduce_roots([Mod(b*a,q) for a in vecs],q)]))
            if newnorm < _norm:
                _norm, scale = newnorm, b
        return 2*_norm/q, scale


    def simulated_run(self,q, scale =1,numsamples = 30):
        """
        we run a simulation to test the quality.
        and scaling.
        """
        vecs = self.vecs_modq(q)
        max_err = 0
        count_zero = 0
        for i in range(numsamples):
            e = self.__call__()
            b = sum([Mod(a*e,q) for a,e in zip(vecs,e)])*scale
            b = ZZ(b) if b <= q//2 else ZZ(b) - q
            max_err = max(max_err, RR(abs(b)))
            # print abs(b)
        ratio  = max_err*2 / q

        # print 'number of zeros = %s'%count_zero
        #print 'Done!'
        return ratio


