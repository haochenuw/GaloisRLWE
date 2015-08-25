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

    def __init__(self,m,H,sigma = 1,prec = 100, method = 'BKZ'):
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

        self.TstarA, self.Acan = self.embedding_matrix(prec = self.prec)

        self.Acaninv = (self.Acan)**(-1)

        #self.disc = (self.TstarA).det() #maybe this is faster in computing discriminants.

        #self.adj = RR(abs(self.disc)**(1.0/(self._degree)))

        #self.final_sigma = self.sigma*self.adj
        #self._T = self.lll_transform_matrix()
        #self.Ared = self.TstarA*self._T

        self.D = MyLatticeSampler(self.TstarA, sigma = self.sigma, method = method)
        self.Ared = self.D.B

        self.final_sigma  =self.D.final_sigma
        # gram-schmidt basis and norms.
        #self._G, self.gs_norms = self.compute_G()
        # self._dd_gen = DiscreteGaussianDistributionIntegerSampler(sigma = sigma)

        # # self.D = [ DiscreteGaussianDistributionIntegerSampler(sigma=ss) for ss in self.stds]
        self.secret = self.__call__()


    def __repr__(self):
        return 'RLWE error sampler with m = %s,  H = %s, secret  = %s and sigma = %s'%(self.m, self.H, self.secret, self.final_sigma.n())

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

    def degree_of_prime(self,q):
        """
        return the degree of q in K
        """
        if not q.is_prime():
            raise ValueError('q must be prime')
        return (self.H).multiplicative_order(q)

    def degree_n_primes(self, min_prime, max_prime, n =1):
        """
        return a bunch of primes of degree n in K. When n = 1, this
        is split primes.
        """
        result = []
        for p in primes(min_prime, max_prime):
            try:
                if self.degree_of_prime(p) == n:
                    result.append(p)
            except:
                pass
        return result

    def basis_lengths(self):
        return [self.Ared.column(i).norm() for i in range(self._degree)]

    def embedding_matrix(self, prec = None):
        """
        We are in a simplified situation because the field K is Galois over QQ,
        so it is either totally real or totally complex.
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


    # deprecated.
    def _lll_transform_matrix(self):
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


    def __call__(self,c = None):
        """
        return an integer vector a = (a_c) indexed by the coset reps of self,
        which represents the vector \sum_c a_c \alpha_c
        Use the algorithm of [GPV].
        http://www.cc.gatech.edu/~cpeikert/pubs/trap_lattice.pdf

        If minkowski = True, return the lattice vector in R^n. Otherwise,
        return the coordinate of the vector in terms of the embedding matrix of self.
        """
        return self.D(c = c)[1]

    def babai(self,c):
        return self.D.babai(c)[1]


    # deprecated.
    """
    def __call__(self, c = None, method = 'GPV', reduced = True):
        if method == 'DD':
            return self._call_dd()
        elif method != 'GPV':
            raise ValueError
        else:
            v = 0
            sigma, Ared, G, norms = self.final_sigma, self.Ared, self._G, self.gs_norms
            n = Ared.nrows()
            if c is None:
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
        if reduced:
            return v, vector(zs[::-1])
        else:
            return v, self._T*vector(zs[::-1])
    """

    def _modq_dict(self,q):
        """
        a sanity check of the generators modulo q.
        """
        cc = self.cosets
        vv = self.vec_modq(q)
        return dict(zip(cc,vv))


    @cached_method
    def vec_modq(self,q, reduced = False):
        """
        the basis elements (normal integral basis) modulo q.

        If reduced is true, return the LLL-reduced basis mod q

        v dot Tz = (vT) dot z
        """
        m = self.m
        degree = self.degree_of_prime(q)
        # print 'degree of %s in self = %s'%(q,degree)
        #if Integers(m)(q) not in self.H1:
        #    raise ValueError('q (= %s) is not a split prime'%q)
        v = finite_cyclo_traces(m,q,self.cosets,self.H1, deg = degree) # could be slow
        if not reduced:
            result = vector(v)
        else:
            result =  vector(v)*self._T
        return result

    def _to_ccn(self, lst):
        """
        convert an element in O_K from C^n to Z^n.
        """
        return list(self.Acan*vector(lst))


    def _to_zzn(self,lst):
        """
        the inversion of the above.
        """
        return list(self.Acaninv*vector(lst))


    # methods for modulus switching.
    def _uniform_a(self,q):
        return [ZZ.random_element(q) for _ in range(self._degree)]


    # in order to be fast, we will sacrifice the speed and compute the product

    def _prod(self,lsta, lstb):
        """
        multiplying two field elements, using the canonical embedding
        """
        lsta, lstb = list(lsta), list(lstb)
        lsta_cc, lstb_cc = self._to_ccn(lsta), self._to_ccn(lstb)
        float_result = self._to_zzn([aa*bb for aa, bb in zip(lsta_cc,lstb_cc)])
        return [ZZ(round(tt.real_part())) for tt in float_result]

    # If we do exact field calculation this will be correct but slow, anyway.


    # RLWE methods
    def rlwe_sample(self,q, add_error = True):
        """
        generate an rlwe sample.
        """
        a = self._uniform_a(q)
        s = self.secret
        b = self._prod(a,s)
        if add_error:
            e = self.__call__()
            #verbose('e = %s'%e)
            newb = [bi + ei for bi, ei in zip(b,e)]
        return (a, [Mod(bi,q) for bi in newb])

    def set_sigma(self,newsigma):
        self.D.final_sigma = newsigma

    def set_secret(self, newsecret):
        self.secret = newsecret


    # this modulus switch.
    def modulus_switch(self,oldq, newq, sample, method = 'Babai', newsigma = None):
        """
        switch a sample from an old modulus to a new one.

        method -- 'GPV' or 'Babai'. The rounding method used.

        return two lists.
        """
        TstarA = self.TstarA
        alpha  = QQ(newq/oldq)
        a, b = sample
        alpha_a = [ZZ(ai)*alpha for ai in a]
        alpha_b = [ZZ(bi)*alpha for bi in b]


        if method == 'Babai':
            round_alpha_a = list(self.babai(TstarA*vector(alpha_a))) # an approximation of scaled_a.
            round_alpha_b = list(self.babai(TstarA*vector(alpha_b)))
        else:
            raise NotImplementedError
        # switch back.


        bprime_lst = _my_list_diff(alpha_b, round_alpha_b)
        #print 'bprime = %s'%bprime_lst
        return (round_alpha_a, [Mod(bi,newq) for bi in round_alpha_b])



    # methods for simulating attacks.


    def _map_to_finite_field(self,lst,q, vec):
        """
        as advertised.
        """
        return sum([aa*bb for aa, bb in zip(lst,vec)])


    def chisquare_attack(self,q,samples):
        """
        Perform a chisquared attack.
        """
        verbose('q = %s'%q)
        vec = self.vec_modq(q)
        s = self.secret
        sbar = self._map_to_finite_field(s, q, vec)
        verbose('sbar = %s'%sbar)

        verbose('vec = %s'%vec)

        degq = self.degree_of_prime(q)

        FF = sbar.parent()
        verbose('finite field  =  %s'%FF)

        errors_dict = dict([(cc,0) for cc in FF])

        verbose('mapping samples to finite field...')
        for a,b in samples:
            abar, bbar = self._map_to_finite_field(a, q, vec), self._map_to_finite_field(b, q, vec)
            ebar = bbar  - sbar*abar
            errors_dict[ebar] += 1

        bins = selecting_bins(q, degq, len(samples))
        print 'number of bins = %s'%bins
        return chisquare_test(errors_dict, bins = bins, std_multiplier = 2)
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
        vec = self.vec_modq(q)
        _norm, scale = RR(sqrt(sum([ZZ(a)**2 for a in reduce_roots(ve,q)]))), 1
        for b in range(2,q//2):
            newnorm =RR(sqrt(sum([ZZ(c)**2 for c in reduce_roots([b*a for a in vec],q)])))
            if newnorm < _norm:
                _norm, scale = newnorm, b
        return 2*_norm/q, scale

    def loo_quality(self,q):
        vec = self.vec_modq(q)
        _norm, scale = RR(max([abs(a) for a in reduce_roots(vec,q)])), 1
        for b in range(2,q//2):
            newnorm =RR(max([abs(c) for c in reduce_roots([b*a for a in vec],q)]))
            if newnorm < _norm:
                _norm, scale = newnorm, b
        return 2*_norm/q, scale


    def lone_quality(self,q):
        """
        We use l
        """
        vec = self.vec_modq(q)
        _norm, scale = RR(sum([abs(ZZ(a)) for a in reduce_roots(vec,q)])), 1
        for b in range(2,q//2+1):
            newnorm = RR(sum([abs(ZZ(c)) for c in reduce_roots([Mod(b*a,q) for a in vec],q)]))
            if newnorm < _norm:
                _norm, scale = newnorm, b
        return 2*_norm/q, scale


    def simulated_run(self,q, scale =1,numsamples = 30):
        """
        we run a simulation to test the quality.
        and scaling.
        """
        vec = self.vec_modq(q)
        max_err = 0
        count_zero = 0
        for i in range(numsamples):
            e = self.__call__()
            b = sum([Mod(a*e,q) for a,e in zip(vec,e)])*scale
            b = ZZ(b) if b <= q//2 else ZZ(b) - q
            max_err = max(max_err, RR(abs(b)))
            # print abs(b)
        ratio  = max_err*2 / q

        # print 'number of zeros = %s'%count_zero
        #print 'Done!'
        return ratio


