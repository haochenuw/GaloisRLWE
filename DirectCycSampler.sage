from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler


class DirectCycSampler:
    """
    for sampling (coefficient-wise) from any cyclotomic field
    """
    def __init__(self,m, sigma = 1):
        self.m = m
        self.sigma = sigma
        #self.D = DiscreteGaussianDistributionIntegerSampler(sigma = sigma)
        self.Gaussian = RealDistribution('gaussian', sigma)
        self.f = cyclotomic_polynomial(m)
        self.n = self.f.degree()
        K.<z> = CyclotomicField(m)
        self.K = K
        self.z = K.gen()
        self.secret = self._to_field(self.__call__())
        self.ff = self.z.coordinates_in_terms_of_powers()
        #self._embedding_matrix = general_embedding_matrix(z, K, prec = 200)
        #self.DD = MyLatticeSampler(self._embedding_matrix)

    def degree_n_primes(self,min_prime,max_prime, n = 1):
        result = []
        m = self.m
        for p in primes(min_prime, max_prime):
            try:
                if (Integers(m)(p)).multiplicative_order() == n:
                    result.append(p)
            except:
                pass
        return result

    def __repr__(self):
        return'RLWE cyclotomic sampler with m = %s, sigma = %s, and secret = %s'%(self.m, self.sigma, self.secret)

    @cached_method
    def _a_root_mod_q(self,q):
        deg  = self.degree_of_prime(q)
        if deg  > 1:
            F.<alpha> = GF(q^deg, impl = 'pari_ffelt')
        else:
            F = GF(q)
        aa = F[x](self.f).roots(multiplicities=False)[0]
        print 'found %s as a root mod %s'%(aa,q)
        return aa


    def vecs_modq(self,q):
        a = self._a_root_mod_q(q)
        return [a**i for i in range(self.n)]

    def _map_to_fq(self,lst,q):
        """
        maps a list (or a field element) to the finite field Fq
        """
        if not isinstance(lst, list):
            lst = self._to_vec(lst)
        aa = self._a_root_mod_q(q)
        F = aa.parent()
        return sum([F(lst[i])*aa**i for i in range(len(lst))])


    def __call__(self):
        """
        sampling done using [DD].
        """
        v0 = [self.Gaussian.get_random_element() for _ in range(self.n + 1)]
        v1 = [a - v0[-1] for a in v0[:-1]]
        return [ZZ(round(float(vv))) for vv in v1]



    def _to_field(self,lst):
        """
        convert from a list to a field.
        """
        return sum([ZZ(lst[i])*self.z**(i+1) for i in range(len(lst))])

    def _to_vec(self, elt):
        return self.ff(elt)

    def degree_of_prime(self,q):
        try:
            return (Integers(self.m)(q)).multiplicative_order()
        except:
            return ZZ(log(self.K.prime_above(q).norm(),q))

    def _uniform_a(self,q):
        return [ZZ.random_element(q) for _ in range(self.n)]

    def light_elos(self,q):
        pass

    #def set_sigma(self,newsigma):
    #    self.sigma = newsigma
    #    self.D = DiscreteGaussianDistributionIntegerSampler(sigma = newsigma)


    def set_secret(self, newsecretvec):
        self.secret = self._to_field(newsecretvec)

    def rlwe_sample(self,q, add_error = True):
        """
        generate an rlwe sample
        """
        a = self._uniform_a(q)
        s = self.secret
        b = self._to_field(a)*s
        if add_error:
            e = self._to_field(self.__call__())
            b += e
        bvec = self._to_vec(b)
        return (a, [Mod(bi,q) for bi in bvec])

    def modulus_switch(self,oldq, newq, sample):
        """
        switch a sample from an old modulus to a new one.
        will use MyLatticeSampler for the rounding.
        also I will use embedding_matrix. I think.

        return two lists
        """
        A = self._embedding_matrix
        alpha  = QQ(newq/oldq)
        a, b = sample
        alpha_a = [ZZ(ai)*alpha for ai in a]
        alpha_b = [ZZ(bi)*alpha for bi in b]
        #print 'alphab = %s'%alpha_b
        DD = self.DD
        round_alpha_a = list(DD.babai(c = A*vector(alpha_a))[1]) # an approximation of scaled_a.
        round_alpha_b = list(DD.babai(c = A*vector(alpha_b))[1])

        bprime_lst = _my_list_diff(alpha_b, round_alpha_b)
        #print 'bprime = %s'%bprime_lst
        return [round_alpha_a, [Mod(bi,newq) for bi in round_alpha_b]]


    def elos_attack(self,q, samples, maxRatio = 2):
        """
        the elos attack
        """
        print 'q = %s'%q
        s = self.secret
        print 's = %s'%s
        sbar = self._map_to_fq(s, q)
        print 'sbar = %s'%sbar
        F = sbar.parent()
        reduced_samples = []
        for a,b in samples:
            abar, bbar = self._map_to_fq(a, q), self._map_to_fq(b, q)
            reduced_samples.append((abar,bbar))
        G = []
        for sguess in F:
            # print 'sguess = %s'%sguess
            good = True
            for abar, bbar in reduced_samples:
                ebar = bbar - abar*sguess
                ebarReduced = ZZ(ebar) if ZZ(ebar) < q//2 else ZZ(ebar) - q
                if float(q/(2*(ebarReduced))) < maxRatio:
                    good = False
                    break
            if good:
                G.append(sguess)
        if len(G) == 0:
            return 'not rlwe'
        elif len(G) > 1:
            return 'not enough samples'
        else:
            if G[0] == sbar:
                return 'success'
            else:
                return 'failed.'

                    # error too large

             # run a chisquare test


    def chisquare_attack(self,q,samples):
        """
        Note that this only works for one
        """
        print 'q = %s'%q
        s = self.secret
        sbar = self._map_to_fq(s, q)
        print 'sbar = %s'%sbar
        F = sbar.parent()
        deg = self.degree_of_prime(q)

        a_dict = dict([(cc,0) for cc in F])
        #b_dict = dict([(cc,0) for cc in F])
        reduced_samples = []
        for a,b in samples:
            abar, bbar = self._map_to_fq(a, q), self._map_to_fq(b, q)
            reduced_samples.append((abar,bbar))
            a_dict[abar] += 1

        bins = selecting_bins(q,deg, len(samples))
        print 'bins = %s'%bins

        result = []
        for sguess in F:
            errors_dict = dict([(cc,0) for cc in F])
            print 'chisquare value for sguess = %s'%sguess
            for abar, bbar in reduced_samples:
                ebar = bbar - abar*sguess
                #ebar = abar*sguess
                errors_dict[ebar] += 1
            _,chisquare = chisquare_test(errors_dict, bins = bins, std_multiplier =2)
            result.append((sguess,chisquare))


        print 'chisquare for a:'
        _,_ = chisquare_test(a_dict, bins = bins, std_multiplier =2)
        #print 'chisquare of b:'
        #_,_ = chisquare_test(b_dict, bins = bins, std_multiplier =2)
        #print
        return result
        # return chisquare_test(errors_dict, bins = bins, std_multiplier =2)[0]





