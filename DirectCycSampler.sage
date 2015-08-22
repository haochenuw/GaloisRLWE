from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler


class DirectCycSampler:
    """
    for sampling (coefficient-wise) from any cyclotomic field
    """
    def __init__(self,m, sigma = 1):
        self.m = m
        self.sigma = sigma
        self.D = DiscreteGaussianDistributionIntegerSampler(sigma = sigma)
        self.f = cyclotomic_polynomial(m)
        self.n = self.f.degree()
        K.<z> = CyclotomicField(m)
        self.K = K
        self.z = K.gen()
        self.secret = self._to_field(self.__call__())
        self.ff = self.z.coordinates_in_terms_of_powers()
        self.A = general_embedding_matrix(z, K, prec = 200)
        self.DD = MyLatticeSampler(self.A)

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

    def _a_root_mod_q(self,q):
        deg  = self.degree_of_prime(q)
        if deg  > 1:
            F.<alpha> = GF(q^deg, impl = 'pari_ffelt')
        else:
            F = GF(q)
        return  F[x](self.f).roots(multiplicities=False)[0]

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
        return sum([F(c[i])*aa**i for i in range(len(lst))])


    def __call__(self):
        return [self.D() for _ in range(self.n)]

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
            raise ValueError('q must be unramified in self.')

    def _uniform_a(self,q):
        return self._to_field([ZZ.random_element(q) for _ in range(self.n)])

    def set_sigma(self,newsigma):
        self.sigma = newsigma

    def set_secret(self, newsecretvec):
        self.secret = self._to_field(newsecretvec)

    def rlwe_sample(self,q, add_error = True):
        """
        generate an rlwe sample
        """
        a = self._uniform_a(q)
        avec = self._to_vec(a)
        s = self.secret
        b = a*s
        if add_error:
            e = self._to_field(self.__call__())
            b += e
        bvec = self._to_vec(b)
        return (avec, [Mod(bi,q) for bi in bvec])

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
        alpha_b = [ZZ(ai)*alpha for ai in a]
        round_alpha_a = list(D.babai(c = A*vector(alpha_a))[1]) # an approximation of scaled_a.
        round_alpha_b = list(D.babai(c = A*vector(alpha_b))[1])
        return [round_alpha_a, [Mod(bi,newq) for newq in round_alpha_b]]

    def elos_chisquare_attack(self,samples,q):
        """
        """
        DD = self.DD
        s = self.secret
        print 's = %s'%s
        for a, b in samples:





