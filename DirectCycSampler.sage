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

    def degree_n_primes(self,min_prime,max_prime, n = 1):
        result = []
        for p in primes(min_prime, max_prime):
            try:
                if (Integers(m)(p)).multiplicative_order() == n:
                    result.append(p)
            except:
                pass
        return result

    def __repr__(self):
        return 'RLWE cyclotomic sampler with m = %s and sigma = %s'%(self.m, self.sigma)

    def vecs_modq(self,q):
        return GF(q)[x](self.f).roots(multiplicities=False)

    def __call__(self):
        return [self.D() for _ in range(self.n)]

    def degree_of_prime(self,q):
        try:
            return (Integers(self.m)(q)).multiplicative_order()
        except:
            raise ValueError('q must be unramified in self.')



