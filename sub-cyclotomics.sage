def AllSubgroups(m):
    """
    input: a positive integer m
    output: all subgroups of Z/mZ^*, given in terms of SubgroupModm instances.
    """
    U = Integers(m).unit_group()
    values = U.gens_values()
    generators = [str(a) for a in U.gens()]
    dict_use = dict(zip(generators, values));
    result = []
    for H in U.subgroups():
        order = H.order()
        gens_modm = [eval(str(a).replace('^','**'),dict_use) for a in H.gens()]
        result.append(SubgroupModm(m,gens_modm))
    return result



class SubgroupModm:
    """
    a subgroup of (Z/mZ)^*
    """

    def __init__(self,m, gens):
        self.m = m
        self.Zm = Integers(m)

        newgens = []
        for a in gens:
            a = self.Zm(a)
            if not a.is_unit():
                raise ValueError('the generator %s must be a unit in the ambient group.'%a)
            newgens.append(a)

        self.gens = newgens

        self.H1 = None
        self.H1 = self.compute_elements() # long time
        self.order = len(self.H1)

    def __repr__(self):
        return "subgroup of (Z/%sZ)^* of order %s generated by %s"%(self.m, self.order, self.gens)

    def is_totally_real(self):
        """
        The fixed field Q(zeta_m)^H is totally real if and only if -1 mod m \in H.
        """
        return self.Zm(-1) in  self.compute_elements()

    def compute_elements(self):
        """
        core function. Gives all the group elements
        """
        if self.H1 is not None:
            return self.H1
        gens = self.gens
        result = [self.Zm(1)]
        for gen in gens:
            if gens != self.Zm(1):
                order = gen.multiplicative_order()
                pows = [gen**j for j in range(order)]
                result = [a*b for a in result for b in pows]
        return result

    def cosets(self):
        """
        another core function, assuming we have elements, this shouldn't be hard.
        """
        Zm = self.Zm
        elts = self.H1
        m = self.m
        result =[]
        explored = []

        for a in range(m):
            if gcd(a,m) == 1 and a not in explored:
                for h in elts:
                    explored.append(h*a)
                result.append(Zm(a))
        # print 'result = %s'%result
        if euler_phi(m) != len(result)*len(elts):
            raise ValueError
        return result

    def _check_cosets(self):
        """
        sanity check that the cosets has been computed correctly.
        """
        H1 = self.H1
        cosets = self.cosets
        from itertools import combinations
        return not any([c[1]*c[0]**(-1) in H1 for c in combinations(cosets, 2)])


    def __hash__(self):
        return hash((self.m,tuple(self.gens)))


    def _associated_characters(self):
        """
        Definition: a Dirichlet character chi of modulus m is associated to
        a subgroup H <= Z/mZ)^* if chi|_H = 1.

        return all the associated characters of self.
        """
        m, Zm = self.m, self.Zm
        G = DirichletGroup(m)
        H1 = Set(self.compute_elements())

        result =[]
        for chi in G:
            ker_chi = Set([Zm(a) for a in chi.kernel()]) # a list of integers
            if H1.issubset(ker_chi):
                result.append(chi)
        return result

    def discriminant(self):
        """
        return, up to sign, the discriminant of the fixed field of self as a subfield of Q(zeta_m).
        """
        return prod([chi.conductor() for chi in self._associated_characters()])


#def split_primes(K,dK, max_prime = 10):
#    return [q for q in primes(max_prime) if K.prime_above(q).norm() == q and Mod(dK,q) != 0]

def split_primes_new(m, H, max_prime = 10):
    Zm = Integers(m)
    return [p for p in primes(3, max_prime) if Zm(p) in [Zm(h) for h in H]]

def min_norm(q,reduced_vecs,variances):
    scaled_roots  =[a for a in reduced_vecs]
    # print scaled_roots
    red_roots =[a.lift() if a.lift() < q//2 else a.lift() - q for a in scaled_roots]
    _norm = RR(sqrt(sum([abs(a)**2*abs(s) for a,s in zip(red_roots,variances)])))
    min_norm, scale = _norm, 1
    for b in range(2,q//2):
        b = GF(q)(b)
        scaled_roots  =[b*a for a in reduced_vecs]
        # print scaled_roots
        red_roots =[a.lift() if a.lift() < q//2 else a.lift() - q for a in scaled_roots]
        _norm = RR(sqrt(sum([abs(a)**2*abs(s) for a,s in zip(red_roots,variances)])))
        if  _norm > min_norm:
            scale = b
            min_norm = _norm
    return min_norm, scale

def adj_factor(m,H,dK):
    """
    is this right? we need to multiply by the
    """
    extra_factor = RR(sqrt(2))
    # I think this new factor should be here, but I'll have to double check.
    # The safest thing to do is to run both programs with primes.
    # OKay, it is even better: I claim that the diagonal terms are
    # (hm + h^2 - h) (of course, when we do not do the 'pruning' of the multi-set).
    dprime = len(H)
    n = euler_phi(m)//dprime
    if Integers(m)(-1) in H: # equivalent to totally real
        extra_factor = 1
    return extra_factor*RR(sqrt(m*dprime))/RR((abs(dK))**(1/(2*n)))


def good_basis_composite(m, H):
    """
    First step to glory.

    Take a composite number m (we assume that it is odd and square free)
    return a list of algebraic integers in the field, such that...

    H -- a list of integers representing a subgroup of (Z/m)^*
    """
    Zm = Integers(m)
    result = {}
    #print 'H = %s'%H
    explored = Set([])
    for a in range(m): # not good, since we are still going through all of m.
        # print 'a = %s'%a
        if Zm(a) not in explored:
            reps = tuple([Zm(j*a) for j in H])
            reps_wo_duplicates = Set(reps)
            rep = reps[0] # we only need one representative!
            if rep in result:
                result[rep] += len(reps_wo_duplicates)
            else:
                result[rep] = len(reps_wo_duplicates)
            explored = explored.union(reps_wo_duplicates)

    # print 'reps for good basis = %s'%result

    #conso_result = {}
    #for a in result.keys():
    #    conso_result[gcd(a,m)] = result[a]
    #print 'consolidate all Galois conjugates'
    #print 'reps for good basis = %s'%conso_result
    return result


def finite_cyclo_traces(m,q,ilst,H):
    """
    trace of zeta_m^i under h in the finite field F_q^r.
    """
    if not q.is_prime():
        raise ValueError
    Zm = Integers(m)
    if Mod(q,m) not in H:
        raise ValueError('we require q mod m to be in H.')
    r = Zm(q).multiplicative_order()
    if r == 1:
        F.<alpha> = GF(q)
    else:
        F.<alpha> = GF(q^r, impl = 'pari_ffelt')
    g = F.multiplicative_generator()
    d = ZZ((q^r-1)/m)
    f_zetam = g^d
    Fq = GF(q)
    return [sum([f_zetam**(i*ZZ(h))for h in H]) for i in ilst]

def construct_K(reps_dict,L,H):
    """
    hopefully faster than the one downstairs
    going to return the minpolys and the field
    which is all we need...?
    """
    dprime = len(H)
    z = L.gen()
    for rep in reps_dict.keys():
        elt = sum([z**(ZZ(rep*h)) for h in H])
        felt = elt.minpoly()
        if felt.degree() == n // dprime:
            return L.subfield(elt)[0]
    return None
        # minpoly_dict.append(felt,mult)
    #return minpoly_dict, K

"""
def construct_basis_in_K(reps_dict, L,K):
    result = {}
    for rep, mult in reps_dict.items():
        elt = (z**rep).trace(K)
        for galelt in elt.galois_conjugates(K):
            result[galelt] = mult
    return result
"""


def search_for_q(v, K, qlst, h = 1):
    """
    Given a number field {K} and a nice basis {v} of O_K,
    search for primes q such that the error could be small
    modulo q.

    get the one that has the smallest 4*norm / q

    h -- the extra term of variance for 1.
    """
    result = []
    n = K.degree()
    #MAX_ratio = 0
    #MAX_q = 0
    Opt_ratio = 10
    Opt_q = 0
    #Scale = 0
    for q in qlst:
        if q.is_prime():
            overalltimer = cputime()

            # these steps can be very slow
            fracq = K.prime_above(q)
            Fq = K.residue_field(fracq)
            _roots  = reduce_roots([Fq(a) for a in v],q)
            # these are probably slow now I have a way to deal with them.
            #max_ratio = 0
            #scale = 0
            l1_norm=  RR(sum([abs(ZZ(a)) for a in _roots])) + RR(sqrt(h))

            l2_norm = RR(sqrt(sum([ZZ(a)**2 for a in _roots])+h))
            #loo_norm = max([abs(a) for a in reduce_roots(_roots,q)])
            ratio  = (4*l2_norm)/q
            #for b in range(1, q//2):
            #    b = Fq(b)
            #    scaled_roots  =[b*a for a in _roots]
            #    red_roots =[a.lift() if a.lift() < q//2 else a.lift() - q for a in scaled_roots]
            #    bz = b.lift() if b.lift() < q//2 else b.lift() - q
            #    if h is not None:
            #        red_norm = RR(sqrt(sum([a**2 for a in red_roots]) + bz**2*h))
            #    else:
            #        # red_norm = RR(sqrt(sum([a**2 for a in red_roots])))
            #        red_norm = max([abs(a) for a in red_roots]) # use l infinity norm, see if it is better.
            #    ratio = q / (4*red_norm)
            #    if ratio > max_ratio:
            #        scale = b
            #        best_roots = [bz] + red_roots
            #    max_ratio = max(ratio, max_ratio)
            #print'q = %s, ratio = %s'%(q, max_ratio)
            #print 'best_roots = %s'%best_roots
            #print 'time: %s'%cputime(overalltimer)
            #print '*************'
        else:
            # composite q. what we do here is consider each prime factor.
            R = Integers(q)
            vq = q.prime_divisors()
            vfracq = [K.prime_above(qi) for qi in vq]
            assert all([fracqi.norm() == qi for (fracqi, qi) in zip(vfracq,vq)]) # check that every prime dividing q splits

            # vF is the list of residue fields
            vF = [K.residue_field(fracqi) for fracqi in vfracq]
            _roots = []
            for a in v:
                reduced_alst = [F(a).lift() for F in vF]
                aq = CRT_list(reduced_alst, vq)
                _roots.append(R(aq))
            # We computed the image of the roots in Z/q using Chinese remainder theorem.
            print 'composite roots = %s'%_roots
            # Next we do the scaling. Maybe this part can be merged? I don't care for now.
            max_ratio = 0
            for b in range(q):
                if gcd(b,q) == 1:
                    b = R(b)
                    scaled_roots  =[b*a for a in _roots]
                    red_roots =[a.lift() if a.lift() < q//2 else a.lift() - q for a in scaled_roots]
                    if h is not None:
                        bz = b.lift() if b.lift() < q//2 else b.lift() - q
                        red_norm = RR(sqrt(sum([a**2 for a in red_roots]) + bz**2*h))
                    else:
                        red_norm = RR(sqrt(sum([a**2 for a in red_roots])))
                    ratio = q / red_norm
                    #if ratio > 2:
                    #    print [bz] + red_roots, ratio
                    max_ratio = max(ratio, max_ratio)
        # print'q = %s, ratio = %s'%(q, ratio)
        #print '*************'
        if ratio < Opt_ratio:
            Opt_ratio, Opt_q = ratio, q

    return Opt_q, Opt_ratio


def loo_norm(vecs,q):
    Fq = GF(q)
    scale = 1
    vecs = reduce_roots(vecs,q)
    min_norm = max([abs(a) for a in vecs]) # use l infinity norm, see if it is better.
    min_roots = vecs
    for b in range(2, q//2):
        b = Fq(b)
        scaled_roots  = reduce_roots([b*a for a in vecs],q)
        red_norm = max([abs(a) for a in scaled_roots]) # use l infinity norm, see if it is better.
        if min_norm > red_norm:
            min_norm, scale, min_roots = red_norm, b, scaled_roots
    return min_norm, scale, min_roots
