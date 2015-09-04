class SubgroupModm:
    """
    a subgroup of (Z/mZ)^*
    """

    def __init__(self,m, gens, elements = None):
        self.m = m
        self.phim = euler_phi(m)
        self.Zm = Integers(m)

        newgens = []
        for a in gens:
            a = self.Zm(a)
            if not a.is_unit():
                raise ValueError('the generator %s must be a unit in the ambient group.'%a)
            newgens.append(a)

        self.gens = newgens

        if elements is None:
            print 'computing group elements...'
            t = cputime()
            self.H1 = self.compute_elements() # long time
            print 'Time = %s'%cputime(t)
            sys.stdout.flush()
        else:
            self.H1 = elements

        self.order = len(self.H1)
        print 'group order = %s'%self.order
        sys.stdout.flush()

        self._degree = ZZ(self.phim // self.order)

        print 'computing coset representatives...'
        t = cputime()
        self.cosets = self.cosets()
        print 'Time = %s'%cputime(t)
        sys.stdout.flush()

        self._is_totally_real = self.is_totally_real()

        if not self._is_totally_real:
            merged_cosets = []
            for c in self.cosets:
                if not any([-c/d in self.H1 for d in merged_cosets]):
                    merged_cosets.append(ZZ(c))
            newcosets = merged_cosets + [-a for a in merged_cosets]
            self.cosets = newcosets

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
        gens = self.gens
        result = [self.Zm(1)]
        for gen in gens:
            if gens != self.Zm(1):
                order = gen.multiplicative_order()
                pows = [gen**j for j in range(order)]
                result = set([a*b for a in result for b in pows])
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
            if euler_phi(m) ==  len(result)*len(elts): # already have enough cosets
                return result

    @cached_method
    def coset(self, a):
        """
        elt -- an integer
        returns the coset representative for this element
        """
        Zm = self.Zm
        for bb in self.cosets:
            if Zm(a)/Zm(bb) in set(self.H1):
                return bb
        raise ValueError('did not find a coset.')

    def extension_degree(self,vec):
        """
        vec -- a vector indexed by cosets of self, representing an element z in K.
        return the degree of the extension QQ(z)/QQ.
        """
        try:
            vec = list(vec)
        except:
            raise ValueError('input can not be turned into a list. Please debug.')
        C = self.cosets
        ele_dict = dict([(a,b) for a,b in zip(C,vec) if b != 0])
        fixGpLen = 0
        for ll in C:
            fixed = True
            for a in ele_dict.keys():
                lla = self.coset(ll*a)
                try:
                    coef = ele_dict[lla]
                except:
                    fixed = False
                    break
                if coef != ele_dict[a]:
                    fixed = False
                    break
            if fixed:
                fixGpLen += 1
        return self._degree // fixGpLen


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

    def multiplicative_order(self, a):
        """
        return the multiplicative order of [a] in the quotien group G/H
        """
        m = self.m
        Zm = self.Zm
        if gcd(m,a) != 1:
            raise ValueError
        a = Zm(a)
        o = self._degree
        for dd in o.divisors()[:-1]:
            if a**dd in self.H1:
                return dd
        return o
        # raise ValueError('did not find multiplicative order? please debug')

    def discriminant(self):
        """
        return, up to sign, the discriminant of the fixed field of self as a subfield of Q(zeta_m).
        """
        return prod([chi.conductor() for chi in self._associated_characters()])

    def intersection(self, other):
        """
        intersection of two subgroups of the same m.
        """
        if self.m != other.m:
            raise ValueError('the underlying m of self and other must be same.')
        H1 = self.H1
        H1other = other.H1
        Hnew = Set(H1).intersection(Set(H1other))
        print 'size of intersection = %s'%len(Hnew)
        Hnew_reduced = _reduce_gens(self.m,Hnew)
        print 'reduced gens for intersection  = %s'%Hnew_reduced
        sys.stdout.flush()
        return SubgroupModm(self.m, Hnew_reduced, elements = Hnew)

def _reduce_gens(m,H1):
    """
    given a full group, get a short list of generators.
    """
    Zm = Integers(m)
    gens = set([])
    gensSpan = set([Zm(1)])
    for a in H1:
        if Zm(a) not in gensSpan:
            #print 'adding %s to the set of generators'%a
            sys.stdout.flush()
            ordera = Zm(a).multiplicative_order()
            #print 'order of a = %s'%ordera
            alst  = [Zm(a)**j for j in range(1, ordera)]
            newelts = set([cc*aa for cc in gensSpan for aa in alst])
            gensSpan  |=  newelts
            gens.add(a)
        #else:
            #print 'already in the span'
        #print 'length of span = %s'%len(gensSpan)
        if len(gensSpan) == len(H1):
            # found enough generators.
            return list(gens)
    raise ValueError('did not find enough generators.')



