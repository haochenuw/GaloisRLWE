def t_matrix(n, prec = 53):
    C = ComplexField(prec)
    mm = ZZ(n//2)
    eyem = Matrix.identity(mm)
    A = eyem.stack(eyem);
    I = C((0,1))
    B = (eyem*C(I)).stack(-eyem*C(I));
    T = A.change_ring(C).augment(B);
    return C(1/sqrt(2))*T


def _real_part(mat):
    """
    real part of a matrix.
    """
    return Matrix([[matentry.real_part() for matentry in matrixrow] for matrixrow in list(mat)])


from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler


class SubCyclotomicsRLWEErrorSampler:
    """
    We write our own sampler. Hahaha.
    """

    def __init__(self,m,H,sigma,prec, disc = None):
        """
        require: m must be square free.

        disc: the discriminant of Q(zeta_m)^H. We pass it
        as an optional parameter, since when the order of H is
        large, the computation could be very slow.
        """

        self.m = m
        self.H = H
        self.H1 =  H.compute_elements()
        self.cosets = H.cosets()
        self.sigma = sigma
        self.prec = prec



        #compute adjustment factors.

        if disc is None:
            self.disc = H.discriminant() # could be very slow.
        else:
            self.disc = disc

        self.degree = euler_phi(m) // len(self.H1)

        self._is_totally_real = H.is_totally_real()


        if not self._is_totally_real:
            merged_cosets = []
            for c in self.cosets:
                if not any([-c/d in self.H1 for d in merged_cosets]):
                    merged_cosets.append(ZZ(c))
            newcosets = merged_cosets + [-a for a in merged_cosets]
            self.cosets = newcosets

        print 'coset = %s'%self.cosets


        self.adj = RR(self.disc**(1/(2*self.degree)))

        # somehow this factor is not needed, will investigate.
        #if not self._is_totally_real:
        #    self.adj /=  RR(sqrt(2))

        self.final_sigma = self.sigma*self.adj

        


        self.TstarA, self.A = self.embedding_matrix(prec = self.prec)


        #self._G = self.compute_G()


        # self._G = self.compute_G()


        # self.stds = [self.final_sigma/ss for ss in self.compute_stds(prec = self.prec)]

        #self._G = self.produce_norms(self.A) # well, need to write this. Also need to make it work
        # #or arbitrary precision.

        # # self.D = [ DiscreteGaussianDistributionIntegerSampler(sigma=ss) for ss in self.stds]


    def __repr__(self):
        return 'RLWE error sampler with H = %s'%self.H

    def compute_G(self, prec = 53):
        B = self.TstarA
        n = self.degree
        n = self.degree
        from mpmath import *
        mp.dps = prec // 2
        BB = mp.matrix([list(w) for w in list(B)])
        Q,R = qr(BB) # QR decomposition
        M = mp.matrix([list(Q.column(i)*R[i,i]) for i in range(n)]);
        M_sage = Matrix([[RealField(prec)(M[i,j]) for i in range(n)] for j in range(n)])
        return M_sage # vectors are columns

        # print 'gram-schmidt vector lengths = %s'%v
        return v

    def embedding_matrix(self, prec = 53):
        """
        to-do: need to separate totally real case from totally complex
        """
        H1 = self.H1
        C = ComplexField(prec)
        zetam = C.zeta(m)
        cosets = self.cosets
        #print 'cosets = %s'%cosets
        n = self.degree

        if self._is_totally_real:
            #print  'totally real case'
            A = _real_part(Matrix([[sum([zetam**(ZZ(l*k*h)) for h in H1]) for l in cosets] for k in cosets]))
            return A,A
        else:
            # we need to merge a with -a mod h.

            #print 'reordered cosets = %s'%newcosets
            # need to rearrange the columns so that the j + n/2 th row is the complex conjugate
            # of the j th row.
            # condense(cosets)

            T = t_matrix(n,prec = prec)
            A = Matrix([[sum([zetam**(ZZ(l*k*h)) for h in H1]) for l in cosets] for k in cosets])
            # print A0.det(), T.det()
            return _real_part(T.conjugate_transpose()*A),A

    def inverse_of_embedding_matrix(self,**kwds):
        """
        maybe I need this, to convert back to a number field element.
        """


    def _gram_schmidt_vectors(self):
        """
        well, this can't be done using a sage matrix, so maybe I have to modify it anyway.
        but if I just need the norm, then maybe its not too bad
        """
        # I want to make it arbitrary precision.
        pass

    def coset_reps(self):
        """
        I need this for representing the basis vectors. Each coset rep c
        represents the element \alpha_c =  \sum_{h \in H} \zeta_m^{ch}.
        """
        return self.cosets

    def __call__(self, minkowski = False):
        """
        return an integer vector a = (a_c) indexed by the coset reps of self,
        which represents the vector \sum_c a_c \alpha_c
        Use the algorithm of [GPV].
        http://www.cc.gatech.edu/~cpeikert/pubs/trap_lattice.pdf
        """
        v = 0
        # also need matrix. Okay.
        sigma, A, G = self.final_sigma, self.TstarA, self._G
        n = A.nrows()
        c = zero_vector(n)

        zs = []
        for i in range(n)[::-1]:
            b_ = G.column(i)
            c_ = c.dot_product(b_) / b_.dot_product(b_)
            sigma_ = sigma/b_.norm()
            assert(sigma_ > 0)
            z = DiscreteGaussianDistributionIntegerSampler(sigma=sigma_, c=c_, algorithm="uniform+online")()
            c = c - z*A.column(i)
            v = v + z*A.column(i)
            zs.append(z)
        if minkowski:
            return v
        else:
            return zs


    def vecs_modq(self,q):
        m = self.m
        if Integers(m)(q) not in self.H1:
            raise ValueError('q (= %s) is not a split prime'%q)
        return finite_cyclo_traces(m,q,self.cosets,self.H1)

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

    ##def peikert_quality(self,q):
    ###    min_norm, scale = self.min_vecs_norm(q, stds = self.stds)
    #
    #    print 'scale = %s'%scale
    ##    return RR(min_norm*4)/ RR(q)

    def simulated_run(self,q, scale =1):
        """
        we run a simulation to test the quality. Remember to set tolerance.
        """
        vecs = self.vecs_modq(q)
        #print 'scale = %s'%scale
        max_err = 0
        count_zero = 0
        for i in range(q):
            e = self.__call__()
            b = sum([a*e for a,e in zip(vecs,e)])
            b = ZZ(b) if b < q//2 else q - ZZ(b)
            if b ==0:
                count_zero +=1
            max_err = max(max_err, RR(abs(b)))
            print abs(b)
        ratio  = max_err*2 / q

        print 'number of zeros = %s'%count_zero
        #print 'Done!'
        return ratio


