import mpmath as mp
mp.dps = 200; mp.pretty = False
from scipy import linalg
import numpy as np

def search_for_q(v, K, qlst, h = None):
    """
    Given a number field {K} and a nice basis {v} of O_K,
    search for primes q such that the error could be small
    modulo q.

    h -- the extra term of variance for 1.
    """
    result = []
    n = K.degree()
    #if h is not None:
    #    print 'h = %s'%h
    MAX_ratio = 0
    MAX_q = 0
    for q in qlst:
        if q.is_prime():
            # print 'q = %s'%q
            fracq = K.prime_above(q)
            if fracq.norm() != q:
                raise ValueError # a prime of degree 1
            Fq = K.residue_field(fracq)
            _roots  = [Fq(a) for a in v]
            max_ratio = 0
            for b in range(1, q//2):
            # scale with b
                # print 'b = %s'%b
                scaled_roots  =[Fq(b)*a for a in _roots]
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
            print'q = %s, ratio = %s'%(q, max_ratio)
            print '*************'
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
            print'q = %s, ratio = %s'%(q, max_ratio)
            print '*************'
        if max_ratio > MAX_ratio:
            MAX_ratio = max(max_ratio, MAX_ratio)
            MAX_q = q

    print 'q = %s, ratio = %s'%(MAX_q, MAX_ratio)


def nice_basis_twopower(k):
    m = 2**k
    L.<z> = CyclotomicField(m)
    a = z+z^(-1)
    K.<c> = NumberField(a.minpoly())
    result = [1,c]
    j = 2
    while j < m//4:
        result.append(c**j - sum([result[j-2*t]*binomial(j,t) for t in range(1, j//2+1)]))
        j += 1
    return result, K

def another_prime_search(f, _range = 10, _min = 1e4, _max = 1e8):
    """
    This works the root of small order attack, when the field is monogenic.
    """
    f = QQ[x](f)
    result = []
    for b in range(-_range, _range):
        fb = ZZ(f(x=b))
        for q in fb.prime_factors():
            roots = GF(q)[x](f).roots(multiplicities=False)
            for root in roots:
                if root != GF(q)(0):
                    order = root.multiplicative_order()
                    if order <= 10 and ((root,q),order) not in result:
                        result.append(((root,q),order))
    return result


def search(fval, K, _range = 10):
    global f, n, _lam, sigma
    f = QQ[x](fval)
    n = f.degree()
    K = NumberField(f, names ='a')
    sigma = 1
    lam = mp.mpf(abs(K.discriminant())**(1/(2*n)))*mp.mpf(mp.sqrt(2*mp.mpf(pi)))
    print 'f = %s'%f
    print 'lam = %s'%lam

    print 'searching for pairs (b,q) with f(b) = 0 mod q...'
    bqlst = prime_search(f,_range = _range)
    print 'got a list of %s candidates.'%len(bqlst)
    print 'final good pairs are: %s'%sift(bqlst)


def prime_search(f, _range = 10,_min = 1e4, _max = 1e9):
    """
    """
    f = QQ[x](f)
    result = []
    for b in range(-_range, _range):
        fb = ZZ(f(x=b))
        qlist = fb.prime_factors()
        for q in qlist:
            if _min < q and q < _max:
                result.append((b,q))
    return result

def sift(lstbq,K, prec = 500):
    global f, lam, sigma
    lam = mp.mpf(abs(K.discriminant())**(1/(2*n)))*mp.mpf(mp.sqrt(2*mp.mpf(pi)))
    result = []
    for b, q in lstbq:
        print 'dealing with the pair (b,q) = (%s,%s)...'%(b,q)
        g = QQ[x](f(x = x-b))
        Kg.<a> = NumberField(g);
        vg = Kg.integral_basis()
        allbutone = check_zero(vg,Kg)
        if not allbutone:
            print 'not all but one basis vector reduce to 0, move to the next candidate'
            continue
        else:
            m1 = [[phi(c) for phi in Kg.complex_embeddings(prec = prec)] for c in vg] # expensive step
            Ng = mp.matrix(m1).transpose()
            print 'computing the SVD...'
            U,s,V = mp.svd(Ng) # expensive step
            print 'SVD done.'
            V0 = V.transpose_conj().column(0)
            print 'Singular values are %s'%s
            print 'V0 = %s with norm %s'%(V0, mp.norm(V0))
            rho = mp.norm(mp.matrix([V0[i,0]*s[i,0]**(-1) for i in range(V0.rows)]))
            print 'Norm is %s'%rho
            quotient = q/(rho*lam*sigma)
            if quotient > 4:
                print 'this pair (%s,%s) is a good pair with quotient = %s'%(b,q,quotient)
                result.append((b,q,quotient))
            else:
                print 'quotient too small: moving forward...'
            print '**********separator**********'
    print 'Done!'
    return result

def check_zero(vg, Kg):
    """
    check if almost vectors reduce to 0 modulo a
    EXAMPLE::

    """
    a = Kg.gen()
    coord = a.coordinates_in_terms_of_powers()
    vpol = [QQ[x](coord(b)) for b in vg];
    result = [pol(x=0) == 0 for pol in vpol]
    return result == [False] + [True]*(len(result)-1)


def astar_a(v,K, totally_real = False):
    """
    return the matrix A*A, where A is the embedding matrix of v, i.e.,
    A[i,j] = \sigma_i(v_j).
    Note we use trace, so the computation is exact.
    """
    import mpmath as mp
    n = len(v)
    if totally_real:
        return Matrix([[(v[i]*v[j]).trace() for j in range(n)] for i in range(n)])
    else:
        # okay, forget about A^*A: just do A^T A instead. No: this is terrible. Don't do that.
        ## phi = find_conjugation(K)
        #print 'phi = %s'%phi
        A = embedding_matrix(v,K,prec = 100)
        return mp.matrix(list(A.conjugate_transpose()*A))

def t_matrix(n, prec = 53):
    C = ComplexField(prec)
    m = ZZ(n//2)
    eyem = Matrix.identity(m)
    A = eyem.stack(eyem);
    I = C((0,1))
    B = (eyem*C(I)).stack(-eyem*C(I));
    T = A.change_ring(C).augment(B);
    return C(1/sqrt(2))*T

def embedding_matrix(v, K, totally_real = False, real = False, prec = 53, test = False):
    """
    Note: this is assuming Galois, so either totally real or totally complex.
    In general, this function needs to be rewritten, To incoporate all real or
    complex embeddings. 
    input:
    v -- a list of algebraic numbers
    K -- the field of definition. We require K to be a Galois extension.
    output:
    M_v -- the embedding where the j-th column is the j-th complex embedding of v

    Example::
    sage: K.<a> = NumberField(x^2+1)
    sage: v = K.integral_basis();
    sage: embedding_matrix(v, real = True)
    array([[ 1.+0.j,  0.+0.j],
       [ 0.+0.j, -1.-0.j]])
    """


    C = ComplexField(prec)
    n = len(v)

    if totally_real:
        return Matrix(n,n,[[phi(a) for phi in K.real_embeddings(prec = prec)] for a in v]).transpose()
    if real:
        M = Matrix(C,n,n//2, [a.complex_embeddings(prec = prec)[::2] for a in v])
        MM = ((M + M.conjugate())/2).augment((M-M.conjugate())/(2*C(I)))
        return MM.transpose()
    else:
        T = t_matrix(n, prec = prec)
        #print 'T = %s'%T
        M = Matrix(C,n,n//2, [a.complex_embeddings(prec = prec)[::2] for a in v]).transpose()
        A = M.stack(M.conjugate())
    if test:
        return A
    else:
        return _real_part(T^(-1)*A)

def _real_part(mat):
    """
    """
    return Matrix([[entry.real_part() for entry in row] for row in list(mat)])

def find_conjugation(K):
    """
    find the complex conjugation automorphism of any
    totally complex Galois extension.
    """
    c = K.gen()
    n = K.degree()
    _absdiffs = []
    import mpmath as mp
    for phi in K.automorphisms():
        if phi(c) != c:
            # if mp.mpf(abs((phi(c) - c).complex_embedding(prec = 100).real_part())).ae(mp.mpf(0)):
            diff =phi(c).complex_embedding(prec = 10*n) - c.complex_embedding(prec=10*n).conjugate()
            _absdiffs.append(abs(diff))
            if mp.mpc(diff).ae(mp.mpc(0),1e-3):
                print 'diff = %s'%diff
                return phi
    raise ValueError('did not find complex conjugation, min diff = %s'%min(_absdiffs))

def length_of_product(lst1,lst2):
    """
    given two lists of numbers ai and bi, return
    L = \sqrt{\sum |a_ib_i|^2}.
    """
    lst1 = [CC(a) for a in lst1]
    lst2 = [CC(b) for b in lst2]
    return sqrt(sum([(abs(a*b))**2 for a,b in zip(lst1,lst2)]))


def comp_norm(ss,u,prec = 100):
    """
    given two vectors [a], [b] of the same length and real entries.
    return the norm of the component-wise product vector (ab).
    """
    RF = RealField(prec)
    ss, u = [RF(a) for a in list(ss)],[RF(a) for a in list(u)]
    return RF(sqrt(sum([(a*b)^2 for a,b in zip(ss,u)])))
