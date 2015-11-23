# some misc functions for use.

from sage.libs.fplll.fplll import FP_LLL

import sys

######################################
# some other statistical uniform tests.
######################################
def finitefield_smallness_test(samples, probThreshold = 1e-5, tolerance = 0.5):
    F = samples[0].parent()
    q = F.characteristic()
    degF = F.degree()
    if degF > 1:
        raise NotImplementedError
    numsamples = len(samples)
    ratio = float(numsamples/q)
    esmall = ratio*(floor(float(q/2*tolerance))*2 + 1)
    elarge = numsamples - esmall
    nsmall = 0
    for e in samples:
        ez = ZZ(e)
        if ez > q/2:
            ez -= q
        if abs(ez) < q/2*tolerance:
            nsmall += 1
    nlarge = numsamples - nsmall
    print 'esmall, elarge =  %s, %s'%(esmall, elarge)
    print 'nsmall, nlarge = %s, %s'%(nsmall, nlarge)
    if min(esmall, elarge) < 5:
        raise ValueError('samples size too small.')
    # Now we have two bins, we do a very tiny chisquare test.
    chisquare = (nsmall - esmall )^2/esmall + (nlarge - elarge)^2/elarge
    T = RealDistribution('chisquared', 1)
    print('chisquare = %s'%chisquare)
    prob = 1 - T.cum_distribution_function(chisquare)
    if prob < probThreshold:
        verbose('non-uniform')
        return False
    else:
        verbose('uniform')
        return True

def subfield_unifrom_test(samples, probThreshold = 1e-5):
    """
    we assume that the samples are from a finite field.
    we separate the ones that are from a proper subfield.
    """
    F = samples[0].parent()
    q = F.characteristic()
    degF = F.degree()
    numsamples = len(samples)
    eltsWithFullDegree = elts_of_full_degree(q,degF)
    nSmall = 0
    nLarge = 0
    for aa in samples:
        if aa.minpoly().degree() < degF:
            nSmall +=1
        else:
            nLarge +=1
    card = q**degF
    eLarge= float(eltsWithFullDegree/card*numsamples)
    eSmall= numsamples - eLarge
    verbose('eSmall, eLarge = %s,%s'%(eSmall, eLarge))
    verbose('nSmall, nLarge = %s,%s'%(nSmall, nLarge))
    if min(eSmall, eLarge) < 5:
        raise ValueError('samples size too small.')
    # Now we have two bins, we do a very tiny chisquare test.
    chisquare = (nSmall - eSmall )^2/eSmall + (nLarge - eLarge)^2/eLarge
    T = RealDistribution('chisquared', 1)
    verbose('chisquare = %s'%chisquare)
    prob = 1 - T.cum_distribution_function(chisquare)
    if prob < probThreshold:
        verbose('non-uniform')
        return False
    else:
        verbose('uniform')
        return True

def elts_of_full_degree(q,n):
    """
    number of elements of F_q^n that do not lie in proper subfields
    """
    return sum([q**(n//d)*moebius(d) for d in n.divisors()])




def adaptive_test(samples, threshold = 100):
    ambientSet = samples[0].parent()
    try:
        numkeys = len(ambientSet)
    except:
        raise ValueError('')
    numsamples = len(samples)
    cutOff = numsamples // 2
    training = samples[: cutOff]
    learning = samples[cutOff:]
    groups = dict([(0,[]), (1,[])])
    _dict = dict([(a,0) for a in ambientSet])
    return
    # then we do the grouping.


def clash_uniform_test(samples, threshold = 1):
    _dict = {}
    try:
        m = samples[0].parent().cardinality()
    except:
        raise ValueError('Parent problems. %s'%samples[0].parent())
    N = len(samples)
    for s in samples:
        if s in _dict:
            _dict[s] += 1
        else:
            _dict[s] = 1
    K1 = sum([a for a in _dict.values() if a == 1])
    K1u = float(N*((m-1)/m)**(N-1))
    stat = K1u - K1
    print 'stat = %s'%stat
    if stat >= threshold:
        print 'non-uniform'
        return False
    else:
        print 'uniform'
        return True

def selecting_bins(q, r, numsamples):
    """
    finite field of size q**r.
    """
    v = (q**r - 1).divisors()
    ideal_bins = numsamples // 5
    for i in range(len(v)):
        if v[i] > ideal_bins:
            return v[i-1]
    return q**r

def _my_dot_product(lst1,lst2):
    return sum([a*b for a,b in zip(lst1,lst2)])

def _my_list_diff(lst1,lst2):
    return [a-b for a,b in zip(lst1,lst2)]

def even_rounding(q,p,e):
    """
    e -- an integer between 0 and q-1
    we do the probablistic rounding such that U(F_q) is mapped to U(F_p)
    """
    try:
        e = ZZ(e)
    except:
        raise ValueError('e must be an integer between 0 and q-1.')
    r = ZZ(q//p)
    if e < p*r and e >= 0:
        return ZZ(Mod(e,p))
    else:
        return ZZ.random_element(0,p)

def basis_transform_matrix(v,K):
    """
    Go back and forth between integral basis and power basis for number fields.

    Given a basis v of K/QQ , where K = QQ(a), return the matrix T such that
    (v_1,\cdots, v_n)T = (1,a, ..., a^{n-1})
    """
    n = K.degree()
    if len(v) != n:
        raise ValueError
    a = K.gen()
    ff = a.coordinates_in_terms_of_powers()
    return Matrix(QQ,n,n,[ff(vi) for vi in v])**(-1)


#### part 2: uniformity test ####


def generate_uniform_samples(q,r, numsamples = 2000):
    if r > 1:
        F.<alpha> = GF(q**r, impl = 'pari_ffelt')
    else:
        F.<alpha> = GF(q)
    print 'generating uniform distribution for a comparison.'

    from itertools import product
    print 'number of samples used = %s'%numsamples
    sys.stdout.flush()

    _dict = dict([(t,0) for t in F])
    alpha = F.gen()
    result = []
    for i in range(numsamples):
        e = sum([ZZ.random_element(0,q)*alpha**ii for ii in range(r)])
        result.append(e)
        if Mod(i, 5000) == 0 and i > 0:
            print '%s samples done.'%i
            print('e = %s'%e)
            sys.stdout.flush()
    return result


# generic chi-square test.
def chisquare_test(hist_dict,bins = None,return_dict = False, alpha = None):
    """
    well, somehow divide the distribution into bins.
    """
    if isinstance(hist_dict,list):
        F = hist_dict[0].parent()
        _dict = dict([(aa, 0) for aa in F])
        for ss in hist_dict:
            _dict[ss] += 1
        hist_dict = _dict
    numkeys = len(hist_dict.keys())
    numsamples = sum(hist_dict.values())
    if bins is None:
        bins = numkeys
    else:
        bins = ZZ(min(bins, numkeys))
    print 'number of bins used = %s'%bins
    if alpha is None:
        alpha = 1 - float(1/(10*bins))
    #print 'alpha = %s'%alpha
    #print 'number of keys = %s'%numkeys
    E = float(numsamples/bins)
    #print 'Expected balls in each bin = %s'%E

    if bins <= 0:
        raise ValueError('number of bins must be positive.')
    if E < 0.1:
        raise ValueError('expected value in each bin is too small.')
    newdict = dict([(a,0) for a in range(bins)])
    quo = ZZ(numkeys//bins)
    rem = numkeys - quo*bins
    keys = hist_dict.keys()
    for i in range(quo*bins):
        newdict[ZZ(Mod(i,bins))] += hist_dict[keys[i]]
    for j in range(quo*bins, numkeys):
        newdict[ZZ.random_element(0,bins)] += hist_dict[keys[j]]

    chisquare = float(sum([(t-E)**2 for t in newdict.values()])/E);
    mu = bins-1
    sigma = float(sqrt(2*bins-2))

    T = RealDistribution('chisquared', mu)
    print('chisquare = %s'%chisquare)
    prob = T.cum_distribution_function(chisquare)
    #if chisquare < mu - mm*sigma or chisquare > mu + mm*sigma:
    if prob > alpha:
        print 'non-uniform'
        uniform =  False
    else:
        print 'uniform'
        uniform = True
    if return_dict:
        return uniform, newdict
    else:
        return uniform


def _another_t_matrix(s, prec = 100):
    """
    r, s -- the signature of the number field
    """
    C = ComplexField(prec)
    eyer = Matrix.identity(s)

    return C(1/2)*block_matrix([[eyer, eyer], [-C(I)*eyer, C(I)*eyer]])

    # small_block = Matrix(C, [[1/2,1/2], [1/(2*C(I)), - 1/(2*C(I))]])

    # return block_diagonal_matrix([eyer] + [small_block for _ in range(s)])

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




def general_embedding_matrix(v, K, prec = 100, canonical = False):
    """
    The embedding matrix of any vector over a number field K.
    Suppose sigma_1, ..., sigma_r and tau_1,..., tau_s are the
    real and complex places of K, return the matrix A
    whose i-th column is sigma_1(v[i]), ..., sigma_r(v[i]), Re(tau_1), Im(tau_1), ..., Re(tau_s), Im(tau_s).

    If the input is a single element a , then return the embedding matrix for power basis
    [1,a,..., a^n-1], where n = [K:QQ].

    canonical: returns the canonical embedding matrix [\sigma_i(a_j)] (no attempt to return a real matrix).
    """
    n = K.degree()
    a = K.gen()
    if not isinstance(v,list):
        v = [v**i for i in range(n)]


    all_embs = K.complex_embeddings(prec = prec)
    C = ComplexField(prec)
    # if it is canonical embedding, then we return a complex matrix, which is Okay.
    if canonical:
        return Matrix(C,n,n, [phi(a) for phi in all_embs for a in v]).transpose()

    rr_embs = [phi for phi in all_embs if phi(a).imag() == 0]
    cc_embs = [phi for phi in all_embs if phi(a).imag() >0]

    r,s = len(rr_embs), len(cc_embs)
    if r+2*s != n:
        raise ValueError

    A2 = Matrix(C,n,s, [[phi(a) for phi in cc_embs] for a in v]).transpose()
    A2 = A2.stack(A2.conjugate())
    t = _another_t_matrix(s,prec = prec) # t is 2s*2s matrix, A2 is 2s*n.
    B2 = _real_part(t*A2)
    if r > 0:
        A1 = Matrix(n,r,[[phi(a) for phi in rr_embs] for a in v]).transpose()
        if s > 0:
            return _real_part(A1).stack(B2)
        else:
            return _real_part(A1)
    else: return B2
        # somehow this is not a Vandermont matrix. So I should do something...


def chisquare_quality(delta,N, c = 5):
    """
    c -- the ratio of sample to ambient set size.
    c = M/N.
    delta -- the statistical distance.

    Used in chisquare.pdf
    """
    T = RealDistribution('gaussian', 1)
    w= T.cum_distribution_function_inv(1 - 1/(N*20))
    _lambda = 4*c*N*delta**2
    #print 'w =  %s'%w
    ss = RR(w*sqrt(2*(N-1)) - _lambda)/RR(sqrt(2*(N-1) + 4*_lambda))
    #print 'ss = %s'%ss
    return 0.904*(1 - T.cum_distribution_function(ss))


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

def finite_field_coerce(F,f,c):
    """
    f: F -> K some embeddin.
    c: some element of K which is in f(F).
    returns the unique element b in F such that
    f(b) = c.
    """
    pp = c.minpoly()
    for b, _ in F[x](pp).roots():
        if f(b) == c:
            return b
    raise ValueError('no roots found.')

def finite_cyclo_traces(m,q,ilst,H,deg =1):
    """
    trace of zeta_m^i under H in the finite field F_q^r.
    where r is the order of q in Zm^* (by CFT).

    """
    if not q.is_prime():
        raise ValueError

    if deg > 1:
        F.<alpha> = GF(q^deg, impl = 'pari_ffelt')
    else:
        F.<alpha> = GF(q)
    Zm = Integers(m)
    #if Mod(q,m) not in H:
    #    raise ValueError('we require q mod m to be in H.')
    r = Zm(q).multiplicative_order()

    if Mod(r,deg):
        raise ValueError('wtf?')
    rel_deg = r// deg
    K.<beta>, f = F.extension(rel_deg, map = True)
    g = K.multiplicative_generator()
    d = ZZ((q^r-1)/m)
    f_zetam = g^d

    return [finite_field_coerce(F,f,(sum([f_zetam**(i*ZZ(h))for h in H]))) for i in ilst]



def simulate_run(D,q,vecs, numsamples):
    """
    D -- a lattice sampler
    q -- a prime
    vecs -- the basis vectors of D modulo q
    """
    Fq = GF(q)
    vecs = [Fq(a) for a in vecs]

    errs = []
    max_err = 0
    #count_zero = 0
    for i in range(numsamples):
        e = D()
        b = sum([a*e for a,e in zip(vecs,e)])
        errs.append(b)
        #b = ZZ(b) if b < q//2 else q - ZZ(b)

        #if b ==0:
        #    count_zero +=1
        #print abs(b)
    errs = Set(errs)

    rederrs = reduce_roots(list(errs),q)
    print 'max err = %s'%max([abs(a) for a in rederrs])
    # print 'errs = %s'%errs

    ratio  = len(errs) / q

    #print 'number of zeros = %s'%count_zero
    #print 'Done!'
    return ratio


def reduce_roots(v,q):
    return  [ZZ(a) if ZZ(a) <=  q//2 else ZZ(a) - q for a in v]


#############################
# some ELOS functions
#############################

# this is evaluation at 1.


def q_over_r_test(D, q, r = 2.5, numsamples = 20, numtrials =10):
    for i in range(numtrials):
        failed = False
        for j in range(numsamples):
            e =  sum(D())
            ez = ZZ(e)
            red_ez = ez if ez < q//2 else ez - q
            if abs(red_ez) > RR(q / r):
                print 'Failed', red_ez
                failed = True
                break
        if not failed:
            print 'Succeeded'

def q_over_r_test_uniform(q, r = 2.5, numsamples = 20, numtrials =10):
    for i in range(numtrials):
        failed = False
        for j in range(numsamples):
            e =  ZZ.random_element(-q//2,q//2+1, distribution = 'uniform')
            if abs(e) > RR(q / r):
                print 'Failed', e
                failed = True
                break
        if not failed:
            print 'Succeeded'