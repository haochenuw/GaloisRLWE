# some misc functions for use.

from sage.libs.fplll.fplll import FP_LLL

import sys

def _my_dot_product(lst1,lst2):
    return sum([a*b for a,b in zip(lst1,lst2)])

def _my_list_diff(lst1,lst2):
    return [a-b for a,b in zip(lst1,lst2)]


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



def test_elos_uniform_with_samples(errors, vq, q, bins = None, std_multiplier = 3):
    """
    a uniform test when we lready has samples. I mean errors, we assume that
    the errors lie in some finite field F_q^r.

    """
    F = vq[0].parent()
    r = F.degree()
    if bins is None:
        bins = ZZ(len(errors)//5)
    smallbins = ZZ(floor(bins**(1/r)))
    print 'small bins = %s'%smallbins
    print 'std multiplier = %s'%std_multiplier
    bins = smallbins**r
    print 'degree of finite field = %s'%r
    from itertools import product

    print 'degree of freedom = %s'%(bins-1)
    numsamples = len(errors)
    print 'number of samples used = %s'%numsamples
    sys.stdout.flush()

    _dict = dict([(tuple(a),0) for a in product(range(smallbins), repeat = r)])


    for i in range(numsamples):
        if Mod(i, 10000) == 0 and i > 0:
            print '10000 samples done'
            sys.stdout.flush()
        error = errors[i]
        e = F(sum([a*b for a,b in zip(error,vq)]))
        e_polylst = [Mod(tt, q) for tt in list(e.polynomial())]
        # pad with zero
        lene = len(e_polylst)
        if lene < r:
            e_polylst += [0 for _ in range(r-lene)]
        verbose('e_poly = %s'%e_polylst)
        _key = tuple([floor(QQ(smallbins/q)*ZZ(e)) for e in e_polylst])
        #verbose('key = %s'%key)
        _dict[_key] += 1

    E = float(numsamples/bins)
    chisquare = float(sum([(t-E)**2 for t in _dict.values()])/E);
    mu = bins-1
    sigma = float(sqrt(2*bins-2))

    # print 'dictionary = %s'%_dict
    print 'chisquare value = %s'%chisquare
    mm = std_multiplier
    if chisquare < mu - mm*sigma or chisquare > mu + mm*sigma:
        print 'non-uniform'
        return True
    else:
        print 'uniform'
        return False


def test_elos_uniform(D,vq,q, numsamples = None, bins = None, sanity_check = False):
    """
    D -- a lattice sampler
    vq -- roots mod q
    q -- the modulus
    """
    if bins is None:
        bins = abs(ZZ(q//100)) + 2

    print 'degree of freedom = %s'%(bins-1)
    if numsamples is None:
        numsamples = bins*5
    print 'number of samples used = %s'%numsamples
    sys.stdout.flush()
    _dict = dict([(a,0) for a in range(bins)])
    for i in range(numsamples):
        if Mod(i, 500) == 0:
            print '500 samples generated'
            sys.stdout.flush()
        if sanity_check:
            e = ZZ.random_element(0,q)
        else:
            error = D()
            verbose('error = %s'%error[1])
            verbose('l_2 norm of error = %s'%error[0].norm())
            e = Mod(sum([a*b for a,b in zip(error[1],vq)]),q)
        _dict[floor(QQ(bins/q)*ZZ(e))] += 1


    E = float(numsamples/bins)
    chisquare = float(sum([(t-E)**2 for t in _dict.values()])/E);
    # T = RealDistribution('chisquared', bins-1);
    mu = bins-1
    sigma = float(sqrt(2*bins-2))

    #print 'dictionary = %s'%_dict
    print 'chisquare value = %s'%chisquare
    if chisquare < mu - 2.5*sigma or chisquare > mu + 2.5*sigma:
        print 'non-uniform'
        return True
    else:
        print 'uniform'
        return False





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

    """
    else:
        for i in range(n):
            vi = v[i]
            coli_real = [sigma(vi) for sigma in real_places]
            coli_complex = []
            for tau in complex_places:
                tauvi = tau(v[i])
                coli_complex += [tauvi.real_part(),tauvi.imag_part()]
            coli = coli_real+coli_complex
            result.append(coli)
    return _real_part(t*Matrix(n,n,result).transpose())
    """

def simulate_uniform(q, samples = 20):
    """
    simulate maximum error for uniform distribution.
    """
    max_abs = 0
    for i in range(samples):
        a = ZZ.random_element(-(q//2),q//2+1,distribution = 'uniform')
        max_abs = max(abs(a), max_abs)
    #print 'max = %s'%max_abs
    return RR(max_abs*2/q)


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