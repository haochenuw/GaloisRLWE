# Some functions for modulus switching purposes.

import sys

def l1_norm(a):
    a = list(a)
    return sum([abs(t) for t in a])

def loo_norm(a):
    a = list(a)
    return max([abs(t) for t in a])

def general_test_modswitch(D, v, p,q, trans_mat, method = 'maxerr', numsamples = None, _eval_value = 1, secret = None, norm_bound = None):
    """
    D -- a MyLatticeSampler instance.

    For explanation, see ...

    There is a way to do this minkowski-freely.

    """
    #if secretNormBound is None:
    #    secretNormBound = float(p*2/D._degree)

    if method == 'maxerr':
        numsamples = 30
    elif method == 'chisquare':
        bins = p // 5 + 2
        numsamples = bins*10
    print 'number of samples = %s'%numsamples
    n = D._degree
    A = D.A
    colsum_vec = D.col_sum()


    if secret is None:
        s = D()[1]
    else:
        s = vector(secret)

    # compute and print the infinity norm of secret


    s_ = A*vector(s)
    s_lst = list(s)
    s_field = sum([si*vi for si, vi in zip(s_lst,v)]) # the field element
    # print 's = %s'%s_field
    l1s = l1_norm(s_)
    print '|s|_oo = %s'%loo_norm(s_).n()

    print '|s|_1 = %s'%l1s.n()
    print '|s|_2 = %s'%(s_.norm().n())

    alpha = QQ(p/q)

    errors = []
    beta = v[0].parent().gen()
    ff = beta.coordinates_in_terms_of_powers()

    TT = trans_mat
    jj, count = 0,0
    while jj < numsamples: # baby case: zero error.
        count += 1
        verbose('already have %s samples'%jj)
        verbose('processing sample %s...'%count)

        # get a uniform a, scale it, and produce the difference
        a_lst = [ZZ.random_element(0,q) for i in range(n)]
        a_ = A*vector(a_lst)

        a_field = _my_dot_product(a_lst,v) # the field element

        as_field = a_field*s_field
        as_lst = list(vector(ff(as_field))*TT)
        e_lst = list(zero_vector(n))
        b_lst = [ZZ(Mod(u1+u2,q)) for u1, u2 in zip(as_lst,e_lst)]

        alpha_a = [ai*alpha for ai in a_lst]


        rounded_alpha_a = list(D.babai(c = A*vector(alpha_a))[1] )  # an approximation of scaled_a.

        aprime_lst =  _my_list_diff(alpha_a, rounded_alpha_a)# a short vector in 1/q O_K.
        aprime_field = _my_dot_product(aprime_lst,v)# a short vector in 1/q O_K.

        # Select samples
        # Maybe l_1 is better?
        aprime_ = A*vector(aprime_lst)
        normaprime = aprime_.norm()
        verbose('norm of a\'= %s'%float(normaprime))

        if norm_bound is not None and normaprime > norm_bound:
            continue
            # move on to the next sample
        else:
            jj +=1

        # round b
        alpha_b = [ai*alpha for ai in b_lst]
        rounded_alpha_b = list(D.babai(c = A*vector(alpha_b))[1]) # an approximation of scaled_a.
        bprime_lst = _my_list_diff(alpha_b, rounded_alpha_b)
        bprime_field = _my_dot_product(bprime_lst,v)# a short vector in 1/q O_K.
        #verbose('norm of b\'= %s'%float((A*vector(bprime_lst)).norm()))

        # compute e
        eprime_field = bprime_field - aprime_field*s_field
        eprime_lst =  list(vector(ff(eprime_field))*TT)

        # verbose('l2 norm of e\'= %s'%float((A*vector(eprime_lst)).norm()))


        eprime_lst = [ZZ(ei) for ei in eprime_lst]
        # naive_eprime_lst = [ZZ(ei) for ei in naive_eprime_lst]
        # evaluation attack
        eprime_at_one = Mod(ZZ[x](eprime_lst)(x = _eval_value), p)
        verbose('e\'(1) = %s'%eprime_at_one)
        errors.append(eprime_at_one)


    # replace this by a chi squere test.
    print 'samples actually used = %s'%count
    errors_red = [ZZ(e) if e <= p//2 else ZZ(e) - p for e in errors]

    if method == 'maxerr':
        return RDF(max([abs(e) for e in errors_red])*2/p)
    elif method == 'chisquare':
        return chi_square_test(errors_red,bins,p)
    else:
        raise NotImplementedError


def chi_square_test(elst,bins,p):
    """
    p -- the modulus
    bins - number of bins
    elst -- the elements list
    """
    numsamples = len(elst)
    _dict = dict([(a,0) for a in range(bins)])
    elst = [ZZ(Mod(ZZ(e),p)) for e in elst]
    for e in elst:
        _dict[floor(QQ(bins/p)*e)] += 1
    E = float(numsamples/bins)
    chisquare = float(sum([(t-E)**2 for t in _dict.values()])/E);
    mu = bins-1
    sigma = float(sqrt(2*bins-2))
    print 'degree of freedome = %s'%mu
    print 'chisquare value = %s'%chisquare
    if chisquare < mu - 2.5*sigma or chisquare > mu + 2.5*sigma:
        print 'non-uniform'
        return True
    else:
        print 'uniform'
        return False

def test_modswitch(D,f,p,q, numsamples = 30, _eval_value = 1, secret = None, naive = False, babai = False, norm_bound = 3.0):
    """
    D -- a MyLatticeSampler instance.

    For explanation, see ...

    There is a way to do this minkowski-freely, but that may not be what we want to do.
    """
    n = D._degree
    P.<x> = PolynomialRing(QQ)
    PQ.<y> = PolynomialQuotientRing(P,P(f))
    # print PQ
    # generate a secret.
    if secret is None:
        s = PQ(list(D()))
    else:
        s = PQ(secret)
    print 's = %s'%s
    print

    A = D.A

    alpha = QQ(p/q)

    errors = []
    i = 0
    count = 0
    while i < numsamples : # baby case: zero error.
        verbose('processing sample %s...'%count)
        count += 1
        sys.stdout.flush()

        # get a uniform a, scale it.
        a = PQ([ZZ.random_element(0,q) for i in range(n)])
        b = PQ([ZZ(Mod(bi,q)) for bi in list(a*s)])
        alpha_a = a*alpha
        alpha_b = b*alpha

        # round it and compute the difference
        if babai:
            rounded_alpha_a = PQ(list(D.babai(c = A*vector(alpha_a))[1])) # an approximation of scaled_a.
            rounded_alpha_b = PQ(list(D.babai(c = A*vector(alpha_b))[1]))
        else:
            rounded_alpha_a = PQ(list(D(c = A*vector(alpha_a))[1])) # an approximation of scaled_a.
            rounded_alpha_b = PQ(list(D(c = A*vector(alpha_b))[1]))

        #verbose('rounded alphaa = %s'%rounded_alpha_a)

        aprime = alpha_a - rounded_alpha_a# a short vector in 1/q O_K.
        bprime = alpha_b - rounded_alpha_b # a short vector in 1/q O_K.

        #verbose('a\'= %s'%aprime)

        normaprime = (A*vector(aprime)).norm()
        if normaprime > norm_bound:
            continue
        i += 1
        verbose('norm of a\'= %s'%float(normaprime))

        # do a naive rounding (result is not good, which is good! since we don't

        #naive_rounded_alpha_a = PQ([round(ai) for ai in list(alpha_a)])
        #naive_aprime = alpha_a - naive_rounded_alpha_a# a short vector in 1/q O_K.
        #verbose('naive a\' = %s'%naive_aprime)

        #naive_rounded_alpha_b = PQ([round(ai) for ai in list(alpha_b)])
        # naive_bprime = alpha_b - naive_rounded_alpha_b# a short vector in 1/q O_K.

        # we compute the rounding error.
        eprime = bprime - aprime*s
        eprime = PQ([ZZ(Mod(ei,p)) for ei in list(eprime)])

        # naive_eprime = naive_bprime - naive_aprime*s
        verbose('e\' = %s'%list(eprime))
        # verbose('norm of e\'= %s'%float((A*vector(eprime)).norm()))

        eprime_lst = list(eprime)
        #naive_eprime_lst = list(naive_eprime)

        try:
            eprime_lst = [ZZ(ei) for ei in eprime_lst]
            #naive_eprime_lst = [ZZ(ei) for ei in naive_eprime_lst]

            # evaluation attack
            eprime_at_one = Mod(ZZ[x](eprime_lst)(x = _eval_value), p)
            verbose('e\'(1) = %s'%eprime_at_one)

            #naive_eprime_at_one = Mod(ZZ[x](naive_eprime_lst)(x = _eval_value), p)
            if not naive:
                errors.append(eprime_at_one)
            else:
                errors.append(naive_eprime_at_one)
        except:
            # got non-integral results, debug
            raise ValueError('not all e_i are integers. Please debug.')

    errors_red = [ZZ(e) if e <= p//2 else ZZ(e) - p for e in errors]
    return RDF(max([abs(e) for e in errors_red])*2/p)


def embedded_l_oo_norm(a):
    return max([abs(ac) for ac in a.complex_embeddings()])



def random_elt(n, q):
    return ([ZZ.random_element(0,q) for _ in range(n)])


def gen_sample(secret):
    e = O(D())
    a = random_elt(q)
    print 'a = %s'%a
    print 'b = %s'%(a*secret + e)
    return a,b


def naive_round(lst):
    naive_round = [round(bi) for bi in lst];
    frac_part = [s-t for s,t in zip(lst,naive_round)]
    return frac_part, naive_round

def l2_norm_cyclotomic(lst,p):
    assert len(lst) == p-1
    return sqrt(RR(p*sum([a**2 for a in lst])) - RR(sum(lst))**2)

# K.<zeta> = CyclotomicField(p);

def _to_field(quotientringelt,K):
    zeta = K.gen()
    g = quotientringelt
    return sage_eval(str(g), locals ={'z':zeta})