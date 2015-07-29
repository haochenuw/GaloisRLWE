# This is the main code for searching for subfields of odd squarefree cyclotomics


# Composite Cyclotomics


import sys

load('mega.sage')

# these are the parameters that can change
m = 13*11*3*5; min_degree = 10; max_degree = 10; max_modulus = 1000;

phim = euler_phi(m)

print 'm, degree, mindegree, max_modulus = (%s, %s, %s, %s)'%(m,phim, min_degree, max_modulus)

Hlst  = [H for H in AllSubgroups(m) if phim//H.order >= min_degree and phim//H.order <= max_degree]


print 'Found %s groups to work on.'%len(Hlst)
print
sys.stdout.flush()


for H in Hlst:

    overall_timer = cputime() #overalltimeer

    # Step 1: compute the elements of H.
    t = cputime()
    print 'Working with H = %s'%H

    n = phim // H.order
    print 'Number field has degree %s'%n

    H1 = H.compute_elements();
    cosets = H.cosets()

    # Step 2: compute representatives for 'good basis', and the ambient cyclotomic field.
    print 'computing a good basis...'
    rep_dict = good_basis_composite(m,H1);
    print 'basis reps computed'
    print 'First step : %s'%cputime(t)
    sys.stdout.flush()

    # Step 3: get the minpolys of the good basis vectors.
    #t = cputime()
    # this part can be optimized? how do we guarantee a generator?

    # where do we need this K? except for discriminant compuation
    #K = construct_K(rep_dict,L,H1)
    #print 'field generation : %s'%cputime(t) # very slow: so we do not want this
    #sys.stdout.flush()

    # this is commented out because field generation is too slow at this point so we just wanted
    # to skip it first: luckily we do not need to know its defining polynomial for the analysis.
    #print 'Found number field of degree %s and defining polynomial %s.'%(K.degree(), K.defining_polynomial())
    ##print
    #sys.stdout.flush()
    dK = H.discriminant() # Because we can compute the field discriminant just from H using Hasse's theorem.

    # Step 4: compute the adjustment factor.
    #adjK = adj_factor(m, H1,dK)
    #print 'adjustment factor = %s'%adjK
    #print
    adjK = RR(abs(dK)**(1/(2*n)))
    print 'adjustment factor = %s'%adjK
    sys.stdout.flush()


    # Step 5: generating stuff (very slow. Want to avoid.)
    #t = cputime()
    #elt_dict = construct_basis_in_K(res,L,K)
    #basis_elts, basis_multiplicities = elt_dict.keys(),elt_dict.values()

    #basis_elts = [(z**rep).trace(K) for rep in res.keys()] # this operation is very slow. Is it now faster?
    #basis_multiplicities = res.values()

    #print 'elements coercion : %s'%cputime(t)
    #print 'multiplicities = %s'%basis_multiplicities
    #sys.stdout.flush()


    # Step 5.5: sanity check.
    #variances = rep_dict.values()
    #stds = [RR(sqrt(s)) for s in variances]
    #assert sum(variances) == m
    #sys.stdout.flush()

    # Step 6: generate a list of moduli that splits in K.
    qs = split_primes_new(m, H1, max_prime = max_modulus)
    print 'modulus candidates are : %s'%qs
    print
    sys.stdout.flush()

    # Step 7: loop over possible moduli to find a good ratio.
    dprime = len(H1)
    #scale = 0
    opt_q, min_ratio = 0, 1
    for q in qs:
        #print 'q = %s'%q
        t = cputime()
        # important!!!
        #reduced_vecs = finite_cyclo_traces(m,q,rep_dict.keys(), H1) # this may again be slow.
        reduced_vecs = finite_cyclo_traces(m,q, cosets, H1)
        # print 'the reduced vectors are %s'%reduced_vecs

        #okay, now I do not scale

        #_norm, scale, _min_roots = loo_norm(reduced_vecs,q)
        # print 'loo norm is %s with min_roots = %s'%(_norm, _min_roots)
        _norm = max([abs(a) for a in reduce_roots(reduced_vecs,q)])
        raw_ratio = _norm*2/q
        ratio = _norm*4/(q*adjK)
        # print 'q = %s, unadjusted ratio = %s'%(q, raw_ratio)
        #print cputime(t)
        sys.stdout.flush()

        if ratio < min_ratio:
            opt_q, min_ratio, scale = q, ratio, scale
    print '*'*25
    print 'n = %s, q = %s, adjusted ratio = %s, scale = %s'%(n, opt_q, min_ratio, scale)
    print 'time = %s'%cputime(overall_timer)
    print '*'*25
    sys.stdout.flush()

print 'Done!'



