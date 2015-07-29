# This is specifically searching for subfields of prime cyclotomic fields.
# Why is it not completely covered by the sub-cyclotomics functionality?
# One reason is we can specify the degree and search for primes.
# so it is a different story.
# Also we do not need the complication of finding a nice basis (it is just there, by simply taking a trace).

# Prime Cyclotomics

import sys

load('mega.sage')

d = 8; max_modulus = 5000 # the maximum value of q we consider.

plst = [e for e in primes(2*d+1, 5000) if Mod(e,d) == 1]; # the candidates p.

print 'the parameters are: d = %s, max_modulus = %s'%(d,max_modulus)
print

print 'We have %s prime candidates = %s'%(len(plst), plst)
sys.stdout.flush()

for p in plst:
    print 'working with p = %s'%p
    print
    sys.stdout.flush()

    L.<z> = CyclotomicField(p); n = L.degree()
    a = GF(p).multiplicative_generator(); a
    dprime = n//d;
    print 'dprime = %s'%dprime
    print
    t = cputime()
    c = sum([z**(b) for b in [a**(i*d) for i in range(dprime)]]); f = c.minpoly();
    print 'field generation took %s seconds.'%cputime(t)
    print 'f = %s'%f
    print
    sys.stdout.flush()

    K.<y> = NumberField(f)
    y = K.gen(); v = y.galois_conjugates(K)

    dK = K.disc(); nK = K.degree()
    adjK = RR(sqrt(p))/RR(abs(dK)**(1/(2*nK))) # changed! Well I don't know if this helps, but it should!
    #if Mod(dprime,2) != 0: # totally complex
    #    adjK *= RR(sqrt(2))
    print 'adjustment = %s'%adjK
    print

    t = cputime()
    # this replaces the split prime functionality (the split prime could also be easier, but testing if q belongs to H)
    qlst = [q for q in primes(max_modulus) if gcd(q,p) == 1 and GF(p)(q).multiplicative_order().divides(dprime)]
    # print 'moduli generation took %s seconds.'%cputime(t)
    sys.stdout.flush()

    opt_q, opt_ratio = search_for_q(v,K,qlst, h = dprime)
    print 'The best (unadjusted) adjusted ratio is (%s) %s, with q = %s'%(opt_ratio, opt_ratio/adjK, opt_q)
    print '*'*20
    sys.stdout.flush()

print 'Done!'