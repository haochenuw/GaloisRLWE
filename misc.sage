load('sampling-subcyclotomics.sage')

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
    return  [ZZ(a) if ZZ(a) < q//2 else ZZ(a) - q for a in v]
