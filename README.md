# GaloisRLWE

These code assists a search for vulnerable Galois RLWE instances.

#### sub-cyclotomics.sage

#### MyLatticeSampler.sage

- Modifies the current implementation of discrete Gaussian lattice sampler in sage (bug fixes + optimization).

- Can sample from a discrete Gaussian on a lattice in R^n, given as the span of the columns of a n-by-n matrix.

- Can perform Babai's nearest plane algorithm.

- Preprocessing: we do a LLL/BKZ reduction on the input, to optimize the above two functions.


#### SubCyCSampler.sage
- Provides a class representing a RLWE instance from a subfield of the m-th cyclotomic field, where m is odd and square free.

- Can generate RLWE samples.

- Can compute primes of a given degree.

- Compute the canonical normal integral basis of O_K, and its embedding matrix (real and complex).

- Can compute the numerical discriminant of the ring of integers O_K.

- Given a prime q, can compute the image of the canonical basis modulo a prime ideal above q.

- Depends on: MyLatticeSampler.sage, misc.sage

#### misc.sage

- Chi-square test for samples in some finite field.

- Small field chi-square test.

- Generate uniform samples in finite fields.

- Some other miscellaneous functions.

- Depends on: None.


