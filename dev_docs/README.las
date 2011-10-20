This is the README file for the lattice siever (in construction).

The lattice siever (las) takes the following arguments/options:

-checknorms    we want to completely factor the leftover norms
-I i           the sieve region used is of width 2^i. The larger i is,
               the more relations we get by special-q, but the harder it is.
               For gnfs numbers of 120-130 digits, one might want i=12, which
               is the default. For 155 digits, one may use i=13 or i=14.
-poly xxx.poly use polynomial xxx.poly (same format as the line siever)
-fb xxx.roots  use factor base xxx.roots (same format as the line siever)
-q0 nnn        left bound of special-q range
-q1 nnn        right bound of special-q range (the last special-q is q1-1,
               assuming q1-1 is prime); if missing, only q0 is sieved
-rho r         Assuming q0 is prime, sieve only the ideal (q0,r). This is
               mainly for debugging purposes.

The special-q range [q0, q1[ is usually taken in the large prime region
[alim+1, 2^lpba[ from the xxx.poly file, but this is not required. For
example GGNFS starts from alim/2. However we need q1 <= 2^lpba.
