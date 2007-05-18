Compile with

gcc -DWANT_ASSERT -I ../ -O0 -g -o sieve -Wall -Wextra sieve.c -lm -lgmp

(Note: linking GMP statically with /usr/lib/libgmp.a for example may yield a
speedup up to 7%.)

Make a factor base with, for example,

echo "rootfind(5017309194362523*x^4 -1406293661386525*x^3 -1131155401311965*x^2 +4737694118287353*x -3415040824020545, 1000000)" | gp -q -p 1000000 rootfind.gp >c80.roots.txt

Run the siever with

./sieve -poly c80.poly -fb c80.roots.txt -1000000 1000000 2001 2101

where c80.poly is the polynomial.
