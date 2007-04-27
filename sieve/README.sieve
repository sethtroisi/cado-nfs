Compile with

gcc -std=c99 -O0 -g -o sieve -Wall -Wextra sieve.c -lm

Make a factor base with, for example,

echo "rootfind(x^6-3, 1000)" | gp -q rootfind.gp >roots.3_317.txt

Run the siever with

./sieve -v -fb roots.3_317.txt -thres 30 -leading 840 -1000 1000 210 220

This will output the values of $a, b$, -1000 <= a <= 1000, 210 <= b <= 220,
where the sum of logs of primes in the factor base is >= 30.
