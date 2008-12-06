s:=1.2;            
p:=10;
while p lt 2^64 do
    for i in [1..20] do
        p:=NextPrime(p);
        d:=Random([2..7]);
        coeffs:=[Random(GF(p)):i in [0..d]];
        F:=PolynomialRing(GF(p))!coeffs;
        printf "in %o", p;
        for c in Reverse(coeffs) do printf " %o", c; end for;
        printf "\n";
        r:=Sort([x[1]: x in Roots(F) ]);
        printf "out";
        for c in r do printf " %o", c; end for;
        printf "\n";
    end for;
    // if ever someday we venture beyond 2^64, then we'd better speed up the
    // growth of coefficients.
    if p gt 2^64 then s:=5; end if;
    p := Ceiling(p*s);
end while;
