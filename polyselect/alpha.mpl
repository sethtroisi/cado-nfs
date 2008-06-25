skew := proc(f, x, s0) local absf, g, i, S;
   absf := add(abs(coeff(f,x,i))*x^i, i=0..degree(f,x));
   g := expand(subs(x=S, absf)/S^(degree(f)/2));
   g := diff(g,S);
   fsolve(g, S=s0)
end:

# takes also into account the linear polynomial x-m
skew2 := proc(f, x, m, s0) local absf, g, i, S;
   absf := add(abs(coeff(f,x,i))*x^i, i=0..degree(f,x));
   g := expand(subs(x=S, absf)/S^(degree(f)/2))*(m/sqrt(S)+sqrt(S));
   g := diff(g,S);
   fsolve(g, S=s0)
end:

get_alpha := proc(f, B) local s, p, e, q, disc;
   s := 0;
   p := 2;
   disc := discrim(f,x);
   while p <= B do
      if disc mod p = 0 then
         # e := est_valuation (f, p, 1, 99, 1, 99);
         e := val (f, p);
#         lprint(p," divides disc(f)", evalf(e,3));
         s := s + evalf((1/(p-1)-e)*log(p));
      else
         q := 0;
         for e in Roots(f) mod p do
            q:=q+1
         od;
         s := s + evalf((1-q*p/(p+1))*log(p)/(p-1))
      fi;
      p := nextprime(p);
   od;
   s;
end:

# returns the smallest value of alpha for p1 <= p <= p2,
# with a degree-d polynomial
min_alpha := proc(d,p1,p2) local p, s;
   s := 0;
   p := nextprime(p1-1);
   while p <= p2 do
      s := s + evalf((1-min(d,p)*p/(p+1))*log(p)/(p-1));
      p := nextprime(p);
   od;
   s
end:

# p-valuation of integer N
valuation := proc(N, p) local v;
   if N=0 then ERROR("input is zero") fi;
   for v from 0 while N mod p^(v+1) = 0 do od;
   v
end:

# estimate average p-valuation of polynomial f in x, from x0 to x1
est_valuation := proc(f0, p, x0, x1, y0, y1) local s, v, w, f, n;
   f := expand(subs(x=x/y,f0)*y^degree(f0));
   s := 0;
   n := 0;
   for v from x0 to x1 do for w from y0 to y1 do
      if igcd(v,w)=1 then
         s:=s+valuation(subs(x=v,y=w,f),p);
         n:=n+1;
      fi
   od od;
   1.0*s/n
end:

# resultant(P,Q,x), with computations with rationals
res := proc(P,Q,x) local q, R;
   q := lcoeff(Q,x);
   if degree(Q)=0 then q^degree(P,x)
   elif degree(Q)>degree(P) then procname(Q,P,x)
   else R:=rem(P,Q,x); q^(degree(P,x)-degree(R,x))*procname(Q,R,x)
   fi
end:

# fraction-free resultant computation
res2 := proc(P0,Q0,x) local P, Q, q, m, d, s;
   P:=P0;
   Q:=Q0;
   if degree(Q,x)>degree(P,x) then procname(Q,P,x)
   else # deg(P) >= deg(Q)
      m := 1;
      d := 1;
      while degree(Q,x)>0 do
         q := lcoeff(Q,x);
         s := degree(P); # multiply by q^deg(P)
         while degree(P,x)>=degree(Q,x) do
            P:=q*P;
            s := s - degree(Q); # divide by q^deg(Q)
            P:=expand(P-lcoeff(P,x)/lcoeff(Q,x)*x^(degree(P,x)-degree(Q,x))*Q);
         od;
         s := s - degree(P); # divide by q^deg(P mod Q)
         if s>0 then m := m*q^s else d:=d*q^(-s) fi;
         P,Q:=Q,P;
      od;
      m/d*Q^degree(P,x);
   fi;
end:

# the following is from Guillaume Hanrot
# determines average valuation for a prime dividing disc(P)
val := proc(P,p)
   (val0(P, p) * p + val0(expand(subs(x=1/(p*x),P)*(p*x)^degree(P)), p))/(p+1)
end:

val0 := proc(P0,p) local v, r, P, Q, P2;
   v := valuation (icontent(P0), p);
   P := P0/p^v;
   Q := diff(P,x);
   for r in Roots(P) mod p do
      if subs(x=r[1],Q) mod p <> 0 then v := v + 1/(p-1)
      else
         P2 := expand(subs(x=p*x+r[1], P));
         v := v + procname(P2, p)/p
      fi
   od;
   v
end:

##############################################################################

# the following procedure decomposes N in base (m,p)
# cf Lemma 2.1 of Kleinjung's paper "On polynomial selection for the
# general number field sieve"
# Example:
# N:=10941738641570527421809707322040357612003732945449205990913842131476349984288934784717997257891267332497625752899781833797076537244027146743531593354333897
# Lemma21(N, 5, 102406, 1197773395291, 639369899891975567556774664501);
Lemma21 := proc(N, d, ad, p, m) local a, r, i, mu;
   a := table();
   r := N;
   a[d] := ad;
   for i from d-1 by -1 to 0 do
      r := (r - a[i+1]*m^(i+1))/p;
      if not type(r, integer) then ERROR("r is not an integer") fi;
      # find -p*m^i/2 <= mu < p*m^i/2 such that r + mu is divisible by m^i
      # and mu = 0 (mod p), thus mu=lambda*p where lambda = -r/p mod m^i
      mu := (-r/p mod (m^i));
      if mu >= m^i/2 then mu:=mu-m^i fi;
      mu := mu*p;
      a[i] := (r + mu) / m^i;
      if not type(a[i], integer) then ERROR("a[i] is not an integer") fi;
   od;
   add(a[i]*x^i, i=0..d)
end:

# compute the sup-norm as in Definition 3.1 of Kleinjung's paper
norme := proc(f) local d, i, j, k, s, ai, ok, n, minn, mins;
   minn := infinity;
   d := degree(f);
   for i from 0 to d do
      ai := abs(coeff(f,x,i));
      if ai=0 then next fi;
      for j from i+1 to d do
         # solve |a_i|*s^(i-d/2) = |a_j|*s^(j-d/2)
         s := evalf(abs(coeff(f,x,j)/ai)^(1/(i-j)));
         n := evalf(ai*s^(i-d/2));
         ok := true;
         for k in {$0..d} minus {i,j} while ok do
            ok := evalb(abs(coeff(f,x,k))*s^(k-d/2) < n);
         od;
         if ok and n<minn then minn:=n; mins:=s fi;
      od
   od;
   mins, minn
end:

sup_norm := proc(f, s) local i;
   seq(evalf(abs(coeff(f,x,i))*s^(i-degree(f)/2)), i=0..degree(f))
end:
