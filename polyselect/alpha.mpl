get_alpha := proc(f, B) local s, p, e, q, disc;
   s := 0;
   p := 2;
   disc := discrim(f,x);
   while p <= B do
      if disc mod p = 0 then
         lprint(p," divides disc(f)");
         e := est_valuation (f, p, 0, 999);
         s := s + evalf((1/(p-1)-e)*log(p));
      else
         q := 0;
         for e in Roots(f) mod p do
            q:=q+1;
         od;
         s := s + evalf((1-q*p/(p+1))*log(p)/(p-1))
      fi;
      p := nextprime(p);
   od;
   s;
end:

# p-valuation of integer N
valuation := proc(N, p) local v;
   if N=0 then ERROR("input is zero") fi;
   for v from 0 while N mod p^(v+1) = 0 do od;
   v
end:

# estimate average p-valuation of polynomial f in x, from x0 to x1
est_valuation := proc(f, p, x0, x1) local s, v;
   s := 0;
   for v from x0 to x1 do s:=s+valuation(subs(x=v,f),p) od;
   1.0*s/(x1+1-x0)
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


