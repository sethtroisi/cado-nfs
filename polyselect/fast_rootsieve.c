#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gmp.h"
#include "gmp_aux.h"
#include "cado.h"
#include "utils.h"
#include "aux.h"
#include "mod_ul_default.h"
#include "mod_ul.h"
#include "fast_rootsieve.h"
#include "macros.h" /* for ASSERT_ALWAYS */


/* The rootsieve : we compute the root property for many polynomials at once.
   We suppose the bounds U and V given.

________            _________________                  
___  __ \_____________  /__  ___/__(_)_______   ______ 
__  /_/ /  __ \  __ \  __/____ \__  /_  _ \_ | / /  _ \
_  _, _// /_/ / /_/ / /_ ____/ /_  / /  __/_ |/ //  __/
/_/ |_| \____/\____/\__/ /____/ /_/  \___/_____/ \___/ 


Here we import Emmanuel's Fast Rootsieve.

 */
/*
int main() {
  // We create example polynomials.
  mpz_poly f,g;
  int degf, degg;
  mpz_t * fco, * gco;   
  int i;
  
  degf = 0;
  degg = 1;

  fco = malloc((degf+1)*sizeof(mpz_t));
  gco = malloc((degg+1)*sizeof(mpz_t));

  mpz_init_set_str(fco[0],"-277565266791543881995216199713801103343120", 10);
  mpz_init_set_str(fco[1],"-18185779352088594356726018862434803054", 10);
  mpz_init_set_str(fco[2],"6525437261935989397109667371894785", 10);
  mpz_init_set_str(fco[3],"-46477854471727854271772677450", 10);
  mpz_init_set_str(fco[4],"-5006815697800138351796828", 10);
  mpz_init_set_str(fco[5],"1276509360768321888", 10);
  mpz_init_set_str(fco[6],"265482057982680", 10);

  mpz_init_set_str(gco[0],"-1291187456580021223163547791574810881",10);
  mpz_init_set_str(gco[1],"34661003550492501851445829",10);  
  
  mpz_poly_init(f,degf);
  mpz_poly_init(g,degg);

  //mpz_poly_set(f,fco);
  //mpz_poly_set(g,gco);

  for( i=0 ; i<= degf; i++)
    mpz_clear(fco[i]);
  
  for( i=0 ; i<= degg; i++)
    mpz_clear(gco[i]);

  free(fco);
  free(gco);

  mpz_poly_clear(f);
  mpz_poly_clear(g);

}*/

int main() {
  int i;	
  mpz_poly g;
  int degg;
  mpz_t * gco;
  
  degg = 1;
  //gco = malloc((degg+1)*sizeof(mpz_t));

  //mpz_init_set_str(gco[0],"-1291187456580021223163547791574810881",10);
  //mpz_init_set_str(gco[1],"34661003550492501851445829",10);

  mpz_poly_init(&g,degg);

  //for( i=0 ; i<= degg; i++)
  //  mpz_clear(gco[i]);


  mpz_poly_clear(g);
  
  //free(gco);


}


void rootsieve(mpz_poly * f,mpz_poly * g, int deg_f, int deg_g, long U, long V, long prime_bound ) {  
  
  // Counters
  long i;
  long p;

  // We prepare the dictionary
  rootsieve_dictionary rdict;
  rdict.B = prime_bound;
  rdict.global_offset = 0;
  rdict.umax = U;
  rdict.vmax = V;
  rdict.space = U*V;  
  rdict.sarr = malloc(rdict.space*sizeof(double));  

  for ( i=0 ; i<rdict.space ; i++ )
    rdict.sarr[i] = 0.0;  

  for ( p = 2 ; p <= prime_bound ; p += 1 + (1 & p))
    if (isprime(p))    
      rdict.global_offset += log(p)/(p-1);

  // This gives the exceptional values for u (once divided by g(l)^2)
  rdict.fdg_gdf = mpz_poly_sum(mpz_poly_prod(f,mpz_poly_derivative(g)),mpz_poly_prod(g,mpz_poly_derivative(f)));  
 
  for ( p = 2 ; p <= prime_bound ; p += 1 + (1 & p))
    if (isprime(p))         
      fill_alpha(rdict,p);  

  
  free(rdict.sarr);
  
}


int fill_alpha(rootsieve_dictionary rdict, long p) {
  

  //"""Fast root sieve modulo p and its powers. The rdict argument must
   // have been setup beforehand"""

  // Counters
  int i;
  long hits;

  long u0,v0,l0;
  double scale;
  int ld,m;
  mpz_poly ff,gg,fdg_gdf;

  ff = rdict.f;
  gg = rdict.g;
  scale = (log(p))*(double)p/(double)(p+1);
  
  l0 = 0;
  ld = m = 0;  

  rdict.cutoff = 1.0e-5;

  rdict.maxpow = (long) (16*log(2)/log(p));
  rdict.ppow   = malloc((rdict.maxpow+1)*sizeof(long));

  rdict.ppow[0] = 1;
  for( i=1 ; i<=rdict.maxpow; i++) 
    rdict.ppow[i] = rdict.ppow[i-1]*p;  

  rdict.gmodp = mpz_poly_coeff_reduction(rdict.g,rdict.ppow[1]);

  fdg_gdf = compose_reduce(rdict.fdg_gdf,l0,rdict,rdict.ppow[rdict.maxpow]);
  ff = compose_reduce(ff,l0,rdict,rdict.ppow[rdict.maxpow]);
  gg = compose_reduce(gg,l0,rdict,rdict.ppow[rdict.maxpow]);
  rdict.l0 = l0;
  hits = rotation_inner(rdict,p,ff,gg,fdg_gdf,u0,v0,l0,ld,m,0,scale,mpz_poly_identity());  
  return hits;  
}


long rotation_inner(rootsieve_dictionary rdict, long p,mpz_poly ff,mpz_poly gg,mpz_poly fdg_gdf, long u0, long v0, long l0,int ld,int m,long twist_v1,double scale,mpz_poly dphi) { 
  
  // Use your precomputations.
  // Update your precomputations

  //if len(hist) > 6:
  //      nhist=copy(hist)
  //      nhist.append([dphi,u0,v0,dphi(0),ld,m])
  //      print "Found tree of height %d" % len(nhist)
  //      for hi in nhist:
  //          print hi
  //

  if (scale < rdict.cutoff) 
    return 0;

  // Hit counter
  long thits,hits;
  thits = 0;
  hits = 0;

  // Abbreviations
  long igl;
  long pm,pm1,pmax;
  
  pm = rdict.ppow[m];
  pm1= rdict.ppow[m+1];
  pmax = rdict.ppow[rdict.maxpow];
 
  
  double scale1,scale2,scale3;
  int look_many_roots;
  mpz_poly minus_f_over_g_modp;
  
  // gg is never scaled down, and it is always evaluated at points where
  // the value is not zero mod p. So the valuation of g mod p never
  // vanishes.
  
  // This is always useful.  

  scale1 = scale/(p-1);
  scale2 = scale/p;
  scale3 = scale1-scale2;

  minus_f_over_g_modp = mpz_poly_coeff_reduction(ff,(unsigned long)rdict.ppow[1]);

  if (m > 0) {
    l0 = rdict.l0;    
    igl = rdict.igl;
    minus_f_over_g_modp = mpz_poly_coeff_product_si(minus_f_over_g_modp,igl);
    minus_f_over_g_modp = mpz_poly_coeff_reduction(minus_f_over_g_modp,(unsigned long)rdict.ppow[1]);
  }
  
  // Deciding whether we're in the general case or not has some
  // subtleties, unfortunately.
  look_many_roots = (!mpz_poly_is_constant(minus_f_over_g_modp)) ||  (m == 0 &&  ! mpz_poly_is_constant(rdict.gmodp));
  
  //nhist = copy(hist);
    
  //if (!look_many_roots)
  //  nhist.append([dphi,u0,v0,l0+p*dphi(0),ld,m]);
  
  long dv0,u1,v1;
  mpz_poly nfdg_gdf,nff,ngg,ndphi,tmp_poly,twist_f;
  long du,du0,new_twist_v1;

  if (look_many_roots) {
    // Then it's the typical case, as encountered for instance at the
    // beginning of the root sieve.
    // Each possible l value must be tried, and will give rise to
    // potentially intersecting lattices. 
    
    long l;
    
    for ( l=0 ; l<p ; l++) {
      
      //nhist=copy(hist);
      //nhist.append([dphi,u0,v0,l0+p*dphi(l),ld,m]);
      long minus_f_over_gmodp_l,gvl,igl;
      mpz_t mp_tmp;
      mpz_init(mp_tmp);
      mpz_poly_eval_si(mp_tmp, minus_f_over_g_modp,l);
      minus_f_over_gmodp_l = mpz_fdiv_ui(mp_tmp,(unsigned long)p);     
      
      
      if (m==0) {
	
	// evaluation //gvl = rdict.gmodp(l);
	mpz_poly_eval_si(mp_tmp, rdict.gmodp,l);
	gvl = (long)mpz_fdiv_ui(mp_tmp,(unsigned long)p);	
	
	if (gvl == 0)
	  continue;

	igl = modular_inverse(gvl,p);
	rdict.l0 = l;
	rdict.igl = igl;
	twist_v1 = -rdict.l0;
	minus_f_over_gmodp_l *= -igl;

      }
                 
      
      dv0 = minus_f_over_gmodp_l;
      u1  = u0;
      v1  = v0;
      
      
      hits = light(rdict.sarr, u1,v1, pm, pm1, twist_v1, rdict.umax,rdict.vmax,-scale1);
      if (hits == 0) 
	continue;
      // print "main, level %d : %d hits" % (m, hits)
      thits += hits;
      // We're optimistic here ; we're counting on the fact that we
      // expect the derivative not to cancel.
            
      // However, there is the possibility that we reach the
      // cancellation point. This only once per (p,p) square,
      // meaning that u,v are constrained.

      // Try to guess the right value for du.
      // Rp = Integers(ppow[m-ld+1])
      // RpP=Rp['w']
      // rhs = Z(Z((RpP(fdg_gdf) - u0*RpP(gg)^2)(l))/ppow[m-ld])
      // Slightly faster when not doing reductions while evaluating.
      // Most probably due to the cost of computing remainders all
      // the way to the final value, which ends up being more
      // expensive than it should.

      // TODO: There has to be a way to remove valuations very early
      // on here. In fact, my guess is that it could very probably
      // be done every time m-ld increases...
      
      
      //   This can be simplified.
      //      rhs = Z(K((fdg_gdf(l)-u0*gg(l)^2)/ppow[m-ld]));
      mpz_t rhs,aux;
      mpz_init(rhs);
      mpz_init(aux);
      mpz_poly_eval_si(rhs,fdg_gdf,l);
      mpz_poly_eval_si(aux,gg,l);
      mpz_mul(aux,aux,aux);
      mpz_mul_si(aux,aux,u0);      
      mpz_sub(rhs,rhs,aux);
      mpz_mod_ui(rhs,rhs,(unsigned long)p);
      ASSERT_ALWAYS(mpz_divisible_ui_p(rhs,(unsigned long)rdict.ppow[m-ld]));
      mpz_fdiv_q_ui(rhs,rhs,(unsigned long)rdict.ppow[m-ld]);
      mpz_clear(aux);
      
      
      if ((ld > 0) && (mpz_sgn(rhs)!= 0))
	continue;
      
      long nspots;

      nfdg_gdf = compose_reduce(fdg_gdf,l,rdict,pmax);
      nff = compose_reduce(ff,l,rdict,pmax);
      ngg = compose_reduce(gg,l,rdict,pmax);      
      
      
      if (m == 0) {
	// XXX This is fairly bizarre. Results are correct with
	// this, but I don't see at the moment why we don't have
	// psi instead of x here.
	ndphi = mpz_poly_identity();
      }
      else {
	ndphi = compose_psi(dphi,l,rdict);
      }

      tmp_poly = mpz_poly_sum(nff,mpz_poly_coeff_product_si(ngg,dv0));	
      nff = mpz_poly_coeff_reduction(mpz_poly_coeff_division(tmp_poly,(unsigned long)p),(unsigned long)p);
      twist_f = mpz_poly_prod(ndphi,ngg);		  

      if (ld > 0)
	if (mpz_sgn(rhs) == 0)
	  // Then all u's will do. Note that this case happens
	  // only rather deep inside the recursion. So indeed
	  // we're going to re-hit the places which have been
	  // hit before, but that's not too worrying.
	  nspots=p;
	else {	    
	  
	  du0 = (igl*igl*mpz_get_si(rhs)) % p;
	  u1 += du0 * pm;
	  v1 += du0 * twist_v1;
	  nff = mpz_poly_sum(nff,mpz_poly_coeff_product_si(twist_f,du0));	    
	  nspots=1;	    
	  new_twist_v1 = p * twist_v1;
	    
	  for ( du=0 ; du<nspots; du++) {
	    hits = light(rdict.sarr, u1,v1, pm1, pm1, 0, rdict.umax, rdict.vmax, scale3);
	    // print "main/excep, level %d : %d hits" % (m, hits)
	    thits += hits;
	    // scale3 is scale/(p-1)-scale/p ; hence it brings
	    // back the cells to the contribution -scale/p, which is the
	    // contribution from a _single_ zero at this location, not
	    // liftable because of ramification. Later recursive calls
	    // will investigate the possibility that despite the multiple
	    // root mod p, we still get roots at higher precision.
	    if (hits >0) 
	      thits += rotation_inner(rdict,p,nff,ngg,nfdg_gdf,u1,v1,l0,ld+1,m+1,new_twist_v1,scale2,ndphi);
			
	    nff = mpz_poly_sum(nff,twist_f);
	    u1 += pm;
	    v1 += twist_v1;
	  }
	    
	}
      
    }
  }
  else {
    // f mod p is a constant polynomial.  This means that our
    // expression has an extra root only for specific u,v pairs.
    // Basically, we have a linear term in u and v which must be
    // cancelled.
    
    dv0 = mpz_poly_constant_coeff_si(minus_f_over_g_modp);
    u1  = u0;
    v1  = v0 + pm * dv0;
    hits=light(rdict.sarr, u1, v1, pm, pm1, twist_v1, rdict.umax, rdict.vmax, -scale);
    if (hits == 0)
      return thits;
    // print "secondary, level %d : %d hits" % (m, hits)
    
    thits += hits;
    tmp_poly = mpz_poly_sum(ff,mpz_poly_coeff_product_si(gg,dv0));	
    nff = mpz_poly_coeff_reduction(mpz_poly_coeff_division(tmp_poly,(unsigned long)p),(unsigned long)p);
    twist_f = mpz_poly_prod(dphi,gg);
    new_twist_v1 = p * twist_v1;

    for ( du=0 ; du<p ; du++) {
      hits = rotation_inner(rdict,p,nff,gg,fdg_gdf,u1,v1,l0,ld,m+1,new_twist_v1,scale,dphi);
      // print "secondary/excep, level %d : %d hits" % (m, hits)
      thits += hits;
      u1 += pm;
      v1 += twist_v1;      
      nff = mpz_poly_sum(nff,twist_f);
    }
  }
  return thits;
}
  

int light (double * alpha, long u0, long v0, long um, long vm, long us, long vs, long skew, double contribution) {
  
  // The lattice has its origin in (u0,v0), and
  // it is generated by (us,0),(skew,vs).
  long u,v;
  long pos,dpos,oversize,cq,cr,carriage_return;
  long hits;

  u = u0;
  v = v0;
  pos = u*vm;
  dpos= us*vm;
  v = v % vs;
  hits=0;
  oversize = vm + skew;
  cq = oversize / vs;
  cr = oversize % vs;
  carriage_return=oversize-cr;
  v = v % vs;
  
  while (u < um) {
    // reduce within the row, since we know that b[0] is zero
    // v = v % vs
    while (v < vm) {
      alpha[pos + v] += contribution;
      hits+=1;
      v += vs;
      u += us;
      v += skew;
      pos += dpos;
      v -= carriage_return;
      //assert cr <= v and v < cr + vs;
      if (v >= vs) 
	v -= vs;
    }
  }

  return hits;

  }




 long modular_inverse( long r,  long N ) {

   //ASSERT_ALWAYS(gcd(r,N)==1);
  residue_t res,inv;
  modulus_t modu;

   long result;

  modul_initmod_ul(modu,N);
  modul_init(res,modu);
  modul_init(inv,modu);
  modul_set_ul_reduced(res,r,modu);
  modul_inv(inv,res,modu);
  result = modul_get_ul(inv,modu);
  modul_clear(inv,modu);
  modul_clear(res,modu);
  modul_clearmod(modu);

  return result;
}

void extend_ppow_list(rootsieve_dictionary rdict,int l, long p) {

  rdict.ppow = realloc(rdict.ppow,l*sizeof(long));
  
  if (l<=rdict.maxpow)
    return;

  int i;
  for( i=rdict.maxpow+1 ; i<=l; i++) {
    rdict.ppow[i] = rdict.ppow[i-1]*p;    
  }

  rdict.maxpow = l;

}


/*
 _  _ |   _  _  _ _ . _ | _
|_)(_)|\/| |(_)| | ||(_||_\
|      /                   

 */

 void mpz_poly_eval_si(mpz_t r, mpz_poly * p,  long a) {

    int i;

    mpz_set(r, p->f[p->deg]);

    for (i = p->deg - 1; i >= 0; i--) {
	mpz_mul_si(r, r, a);
	mpz_add(r, r, p->f[i]);
    }
}

void mpz_poly_init(mpz_poly * p, int deg) {

  int i;

  p->deg = deg;
  p->f = malloc((deg+1)*sizeof(mpz_t));

  for( i=0 ; i<=p->deg ; i++) 
     mpz_init(p->f[i]);    
  
}

void mpz_poly_clear(mpz_poly * p) {

  int i;

  for( i=0 ; i<=p->deg ; i++) 
     mpz_clear(p->f[i]);   
  
  free(p->f);
  //mpz_out_str(stdout,10,p.f[0]);
  //printf("\n");
}

void mpz_poly_trim(mpz_poly * p) {
   int i;
   for (i = p->deg; i >= 1; i--)
     if(mpz_sgn(p->f[i])!=0)
       break;

   p->deg = i;   
     
}

 void  mpz_poly_set_linear_si(mpz_poly * p,long l, long p) {
   mpz_set_si(p->f[0],l);
   mpz_set_si(p->f[1],p);
 }

 mpz_poly mpz_poly_identity() {
   mpz_poly id;
   mpz_poly_init(&id,1);
   mpz_set_si(id.f[1],1);
   return id;
 }

 void mpz_poly_set(mpz_poly * p,mpz_t * clist) {

   // We suppose clist has length at least p->deg+1
   int i;

   for( i=0 ; i<=p->deg ; i++)
     mpz_set(p->f[i],clist[i]);
 }

 void mpz_poly_init_set(mpz_poly * f,int deg,mpz_t * clist) {

   // We suppose clist has length at least p->deg+1
   int i;
   
   p->deg = deg;
   p->f = malloc((deg+1)*sizeof(mpz_t));

   for( i=0 ; i<=p->deg ; i++) {
     mpz_init(p->f[i]);
     mpz_set(p->f[i],clist[i]);
   }
 }

// Compute the derivative

mpz_poly mpz_poly_derivative(mpz_poly f) {

  mpz_poly diff;  
   int n; // Loop counter
    
  if (f.deg<=0) {
    mpz_poly_init(diff,0);        
    return diff;
  }
  
  mpz_poly_init(diff,f.deg-1);
  
  for( n=0 ; n<=f.deg-1 ; n++ )
    mpz_mul_si(diff.f[n],f.f[n+1],n+1);
  
  return diff;
  
}

// Reduction mod m of coeffs St: 2h12
mpz_poly mpz_poly_coeff_reduction(mpz_poly f, unsigned long m) {

  mpz_poly reduced;
   int n;
  
  mpz_poly_init(reduced,f.deg);

  for ( n=0 ; n<f.deg+1 ; n++)
    mpz_mod_ui(reduced.f[n],f.f[n],m);

  // Just in case some leading coeffs were multiples of m
  mpz_poly_trim(reduced);

  return reduced;
    
}
// End : 2h21 - 9 minutes

// Losing time : 5 minutes

// Division of every coeff by something : Start 2h27
mpz_poly mpz_poly_coeff_division(mpz_poly f, unsigned long m) {

  mpz_poly divised;
  int n;
  
  mpz_poly_init(divised,f.deg);

  for ( n=0 ; n<f.deg+1 ; n++) {
    ASSERT_ALWAYS(mpz_divisible_ui_p(f.f[n],m));
    mpz_fdiv_q_ui(divised.f[n],f.f[n],m);
  }      

  return divised;
    
} 
// End 2h47 - 20 min.

mpz_poly mpz_poly_coeff_product_si(mpz_poly f,  long m) {

  mpz_poly multiplied;
   int n;
  
  mpz_poly_init(multiplied,f.deg);

  for ( n=0 ; n<f.deg+1 ; n++)     
    mpz_mul_si(multiplied.f[n],f.f[n],m);
  
  return multiplied;
    
} 

// Sum of polynomials 2h47

mpz_poly mpz_poly_sum(mpz_poly f, mpz_poly g) {
  if (f.deg<=g.deg) 
    return mpz_poly_sum_sorted(f,g);
  else 
    return mpz_poly_sum_sorted(g,f);
}

// We assume the degree of g is less than the degree of f.
mpz_poly mpz_poly_sum_sorted(mpz_poly f, mpz_poly g) {
    
  mpz_poly s;
  int n;

  mpz_poly_init(s,g.deg);  
  
  for( n=0 ; n<=f.deg ; n++) {
    mpz_init(s.f[n]);
    mpz_add(s.f[n],f.f[n],g.f[n]);
  }
  
  for( n=f.deg+1 ; n<=g.deg ; n++) {
    mpz_init(s.f[n]);
    mpz_set(s.f[n],g.f[n]);
  }

  //mpz_poly_trim(s);

  return s;

} // End 3h34 - 47 minutes

// Product of polynomials : Start 3h34
mpz_poly mpz_poly_prod(mpz_poly f, mpz_poly g) {
  
  mpz_poly prod;
  mpz_t tmp;
  int i,n;
  
  mpz_poly_init(prod,f.deg+g.deg);

  for( n=0 ; n<=prod.deg ; n++) {
    for ( i=0 ; i<=n ; i++) {
      mpz_mul(tmp,f.f[i],g.f[n-i]);
      mpz_add(prod.f[n],prod.f[n],tmp);
    }
  }

  return prod;
  
} // End 3h53 -- 20 minutes

 mpz_poly compose_psi(mpz_poly f,long l,rootsieve_dictionary rdict) {

   mpz_poly compo;
   long i,k;
   mpz_t tmpsum,tmpterm,tmplpow;
  
   mpz_init(tmpsum);
   mpz_init(tmpterm);
   mpz_init(tmplpow);    

   mpz_poly_init(compo,f.deg);

   for( k=0 ; k<=compo.deg ; k++) {

     mpz_set_si(tmplpow,1);
     mpz_set_si(tmpsum,0);
       
     for( i=k ; i<compo.deg ; i++) {
       // Compute the thing
       mpz_bin_uiui(tmpterm,(unsigned long)i,(unsigned long)k);
       mpz_mul(tmpterm,tmpterm,f.f[i]);
       mpz_mul_ui(tmpterm,tmpterm,(unsigned long)rdict.ppow[k]);
       mpz_mul(tmpterm,tmpterm,tmplpow);
      
       // Update pows
       mpz_mul_si(tmplpow,tmplpow,l);      
      
       // Add to the sum
       mpz_add(tmpsum,tmpsum,tmpterm);      
     }

     mpz_set(compo.f[k],tmpsum);

   }

   return compo;
 }

mpz_poly compose_reduce(mpz_poly f,long l,rootsieve_dictionary rdict, unsigned long m)  {
   
  return mpz_poly_coeff_reduction(compose_psi(f,l,rdict),m);
  
} 

// Constant test
int mpz_poly_is_constant(mpz_poly f) {
  return (f.deg==0);
}

void mpz_poly_constant_coeff(mpz_t op,mpz_poly f) {  
  mpz_set(op,f.f[0]);
}

long mpz_poly_constant_coeff_si(mpz_poly f) {  
  ASSERT_ALWAYS(mpz_fits_slong_p(f.f[0]));
  return mpz_get_si(f.f[0]);
}
