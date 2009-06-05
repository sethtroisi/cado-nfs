#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gmp.h"
#include "gmp_aux.h"
#include "cado.h"
#include "utils.h"
#include "auxiliary.h"
#include "mod_ul_default.h"
#include "mod_ul.h"
#include "fast_rootsieve.h"
#include "macros.h" /* for ASSERT_ALWAYS */
//#define LIGHT_MESSAGES

/* Rootsieve : we compute the root property for many polynomials at once.
   The code is a C port of Emmanuel Thomé's code in Sage.
*/


int main() {
  // This main function just runs an example.
  poly_t f,g,h;
  int degf, degg;
  mpz_t * fco, * gco;   
  int i;
  
  degf = 6;
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
  
  poly_alloc(f,degf);
  poly_alloc(g,degg);

  poly_set(f,fco,degf);
  poly_set(g,gco,degg);

  for( i=0 ; i<= degf; i++)
    mpz_clear(fco[i]);
  
  for( i=0 ; i<= degg; i++)
    mpz_clear(gco[i]);

  free(fco);
  free(gco);

  

  //printf("%d\n",g->deg);
  poly_alloc(h,f->deg);
  rootsieve(h,f,g,300);
  printf("Polynomial f :\n");
  poly_print(f);
  printf("\n");
  printf("Polynomial g :\n");
  poly_print(g);
  printf("\n");
  printf("Rotated polynomial :\n");
  poly_print(h);
  printf("\n");
  

  poly_free(f);
  poly_free(g);
  poly_free(h);

}


/*
int main() {
  int i,j;
  double * alpha;
  long U0,U1,V0,V1;

  U0 = 0;
  U1 =  30;
  V0 =  0;
  V1 =  30;

  alpha = malloc((U1-U0+1)*(V1-V0+1)*sizeof(double));
  
  for ( i=0 ; i<(U1-U0+1)*(V1-V0+1) ; i++ )
    alpha[i] = 0.0;
  
  light_rectangle(alpha,3,5,7,2,3,U0,U1,V0,V1,1.0);

  for ( j=V1-V0 ; j>=0 ; j-- ) {
    for ( i=0 ; i<(U1-U0+1) ; i++ )
      printf("%d ",(int)alpha[i +(U1-U0+1)*j]);  
    printf("\n");
  }

  
  } Light rectangle test */



/* h must be allocated : it must have f->deg+1 slots. */
void rootsieve(poly_t h,poly_t f,poly_t g, unsigned long prime_bound ) {  // Ready
  
  // Counters
  long i,j;
  unsigned long p;

  // The bounds
  long U0,U1,V0,V1;

  // We compute the bounds of the rectangle limiting the rotation coefficients
  rotate_bounds(f->coeff, f->deg, g->coeff[1], g->coeff[0], &V0, &V1, &U0, &U1, 1);
  
  // We prepare the dictionary
  rootsieve_dictionary rdict;

  // The rectangle
  rdict->U0 = U0;
  rdict->U1 = U1;
  rdict->V0 = V0;
  rdict->V1 = V1;
  rdict->space = (U1-U0+1)*(V1-V0+1);  
  rdict->sarr = malloc(rdict->space*sizeof(double));  
  
  for ( i=0 ; i<rdict->space ; i++ )
    rdict->sarr[i] = 0.0;  

  // The rectanglewise upper bound for the powers of prime p.
  rdict->m0 = (int)floor(log(MAX(MAX(MAX(labs(U0),labs(U1)),labs(V0)),labs(V1)))/log(p));

  // Other bounds
  rdict->B = prime_bound;
  rdict->global_offset = 0;
  rdict->cutoff = 1.0e-5;

  for ( p = 2 ; p <= prime_bound ; p += 1 + (1 & p))
    if (isprime(p))    
      rdict->global_offset += log(p)/(p-1);

  // This gives the exceptional values for u (once divided by g(l)^2)
  poly_t dg,df,aux;  
  poly_alloc(dg,g->deg-1);  
  poly_alloc(df,f->deg-1);  
  poly_alloc(aux,-1);  
  poly_alloc(rdict->fdg_gdf,-1);  
  poly_derivative(dg,g);
  poly_derivative(df,f);
  poly_mul(rdict->fdg_gdf,f,dg);
  poly_mul(aux,g,df);  
  poly_sub(rdict->fdg_gdf,rdict->fdg_gdf,aux);
  poly_free(aux);
  poly_free(dg);
  poly_free(df);
  
  poly_alloc(rdict->f,f->deg);
  poly_copy(rdict->f,f);
  poly_alloc(rdict->g,g->deg);
  poly_copy(rdict->g,g);
  
  
  for ( p = 2 ; p <= prime_bound ; p += 1 + (1 & p))
    if (isprime(p))   {
      printf("Filling alpha for prime: %u\n",p);
      fill_alpha(rdict,p);        
    }


  // Finding the best polynomial
  double minimum;
  int u,v;
  int best_u,best_v;

  // Assignment
  minimum = DBL_MAX;
  best_u = 0;
  best_v = 0;

  // Finding the minimum
  for ( j=V1-V0 ; j>=0 ; j-- ) {
    for ( i=0 ; i<(U1-U0+1) ; i++ )
      if(rdict->sarr[i +(U1-U0+1)*j] < minimum) {
	best_u = U0 + i;
	best_v = V0 + j;
	minimum = rdict->sarr[i +(U1-U0+1)*j];
      }      
  }



  
  // Packing the selected rotated polynomial h
  poly_t tmpoly;
  mpz_t tmpmp;
  poly_alloc(tmpoly,1);
  mpz_init_set_si(tmpmp,best_v);
  poly_setcoeff(tmpoly,0,tmpmp);
  mpz_set_si(tmpmp,best_u);
  poly_setcoeff(tmpoly,1,tmpmp);
  
  poly_mul(h,tmpoly,g);
  poly_add(h,h,f);

  printf("Best alpha value : %f\n",minimum);
  printf("The pair u,v : %ld %ld\n",best_u,best_v);
  
  // Cleaning
  mpz_clear(tmpmp);
  poly_free(tmpoly);

  // Dictionary cleaning
  free(rdict->sarr);

  poly_free(rdict->f);
  poly_free(rdict->g);
  poly_free(rdict->fdg_gdf);    
    
}


long fill_alpha(rootsieve_dictionary rdict,unsigned long p) { // Ready
  

  //"""Fast root sieve modulo p and its powers. The rdict argument must
  // have been setup beforehand"""

  // Counters
  int i;
  long hits;

  long u0,v0,l0;
  double scale;
  int ld,m;
  poly_t ff,gg,fdg_gdf;
  
  poly_alloc(ff,rdict->f->deg);
  poly_copy(ff,rdict->f);
  poly_alloc(gg,rdict->g->deg);
  poly_copy(gg,rdict->g);
  poly_alloc(fdg_gdf,rdict->fdg_gdf->deg);  
  scale = (log(p))*(double)p/(double)(p+1);
  
  l0 = 0;
  ld = m = 0;  
  // In theory, one does not need to initialize these values,
  // we do it now to control warnings.
  u0 = 0; 
  v0 = 0; 

  rdict->maxpow = (long) (31*log(2)/log(p));
  rdict->ppow   = malloc((rdict->maxpow+1)*sizeof(mpz_t));

  mpz_init_set_ui(rdict->ppow[0],1);
  for( i=1 ; i<=rdict->maxpow; i++) {
    mpz_init(rdict->ppow[i]);
    mpz_mul_ui(rdict->ppow[i],rdict->ppow[i-1],p);
  }
  
  poly_alloc(rdict->gmodp,rdict->g->deg);
  poly_reduce_mod_mpz(rdict->gmodp,rdict->g,rdict->ppow[1]);

  // You must get the universe right.
  compose_reduce(fdg_gdf,rdict->fdg_gdf,l0,rdict,rdict->ppow[rdict->maxpow]);
  compose_reduce(ff,ff,l0,rdict,rdict->ppow[rdict->maxpow]);
  compose_reduce(gg,gg,l0,rdict,rdict->ppow[rdict->maxpow]);  
  rdict->l0 = l0;  

  poly_t id;
  poly_alloc_identity(id);  
  
  hits = rotation_inner(rdict,p,ff,gg,fdg_gdf,u0,v0,l0,ld,m,0,scale,id); 

  poly_free(ff);
  poly_free(gg);
  poly_free(fdg_gdf);
  poly_free(id);

  for( i=0 ; i<=rdict->maxpow; ++i)
    mpz_clear(rdict->ppow[i]);  
  free(rdict->ppow); 

  poly_free(rdict->gmodp);

  return hits;  
}


long rotation_inner(rootsieve_dictionary rdict, unsigned long p,poly_t ff,poly_t gg,poly_t fdg_gdf, long u0, long v0, long l0,int ld,int m,long twist_v1,double scale,poly_t dphi) {
  
  // Use your precomputations.
  // Update your precomputations

  //if len(hist) > 6:
  //      nhist=copy(hist)
  //      nhist.append([dphi,u0,v0,dphi(0),ld,m])
  //      print "Found tree of height %d" % len(nhist)
  //      for hi in nhist:
  //          print hi
  //

  if (scale < rdict->cutoff) {    
    printf("Cutoff aborted\n");
    return 0;
  }

  // Hit counter
  long thits,hits;
  thits = 0;
  hits = 0;

  // Abbreviations
  long igl;
  long pm,pm1,pmax;
    
  pm    = mpz_get_si(rdict->ppow[m]);
  pm1   = mpz_get_si(rdict->ppow[m+1]);
  pmax  = mpz_get_si(rdict->ppow[rdict->maxpow]);  
 
  
  double scale1,scale2,scale3;
  int look_many_roots;
  poly_t minus_f_over_g_modp;
  
  // gg is never scaled down, and it is always evaluated at points where
  // the value is not zero mod p. So the valuation of g mod p never
  // vanishes.
  
  // This is always useful.  

  scale1 = scale/(p-1);
  scale2 = scale/p;
  scale3 = scale1-scale2;  
  
  //printf("%d",ff->deg);
  //poly_print(ff);
  //printf("\n");
  poly_alloc(minus_f_over_g_modp, ff->deg);
  poly_reduce_mod_mpz(minus_f_over_g_modp, ff,rdict->ppow[1]);

  if (m > 0) {
    l0 = rdict->l0;    
    igl = rdict->igl;
    mpz_t aux;
    mpz_init_set_si(aux,igl);
    poly_mul_mpz(minus_f_over_g_modp,minus_f_over_g_modp,aux);
    mpz_clear(aux);
    poly_reduce_mod_mpz(minus_f_over_g_modp,minus_f_over_g_modp,rdict->ppow[1]);    
  }
  

  // Deciding whether we're in the general case or not has some
  // subtleties, unfortunately.
  look_many_roots = (!poly_is_constant(minus_f_over_g_modp)) ||  (m == 0 &&  ! poly_is_constant(rdict->gmodp));
  
  //nhist = copy(hist);
    
  //if (!look_many_roots)
  //  nhist.append([dphi,u0,v0,l0+p*dphi(0),ld,m]);
  
  long dv0,u1,v1;
  poly_t nfdg_gdf,nff,ngg,ndphi,tmp_poly,twist_f;
  long du,du0,new_twist_v1;

  if (look_many_roots) {
    // Then it's the typical case, as encountered for instance at the
    // beginning of the root sieve.
    // Each possible l value must be tried, and will give rise to
    // potentially intersecting lattices. 
    
    long l;
    poly_alloc(nfdg_gdf,fdg_gdf->deg);
    poly_alloc(nff,ff->deg);
    //printf("Degre gg : %d\n",gg->deg);
    poly_alloc(ngg,gg->deg);
    
    for ( l=0 ; l<p ; l++) {
      
      //nhist=copy(hist);
      //nhist.append([dphi,u0,v0,l0+p*dphi(l),ld,m]);
      long minus_f_over_gmodp_l,gvl,igl;
      mpz_t mp_tmp, mpl;
      mpz_init(mp_tmp);
      mpz_init_set_si(mpl,l);
			
      poly_eval_mod_mpz(mp_tmp, minus_f_over_g_modp,mpl,rdict->ppow[1]);
      poly_free(minus_f_over_g_modp);
      minus_f_over_gmodp_l = mpz_get_si(mp_tmp);     
				      
  
      if (m==0) {

	// evaluation //gvl = rdict->gmodp(l);
	poly_eval_mod_mpz(mp_tmp,rdict->gmodp,mpl,rdict->ppow[1]);
	gvl = mpz_get_si(mp_tmp);     

	if (gvl == 0)
	  continue;
					
	igl = (long)modular_inverse((unsigned long)gvl,p);
	rdict->l0 = l;
	rdict->igl = igl;
	twist_v1 = -rdict->l0;
	// The next line divides f(l) by g(l) mod p
	minus_f_over_gmodp_l = (long)modular_product((unsigned long)((minus_f_over_gmodp_l==0)?(minus_f_over_gmodp_l):(p-minus_f_over_gmodp_l)),(unsigned long)igl,p);

      }

      mpz_clear(mp_tmp);

      
      dv0 = minus_f_over_gmodp_l;
      u1  = u0;
      v1  = v0 + pm*dv0;
#ifdef LIGHT_MESSAGES
      printf("turn on the light 1\n");
#endif
      hits = light_rectangle(rdict->sarr, u1,v1, pm, pm1, twist_v1, rdict->U0,rdict->U1,rdict->V0,rdict->V1,-scale1);     
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
      poly_eval(rhs,fdg_gdf,mpl);
      poly_eval(aux,gg,mpl);
      mpz_mul(aux,aux,aux);
      mpz_mul_si(aux,aux,u0);      
      mpz_sub(rhs,rhs,aux);
      mpz_fdiv_q(rhs,rhs,rdict->ppow[m-ld]);
      mpz_clear(aux);
      mpz_clear(mpl);
      
      
      if ((ld > 0) && (mpz_sgn(rhs)!= 0))
	continue;
      
      long nspots;      

      compose_reduce(nfdg_gdf,fdg_gdf,l,rdict,rdict->ppow[rdict->maxpow]);
      compose_reduce(nff,ff,l,rdict,rdict->ppow[rdict->maxpow]);
      compose_reduce(ngg,gg,l,rdict,rdict->ppow[rdict->maxpow]);  


      if (m == 0) {
	// XXX This is fairly bizarre. Results are correct with
	// this, but I don't see at the moment why we don't have
	// psi instead of x here.
	poly_alloc_identity(ndphi);
      }
      else {
	poly_alloc(ndphi,dphi->deg);
	compose_psi(ndphi,dphi,l,rdict);
      }

      poly_alloc(tmp_poly,ngg->deg);
      poly_mul_ui(tmp_poly,ngg,(unsigned long)dv0);		
      poly_add(tmp_poly,nff,tmp_poly);
      poly_div_ui_mod_ui(nff,tmp_poly,p,p);
      poly_alloc(twist_f,-1);
      poly_mul(twist_f,ndphi,ngg);				

      if (ld > 0) {
	if (mpz_sgn(rhs) == 0)
	  // Then all u's will do. Note that this case happens
	  // only rather deep inside the recursion. So indeed
	  // we're going to re-hit the places which have been
	  // hit before, but that's not too worrying.
	  nspots=p;
	
      }
      else {	
	du0 = (igl*igl*mpz_get_si(rhs)) % p;
	u1 += du0 * pm;
	v1 += du0 * twist_v1;
	poly_mul_ui(tmp_poly,twist_f,(unsigned long)du0);
	poly_add(nff,nff,tmp_poly);
	nspots=1;
      }	    
		
      mpz_clear(rhs);      
      new_twist_v1 = p * twist_v1;

      for ( du=0 ; du<nspots; du++) {
#ifdef LIGHT_MESSAGES
	printf("turn on the light 2\n");
#endif
	hits = light_rectangle(rdict->sarr, u1,v1, pm1, pm1, 0, rdict->U0, rdict->U1, rdict->V0, rdict->V1, scale3);
	
	// print "main/excep, level %d : %d hits" % (m, hits)
	thits += hits;
	// scale3 is scale/(p-1)-scale/p ; hence it brings
	// back the cells to the contribution -scale/p, which is the
	// contribution from a _single_ zero at this location, not
	// liftable because of ramification. Later recursive calls
	// will investigate the possibility that despite the multiple
	// root mod p, we still get roots at higher precision.
	if (hits >0) {	      
	  
	  if (m<rdict->m0 && m<rdict->maxpow) 
	    thits += rotation_inner(rdict,p,nff,ngg,nfdg_gdf,u1,v1,l0,ld+1,m+1,new_twist_v1,scale2,ndphi);	      	    
	  else {
	    if(m == rdict->m0)
	      printf("Rectangle aborted ");
	    
	    if(m == rdict->maxpow)
	      printf("Power aborted");

	    printf("\n");
	  }	    
	}

	poly_add(nff,nff,twist_f);		
	u1 += pm;
	v1 += twist_v1;
      }	
    }
    poly_free(nfdg_gdf);
    poly_free(nff);
    poly_free(ngg);
    poly_free(ndphi);
    poly_free(twist_f);
  }
  else {
    // f mod p is a constant polynomial.  This means that our
    // expression has an extra root only for specific u,v pairs.
    // Basically, we have a linear term in u and v which must be
    // cancelled.
    mpz_t mp_tmp;
    mpz_init(mp_tmp);
    poly_getcoeff(mp_tmp,0,minus_f_over_g_modp);
    poly_free(minus_f_over_g_modp);
    dv0 = mpz_get_si(mp_tmp);
    mpz_clear(mp_tmp);
    u1  = u0;
    v1  = v0 + pm * dv0;
#ifdef LIGHT_MESSAGES
    printf("turn on the light 3\n");
#endif
    hits=light_rectangle(rdict->sarr, u1, v1, pm, pm1, twist_v1, rdict->U0, rdict->U1, rdict->V0, rdict->V1, -scale);
    if (hits == 0)
      return thits;
    // print "secondary, level %d : %d hits" % (m, hits)

    thits += hits;

    poly_alloc(nff,ff->deg);
    poly_alloc(tmp_poly,gg->deg);
    poly_mul_ui(tmp_poly,gg,(unsigned long)dv0);
    poly_add(tmp_poly,ff,tmp_poly);
    poly_div_ui_mod_ui(nff,tmp_poly,p,p);
    poly_free(tmp_poly);
    poly_mul(twist_f,dphi,gg);				
    new_twist_v1 = p * twist_v1;
    
    for ( du=0 ; du<p ; du++) {
      if (m<rdict->m0 && m<rdict->maxpow) 
	hits = rotation_inner(rdict,p,nff,gg,fdg_gdf,u1,v1,l0,ld,m+1,new_twist_v1,scale,dphi);
      else {
	if(m == rdict->m0)
	  printf("Rectangle aborted ");
	
	if(m == rdict->maxpow)
	  printf("Power aborted");
	
	printf("\n");
      }	    
      // print "secondary/excep, level %d : %d hits" % (m, hits)
      thits += hits;
      u1 += pm;
      v1 += twist_v1;      
      poly_add(nff,nff,twist_f);
    }
    
    poly_free(nff);
  }

  return thits;
}


long light_rectangle (double * alpha, long u0, long v0, long us, long vs, long skew, long U0, long U1, long V0, long V1, double contribution) {
  
  // The lattice has its origin in (u0,v0), and it is generated by (us,0),(skew,vs).
  // We sieve the u,v rectangle of lower left coordinates (U0,V0) and upper right 
  // coordinates (U1,V1). We start from the lower left, and we sieve first from left to
  // right, and then from the lower to the upper end.
  long u,v;
  long q;
  long pos;
  long hits;

  hits=0;

  // We will start at the lower left lattice point
  
  q = (long)floor((v0-V0)/vs);
  v = v0-q*vs;
  u = u0-q*skew;
  q = (long)floor((u-U0)/us);
  u = u - q*us;
  pos = (U1-U0+1)*(v-V0) + (u-U0);
  
  // at this stage, (u,v) is the lower left point of the lattice
  
  while(v<=V1) {
    
    while(u<=U1) {      
      alpha[pos] += contribution;
      u += us;
      pos += us;   
      hits++;
    }
    
    // we do the carriage return
    u += skew;
    v += vs;
    q = (long)floor((u-U0)/us);
    u = u - q*us;
    pos = (U1-U0+1)*(v-V0) + (u-U0);
    
    // now, we are in the left end of the rectangle.
    
  }

  return hits;

}

unsigned long modular_product( unsigned long a, unsigned long b, unsigned long N ) {

  residue_t res,opa,opb;
  unsigned long result;
  modulus_t modu;
	
  modul_initmod_ul(modu,N);
  modul_init(res,modu);
  modul_init(opa,modu);
  modul_init(opb,modu);
  modul_set_ul_reduced(opa,a,modu);
  modul_set_ul_reduced(opb,b,modu);	
  modul_mul(res,opa,opb,modu);
  result = modul_get_ul(res,modu);
  modul_clear(opa,modu);
  modul_clear(opb,modu);
  modul_clear(res,modu);
  modul_clearmod(modu);

  return result;
}


unsigned long modular_inverse( unsigned long r,  unsigned long N ) {

  //ASSERT_ALWAYS(gcd(r,N)==1);
  residue_t res,inv;
  modulus_t modu;

  unsigned long result;

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

// Polynomials
// Ensures that f(X)=X
void poly_alloc_identity(poly_t f) {
  poly_alloc(f,1);
  mpz_set_ui(f->coeff[0],0);
  mpz_set_ui(f->coeff[1],1);
}

// f <- g(psi) 
// f must be allocated (with the degree of g for example), f and g can be the same.
void compose_psi(poly_t f, const poly_t g, unsigned long l, rootsieve_dictionary rdict) {
  
  int d;  
  d = g->deg;

  long i,k;
  mpz_t tmpsum,tmpterm,tmplpow;
  poly_t aux;
  poly_alloc(aux,d);
  
  mpz_init(tmpsum);
  mpz_init(tmpterm);
  mpz_init(tmplpow);    

  for( k=0 ; k<=d ; k++) {

    mpz_set_ui(tmplpow,1);
    mpz_set_ui(tmpsum,0);
       
    for( i=k ; i<d ; i++) {
      // Compute the thing
      mpz_bin_uiui(tmpterm,(unsigned long)i,(unsigned long)k);
      mpz_mul(tmpterm,tmpterm,g->coeff[i]);
      mpz_mul_ui(tmpterm,tmpterm,(unsigned long)rdict->ppow[k]);
      mpz_mul(tmpterm,tmpterm,tmplpow);
      
      // Update pows
      mpz_mul_ui(tmplpow,tmplpow,l);      
      
      // Add to the sum
      mpz_add(tmpsum,tmpsum,tmpterm);      
    }

    poly_setcoeff(aux,i,tmpsum);
  }

  poly_copy(f,aux);

  mpz_clear(tmpsum);
  mpz_clear(tmpterm);
  mpz_clear(tmplpow);  
  poly_free(aux);

}
/* We assume (rdict->delta) > 0. No problem when f and g are the same. */
void compose_reduce(poly_t f, const poly_t g,unsigned long l,rootsieve_dictionary rdict, mpz_t m) {
  poly_t aux;
  poly_alloc(aux,g->deg+1);
  compose_psi(aux,g,l,rdict);
  poly_reduce_mod_mpz(f,aux,m);
  poly_free(aux);
}



