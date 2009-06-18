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
#include "macros.h"
#define VERBOSE 1
#define DEBUG 0
#define PLB 2
#define PUB 2000

/* Rootsieve : we compute the root property for many polynomials at once.
   The code is a C port of Emmanuel Thomé's code in Sage.
*/


int main() {
  // This main function just runs an example.
  poly_t f,g,h;
  int degf, degg;  
  
  degf = 6;
  degg = 1;

  poly_alloc(f,degf);
  poly_alloc(g,degg);

  poly_setcoeff_str(f,6,"265482057982680",10);
  poly_setcoeff_str(f,5,"1276509360768321888",10);
  poly_setcoeff_str(f,4,"-5006815697800138351796828",10);
  poly_setcoeff_str(f,3,"-46477854471727854271772677450",10);
  poly_setcoeff_str(f,2,"6525437261935989397109667371894785",10);
  poly_setcoeff_str(f,1,"-18185779352088594356726018862434803054",10);
  poly_setcoeff_str(f,0,"-277565266791543881995216199713801103343120",10);

  poly_setcoeff_str(g,1,"34661003550492501851445829",10);
  poly_setcoeff_str(g,0,"-1291187456580021223163547791574810881",10);

  poly_alloc(h,f->deg);
  rootsieve(h,f,g,PUB);
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

/* h must be allocated : it must have f->deg+1 slots. */
void rootsieve(poly_t h,poly_t f,poly_t g, unsigned long prime_bound ) {  // Ready
  
  if(DEBUG) {
      printf("Received polynomials f ...\n");
      poly_print(f);
      printf("and g ... \n");
      poly_print(g);
  }

  
  // Counters
  long i,j;
  unsigned long p;
  long hits;

  // The bounds
  long U0,U1,V0,V1;

  // We compute the bounds of the rectangle limiting the rotation coefficients
  if(DEBUG) 
    printf("Computing rotation bounds...\n");
  rotate_bounds(f->coeff, f->deg, g->coeff[1], g->coeff[0], &V0, &V1, &U0, &U1, 0);
  
  // We prepare the dictionary
  rootsieve_dictionary rdict;
  
  // The rectangle
  rdict->U0 = U0;
  rdict->U1 = U1;
  rdict->V0 = V0;
  rdict->V1 = V1;
  rdict->space = (U1-U0+1)*(V1-V0+1);  
  rdict->sarr = malloc(rdict->space*sizeof(double));  

  if(DEBUG) {
    printf("Rotation takes place in the rectangle of lower left corner (%ld,%ld) and upper right corner (%ld,%ld) ...\n",U0,V0,U1,V1);
    printf("Array alpha has length (U1-U0+1)*(V1-V0+1) = %ld ...\n",rdict->space);
  }
  
  for ( i=0 ; i<rdict->space ; i++ )
    rdict->sarr[i] = 0.0;  
  
  if(DEBUG) {
    printf("Alpha array initialized to 0.0 everywhere...\n");
  }

  // Other bounds
  rdict->B = prime_bound;
  rdict->cutoff = 1.0e-5;

  if(DEBUG) {
    printf("Prime upper bound : %d\n",rdict->B);
  }


  // This gives the exceptional values for u (once divided by g(l)^2)
  poly_t dg,df,aux;  
  poly_alloc(dg,g->deg-1);  
  poly_alloc(df,f->deg-1);  
  poly_alloc(aux,-1);  
  poly_alloc(rdict->fdg_gdf,-1);  
  poly_derivative(dg,g);
  poly_derivative(df,f);

  if(DEBUG) {
    printf("Computed derivative of g ...\n");
    poly_print(dg);
    printf("Computed derivative of f ...\n");
    poly_print(df);    
  }

  poly_mul(rdict->fdg_gdf,f,dg);
  poly_mul(aux,g,df);  
  poly_sub(rdict->fdg_gdf,rdict->fdg_gdf,aux);

  if(DEBUG) {
      printf("Computed f*dg - g*df ...\n");
      poly_print(rdict->fdg_gdf);
  }

  poly_free(aux);
  poly_free(dg);
  poly_free(df);
  
  poly_alloc(rdict->f,f->deg);
  poly_copy(rdict->f,f);
  poly_alloc(rdict->g,g->deg);
  poly_copy(rdict->g,g);

  if(DEBUG) {
    printf("Added f and g to the dictionary...\n");
    printf("Now entering prime number loop...\n");
  }

  hits = 0;
  
  for ( p = PLB ; p <= prime_bound ; p += 1 + (1 & p))
    if (isprime(p))   {
      if(DEBUG) {
	printf("Filling alpha for prime: %u\n",p);
      }
      hits += fill_alpha(rdict,p);        
    }

  if(DEBUG) {
    printf("Finished prime number loop...\n");
    printf("Counted %ld hits.\n",hits);
    printf("--------------------------------------------------\n");
    //printf("Resulting alpha array :\n");
  }
  
  //for ( j=V1-V0 ; j>=0 ; j-- ) {
  //  for ( i=0 ; i<(U1-U0+1) ; i++ )
  //    printf("%.3f",rdict->sarr[i +(U1-U0+1)*j]);  
  //  printf("\n");
  //}
  
  if(DEBUG) {
    printf("Finding a best polynomial ...\n");
  }

  // Finding the best polynomial
  double minimum;
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

  if(VERBOSE || DEBUG) {
    printf("Best alpha value : %f\n",minimum);
    printf("The pair u,v : %d %d\n",best_u,best_v);
  }
  
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
  
  if(DEBUG) {    
    printf("Creating polynomials ff and gg ...\n");
  }
  poly_alloc(ff,rdict->f->deg);
  poly_copy(ff,rdict->f);
  
  if(DEBUG) {    
    printf("Polynomial ff ...\n");
    poly_print(ff);
  }
  poly_alloc(gg,rdict->g->deg);
  poly_copy(gg,rdict->g);

  if(DEBUG) {    
    printf("Polynomial gg ...\n");
    poly_print(gg);
  }  

  poly_alloc(fdg_gdf,rdict->fdg_gdf->deg);  
  scale = log(p)*(double)p/(double)(p+1);
  
  // The rectanglewise upper bound for the powers of prime p.  
  rdict->m0 = (int)floor(log(MAX(labs(rdict->U1-rdict->U0),labs(rdict->V1-rdict->V0)))/log(p));

  if(DEBUG) {
    printf("The rectangle upper bound is %d ...\n",rdict->m0);
  }

  l0 = 0;
  ld = m = 0;  
  // In theory, one does not need to initialize these values,
  // we do it now to control warnings.
  u0 = 0; 
  v0 = 0; 
  
  
  // We must be able to compose polynomials with psi(x) = l+p*x, that's the reason of taking the max of the degrees.
  rdict->maxpow = MAX(MAX(rdict->m0 +1,rdict->g->deg),MAX(rdict->f->deg,rdict->fdg_gdf->deg)); 

  ASSERT_ALWAYS(rdict->maxpow >= rdict->m0);

  if(DEBUG) {    
    printf("Maximum power of %u is %d ...\n",p,rdict->maxpow);
  }

  rdict->ppow = malloc((rdict->maxpow+1)*sizeof(mpz_t));

  mpz_init_set_ui(rdict->ppow[0],1);
  for( i=1 ; i<=rdict->maxpow; i++) {
    mpz_init(rdict->ppow[i]);
    mpz_mul_ui(rdict->ppow[i],rdict->ppow[i-1],p);
  }

  if(DEBUG) {    
    printf("Initialized powers of %u ...\n",p);
  }
  
  poly_alloc(rdict->gmodp,rdict->g->deg);
  poly_reduce_mod_mpz(rdict->gmodp,rdict->g,rdict->ppow[1]);
  
  if(DEBUG) {    
    printf("Reduced polynomial g modulo %u, the result is\n",p);
     poly_print(rdict->gmodp);
  }

  if(DEBUG) {    
     printf("Composing polynomials ...\n");
  }
  // You must get the universe right.
  compose_reduce(fdg_gdf,rdict->fdg_gdf,l0,rdict->ppow,rdict->maxpow,rdict->ppow[rdict->maxpow]);

  if(DEBUG) {        
     printf("Composition/reduction of fdg_gdf, result is\n");
     poly_print(fdg_gdf);
     printf("and the original polynomial fdg_gdf is\n");
     poly_print(rdict->fdg_gdf);
  }

  compose_reduce(ff,ff,l0,rdict->ppow,rdict->maxpow,rdict->ppow[rdict->maxpow]);

  if(DEBUG) {    
     printf("Composition/reduction of ff, result is\n");
     poly_print(ff);
  }

  compose_reduce(gg,gg,l0,rdict->ppow,rdict->maxpow,rdict->ppow[rdict->maxpow]);  
  rdict->l0 = l0;  

  if(DEBUG) {    
     printf("Composition/reduction of gg, result is\n");
     poly_print(gg);
  }

  poly_t id;
  poly_alloc(id,1);
  poly_set_identity(id);  
  
  if(DEBUG) {    
     printf("Entering rootsieving recursion ...\n");     
  }

  hits = rotation_inner(rdict,p,ff,gg,fdg_gdf,u0,v0,l0,ld,m,0,scale,id); 

  if(DEBUG) {    
    printf("Finished rootsieving modulo %u and its powers...\n",p);    
  }

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
  
  // Some necessary assertions
  ASSERT_ALWAYS(m>=ld);

  if(DEBUG) {
    printf("(Rotation Inner Start) m=%d lambda=%d\n",m,ld);
  }
  
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
 
  if(DEBUG) {
    printf("Charged powers of %u : %ld, %ld, %ld\n",p,pm,pm1,pmax);
  }
  
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
  
  poly_alloc(minus_f_over_g_modp, ff->deg);
  poly_reduce_mod_mpz(minus_f_over_g_modp, ff,rdict->ppow[1]);  
  
  if(DEBUG) {
    printf("Polynomial -f/g mod %lu :\n",p);
    poly_print(minus_f_over_g_modp);
  }

  if (m > 0) {
    if(DEBUG) {
      printf("Entering positive m section...\n");
    }
    l0 = rdict->l0;    
    igl = rdict->igl;
    if(DEBUG) {
      printf("l0=%ld, 1/g(l) mod p=%ld\n",l0,igl);
      
    }    
    ASSERT_ALWAYS(0<igl && igl<p);
    mpz_t aux;
    mpz_init_set_si(aux,igl);
    poly_mul_mpz(minus_f_over_g_modp,minus_f_over_g_modp,aux);
    mpz_clear(aux);
    poly_reduce_mod_mpz(minus_f_over_g_modp,minus_f_over_g_modp,rdict->ppow[1]);    
    if(DEBUG) {
      printf("The definitive -f/g mod %lu : \n",p);
      poly_print(minus_f_over_g_modp);
    }
  }
  

  // Deciding whether we're in the general case or not has some
  // subtleties, unfortunately.
  if(DEBUG) {
    printf("Computing the \"Possibly multiple root\" condition.  \n");
    printf("The polynomial -f/g mod %lu is%s constant\n",p,(poly_is_constant(minus_f_over_g_modp))?"":" not");
    if(m==0) {
      printf("The polynomial g mod %lu\n",p);
      poly_print(rdict->gmodp);
      printf("is%s constant\n",(poly_is_constant(rdict->gmodp))?"":" not");
    }
  }
  
  // The following must hold, otherwise g is not irreducible
  //ASSERT_ALWAYS(rdict->gmodp->deg >=0);
  if(rdict->gmodp->deg == -1) {

    if(DEBUG) {
      printf("Error : g can be factored by %lu, hence it is not irreducible. Halting.\n",p);
    }

    ASSERT_ALWAYS(thits==0);
    return thits;
  }

  look_many_roots = (!poly_is_constant(minus_f_over_g_modp)) ||  (m == 0 &&  ! poly_is_constant(rdict->gmodp));
  


  if(DEBUG) {
    printf("We%s expect multiple roots.\n",(look_many_roots)?"":" don't");
  }
  //nhist = copy(hist);
    
  //if (!look_many_roots)
  //  nhist.append([dphi,u0,v0,l0+p*dphi(0),ld,m]);
  
  long dv0,u1,v1;
  long du,du0,new_twist_v1;
  poly_t nfdg_gdf,nff,ngg,ndphi,twist_f,tmp_poly;
  mpz_t rhs;
  
  poly_alloc(nfdg_gdf,fdg_gdf->deg);
  poly_alloc(nff,ff->deg);
  poly_alloc(ngg,gg->deg);
  poly_alloc(ndphi,1);
  poly_alloc(twist_f,0);
  poly_alloc(tmp_poly,0);
  mpz_init(rhs);

  if (look_many_roots) {
    // Then it's the typical case, as encountered for instance at the
    // beginning of the root sieve.
    // Each possible l value must be tried, and will give rise to
    // potentially intersecting lattices. 
    
    if(DEBUG) {
      printf("Entering the multiple root section...\n");
    }
    
    long l;

    
    for ( l=0 ; l<p ; l++) {
      
      //nhist=copy(hist);
      //nhist.append([dphi,u0,v0,l0+p*dphi(l),ld,m]);
      long minus_f_over_gmodp_l,gvl;
      mpz_t mp_tmp, mpl;
      mpz_init(mp_tmp);
      mpz_init_set_si(mpl,l);
      
      poly_eval_mod_mpz(mp_tmp, minus_f_over_g_modp,mpl,rdict->ppow[1]);
      minus_f_over_gmodp_l = mpz_get_si(mp_tmp);     

      if(DEBUG) {
	printf("Evaluated -f(%ld)/g(%ld) mod %lu to ",l,l,p);
	mpz_out_str(stdout,10,mp_tmp);
	printf("\n");
      }
  
      if (m==0) {

	// evaluation 
	//gvl = rdict->gmodp(l);
	
	if(DEBUG) {
	  printf("Entered m==0 section ...\n");
	}
	
	poly_eval_mod_mpz(mp_tmp,rdict->gmodp,mpl,rdict->ppow[1]);
	gvl = mpz_get_si(mp_tmp);

	if(DEBUG) {
	  printf("Evaluated g(%ld) mod %lu to %ld ...\n",l,p,gvl);
	}

	if (gvl == 0) {
	  if(DEBUG) {
	    printf("Singular g(%ld) mod %lu aborted...\n",l,p);
	  }
	  mpz_clear(mp_tmp);
	  mpz_clear(mpl);	  
	  continue;
	}
		
	igl = (long)modular_inverse((unsigned long)gvl,p);
	
	rdict->l0 = l;
	rdict->igl = igl;
	twist_v1 = -rdict->l0;

	if(DEBUG) {
	  printf(" 1/g(%ld) mod %lu = %ld, twist_v1 = %ld ...\n",l,p,igl,twist_v1);	  
	}
	
	ASSERT_ALWAYS(0<igl && igl<p);

	// The next line divides f(l) by g(l) mod p
	minus_f_over_gmodp_l = (long)modular_product((unsigned long)((minus_f_over_gmodp_l==0)?(minus_f_over_gmodp_l):(p-minus_f_over_gmodp_l)),(unsigned long)igl,p);
	
	if(DEBUG) {
	  printf("-f(%ld)/g(%ld) mod %lu = %ld ...\n",l,l,p,minus_f_over_gmodp_l);
	}

	ASSERT_ALWAYS(0<=minus_f_over_gmodp_l && minus_f_over_gmodp_l<p);
      }

      mpz_clear(mp_tmp);
      
      dv0 = minus_f_over_gmodp_l;
      u1  = u0;
      v1  = v0 + pm*dv0;

      if(DEBUG) {
	printf("The lattice : dv0 = %ld, u1 = %ld, v1 = %ld...\n",dv0,u1,v1);
      }

      hits = light_rectangle(rdict->sarr, u1,v1, pm, pm1, twist_v1, rdict->U0,rdict->U1,rdict->V0,rdict->V1,-scale1);

      if (hits == 0) {
	if(DEBUG) {
	  printf("No hits, continuing to next iteration ...\n");
	}		
	mpz_clear(mpl);	
	continue;	
      }

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
     
      mpz_t aux;      
      mpz_init(aux);
      poly_eval(rhs,fdg_gdf,mpl);
      poly_eval(aux,gg,mpl);
      mpz_mul(aux,aux,aux);
      mpz_mul_si(aux,aux,u0);      
      mpz_sub(rhs,rhs,aux);

      ASSERT_ALWAYS(m>=ld);

      mpz_fdiv_q(rhs,rhs,rdict->ppow[m-ld]);
      mpz_clear(aux);
      mpz_clear(mpl);     
      
      
      if ((ld > 0) && (mpz_sgn(rhs)!= 0)) {

	if(DEBUG) {
	  printf("Root is simple, aborting the anticontribution section...\n");
	}
	continue;
	
      }
	
      long nspots;
      nspots = -1; // This bogus initial value is used for testing that assignment is actually done.

      compose_reduce(nfdg_gdf,fdg_gdf,l,rdict->ppow,rdict->maxpow,rdict->ppow[rdict->maxpow]);      
      compose_reduce(nff,ff,l,rdict->ppow,rdict->maxpow,rdict->ppow[rdict->maxpow]);
      compose_reduce(ngg,gg,l,rdict->ppow,rdict->maxpow,rdict->ppow[rdict->maxpow]);  

      if(DEBUG) {
	  printf("Composed reduced polynomials, nfdg_fdg is...\n");
	  poly_print(nfdg_gdf);
	  printf("nff is...\n");
	  poly_print(nff);
	  printf("ngg is...\n");
	  poly_print(ngg);
      }

      
      if (m == 0) {
	// XXX This is fairly bizarre. Results are correct with
	// this, but I don't see at the moment why we don't have
	// psi instead of x here.

	if(DEBUG) {
	  printf("Fairly bizarre (m==0)...\n");
	}
	
	poly_set_identity(ndphi);

	if(DEBUG) {
	  printf("Polynomial ndphi is ...\n");
	  poly_print(ndphi);
	}		
      }
      else {        
	compose_psi(ndphi,dphi,l,rdict->ppow,rdict->maxpow);
	if(DEBUG) {
	  printf("Polynomial ndphi is ...\n");
	  poly_print(ndphi);
	}	
      }
      
      poly_mul_ui(tmp_poly,ngg,(unsigned long)dv0);       
      poly_add(tmp_poly,nff,tmp_poly);      
      poly_div_ui_mod_ui(nff,tmp_poly,p,p);            
      poly_mul(twist_f,ndphi,ngg);				

      if(DEBUG) {
	  printf("Polynomial twist_f is ...\n");
	  poly_print(twist_f);
      }

      if (ld > 0) {
	if(DEBUG) {
	  printf("Entered lambda>0 spots section...\n");
	}
	// Then all u's will do. Note that this case happens
	// only rather deep inside the recursion. So indeed
	// we're going to re-hit the places which have been
	// hit before, but that's not too worrying.

	// This test was redundant, left here for reference.
	ASSERT_ALWAYS(mpz_sgn(rhs) == 0);
	nspots=p;		
      }
      else {	
	if(DEBUG) {
	  printf("Entered lambda=0 spots section...\n");
	  printf("Value of g(l) mod p is %ld ...\n",igl);	
	}
      	
	ASSERT_ALWAYS(igl>0 && igl<p);

	du0 = (igl*igl*mpz_get_si(rhs)) % p;
	u1 += du0 * pm;
	v1 += du0 * twist_v1;
	poly_mul_ui(tmp_poly,twist_f,(unsigned long)du0);
	poly_add(nff,nff,tmp_poly);

	if(DEBUG) {
	  printf("Constructed new f(l+p*x) polynomial :\n");
	  poly_print(nff);
	}
	nspots=1;
      }	          

      if(DEBUG) {
	printf("Number of spots to visit : %ld\n",nspots);	
      }
      
      ASSERT_ALWAYS(nspots>0);      
      new_twist_v1 = p * twist_v1;

      if(DEBUG) {
	printf("The new twist is new_twist_v1=%ld \n",new_twist_v1);
      }

      if(DEBUG) {
	printf("Entering anticontribution loop...\n");
      }
      
      
      for ( du=0 ; du<nspots; du++) {

	if(DEBUG) {
	  printf("Distributing contribution\n");
	}
	hits = light_rectangle(rdict->sarr, u1,v1, pm1, pm1, 0, rdict->U0, rdict->U1, rdict->V0, rdict->V1, scale3);
	
	if(DEBUG) {
	  printf("Had %ld hits...\n",hits);
	}

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
	    if(DEBUG) {
	      if(m == rdict->m0)
		printf("Rectangle aborted ");
	      
	      if(m == rdict->maxpow)
		printf("Power aborted");
	      
	      printf("\n");
	    }	    	    
	  }	    
	}
	else {
	  if(DEBUG) {
	    printf("No recursion because no hits...\n");
	  }	  
	}

	poly_add(nff,nff,twist_f);		
	u1 += pm;
	v1 += twist_v1;
	
	if(DEBUG) {
	  printf("Updated lattice...\n");
	  printf("u1=%ld, v1=%ld, and f(l+p*x) is\n",u1,v1);
	  poly_print(nff);
	}
	}
	
      }
    
    if(DEBUG) {
      printf("Completed multiple root section: exiting...\n");
    }    
  }
  else {
    // f mod p is a constant polynomial.  This means that our
    // expression has an extra root only for specific u,v pairs.
    // Basically, we have a linear term in u and v which must be
    // cancelled.

    if(DEBUG) {
      printf("Entering simple root section...\n");
      poly_t auxpoly;
      poly_alloc(auxpoly,0);
      poly_reduce_mod_mpz(auxpoly,ff,rdict->ppow[1]);
      ASSERT_ALWAYS(poly_is_constant(auxpoly));
      //printf("f mod %lu is supposed to be constant, and the test gives %s \n",p,(poly_is_constant(auxpoly))?"constant":"not constant");
      poly_free(auxpoly);
    }
 
    mpz_t mp_tmp;
    mpz_init(mp_tmp);

    if (m==0) {
      
      // This case is rare in practice, but it appeared in the example in the
      // example given in the file kleinjung.sage, for p=37. 
      //
      // We need to initialize 1/g(l) mod p, knowing that g is constant.     
      // IMPORTANT : we know that this inverse igl will be valid for all values
      // of l, that's why we can compute it here, instead of in the main loop.

      long gvl,igl;             
	
      if(DEBUG) {
	printf("Entered m==0 section ...\n");
      }

      poly_getcoeff(mp_tmp,0,rdict->gmodp);      
      gvl = mpz_get_si(mp_tmp);
      ASSERT_ALWAYS(0<=gvl && gvl<p);

      if(DEBUG) {
	printf("Constant polynomial g mod %lu evaluates to %ld ...\n",p,gvl);
      }
      
      ASSERT_ALWAYS(gvl!=0); // Otherwise g is reducible
		
      igl = (long)modular_inverse((unsigned long)gvl,p);
	
      rdict->igl = igl;

      if(DEBUG) {
	printf(" 1/g mod %lu = %ld, ...\n",p,igl);	  
      }
	
      ASSERT_ALWAYS(0<igl && igl<p);
      
    }
        
    poly_getcoeff(mp_tmp,0,minus_f_over_g_modp);

    if(DEBUG) {
      printf("The constant function -f/g mod %lu equals ",p);
      mpz_out_str(stdout,10,mp_tmp);
      printf("\n");
    }

    dv0 = mpz_get_si(mp_tmp);
    mpz_clear(mp_tmp);
    u1  = u0;
    v1  = v0 + pm * dv0;
    
    if(DEBUG) {
      printf("New lattice : u1=%ld, v1=%ld\n Adding contribution...",u1,v1);
    }
    
    hits=light_rectangle(rdict->sarr, u1, v1, pm, pm1, twist_v1, rdict->U0, rdict->U1, rdict->V0, rdict->V1, -scale);
    
    if(DEBUG) {
      printf("Counted %ld hits in the rectangle.\n",hits);
    }

    if (hits == 0) {
      
      if(DEBUG) {
	printf("No hits, returning total hits and terminating (Rotation Inner)...\n");
      }
                 
      poly_free(minus_f_over_g_modp);
      poly_free(nfdg_gdf);
      poly_free(nff);
      poly_free(ngg);
      poly_free(ndphi);
      poly_free(twist_f);
      poly_free(tmp_poly);
      mpz_clear(rhs);

      return thits;

    }
    // print "secondary, level %d : %d hits" % (m, hits)

    thits += hits;    
            
    poly_mul_ui(tmp_poly,gg,(unsigned long)dv0);
    poly_add(tmp_poly,ff,tmp_poly);
    poly_div_ui_mod_ui(nff,tmp_poly,p,p);    
    poly_mul(twist_f,dphi,gg);				
    new_twist_v1 = p * twist_v1;
    
    if(DEBUG) {
      printf("New twisted polynomial :\n");
      poly_print(twist_f);
      printf("New lattice twist : new_twist_v1=%ld\n",new_twist_v1);
    }
    
    if (m<rdict->m0 && m<rdict->maxpow) {
      for ( du=0 ; du<p ; du++) {      

	if(DEBUG) {
	  printf("Recalling (Rotation Inner)...\n");
	}
	hits = rotation_inner(rdict,p,nff,gg,fdg_gdf,u1,v1,l0,ld,m+1,new_twist_v1,scale,dphi);
	
	// print "secondary/excep, level %d : %d hits" % (m, hits)
	if(DEBUG) {
	  printf("Rotation Inner returned.\n",hits);
	}
	thits += hits;
	u1 += pm;
	v1 += twist_v1;      
	poly_add(nff,nff,twist_f);
	
	if(DEBUG) {
	  printf("New lattice : u1=%ld, v1=%ld, and the twisted polynomial is\n",u1,v1);
	  poly_print(nff);
	}
      }
    }
    else {
	
      if(DEBUG) {	  

	if(m == rdict->m0)
	  printf("Rectangle aborted ");
	  
	if(m == rdict->maxpow)
	  printf("Power aborted");
	  
	printf("\n");
      }	
    }
     
    if(DEBUG) {
      printf("Completed simple root section: exiting...\n");
    } 
    
  }
  poly_free(minus_f_over_g_modp);
  poly_free(nfdg_gdf);
  poly_free(nff);
  poly_free(ngg);
  poly_free(ndphi);
  poly_free(twist_f);
  poly_free(tmp_poly);
  mpz_clear(rhs);

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

  // We ensure that the lattice basis is positively oriented.
  if(us<0)
    us = -us;
  if(vs<0) {
    vs = -vs;
    skew = -skew;
  }
  while(skew<0)
    skew += us;  
  while(v0<V0) {
    v0 += vs;
    u0 += skew;
  }
  while(u0<U0) 
    u0 += us;  
  
  ASSERT_ALWAYS(us>=0 && skew >= 0 && vs>=0 && u0>=U0 && v0>=V0);  

  if(DEBUG) {
    printf("Lighting lattice u0=%ld ,v0=%ld, us=%ld, vs=%ld, skew=%ld \n",u0,v0,us,vs,skew);
  }

  q = (long)floor((v0-V0)/vs);
  v = v0-q*vs;
  u = u0-q*skew;
  q = (long)floor((u-U0)/us);
  u = u - q*us;
  pos = (U1-U0+1)*(v-V0) + (u-U0);

  if(DEBUG)
    printf("We start at point u=%ld,v=%ld \n",u,v);
  
  // at this stage, (u,v) is the lower left point of the lattice
  while(v<=V1) {    
    while(u<=U1) {       
      if(DEBUG) {
	printf("Sieving %ld,%ld, position %ld. \n",u,v,pos);
      }
      ASSERT(pos >= 0);      
      ASSERT(pos <(U1-U0+1)*(V1-V0+1));
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
  
  if(DEBUG)
    printf("There were %ld hits.\n",hits);
  
  
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
void poly_set_identity(poly_t f) {
  cleandeg(f,1);
  mpz_set_ui(f->coeff[0],0);
  mpz_set_ui(f->coeff[1],1);
}

// psi(x) = l + p*x
// f <- g(psi) 
// f must be allocated, f and g can be the same. 
// TODO : Ensure that no constraints on the degree of f are necessary.
void compose_psi(poly_t f, const poly_t g, unsigned long l, mpz_t * ppow, int maxpow) {
  
  int d;  
  d = g->deg;

  // Suspicious case
  if(d<=0) {
    if(DEBUG) {
      printf("(Compose psi) Warning : received constant polynomial of degree %d\n",d);
    }
    poly_copy(f,g);
    return;
  }

  ASSERT_ALWAYS(d<=maxpow); // Otherwise we will get a segfault.
  ASSERT_ALWAYS(g->alloc >= d+1); // Idem
  

  unsigned long i,k;
  mpz_t tmpsum,tmpterm,tmplpow;
  poly_t aux;
  poly_alloc(aux,d);


  mpz_init(tmpsum);
  mpz_init(tmpterm);
  mpz_init(tmplpow);    

  //printf("Degree of polynomial to be reduced is %d ...\n",d);
  //printf("(Compose Psi) Entering compose loop ...\n");

  for( k=0 ; k<=d ; k++) {

    mpz_set_ui(tmplpow,1);
    mpz_set_ui(tmpsum,0);
       
    for( i=k ; i<=d ; i++) {
      // printf("Computing term %lu of the sum...\n",i);
      // Compute the thing
      mpz_bin_uiui(tmpterm,i,k);
      //printf("Before coeff...\n");
      mpz_mul(tmpterm,tmpterm,g->coeff[i]);
      //printf("Before ppow...\n");      
      mpz_mul(tmpterm,tmpterm,ppow[k]);
      //mpz_mul(tmpterm,tmpterm,rdict->ppow[k]);
      //printf("After both\n");
      mpz_mul(tmpterm,tmpterm,tmplpow);
      
      // Update pows
      mpz_mul_ui(tmplpow,tmplpow,l);      
      
      // Add to the sum
      mpz_add(tmpsum,tmpsum,tmpterm);      
    }
    //printf("Computed the sum...\n");
    //printf("(Compose Psi) Coeff %lu of aux is \n",k);    
    //mpz_out_str(stdout,10,tmpsum);
    //printf("\n");
    poly_setcoeff(aux,k,tmpsum);
  }

  //printf("(Compose Psi) The aux polynomial, before copying into f :\n");
  //poly_print(aux);
  poly_copy(f,aux);

  //printf("(Compose Psi) The output's f polynomial :\n");
  //poly_print(f);

  mpz_clear(tmpsum);
  mpz_clear(tmpterm);
  mpz_clear(tmplpow);  
  poly_free(aux);

}
/* We assume (rdict->delta) > 0. No problem when f and g are the same. */
void compose_reduce(poly_t f, const poly_t g,unsigned long l,mpz_t * ppow, int maxpow, mpz_t m) {
  poly_t aux;
  poly_alloc(aux,g->deg);
  compose_psi(aux,g,l,ppow,maxpow);
  poly_reduce_mod_mpz(f,aux,m);
  //printf("(Compose reduce) The original and reduced polynomials are\n");
  //poly_print(g);
  //printf("(Compose reduce) and\n");
  //poly_print(f);
  poly_free(aux);
}



