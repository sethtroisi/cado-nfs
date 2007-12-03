#ifndef BW_LVBLOCK_STEPS_H_
#define BW_LVBLOCK_STEPS_H_
/* l : traveling from 0 to nbys-1 */
/* k : traveling from 0 to bw_longsize */
/* j : traveling from 0 to ncols */

#define bw_lvblock_step_p00(x)  x+=bw_longsize*nbys
#define bw_lvblock_step_n00(x,n)  x+=n*bw_longsize*nbys
#define bw_lvblock_step_z00(x)  x-=ncols*bw_longsize*nbys

#ifdef BLOCKS_TOGETHER
#define bw_lvblock_step_nnn(x,j,k,l)	x+=l+nbys*(k+bw_longsize*j)
/* travel through l,k,j in this order */
#define bw_lvblock_step_00p(x)	x++				/* nnn(0,0,1) */
#define bw_lvblock_step_0pz(x)					/* nnn(0,1,-) */
#define bw_lvblock_step_pz0(x)					/* nnn(1,-,0) */

/* travel through k,l,j */
#define bw_lvblock_step_0p0(x)	x+=nbys				/* nnn(0,1,0) */
#define bw_lvblock_step_0n0(x,n) x+=n*nbys			/* nnn(0,n,0) */
#define bw_lvblock_step_0z0(x)	x-=nbys*bw_longsize		/* nnn(0,-,0) */
#define bw_lvblock_step_0zp(x)	x+=1-bw_longsize*nbys 		/* nnn(0,-,1) */
#define bw_lvblock_step_p0z(x)	x+=(bw_longsize-1)*nbys		/* nnn(1,0,-) */
#define bw_lvblock_step_00z(x)	x-=nbys				/* nnn(0,0,-) */

#else
#define bw_lvblock_step_nnn(x,j,k,l)	x+=k+bw_longsize*(l+nbys*j)
/* travel through k,l,j */
#define bw_lvblock_step_0p0(x)	x++				/* nnn(0,1,0) */
#define bw_lvblock_step_0n0(x,n) x+=n				/* nnn(0,1,0) */
#define bw_lvblock_step_0z0(x)	x-=bw_longsize			/* nnn(0,-,0) */
#define bw_lvblock_step_0zp(x)					/* nnn(0,-,1) */
#define bw_lvblock_step_p0z(x)					/* nnn(1,0,-) */
#define bw_lvblock_step_00z(x)	x-=nbys*bw_longsize		/* nnn(0,0,-) */

/* travel through l,k,j */
#define bw_lvblock_step_00p(x)	x+=bw_longsize			/* nnn(0,0,1) */
#define bw_lvblock_step_0pz(x)	x+=1-bw_longsize*nbys		/* nnn(0,1,-) */
#define bw_lvblock_step_pz0(x)	x+=(nbys-1)*bw_longsize		/* nnn(1,-,0) */
#endif

#endif /* BW_LVBLOCK_STEPS_H_ */
