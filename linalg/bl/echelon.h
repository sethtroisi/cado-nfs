extern void ConstructPerm(unsigned long n, unsigned long US[n + 1],
			  unsigned long **p);
extern void Construction_S_Winv(unsigned long n, unsigned long **m,
				unsigned long USet[n + 1], unsigned long **S,
				unsigned long Set0[1], unsigned long **Winv,
				unsigned long **SGPerm);
extern void ConstructPermBit(unsigned long n, unsigned long *US,
			     unsigned long *p);
void Construction_S_WinvBit(unsigned long n, unsigned long *m,
			    unsigned long *USet, unsigned long *S,
			    unsigned long *Set0, unsigned long *Winv,
			    unsigned long *SGPerm);
