#ifndef PAD_H_
#define PAD_H_

#ifndef PAD
#define PAD_(X,Y) X ## _ ## Y
#define PAD(X,Y) PAD_(X,Y)
#endif
#ifdef  P
#undef P
#endif

#endif	/* PAD_H_ */
