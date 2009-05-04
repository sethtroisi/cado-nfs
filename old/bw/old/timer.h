#ifndef TIMER_H_
#define TIMER_H_

#ifdef	__cplusplus
extern "C" {
#endif

#define TIMER_SET	1
#define TIMER_ASK	2
#define TIMER_ASKSET	(TIMER_ASK | TIMER_SET)
#define TIMER_MTN(x)	(((x)-1)<<2)
#define TIMER_MT	TIMER_MTN(1024)
#define TIMER_NT(x)	(1+((x)>>2))

extern double timer(int);
extern double timer_r(struct timeval *,int);
extern double wall_clock(int);
extern double wall_clock_r(struct timeval *, int);

#ifdef	__cplusplus
}
#endif

#endif	/* TIMER_H_ */
