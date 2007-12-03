#ifndef TICKS_HPP_
#define TICKS_HPP_

#ifdef	__linux__
extern unsigned long linux_tid();	/* This is only a system call */
#endif

/* give the number of running threads as an argument */
/* This will bark if a thread is not running presently */
extern double thread_ticks(int = 1);

extern double oncpu_ticks();

extern double wallclock_ticks();

#endif	/* TICKS_HPP_ */
