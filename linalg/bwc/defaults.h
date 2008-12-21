#ifndef DEFAULTS_H_
#define DEFAULTS_H_

/* All the defaults in this file can be tuned by giving the proper
 * arguments to prep (for check_interval) and krylov (for the other two).
 */

/* This is the periodicity for checking dotproducts. Set to 0 to disable
 * checking altogether. */
#define CHECK_INTERVAL_DEFAULT  100

/* This is the periodicity for gathering and flushing the A files. This
 * may become unsignificant if v_save_interval is below that. The A files are
 * written in append mode only, so the amount per extra iteration is tiny
 * (m*n times the size of a scalar).
 */
#define A_FLUSH_INTERVAL_DEFAULT        100

/* This is the periodicity for writing out the intermediate V files. This
 * files are somewhat big, so the interval should not be too short. Since
 * V files are rotated on each write, there is no impact on the total
 * disk usage, though.  This periodicity obviously affects the
 * checkpointing capability of the program.
 */
#define V_SAVE_INTERVAL_DEFAULT        1000


#endif	/* DEFAULTS_H_ */
