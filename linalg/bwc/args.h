#ifndef ARGS_H_
#define ARGS_H_

#ifdef __cplusplus
extern "C" {
#endif

extern int argparse(int * p_argc, char *** p_argv);
extern int finish();

#ifdef __cplusplus
}
#endif

#endif	/* ARGS_H_ */
