#!/bin/sh

# a.out is the c program below. See also taskset.
# This silly script clearly lacks integrating ;-))


bind_threads() {
    p=$1
    shift
    pid=`pidof $p`
    ls /proc/$pid/task | grep -v $pid | while read tpid ; do
        ~/a.out $tpid $1
        shift
    done
}

# Usage: ./bind_threads bw-krylov-mt 1 2 5 6

bind_threads "$@"


# #define _GNU_SOURCE
# #include <sched.h>
# #include <stdio.h>
# #include <stdlib.h>
# #include <string.h>
# #include <unistd.h>
# #include <errno.h>
# 
# int main(int argc, char * argv[])
# {
#     cpu_set_t mask[1];
#     CPU_ZERO(mask);
# 
#     pid_t pid = atoi(argv[1]);
# 
#     for(int i = 2 ; i < argc ; i++) {
#         CPU_SET(atoi(argv[i]), mask);
#     }
# 
#     int rc = sched_setaffinity(pid, 8, mask);
# 
#     if (rc < 0) {
#         fprintf(stderr, "%s\n", strerror(errno));
#         exit(1);
#     }
#     return 0;
# }
