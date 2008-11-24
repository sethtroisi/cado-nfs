#!/bin/sh

# This script can be used to run an mpi program with one xterm per job,
# and wraps the program as a child of some ``controlling application''.
# The controlling application is specified as the --tool option, and may
# be gdb (default), or valgrind (or anything you want, provided you
# modify the script accordingly).

# For a standard job, one would change ./a.out <params> to gdb --args ./a.out
# <params>, or valgrind ./a.out <params>.
#
# Here the pattern goes as the typical command:
#       mpirun -np NN <mpi args> ./a.out <args>
# being changed into
#       mpirun -np NN <mpi args> ./tool.sh --tool gdb ./a.out <args>
# or
#       mpirun -np NN <mpi args> ./tool.sh --tool valgrind ./a.out <args>

# Note that per se, this script does not deal with mpi, so it's a generic
# way of running a program controlled by gdb/valgrind within a separate
# xterm...

tool=run

if [ $1 = "--inside" ] ; then
    # --inside is an internal parameter, do not use.
    shift
else
    xterm -e $0 --inside "$@"
    exit 0
fi

if [ $1 = "--tool" ] ; then
    shift
    tool=$1
    shift
fi

env | grep OMPI

CAN_PRINT=1
export CAN_PRINT

if [ $tool = "gdb" ] ; then
    gdb --args "$@"
    exit 0
elif [ $tool = "xgdb" ] ; then
    gdb -ex run --args "$@"
    exit 0
elif [ $tool = "memcheck" ] || [ $tool = "valgrind" ] ; then
    valgrind $@
    echo 'press ctrl-D'
    cat
elif [ $tool = "run" ] ; then
    $@
    echo 'press ctrl-D'
    cat
else
    echo "Unknown tool $tool, sorry" >&2
    echo 'press ctrl-D'
    cat
fi
