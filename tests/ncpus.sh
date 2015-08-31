#!/usr/bin/env bash

set -ex

if [ "$NCPUS_FAKE" ] ; then
    echo $NCPUS_FAKE
    exit 0
fi

if [ -f /proc/cpuinfo ] ; then
    # grep -c fails with no match, while |wc -l is happy
    nphysical=$(sort -u /proc/cpuinfo  | grep '^physical' | wc -l)
    if [ "$nphysical" -eq 0 ] ; then
        grep -c ^processor /proc/cpuinfo
        exit 0
    fi
    variants="$(grep '^cpu cores' /proc/cpuinfo | uniq | wc -l)"
    if [ "$variants" = 0 ] ; then
        grep -c ^processor /proc/cpuinfo
        exit 0
    fi
    if [ "$variants" != 1 ] ; then
        echo "inhomogeneous platform ?" >&2
        exit 1
    fi
    cores_per_cpu=$(grep '^cpu cores' /proc/cpuinfo | head -1 | cut -d: -f2)
    echo $((nphysical*cores_per_cpu))
elif [ "$(uname -s)" = Darwin ] ; then
    # does this count hyperthreading or not ?
    sysctl -n hw.ncpu
elif [ "$(uname -s)" = OpenBSD ] ; then
    # does this count hyperthreading or not ?
    sysctl -n hw.ncpu
elif [ "$(uname -s)" = MINGW32_NT-6.1 ] ; then
    # not clear whether it's physical or logical.
    # wmic cpu get Caption | tail -n +2 | grep -c .
    # anyway we don't believe mingw will support multithreading well. For
    # sure make -j2 seems to just not work...
    echo 1
else
    # this would work as well on linux and darwin, but pretty surely does
    # not count hyperthreading. Does not work on openbsd5.3
    getconf _NPROCESSORS_ONLN
fi

