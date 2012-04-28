#!/usr/bin/env sage 

import sys

load tools.sage
load readparam.sage
load makefb0.sage

if len(sys.argv) != 2:
    print "Usage: %s <paramfile>"%sys.argv[0]
    sys.exit(1)

print "Constructing factor base for side 0...\n"
makefb(sys.argv[1], 0)
print "Constructing factor base for side 1...\n"
makefb(sys.argv[1], 1)

