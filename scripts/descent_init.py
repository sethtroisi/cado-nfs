#!/usr/bin/env python3
#
# Example for the p59 in tests/test_full_p59
# scripts/descent_init.py
#   --p 43341748620473677010074177283795146221310971425909898235183
#   --target 15384226356769205532939362866999574621778419520781718884300
#   --cadobindir /tmp/cado-nfs-build-master
#   --mfb 40 --lpb 25 --lim 1000000 --tkewness 8388608 --I 10
#   --ncurves 20

import os
import subprocess
import sys
import argparse
import re
import math
import tempfile
import shutil

def Myxgcd(a, b, T):
    ainit = a
    binit = b
    bound = math.floor(math.sqrt(b*T))
    x = 0
    lastx = 1
    y = 1
    lasty = 0
    while abs(b) > bound:
        q = a // b
        r = a % b
        a = b
        b = r
        newx = lastx - q*x
        lastx = x
        x = newx
        newy = lasty - q*y
        lasty = y
        y = newy
    return [ [ b, x ], [ a, lastx ] ]

def get_root_modp(polyfile, p):
    G = [0, 0]
    with open(polyfile, "r") as file:
        for line in file:
            if line[0] == 'Y':
                ll = line.split()
                if line[1] == '0':
                    G[0] = int(ll[1])
                else:
                    assert line[1] == '1'
                    G[1] = int(ll[1])
    assert G[0] != 0 and G[1] != 0
    invg1 = pow(G[1], p-2, p)
    r = (-G[0]*invg1) % p
    assert (G[0] + r*G[1]) % p == 0
    return r

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Descent initialization for DLP")
    # Required
    parser.add_argument("--p", help="Prime defining GF(p) in which we work",
            type=int, required=True)
    parser.add_argument("--target", help="Element whose DL is wanted",
            type=int, required=True)
    parser.add_argument("--cadobindir", help="Cado build directory",
            required=True, type=str)
    # Optional
    parser.add_argument("--tkewness", help="Tkewness",
            type=int, default=2**30)
    parser.add_argument("--workdir", help="Temporary working directory")
    parser.add_argument("--lim", help="Factor base bound", default=2**26)
    parser.add_argument("--lpb", help="Large prime bound", default=64)
    parser.add_argument("--mfb", help="Cofactor bound", default=100)
    parser.add_argument("--ncurves", help="ECM effort in cofac",
            default=80)
    parser.add_argument("--I", help="Sieving range", default=14)
    parser.add_argument("--nthreads", help="How many threads to use", default=4)
    parser.add_argument("--todofile", help="Output primes to this toto file",
            type=str)
    parser.add_argument("--poly", help="Polynomial file", type=str)
    parser.add_argument("--finallpb", help="Primes smaller than this bound are not printed in the todo file", type=int)
    args = parser.parse_args()
    if args.todofile != None and args.poly == None:
        print("--poly option is required if you give a todofile")
        parser.print_help()
        sys.exit(1)
    p = args.p
    z = args.target
    Tkew = args.tkewness
    gg = Myxgcd(z, p, Tkew)
    wdir = args.workdir
    if wdir == None:
        wdir = tempfile.mkdtemp(dir="/tmp")

    print ("Working directory is: " + wdir)
    polyfilename = wdir + "/desc.poly"
    polyfile = open(polyfilename, 'w')
    polyfile.write("n: " + str(p) + "\n")
    polyfile.write("skew: 1\n")
    polyfile.write("c1: " + str(gg[0][0]) + "\n")
    polyfile.write("c0: " + str(gg[1][0]) + "\n")
    polyfile.write("Y1: " + str(gg[0][1]) + "\n")
    polyfile.write("Y0: " + str(gg[1][1]) + "\n")
    polyfile.close()

    print ("Creating factor bases.")
    for side in range(2):
        subprocess.check_call([args.cadobindir + "/sieve/makefb",
            "-t", str(args.nthreads),
            "-poly", polyfilename,
            "-alim", str(args.lim),
            "-side", str(side),
            "-out", wdir + str("/desc.fb") + str(side)],
            stdout=subprocess.DEVNULL)
    
    print ("Sieving.")
    q0 = Tkew
    q1 = Tkew + 100000
    relfilename = wdir + "/desc.rels"
    subprocess.check_call([args.cadobindir + "/sieve/las",
        "-poly", polyfilename,
        "-fb0", wdir + str("/desc.fb") + str(0),
        "-fb1", wdir + str("/desc.fb") + str(1),
        "-lim0", str(args.lim),
        "-lim1", str(args.lim),
        "-lpb0", str(args.lpb),
        "-lpb1", str(args.lpb),
        "-mfb0", str(args.mfb),
        "-mfb1", str(args.mfb),
        "-ncurves0", str(args.ncurves),
        "-ncurves1", str(args.ncurves),
        "-I", str(args.I),
        "-q0", str(q0),
        "-q1", str(q1),
        "-exit-early",
        "-out", relfilename
        ])

    rel = None
    with open(relfilename, 'r') as relfile:
        for line in relfile:
            if line[0] == '#':
                continue;
            else:
                rel = line
                break
    if rel == None:
        raise NameError("No relation found!")
    rel = rel.split(':')
    ab = rel[0].split(',')
    a = int(ab[0])
    b = int(ab[1])
    Num = a*gg[0][0] + b*gg[1][0]
    Den = a*gg[0][1] + b*gg[1][1]
    assert (z*Den-Num) % p == 0

    factNum = [ int(x, 16) for x in rel[2].split(',') ]
    factDen = [ int(x, 16) for x in rel[1].split(',') ]
    NN = 1
    for x in factNum:
        NN = NN*x
    assert (NN == abs(Num))

    DD = 1
    for x in factDen:
        DD = DD*x
    assert (DD == abs(Den))
    print (Num, Den, factNum, factDen)

    if not args.todofile == None:
        finallpb = 20
        if args.finallpb != None:
            finallpb = args.finallpb
        with open(args.todofile, "w") as file:
            for p in factNum + factDen:
                if p > 1<<finallpb:
                    r = get_root_modp(args.poly, p)
                    file.write("0 " + str(p) + " " + str(r) + "\n")

    shutil.rmtree(wdir)

