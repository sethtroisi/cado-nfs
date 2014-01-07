#!/usr/bin/env python3

# Generates gnuplot commands to plot the sieve region in the a,b-plane, and the relations found.
# One .png file is generated for each special-q in the siever output file.
# Takes one command line parameter with the siever output file; if not specified, reads from stdin

import sys
import re

pattern1 = re.compile(r"# Sieving algebraic q=(\d+); rho=(\d+); a0=(-?\d+); b0=(-?\d+); a1=(-?\d+); b1=(-?\d+);")
pattern2 = re.compile(r"# (\d+) relation\(s\) for algebraic")
pattern3 = re.compile(r"(-?\d+),(\d+):")
pattern4 = re.compile(r"# I=(\d+); J=(\d+)")

def print_gnuplot(q, r, a0, b0, a1, b1, i, j, rels, outputfilename):
  # First corner: (-I,0) -> (-I * a0, -I * b0)
  c1x = -I//2 * a0
  c1y = -I//2 * b0
  # Second corner: (-I,J) -> (-I * a0 + J * a1, -I * b0 + J * b1)
  c2x = -I//2 * a0 + J * a1
  c2y = -I//2 * b0 + J * b1
  # Third corner: (I,J) -> (I * a0 + J * a1, I * b0 + J * b1)
  c3x = I//2 * a0 + J * a1
  c3y = I//2 * b0 + J * b1
  # Fourth corner: (I,0) -> (I * a0, I * b0)
  c4x = I//2 * a0
  c4y = I//2 * b0

  xmin=min(c1x, c2x, c3x, c4x)
  xmax=max(c1x, c2x, c3x, c4x)
  ymin=min(c1y, c2y, c3y, c4y)
  ymax=max(c1y, c2y, c3y, c4y)

  if ymin < 0:
    xmin = min(xmin, -xmax)
    xmax = max(xmax, -xmin)
    ymax = max(ymax, -ymin)

  print('set terminal png size 1200,900 enhanced font "Helvetica,20"')
  print("set output 'sieve_region.%d.%d.png'" % (q, r))

  print("set xrange [%d:%d]" % (xmin, xmax))
  print("set yrange [%d:%d]" % (ymin, ymax))
  print("set xzeroaxis")
  print("set yzeroaxis")
  print('set title "q=%d, r=%d, I=%d, J=%d, yield: %d"' % (q, r, i, j, rels))
  print("set object 1 polygon from %d,%d to %d,%d to %d,%d to %d,%d to %d,%d"
        % (c1x, c1y, c2x, c2y, c3x, c3y, c4x, c4y, c1x, c1y))
  print('set label 1 "(-I/2,0)" at %d,%d' % (c1x, c1y))
  print('set label 2 "(-I/2,J)" at %d,%d' % (c2x, c2y))
  print('set label 3 "(I/2,J)" at %d,%d' % (c3x, c3y))
  print('set label 4 "(I/2,0)" at %d,%d' % (c4x, c4y))

  print('plot "%s"' % outputfilename)


OUTPUTFILE = None
OUTPUTFILENAME = None

I=2048
J=1024

if len(sys.argv) > 1:
  INPUTFILE = open(sys.argv[1])
else:
  INPUTFILE = sys.stdin

for line in INPUTFILE:
  match = pattern1.match(line)
  if match:
    (q, r, a0, b0, a1, b1) = map(int, match.groups())
    OUTPUTFILENAME = "relations_ab.%d.%d" % (q, r)
    OUTPUTFILE = open(OUTPUTFILENAME, "w")
    continue
  match = pattern2.match(line)
  if match:
    rels = int(match.group(1))
    if OUTPUTFILE:
      OUTPUTFILE.close()
    print_gnuplot(q, r, a0, b0, a1, b1, I, J, rels, OUTPUTFILENAME)
    continue
  match = pattern3.match(line)
  if match:
    (a, b) = map(int, match.groups())
    OUTPUTFILE.write("%d %d\n" % (a,b))
    continue
  match = pattern4.match(line)
  if match:
    (I, J) = map(int, match.groups())

if OUTPUTFILE:
  OUTPUTFILE.close()
