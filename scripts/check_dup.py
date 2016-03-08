#!/usr/bin/env python3
import sys
import subprocess
import re

sq_pattern = re.compile(r"# Sieving algebraic q=(\d+); rho=(\d+);")
dupecheck_pattern = re.compile(r"# DUPECHECK Checking if relation \(a,b\) = \((-?\d+),(\d+)\) is a dupe of sieving special-q -q0 (\d+) -rho (\d+)")
isdupe_pattern = re.compile("# DUPECHECK relation is probably a dupe")
notdupe_pattern = re.compile("# DUPECHECK relation is probably not a dupe")
oldIJ_pattern = re.compile(r"# DUPECHECK oldI = \d+, I = (\d+), oldJ = \d+, J = (\d+)")
old_ij_pattern = re.compile(r"# DUPECHECK relation had i=(-?\d+), j=(\d+), remaining lognorms (\d+), (\d+)")
IJ_pattern = re.compile(r"# I=(\d+); J=(\d+)")

def run_las(a, b, q, rho, cmdline, I=None, J=None, ij_lognorm=None):
  sys.stdout.write("Re-running special-q %d,%d, looking for relation %d,%d" % (q, rho, a, b))
  if ij_lognorm is not None:
    sys.stdout.write(" with i,j,lognorm0,lognorm1 = %d,%d,%d,%d\n" % tuple(ij_lognorm))
  else:
    sys.stdout.write("\n")
    
  assert not "-q0" in cmdline
  assert not "-q1" in cmdline
  assert not "-rho" in cmdline
  my_cmdline = cmdline + ["-q0", str(q), "-rho", str(rho)]
  # print("Running: %s" % " ".join(my_cmdline))
  p = subprocess.Popen(my_cmdline, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  (stdout, stderr) = p.communicate()
  stdout = stdout.decode("ascii")
  # print("Stdout: %s" % stdout)

  # Verify that the I,J-values agree
  match = IJ_pattern.search(stdout, re.MULTILINE)
  if match and I is not None and J is not None:
    (dupe_I, dupe_J) = map(int, match.groups())
    # print("Found dupe_I=%d, dupe_J=%d" % (dupe_I, dupe_J))
    if I != dupe_I or J != dupe_J:
      print("Error, re-run used different I,J = %d,%d, dupecheck used %d,%d" %
            (dupe_I, dupe_J, I, J))

  # Check that the i,j-coordinates and remaining lognorms agree
  if ij_lognorm is not None:
    pattern = r"# i=(-?\d+), j=(\d+), lognorms = (\d+), (\d+)\s%d,%d:" % (a,b)
    match = re.search(pattern, stdout, re.MULTILINE)
    if match:
      new_ij_lognorm = tuple(map(int, match.groups()))
      ij = tuple(ij_lognorm[0:2])
      new_ij = new_ij_lognorm[0:2]
      if ij != new_ij:
        print("Error, re-run had different i,j = %d,%d, dupecheck used %d,%d"
              % (new_ij + ij))
      lognorm = tuple(ij_lognorm[2:4])
      new_lognorm = new_ij_lognorm[2:4]
      if lognorm != new_lognorm:
        print("Warning, re-run had different lognorms: %d,%d, dupecheck used %d,%d"
              % (new_lognorm + lognorm))
    else:
      print("Warning, no i,j,lognorm0,lognorm1 = %d,%d,%d,%d received from re-run")
  
  found = bool(re.search("^%s,%s:" % (a,b), stdout, re.MULTILINE))
  return found

def pop_param(arr, param, must_occur=False, has_value=True):
  if param in arr:
    i = arr.index(param)
    if has_value:
      v = arr[i+1]
      del(arr[i:i+2])
      return v
    else:
      del(arr[i])
      return True
  elif must_occur:
    raise Exception("Required parameter %s did not occur" % param)
  else:
    return None

for filename in sys.argv[1:]:
  isdupe = None
  I = None
  J = None
  old_ij_lognorm = None
  f = open(filename, "r")
  cmdline = f.readline().split()
  if not cmdline[0] == "#" or not re.match("\(.*\)$", cmdline[1]):
    print("Command line not found in first line of %s" % filename)
    sys.exit(1)
  cmdline = cmdline[2:]

  q0 = int(pop_param(cmdline, "-q0", must_occur=True))
  q1 = pop_param(cmdline, "-q1")
  q1 = None if q1 is None else int(q1)
  rho = pop_param(cmdline, "-rho")
  rho = None if rho is None else int(rho)
  pop_param(cmdline, "-out")
  pop_param(cmdline, "-dup", has_value=False)
  while pop_param(cmdline, "-v", has_value=False):
    pass
  cmdline += ["-v"] * 2
  
  print("Using command line template: %s" % " ".join(cmdline))
  
  for line in f:
    match = sq_pattern.match(line)
    if match:
      cur_q, cur_rho = map(int, match.groups())
    match = dupecheck_pattern.match(line)
    if match:
      # print(line)
      assert isdupe is None
      (a, b, q, rho) = map(int, match.groups())
    
    match = oldIJ_pattern.match(line)
    if match:
      I, J = map(int, match.groups())
      # print("Found: I=%d, J=%d" % (I, J))

    match = old_ij_pattern.match(line)
    if match:
      old_ij_lognorm = list(map(int, match.groups()))
    
    match = isdupe_pattern.match(line)
    if match:
      isdupe = True
    
    match = notdupe_pattern.match(line)
    if match:
      isdupe = False
    if not isdupe is None:
      if old_ij_lognorm is None:
          print("Warning, no i,j,lognorm0,lognorm1 = %d,%d,%d,%d received from dupecheck")
      found = run_las(a, b, q, rho, cmdline, I, J, old_ij_lognorm)
      errormsg = "" if isdupe == found else "Error: "
      print("%sRelation %s,%s in q,rho=%s,%s was %sconsidered a dupe and is %sfound by q=%s, rho=%s" % 
            (errormsg, a, b, cur_q, cur_rho, "" if isdupe else "not ", "" if found else "not ", q, rho));
      isdupe = None
      I = None
      J = None
      old_ij_lognorm = None
