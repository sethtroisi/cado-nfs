import sys
import subprocess
import re

sq_pattern = re.compile(r"# Sieving algebraic q=(\d+); rho=(\d+);")
dupecheck_pattern = re.compile(r"# DUPECHECK Checking if relation \(a,b\) = \((-?\d+),(\d+)\) is a dupe of sieving special-q q=(\d+); rho=(\d+)")
isdupe_pattern = re.compile("# DUPECHECK relation is probably a dupe")
notdupe_pattern = re.compile("# DUPECHECK relation is probably not a dupe")

def run_las(a, b, q, rho, cmdline):
  assert not "-q0" in cmdline
  assert not "-q1" in cmdline
  assert not "-rho" in cmdline
  my_cmdline = cmdline + ["-q0", str(q), "-rho", str(rho)]
  # print("Running: %s" % " ".join(my_cmdline))
  p = subprocess.Popen(my_cmdline, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  (stdout, stderr) = p.communicate()
  stdout = stdout.decode("ascii")
  # print("Stdout: %s" % stdout)
  found = bool(re.search("^%s,%s:" % (a,b), stdout, re.MULTILINE))
  return found

def pop_param(arr, param, must_occur=False):
  if param in arr:
    i = arr.index(param)
    v = arr[i+1]
    del(arr[i:i+2])
    return v
  elif must_occur:
    raise Exception("Required parameter %s did not occur" % param)

for filename in sys.argv[1:]:
  isdupe = None
  f = open(filename, "r")
  cmdline = f.readline().split()
  if not cmdline[0] == "#" or not re.match("\(.*\)$", cmdline[1]):
    print("Command line not found in first line of %s" % filename)
    sys.exit(1)
  cmdline = cmdline[2:]

  q0 = int(pop_param(cmdline, "-q0", True))
  q1 = pop_param(cmdline, "-q1")
  q1 = None if q1 is None else int(q1)
  rho = pop_param(cmdline, "-rho")
  rho = None if rho is None else int(rho)
  pop_param(cmdline, "-out")
  
  for line in f:
    match = sq_pattern.match(line)
    if match:
      cur_q, cur_rho = match.groups()
    match = dupecheck_pattern.match(line)
    if match:
      # print(line)
      assert isdupe is None
      (a, b, q, rho) = match.groups()
    
    match = isdupe_pattern.match(line)
    if match:
      isdupe = True
    
    match = notdupe_pattern.match(line)
    if match:
      isdupe = False
    if not isdupe is None:
      found = run_las(a, b, q, rho, cmdline)
      errormsg = "" if isdupe == found else "Error: "
      print("%sRelation %s,%s in q,rho=%s,%s was %sconsidered a dupe and is %sfound by q=%s, rho=%s" % 
            (errormsg, a, b, cur_q, cur_rho, "" if isdupe else "not ", "" if found else "not ", q, rho));
      isdupe = None
