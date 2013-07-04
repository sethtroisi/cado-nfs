#!/usr/bin/env python
import sys

for line in sys.stdin:
  line2 = line.split('#',1)[0].strip()
  if line2:
    print(line2)
