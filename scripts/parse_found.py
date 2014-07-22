#!/usr/bin/env python3

import re
import sys
import math
import datetime

# Relations counts received within each interval are added up.
# Define here how long an interval is.
SECONDS_PER_INTERVAL = 7200

# Format required to convert logger messages into datetime objects
DATE_FMT = "%Y-%m-%d %H:%M:%S,%f"
# The time stamp of when we started sieving RSA220
START_TIME = datetime.datetime.strptime("2014-05-17 22:10:24,702", DATE_FMT)

# Regex to parse relation counts as printed in the log. Example:
# "PID13816 2014-05-27 13:23:24,843 Info:Lattice Sieving: Found 7819 relations in '/localdisk/kruppaal/factoring/rsa220/rsa220.upload/rsa220.sieving.985030000-985040000.rlat6a.gz', total is now 156535983/1500000000"
FOUND_RE = re.compile(r"PID\d+ (\d+-\d+-\d+ \d+:\d+:\d+,\d+) Info:Lattice Sieving: Found (\d+) relations in '\S+', total is now (\d+)/\d+")

sums = {}

for line in sys.stdin:
  # print(line.strip())
    match = FOUND_RE.match(line)
    if match:
        (timestamp, new_rels, rels) = match.groups()
        time = datetime.datetime.strptime(timestamp, DATE_FMT)
        sec = (time - START_TIME).total_seconds()
        if False:
            print("%s %s" % (sec, rels))
        else:
            interval = math.floor(sec / SECONDS_PER_INTERVAL)
            sums[interval] = sums.get(interval, 0) + int(new_rels)

cumulative_value = 0
for i in range(max(sums)): # don't do +1 to avoid partial interval
    value = sums[i] if i in sums else 0
    cumulative_value += value
    print("%d %d %d" % (i, value, cumulative_value))
