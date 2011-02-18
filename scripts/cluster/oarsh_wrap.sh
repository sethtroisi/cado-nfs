#!/bin/sh

exec oarsh -x -q -oUserKnownHostsFile=/dev/null  "$@"
