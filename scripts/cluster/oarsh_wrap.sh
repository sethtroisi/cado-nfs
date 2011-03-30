#!/bin/sh

unset DISPLAY
exec oarsh -x -q -oUserKnownHostsFile=/dev/null  "$@"
