#!/usr/bin/env bash

SCRIPTNAME="$0"

function usage() {
  echo "Usage: $SCRIPTNAME [-s <n>] [-p <path>] <parameters for wuclient2.py>"
  echo "Mandatory parameters for wuclient2.py: --server=<URL>"
  echo "Optional parameters:"
  echo "-s <n>: Start <n> client scripts on each host. Default: 1"
  echo "-p <path>: Look for wuclient2.py script in directory specified by <path>."
  echo "           Default: ./"
}


if [ -z "$OAR_NODE_FILE" ]
then
  echo "OAR_NODE_FILE shell environment variable not set. $SCRIPTNAME is meant to be run inside an OAR job." >&2
  usage
  exit 1
fi

if [ ! -f "$OAR_NODE_FILE" ]
then
  echo "Variable OAR_NODE_FILE is set, but file $OAR_NODE_FILE does not exist. That's odd." >&2
  usage
  exit 1
fi

# If an -s parameter is given, start that many clients per node

NRCLIENTS=1
if [ "$1" = "-s" ]
then
  NRCLIENTS="$2"
  shift 2
fi

SCRIPTPATH="."
if [ "$1" = "-p" ]
then
  SCRIPTPATH="$2"
  shift 2
fi

HAVE_SERVER=false
for ARG in "$@"
do
  if [[ "$ARG" =~ "--server" ]]
  then
    HAVE_SERVER=true
  fi
done
if ! $HAVE_SERVER
then
  echo "You must specify the URL of the workunit server with --server=<URL>"
  usage
  exit 1
fi


NODES=`cat "$OAR_NODE_FILE" | uniq | tr '\n' ' '`

echo "Starting $NRCLIENTS clients each on $NODES with parameters: $@"

for NODE in `cat "$OAR_NODE_FILE" | uniq`
do 
    for I in `seq 1 "$NRCLIENTS"`
    do 
        oarsh "$NODE" "$SCRIPTPATH/wuclient2.py" -d "$@"
    done
done
