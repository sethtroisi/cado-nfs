#!/bin/bash

# This script is executed by the jenkins connection. It has to do the
# following things.

# - Use the JOB_NAME to decide exactly which test is to be run.
# - check the persistent test database for the test, and see if one is
#   running or pending. If yes, hold.
# - submit one OAR job which will compile the revision with the
#   appropriate software and compilation flags
# - submit several OAR jobs which will perform the required checks (those
#   OAR jobs may require various numbers of nodes)

# As long as the script detects that a previous script has not completed,
# it will stand by, and the jenkins master will see this as just taking
# long to complete. So this means that the timeout value has to be
# adjusted accordingly, and that the number of nodes required by checks
# should be kept minimal.

set -e

# cd to where we are.
cd "`dirname $0`"

if ! [ "$JOB_NAME" ] ; then
    echo "This script must be run with \$JOB_NAME set" >&2
    exit 1
fi

if ! [ "$BUILD_TAG" ] ; then
    BUILD_TAG="$JOB_NAME-test-`date +%Y%m%d%H%M%S`"
fi

while [ $# -gt 0 ] ; do
    if false ; then
        shift
    else
        echo "Unexpected argument $1" >&2
        exit 1
    fi
done

SRC="$HOME/jenkins/workspace/$JOB_NAME"
CONTROL="$PWD/$JOB_NAME"
LOGFILEDIR="$CONTROL/logfiles"
JOBDIR="$CONTROL/jobs.d"

lockfile="$CONTROL/pending"
deadlockfile="$LOGFILEDIR/$BUILD_TAG.log"

mkdir -p "$LOGFILEDIR" "$JOBDIR" "$CONTROL" || :
if ! [ -f "$SRC/local.sh" ] && [ -f "$CONTROL"/local.sh ] ; then
    cp -f "$CONTROL"/local.sh "$SRC"/local.sh
fi
if ! [ -f "$SRC/local.sh" ] && [ -f local.sh.$JOB_NAME ] ; then
    cp -f local.sh.$JOB_NAME "$SRC"/local.sh
fi
if ! [ -f "$SRC/local.sh" ] && [ -f local.sh ] ; then
    cp -f local.sh "$SRC"/local.sh
fi
if ! [ -f "$JOBDIR/00build" ] ; then
    echo "Creating default job file $JOBDIR/00build" >&2
    cat > "$JOBDIR/00build" <<-'EOF'
#!/bin/bash
#OAR -l /nodes=1,walltime=1

set -e
set -x
make cmake
make gf2x-build
make -j$(grep -c ^processor /proc/cpuinfo)
EOF
    chmod 755 "$JOBDIR/00build"
fi

cd "$SRC"


if [ -e "$lockfile" ] ; then
    echo "File $lockfile indicates lock, owner: $(head -1 $lockfile)" >&2
    while [ -e "$lockfile" ] ; do sleep 10 ; echo -n . >&2 ; done ; echo >&2
else
    echo "`date` $BUILD_TAG" > $lockfile
fi

# ok, from this point on, we should not rely simply on -e, since
# otherwise we'll quite likely leave something lacking cleanup every once
# in a while...

set +e

progress() {
    echo "$@" >&2
    echo "`date` $@" >> $lockfile
}

submit_job() {
    eval $(oarsub -S -n $BUILD_TAG -O "$LOGFILEDIR"/'%jobname%.OAR.%jobid%.out' -E "$LOGFILEDIR"/'%jobname%.OAR.%jobid%.err' "$@" | tee /dev/stderr | grep OAR_JOB_ID)
    if ! [ "$OAR_JOB_ID" ] ; then
        progress "oarsub failed !"
        progress "command line was: oarsub -S -n $BUILD_TAG -O "$LOGFILEDIR"/'%jobname%.OAR.%jobid%.out' -E "$LOGFILEDIR"/'%jobname%.OAR.%jobid%.err' $@"
        mv $lockfile $deadlockfile
        exit 1
    fi
    progress "submitted job $OAR_JOB_ID: $@"
    progress "stdout in $LOGFILEDIR/$BUILD_TAG.OAR.$OAR_JOB_ID.out"
    progress "stderr in $LOGFILEDIR/$BUILD_TAG.OAR.$OAR_JOB_ID.err"
    echo $OAR_JOB_ID
}

submit_jobs_and_wait_for_completion()
{
    # It's difficult because we need to make the thing synchronous.
    d=`mktemp -d /tmp/job.XXXXXXXX`
    cat > "$d/notify.sh" <<EOF
#!/bin/sh
job=\$1; shift
name=\$1; shift
tag=\$1
echo "\$@" >> $d/\$job-\$name
echo "\$@" >> $d/\$job-\$tag
touch $d/\$job-\$tag
EOF
    chmod 755 "$d/notify.sh"

    subjobs=()
    for jobfile in "$@" ; do
        x="$JOBDIR/$jobfile"
        jobid=$(submit_job --notify "exec:$d/notify.sh" "$x")
        if ! [ "$jobid" ] ; then
            progress "oarsub failed !"
            if [ "${#subjobs[@]}" -gt 0 ] ; then
                oardel "${subjobs[@]}"
            fi
            mv $lockfile $deadlockfile
            exit 1
        fi
        subjobs=("${subjobs[@]}" $jobid)
    done

    donejobs=()

    status=0
    while [ "${#donejobs[@]}" -lt "${#subjobs[@]}" ] ; do
        seen=
        for jobid in "${subjobs[@]}" ; do
            if [ -f $d/$jobid-END ] ; then
                progress "successful job $jobid: $@ [`cat $d/$jobid-END`]"
                donejobs=("${donejobs[@]}" $jobid)
                seen=1
            elif [ -f $d/$jobid-ERROR ] ; then
                progress "FAILED job $jobid: $@ [`cat $d/$jobid-ERROR`]"
                donejobs=("${donejobs[@]}" $jobid)
                status=1
                seen=1
            fi
        done
        if ! [ "$seen" ] ; then
            sleep 10
        fi
    done
    rm -rf $d
    return $status
}

# 00build is special, we want to synchronize after it completes.

if [ -f "$JOBDIR/00build" ] ; then
    if ! submit_jobs_and_wait_for_completion 00build ; then
        mv $lockfile $deadlockfile
        exit 1
    fi
fi

jobfiles=($(ls $JOBDIR | grep -v '~$' | grep -v 00build))

progress "jobs to run: ${jobfiles[@]}"

if ! submit_jobs_and_wait_for_completion "${jobfiles[@]}" ; then
    mv $lockfile $deadlockfile
    exit 1
fi

mv $lockfile $deadlockfile
exit 0
