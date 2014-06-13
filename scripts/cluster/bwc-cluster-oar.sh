#!/usr/bin/env bash

notify() {
    LD_LIBRARY_PATH= \
    /usr/bin/curl -kni \
        https://api.grid5000.fr/sid/notifications        \
        -XPOST  \
        -d "to[]=xmpp:$USER@jabber.grid5000.fr&body=$*"
}

if [[ "$*" =~ --debug ]] ; then
    debug=1
fi

command_line="$0 $*"
ifdebug() { : ; }

if [ "$debug" ] ; then
    ifdebug() { eval "$@" ; }
fi

ifdebug date
ifdebug echo "$0 $@"

if ! [ "$CADO_NFS_AVOID_TRAMPOLINE" ] ; then
# Avoid NFS stale file handle errors !
kid=/tmp/bwc-cluster-oar.sh.$$
if [ "$0" != "bash" ] && [ "$0" != "-bash" ] ; then
if [ "$0" != "$kid" ] ; then
    ifdebug echo "script copied to $kid"
    cp "$0" "$kid"
    export CADO_NFS_SOURCE="$(cd `dirname $0`/../.. ; pwd)"
    export CADO_NFS_AVOID_TRAMPOLINE=1
    exec "$kid" "$@" < /dev/null
fi
fi
fi

if [[ "$TERM" =~ xterm ]] ; then
BS="[01m"
BE="[00m"
fi

usage_terse() {
    cat >&2 <<EOF
${BS}Usage: $0 [options]${BE}
${BS}Recognized options are:${BE}

--mksol                     trigger mksol mode (default is krylov)
-c <configfile>             read alternate config file
--slot <k>                  work with vector (slot, see --help) number k
--slot-length <n>           slot length (see --help)
<name>=<value>              set the corresponding shell variable
                            (${BS}USE THIS${BE} to set mpi and thr values)
--help                      show extended usage info
EOF
}

usage_additional() {
    cat >&2 <<EOF
${BS}Other possibly useful options:${BE}

-i <interval>               sets the interval value for checkpointing
-s <server>                 adds http://<server>/ and
                                 rsync://<server>:8873/ to the server list
-n <nodefile>               use alternate node file
--debug                     activate debug mode
--cado-nfs-source <path>    use alternate cado-nfs source path
--backlog-period <delay>    periodicity of the background backup process
--local                     work locally
--oar                       use OAR scheduler (should be auto-detected)
EOF
}

usage_extended() {
    usage_terse
    echo >&2
    usage_additional
    echo >&2
    cat >&2 <<EOF
${BS}Notes on config files:${BE}

config files are read in the following order:

Before command line parsing:

    $HOME/.cado-nfs-bwc-cluster.conf
    $CADO_NFS_SOURCE/cado-nfs-bwc-cluster.conf

    (the latter is however done after a first pass which checks for
    --cado-nfs-source arguments)

Within command line parsing:

    any file specified with -c

config files are evaluated, and thus may override the variables set by
the preceding command line switches. In particular, if a config files
redefines http_servers and/or rsync_servers, any preceding -s option is
ignored.

${BS}Notes on slots:${BE}

slots are vector numbers for krylov, but for mksol they define both a
vector number and starting iteration. In the latter case, --slot-length
defines the relevant number of iterations to the next slot. The following
formats are valid:
    0       vector 0, from iteration 0 (for krylov)
    0,0     vector 0, from iteration 0 (for mksol)
    0,1000  vector 0, from iteration 1000 (for mksol)
EOF

}


### Install basic variables which are needed in config files.

: ${CADO_NFS_SOURCE:=$(cd `dirname $0`/../.. ; pwd)}
cluster=`hostname --short | cut -d- -f1`


### Learn how to read config files, and read the default ones.

read_config_file() {
    config_file=$1
    if [ $# = 2 ] ; then
        config_file=$2
    fi
    if [ -r $config_file ] ; then
        ifdebug echo "Reading config file $config_file" >&2
        source $config_file
    elif [ $# = 2 ] ; then
        echo "$config_file: cannot read file" >&2
        exit 1
    fi
}

req_arg() {
    arg=$1
    req=$2
    n=$3
    if [ "$n" -lt "$req" ] ; then
        echo "ERROR: $arg requires $req argument(s)" >&2
        usage_terse
        exit 1
    fi
}


read_config_file $HOME/.cado-nfs-bwc-cluster.conf

### Parse command-line once so as to decide whether CADO_NFS_SOURCE needs
### to be overridden.

args=("$@")
for i in `seq 0 $(($#-1))` ; do
    arg="${args[$i]}"
    if [ "$arg" != "--cado-nfs-source" ] ; then
        continue
    fi
    req_arg $arg 1 $(($#-$i-1))
    CADO_NFS_SOURCE="${args[$(($i+1))]}"
done

read_config_file $CADO_NFS_SOURCE/cado-nfs-bwc-cluster.conf

### Parse command line

while [ "$#" -gt 0 ] ; do
    arg="$1"
    shift
    if [[ "$arg" =~ [a-z]+=.* ]] ; then
        eval "$arg"
    elif [ "$arg" = "-i" ] ; then
        req_arg $arg 1 $#
        interval="$1"
        shift
    elif [ "$arg" = "-s" ] ; then
        req_arg $arg 1 $#
        rsync_servers=(${rsync_servers[@]} rsync://${1}:8873/)
        http_servers=(${http_servers[@]} http://${1}/)
        shift
    elif [ "$arg" = "-n" ] ; then
        req_arg $arg 1 $#
        nodefile="$1"
        shift
    elif [ "$arg" = "--cado-nfs-source" ] ; then
        req_arg $arg 1 $#
        CADO_NFS_SOURCE="$1"
        shift
    elif [ "$arg" = "--backlog-period" ] ; then
        req_arg $arg 1 $#
        backlog_period="$1"
        shift
    elif [ "$arg" = "--slot" ] ; then
        req_arg $arg 1 $#
        slot="$1"
        shift
    elif [ "$arg" = "-c" ] ; then
        req_arg $arg 1 $#
        config_file="$1"
        read_config_file -r $config_file
        shift
    elif [ "$arg" = "--slot-length" ] ; then
        req_arg $arg 1 $#
        slotlength="$1"
        shift
    elif [ "$arg" = "--no-rebuild" ] ; then no_rebuild=1
    elif [ "$arg" = "--insist" ] ; then insist=1
    elif [ "$arg" = "--local" ] ; then locally=1
    elif [ "$arg" = "--mksol" ] ; then mksol=1
    elif [ "$arg" = "--debug" ] ; then debug=1
    elif [ "$arg" = "--oar" ] ; then
        nodefile=$OAR_NODEFILE
        SSH="$CADO_NFS_SOURCE/scripts/cluster/oarsh_wrap.sh"
        # TODO: other schemes may be used.
    elif [ "$arg" = "--help" ] ; then usage_extended; exit 0
    else
        echo "Unexpected argument: $arg" >&2
        exit 1
    fi
done

### Start of the real thing. First we define the functions used in the
### script.

echo "Command line: $command_line"
echo "Started: `date`"

ifdebug set -x

set -e

fail() { echo -e "$@" >&2 ; exit 1 ; }

### Start by recognizing the scheduler environment. Use it do define some
### variables. This includes SSH, although we override it only if it has
### not been defined yet.

autodetect_scheduler_and_set_dependent_vars() {
    if [ "$locally" ] ; then
        echo "Working locally"
        totalsize=1
        nodes=(`hostname --short`)
        : ${SSH=ssh}
    elif [ "$nodefile" ] ; then
        totalsize=`wc -l $nodefile | cut -d\  -f1`
        nodes=(`uniq $nodefile`)
        : ${SSH=ssh}
    elif [ "$OAR_NODEFILE" ] ; then
        echo "auto-detected OAR"
        nodefile=$OAR_NODEFILE
        totalsize=`wc -l $nodefile | cut -d\  -f1`
        nodes=(`uniq $nodefile`)
        : ${SSH="$CADO_NFS_SOURCE/scripts/cluster/oarsh_wrap.sh"}
    else
        fail "Missing nodefile, or missing --local option"
    fi
    nnodes=${#nodes[@]}
    leader=${nodes[0]}
}


### The conf files have to set some variables. Check them. Also enforce
### some policies

check_mandatory_definitions() {
    missing=""
    for v in    matrix_base server_subdir bw_m bw_n interval    \
        mpi thr                                                 \
        backlog_period slotlength slot http_servers             \
        rsync_servers build_tree ;                              \
    do
        if ! [ "$(eval "echo \$$v")" ] ; then
            missing="$missing $v"
        fi
    done
    if [ "$missing" ] ; then
        fail "The following variables MUST be set by config files"      \
            " or the command line:"     \
            "\n\t$missing"
    fi
    # http server is in reality anything which can be accessed by the curl
    # library. Unfortunately, this does not include rsync URIs, which are our
    # primary choice for the rest...

    if ! [[ $server_subdir =~ /$ ]] ; then
        server_subdir="$server_subdir/"
    fi

    # libraries may be installed on the local system anyway, so it's not
    # technically mandatory to define these.
    if [ "$GMP" ] && ! [ -d "$GMP" ] ; then fail "GMP directory $GMP not found" ; fi
    if [ "$CURL" ] && [ "$CURL" != 1 ] && ! [ -d "$CURL" ] ; then fail "CURL directory $CURL not found" ; fi
    # It's ok to not have $MPI, for local work.
    if [ "$MPI" ] && [ "$MPI" != 1 ] && ! [ -d "$MPI" ] ; then fail "MPI directory $MPI not found" ; fi
}

check_servers() {
    while true ; do
        echo "Checking rsync servers:"
        nok=0
        for rsync_server in ${rsync_servers[@]} ; do
            if rsync $rsync_server$server_subdir/ >/dev/null 2>&1 ; then
                echo -e "\t$rsync_server: alive"
                let nok+=1
            else
                echo -e "\t$rsync_server: ${BS}unreachale${BE}"
            fi
        done
        if [ "$nok" ] ; then
            break
        else
            if [ "$insist" ] ; then
                echo "All rsync servers failed. Trying again in 30 seconds" >&2
                sleep 30
            else
                 cat >&2 <<EOF
All rsync servers failed.
Arrange so that one is started on the relevant servers with the command line:
    ${BS}rsync --daemon --config=rsyncd.conf${BE}
Example config file:
    port=8873
    use chroot = no
    max connections = 8
    lock file = /home/nancy/ethome/cado/rsyncd.lock

    [cado]
          path = /srv/cado
          comment = cado
          read only = no
EOF
                exit 1
            fi
        fi
    done
}

get_possible_bcodes() {
    split="$1"
    while true ; do
        res=()
        for rsync_server in ${rsync_servers[@]} ; do
            x=$( (rsync  --include="$matrix_base.$split.????????.bin" --exclude='*' $rsync_server$server_subdir/$split/ || :) | perl -ne "m,$matrix_base\.$split\.([\da-f]+)\.bin, && print qq{\$1\n};")
            if [ "$x" ] ; then
                res=(${res[@]} $x)
            fi
        done
        if [ "$res" ] ; then break ; fi
        # If rebuilding is allowed, then we may exit from here.
        if ! [ "$no_rebuild" ] ; then break ; fi
        if [ "$insist" ] ; then
            echo "All rsync servers failed. Trying again in 30 seconds" >&2
            sleep 30
        else
            break
        fi
    done
    possible_bcodes=(`echo "${res[@]}" | xargs -n 1 echo | sort -u`)
    if [ ${#possible_bcodes[@]} -eq 0 ] ; then
        echo "No known balancing for $split" >&2
        return 1
    elif [ ${#possible_bcodes[@]} -gt 1 ] ; then
        echo "Several known possible balancings for $split: ${possible_bcodes[*]}" >&2
        return 1
    else
        bcode=${possible_bcodes[0]}
        echo "Automatically selecting $bcode as only bcode already configured for $split"
    fi
    return 0
}



# TODO: presently, when the bcode is not known, then the balancing is
# computed, which also has the effect of re-computing the dispatch. Of
# course, the thing to do is to update this file to improve its
# knowledge.
pick_cluster_configuration() {
    if [ "$bcode" ] ; then
        return
    else
        # Try to list the known bcodes, if any.
        if get_possible_bcodes ${nh}x${nv} ; then
            return
        else
            # Then it's as if we had nothing. bcode is not
            # determined, we'll have to recompute the balancing
            if [ ${#possible_bcodes[@]} -gt 1 ] ; then
                echo "Cannot decide which bcode to take" >&2
                exit 1
            fi
        fi
    fi
}

late_variables() {
    jh=`echo $mpi | cut -dx -f1`
    jv=`echo $mpi | cut -dx -f2`
    th=`echo $thr | cut -dx -f1`
    tv=`echo $thr | cut -dx -f2`
    nh=$((jh*th))
    nv=$((jv*tv))
}

display_all_variables() {
    env | grep OAR || :
    for v in CC CFLAGS CURL CXXFLAGS GMP LD_LIBRARY_PATH MPI    \
        MV2_ENABLE_AFFINITY PATH PTHREADS SSH backlog_period    \
        bcode build_tree bw_m bw_n cluster force_build_tree     \
        http_servers interval jh jv leader matrix_base mpi       \
        mpirun_extra needs_mpd nh nnodes nodes nv               \
        rsync_servers server_subdir slot slotlength src th       \
        thr totalsize tv wdir ; \
    do
        if [ "$(eval "echo \$$v")" ] ; then
            eval "echo $v=\$$v"
        fi
    done
}

mkdir_everywhere() {
    # $build_tree/linalg/bwc/bwc.pl shell command="mkdir -p $wdir" mpi=$mpi thr=$thr matrix=$http_server$server_subdir$matrix_base.bin wdir=$wdir m=$bw_m n=$bw_n
    for m in "${nodes[@]}" ; do $SSH $m mkdir -p $wdir & : ; done 
    wait
}

build() {
    # if [ -f .git ] ; then git reset --hard HEAD ; fi
    if [ -d "$build_tree" ] ; then
        cd "$build_tree"
    else
        mkdir -p "$build_tree"
        cd "$build_tree"
        cmake "$CADO_NFS_SOURCE"
    fi
    nprocs=`grep -cw ^processor /proc/cpuinfo`
    make -j $nprocs
    make -j $nprocs shell
    cd "$CADO_NFS_SOURCE"
}

start_local_rsync_server() {
    # FIXME: This check is local, and we want something global. What if
    # the leader has /tmp/rsyncd.conf, and not the other nodes ???
    # if ! [ -f /tmp/rsyncd.conf ] ; then
    # actually the rsync daemon does not care if it gets re-run. And if
    # conf file changes, it gets reread anyway.
        cat > /tmp/rsyncd.conf <<EOF
port = 8873
use chroot = no
max connections = 8
lock file = /tmp/rsyncd.lock

[tmp]
      path = $wdir
      comment = cado
      read only = yes
EOF
        for m in "${nodes[@]}" ; do rsync -e $SSH /tmp/rsyncd.conf ${m}:/tmp/rsyncd.conf & : ; done 
        wait
        for m in "${nodes[@]}" ; do $SSH -n $m rsync --daemon --config /tmp/rsyncd.conf & : ; done 
        wait
    # fi
}

http_head_somewhere() {
    if [ "$1" = "-r" ] ; then
        required=1
        shift
    fi
    what="$1"
    for http_server in ${http_servers[@]} ; do
        echo "Checking for $http_server$server_subdir$what"
        output=$(HEAD $http_server$server_subdir$what)
        if echo "$output" | grep -q '^200 OK' ; then
            echo -e "\tok"
            return 0
        fi
    done
    if [ "$required" ] ; then
        echo "$what found on none of ${http_servers[@]}" >&2
        exit 1
    fi
    return 1
}

######################################################################

### Go !

ln -sf OAR.$OAR_JOB_NAME.$OAR_JOBID.stdout OAR.$OAR_JOB_NAME.current.stdout

autodetect_scheduler_and_set_dependent_vars

notify "`date`: start $OAR_JOBID $cluster/$nnodes mpi=$mpi thr=$thr, leader=$leader"

echo "Running on $cluster, $nnodes nodes, leader is $leader"

check_mandatory_definitions

http_head_somewhere -r "$matrix_base.bin"
ncoeffs=$((`echo "$output" | awk '/^Content-Length/ { print $2 }'`/4))
http_server_with_matrix=$http_server
http_head_somewhere "$matrix_base.rw.bin"
nrows=$((`echo "$output" | awk '/^Content-Length/ { print $2 }'`/4))
http_head_somewhere "$matrix_base.cw.bin"
ncols=$((`echo "$output" | awk '/^Content-Length/ { print $2 }'`/4))
let ncoeffs-=$nrows

echo "$matrix_base.bin: $nrows rows $ncols cols $ncoeffs coeffs"

check_servers

# set_build_variables
cd "$CADO_NFS_SOURCE"
build
# set_runtime_variables
if [ "$mpi" ] && [ "$thr" ] ; then
    late_variables
    pick_cluster_configuration
    echo "bcode is $bcode"
else
    echo "Define mpi= and thr= first" >&2
    exit 1
fi
display_all_variables
mkdir_everywhere
start_local_rsync_server

fetch="$CADO_NFS_SOURCE/scripts/cluster/bwc-fetchfiles.pl --localdir $wdir -v"
for rsync_server in ${rsync_servers[@]} ; do
    fetch="$fetch -s $rsync_server$server_subdir"
done
if $fetch --on $leader "md5s" ; then
    fetch="$fetch --md5db $wdir/md5s"
else
    echo "md5 db not found on server"
fi

if [ "$force_rebuild" ] ; then
    bcode=todo
fi

if [ "$bcode" = todo ] ; then
    unset bcode
fi

if [ "$no_rebuild" ] && ! [ "$bcode" ] ; then
    echo "${BS}ERROR: Need to rebuild, although rebuilding is forbidden${BE}" >&2
    exit 1
fi

if ! [ "$bcode" ] ; then
    creating=1
    echo "${BS}Computing a balancing permutation${BE}"
    $fetch --on $leader "$matrix_base.rw.bin $matrix_base.cw.bin"
    find $wdir -name ${matrix_base}.${nh}x${nv}.????????.bin | xargs rm -vf
    $build_tree/linalg/bwc/mf_bal       \
        --shuffled-product              \
        --display-correlation           \
        --mfile $wdir/$matrix_base.bin  \
        --out $wdir/                    \
        $nh $nv
    bfile=`find $wdir -name ${matrix_base}.${nh}x${nv}.????????.bin -printf '%f\n'`
    bcode=$bfile
    bcode=${bcode##${matrix_base}.${nh}x${nv}.}
    bcode=${bcode%%.bin}
    for m in "${nodes[@]}" ; do $SSH $m find $wdir -name ${matrix_base}.${nh}x${nv}.${bcode}\*bucketT.bin \| xargs rm -vf & : ; done
    wait
    $fetch --put --on $leader "${nh}x${nv}/$bfile"
fi
bfile="${matrix_base}.${nh}x${nv}.${bcode}.bin"
echo "need to get $bfile"

bwc_common_args="nullspace=left matrix=$http_server_with_matrix$server_subdir$matrix_base.bin wdir=$wdir m=$bw_m n=$bw_n balancing=$bfile mpi=$mpi thr=$thr"

bwc() {
    cmd="$1"
    shift
    $build_tree/linalg/bwc/bwc.pl $cmd $bwc_common_args "$@"
}

sed_on_cachelist() {
    perl -ne "/^get-cache/ or die; chomp(\$_); my @x=split ' ', \$_; shift @x; my \$h = shift @x; @x=map { qq{${nh}x${nv}/\$_}; } @x; unshift @x, \$h; print join(' ', @x), qq{\n};" -i $wdir/cachelist
}

if [ "$creating" ] ; then
    echo "Running original dispatch"
    bwc dispatch sanity_check_vector=H1 save_submatrices=1 sequential_cache_build=1
    bwc dispatch export_cachelist=$wdir/cachelist
    sed_on_cachelist
    $fetch --put --dispatch-list $wdir/cachelist
else
    echo "Fetching cache files from remote"
    $fetch --on $leader "${nh}x${nv}/$bfile"
    bwc dispatch export_cachelist=$wdir/cachelist
    sed_on_cachelist
    $fetch --dispatch-list $wdir/cachelist
fi

# time 
expected_size_V=$((nrows * (64/8)))
expected_size_S=$((nrows * ($bw_n/8)))
expected_size_A_onestep=$((bw_m * 64 / 8))

save_once() {
    # For save_once, it's normal to have failures !
    set +x
    set +e
    # we're still alive: propose ourselves as exporters of the matrix
    # files...
    matrixfiles=(`cd $wdir ; ls ${matrix_base}.*.bin`)
    for f in "${matrixfiles[@]}" ; do
        ln -sf "rsync://`hostname`:8873/tmp" $HOME/.cache/$f
    done
    localfiles=(`cd $wdir ; ls [AVS]*`)
    process=()
    for f in "${localfiles[@]}" ; do
        if [ -f "$wdir/done.$f" ] ; then
            # echo "$f already saved"
            continue
        fi

        if [[ $f =~ ^A[0-9]+-[0-9]+\.([0-9]+)-([0-9]+) ]] ; then
            i0=${BASH_REMATCH[1]}
            i1=${BASH_REMATCH[2]}
            exp=$((expected_size_A_onestep*(i1-i0)))
        elif [[ $f =~ ^V ]] ; then
            exp=$expected_size_V
        elif [[ $f =~ ^S ]] ; then
            exp=$expected_size_S
        else
            echo "What is $f ???" >&2
            continue
        fi
        bytes=`du --bytes --apparent-size $wdir/$f | awk '// { print $1 }'`
        if [ "$bytes" != "$exp" ] ; then
            echo "$f has $bytes != $exp bytes, not saving"
            continue
        fi
        # echo "will save $f"
        process=("${process[@]}" "$slot/$f")
    done
    if [ ${#process[@]} -gt 0 ] ; then
        $fetch --put --on $leader "${process[*]}"
    fi
}

backlog() {
    # For save_once, it's normal to have failures !
    set +x
    set +e
    while true ; do
        date
        save_once
        sleep $((RANDOM % $backlog_period))
    done < /dev/null
}

get_last() {
    myslot="$1"
    base="$1/" ; shift
    if [ "$mksol" ] && [[ "$myslot" =~ ^([0-9]+),([0-9]+)$ ]] ; then
        n1=${BASH_REMATCH[1]}
        n2=$((n1+64))
        last=${BASH_REMATCH[2]}
    elif ! [ "$mksol" ] && [[ "$myslot" =~ ^([0-9]+)$ ]] ; then
        n1=${BASH_REMATCH[1]}
        n2=$((n1+64))
        last=0
    else
        echo "$1: $!" >&2 ; exit 1
    fi
    for rsync_server in ${rsync_servers[@]} ; do
        x=$(rsync --include='V*' --exclude='*' $rsync_server$server_subdir$base | perl -ne 'm{\.(\d+)$} && print "$1\n";' | sort -n | tail -1)
        : ${x:=0}
        if [ "$x" -gt "$last" ] ; then
            last=$x
        fi
    done
    echo $last
}

# Caveat: we have to make sure to move V files before starting mksol.
# This script has no provision for doing this automatically. If it isn't
# done, we'll do bogus computations.
start=`get_last $slot`

# This could be something to consider including
if ! $fetch --on $leader "X $slot/V0-64.${start}" ; then
    if [ "$start" != 0 ] || [ "$slot" != 0 ] ; then
        echo "not finding X where start > 0 looks like a bug" >&2
        exit 1
    fi
    if ! $fetch --on $leader "X Y.0" ; then
        echo "init vectors not found, re-running prep"
        bwc prep interval=${interval} 2>&1
        $fetch --put --on $leader "X Y.0"
    fi
    bwc :ysplit interval=${interval} 2>&1
    $fetch --put --on $leader "$slot/V0-64.${start}"
fi

checkpoint_mode="keep_rolling_checkpoints=4 keep_checkpoints_younger_than=$((backlog_period*10))"

# Either we prepare re-computation of the check vector, or we defer its
# computation.
if ! $fetch --on $leader "C.$interval" ; then
    echo "Check vector C.${interval} missing. Re-creating"
    bwc secure interval=${interval} 2>&1
    $fetch --put --on $leader "C.$interval"
fi
if ! [ -f "$wdir/C.$interval" ] ; then
    echo "Check vector C.${interval} missing. Keeping all V files"
    checkpoint_mode="skip_online_checks=1"
fi

# Wakes up bwc.pl with SIGHUP if jobs should be restarted
watch_server_progress() {
    new_best=$start
    while sleep 300 ; do
        my_best=$(cd $wdir ; ls V* | perl -ne 'm{\.(\d+)$} && print "$1\n";' | sort -n | tail -1)
        tmp_best=`get_last $slot`
        if [ "$tmp_best" -gt "$new_best" ] ; then
            new_best=$tmp_best
            echo "`date` *** Central data reaches iteration $tmp_best (here: $my_best)"
        fi
        if [ $my_best -lt $((new_best - 2*interval)) ] ; then
            echo "`date` *** Central data has already reached iteration $new_best, restarting ***"
            $fetch --on $leader "X $slot/V0-64.${new_best}"
            pkill -USR1 -f bwc.pl
            pkill -USR1 -f mpiexec
        fi
    done
}

# Now that we fork, we avoid -e, since it does not dismiss children !
set +e
backlog &
backlog_pid=$!

watch_server_progress &
watch_server_progress_pid=$!

 
while true ; do
    start=`get_last $slot`
    if [ "$mksol" ] ; then
        if [[ "$slot" =~ ^([0-9]+),([0-9]+)$ ]] ; then
            n1=${BASH_REMATCH[1]}
            n2=$((n1+64))
            slotbase=${BASH_REMATCH[2]}
            slotend=$((slotbase + slotlength))
        else
            echo "$slot: bad slot" >&2 ; exit 1
        fi
        # This is fairly stupid. F contains too much information in reality
        $fetch --on $leader "F${n1}-${n2}"
        bwc mksol interval=${interval} $checkpoint_mode start=${start} end=$slotend 2>&1
    else
        bwc krylov interval=${interval} $checkpoint_mode start=${start} 2>&1
    fi
    rc=$?
    if [ $rc -ge 128 ] && ( [ $((rc & 127)) = 1 ] ||  [ $((rc & 127)) = 10 ] ) ; then
        echo "*** bwc process aborted with SIGHUP/SIGUSR1, resuming ***"
        pkill -9 -f mpiexec
    else
        break
    fi
done

kill -9 $backlog_pid
kill -9 $watch_server_progress_pid

save_once
echo bye
