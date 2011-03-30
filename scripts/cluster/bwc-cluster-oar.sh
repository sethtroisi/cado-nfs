#!/bin/bash

date
echo "$0 $@"

# Avoid NFS stale file handle errors !
kid=/tmp/bwc-cluster-oar.sh.$$
if [ "$0" != "bash" ] && [ "$0" != "-bash" ] ; then
if [ "$0" != "$kid" ] ; then
    echo "script copied to $kid"
    cp "$0" "$kid"
    mat="$1"
    shift
    exec "$kid" "$mat" "--cado-source" "$(cd `dirname $0`/../.. ; pwd)" "$@" < /dev/null
fi
fi

bw_m=64
bw_n=64
interval=100
backlog_period=120
# http server is in reality anything which can be accessed by the curl
# library. Unfortunately, this does not include rsync URIs, which are our
# primary choice for the rest...
# rsync_server=rsync://rsa768.rennes.grid5000.fr:8873/
# http_server=http://rsa768.rennes.grid5000.fr/
rsync_servers=()
http_servers=()

matrix_base="$1"
shift
if [ "$matrix_base" = "" ] ; then
    echo "Usage: $0 [matrix base name]" >&2
    exit 1
fi

# This is required to end with a /
server_subdir=cado/$matrix_base/

SSH=ssh

export SSH

req_arg() {
    arg=$1
    req=$2
    n=$3
    if [ "$n" -lt "$req" ] ; then
        echo "$arg requires $req argument(s)" >&2
        exit 1
    fi
}

while [ "$#" -gt 0 ] ; do
    arg="$1"
    shift
    if [ "$arg" = "-i" ] ; then
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
    elif [ "$arg" = "--cado-source" ] ; then
        req_arg $arg 1 $#
        cado_source="$1"
        shift
    elif [ "$arg" = "--backlog-period" ] ; then
        req_arg $arg 1 $#
        backlog_period="$1"
        shift
    elif [ "$arg" = "--local" ] ; then locally=1
    elif [ "$arg" = "--mksol" ] ; then mksol=1
    elif [ "$arg" = "--oar" ] ; then
        nodefile=$OAR_NODEFILE
        SSH="$cado_source/scripts/cluster/oarsh_wrap.sh"
        # TODO: other schemes may be used.
    else
        echo "Unexpected argument: $arg" >&2
        exit 1
    fi
done

set -e
set -x

cluster=`hostname --short | cut -d- -f1`

if [ "$locally" ] ; then
    echo "Working locally"
    totalsize=1
    nodes=(`hostname --short`)
elif [ "$nodefile" ] ; then
    totalsize=`wc -l $nodefile | cut -d\  -f1`
    nodes=(`uniq $nodefile`)
elif [ "$OAR_NODEFILE" ] ; then
    echo "auto-detected OAR"
    nodefile=$OAR_NODEFILE
    totalsize=`wc -l $nodefile | cut -d\  -f1`
    nodes=(`uniq $nodefile`)
    SSH="$cado_source/scripts/cluster/oarsh_wrap.sh"
else
    echo "Missing nodefile, or --local option" >&2
    exit 1
fi

nnodes=${#nodes[@]}
leader=${nodes[0]}
    
fail() { echo "$@" >&2 ; exit 1 ; }

set_build_variables() {
    export PATH=$HOME/Packages/cmake-2.6.4/bin/:$PATH
    export build_tree=./prod/$cluster

    # defaults
    export GMP="$HOME/Packages/gmp-5.0.1"
    export CURL="$HOME/Packages/curl-7.21.1"
    export PTHREADS=1
    export CFLAGS="-O3 -DNDEBUG"
    export CXXFLAGS="$CFLAGS"


    case "$cluster" in
        edel|graphene|griffon|parapide|parapluie|genepi)
            export MPI=$HOME/Packages/mvapich2-1.6/ 
            if [ -d "$HOME/Packages/ib/lib" ] ; then
                export LD_LIBRARY_PATH=$HOME/Packages/ib/lib:$LD_LIBRARY_PATH
            fi
            ;;
        chinqchint)
            export MPI=$HOME/Packages/mpich2-1.3.2p1/
            needs_mpd=1
            ;;
        truffe)
            unset GMP
            unset CURL
            unset MPI
            ;;
        tiramisu)
            MPI=1
            GMP=/users/caramel/logiciels/gmp-5.0.1-patched/core2
            CURL=1
            ;;
        *)
            if ! [ "$locally" ] ; then
                echo "Cluster $cluster not configured" >&2
                exit 1
            fi
            ;;
    esac

    # libraries may be installed on the local system anyway
    if [ "$GMP" ] && ! [ -d "$GMP" ] ; then fail "GMP directory $GMP not found" ; fi
    if [ "$CURL" ] && [ "$CURL" != 1 ] && ! [ -d "$CURL" ] ; then fail "CURL directory $CURL not found" ; fi
    # It's ok to not have $MPI, for local work.
    if [ "$MPI" ] && [ "$MPI" != 1 ] && ! [ -d "$MPI" ] ; then fail "MPI directory $MPI not found" ; fi
}

set_runtime_variables() {
    if [[ $MPI =~ mvapich2 ]] ; then export MV2_ENABLE_AFFINITY=0 ; fi
    if [[ $MPI =~ openmpi ]] ; then mpirun_extra="--mca plm_rsh_agent $SSH" ; fi
    export wdir=/tmp/bwc.$cluster
}

# TODO: presently, when the bcode is not known, then the balancing is
# computed, which also has the effect of re-computing the dispatch. Of
# course, the thing to do is to update this file to improve its
# knowledge.
pick_cluster_configuration() {
    if [ "$mpi" ] && [ "$thr" ] && [ "$bcode" ] ; then
        return
    fi
    case "$matrix_base/$cluster/$nnodes" in
        example/truffe/1) mpi=1x1 thr=2x2 bcode=78713a01;;
        # c72/*/1) mpi=1x1 thr=1x1 bcode=26082101;;
        c72/*/1) mpi=1x1 thr=2x2 bcode=5aae0b01;;
        snfs247.small/edel/30)    mpi=5x6 thr=4x2 bcode=ca4a7001;; # old 0.40
        snfs247.small/genepi/12)    mpi=3x4 thr=4x2 bcode=0a426301;;
        snfs247.small/chinqchint/12)    mpi=3x4 thr=4x2 bcode=0a426301;; # new 1.80 old 1.73 mpich2-1.3.2p1
        womack/chinqchint/12)    mpi=3x4 thr=4x2 bcode=todo;;
        snfs247.small/chinqchint/30)    mpi=5x6 thr=4x2 bcode=todo;;
        snfs247.small/parapluie/20) mpi=5x4 thr=4x3 bcode=ca4a7001;; # new 0.72 old:0.72
        snfs247.small/parapluie/10) mpi=5x2 thr=4x6 bcode=todo;;
        # parapluie/12) mpi=4x3 thr=5x4 bcode=todo;; # old 0.92
        # parapluie/12)    mpi=2x6 thr=6x2 bcode=todo;;
        # parapluie/12)    mpi=3x4 thr=4x3 bcode=todo;;
        # parapluie/12)    mpi=1x12 thr=12x1 bcode=todo;;
        snfs247.small/parapluie/6)    mpi=2x3 thr=6x4 bcode=todo;;
        snfs247.small/parapluie/8)    mpi=1x8 thr=8x1 bcode=todo;;
        # graphene/60)    mpi=5x12 thr=4x1 bcode=todo;;       # old 0.36
        snfs247.small/graphene/60)    mpi=10x6 thr=2x2 bcode=ca4a7001;; # old 0.34
        womack/graphene/30)    mpi=5x6 thr=2x2 bcode=65bc0c01;;
        snfs247.small/griffon/30)    mpi=5x6 thr=4x2 bcode=ca4a7001;; # new 0.78 old 0.81
        snfs247.small/griffon/32)    mpi=4x8 thr=4x2 bcode=cf486301;; # new 0.72
        snfs247.small/graphene/16)    mpi=4x4 thr=2x2 bcode=50e16501;; # new 0.90
        snfs247.small/griffon/60)    mpi=10x6 thr=2x2 bcode=ca4a7001;; # new 0.56
        *)
            if ! ( [ "$mpi" ] && [ "$thr" ] ) ; then
                echo "Configuration $matrix_base/$cluster/$nnodes not configured" >&2
                exit 1
            fi
            ;;
    esac
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
    for v in wdir leader cluster totalsize nodes nnodes src PATH force_build_tree GMP CURL MPI PTHREADS CFLAGS CXXFLAGS MPI needs_mpd LD_LIBRARY_PATH SSH MV2_ENABLE_AFFINITY mpirun_extra mpi thr bcode jh jv th tv nh nv build_tree CFLAGS CC ; do
        eval "echo $v=\$$v"
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
        cmake "$cado_source"
    fi
    make -j 32
    make -j 32 shell
    cd "$cado_source"
}

start_local_rsync_server() {
    # FIXME: This check is local, and we want something global. What if
    # the leader has /tmp/rsyncd.conf, and not the other nodes ???
    if ! [ -f /tmp/rsyncd.conf ] ; then
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
    fi
}

http_get_somewhere() {
    if [ "$1" = "-r" ] ; then
        reuired=1
        shift
    fi
    what="$1"
    for http_server in ${http_servers[@]} ; do
        echo "Checking for $http_server$server_subdir$what"
        output=$(HEAD $http_server$server_subdir$what)
        if echo "$output" | grep -q '^200 OK' ; then
            echo "ok"
            return 0
        fi
    done
    if [ "$required" ] ; then
        echo "$what found on none of ${http_servers[@]}" >&2
        exit 1
    fi
    return 1
}

http_get_somewhere -r "$matrix_base.bin"
ncoeffs=$((`echo "$output" | awk '/^Content-Length/ { print $2 }'`/4))
http_get_somewhere "$matrix_base.rw.bin"
nrows=$((`echo "$output" | awk '/^Content-Length/ { print $2 }'`/4))
http_get_somewhere "$matrix_base.cw.bin"
ncols=$((`echo "$output" | awk '/^Content-Length/ { print $2 }'`/4))
let ncoeffs-=$nrows

echo "$matrix_base.bin: $nrows rows $ncols cols $ncoeffs coeffs"

set_build_variables
cd "$cado_source"
build
set_runtime_variables
pick_cluster_configuration
late_variables
display_all_variables
mkdir_everywhere
start_local_rsync_server

fetch="$cado_source/scripts/cluster/bwc-fetchfiles.pl --localdir $wdir"
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

if ! [ "$bcode" ] ; then
    creating=1
    echo "Computing a balancing permutation"
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

bwc_common_args="nullspace=left matrix=$http_server$server_subdir$matrix_base.bin wdir=$wdir m=$bw_m n=$bw_n balancing=$bfile mpi=$mpi thr=$thr"

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
    bwc dispatch sanity_check_vector=H1
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
        process=("${process[@]}" "$f")
    done
    if [ ${#process} -gt 0 ] ; then
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
    x=$(rsync  --include='V*' --exclude='*' $rsync_server$server_subdir | perl -ne '/V0-64\.(\d+)$/ && print "$1\n";' | sort -n | tail -1)
    : ${x:=0}
    echo $x
}

# Caveat: we have to make sure to move V files before starting mksol.
# This script has no provision for doing this automatically. If it isn't
# done, we'll do bogus computations.
start=`get_last`

# This could be something to consider including
if ! $fetch --on $leader "X V0-64.${start}" ; then
    if [ "$start" != 0 ] ; then
        echo "not finding X where start > 0 looks like a bug" >&2
        exit 1
    fi
    if ! $fetch --on $leader "X Y.0" ; then
        echo "init vectors not found, re-running prep"
        bwc prep interval=${interval} 2>&1
        $fetch --put --on $leader "X Y.0"
    fi
    bwc :ysplit interval=${interval} 2>&1
    $fetch --put --on $leader "V0-64.${start}"
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

# Now that we fork, we avoid -e, since it does not dismiss children !
set +e
backlog &
backlog_pid=$!
    
if [ "$mksol" ] ; then
    $fetch --on $leader "F0-64"
    bwc mksol interval=${interval} $checkpoint_mode start=${start} 2>&1
else
    bwc krylov interval=${interval} $checkpoint_mode start=${start} 2>&1
fi

kill -9 $backlog_pid
save_once
echo bye
