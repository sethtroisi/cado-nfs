#!/bin/bash

matrix_base=snfs247.small

echo "$0 $@"

date

# Avoid NFS stale file handle errors !
kid=/tmp/bwc-cluster-oar.sh.$$
if [ "$0" != "bash" ] && [ "$0" != "-bash" ] ; then
if [ "$0" != "$kid" ] ; then
    echo "script copied to $kid"
    cp "$0" "$kid"
    exec "$kid" "$@" < /dev/null
fi
fi

cluster=`hostname --short | cut -d- -f1`
totalsize=`wc -l $OAR_NODEFILE | cut -d\  -f1`
nodes=(`uniq $OAR_NODEFILE`)
nnodes=${#nodes[@]}
leader=${nodes[0]}

set_build_variables() {
    src=$HOME/cado

    export PATH=$HOME/Packages/cmake-2.6.4/bin/:$PATH
    export build_tree=./prod/$cluster

    # defaults
    export GMP="$HOME/Packages/gmp-5.0.1"
    export CURL="$HOME/Packages/curl-7.21.1"
    export MPI="$HOME/Packages/mvapich2-1.6rc2"
    export PTHREADS=1
    export CFLAGS="-O3 -DNDEBUG"
    export CXXFLAGS="$CFLAGS"


    case "$cluster" in
        edel|graphene|griffon|parapide|parapluie)
            export MPI=$HOME/Packages/mvapich2-1.6rc2/ 
            needs_mpd=1
            if [ -d "$HOME/Packages/ib/lib" ] ; then
                export LD_LIBRARY_PATH=$HOME/Packages/ib/lib:$LD_LIBRARY_PATH
            fi
            ;;
        *)
            echo "Cluster $cluster not configured" >&2
            exit 1
            ;;
    esac

    if ! [ -d "$GMP" ] ; then echo "GMP directory $GMP not found" ; fi
    if ! [ -d "$CURL" ] ; then echo "CURL directory $CURL not found" ; fi
    if ! [ -d "$MPI" ] ; then echo "MPI directory $MPI not found" ; fi
}

set_runtime_variables() {
    export SSH="$src/scripts/cluster/oarsh_wrap.sh"
    if [[ $MPI =~ mvapich2 ]] ; then export MV2_ENABLE_AFFINITY=0 ; fi
    if [[ $MPI =~ openmpi ]] ; then mpirun_extra="--mca plm_rsh_agent $SSH" ; fi
    export wdir=/tmp/bwc.$cluster
}

pick_cluster_configuration() {
    case "$cluster/$nnodes" in
        edel/30)    mpi=5x6 thr=4x2 bcode=8e0d9701;;
        parapluie/20) mpi=5x4 thr=4x3 bcode=8e0d9701;;
        # parapluie/12)    mpi=2x6 thr=6x2 bcode=32d5f701;;
        # parapluie/12)    mpi=3x4 thr=4x3 bcode=32d5f701;;
        parapluie/12)    mpi=1x12 thr=12x1 bcode=32d5f701;;
        parapluie/6)    mpi=2x3 thr=6x4 bcode=32d5f701;;
        parapluie/8)    mpi=1x8 thr=8x1 bcode=67030201;;
        *)
            echo "Configuration $cluster/$nnodes not configured" >&2
            exit 1
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
    env | grep OAR
    for v in wdir leader cluster totalsize nodes nnodes src PATH force_build_tree GMP CURL MPI PTHREADS CFLAGS CXXFLAGS MPI needs_mpd LD_LIBRARY_PATH SSH MV2_ENABLE_AFFINITY mpirun_extra mpi thr bcode jh jv th tv nh nv build_tree CFLAGS CC ; do
        eval "echo $v=\$$v"
    done
}

mkdir_everywhere() {
    for m in "${nodes[@]}" ; do $SSH $m mkdir -p $wdir & : ; done 
    wait
}

build() {
    if [ -f .git ] ; then
        git reset --hard HEAD
    fi
    if [ -d "$build_tree" ] ; then
        cd "$build_tree"
    else
        mkdir -p "$build_tree"
        cd "$build_tree"
        cmake "$src"
    fi
    make -j 32
    cd "$src"
}

set_build_variables
cd "$src"
build
set_runtime_variables
pick_cluster_configuration
late_variables
display_all_variables
mkdir_everywhere

bfile="${matrix_base}.${nh}x${nv}.${bcode}.bin"
echo "need to get $bfile"

fetch="$src/scripts/cluster/bwc-fetchfiles.pl --localdir $wdir"

$fetch --fetch $leader "${nh}x${nv}/$bfile"

bwc_common_args="nullspace=left matrix=$matrix_base.bin wdir=$wdir mn=64 balancing=$wdir/$bfile mpi=$mpi thr=$thr"

bwc() {
    cmd="$1"
    shift
    $build_tree/linalg/bwc/bwc.pl $cmd $bwc_common_args "$@"
}

bwc u64_dispatch export_cachelist=$wdir/cachelist
perl -ne "/^get-cache/ or die; chomp(\$_); my @x=split ' ', \$_; shift @x; my \$h = shift @x; @x=map { qq{${nh}x${nv}/\$_}; } @x; unshift @x, \$h; print join(' ', @x), qq{\n};" -i $wdir/cachelist

if ! [ -f /tmp/rsyncd.conf ] ; then
    cat > /tmp/rsyncd.conf <<EOF
port=8873
use chroot = no
max connections = 8
lock file = /tmp/rsyncd.lock

[tmp]
      path = $wdir
      comment = cado
      read only = yes
EOF
    $MPI/bin/mpiexec -n $nnodes $build_tree/linalg/bwc/bcast-file /tmp/rsyncd.conf
    for m in "${nodes[@]}" ; do $SSH -n $m rsync --daemon --config /tmp/rsyncd.conf & : ; done 
    wait
fi

$fetch --dispatch-list $wdir/cachelist

backlog() {
while true ; do
    date
    rsync --port=8873 -av ${wdir}/[AVS]* rsa768.rennes.grid5000.fr::cado/${bcode}/
    sleep $((RANDOM % 900))
done < /dev/null
}

backlog &

backlog_pid=$!

get_last() {
    tmp_bcode="$1"
    x=$(rsync  --port=8873 --include='V*' --exclude='*' root@rsa768.rennes.grid5000.fr::cado/${tmp_bcode}/ | perl -ne '/V0-64\.(\d+)\./ && print "$1\n";' | sort -n | tail -1)
    : ${x:=0}
    echo $x
}


if [ "$mksol" ] ; then
    rsync -aP --port=8873 rsa768.rennes.grid5000.fr::cado/F0-64 $wdir/
fi

# Make sure to move V files before starting mksol
start=`get_last $bcode`
interval=4000

rsync -aP --port=8873 root@rsa768.rennes.grid5000.fr::cado/${bcode}/{V0-64.${start}.${bcode},C.${interval}.${bcode},X.${bcode}} /tmp/bwc.$cluster/

# This could be something to consider including
# +bwc u64n_prep interval=${interval} 2>&1
# +bwc :ysplit interval=${interval} 2>&1

if [ "$mksol" ] ; then
    bwc u64_mksol interval=${interval} keep_rolling_checkpoints=4 start=${start} 2>&1
else
    bwc u64_krylov interval=${interval} keep_rolling_checkpoints=4 start=${start} 2>&1
fi

echo bye
kill -9 $backlog_pid
