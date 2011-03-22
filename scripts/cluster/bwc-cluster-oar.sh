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
    export MPI="$HOME/Packages/mvapich2-1.6"
    export PTHREADS=1
    export CFLAGS="-O3 -DNDEBUG"
    export CXXFLAGS="$CFLAGS"


    case "$cluster" in
        edel|graphene|griffon|parapide|parapluie)
            export MPI=$HOME/Packages/mvapich2-1.6/ 
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
        edel/30)    mpi=5x6 thr=4x2 bcode=8e0d9701;; # 0.40
        chinqchint/12)    mpi=3x4 thr=4x2 bcode=1cb0aa01;; # 1.73 mpich2-1.3.2p1
        chinqchint/30)    mpi=5x6 thr=4x2 bcode=8e0d9701;;
        parapluie/20) mpi=5x4 thr=4x3 bcode=8e0d9701;; # 0.72
        parapluie/10) mpi=5x2 thr=4x6 bcode=8e0d9701;;
        # parapluie/12) mpi=4x3 thr=5x4 bcode=8e0d9701;; # 0.92
        # parapluie/12)    mpi=2x6 thr=6x2 bcode=32d5f701;;
        # parapluie/12)    mpi=3x4 thr=4x3 bcode=32d5f701;;
        # parapluie/12)    mpi=1x12 thr=12x1 bcode=32d5f701;;
        parapluie/6)    mpi=2x3 thr=6x4 bcode=32d5f701;;
        parapluie/8)    mpi=1x8 thr=8x1 bcode=67030201;;
        # graphene/60)    mpi=5x12 thr=4x1 bcode=8e0d9701;;       # 0.36
        graphene/60)    mpi=10x6 thr=2x2 bcode=8e0d9701;;       # 0.34
        griffon/30)    mpi=5x6 thr=4x2 bcode=8e0d9701;; # 0.81
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

bwc_common_args="nullspace=left matrix=http://rsa768.rennes.grid5000.fr/cado/$matrix_base.bin wdir=$wdir mn=64 balancing=$wdir/$bfile mpi=$mpi thr=$thr"

bwc() {
    cmd="$1"
    shift
    $build_tree/linalg/bwc/bwc.pl $cmd $bwc_common_args "$@"
}

save_once() {
    rsync --port=8873 -av ${wdir}/[AVS]* rsa768.rennes.grid5000.fr::cado/
}

backlog() {
while true ; do
    date
    rsync --port=8873 -av ${wdir}/[AVS]* rsa768.rennes.grid5000.fr::cado/
    sleep $((RANDOM % 900))
done < /dev/null
}

fetch="$src/scripts/cluster/bwc-fetchfiles.pl --localdir $wdir"

$fetch --fetch $leader "md5s"

fetch="$fetch --md5db $wdir/md5s"

$fetch --fetch $leader "${nh}x${nv}/$bfile"

bwc dispatch export_cachelist=$wdir/cachelist
perl -ne "/^get-cache/ or die; chomp(\$_); my @x=split ' ', \$_; shift @x; my \$h = shift @x; @x=map { qq{${nh}x${nv}/\$_}; } @x; unshift @x, \$h; print join(' ', @x), qq{\n};" -i $wdir/cachelist

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

$fetch --dispatch-list $wdir/cachelist


backlog &

backlog_pid=$!

get_last() {
    x=$(rsync  --port=8873 --include='V*' --exclude='*' rsa768.rennes.grid5000.fr::cado/ | perl -ne '/V0-64\.(\d+)$/ && print "$1\n";' | sort -n | tail -1)
    : ${x:=0}
    echo $x
}


# Make sure to move V files before starting mksol
start=`get_last`
interval=4000

rsync -aP --port=8873 root@rsa768.rennes.grid5000.fr::cado/{V0-64.${start},C.${interval},X} /tmp/bwc.$cluster/

checkpoint_mode="keep_rolling_checkpoints=4"

if ! [ -f "$wdir/C.${interval}" ] ; then
    echo "Check vector C.${interval} missing. Keeping all V files"
    checkpoint_mode="skip_online_checks=1"
fi

    
if [ "$mksol" ] ; then
    rsync -aP --port=8873 rsa768.rennes.grid5000.fr::cado/F0-64 $wdir/

    bwc mksol interval=${interval} $checkpoint_mode start=${start} 2>&1
else
    # This could be something to consider including
    if ! [ -f "$wdir/X" ] ; then
        bwc prep interval=${interval} 2>&1
        bwc :ysplit interval=${interval} 2>&1
    fi

    bwc krylov interval=${interval} $checkpoint_mode start=${start} 2>&1
fi

kill -9 $backlog_pid
save_once
echo bye
