#!/usr/bin/env bash

set -e
# set -x

# This script is intended *for testing only*. It's used for, e.g.,
# coverage tests. This even goes with dumping all intermediary data to
# magma-parseable format. So don't try to use this script for real
# examples. Real examples require the user to craft his own bwc.pl
# command line with care (command-lines as created by this tool might be
# a source of inspiration, though).

while [ $# -gt 0 ] ; do
    a="$1"
    shift
    if [ "$a" = "--" ] ; then
        break
    else
        eval "$a"
    fi
done

pass_bwcpl_args=("$@")

# various configuration variables. environment can be used to override them
: ${scriptpath=$0}
: ${m=8}
: ${n=4}
: ${prime=4148386731260605647525186547488842396461625774241327567978137}
: ${mpi=1x1}
: ${thr=2x2}
# Set the "matrix" variable in order to work on a real matrix.
: ${matrix=}
# Set the "bindir" variable to use pre-built binaries (must point to the
# directories with the bwc binaries)
: ${bindir=}

# This is related to the random matrix which gets automatically created
# if no matrix was given on the cmdline.
: ${random_matrix_size=1000}
# there is no random_matrix_coeffs_per_row ; see random_matrix.c
: ${random_matrix_maxcoeff=10}
: ${random_matrix_minkernel=10}
: ${mats=$HOME/Local/mats}
: ${matsfallback=/local/rsa768/mats}
: ${pre_wipe=}
: ${seed=$RANDOM}

wordsize=64
# XXX note that $wdir is wiped out by this script !
: ${wdir=/tmp/bwcp}
: ${shuffle=1}

# By default we don't enable magma. Just say magma=magma (or
# magma=/path/to/magma) to get (very expensive) magma checks.
: ${magma=}

# This has to be provided as auxiliary data for the GF(p) case (or we'll
# do without any rhs whatsoever).
: ${rhsfile=}
: ${nullspace=right}
: ${interval=50}
: ${mm_impl=basicp}

if ! type -p seq >/dev/null ; then
    seq() {
        first="$1"
        shift
        if [ "$2" ] ; then
            incr="$1"
            shift
        else
            incr=1
        fi
        last="$1"
        while [ "$first" -le "$last" ] ; do
            echo "$first"
            let first=first+incr
        done
    }
fi

usage() {
    cat >&2 <<-EOF
    Usage: $scriptpath [param1=value param2=value2 ...]
    This runs a complete block Wiedemann algorithm, given the parameters
    specified by various environment variables (the command line can also
    be used to set those variables).
EOF
    exit 1
}

argument_checking() {
    case "$nullspace" in
        left|right) ;;
        *) echo "\$nullspace must be left or right" >&2; usage;;
    esac
    if [ "$prime" != 2 ] ; then
        if [ "$rhsfile" ] && ! [ "$matrix" ] ; then
            echo "Please specify \$rhsfile only with \$matrix" >&2
            usage
        fi
        if [ "$matrix" ] && [ "$nrhs" ] ; then
            echo "Please specify \$nrhs only without \$matrix (so that a random problem gets generated)" >&2
            usage
        fi
    fi
}

derived_variables() {
    if ! [ -d $mats ] ; then mats=$matsfallback; fi
    if [ "$prime" = 2 ] ; then
        splitwidth=64
    else
        splitwidth=1
    fi
}

prepare_wdir() {
    if [ -d $wdir ] ; then
        if [ "$pre_wipe" ] ; then
            rm -rf $wdir 2>/dev/null
        else
            echo "Won't wipe $wdir unless \$pre_wipe is set" >&2
            exit 1
        fi
    fi
    mkdir $wdir
}


create_test_matrix_if_needed() {
    if [ "$matrix" ] ; then
        if ! [ -e "$matrix" ] ; then
            echo "matrix=$matrix inaccessible" >&2
            usage
        fi
        # Get absolute path.
        matrix=$(readlink /proc/self/fd/99 99< $matrix)
        if [ -e "$rhsfile" ] ; then
            rhsfile=$(readlink /proc/self/fd/99 99< $rhsfile)
        fi
        return
    fi

    # It's better to look for a kernel which is not trivial. Thus
    # specifying --kright for random generation is a good move prior to
    # running this script for nullspace=right
    kside="--k$nullspace"

    # We only care the binary matrix, really. Nevertheless, the random
    # matrix is created as text, and later transformed to binary format.
    # We also create the auxiliary .rw and .cw files.

    # We're not setting the density, as there is an automatic setting in
    # random_matrix for that (not mandatory though, since
    # nrows,ncols,density, may be specified in full).

    rmargs=(${random_matrix_size} -s $seed)
    if [ "$prime" = 2 ] ; then
        basename=$mats/t${random_matrix_size}
        matrix="$basename.matrix.bin"
        rmargs=("${rmargs[@]}"  --k$nullspace ${random_matrix_minkernel})
    elif ! [ "$nrhs" ] ; then
        basename=$mats/t${random_matrix_size}p
        matrix="$basename.matrix.bin"
        rmargs=("${rmargs[@]}" --k$nullspace ${random_matrix_minkernel})
        rmargs=("${rmargs[@]}" -c ${random_matrix_maxcoeff})
    else
        basename=$mats/t${random_matrix_size}p+${nrhs}
        density=`echo "l($random_matrix_size)^2/2" | bc -l | cut -d. -f1`
        if [ "$density" -lt 12 ] ; then density=12; fi
        rmargs=("${rmargs[@]}" -d $density)
        matrix="$basename.matrix.bin"
        rhsfile="$basename.rhs.bin"
        rmargs=("${rmargs[@]}" -c ${random_matrix_maxcoeff})
        rmargs=("${rmargs[@]}" rhs="$nrhs,$prime,$rhsfile")
    fi
    rmargs=("${rmargs[@]}" --freq --binary --output "$matrix")
    ${bindir}/random_matrix "${rmargs[@]}"
    rwfile=${matrix%%bin}rw.bin
    cwfile=${matrix%%bin}cw.bin
    ncols=$((`wc -c < $cwfile` / 4))
    nrows=$((`wc -c < $rwfile` / 4))
}

# This is only useful in the situation where an input matrix has been
# provided by the user (which may be ither in text or binary format). In
# this case, we must make sure that we have the .cw and .rw files too.
create_auxiliary_weight_files() {
    if [ "$prime" != 2 ] ; then withcoeffs=--withcoeffs ; fi
    case "$matrix" in
        *.txt)
            if [ "$rhsfile" ] ; then
                # It's really a hassle to keep the conversion code (which
                # existed until dad7019)
                echo "Please supply $rhs as a *binary* file, please\n" >&2
                exit 1
            fi
            matrix_txt="$matrix"
            matrix=${matrix%%txt}bin
            rwfile=${matrix%%bin}rw.bin
            cwfile=${matrix%%bin}cw.bin
            if [ "$matrix" -nt "$matrix_txt" ] && [ "$rwfile" -nt "$matrix_txt" ] && [ "$cwfile" -nt "$matrix_txt" ] ; then
                echo "Taking existing $mfile, $rwfile, $cwfile as accompanying $matrix_txt"
            else
                echo "Creating files $matrix, $rwfile, $cwfile from $matrix_txt"
                $bindir/mf_scan  --ascii-in $withcoeffs --mfile $matrix_txt  --freq --binary-out --ofile $matrix
            fi
            ;;
        *.bin)
            rwfile=${matrix%%bin}rw.bin
            cwfile=${matrix%%bin}cw.bin
            if [ "$rwfile" -nt "$matrix" ] && [ "$cwfile" -nt "$matrix" ] ; then
                echo "Taking existing $rwfile, $cwfile as accompanying $matrix"
            else
                $bindir/mf_scan  --binary-in $withcoeffs --mfile $matrix --freq
            fi
            ;;
    esac
    ncols=$((`wc -c < $cwfile` / 4))
    nrows=$((`wc -c < $rwfile` / 4))
}

prepare_common_arguments() {
    # This sets the common arguments for all bwc binaries
    read -s -r -d '' common <<-EOF
        matrix=$matrix
        mpi=$mpi
        thr=$thr
        m=$m
        n=$n
        wdir=$wdir
        prime=$prime
        mm_impl=$mm_impl
        nullspace=$nullspace

EOF
}

argument_checking
derived_variables
prepare_wdir

if ! [ "$matrix" ] ; then
    # This also sets rwfile cwfile nrows ncols
    create_test_matrix_if_needed
else
    # This also sets rwfile cwfile nrows ncols
    create_auxiliary_weight_files
fi

prepare_common_arguments

if [ "$rhsfile" ] ; then
    set $common rhs=$rhsfile nrhs=$nrhs
else
    set $common
fi

if [ "$magma" ] ; then
    echo "### Enabling magma checking ###"
    interval=1
    set "$@" save_submatrices=1 interval=1
fi

$bindir/bwc.pl :complete "$@" "${pass_bwcpl_args[@]}"

if ! [ "$magma" ] ; then
    exit 0
fi

set +x

mdir=$wdir

if ! [[ "$mpi,$thr" =~ ^([0-9]*)x([0-9]*),([0-9]*)x([0-9]*)$ ]] ; then
    echo "bad format for mpi=$mpi and thr=$thr" >&2
    exit 1
fi

Nh=$((${BASH_REMATCH[1]}*${BASH_REMATCH[3]}))
Nv=$((${BASH_REMATCH[2]}*${BASH_REMATCH[4]}))

split=${Nh}x${Nv}

echo find $wdir -name "`basename $matrix .bin`.$split.????????.bin"
bfile=$(find $wdir -name "`basename $matrix .bin`.$split.????????.bin")
if ! [[ "$bfile" =~ \.[0-9a-f]{8}\.bin$ ]] ; then
    echo "Found no bfile, or nothing satisfactory ($bfile)" >&2
    exit 1
fi

cmd=`dirname $0`/convert_magma.pl


# This is for the **unbalanced** matrix !!

echo "m:=$m;n:=$n;interval:=$interval;nrhs:=$nrhs;" > $mdir/mn.m

echo "Saving matrix to magma format"
$cmd weights < $rwfile > $mdir/rw.m
$cmd weights < $cwfile > $mdir/cw.m
$cmd bpmatrix < $matrix > $mdir/t.m
$cmd balancing < "$bfile" > $mdir/b.m

placemats() {
    if [ "$nullspace" = left ] ; then
        transpose_if_left="Transpose"
    else
        transpose_if_left=""
    fi
    cat <<-EOF
        p:=$prime;
        nullspace:="$nullspace";
        xtr:=func<x|$transpose_if_left(x)>;
        M:=Matrix(GF(p),Matrix (M));
        nr:=Nrows(M);
        nc:=Ncols(M);
        nh:=$Nh;
        nv:=$Nv;
        nr:=nh*nv*Ceiling(Maximum(nr, nc)/(nh*nv));
        nc:=nr;
        x:=Matrix(GF(p),nr,nc,[]);InsertBlock(~x,M,1,1);M:=x;
        nrp:=nv*Ceiling (nr/(nh*nv));
        ncp:=nh*Ceiling (nc/(nh*nv));
        Mt:=Matrix(GF(p),nh*nrp,nv*ncp,[]);
EOF
    for i in `seq 0 $((Nh-1))` ; do
        cat <<-EOF
            nr$i:=nrp; /* nr div nh + ($i lt nr mod nh select 1 else 0); */
            snr$i:=$i*nrp; /* $i*(nr div nh) + Min($i, nr mod nh); */
EOF
        for j in `seq 0 $((Nv-1))` ; do
            $cmd bpmatrix < ${bfile%%.bin}.h$i.v$j> $mdir/t$i$j.m
            cat <<-EOF
                nc$j:=ncp; /*  div nv + ($j lt nc mod nv select 1 else 0); */
                snc$j:=$j*ncp; /* (nc div nv) + Min($j, nc mod nv); */
                load "$mdir/t$i$j.m";
                M$i$j:=Matrix(GF(p),Matrix(var));
                x:=RMatrixSpace(GF(p),nr$i,nc$j)!0;
                InsertBlock(~x,$transpose_if_left(M$i$j),1,1);
                M$i$j:=x;
                InsertBlock(~Mt,M$i$j,1+snr$i,1+snc$j);
EOF
        done
    done
    echo "mlist:=["
    for i in `seq 0 $((Nh-1))` ; do
        if [ "$i" != 0 ] ; then echo "," ; fi
        echo -n "["
        for j in `seq 0 $((Nv-1))` ; do
            if [ "$j" != 0 ] ; then echo -n ", " ; fi
            echo -n "M$i$j";
        done
        echo -n "]"
    done
    echo "];"
}

placemats > $mdir/placemats.m

echo "Saving vectors to magma format"
$cmd x $wdir/X > $mdir/x.m

print_all_rhs() {
    echo "RHS:=[VectorSpace(GF(p), nr_orig)|];"
    for j in `seq 0 $splitwidth $((nrhs-1))` ; do
        let j1=j+splitwidth
        $cmd spvector64 < $wdir/V${j}-${j1}.0 > $mdir/V${j}.0.m
        cat <<-EOF
            load "$mdir/V${j}.0.m";
            Append(~RHS, g(var));
EOF
    done
}
print_all_rhs > $mdir/R.m

print_all_v() {
    echo "$2:=[VectorSpace(GF(p), nr_orig)|];" 
    for j in `seq 0 $splitwidth $((n-1))` ; do
        let j1=j+splitwidth
        $cmd spvector64 < $wdir/V${j}-${j1}.$1 > $mdir/V${j}.$1.m
        cat <<-EOF
            load "$mdir/V${j}.${1}.m";
            Append(~$2, g(var));
EOF
    done
}
print_all_v 0 Y > $mdir/Y0.m

# When magma debug has been activated, we normally have the first 10
# iterates of each sequence.
echo "VV:=[];" > $mdir/VV.m
for i in {0..10} ; do
    print_all_v $i V_$i > $mdir/V.$i.m
    cat >> $mdir/VV.m <<-EOF
    load "$mdir/V.$i.m";
    Append(~VV, V_$i);
EOF
done

echo "Saving check vectors to magma format"
$cmd spvector64 < $wdir/C.0 > $mdir/C0.m
$cmd spvector64 < $wdir/C.1 > $mdir/C1.m
$cmd spvector64 < $wdir/C.$interval > $mdir/C$interval.m


echo "Saving krylov sequence to magma format"
afile=$(basename `find $wdir -name 'A*[0-9]'`)
$cmd spvector64 < $wdir/$afile > $mdir/A.m



echo "Saving linear generator to magma format"
$cmd spvector64 < $wdir/$afile.gen > $mdir/F.m
if [ -f $wdir/$afile.gen.rhs ] ; then
$cmd spvector64 < $wdir/$afile.gen.rhs > $mdir/rhscoeffs.m
else
    # We must create an empty file, because otherwise magma won't accept
    # "load" being conditional...
    echo > $mdir/rhscoeffs.m
fi

print_all_Ffiles() {
    echo "Fchunks:=[];";
    for s in `seq 0 $splitwidth $((n-1))` ; do
        let s1=s+splitwidth
        echo "Fs:=[];"
        for j in `seq 0 $splitwidth $((n-1))` ; do
            let j1=j+splitwidth
            $cmd spvector64 < $wdir/F.sols${s}-${s1}.${j}-${j1} > $mdir/F.$s.$j.m
            cat <<-EOF
                load "$mdir/F.$s.$j.m";
                Append(~Fs, g(var));
EOF
        done
        echo "Append(~Fchunks, Fs);"
    done
}
print_all_Ffiles > $mdir/Fchunks.m

# Ffiles=(`ls $wdir | perl -ne '/^F.sols\d+-\d+.\d+-\d+$/ && print;' | sort -n`)



echo "Saving mksol data to magma format"


echo "vars:=[];" > $mdir/S.m
k=$interval;
while true ; do
    for i in `seq 0 $((nrhs-1))` ; do
        let i1=i+1
        for j in `seq 0 $((n-1))` ; do
            let j1=j+1
            if ! [ -f $wdir/S.sols${i}-${i1}.${j}-${j1}.$k ] ; then
                if [ $i = 0 ] && [ $j = 0 ] ; then
                    echo "nblocks:=$((k/interval-1));" >> $mdir/S.m
                    break
                else
                    echo "Weird. Short of S files not at i,j=0 ?" >&2
                    exit 1
                fi
            fi
            $cmd spvector64 < $wdir/S.sols${i}-${i1}.${j}-${j1}.$k > $mdir/S${i},${j}.${k}.m
            echo "load \"$mdir/S${i},${j}.${k}.m\"; Append(~vars, var);" >> $mdir/S.m
        done
        if [ $i = 0 ] && [ $j = 0 ] ; then break ; fi
    done
    if [ $i = 0 ] && [ $j = 0 ] ; then break ; fi
    k=$((k+interval))
done

echo "Saving gather data to magma format"
echo "vars:=[];" > $mdir/K.m
(cd $wdir ; find  -name K.sols\*[0-9]) | while read f ; do
$cmd spvector64 < $wdir/$f > $mdir/$f.m
echo "load \"$mdir/$f.m\"; Append(~vars, var);" >> $mdir/K.m
done

echo "Running magma verification script"

s=$0
s=${s%%sh}m
cd $wdir
# magma does not exit with a useful return code, so we have to grep its
# output.
$magma -b $s < /dev/null | tee magma.out | grep -v '^Loading'
if egrep -q "(Runtime error|Assertion failed)" magma.out ; then
    exit 1
fi
