#!/usr/bin/env bash

set -e
set -x

# This script is intended *for testing only*. It's used for, e.g.,
# coverage tests. This even goes with dumping all intermediary data to
# magma-parseable format. So don't try to use this script for real
# examples. Real examples require the user to craft his own bwc.pl
# command line with care (command-lines as created by this tool might be
# a source of inspiration, though).

if [ "$*" ] ; then eval "$*" ; fi

# various configuration variables. environment can be used to override them
: ${scriptpath=$0}
: ${m=8}
: ${n=4}
: ${prime=4148386731260605647525186547488842396461625774241327567978137}
: ${Mh=1}
: ${Mv=1}
: ${Th=2}
: ${Tv=2}
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
: ${wipe=}
: ${seed=}

wordsize=64
# XXX note that $wdir is wiped out by this script !
: ${wdir=/tmp/bwcp}
: ${buildopts="MPI=1 DEBUG=1"}
: ${shuffle=1}
: ${nomagma=}

# This has to be provided as auxiliary data for the GF(p) case (or we'll
# do without any rhs whatsoever).
: ${rhs=}
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
    if ! [ "$prime" ] ; then
        prime=2
    fi
    if [ "$prime" = 2 ] ; then
        if [ $((m % 64)) != 0 ] || [ $((n % 64)) != 0 ] ; then
            echo "\$m and \$n must be multiples of 64" >&2
            usage
        fi
        if [ "$rhs" ] || [ "$nrhs" ] ; then
            echo "\$rhs and \$nrhs are unsupported in the GF(2) case" >&2
            usage
        fi
    else
        if [ "$rhs" ] && ! [ "$matrix" ] ; then
            echo "Please specify \$rhs only with \$matrix" >&2
            usage
        fi
        if [ "$matrix" ] && [ "$nrhs" ] ; then
            echo "Please specify \$nrhs only without \$matrix (so that a random problem gets generated)" >&2
            usage
        fi
    fi
}

derived_variables() {
    # Set here the variables which merit global scope, as they have
    # global relevance.
    Nh=$((Mh*Th))
    Nv=$((Mv*Tv))
    mpi=${Mh}x${Mv}
    thr=${Th}x${Tv}

    mpi_njobs_lingen=$((Nh*Nv))
    mpi_njobs_other=$((Mh*Mv))

    nwords=$(echo "1+l($prime)/l(2)/$wordsize"| bc -l | cut -d. -f1)
    n32bit=$(echo "1+l($prime)/l(2)/32"| bc -l | cut -d. -f1)
    bits_per_coeff=$((nwords*$wordsize))

    : ${plingen_program:=plingen_p_${nwords}}

    top=`dirname $0`/../..

    if ! [ -d $mats ] ; then mats=$matsfallback; fi
    if [ "$prime" = 2 ] ; then
        splitwidth=64
    else
        splitwidth=1
    fi

    sequences=()
    all_splits=0
    for j in `seq 0 $splitwidth $((n-1))`; do
        sequences=("${sequences[@]}" $j..$((j+splitwidth)))
        all_splits=$all_splits,$((j+splitwidth))
    done

    if [ "$bindir" ] ; then
        mpi_bindir=$(perl -ne '/HAVE_MPI\s*"(.*)"\s*$/ && print "$1\n";' $bindir/../../cado_mpi_config.h)
        if ! [ "$mpi_bindir" ] ; then
            disable_parallel_plingen=1
        elif ! [ -d "$mpi_bindir" ] ; then
            echo "HAVE_MPI variable in cado_mpi_config.h ($mpi_bindir) does not point to a directory" >&2
            exit 1
        fi
    else
        eval "export $buildopts"
        if ! [ "$MPI" ] || [ "$MPI" = 0 ] ; then
            disable_parallel_plingen=1
        elif [ "$MPI" = 1 ] ; then
            :
        elif [ -d "$MPI" ] ; then
            mpi_bindir="$MPI"
        else
            echo "MPI environment variable ($MPI) does not point to a directory (nor is it 0 or 1)" >&2
            exit 1
        fi
    fi

    if [ "$disable_parallel_plingen" ] ; then
        mpi_njobs_lingen=1
    fi
    if [ $mpi_njobs_lingen -gt 1 ] ; then
        if ! [ -d "$mpi_bindir" ] ; then
            echo "We need mpi_bindir to run parallel plingen" >&2
            exit 1
        fi
    fi
    if ! [ "$mpi_bindir" ] || [ -d "$mpi_bindir" ] ; then
        :
    else
        echo "mpi_bindir must be either empty or point to a directory" >&2
        exit 1
    fi

    echo "Working on ${#sequences[@]} sequences"

    if ! [ "$seed" ] ; then
        seed=$RANDOM
        echo "Setting seed to $seed"
    fi
}

optionally_build() {
    # build the software. This is optional, really. The only thing is that
    # it's so easy to mess with the location of the binaries that building it
    # by ourselves looks like a reasonable thing to do, for this illustrative
    # script.

    if ! [ "$bindir" ] ; then
        eval "export $buildopts"
        make -s -C $top -j 4
        make -s -C $top -j 4 $plingen_program
        eval `make -s -C $top show`
        bindir=$top/$build_tree/linalg/bwc
    fi
}

prepare_wdir() {
    if [ -d $wdir ] ; then
        if [ "$wipe" ] ; then
            rm -rf $wdir 2>/dev/null
        else
            echo "Won't wipe $wdir unless \$wipe is set" >&2
            exit 1
        fi
    fi
    mkdir $wdir
    echo $seed > $wdir/seed.txt
}


create_test_matrix_if_needed() {
    if [ "$matrix" ] ; then
        if ! [ -e "$matrix" ] ; then
            echo "matrix=$matrix inaccessible" >&2
            usage
        fi
        # Get absolute path.
        matrix=$(readlink /proc/self/fd/99 99< $matrix)
        if [ -e "$rhs" ] ; then
            rhs=$(readlink /proc/self/fd/99 99< $rhs)
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
    if [ "$prime" = 2 ] ; then
        $bindir/random_matrix  ${random_matrix_size} -s $seed --k$nullspace ${random_matrix_minkernel} > $mats/t${random_matrix_size}.txt
        matrix=$mats/t${random_matrix_size}.txt
    elif ! [ "$nrhs" ] ; then
        $bindir/random_matrix  ${random_matrix_size} -s $seed -c ${random_matrix_maxcoeff} --k$nullspace ${random_matrix_minkernel} > $mats/t${random_matrix_size}p.txt
        matrix=$mats/t${random_matrix_size}p.txt
    else
        # We use specialized magma code to create the example.
        basename=$mats/t${random_matrix_size}p+${nrhs}
        # I'm adding a density info here, because I've too often stumbled
        # across problems.
        density=`echo "l($random_matrix_size)^2/2" | bc -l | cut -d. -f1`
        if [ "$density" -lt 12 ] ; then density=12; fi
        $bindir/random_matrix $random_matrix_size -s $seed -d $density -c ${random_matrix_maxcoeff}  rhs=$nrhs,$prime,"${basename}".rhs.txt > "${basename}".matrix.txt
        # TODO: This will probably have to be tweaked, as we rather
        # expect the matrix to come with some expected excess.
        matrix="$basename.matrix.txt"
        rhs="$basename.rhs.txt"
    fi
}

create_auxiliary_weight_files() {
    if [ "$prime" != 2 ] ; then withcoeffs=--withcoeffs ; fi
    case "$matrix" in
        *.txt)
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
                $bindir/mf_scan  --binary-in $withcoeffs --mfile $matrix --binary-out --freq
            fi
            ;;
    esac
    ncols=$((`wc -c < $cwfile` / 4))
    nrows=$((`wc -c < $rwfile` / 4))
}



create_balancing_file_based_on_mesh_dimensions() {
    if [ "$shuffle" = 1 ] ; then
        shuffle_option=--shuffled-product
    else
        # This removes the de-correlating permutation.
        shuffle_option="noshuffle=1"
    fi

    if [ "$prime" != 2 ] ; then withcoeffs=--withcoeffs ; fi

    $bindir/mf_bal $shuffle_option mfile=$matrix $Nh $Nv out=$wdir/ $withcoeffs

    if [ $(ls $wdir/`basename $matrix .bin`.${Nh}x${Nv}.*.bin | wc -l) != 1 ] ; then
        echo "Weird -- should have only one balancing file as output." >&2
        exit 1
    fi
    bfile=$(ls $wdir/`basename $matrix .bin`.${Nh}x${Nv}.????????.bin)
    echo "Using balancing file $bfile"

    # The checksum file is used later on for the magma code.
    checksum=${bfile#$wdir/`basename $matrix .bin`.${Nh}x${Nv}.}
    checksum=`basename $checksum .bin`
}

prepare_common_arguments() {
    # This sets the common arguments for all bwc binaries
    read -s -r -d '' common <<-EOF
        matrix=$matrix
        mpi=$mpi
        thr=$thr
        balancing=$bfile
        m=$m
        n=$n
        wdir=$wdir
        prime=$prime
        mm_impl=$mm_impl
        nullspace=$nullspace

EOF
}

convert_rhs_text_to_matrix() {
    if ! [ -r "$rhs" ] ; then
        echo "$rhs unreadable" >&1
        usage
    fi
    (set +x
    rhs="$1"
    set `head -2 $rhs | tail -n 1`
    nrhs=$#
    set `head -1 $rhs`
    if [ $# = 1 ] ; then set "$@" $nrhs ; fi
    echo "$@"
    tail -n +2 $rhs | tr ':' ' ' | while read a; do
        set $a
        echo -n "$#"
        j=0
        while [ $# -gt 0 ] ; do
            echo -n " $j:$1" 
            let j+=1
            shift
        done
        echo
    done)
}

create_binary_rhs_matrix() {
    echo "Running perl now" >&2
    # This is really done in perl.
    read -s -r -d '' code <<-'EOF'
        use warnings;
        use strict;
        use Fcntl;
        print STDERR "Create binary RHS started\n";
        sysopen(STDIN,'&STDIN',O_RDONLY);
        binmode(STDIN, ':bytes');
        my ($nrhs,$n32b,$prefix,$group) = @ARGV;
        $group = 1 unless $group;
        my @fds;
        for(my $i = 0 ; $i < $nrhs; $i+=$group) {
            my $j = $i + $group;
            my $fh;
            my $fname=$prefix . "${i}-${j}.0";
            print STDERR "Opening $fname as file number $i\n";
            sysopen($fh, $fname, O_CREAT | O_WRONLY) or die "$fname: $!";
            push @fds, $fh;
        }
        my $rowidx=0;
        while(sysread(STDIN, my $x, 4)) {
            $rowidx++;
            my $v = unpack("L", $x);
            die unless $v == $nrhs;
            for(my $i = 0 ; $i < $nrhs; $i+=$group) {
                for(my $k = 0 ; $k < $group ; $k++) {
                    sysread(STDIN, $x, 4) or die;
                    $x=unpack("L", $x);
                    die unless $x == $i + $k;
                    my $todo = 4 * $n32b;
                    for(my $done = 0; $done < $todo; ) {
                        my $delta = sysread(STDIN, $x, $todo - $done);
                        syswrite($fds[$i], $x, $delta);
                        $done += $delta;
                    }
                    # Do this because we expect 64-bit data !
                    if ($n32b & 1) {
                        $x="\0\0\0\0";
                        syswrite($fds[$i], $x, 4);
                    }
                }
            }
        }
        for(my $i = 0 ; $i < $nrhs; $i+=$group) {
            close($fds[$i]);
        }
        print STDERR "Create binary RHS ok\n";

EOF
    set -e
    perl -e "$code" "$@"
}

argument_checking
derived_variables
optionally_build
prepare_wdir

create_test_matrix_if_needed
# This also sets rwfile cwfile nrows ncols
create_auxiliary_weight_files
create_balancing_file_based_on_mesh_dimensions

prepare_common_arguments

if [ "$rhs" ] ; then
    # http://stackoverflow.com/questions/4489139/bash-process-substitution-and-syncing
    # $bindir/mf_scan  --ascii-in --with-long-coeffs $n32bit --mfile <(convert_rhs_text_to_matrix $rhs)  --binary-out --ofile >(create_binary_rhs_matrix $nrhs $n32bit ${rhs}.bin $nrhs)
    # fortunately mf_scan is stdout-silent...
    $bindir/mf_scan  --ascii-in --with-long-coeffs $n32bit --mfile <(convert_rhs_text_to_matrix $rhs)  --binary-out --ofile - | create_binary_rhs_matrix $nrhs $n32bit ${rhs}.bin $nrhs
    rhsbin=`echo $rhs | sed -e s/txt/bin/`
    if [ "$rhsbin" == "$rhs" ] ; then rhsbin=$rhs.bin ; fi
    mv ${rhs}.bin0-$nrhs.0 $rhsbin
    set $common rhs=$rhsbin nrhs=$nrhs
else
    set $common
fi
$bindir/bwc.pl :complete "$@"

if [ "$nomagma" ] ; then
    exit 0
fi

## Now take solution 0, for instance.
#j0=0
#while [ $j0 -lt $n ] ; do
#    let j1=$j0+1
#    ln -s F.sols0-1.${j0}-${j1} $wdir/F${j0}-${j1}
#    j0=$j1
#done


#j0=0
#while [ $j0 -lt $n ] ; do
#    let j1=$j0+1
#    $bindir/bwc.pl mksol  $common interval=$interval ys=$j0..$j1 nsolvecs=$nrhs
#    j0=$j1
#done

set +x

# [ "$?" = 0 ] && $bindir/bwc.pl acollect    $common -- --remove-old
# 
# 
# 
# # set +e
# # $bindir/bench  --nmax 1000 --prime $prime --nbys 1 --impl basicp  -- $matrix
 
mdir=$wdir

split=${Nh}x${Nv}
b=`basename $matrix .bin`
c=$b.$split.$checksum

cmd=`dirname $0`/convert_magma.pl


# This is for the **unbalanced** matrix !!

echo "m:=$m;n:=$n;" > $mdir/mn.m

echo "Saving matrix to magma format"
$cmd weights < $rwfile > $mdir/rw.m
$cmd weights < $cwfile > $mdir/cw.m
$cmd bpmatrix < $matrix > $mdir/t.m
$cmd balancing < $wdir/$c.bin > $mdir/b.m

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
            $cmd bpmatrix < $wdir/$c.h$i.v$j> $mdir/t$i$j.m
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
$cmd spvector64 < $wdir/$afile > $mdir/A.m



echo "Saving linear generator to magma format"
$cmd spvector64 < $wdir/$afile.gen > $mdir/F.m
if [ -f $wdir/$afile.gen.rhs ] ; then
$cmd spvector64 < $wdir/$afile.gen.rhs > $mdir/rhscoeffs.m
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



echo "Saving mksol and gather data to magma format"


if [ -f $wdir/K.sols0-1.0 ] ; then
$cmd spvector64 < $wdir/K.sols0-1.0 > $mdir/K.sols0-1.0.m
fi

echo "vars:=[];" > $mdir/S.m

for i in `seq 0 $((nrhs-1))` ; do
    let i1=i+1
    for j in `seq 0 $((n-1))` ; do
        let j1=j+1
        $cmd spvector64 < $wdir/S.sols${i}-${i1}.${j}-${j1}.$interval > $mdir/S${i},${j}.m
        echo "load \"$mdir/S${i},${j}.m\"; Append(~vars, var);" >> $mdir/S.m
    done
done

