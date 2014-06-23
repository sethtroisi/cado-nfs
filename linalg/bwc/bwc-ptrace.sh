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
# do without any SM whatsoever).
: ${sm=}
: ${nullspace=right}
: ${interval=50}
: ${mm_impl=basicp}


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
        if [ "$sm" ] || [ "$nsm" ] ; then
            echo "\$sm and \$nsm are meaningless in the GF(2) case" >&2
            usage
        fi
    else
        if [ "$sm" ] && ! [ "$matrix" ] ; then
            echo "Please specify \$sm only with \$matrix" >&2
            usage
        fi
        if [ "$matrix" ] && [ "$nsm" ] ; then
            echo "Please specify \$nsm only without \$matrix (so that a random matrix gets generated)" >&2
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

    top=`dirname $0`/../..

    if ! [ -d $mats ] ; then mats=$matsfallback; fi
    if [ "$prime" = 2 ] ; then
        splitwidth=64
    else
        splitwidth=1
    fi

    sequences=()
    all_splits=0
    j0=0
    while [ $j0 -lt $n ] ; do
        sequences=("${sequences[@]}" $j0..$((j0+splitwidth)))
        let j0=$j0+$splitwidth
        all_splits=$all_splits,$j0
    done

    if [ "$bindir" ] ; then
        if grep -q 'undef HAVE_MPI' $bindir/../../cado_mpi_config.h ; then
            disable_parallel_plingen=1
        fi
    else
        eval "export $buildopts"
        if ! [ "$MPI" ] || [ "$MPI" = 0 ] ; then
            disable_parallel_plingen=1
        fi
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
        make -s -C $top -j 4 plingen_p_$nwords
        # make -s -C $top -j 4 plingen_pz
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

create_example_dlp_matrix_with_sm() {
    # This creates the four files $basename.{matrix,sm}{,2}.txt
    # magma is required to create the test data.
    read -s -r -d '' code <<-EOF
        M:=Matrix(Integers(),N+nsm,N,[]);
        for i in [1..N+nsm] do for j in [1..12] do k:=Random([1..N]); M[i,k]+:=Random(5)-2; end for; end for;
        F:=Open(basename cat ".matrix.txt", "w");
        fprintf F, "%o %o\n", N, N;
        for i in [1..N] do
            nz:=[j:j in [1..N]|M[i,j] ne 0];
            fprintf F, "%o %o\n", #nz, Join([Sprintf("%o:%o", j-1,M[i,j]):j in nz], " ");
        end for;
        delete F;

        F:=Open(basename cat ".matrix2.txt", "w");
        fprintf F, "%o %o\n", nsm, N;
        for i in [N+1..N+nsm] do
            nz:=[j:j in [1..N]|M[i,j] ne 0];
            fprintf F, "%o %o\n", #nz, Join([Sprintf("%o:%o", j-1,M[i,j]):j in nz], " ");
        end for;
        delete F;

        sm:=[Random(VectorSpace(GF(p),N+nsm)): j in [1..nsm-1]];
        /* adjust last SM so that there is a kernel mod p */
        v0:=Vector(GF(p),[1..N]);

        v:=v0*Transpose(Matrix(GF(p),M));
        for j in [1..nsm-1] do
            v+:=(N+j)*sm[j];
        end for;
        Append(~sm, v/-(N+nsm));

        assert Rank(Matrix(GF(p),M)) eq N;
        Mfull:=HorizontalJoin(Matrix(GF(p),M),Transpose(Matrix(sm)));
        assert IsZero(Vector(GF(p),[1..N+nsm])*Transpose(Mfull));

        F:=Open(basename cat ".sm.txt", "w");
        fprintf F, "%o\n", N;
        for i in [1..N] do
            fprintf F, "%o\n", Join([IntegerToString(Integers()!sm[j][i]):j in [1..nsm]]," ");
        end for;
        delete F;

        F:=Open(basename cat ".sm2.txt", "w");
        fprintf F, "%o\n", nsm;
        for i in [N+1..N+nsm] do
            fprintf F, "%o\n", Join([IntegerToString(Integers()!sm[j][i]):j in [1..nsm]]," ");
        end for;
        delete F;

EOF
    (
    echo "N:=$random_matrix_size; p:=$prime; nsm:=$nsm;"
    echo "basename:=\"$basename\";"
    echo "$code"
    ) | magma -b
}

create_test_matrix_if_needed() {
    if [ "$matrix" ] ; then
        if ! [ -e "$matrix" ] ; then
            echo "matrix=$matrix inaccessible" >&2
            usage
        fi
        # Get absolute path.
        matrix=$(readlink /proc/self/fd/99 99< $matrix)
        if [ -e "$sm" ] ; then
            sm=$(readlink /proc/self/fd/99 99< $sm)
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
        basename=$mats/t${random_matrix_size}p+${nsm}sm
        create_example_dlp_matrix_with_sm
        # TODO: This will probably have to be tweaked, as we rather
        # expect the matrix to come with some expected excess.
        matrix="$basename.matrix.txt"
        sm="$basename.sm.txt"
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
        # I have the impression that it's better to leave this here !
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

dispatch_matrix() {
    # ys=0..$splitwidth here is really a hack. It merely has to match the
    # version which is used in production.
    set $common ys=0..$splitwidth
    submats="save_submatrices=1"
    if ! [ "$nomagma" ] ; then set "$@" save_submatrices=1 ; fi
    if [ "$prime" = 2 ] ; then set "$@" sanity_check_vector=H1 ; fi
    $bindir/bwc.pl dispatch "$@"

}


convert_sm_text_to_matrix() {
    if ! [ -r "$sm" ] ; then
        echo "$sm unreadable" >&1
        usage
    fi
    (set +x
    sm="$1"
    set `head -2 $sm | tail -n 1`
    nsm=$#
    echo "`head -1 $sm` $nsm" 
    tail -n +2 $sm | while read a; do
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

create_binary_sm_matrix() {
    echo "Running perl now" >&2
    # This is really done in perl.
    set +e
    read -s -r -d '' code <<-'EOF'
        use warnings;
        use strict;
        use Fcntl;
        print STDERR "Create binary SM started\n";
        sysopen(STDIN,'&STDIN',O_RDONLY);
        binmode(STDIN, ':bytes');
        my ($nsm,$n32b,$prefix) = @ARGV;
        my @fds;
        for(my $i = 0 ; $i < $nsm; $i++) {
            my $j = $i + 1;
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
            die unless $v == $nsm;
            for(my $i = 0 ; $i < $nsm; $i++) {
                sysread(STDIN, $x, 4) or die;
                $x=unpack("L", $x);
                die unless $x == $i;
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
        for(my $i = 0 ; $i < $nsm; $i++) {
            close($fds[$i]);
        }
        print STDERR "Create binary SM ok\n";

EOF
    perl -e "$code" "$@"
}


do_prep_step_gfp() {
    if [ "$sm" ] ; then
        set `head -2 $sm | tail -n 1`
        nsm=$#
        if [ "$n" -lt "$nsm" ] ; then
            echo "Error: $sm has $nsm SMs, so we need a block size n>=$nsm (currently have n=$n)" >&2
            exit 1
        fi
        echo "Using $nsm vectors from SM matrix"
        # We create a binary matrix, this will ease things significantly.
        $bindir/mf_scan  --ascii-in --with-long-coeffs $n32bit --mfile <(convert_sm_text_to_matrix $sm)  --binary-out --ofile >(create_binary_sm_matrix $nsm $n32bit $wdir/V)
    else
        nsm=0
    fi
    j0=$nsm

    if [ $n -eq 1 ] && [ $j0 -eq 0 ] ; then
        # This should be the fallback to the "normal" case...
        $bindir/bwc.pl prep   $common
        ln -s Y.0 $wdir/V0-1.0
    else
        # prep won't work. Let's be stupid. We use the SM provided, if any,
        # and later use auto-generated vectors. We will shift the things a
        # bit in lingen so as to get a proper generators (we'll do A(X) div X
        # for the columns corresponding to our random vectors here).
        if [ "$nullspace" = right ] ; then
            nbytes=$((bits_per_coeff/8 * $ncols))
        else
            echo "Untested" >&2
            exit 1
            nbytes=$((bits_per_coeff/8 * $nrows))
        fi
        if [ "$nsm" -gt 0 ] ; then
            echo "Padding $nsm vectors from SM matrix with $((n-nsm)) random vectors"
        fi
        while [ $j0 -lt $n ] ; do
            let j1=j0+1
            echo "Generating $wdir/V${j0}-${j1}.0"
            dd if=/dev/urandom bs=1 count=$nbytes of=$wdir/V${j0}-${j1}.0
            j0=$j1
        done
        (echo 1 ; seq 0 $((m-1))) > $wdir/X
        # TODO: create Y, too.
        # XXX Does it make any sense, after all ??
    fi
}


do_prep_step_gf2() {
    $bindir/bwc.pl prep $common seed=$seed
    $bindir/bwc.pl :ysplit $common splits=$all_splits
}

do_prep_step() {
    if [ "$prime" = 2 ] ; then
        do_prep_step_gf2
    else
        # This merges prep + ysplit, effectively removing the need for a Y
        # vector.
        do_prep_step_gfp
    fi
}

do_limited_krylov_step_just_for_magma() {
    # Just for playing. The output of this is used later on for checking with
    # magma.
    if ! [ "$nomagma" ] ; then
        $bindir/bwc.pl secure  $common interval=1
        for seq in "${sequences[@]}" ; do
            $bindir/bwc.pl krylov  $common interval=1 end=$interval ys=$seq skip_online_checks=1
        done
        rm -f $wdir/A*
    fi
}

do_krylov_step() {
    # Now do krylov for real.
    $bindir/bwc.pl secure  $common interval=$interval
    for seq in "${sequences[@]}" ; do
        $bindir/bwc.pl krylov  $common interval=$interval ys=$seq
    done
}

do_lingen_step_gfp() {
    afile=$($bindir/acollect wdir=$wdir m=$m n=$n bits-per-coeff=$bits_per_coeff --remove-old | tail -1)

    if [ $mpi_njobs_lingen -gt 1 ] ; then
        precmd="mpirun -n $mpi_njobs_lingen"
    fi
    if ! [ "$disable_parallel_plingen" ] ; then
        $precmd $bindir/plingen_p_$nwords lingen-mpi-threshold=10000 lingen-threshold=10 m=$m n=$n wdir=$wdir prime=$prime afile=$afile mpi=$mpi thr=$thr nsm=$nsm
    else
        $bindir/plingen_p_$nwords lingen-mpi-threshold=10000 lingen-threshold=10 m=$m n=$n wdir=$wdir prime=$prime afile=$afile nsm=$nsm
    fi

    ln $wdir/$afile.gen $wdir/F

    # This does the splitting as documented in mksol.c, e.g. with F on disk
    # stored as the transpose of the reversal of the F in A*F=G+O(X^t).

    # TODO: We should get rid of splitting altogether, it's ridiculous.
    # Better do everything within lingen/plingen.

    # Note that the transpose was absent in commits 4f7d835 and earlier.
    $bindir/split wdir=$wdir m=$m n=$n splits=$all_splits     \
        ifile=F ofile-fmt=F.%u-%u --binary-ratio $((bits_per_coeff/8))/1
    for seq in "${sequences[@]}" ; do
        dseq=$(echo $seq | sed -e 's/\.\./-/g')
        $bindir/split wdir=$wdir m=$m n=$n              \
            splits=$all_splits                          \
            ifile=F.$dseq ofile-fmt=F.sols%u-%u.$dseq   \
            --binary-ratio $((bits_per_coeff/8))/1
    done

}

do_lingen_step_gf2() {
    $bindir/bwc.pl acollect    $common -- --remove-old
    $bindir/bwc.pl lingen      $common lingen_threshold=64
    # This keeps *all* solutions, but splits the data into several
    # chunks, each appropriate for one of the sequences.
    $bindir/bwc.pl :fsplit     $common splits=$all_splits
    # We now also need a second split. But it's significantly more
    # difficult. We have 64 coeffs to skip.
    #
    # TODO: We should get rid of splitting altogether, it's ridiculous.
    # Better do everything within lingen/plingen.
    if [ "${#sequences[@]}" -gt 1 ] ; then
        for seq in "${sequences[@]}" ; do
            dseq=$(echo $seq | sed -e 's/\.\./-/g')
            $bindir/split wdir=$wdir m=$m n=$n              \
                splits=$all_splits                          \
                ifile=F.sols0-$n.$dseq ofile-fmt=F.sols%u-%u.$dseq   \
                --binary-ratio $splitwidth/8
        done
    fi
#    for seq in "${sequences[@]}" ; do
#        dseq=$(echo $seq | sed -e 's/\.\./-/g')
#        ln -s F$dseq $wdir/F.sols0-$n.$dseq
#    done
}

do_lingen_step() {
    if [ "$prime" = 2 ] ; then
        do_lingen_step_gf2
    else
        do_lingen_step_gfp
    fi
}

do_mksol_step() {
    set $common interval=$interval
    if [ "$nsm" ] && [ "$nsm" -gt 0 ] ; then
        set "$@" nsolvecs=$nsm
    else
        set "$@" nsolvecs=$splitwidth
    fi
    for seq in "${sequences[@]}" ; do
        $bindir/bwc.pl mksol ys=$seq "$@"
    done
}

do_gather_step() {
    set $common interval=$interval nsolvecs=$splitwidth
    if [ "$nsm" ] && [ "$nsm" -gt 0 ] ; then
        if ! [ -f $wdir/$afile.gen.rhs ] ; then
            echo "sm-outside-matrix requires the file $wdir/$afile.gen.rhs to be computed by plingen" >&2
            exit 1
        fi
        set "$@" nsm=$nsm inhomogeneous_rhs=$wdir/$afile.gen.rhs
    fi
    set +e
    $bindir/bwc.pl gather "$@"
    rc=$?
    set +e
    if [ "$rc" != 0 ] ; then
        if [ "$nsm" ] && [ "$nsm" -gt 0 ] ; then
            echo "Accepting failure within \"gather\" for now" >&2
        else
            echo "\"$bindir/bwc.pl gather\" failed, exit code $?" >&2
            exit 1
        fi
    fi
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
dispatch_matrix

do_prep_step

do_limited_krylov_step_just_for_magma

do_krylov_step

do_lingen_step

do_mksol_step

do_gather_step

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
#    $bindir/bwc.pl mksol  $common interval=$interval ys=$j0..$j1 nsolvecs=$nsm
#    j0=$j1
#done















set +x

# [ "$?" = 0 ] && $bindir/bwc.pl acollect    $common -- --remove-old
# 
# 
# 
# # set +e
# # $bindir/bench  --nmax 1000 --prime $prime --nbys 1 --impl basicp  -- $matrix
# 
mdir=$wdir

split=${Nh}x${Nv}
b=`basename $matrix .bin`
c=$b.$split.$checksum

if [ "$nomagma" ] ; then
    exit 0
fi

cmd=`dirname $0`/convert_magma.pl

$cmd weights < $rwfile > $mdir/rw.m
$cmd weights < $cwfile > $mdir/cw.m

# This is for the **unbalanced** matrix !!

echo "m:=$m;n:=$n;" > $mdir/mn.m

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


if [ -f  $wdir/Y.0  ] ; then
    $cmd spvector64 < $wdir/Y.0 > $mdir/Y0.m
else
    :
    # otherwise we're probably in the case where n>1
fi
for j in `seq 0 $((n-1))` ; do
    let j1=$j+1
    $cmd spvector64 < $wdir/V${j}-${j1}.0 > $mdir/V${j}.0.m
    $cmd spvector64 < $wdir/V${j}-${j1}.1 > $mdir/V${j}.1.m
    $cmd spvector64 < $wdir/V${j}-${j1}.2 > $mdir/V${j}.2.m
    $cmd spvector64 < $wdir/V${j}-${j1}.3 > $mdir/V${j}.3.m
    $cmd spvector64 < $wdir/V${j}-${j1}.4 > $mdir/V${j}.4.m
    $cmd spvector64 < $wdir/V${j}-${j1}.5 > $mdir/V${j}.5.m
    $cmd spvector64 < $wdir/V${j}-${j1}.6 > $mdir/V${j}.6.m
    $cmd spvector64 < $wdir/V${j}-${j1}.7 > $mdir/V${j}.7.m
    $cmd spvector64 < $wdir/V${j}-${j1}.8 > $mdir/V${j}.8.m
    $cmd spvector64 < $wdir/V${j}-${j1}.9 > $mdir/V${j}.9.m
    $cmd spvector64 < $wdir/V${j}-${j1}.10 > $mdir/V${j}.10.m
done
$cmd spvector64 < $wdir/C.0 > $mdir/C0.m
$cmd spvector64 < $wdir/C.1 > $mdir/C1.m
$cmd spvector64 < $wdir/C.$interval > $mdir/C$interval.m
$cmd x $wdir/X > $mdir/x.m
$cmd spvector64 < $wdir/$afile > $mdir/A.m
$cmd spvector64 < $wdir/$afile.gen > $mdir/F.m
if [ -f $wdir/$afile.gen.rhs ] ; then
$cmd spvector64 < $wdir/$afile.gen.rhs > $mdir/rhs.m
fi
if [ -f $wdir/K.sols0-1.0 ] ; then
$cmd spvector64 < $wdir/K.sols0-1.0 > $mdir/K.sols0-1.0.m
fi

echo "vars:=[];" > $mdir/S.m

for i in `seq 0 $((nsm-1))` ; do
    let i1=i+1
    for j in `seq 0 $((n-1))` ; do
        let j1=j+1
        $cmd spvector64 < $wdir/S.sols${i}-${i1}.${j}-${j1}.$interval > $mdir/S${i},${j}.m
        echo "load \"$mdir/S${i},${j}.m\"; Append(~vars, var);" >> $mdir/S.m
    done
done

