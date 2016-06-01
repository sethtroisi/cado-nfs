#!/bin/bash

nfields=10
primelist="2 3 5 7 11 13 17 19 23"

while [ $# -gt 0 ] ; do
    if [ "$1" = "--binary" ] ; then
        binary=$2
        shift
    elif [ "$1" = "--magma" ] ; then
        magma=$2
        shift
    elif [ "$1" = "--fields" ] ; then
        nfields=$2
        shift
    elif [ "$1" = "--primelist" ] ; then
        primelist=$2
        shift
    else
        echo "unexpected arg: $1" >&2
        exit 1
    fi
    shift
done

if ! [ "$binary" ] || ! [ -x "$binary" ] ; then
    echo "binary=$binary: not found or not executable" >&2
    exit 1
fi
if ! [ "$magma" ] || ! [ -x "$magma" ] ; then
    echo "magma=$magma: not found or not executable" >&2
    exit 1
fi

input_file_for_binary() {
    d=$1
    shift
    echo $d
    echo "$@"
    for i in `seq 1 $d` ; do
        echo -n 1
        for j in `seq 1 $d` ; do
            [ $i != $j ]
            echo -n " $?"
        done
        echo
    done
}

read -s -r -d '' perl_code <<-'EOF'
    while (<>) {
        /the maximal order/ && do { $x=1; next; };
        s/^\s*//;
        print if $x;
    }

EOF

magma_code() {
cat <<-'EOF'
    f:=Polynomial(StringToIntegerSequence(COEFFS));
    n:=Degree(f);
    K<alpha>:=NumberField(f);
    basis:=Matrix(Rationals(),n,n,[StringToRational(x):x in Split(RESULT," ")]);
    O:=Order([LeadingCoefficient(f)*alpha]);
    try
        myO:=Order([FieldOfFractions(O)!Eltseq(r):r in Rows(basis)]);
        test := Gcd(Index(pMaximalOrder(myO,p),myO), p) eq 1;
        printf "ok p:=%o; d:=%o; COEFFS:=\"%o\";\n", p, n, COEFFS;
    catch e
        printf "FAILED p:=%o; d:=%o; COEFFS:=\"%o\"; RESULT:=\"%o\";\n", p, n, COEFFS, RESULT;
    end try;
    quit;
EOF
}


computed_basis() {
    $binary "$@" | perl -e "$perl_code" | xargs echo
}

primes=($(echo $primelist))

for n in `seq 1 $nfields` ; do
    d=$((RANDOM % 10))
    let d+=2
    coeffs=()
    x=()
    for i in `seq 0 $d` ; do
        if [ $i -lt $((d+1)) ] ; then
            x=("${x[@]}" $((RANDOM-10000)))
        else
            while true ; do
                lc=$RANDOM
                if [ $lc != 0 ] ; then break; fi
            done
            x=("${x[@]}" $lc)
        fi
    done
    set -- "${x[@]}"

    for p in "${primes[@]}" ; do
        res=$(computed_basis <(input_file_for_binary $d "${x[@]}") $p)
        set -- "${x[@]}"
        displayed=$($magma -b <(echo "COEFFS:=\"$*\"; RESULT:=\"$res\"; p:=$p;") <(magma_code))
        echo $displayed
        if echo $displayed | grep -q FAILED ; then
            exit 1
        fi
    done
done
