#!/bin/sh

CADO_NFS_SOURCE_DIR=$1
CADO_NFS_BINARY_DIR=$2

t=`mktemp -d /tmp/cado-check.XXXXXXX`

POLYFILE=$t/p3dd15-f4g3-GJL-1.poly
PARAMFILE=$t/p3dd15-f4g3-GJL-1.params
WDIR=$t/p3dd15

mkdir $WDIR

cat > $POLYFILE <<EOF
n: 100000000010189
# ell: 10000000002037900000103825911
skew: 1.000
Y0: 23501753462
Y1: 17389930810
Y2: -9893995512
Y3: 1259659241
c0: 1
c1: 3
c2: 4
c3: -2
c4: 1
# p = n = 5 mod 24
# f := x^4 - 2*x^3 + 4*x^2 + 3*x + 1;
# sgn(f) = (0, 2) => rk = 1
# smfexp:=100000000040758000006229554630423173648064579819794979920; #
# ell^2-1
# alpha(f) = 1.31698
# Murphy E (f) = 1.86762e-3
# g := 1259659241*x^3 - 9893995512*x^2 + 17389930810*x + 23501753462;
# sgn(g) = (1, 1) => rk = 1
# smgexp:=1000000000611370000155738865621158682343169972959968560744281687967628329623524533030;
# # ell^3-1
# alpha(g) = 1.49903
# Murphy E (g) = 5.46052e-4
# Murphy E (f, g) = 2.28411e-7
# log_2(Nf) = 2
# log_2(Ng) = 35
# varphi = x^3 + 36917133770091*x^2 + 3855640591067*x + 17939575058717
# varphi = gcd(f, g) mod p.
EOF


cat > $PARAMFILE <<EOF
name = p3dd15-f4g3-GJL-1
dlp = true
N = 100000000010189
gorder = 10000000002037900000103825911

slaves.nrclients = 2
tasks.threads = 2
tasks.execpath = $CADO_NFS_BINARY_DIR
slaves.scriptpath = $CADO_NFS_SOURCE_DIR/scripts/cadofactor
tasks.workdir = $WDIR
slaves.basepath= $WDIR/client
slaves.hostnames = localhost

tasks.polyselect.import = $POLYFILE

# this should be uncommented for nominal execution, since the computations
# involve 1 unit on side 1, 1 SM on side 0
tasks.explicit_units1 = true 
tasks.addfullcol = true

tasks.I = 13
tasks.polyselect.degree = 4
tasks.polyselect.admax = 0
tasks.polyselect.incr = 60
tasks.polyselect.adrange = 500
tasks.polyselect.P = 420
tasks.polyselect.nq = 1000

rlim = 50000
alim = 50000
lpbr = 17
lpba = 17
tasks.sieve.mfbr = 17
tasks.sieve.mfba = 17
tasks.sieve.qrange = 2000
tasks.sieve.rels_wanted = 30000
tasks.reconstructlog.partial = false
checkdlp = false
EOF

cleanup() {
    if ! [ "$CADO_DEBUG" ] ; then
        rm -rf $t
    else
        echo "(debug mode, temporary files are kept in $t)"
    fi
}

${CADO_NFS_SOURCE_DIR}/scripts/cadofactor/cadofactor.py $PARAMFILE && cleanup
