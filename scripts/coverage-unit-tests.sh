#!/bin/bash

set -e
export DIR=`mktemp -d /tmp/cado-cov.XXXXXXXXXXXX`
export CADO_DEBUG=1
export COV=1
F=`mktemp /tmp/cado-cov.local.XXXXXXXXX.sh`
cat > $F <<EOF
CFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"
CXXFLAGS="-O0 -g -fprofile-arcs -ftest-coverage"
EOF
export LOCALFILE=$F

if [ -x "$HOME/cmake/bin/cmake" ] ; then
    echo "Using $HOME/cmake/bin for cmake search path"
    PATH="$HOME/cmake/bin:$PATH"
fi

# Add $HOME/bin unconditionally, for convenience. It's easy to install
# e.g. git in some directory under $HOME, and then have a link from
# there.
PATH="$HOME/bin:$PATH"

cd $DIR

git clone git://scm.gforge.inria.fr/cado-nfs/cado-nfs.git src
SOURCETREE=src

cp $LOCALFILE "$SOURCETREE/local.sh"

CADO_DIST_ARCHIVE_NAME=cado-nfs
export CADO_DIST_ARCHIVE_NAME

# first build a tarball
make -C $SOURCETREE dist

# now extract the tarball and test it
tar zxf $SOURCETREE/$CADO_DIST_ARCHIVE_NAME.tar.gz

if [ "$LOCALFILE" ] ; then
    cp -f $LOCALFILE $CADO_DIST_ARCHIVE_NAME/local.sh
fi


cd $CADO_DIST_ARCHIVE_NAME

echo "Starting compilation at: `date`"
# This is a kludge.
touch files.dist
make cmake
# make -j8
echo "Starting tests at: `date`"
cat > $DIR/epilog.html <<EOF
<p style="text-align: center;">See also the <a href="./@basedir@/make-test.txt">test results</a></p>
EOF
maketest() {
export CTEST_OUTPUT_ON_FAILURE=yes
if make -j8 test $EXTRA_MAKETEST_ARGS; then
    :
else
    cat > $DIR/epilog.html <<EOF
<p style="text-align: center;">Warning: some tests have failed, see <a href="./@basedir@/make-test.txt">test results</a></p>
EOF
fi
}

maketest | tee $DIR/make-test.txt

# find . -name 'test_*.gcda' -print0 | xargs -0 -r rm
find $(find $DIR/$CADO_DIST_ARCHIVE_NAME/build -type d -name tests) -name '*.gcda' -print0 | xargs -0 -r rm


cd $DIR/$CADO_DIST_ARCHIVE_NAME
geninfo --no-checksum --ignore-errors gcov,source -q --output-filename $DIR/cado-nfs.info  ./ --no-external
rm -rf ~/.webdir/cado-unit-tests/ || :
/bin/cp -pf $(which genhtml) $DIR/genhtml
ex $DIR/genhtml <<EOF
/sub get_date_string()$
/{
mark a
/}
mark b
'a,'b c
{
       return scalar localtime;
}
.
wq
EOF
chmod 755 $DIR/genhtml
$DIR/genhtml --html-epilog $DIR/epilog.html  -o ~/.webdir/cado-unit-tests/ $DIR/cado-nfs.info
cp $DIR/make-test.txt ~/.webdir/cado-unit-tests
rm -rf $DIR $F

