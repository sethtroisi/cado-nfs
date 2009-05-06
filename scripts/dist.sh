#!/bin/sh

# Usage: scripts/dist.sh <package basename>
# This creates <package basename>.tar.gz

set -e

pkg=$1
shift

t=`mktemp -d /tmp/XXXXXXXX`

if ! scripts/check_file_lists.pl ; then
    echo "Refusing to make the distribution. Sorry" >&2
    exit 1
fi

read_files_dist() {
    # This command reads the files.dist list. All nodes ending in / imply
    # that the whole subdirectory is matched.
    grep '^[^#]' files.dist | while read x ; do
        case $x in
            */) find $x ;;
            *) echo $x;;
        esac
    done | sed -e 's+//*+/+g' -e 's+^+/+'
    # On some machines, find foo/ will produce foo//file, which does not
    # play well with rsync.
}

read_files_dist | rsync -a     \
                    --include-from=-    \
                    --include='*/'      \
                    --exclude='*'       \
                    --delete-excluded   \
                    --prune-empty-dirs  \
                    ./ $t/$pkg/

(cd $t ; tar cf - $pkg) | gzip -9 > $pkg.tar.gz

echo "Built $pkg.tar.gz"
rm -rf $t
