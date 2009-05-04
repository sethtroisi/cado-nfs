#!/bin/sh

src="$1"
dst="$2"

copy="Copyright (C) 1999--`date +%Y` Emmanuel Thom'e --- see LICENSE file"

b=`basename $src`
if [ -d "$dst" ] ; then
	dst="$dst/$b"
fi

case "$b" in
	*.[ch]) (echo "/* $copy */" ; echo ; cat $src) > $dst ;;
	*.[ch]pp) (echo "/* $copy */" ; echo ; cat $src) > $dst ;;
	*.m4) (echo "dnl $copy" ; echo "dnl" ; cat $src) > $dst ;;
	*) echo "$b : not handled" ; exit 1;;
esac

