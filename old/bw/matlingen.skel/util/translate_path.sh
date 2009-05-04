#!/bin/bash

up_path() {
	X="`pwd`"
	Y="`echo $1 | sed -e 's,/\.$,,g'`"
	Z=""
	let spin=0
	while [ "$Y" != . ] ; do
		Z="`basename $X`/$Z"
		Y="`dirname $Y`"
		X="`dirname $X`"
		let spin+=1
		if [ $spin -ge 50 ] ; then
			echo "Bug in util/translate_path.sh !" >&2
			exit 1
		fi
	done
	Z="`echo $Z | sed -e s,/\$,,g`"
	if [ ! -n "$Z" ] ; then
		Z=.
	fi
	echo "$Z"
}

Z="$1/$(up_path "$2")"
echo "$Z" | sed -e 's,//*,/,g' -e 's,\./,,g' -e 's,/\.,,g'
