#!/bin/sh

input="$1"

if [ ! -f "$input" ] ; then echo "$1 : no such file" ; exit 1 ; fi

pol=`< "$input" tr -d c:  | awk ' /^[0-9]/ { print $2 "*x^" $1,"+"; }' | xargs echo | sed -e 's/ *+$//g'`

alim=`< $input awk '/^alim:/ { print $2; }'`

echo "rootfind($pol, $alim)" | gp -q -p $alim rootfind.gp
