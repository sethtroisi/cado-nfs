#!/usr/bin/env bash

OUTPUT_FILE="$1"

if ! [ -d .git ] ; then exit 0 ; fi

function list_modified() {
  echo '#define CADO_MODIFIED_FILES "\'
  git status --porcelain -uno | while read STATUS FILE
  do
    if [[ "$STATUS" = M  ||  "$STATUS" = A ]]
    then
      SHA1="`sha1sum "$FILE" | cut -d " " -f 1`"
      echo "# $STATUS" "$FILE" "$SHA1"'\n\'
    else
      echo "# $STATUS" "$FILE"'\n\'
    fi
  done
  echo '"'
}

TEMPFILE="`mktemp /tmp/XXXXXXX`"
list_modified >> "$TEMPFILE"
if ! [ -f "$OUTPUT_FILE" ] || ! cmp -s "$TEMPFILE" "$OUTPUT_FILE"
then
  rm -f "$OUTPUT_FILE"
  cp "$TEMPFILE" "$OUTPUT_FILE"
fi
rm -f "$TEMPFILE"
