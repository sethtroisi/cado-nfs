#!/usr/bin/env bash

OUTPUT_FILE="$1"

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

list_modified >| "$OUTPUT_FILE"
