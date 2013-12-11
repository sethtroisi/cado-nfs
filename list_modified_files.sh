#!/usr/bin/env bash

echo '#define CADO_MODIFIED_FILES "\'
git status --porcelain -uno | while read STATUS FILE
do
  SHA1="`sha1sum "$FILE" | cut -d " " -f 1`"
  echo "# $STATUS" "$FILE" "$SHA1" '\n\'
done
echo '"'
