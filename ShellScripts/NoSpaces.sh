#!/bin/sh
#
# Remove all spaces from the files in the current dir
for f in *;
do
  mv "$f" "${f//[[:space:]]}"
done
