#!/bin/bash
#
# Remove all spaces from the files in the current dir
for i in * ; do
  if [ "$i" != ${i//[[:space:]]} ] ;
  then
    mv "$i" "${i//[[:space:]]}"
  fi
done
