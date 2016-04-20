#!/bin/sh

# Names of new files
name=""

# Rename files s1-19, clipping off first page
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
do
  pdftk s$i.pdf cat ~1 output a$i.pdf
  name="$name a$i.pdf"
done

# Combine files
pdftk $name cat output Combined.pdf
