#! /bin/sh

# Extracts text from docx file

# This line from http://woozle.org/~neale/papers/docx.html
unzip -qc "$1" word/document.xml | sed 's#</w:p>#\n\n#g;s#<[^>]*>##g' > $1.txt

# Convert form dos and change the file extension
for file in *.docx.txt
do
    dos2unix $file
    mv ${file} ${file%.docx.txt}.txt
done


