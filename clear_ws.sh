#! /usr/bin/bash

file_to_remove=$(ls -t install/*.csv | tail -n +11)

IFS=$'\n'
for line in $file_to_remove
do
    rm "$line"
done
IFS=

rm -rf log/*