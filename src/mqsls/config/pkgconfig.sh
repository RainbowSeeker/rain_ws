#! /usr/bin/bash
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <model_dir>"
    exit 1
fi
model_dir=$1
python3 $model_dir/rope/gen_rope.py

file_to_remove=$(ls -t install/*.csv | tail -n +11)

IFS=$'\n'
for line in $file_to_remove
do
    rm "$line"
done
IFS=

rm -rf log/*