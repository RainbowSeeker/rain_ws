#! /usr/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <model_dir>"
    exit 1
fi
# generate rope
model_dir=$1
python3 $model_dir/rope/gen_rope.py

# kill old process
pkill -9 --full "gz sim"
pkill -9 "px4"
pkill -9 "gz_plugin_node"
pkill -9 "MicroXRCEAgent"

# rm old files
file_to_remove=$(ls -t install/*.csv | tail -n +11)
IFS=$'\n'
for line in $file_to_remove
do
    rm "$line"
done
IFS=