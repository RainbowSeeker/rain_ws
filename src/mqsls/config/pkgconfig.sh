#! /usr/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <model_dir>"
    exit 1
fi
# generate rope
model_dir=$1
python3 $model_dir/rope/gen_rope.py

# kill old process
function kill_old_process {
    process_name=$1
    pid=$(ps -ef | grep $process_name | grep -v grep | awk '{print $2}')
    if [ -n "$pid" ]; then
        kill -9 $pid
    fi
}
kill_old_process "gz sim"
kill_old_process "px4"

# rm old files
file_to_remove=$(ls -t install/*.csv | tail -n +11)
IFS=$'\n'
for line in $file_to_remove
do
    rm "$line"
done
IFS=