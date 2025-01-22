#! /usr/bin/bash

current_dir=$(dirname $0)
source_dir=~/MQSLS_Model/build/codegen/lib/mso_forces
source_files=$(ls $source_dir | grep -E '*\.(c|h)$')

rm -rf $current_dir/force_opt/*
for file in $source_files
do
    cp $source_dir/$file $current_dir/force_opt/
done