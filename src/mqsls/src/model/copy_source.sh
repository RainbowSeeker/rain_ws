#! /usr/bin/bash

current_dir=$(dirname $0)
source_dir=~/MQSLS_Model/build/payload_controller_ert_rtw
source_files=$(ls $source_dir | grep -E '*\.(cpp|h)$')

rm -rf $current_dir/controller/*
for file in $source_files
do
    cp $source_dir/$file $current_dir/controller/
done