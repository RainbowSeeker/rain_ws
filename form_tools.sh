#! /usr/bin/bash

workdir=$(pwd)
remote_workdir="~/$(basename $workdir)"
source_cmd="source /opt/ros/humble/setup.bash && source $remote_workdir/install/setup.bash"

help() {
    echo "Usage: $0 [options] [arguments]"
    echo "Options:"
    echo "  -h, --help:                 Display this help message."
    echo "  -s, --sync:                 Sync the install directory in all machines."
    echo "  -u, --upload:               Upload all files to target machines except from .gitignore."
    echo "  -p, --prelaunch:            PreLaunch the formation on all machines."
    echo "  -l, --launch:               Launch the formation on all machines."
    echo "  -k, --kill:                 Kill the formation on all machines."
    echo "  -d, --disarm <obj_id>/all:  Disarm the obj_id or all machines."
}

every_run() {
    IFS=$'\n'
    comment=$1
    remote_exec=$2
    upper_cmd=$3
    amc_id=0
    for line in $(cat machines.txt)
    do
        machine=$(echo $line | awk '{print $1}')
        password=$(echo $line | awk '{print $2}')
        amc_id=$[amc_id+1]
        run_cmd=$(eval "echo \"$upper_cmd\"")
        echo "Perform [$comment] on machine [$machine]...$run_cmd"
        if [ "$remote_exec" = true ]; then
            run_cmd="tmux new-session -d -s '$comment' '$run_cmd'"
            sshpass -p "$password" ssh $machine "$run_cmd"
        else
            eval $run_cmd
        fi
    done
    IFS=
}

upload() {
    machine=$(cat machines.txt | awk '{print $1}' | head -n 1)
    password=$(cat machines.txt | awk '{print $2}' | head -n 1)
    run_cmd="rsync -az --delete --exclude-from='.gitignore' -e ssh ~/rain_ws $machine:~/"
    eval sshpass -p "$password" $run_cmd
    echo "Upload completed."
    run_cmd="$source_cmd; cd $remote_workdir && colcon build --packages-select mqsls"
    sshpass -p "$password" ssh $machine "$run_cmd"
}

sync() {
    every_run "sync" false "sshpass -p \$password rsync -az --delete -e ssh $workdir/install/ \$machine:$remote_workdir/install/"
}

prelaunch() {
    every_run "prelaunch" true "$source_cmd; ros2 launch mqsls actual_3dof_mqsls_prelaunch.py amc_id:=\$amc_id hil_enable:=False"
}

launch() {
    every_run "launch" true "$source_cmd; ros2 launch mqsls actual_3dof_mqsls_launch.py amc_id:=\$amc_id"
}

kill() {
    every_run "kill" true "tmux kill-server"
}

disarm() {
    obj_id=$1
    if [ "$obj_id" = "all" ]; then
        every_run "disarm" false "$source_cmd; ros2 run formation commander_node \$amc_id 400 0"
    elif [ "$obj_id" -ge 0 ]; then
        eval "$source_cmd; ros2 run formation commander_node $obj_id 400 0"
    else
        echo "Invalid obj_id."
    fi
}

while [ "$1" != "" ]; do
    case $1 in
        -h | --help )           help
                                exit
                                ;;
        -s | --sync )           sync
                                exit
                                ;;
        -u | --upload )         upload
                                exit
                                ;;
        -p | --prelaunch )      prelaunch
                                exit
                                ;;
        -l | --launch )         launch
                                exit
                                ;;
        -k | --kill )           kill
                                exit
                                ;;
        -d | --disarm )         disarm $2
                                exit
                                ;;                       
        * )                     help
                                exit 1
    esac
    shift
done